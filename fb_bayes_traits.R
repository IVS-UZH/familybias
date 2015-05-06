library(familybias)
library(parallel)
library(devtools)
source_url('https://raw.githubusercontent.com/IVS-UZH/phylo-convert/master/phylo-convert.R')
source_url('https://github.com/IVS-UZH/familybias/raw/master/R/familybias.R') 

# the interface function
familybias_bayestraits <- function(df, family.names, r.name, bayes_traits_cmd ='./BayesTraitsV2', tmpdir='bayes_traits_data', ...) {
  
  buildResults <- function(current.df, varname, fbias.pre.correction, fbias.ml.correction, fbias.mcmc.correction)
  {
    list(
      fbias.pre.correction = fbias.pre.correction,
      fbias.ml.correction = fbias.ml.correction,
      fbias.mcmc.correction = fbias.mcmc.correction
    )
  }
  
  .processVariable <- function()
  {
   # prepare for BayesTraits correction
   bayes_traits_data <- current.df
   names(bayes_traits_data)[match(r.name, names(bayes_traits_data))] <- '.VARIABLE.'
   # prepare the list of units we are interested in processing 
   units <- droplevels(subset(fbias.pre.correction$large.families.estimate, majority.prop<1 | is.na(majority.prop), select=c('family.name','taxonomic.level')))

   # run BayesTraits and gather the stats
   bayes_traits_output <- .runBayesTraitsAndParseResults(bayes_traits_data, units)
 
   # build the override tables for the ML model
   bias_overrides_ml <- lapply(names(bayes_traits_output), function(unit.name) {
  	 # extract the stats
  	 stats <- bayes_traits_output[[unit.name]]
  	 # the family is biased of if the ML p value is <=0.05
  	 # then we take the majority values computed by BayesTraits
  	 if (stats$lr.p.ml <= 0.05) 
  			stats$maj.value.data.ml
  	 else
  			list(mj.prop=NA, mj.val=NA)
   })
   names(bias_overrides_ml) <- names(bayes_traits_output)

   # build the override tables for the MCMC model
   bias_overrides_mcmc <- lapply(names(bayes_traits_output), function(unit.name) {
  	 # extract the stats
  	 stats <- bayes_traits_output[[unit.name]]
  	 # the family is biased of if the bayes factor is >2
  	 # then we take the majority values computed by BayesTraits
  	 if ((stats$lr.mcmc > 2)) 
  		 stats$maj.value.data.mcmc
  	 else
  		 list(mj.prop=NA, mj.val=NA)
   })
   names(bias_overrides_mcmc) <- names(bayes_traits_output)
 
   # run the family bias estimation (with BayesTraits ML correction)
   fbias.ml.correction <- invoke_familybias(df=current.df, bias.override = c_override_lists(args$bias.override, bias_overrides_ml))

   # run the family bias estimation (with BayesTraits MCMC correction)
   fbias.mcmc.correction <- invoke_familybias(df=current.df, bias.override = c_override_lists(args$bias.override, bias_overrides_mcmc))
 
   buildResults(current.df, r.name, fbias.pre.correction, fbias.ml.correction, fbias.mcmc.correction)	
  }
  
  # runs the BayesTraits algorithm for the data
  .runBayesTraitsAndParseResults <- function(bayes_traits_data, units)
  {	
  	# prepare the trees
  	suppressWarnings({
  		trees <- as.phylo(bayes_traits_data, levels=attr(bayes_traits_data, 'taxa'), roots=units)
  		if(!'multiPhylo' %in% class(trees)) trees <- list(trees)
  	})
	
  	# init the path
  	unlink(tmpdir,recursive=T)
  	dir.create(paste(tmpdir, 'results', sep='/'), recursive=T, showWarnings=F)
	
  	# write out the tree and data for every unit
  	# run BayesTraits and parse the input back
  	mclapply(seq_along(trees), function(i)
  	{
  		# extract the appropriate data 
  		current.tree <- trees[[i]]
  		unitname <- current.tree$node.label[1]
  		lowest_taxon <- attr(bayes_traits_data, 'taxa')[[length(attr(bayes_traits_data, 'taxa'))]]
  		current.data <- bayes_traits_data[bayes_traits_data[[lowest_taxon]] %in% current.tree$tip.label, c(lowest_taxon, '.VARIABLE.')]

  		# format the data so that BayesTraits can read it
  		# there can be no node labels in trees
  		current.tree$node.label <- NULL
  		# and the trees should be cleaned up
  		current.tree <- collapse.singles(current.tree) 
  		# recode the variable levels as letters
  		VARIABLE.recoding <- data.frame(original = unique(current.data$.VARIABLE.))
  		VARIABLE.recoding$recoded <- factor(LETTERS[1:nrow(VARIABLE.recoding)])
  	  current.data$.VARIABLE. <- VARIABLE.recoding$recoded[match(current.data$.VARIABLE., VARIABLE.recoding$original)]

      # output the files
      output_file_path <- paste0(tmpdir,'/','UNIT', i)
  		write.nexus(current.tree, file=paste0(output_file_path, ".tree"))
  		write.table(current.data, file=paste0(output_file_path, ".data"), row.names = F, quote = F, col.names = F, sep = "\t")
		
  		# write out the command files
  		cfile <- function(.) paste0(tmpdir, '/commands.UNIT', i, '.', ., '.txt')
  		write(c('1','1','RestrictAll qAB','run'),sep='\n', file=cfile('ER'))
  		write(c('1','2','RestrictAll qAB','run'),sep='\n', file=cfile('ER.MCMC'))
  		write(c('1','1','run'),sep='\n', file=cfile('ARD'))
  		write(c('1','2','run'),sep='\n', file=cfile('ARD.MCMC'))	
   	 
  		# run all the commands and parse the output		
  		rfile <- function(.) paste0(tmpdir, '/results/result.UNIT', i, '.', ., '.txt')
  		output <- list()
  		for(command in c('ER', 'ER.MCMC', 'ARD', 'ARD.MCMC'))
  		{	
  			exe <- paste(bayes_traits_cmd, 
  			             paste0(output_file_path, ".tree"),
  									 paste0(output_file_path, ".data"), 
  									  '<', cfile(command), '>', rfile(command))
  			system(exe)
			
  			tryCatch(
  				output[[command]] <- parseBayesTraitOutput(rfile(command)),
  				error = function(e)
  				{
  					stop('Unable to parse BayesTraits results')
  				})
			
  	  }		
		
  		# helper func to determine the majority value for the BayesTraits results
  		find_majority_value <- function(tab) 
  		{
  			# get the transition probabilities
  			# its the last row of the table and the columns are on form qXY
  			trans_prob <- unlist(tab[nrow(tab), grepl('^q[A-Z][A-Z]$', names(tab)), drop=T])			
  			# we want to look for the highest transition rate
  			max.transition.rate <- names(trans_prob)[which.max(trans_prob)]
			
  			# and the majority value is last letter
  			maj.val <- substr(max.transition.rate, nchar(max.transition.rate), nchar(max.transition.rate))
			
  			# recode it back and return it together with the proportion
  			maj.val.original <- as.character(VARIABLE.recoding$original[match(maj.val, VARIABLE.recoding$recoded)])			
  			list(mj.prop=sum(as.character(current.data$.VARIABLE.) %in% as.character(maj.val))/nrow(current.data), mj.val=maj.val.original)
  		}
		
		
  		# gather the result statistics
  		stats <- list(
  			# 1. the ML likelihood ratio and majority value of the ARD model
  			lr.ml=2*(output$ARD$Lh-output$ER$Lh),
  			maj.value.data.ml = find_majority_value(output$ARD),
  			# 2. the bayes factor of MCMC evaluation and the majority value of the 
  			# ARD.MCMC model
  			lr.mcmc=2*(output$ARD.MCMC$Harmonic.Mean[nrow(output$ARD.MCMC)]-output$ER.MCMC$Harmonic.Mean[nrow(output$ER.MCMC)]),
  			maj.value.data.mcmc = find_majority_value(output$ARD.MCMC)
  		)	
  		stats$lr.p.ml <- 1 - pchisq(stats$lr.ml, 2*choose(length(unique(current.data$.VARIABLE.)),2)-1)
		
		
  		# return the stats
  		stats 
  	}) -> results
	
  	names(results) <- sapply(trees, function(.) .$node.label[1])
	
  	# check if there are any errors
  	errors <- sapply(results, inherits, 'try-error')
  	if(any(errors))
  	{
  		errors <- results[errors]
  		errors <- sapply(errors, function(.) attr(., 'condition')$message)
  		errors <- paste(names(errors), errors, sep=':', collapse='; ')
  		stop(errors)
  	}
	
	
  	results
  }


  parseBayesTraitOutput <- function(path)
  {
  	lines <- readLines(path)
	
  	# extract only the table 
  	# we use a very simple heuristic here: tables contain at least 5 tabs
  	tbl <- lines[grepl('(\t[^\t]*){5,}', lines, perl=T)]
  	tbl <- gsub('[[:space:]]+$', '', tbl)
	
  	read.table(textConnection(tbl), header=T, sep='\t')
  }
  
  
  c_override_lists <- function(x, y) {
    if(is.null(x)) return(y)
    if(is.null(y)) return(x)
      
    for(i in seq_along(y)) x[[names(y)[i]]] <- y[[i]]  
      
    x  
  }

  invoke_familybias <- function(...) {
    # build the familybias call 
    # we remove all arguments that are not relevant to familybias
    # + override other arguments
    fb_args <- args[setdiff(names(args), c('bayes_traits_cmd', 'tmpdir'))] 
      
    do.call(familybias, c_override_lists(fb_args, list(...)))
  }
  
  # CODE STARTS HERE
  
  
  # get the argument list
  args <- as.list(match.call())[-1]

  current.df <- df

	# we need to ensure that there are no duplicate names in the taxonomy. Best do
	# it now, so that things can be consistent. I will just go radical
	# here and rename ALL taxa, this is the simplest way. Also, I will make sure
	# that all taxa names are compatible with whatever BayesTraits wants. Its
	# ugly, but will hopefully allow us to avoid hard to detect bugs later on
	for(l in family.names[-length(family.names)])
	{
		current.df[[l]] <- ifelse(!is.na(current.df[[l]]), paste(current.df[[l]], l), NA)
	}

	# check the tips
	lowest_taxon <- family.names[[length(family.names)]]

	# make sure that there are no NAs on the lowest level
  stopifnot(!any(is.na(current.df[[lowest_taxon]])))

	# similarly, make sure that every tip has a unique name and that there are no
	# spaces in the language names as it appears to throw BayesTraits off	
	current.df[[lowest_taxon]] <- gsub(' +', '', current.df[[lowest_taxon]])
  
  
  
	# chek if there are duplicates on the lowest level, if yes, we want to
	# introduce subsystems
	if(any(duplicated(current.df[[lowest_taxon]])))
	{
		lowest_taxon1 <- paste0(lowest_taxon, '-subsystem')
		current.df[[lowest_taxon1]] <- factor(paste(current.df[[lowest_taxon]], 1:nrow(current.df), sep='-'))
		attr(current.df, 'taxa') <- c(family.names, lowest_taxon1)
	}
	else attr(current.df, 'taxa') <- family.names
    
		
	# the inner processor might crash
	# if it crashes, we simply report an error		
	r <- buildResults(current.df, r.name, NULL, NULL, NULL)
  
  r <- 'stuff'

	result <- tryCatch(
		{
		  # run the family bias estimation (before BayesTraits correction)
      
	    fbias.pre.correction <- invoke_familybias(df=current.df)
	
		  r <- buildResults(current.df, r.name, fbias.pre.correction, NULL, NULL)
		  .processVariable()
		},
		error = function(e)
		{
			r$problems <- as.character(e$message)
		
			r
		})
	
	result	
}



mean.familybias <- function(x, ...)
{
	varfields <- grep("^number\\.", names(x$units), value=T)
	varlevels <- c('diverse', sub("^number\\.( )*", "", varfields))	
	
	# now we want to build the dataframe
	hdr <- expand.grid(rev(c(lapply(x$predictors,  function(pred) unique(as.character(x$units[, pred]))), list(as.character(varlevels)))))
	hdr <- hdr[, rev(1:ncol(hdr)), drop=F]
	names(hdr) <- c(x$predictors, 'majority.response')
	hdr0 <- do.call("paste", c(hdr, sep = "\r"))
	
	
	Freq = numeric(length=nrow(hdr))
	
  
  if(length(x$extrapolations)>0) {
	for(ex in x$extrapolations) 
	{
		ex = do.call("paste", c(ex[, names(hdr), drop=F], sep = "\r"))
		for(i in seq_along(hdr0))
			Freq[i] = Freq[i] + sum(ex %in% hdr0[i])
    
    # cat('---\n')
    # print(ex)
    # print(Freq)
    # stop()
	}
	
	hdr$Freq <- Freq/length(x$extrapolations)
  } else
  {
    ex <- x$large.families.estimate
    
		ex = do.call("paste", c(ex[, names(hdr), drop=F], sep = "\r"))
		for(i in seq_along(hdr0))
			Freq[i] = Freq[i] + sum(ex %in% hdr0[i])
    
    hdr$Freq <- Freq
     
  }
	
	hdr
}


# do not run
data(vprel.g)
head(vprel.g)

taxa.names = c("stock", "mbranch", "sbranch", "ssbranch", "lsbranch", "language")


x <- familybias_bayestraits(df = vprel.g, family.names = taxa.names,
                        r.name = 'DRYREL0', p.names = 'DRYSOV4',
                        extrapolate=F, small.family.size=4, B=50)
#
mean(x$fbias.pre.correction)