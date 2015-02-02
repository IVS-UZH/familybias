multiple.binom.test <- function(x, p.threshold=0.1)
{
	# print(x)
	x <- table(x)
	# print(x)
	#
	# stop()
	
	# do pair-wise binom tests
	p.vals <- sapply(seq_along(x), function(i) binom.test(x[i], sum(x), p=1/length(x), alternative="greater")$p.value)
	names(p.vals) <- names(x)
	
	p.vals <- sort(p.vals, decreasing=F)
	
	# perform Holm-Bonferroni correction
	fits <- p.vals <= sapply(length(x):1, function(k) p.threshold/k)

	
	nr <- which(!fits)[1]
	if(!is.na(nr))
		fits[nr:length(fits)] <- F


	confirmed <- names(fits[fits])	

		
	if(length(confirmed) != 1)
		return(NA) 
	else
 		return(confirmed)	
}


familybias <- function(df, family.names, r.name, response_is_continous = FALSE, p.names=NULL, extrapolate=TRUE, 
                       B=1000, small.family.size=4, diverse.r = "diverse", 
                       bias.test = multiple.binom.test, p.threshold = 0.1,
                       verbose=F, lapplyfunc=lapply, bias.override = list())
{
	# setup the result list
	result <- structure(list(), class="familybias")
	result$predictors <- p.names
	
	response_is_continous <- response_is_continous || !(is.factor(df[[r.name]]) || is.character(df[[r.name]]))
	
	# remove the factors from the data frame
	t <- df
	df <- as.data.frame(as.matrix(df), stringsAsFactors=F)
	names(df) <- names(t)

	# remove all the entries which contain an NA either in the response or the predictor
	# chunk User-level interface:families():4
	df <- df[!apply(
			is.na(df[, c(r.name, p.names), drop=F]), 
			1, 
			any),
	]
	
	
	# ensure that the lowest-level names are unique
	if(any(duplicated(df[, family.names[length(family.names)]])))
	{
		for(nn in unique(df[, family.names[length(family.names)]]))
		if(sum(df[, family.names[length(family.names)]] %in% nn)>1)
			df[df[, family.names[length(family.names)]] %in% nn, family.names[length(family.names)]] <- paste(df[df[, family.names[length(family.names)]] %in% nn, family.names[length(family.names)]], 1:sum(df[, family.names[length(family.names)]] %in% nn), sep=".")

		
		warning(paste("Duplicates found at the lowest taxonomic level '", 
		 		family.names[length(family.names)],
		         "'", sep=""))
	}
	
	
	# cthe response variable may not be a constant (TODO: let the user specify the response levels)
	if(length(unique(df[, r.name]))==1)
		stop(paste("The response variable", r.name, "is a constant!"))
		
	# split the units
	if(verbose) cat("Splitting the data...\n")
	result$units <- split.families(df, family.names, r.name, p.names, verbose, response_is_continous)
	
	
	# estimate the bias of large families
	if(verbose) cat("Estimating the bias of large language units...\n")
	flush.console()
	result$large.families.estimate <- test.families(result$units, small.family.size, diverse.r, bias.test, verbose, lapplyfunc, p.threshold, bias.override)

	#return(result$large.families.estimate)
		
	
	# combine units and large.families.estimate
	result$units.estimate <- 
	{	
		small <- subset(result$units, !family.name %in% result$large.families.estimate$family.name)

		if(nrow(small) ==0) result$large.families.estimate 
		else
		{

			small$majority.response <- NA
			small$majority.prop <- NA
			small$distribution <- "small"

			rbind(result$large.families.estimate, small)
		}
	}	
	
	# cperform extrapoation if required
	if (extrapolate)
		result <- structure(c(
				result,
		 		extrapolate.families(result, p.names, small.family.size, B,  diverse.r, verbose, lapplyfunc)),
			    class="familybias")
	
	# we are done here!
	result
}


split.families <- function(DF, family.names, r.name, p.names, verbose, response_is_continous=FALSE)
{
	# insert dummy names for non-specified levels
	if(length(family.names)>1)
		for(taxon in (length(family.names)-1):1)
			DF[, family.names[taxon]] <-  
			ifelse(is.na(DF[, family.names[taxon]]),
		       	DF[, family.names[taxon+1]], 
		       	DF[, family.names[taxon]])
		
	# split the languages acording to their families (highest taxonomic level)
	units <-split(DF, DF[, family.names[1]], drop=T)
		  
	# this is the recursive splitting function
	# it will attempt to split the languages into taxonomic groups 
	# such that the predictor values are constant within each group
	# if this is impossible, pseudo groups will be postulated
	process.unit <- function(unit)
		# if the predictor values are constant withign the group, we have
		# reacehd the required split depth
		if (nrow(unique(unit[, p.names, drop=F]))==1) 
			list(unit)
		# otherwise, we must split further 
		else	
		{
			# find the highest not-yet-split level
			notsplit <- which(apply(
					unit[, family.names],
					2,
					function(x) length(unique(x)))>1)[1]
					
					
			# further splitting would result in singletons
			# - we create pseudo groups here by splittign the unit accordign to predictor values
			if (length(unique(unit[, family.names[notsplit]])) == nrow(unit))
			{
				pseudogroup.level <- family.names[notsplit-1]
				unit[, pseudogroup.level] <- paste(unit[, pseudogroup.level], "pseudo-group")
			
				pseudogroups <- split(unit, as.list(unit[, p.names, drop=F]), drop=T)
								
				for(i in seq_along(pseudogroups))
				{
					pseudogroups[[i]][, pseudogroup.level] = paste(pseudogroups[[i]][, pseudogroup.level], i, sep='')
				}
				
				return(pseudogroups)
				#return(split(unit, as.list(unit[, p.names, drop=F]), drop=T))
			}
			# otherwise, we should split further, recursively applying the process.unit function 
			else
			return(do.call(c, lapply(
					split(unit, unit[, family.names[notsplit]], drop=T), 
				 	process.unit)))
		} 
	# end process.unit
	
	# if predictors are present, we must go through the recursive splitting procedure	 
	if (!is.null(p.names))	
		units <- do.call(c, lapply(units, process.unit))	
	# otherwise, we need only to split according to teh families (which we have already done)
	
	
	# convert the result to an easy-to-work-with data frame
	# 	one unit at a time
	result <- do.call(rbind, lapply(units, function(unit)
	{
		# determine at what actual taxonomic level the unit was split
		split.level <- which(apply(
				unit[, rev(family.names), drop=F],
				2,
				function(x) length(unique(x)))==1)[1]
		split.level <- rev(family.names)[split.level]
		
		response_is_continous <- F
		if(!response_is_continous)
		{
		# one unit corresponds to a singel row in the resulting data frame
		result <- as.data.frame(c(list(
			# family.name
			as.character(unit[1, split.level]), 
			# taxonomic.level
			split.level), 
			# p.names	 
			lapply(unlist(unit[1, p.names]), as.character),
			# response value counts		 
			lapply(unique(DF[, r.name]),  function(x) sum(unit[, r.name] %in% x)), 
			# total number of languages
			list(nrow(unit))
		))
		# ...and set the column names
		names(result) <- c(
   			"family.name",
  			"taxonomic.level",
   			p.names,
   			paste("number.", unique(DF[, r.name]), sep=""),
   		    "number"
		)
	  }
		else
		{
		# one unit corresponds to a singel row in the resulting data frame
		result <- as.data.frame(c(list(
			# family.name
			as.character(unit[1, split.level]), 
			# taxonomic.level
			split.level), 
			# p.names	 
			lapply(unlist(unit[1, p.names]), as.character),
			# total number of languages
			list(nrow(unit))
		))
		# ...and set the column names
		names(result) <- c(
   			"family.name",
  			"taxonomic.level",
   			p.names,
   		    "number"
		)
	  }	
	
		if (response_is_continous)
			result$raw_responce <- I(list(unit[[r.name]]))
		else
			result$raw_responce <- I(list(factor(unit[[r.name]], levels=unique(DF[[r.name]]))))

		
		result 		
	}))
	
	# clean up the row names (R has a bad habbit of giving them horrible labels)
	rownames(result) <- NULL

	result
}




test.families <- function(DF, small.family.size, diverse.r, bias.test, verbose, lapplyfunc=lapply, p.threshold, bias.override = list())
{
	# find the names of the columns which store the frequency data (of the response variable)
	# varfields <- grep("^number\\.", names(DF), value=T)
	# varlevels <- sub("^number\\.", "", varfields)
	# <----- we move to raw_responce now!

	# we can only test "large" families, so we extract them from the data
	DF <- subset(DF, number>small.family.size)
	
	if(verbose) pb <- txtProgressBar(min=1, max=nrow(DF), initial=1, style=3)
	
	# for each family, perform the bias test
	result <- cbind(DF,  do.call(rbind, lapplyfunc(1:nrow(DF), function(row)
	{
		# test the frequency distribution
		vars <- unlist(DF[row, 'raw_responce'])
		
		family.name <- as.character(DF$family.name[[row]])
				
		if (family.name %in% names(bias.override))
			bias.val <- bias.override[[family.name]]$bias.val
		else
			bias.val <- bias.test(vars, p.threshold = p.threshold)
		
		
		if(verbose) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
		
		# determine the majority (modal) type
		mj.val <- bias.val

		# determine the proportion of the majority type in the family
		# plug in the override
		if (family.name %in% names(bias.override))
		{
			mj.prop <- as.numeric(bias.override[[family.name]]$mj.prop)
			# test if mj.val is in there as well
			if('mj.val' %in% names(bias.override[[family.name]]))
				mj.val <- bias.override[[family.name]]$mj.val
		}
		else
		{
			# if(family.name == 'Semitic stock')
			# {
			# 	print(mj.val)
			# 	print(vars)
			# 	stop()
			# }
			#
			mj.prop <- if (is.na(mj.val)) NA else sum(vars %in% mj.val)/length(vars)
		}
		
		distribution <- if(!is.na(mj.val)) "biased" else "diverse"
		if(is.na(mj.val)) mj.val <- diverse.r
		
		
		if (!is.numeric(mj.prop) || is.nan(mj.prop))
			mj.prop <- NA

		# return it as data frame
		data.frame(
			majority.response=mj.val,
			majority.prop=mj.prop,
			distribution=distribution, 
			row.names=NULL)
	})))
	
	if(verbose) close(pb)
	# we are done here
	result
}




extrapolate.families <- function(dta,  p.names, small.family.size, B, diverse.r, verbose, lapplyfunc=lapply)
{
	if(verbose) cat("Performing small family extrapolation (", B, " simulations)\n", sep="")
	flush.console()
	
	
	
	# find the names of the columns which store the frequency data (of the response variable)
	varfields <- grep("^number\\.", names(dta$units), value=T)
	varlevels <- sub("^number\\.", "", varfields)

	# extract the small units and add columns for estimation parameters
	dta$units <- subset(dta$units, number <= small.family.size)
	if(nrow(dta$units)==0)
	{
		warning("Cannot perform extrapolation: no small families present!\n")
		return(NULL)
	}	
	
	dta$units$majority.response <- diverse.r
	dta$units$majority.prop <- NA
	dta$units$distribution <- NA

	# collapse the estimated large families with (yet to be estimated) small families
	dta <- rbind(dta$units, dta$large.families.estimate)
	
	# and split the data according to predictors
	split.data <- if(is.null(p.names)) 
					list(dta)
		 		 else
					split(dta, 
				  		  as.list(dta[, p.names, drop=F]), 
				  		  drop=T)

	# this function is used to compute the prior probabilities
	compute.prior <- function(taxa.list)
	if(nrow(taxa.list)>0)
		data.frame(
			# prior probability of bias (laplace estimate)
			p.bias = with(taxa.list, 
				(sum(distribution %in% "biased") + 1)/
				(length(distribution) + 2)),
			# and its 95% credibility interval
			p.bias.95.cred = I(list(with(taxa.list,
				round(qbeta(c(0.025,0.975),
			  	sum(distribution %in% "biased") + 1,
			  	sum(!distribution %in% "biased") + 1),
			 	digits=2)))),
			# prior probability of deviation
	p.deviate = if(any(taxa.list$distribution %in% "biased"))
				mean(1 - taxa.list$majority.prop, na.rm=T)	
			  else 0.05
	# no data, just assume some defaults
	) else data.frame(
		p.bias 	= 0.5,
		p.bias.95.cred = I(list(c(0,1))),
		p.deviate = 0.05
	)	
	

	# compute the prior probabilities
	prior.data <- do.call(rbind, lapply(split.data, function(X)
		if(is.null(p.names))
		 	compute.prior(subset(X, number > small.family.size))
		else cbind(
				X[1,p.names, drop=F],
				compute.prior(subset(X, number > small.family.size)))
	))

	rownames(prior.data) <- NULL
	
	if (verbose) pb <- txtProgressBar(min=1, max=B, initial=1, style=3)
	
	# do the simulations B times
	extrapolations <- lapplyfunc(1:B, function(simulation_id)
	{
	if(verbose) setTxtProgressBar(pb, simulation_id)
	
	# which are done separetely for each predictor group (if any)
	do.call(rbind, lapply(seq_along(split.data), function(group_id)
	# and combines the already known estimation parameters for large families
	rbind(subset(split.data[[group_id]], number > small.family.size), 
	# with following simulated estimates for small families
	{
		# extract the small families
		units <- subset(split.data[[group_id]], number <= small.family.size)
		pr.data <- prior.data[group_id, ]


		# no data, nothign to estimate
		if(nrow(units)==0) return(subset(split.data[[group_id]], number > small.family.size))
		
		# first, assign all small families the label "diverse"
		units$distribution <- "diverse"
		units$majority.response <- diverse.r
		
		# assign a part of the data the label biased:representative 
		units$distribution[as.logical(rbinom(nrow(units), 1, pr.data$p.bias))] <-  'biased:representative'
		
		
		# part of these units are deviates
		units$distribution[units$distribution %in% "biased:representative"][as.logical(rbinom(sum(units$distribution %in% "biased:representative"), 1, pr.data$p.deviate))] <-  'biased:deviate'
		
		# compute the majority value
		for(i in grep("biased", units$distribution))
		{
			m <- units[i, varfields]
			m <- varlevels[m == max(m)]
			units$majority.response[i] <- sample(m, 1)
		}
		
		# adjust for deviates
		for(i in grep("deviate", units$distribution))
			units$majority.response[i] <- sample(setdiff(varlevels, units$majority.response[i]), 1)
		
		# compute the majority proportion
		for(i in grep("biased", units$distribution))
		{
			m <- units[i, varfields]
			names(m) <- varlevels
			units$majority.prop[i] <- m[units$majority.response[i]]/sum(m)
		}

		# correct the labels
		units$distribution[grep("biased", units$distribution)] <- "biased"

		# we are done here
		units
	}
	)))})
	if(verbose) close(pb)
	
	# return the estimates
	list(prior=prior.data, extrapolations=extrapolations)
}

# chunk FamilyBias class printing:4
make.freq <- function(X, x, pretty=F)
{
	varfields <- grep("^number\\.", names(x$units), value=T)
	varlevels <- sub("^number\\.( )*", "", varfields)
	
	B <- length(X)
	
	X <- do.call(rbind, X)
	
	# a very weird R bug - randomly inserts initial spaces???
	X$majority.response <- gsub("(^ +)|( +$)", "", X$majority.response)
	
	X <- if(is.null(x$predictors)) 
			list(X) 
		else
	 		split(X, as.list(X[, x$predictors, drop=F]), drop=T)
	
	X <- do.call(rbind, lapply(X, function(XX)
	{ 
	
		varlevels.f <- sapply(varlevels, function(varlevel)
		{
			f <- sum(XX$majority.response %in% varlevel)									
			f <- if(pretty)
			{
				f0 <- c(f/B, f/nrow(XX)*100)
				f0 <- format(f0,digits=3,drop0trailing=T)
				paste(f0[1], " (", paste(f0[2], "%", sep=""), ")", sep="")
			} 
				else
				f/B
				
			f
		})
						
		f <- sum(!XX$majority.response %in% varlevels)
		f <- if(pretty)
		{
			f0 <- c(f/B, f/nrow(XX)*100)
			f0 <- format(f0,digits=3,drop0trailing=T)
			paste(f0[1], " (", paste(f0[2], "%", sep=""), ")", sep="")
		} 
			else
			f/B
		
		m <- data.frame(as.list(c(varlevels.f, f)))
		
		pr <- unique(XX[, x$predictors, drop=F])
		if(nrow(pr)==1)
			 m <- cbind(pr, m)
				
		names(m) <- c(x$predictors, varlevels, "diverse")
		
		
		m
	}))

	X
}

print.familybias <- function(x, ...)
{
# chunk FamilyBias class printing:2
cat(paste("Family bias analysis over", nrow(x$units), "taxa\n"))
# chunk FamilyBias class printing:3
if(!is.null(x$predictors))
{
	cat("--Predictors:\n")
	
	print(unique(x$units[, x$predictors, drop=F]),
		  row.names=F)
}

# chunk FamilyBias class printing:4
print.freq <- function(X)
{
	X <- make.freq(X, x, pretty=T)

	print(format(X),row.names=F)
}

# chunk FamilyBias class printing:5
cat("\n--Large families  estimation results:\n")
print.freq(list(x$large.families.estimate))
# chunk FamilyBias class printing:6
if(!is.null(x$extrapolations))
{
	cat("\n--Prior probabilities:\n")
	print(format(x$prior,digits=3,drop0trailing=T),row.names=F)
		
	spread <- do.call(rbind, 
			     lapply(x$extrapolations, 
						function(X) table(X$majority.response)/nrow(X)))
	cat("\n--Estimation spread\n")
	print(summary(spread))
		
	cat(paste("\n--Probabilistic estimation results (", 
			  length(x$extrapolations), 
			  " simulation",
			  if (length(x$extrapolations)==1) "" else "s",
			  "):\n", sep=""))
	print.freq(x$extrapolations)
}

}

as.data.frame.familybias <- function(x, ...) #res=c("large", "extrapolations", "prior"), pretty=F, ...)
{
	args <- list(...)

	res <- list(...)$res
	if (is.null(res)) res <- 'large'
		
	pretty <- list(...)$pretty
	if (is.null(pretty)) pretty <- F	
	
	
	switch(res,
		large = make.freq(list(x$large.families.estimate), x, pretty=pretty),
		extrapolations = make.freq(x$extrapolations, x, pretty=pretty),
		prior = if(pretty) format(x$prior,digits=3,drop0trailing=T) else x$prior,
		NULL
		)
}


# mean.familybias <- function(x, ...)
# {
# 	varfields <- grep("^number\\.", names(x$units), value=T)
# 	varlevels <- sub("^number\\.( )*", "", varfields)	
# 	
# 	# get the optional order argument
# 	order <- list(...)$order
# 	if (is.null(order)) order <- c('diverse', as.character(varlevels)) 
# 	if (length(setdiff(varlevels, order))!=0)
# 	 stop('Order must include all response variable levels!')
# 	if (!('diverse' %in% order)) order <- c('diverse', order)
# 	
# 	# now we want to build predictor dataframe
# 	hdr <- expand.grid(rev(c(lapply(x$predictors,  function(pred) unique(x$units[, pred])))))
# 	hdr <- hdr[, rev(1:ncol(hdr)), drop=F]
# 	names(hdr) <- x$predictors
# 	hdr0 <- hdr
# 	for(i in 1:ncol(hdr))
# 	{
# 		dd <- as.character(hdr[, i])
# 		r = rle(dd)
# 		pos.first <- sapply(1:length(r$lengths), function(l) if(l==1) 1 else 1+ sum(r$lengths[1:(l-1)]))
# 		dd[-pos.first] <- ""
# 		hdr[, i] <- factor(dd)
# 	}
# 		
# 	dummymat <- as.data.frame(matrix(0, ncol=length(order), nrow=nrow(hdr)))
# 	names(dummymat) <- order
# 	
# 	# collect the data from the extrapolations
# 	tab <- cbind(hdr, dummymat)
# 	for(ex in x$extrapolations)
# 		for(col in order)
# 			for(row in 1:nrow(hdr0))
# 			{
# 				# get the relevant subtable
# 				x0 <- ex
# 				for(pred in names(hdr0))
# 					x0 <- x0[as.character(x0[, pred]) %in% hdr0[row, pred], , drop=F]
# 				tab[row, col] <- tab[row, col] + sum(as.character(x0$majority.response) %in% col)/length(x$extrapolations)	
# 			}
# 	
# 	tab
# }

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
	


	for(ex in x$extrapolations) 
	{
		ex = do.call("paste", c(ex[, names(hdr), drop=F], sep = "\r"))
		for(i in seq_along(hdr0))
			Freq[i] = Freq[i] + sum(ex %in% hdr0[i])
	}
	
	hdr$Freq <- Freq/length(x$extrapolations)
	
	hdr
}
#
#
# ta <- c("stock", "mbranch", "sbranch", "ssbranch", "lsbranch", "language")
# load('data/pro.gender.g.rda')
#
# pro.gender.g$xx <- rnorm(nrow(pro.gender.g))
#
# head(pro.gender.g)
#
# x <- familybias(pro.gender.g, ta, r.name='xx', B = 20, extrapolate=T)
# x
#
# stop()