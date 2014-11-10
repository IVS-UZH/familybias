\name{familybias}
\alias{familybias}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Family bias estimation - familybias()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function estimates biases in the diachronic development of features in language families, with optional extrapolation to small families and isolates.
}
\usage{
familybias(df, family.names, r.name, p.names = NULL, extrapolate = TRUE, 
	B = 1000, small.family.size = 4, diverse.r = "diverse", 
	bias.test = multiple.binom.test, p.threshold = 0.1,
	verbose=F, lapplyfunc=lapply, bias.override = list()) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
An R data frame which contains the language data. 
Each row in the table corresponds to a single language variety.
}
 \item{family.names}{
	A character vector containing the column names of the data frame \code{df}
	which correspond to the taxonomic levels, in order from highest to lowest. The last
			  item of this vector is the column which stores the names of individual language varieties (data points).
}
  \item{r.name}{
	A character string containing the column name which store the values of the response variable. 
This is the variable whose bias we want to characterize using \code{bias.test}.
}
  \item{p.names}{
	A character vector containing the names of the columns which store the predictor variables (conditions) in the data frame \code{df}. An empty value NULL (default) corresponds to a design without predictors. The bias of the
 response variable will be tested conditioned by the predictor values: when one or more predictor
is present, the data is first split into sets such that the predictors are constant within these 
sets. Predictors must be categorical.
}
  \item{extrapolate}{
	A logical value which controls whether or not to perform probabilistic estimation of
	biases in small families and isolates.
}
  \item{B}{
	A number which specifies how many simulations are desired in the probabilistic estimation.
	 Per default, 1000 simulations will be performed.
}
  \item{small.family.size}{
	A number. If a family contains less than this number of individual languages, it is considered too small for estimating biases by the \code{bias.test} and is instead subject to probabilistic estimation (if enabled). Defaults to 4, meaning that 
	language units with 4 or less members are considered too `small' for the \code{bias.test}}
	
  \item{diverse.r}{
The dummy response value by which families with no bias are to be reported (see below). Default is `diverse'. 
}
  \item{bias.test}{a function to be used as a direct bias test. Its single argument is a frequency distribution
of the response variable (a named vector of counts), given predictor level. The function needs to return the value of the variable which represents the bias or NA if there is no bias. The default function is an exact binomial test applied to each pair of levels in the response variable, with a Holm–Bonferroni correction applied in case of multiple pairs. The default \code{p.threshold} is set to 0.1. }

  \item{p.threshold}{The threshold of statistical significance, passed along to the bias.test function}

  \item{verbose}{TRUE if progress messages are desired.}

  \item{lapplyfunc}{A function which will be used for vectorized computation. The arguments and return results of the function must mimic \code{lapply} in \R base.
  This parameter can be used to speed up computation by using a parallelized version of lapply (like the one provided by the \pkg{multicore} package, e.g. \code{lapplyfunc = mclapply}).}

  \item{bias.override}{A list which allows the user to specify custom pre-computed bias estimates for individual families (e.g. using Bayesian methods with continuous-time Markov models as in \code{BayesTraits}). The names in the list should match the family names. Each element of the list should specify two items: \code{bias.val} - the biased category (i.e. the value reported as \code{majority.response} below, or NA if no bias is present) and \code{mj.prop} - the proportion of the biased category in the total counts (\code{as.numeric}). Optionally, the list can also specify \code{mj.val} — for the override of the majority value.}

}
\details{
The function performs family bias estimation as described in Bickel (2011, 2013). Specifically, the language data contained in the input data frame is first split into the largest possibly taxonomic groups (as defined by the \code{family.names} parameter) such that the predictor values are constant within each group. If no predictors are specified, these groups simply correspond to language families (at the maximal taxonomic level given in \code{family.names}).
 % If predictors are specified, these groups are the largest available family branches with uniform predictor values across the board.

For each sufficiently large family (as defined by \code{small.family.size}), a \code{bias.test} is performed. For smaller families,  \code{familybias} can perform an optional probabilistic estimation (controlled by the \code{extrapolate} and \code{B} parameters). Here, the bias estimation results obtained from the larger units are used as an average expectation of the presence of the bias (in any direction) within the smaller units. The probability of bias is estimated via Laplace's rule of succession. For instance, if 70\% of the large families are estimated to have a bias, the extrapolation algorithm will also expect about 70\% of small families to have a bias. Following this, for each small family, the extrapolation algorithm assigns the label `biased' with probability .7 and the label `diverse' with probability .3. 

The members of small families that are estimated to be biased can be representative of the bias or not (i.e. they can be the odd guy(s) out in the family). The probability of this is estimated from the strength of the bias in large families: e.g. if among biased large families, biases tend to be very strong (e.g. on average covering over 90\% of members), we estimate a .9 probability that the members of small families are representative of the larger unknown family from which they derive. We then take the majority value of the small family as the representative value; and a random value in the case of ties. In the small families where members are estimated to represent deviating exceptions, we declare them as survivors of a family that had a bias in an alternative direction (randomly chosen but weighted by the probability of directions given by the general bias estimate and the frequency distribution within the family).

The extrapolation is repeated \code{B} times. The spread within the extrapolation results of individual extrapolation runs can be used to quantify the error of the extrapolation. 

To obtain averaged extrapolation results, use the function \code{mean} on the result of \code{familybias}. The function \code{mean} returns a contingency table in long format (as a \code{dataframe}), with the column \code{Freq} containing the average frequencies. (Compact contingency tables can be obtained by \code{xtabs(Freq~., x)}, where x is the result of applying the function \code{mean}.) 

The package also defines a \code{print} method for the class `familybias', showing a summary of the bias estimation results and of the extrapolation parameters.


}
\value{
The return value is a list of class `familybias'. This list has the following components:

\describe{
\item{units}{ a data frame where each row corresponds to a single genealogical taxon conditioned by the values of the predictor variables (if any are present).}
\item{predictors}{ a character vector with the names of the predictor variables. This component is omitted if no predictors are specified. }
\item{large.families.estimate}{ a data frame with the results of the direct bias estimation. Each row corresponds to a `large' genealogical taxon. }
\item{units.estimate}{same as large.families.estimate, with added rows for small units. These units are not estimated, their distribution is given as `small'.}
\item{extrapolations}{ a list of length \code{B}, containing data frames with the simulated results of probabilistic bias estimation of small families, together with the direct estimation results. This component is omitted if probabilistic bias estimation was not performed (extrapolate parameter set to FALSE)}
\item{prior}{ a list which describes the prior information used for probabilistic bias estimation (i.e. the estimated probabilities of being biased and of being representative of the bias (see `Details'). This component is omitted if probabilistic bias estimation was not performed (extrapolate parameter set to FALSE)}
}

The data frames contained in the components \code{units}, \code{large.families.estimate}, and \code{extrapolations} have the following columns:

\describe{
\item{family.name}{ name of a taxon. The actual taxon is not necessarily a family, it is the largest available taxonomic unit not split by predictor levels. }
		
\item{taxonomic.level}{ the taxonomic level of the taxon. }
		
\item{...}{ values of the predictors (if specified). }
		
\item{number.*}{ the frequency of a particular type (value of the response variable) within the unit. }

\item{number}{ the total number of data points within the taxon. }	
}

In addition, the data frames within the components \code{large.families.estimate} and \code{extrapolations} contain the following columns, related to bias estimation:

\describe{
	\item{majority.response}{ the predominant value of the response variable within the unit. For units which are estimated to be not biased, the value is controlled by the \code{diverse.r} parameter.}
					
	\item{majority.prop}{ proportion of languages within the unit which have the predominant value of the response variable or NA if the unit is estimated to have no bias. }

	\item{distribution}{ the result of unit bias estimation: either `biased' or `diverse'}
	
}

}
\references{
Bickel, B., 2011. Statistical modeling of language universals. \emph{Linguistic Typology} 15, 401 – 414.

Bickel, B., 2013. Distributional biases in language families. In Bickel, B., L. A. Grenoble, D. A. Peterson, & A. Timberlake (eds.) \emph{Language typology and historical contingency}, 415-444, Amsterdam: Benjamins.

}

\author{
Zakharko, Taras & Bickel, Balthasar 

}
\note{
\describe{
\item{Taxonomic levels}{
The genealogical information about the individual languages is stored as a strictly hierarchical sequence of taxa names (e.g. language family name, branch name, language name). As some languages may lack classification at particular taxonomic levels (in particular, language isolates), the language data is allowed to contain `gaps' in the taxonomic description. These gaps appear in the form of NA values. The algorithm creates ad-hoc taxa for such entries, by copying the name of a non-empty lower level taxon. For instance, the language isolate Basque is typically classified on the language level only, and so a data frame containing Basque will contain NA values for entries corresponding to language family and every intermediate branch. In such a case, the algorithm treats the Basque language as a member of a virtual Basque family, Basque branch etc.
}
\item{Pseudo-groups}{
When predictor variables are specified, the splitting algorithm  splits the language data into the largest possible genealogical units such that the distribution of the predictor variables within these units is constant. If no such unit exists above the individual language level (i.e. even the lowest branch has non-uniform predictor values), the algorithm will divide the unit into pseudo-groups, conditioned by the predictor values. See Bickel (2013, section 6) for discussion and justification.
}
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
# Example 1: relative clause position and word order
taxa.names = c("stock", "mbranch", "sbranch", "ssbranch", "lsbranch", "language")

data(vprel.g)

vprel.bias <- familybias(df = vprel.g, family.names = taxa.names, 
                        r.name = 'DRYREL0', p.names = 'DRYSOV4', 
                        extrapolate=T, small.family.size=4, B=50,
                        verbose=T) 	
	
vprel.bias

mean(vprel.bias)

xtabs(Freq~., vprel.bias)

# Example 2: hotbeds of pronominal gender
data(pro.gender.g)
pro.gender.bias <- familybias(df = pro.gender.g, family.names = taxa.names,
                              r.name = 'SIEGEN2', p.names = 'hot', 
                              extrapolate=T, small.family.size=4, 
                              B=50, verbose=T) 	
	
pro.gender.bias

mean(pro.gender.bias)

xtabs(Freq~., pro.gender.bias)

}
}
