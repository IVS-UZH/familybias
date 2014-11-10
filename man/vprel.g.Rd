\name{vprel.g}
\alias{vprel.g}
\docType{data}
\title{
Relative clause position and word order
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Language sample on relative clause position and word order from Dryer (2005a,2005b) linked to the genealogical classification of Nichols & Bickel (2009). 

}
\usage{data(vprel.g)}
\format{
  A data frame containing following information on 506 language varieties:
  \describe{
    \item{\code{stock}}{Language family name}
    \item{\code{mbranch}}{Major branch name}
    \item{\code{sbranch}}{Subbranch name}
    \item{\code{ssbranch}}{Sub-subbranch name}
    \item{\code{lsbranch}}{Lowest subbranch name}
    \item{\code{language}}{Language (dialect/variety) name}

    \item{\code{DRYSOV4}}{word order in respect to the position of the verb, \code{final} vs \code{non_final}}

    \item{\code{DRYREL0}}{order of the relative clause in the relation tot he noun, \code{NRel} vs \code{RelN}}
 }
}
\references{

Nichols, Johanna & Balthasar Bickel, 2009. The AUTOTYP  genealogy and geography database: 2009
release. Electronic database, \url{http://www.uni-leipzig.de/~autotyp}.


Dryer, Matthew S. 2005a. Order of Relative Clause and Noun.
In Haspelmath, Martin, Matthew S. Dryer, David Gil, & Bernard Comrie (eds.) \emph{The world atlas of language structures,}
182–185. Oxford: Oxford University Press.


Dryer, Matthew S. 2005b. Order of Subject, Object and Verb.
In Haspelmath, Martin, Matthew S. Dryer, David Gil, & Bernard Comrie (eds.) \emph{The world atlas of language structures,}
182–185. Oxford: Oxford University Press.
}
\examples{
data(vprel.g)

head(vprel.g)
}
\keyword{datasets}
