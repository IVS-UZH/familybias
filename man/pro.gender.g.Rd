\name{pro.gender.g}
\alias{pro.gender.g}
\docType{data}
\title{
Hotbeds of pronominal gender
}
\usage{data(pro.gender.g)}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Language sample on gender distinction in independent pronouns (Siewierska 2005)
and residence within a typological hotbed (Nichols 1992). 
The taxonomic classification is Nichols & Bickel's (2009).
}
\format{
  A data frame containing following information about 383 language varieties:
  \describe{
    \item{\code{stock}}{Language family name}
    \item{\code{mbranch}}{Major branch name}
    \item{\code{sbranch}}{Subbranch name}
    \item{\code{ssbranch}}{Sub-subbranch name}
    \item{\code{lsbranch}}{Lowest subbranch name}
    \item{\code{language}}{Language (dialect/variety) name}

    \item{\code{SIEGEN2}}{gender distinction in independent personal pronouns: \code{has_gender} \emph{vs.} \code{no_gender}}

    \item{\code{hot}}{whether the language resides within a `hotbed' (according to Nichols 1992): \code{outside hotbeds} \emph{vs.} \code{inside hotbeds}}
 }
}
\references{


Nichols, Johanna & Balthasar Bickel, 2009. The AUTOTYP  genealogy and geography database: 2009
release. Electronic database, \url{http://www.uni-leipzig.de/~autotyp}.

Siewierska, Anna, 2005. Gender distinctions in independent personal pronouns. In Haspelmath,
Martin, Matthew S. Dryer, David Gil, & Bernard Comrie (eds.) \emph{The world atlas of language structures,}
182–185. Oxford: Oxford University Press.

Nichols, Johanna, 1992. \emph{Linguistic diversity in space and time.} Chicago: Thee University of Chicago Press.
}
\examples{
data(vprel.g)

head(vprel.g)
}
\keyword{datasets}