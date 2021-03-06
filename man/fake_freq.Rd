% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_dataset.R
\name{fake_freq}
\alias{fake_freq}
\title{fake haplotype frequencies}
\usage{
fake_freq(nhaps = 100, nsnps = 10)
}
\arguments{
\item{nhaps}{number of haplotypes}

\item{nsnps}{number of snps}
}
\value{
data.frame of 1,2, nhaps x nsamples + frequency column
}
\description{
create a fake haplotype frequency dataset
}
\details{
no attempt is made at biological realism, this is purely for testing code
}
\examples{
freq=fake_freq(nhaps=100,nsnps=5)
dim(freq)
head(freq)
}
\author{
Chris Wallace
}
