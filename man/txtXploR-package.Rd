\name{txtXploR-package}
\alias{txtXploR-package}
\alias{txtXploR}
\docType{package}
\title{
Text Mining, tweet-orient exploration.
}
\description{
Facilitates the use of popular tm package' functions. All functions but one are executed from term-document matrix. Involves very little data manipulation to produce charts and other outputs. Get some tweets, turn to corpus, use the clean function then get a term-document matrix, which is required to use other functions. 
}
\details{
\tabular{ll}{
Package: \tab txtXploR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-25\cr
}
Functions:

\itemize{
    \item acc_fetch
    \item clean_corpus
    \item hier_clust
    \item pam_k
    \item plot_topic
    \item term_network
    \item tag_cloud
    
}

Please see their individual documentation for info.
}
\author{
John Coene

Maintainer: <jcoenep@gmail.com>
}
\references{
\href{http://cran.r-project.org/doc/contrib/Zhao_R_and_data_mining.pdf}{RDataMining} by Yanchang Zhao.
\href{http://socialfunction.wordpress.com/}{SocialFunction()}
}