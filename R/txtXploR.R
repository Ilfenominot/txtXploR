#' clean_corpus
#' 
#' @description
#' Cleans corpus of documents.
#' 
#' @param corpus Required. Corpus of documents to clean.
#' @param stem Optional, boolean. If TRUE stems documents after cleaning. Defaults to FALSE.
#' @param rem_words Optional, additional words to remove from documents.
#' 
#' @details
#' Uses tm_map functions from \href{http://cran.r-project.org/web/packages/tm/tm.pdf}{tm package}.
#' 
#' * Strips white space
#' * To lower case
#' * Removes @@tags
#' * Removes punctuation
#' * Removes numbers
#' * Removes stopwords
#' * To plain text
#' 
#' @return cleaned corpus
#' 
#' @seealso
#' \href{http://cran.r-project.org/web/packages/tm/tm.pdf}{tm package}
#' 
#' @examples
#' 
#' \dontrun{
#'  
#'  library(tm)
#'  
#'  data(crude)
#'  
#'  crude <- clean_corpus(crude, stem = TRUE, rem_words = c("saudi"))
#'  
#'  }
#'  
clean_corpus <- function(corpus, stem = FALSE, rem_words="the") {
  if (class(corpus)[1] != "VCorpus") {
    stop ("Need object of type corpus")
  } else if (missing(corpus)) {
    stop("Missing corpus")
  } else if (!is.character(rem_words)) {
    stop("rem_words of wrong class !character")
  } 
  library(tm)
  corpus <- tm_map(corpus, stripWhitespace)
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, function(x) gsub("@[^ ]*", "", x))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, removeNumbers)
  corpus <- tm_map(corpus, function(x) gsub('http[[:alnum:]]*', '', x))
  corpus <- tm_map(corpus, removeWords, c(stopwords(kind = 'en'), 'amp','via','and','for',
                                          'from', rem_words))
  corpus <- tm_map(corpus, PlainTextDocument) 
  if (stem == TRUE) {
    library("tm")
    dict <- corpus
    corpus <- tm_map(corpus, stemDocument)
    myCorpus <- tm_map(corpus, stemCompletion, dictionary=dict)
  }
  return(corpus)
}

#' hier_clust
#' 
#' @description
#' Hierarchical cluster analysis presented in dendrogram.
#' 
#' @param tdm Required. Term-document matrix.
#' @param clusters Optional. Number of clusters based on distance between terms, defaults to 5.
#' @param sparsity Optional. Maximum allowed sparsity. Between 0 and 1, defaults to 0.98.
#' @param method Optional. Clustering method, defaults to ward.D (see details) - other: single linkage, complete linkage, average linkage, median and centroid
#' 
#' @return Clusters, plot fit.
#' 
#' @details
#' Transforms tdm, removes sparse terms and clusters words. Default method is ward.D which clusters terms according to the varianced.
#' 
#' @seealso
#' \code{\link{pam_k}} to cluster around medoids
#' \code{\link[tm]{TermDocumentMatrix}} TermDocumentMatrix from corpus
#' \code{\link[stats]{hclust}} Clustering methods
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(tm)
#' 
#' data(crude)
#' 
#' tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#' 
#' fit <- hier_clust(tdm, sparsity = 0.97)
#' 
#' plot(fit)
#' 
#' }
hier_clust <- function (tdm, clusters = 5, sparsity = 0.98, method = "ward.D") {
  if(missing(tdm)) {
    stop("Missing TermDocumentMatrix")
  } else if (class(tdm)[1] != "TermDocumentMatrix") {
    stop("object must be TermDocumentMatrix")
  } else if (clusters <= 1) {
    stop ("cluster must be between 2 and 58")
  }
  tdm2 <- removeSparseTerms (tdm, sparse=sparsity)
  m2 <- as.matrix(tdm2)
  dist_matrix <- dist(scale(m2))
  fit <- hclust(dist_matrix, method='ward.D')
  plot(fit)
  rect.hclust (fit, k=clusters)
  return(fit)
}

#' pam_k
#' 
#' @description
#' Like k-means, attempts to break dataset into groups and minimize distance between points.
#' 
#' @param tdm Required. object of class TermDocumentMatrix.
#' @param sparsity Optional. Maximum allowed sparsity. Between 0 and 1, defaults to 0.98.
#' 
#' @return
#' 2D clustering plot - clusters 
#' Silhouette plot - higher the si the better clustered
#' Prints clusters in console
#' 
#' @details
#' Partitions around \href{http://en.wikipedia.org/wiki/K-medoids}{k-medoids} and produces 2D clustering plot as well as silhouette plot to test the fit of clusters.
#' 
#' @seealso
#' \code{\link{hier_clust}} for hierarchical clustering
#' \href{http://cran.r-project.org/web/packages/fpc/fpc.pdf}{fpc} package (dependency)
#' 
#' @examples
#' 
#' \dontrun{
#' library(tm)
#' 
#' data(crude)
#' 
#' tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#' 
#' pam_k(tdm, sparsity = 0.99)
#' }

pam_k <- function(tdm, sparsity = 0.98) {
  if (missing(tdm)) {
    stop("no term-document matrix")
  } else if (class(tdm)[1] != "TermDocumentMatrix") {
    stop("object must be TermDocumentMatrix")
  } 
  tdm <- removeSparseTerms (tdm, sparse=sparsity)
  m <- t(as.matrix(tdm))
  library(fpc)
  pam <- pamk(m)
  c <- pam$nc
  k <- pam$nc
  pam <- pam$pamobject
  for (i in 1:k) {
    cat(paste('cluster ', i,': '))
    cat(colnames(pam$medoids)[which(pam$medoids[i,]==1)],'\n')
  }
  plot <- plot(pam, color=F, labels=4, lines=0, cex=.8, col.clus=1,
               col.p=pam$clustering)
  return(plot)
}

#' plot_topic
#' 
#' @description
#' 
#' Models topics using Latent Dirichlet Allocation and plots against time using qplot.
#' 
#' @return Plots of topics over time (ggplot2)
#' 
#' @param tdm Required. term-document matri; must be of class "TermDocumentMatrix"
#' @param date Required. Vector of dates to plot topics against. POSIXlt class
#' @param topics Optional. Number of topics to find (numeric). Defaults to 5.
#' @param terms Optional. Number of terms to display for each topics; defaults to 4.
#' 
#' @details
#' Uses \href{http://cran.r-project.org/web/packages/topicmodels/topicmodels.pdf}{topicmodels} package as well as \href{http://ggplot2.org/}{ggplot2} package to plot topics over time.
#' Topics are found with the LDA model - \href{http://en.wikipedia.org/wiki/Latent_Dirichlet_allocation}{Latent Dirichlet Allocation}.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(tm)
#' 
#' data(crude)
#' 
#' tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#' x <- as.POSIXlt(c(Sys.Date() - 1, Sys.Date(), Sys.Date() + 1 ))
#' dates <- sample(x, 20, replace =TRUE)
#' 
#' plot_topic(tdm, dates)
#' 
#' }
plot_topic <- function(tdm, date, topics = 5, terms = 4) {
  if (missing(tdm)) {
    stop("Missing term-document matrix")
  } else if (missing(date)) {
    stop("no date vector provided")
  } else if (class(tdm)[1] != "TermDocumentMatrix") {
    stop ("tdm not of TermDocumentMatrix class")
  } else if (class(topics) != "numeric") {
    stop("topics not numeric")
  } else {
    libs <- c("ggplot2", "topicmodels")
    lapply(libs, library, character.only=TRUE)
    dtm <- as.DocumentTermMatrix (tdm)
    rowTotals <- apply(dtm , 1, sum)
    dtm.new   <- dtm[rowTotals> 0, ]
    lda <- LDA(dtm.new, k = topics)
    term <- terms (lda, terms)
    term <- apply(term, MARGIN = 2, paste, collapse = ", ")
    topics <- topics(lda, 1)
    topics <- data.frame(date=date[1:length(topics)], topics)
  }
  qplot(date, ..count.., data = topics, geom = "density" , fill = term[topics], position ="stack")
}

#' tag_cloud
#' 
#' @description creates word-cloud of most frequent terms in term-document matrix fed to the function
#' @param tdm Required. Term-document matrix
#' @param min_freq Optional, defaults to 1. Minimum frequency of term for it to appear in Wordcloud
#' @param scale Optional. A vector of length 2 indicating the range of the size of the words.
#' @param order Optional, defaults to FALSE. Whether the order of the terms in wordcloud are random
#' @param colour Optional. Brewer colour palette, defaults to brewer.pal(8, 'dark2')
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(tm)
#' library(wordcloud)
#' 
#' data(crude)
#' 
#' tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#' 
#' tag_cloud(tdm, min_freq = 15, scale = c(3, .6), order = FALSE, colour = brewer.pal(9,"BuGn"))
#' 
#' }
#' 
tag_cloud <- function(tdm, min_freq = 1, scale = c(4, .5), order = FALSE,
                      colour = brewer.pal(8, 'Dark2')) {
  if (class(tdm)[1] != "TermDocumentMatrix") {
    stop ("object not of TermDocumentMatrix Class")
  } else if (missing(tdm)) {
    stop("no Term-Document Matrix")
  } else {
    m <- as.matrix(tdm)
    v <- sort(rowSums(m), decreasing=T)
    names <- names(v)
    d <- data.frame(words=names, freq=v)
    library(wordcloud) 
  }
  wordcloud(d$words, d$freq, scale=scale, min.freq=min_freq,
            random.order=order, rot.per=.25, color=colour)
}

#' term_bar
#' 
#' @description term frequency using tdm - barplot 
#' 
#' @param tdm Required. Object of class TermDocumentMatrix.
#' @param term_freq Optional, defaults to 1. Minimum frequency of terms to plot
#' @param colour Optional, defaults to heat. Color to be used for bars.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(tm)
#' 
#' data(crude)
#' 
#' tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#' 
#' term_bar(tdm, min.freq=15)
#' 
#' }
term_bar <- function(tdm, term_freq = 1, colour = "heat") {
  if (class(tdm)[1] != "TermDocumentMatrix") {
    stop("object not of TermDocumentMatrix class")
  } else if (missing(tdm)) {
    stop("missing TermDocumetnMatrix")
  } else if (colour == "heat"){
    term_frequency <- rowSums(as.matrix(tdm))
    term_frequency <- subset(term_frequency, term_frequency>=term_freq)
    term_frequency <- sort(term_frequency, decreasing=TRUE)
    colour <- heat.colors(length(term_frequency))
  } else {
    term_frequency <- rowSums(as.matrix(tdm))
    term_frequency <- subset(term_frequency, term_frequency>=term_freq)
    term_frequency <- sort(term_frequency, decreasing=TRUE)
  }
  barplot (term_frequency, las=2, col = colour)
}

#' term_network
#' 
#' @description
#' Creates netowrk of terms object (igraph) from term-document matrix
#' 
#' @param tdm Required. Term-document matrix.
#' @param sparsity Optional. Maximum allowed sparsity. Between 0 and 1, defaults to 0.97.
#' @param weighted Optional boolean. Whether the graph should be weighted. Defaults to TRUE.
#' @param graphml Optional. If TRUE will write graph as graphml file in wd, defaults to FALSE.
#' 
#' @seealso
#' \href{http://cran.r-project.org/web/packages/igraph/igraph.pdf}{igraph}
#' 
#' @examples
#' 
#'  \dontrun{
#'  
#'  library(tm)
#'  
#'  data(crude)
#'  
#'  tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))
#'  
#'  net <- term_network(tdm, sparsity = 0.99, weighted=FALSE)
#'  
#'  tkplot(net)
#'  
#'  }
#' 
term_network <- function(tdm, sparsity = 0.98, weighted = TRUE, graphml = FALSE) {
  if(class(tdm)[1] != "TermDocumentMatrix") {
    stop("object of wrong class - need TermDocumentMatrix")
  } else if (missing(tdm)){
    stop("Missing TermDocumentMatrix")
  } else {
    library("igraph")
    tdm2 <- removeSparseTerms (tdm, sparse=sparsity)
    m2 <- as.matrix(tdm2)
    m2[m2>=1] <- 1
    term_matrix <- m2 %*% t(m2)
    if (weighted == FALSE) {
      weighted = NULL
      g <- graph.adjacency(term_matrix, weighted=weighted, mode="undirected")
      g <- simplify(g)
      V(g)$label <- V(g)$name
      V(g)$degree <- degree(g)
    } else {
      g <- graph.adjacency(term_matrix, weighted=weighted, mode="undirected")
      g <- simplify(g)
      V(g)$label <- V(g)$name 
      egam <- (log(E(g)$weight)+.4) / max(log(E(g)$weight)+.4)
      E(g)$color <- rgb(.5, .5, 0, egam)
      E(g)$width <- egam
      V(g)$degree <- degree(g)
    }
    if (graphml == TRUE) {
      write.graph(g, "term_network.graphml", format='graphml')
      print("saving graph as object and writing term_network.graphml")
    } else {
      print("saving graph as object")
    }
  } 
  return(g)
}