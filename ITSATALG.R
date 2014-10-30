################################################################################
#                                                                              #
#                         _                 _       _                          #
#                      __| |_   _  ___ _ __| | __ _| |__                       #
#                     / _` | | | |/ _ \ '__| |/ _` | '_ \                      #
#                    | (_| | |_| |  __/ |  | | (_| | |_) |                     #
#                     \__,_|\__, |\___|_|  |_|\__,_|_.__/                      #
#                           |___/                                              #
#                                                                              #
################################################################################

# Author: R.J. Dyer
# Email: rjdyer@vcu.edu
# Web: http://dyerlab.com
# This is the complete R script for the analyses conducted for the manuscript
#   "Is there such as thing as landscape genetics?"  The data input for this is 
#   a data frame "df.filtered.rda" that should be in the same directory as this 
#   script.

rm(list=ls())
require(tm)
require(ape)
require(MASS)
library(pvclust)
require(ggplot2)
require(skmeans)
require(SnowballC)
require(wordcloud)
require(RColorBrewer)


# Load in the raw text data frame (has columns with Intro & Methods)
load("df.filtered.rda")

################################################################################
#                                                                              #
#  Translate the raw data from filtered text to vector space format            #
#                                                                              #
################################################################################

filter_and_stem <- function( raw_text, type, filenames, standardize=TRUE, verbose=FALSE ) {
  
  # to Corpus analyses
  corpus <- Corpus( DataframeSource(data.frame(raw_text)))
  corpus <- tm_map( corpus, content_transformer(removePunctuation))
  corpus <- tm_map( corpus, content_transformer(tolower))
  corpus <- tm_map( corpus, content_transformer(removeNumbers ))
  corpus <- tm_map( corpus, content_transformer(stripWhitespace))
  corpus <- tm_map( corpus, content_transformer(
    function(x) { removeWords(x, stopwords(kind="en")) } ))
  corpus <- tm_map( corpus, content_transformer(
    function(x) { removeWords(x, stopwords(kind="SMART")) } ))
  
  # Report Unique Terms
  d <- paste( raw_text, collapse=" ")
  d <- strsplit(d,split = " ")[[1]]
  cat("Unique terms:",length(unique( d )),"\n")
  
  # Stem each document (SnowballC currently mpi-like error)
  for( i in 1:length(corpus) ){
    x <- corpus[[i]]
    y <- stemDocument( x )
    corpus[[i]] <- y
  }
  
  # Turn it into a matrix and keep terms > 5%
  tdm <- TermDocumentMatrix( corpus )  
  m <- as.matrix( tdm )
  
  cat("Stemmed Terms:",nrow(m),"\n")
  n <- m
  n[ n > 0 ] <- 1
  r <- rowSums( n ) 
  m <- m[ r>5 ,]
  cat("Terms 95:", nrow(m),"\n")
  
  # Normalize the vectors to unit length
  r <- matrix( colSums( m ), nrow=nrow(m), ncol=ncol(m))
  m <- m / r
  
  cat("Retained Terms:",nrow(m),"\n")
  data <- data.frame( Group=type, stringsAsFactors = FALSE )
  data <- cbind( data, t(m) )
  
  # Remove terms that are homogeneous within groups
  if( standardize ){
    for( i in seq(ncol(data),2,by=-1)){
      s <- as.numeric( by( data[,i], data$Group, var))
      if( any( s == 0 ) ) {
        data[,i] <- NULL
        if( verbose ) 
          cat("Dropping: ", names(data)[i], paste(s,sep=", "),"\n")
      }
    }
  }
  cat("Variable Terms:",ncol(data)-1,"\n")
  
  rownames(data) <- filenames
  
  return( data )  
}




# Translate raw data
cat("\nFilter & Stem Introduction\n")
df.intro <- filter_and_stem( df$Introduction, df$Type, df$File )

cat("\nFilter & Stem Methods\n")
df.methods <- filter_and_stem( df$Methods, df$Type, df$File )


################################################################################
#                                                                              #
#   Perform discriminant analysis and classify the LG papers into groups       #
#                                                                              #
################################################################################

discriminate_and_classify <- function( df ){
  training <- droplevels( df[ df$Group != "Landscape Genetics", ]  )
  suppressWarnings( fit.lda <- lda( Group ~ ., data=training ) )
  
  # classify
  test <- predict( fit.lda, newdata = df)
  cat("LG Classification\n")
  t <- table(test$class,df$Group)
  print(t)
  
  # test for equality
  cat("Testing for Equality of Allocation for LG papers")
  fit.b <- binom.test( x=t[,2],p=0.5)
  print(fit.b)
  
  # return as data frame
  ret <- data.frame( Type=df$Group, 
                     Predicted=test$class,
                     LDA1=test$x )
  return( ret )
}

# Fit Introduction
cat("\nDiscriminate & Classify Introduction\n")
df.lda.intro <- discriminate_and_classify( df.intro )
df.lda.intro$Section <- "Introduction"

# Fit the Methods
cat("\nDiscriminate & Classify Methods\n")
# 553, 587 & 877 are constant within groups (+1 in df.methods )
df.methods[,554] <- df.methods[,588] <- df.methods[,878] <- NULL
df.lda.methods <- discriminate_and_classify( df.methods )
df.lda.methods$Section <- "Methods"

# Plot both 
d <- rbind( df.lda.intro, df.lda.methods)
p <- ggplot( d, aes(x=LD1, fill=Type, color=Type)) + geom_density(alpha=0.75)
p <- p + facet_grid(Section~.) + theme_gray( base_size=16 ) 
p <- p + xlab("Linear Discriminant Axis 1") + ylab("Manuscript Density")
print(p)






################################################################################
#                                                                              #
#   Differential Stemmed Term Usage LG vs. LE vs. PG for Intro & Methods       #
#                                                                              #
################################################################################

make_relative_frequency_vector <- function( x ){
  x[ x>0] <- 1
  ret <- colSums(x) / nrow(x)
  return(ret)
}

## Introduction
intro.lg <- make_relative_frequency_vector( df.intro[ df.intro$Group=="Landscape Genetics",seq(2,ncol(df.intro))] )
intro.le <- make_relative_frequency_vector( df.intro[ df.intro$Group=="Landscape Ecology", seq(2,ncol(df.intro))] )
intro.pg <- make_relative_frequency_vector( df.intro[ df.intro$Group=="Population Genetics",seq(2,ncol(df.intro))])

lg.le <- intro.lg - intro.le
print("LE-LG: Intro LE Bias")
print(sort(lg.le)[1:10])
print("LE-LG: Intro LG Bias")
print(rev(sort(lg.le))[1:10])

lg.pg <- intro.lg - intro.pg
print("PG-LG: Intro PG Bias")
print(sort(lg.pg)[1:10])
print("PG-LG: Intro LG Bias")
print(rev(sort(lg.pg))[1:10])




## Methods
methods.lg <- make_relative_frequency_vector( df.methods[ df.methods$Group=="Landscape Genetics",seq(2,ncol(df.methods))] )
methods.le <- make_relative_frequency_vector( df.methods[ df.methods$Group=="Landscape Ecology", seq(2,ncol(df.methods))] )
methods.pg <- make_relative_frequency_vector( df.methods[ df.methods$Group=="Population Genetics",seq(2,ncol(df.methods))])

lg.le <- methods.lg - methods.le
print("LE-LG: Methods LE Bias")
print(sort(lg.le)[1:10])
print("LE-LG: Methods LG Bias")
print(rev(sort(lg.le))[1:10])

lg.pg <- methods.lg - methods.pg
print("PG-LG: Methods PG Bias")
print(sort(lg.pg)[1:10])
print("PG-LG: Methods LG Bias")
print(rev(sort(lg.pg))[1:10])




################################################################################
#                                                                              #
#   Conduct Clustering                                                         #
#                                                                              #
################################################################################
x <- df.methods[df.methods$Group=="Landscape Genetics",-1]
D <- dist(x,method = "euclidean")
hc <- hclust(D, method="ward.D" )
labelColors <- c("#6699aa","#ABA380","#ff7f00","#e31a1c")
groups <- cutree(hc, 3)
hcp <- as.phylo(hc)
plot( hcp, type="unrooted",  show.tip.label=FALSE, cex=0.5  )
tiplabels(pch=21, bg=labelColors[groups], cex=1.5)
X.mat <- x

################################################################################
#                                                                              #
#   Develop a wordcloud for two main groups                                    #
#                                                                              #
################################################################################

t <- table(groups)
cat("Manuscripts in each group: \n")
print(t)

pal1 <- colorRampPalette(c("#447799","#6699aa","#C9CFD2"))
pal2 <- colorRampPalette(c("#4D4915","#ABA380","#C9CFD2"))
pal3 <- colorRampPalette(c("#ffeda0","#feb24c","#f03b20"))
numWordsSamp <- 10

grp1 <- x[ groups==1, ]
grp1[ grp1 > 0 ] <- 1
v1 <- sort( colSums(grp1)/nrow(grp1),decreasing=TRUE)
d1 <- data.frame( word=names(v1), freq=v1)
summary(d1)
d1 <- d1[1:numWordsSamp,]
d1$rank <- seq(numWordsSamp,1,by=-1)
wordcloud( d1$word, seq(numWordsSamp,1,by=-1),colors=rev(pal1(10)), scale=c(5,2), min.freq=1)



grp2 <- x[ groups==2, ]
grp2[ grp2 > 0 ] <- 1
v2 <- sort( colSums(grp2)/nrow(grp2),decreasing=TRUE)
d2 <- data.frame( word=names(v2), freq=v2)
summary(d2)
d2 <- d2[ 1:numWordsSamp,]
d2$rank <- seq(numWordsSamp,1,by=-1)
wordcloud( d2$word,d2$rank,colors=rev(pal2(10)), scale=c(5,2), min.freq=1)


grp3 <- x[ groups==3,]
v3 <- sort( colSums(grp3)/nrow(grp3),decreasing=TRUE)
d3 <- data.frame( words=names(v3),freq=v3)
d3 <- d3[ d3$freq > 0,]
d3 <- d3[1:numWordsSamp,]
d3$rank <- seq(numWordsSamp,1,by=-1)
wordcloud( d3$word, d3$rank, colors=pal3(10), scale=c(5,2), min.freq=1)

grp.diff <- abs(colSums(grp1)/nrow(grp1) - colSums(grp2)/nrow(grp2) )
v.diff <- sort( grp.diff, decreasing=TRUE)
d.diff <- data.frame( words=names(v.diff), freq=v.diff )
summary(d.diff)


d.diff.s <- d.diff[ d.diff$freq>=0.10,]
d.diff.s$freq <- d.diff.s$freq / sum(d.diff.s$freq)
wordcloud( d.diff.s$word,
           d.diff.s$freq
           )





















################################################################################
#                                                                              #
#   Develop a wordcloud for two main groups                                    #
#                                                                              #
################################################################################

if(0) {
  
  X <- t(as.matrix(x))
  suppressMessages(bootstrap.result <- pvclust( X, method.hclust = "ward", method.dist="euclidean", nboot=10000))
  save(bootstrap.result,file="clustering.bootstrap.result.rda")  
  print("Finished Bootstrap") 
} 

# plot the supported clusters
load("clustering.bootstrap.result.rda")
source("PlotBranchByTrait.R")
edges <- hcp$edge
val <- rep(0,nrow(edges))
colors <- rep("black",nrow(edges))  
groups <- rep(0,length(hcp$tip.label))
clusters <- cutree(hc, 3)
r <- pvpick(bootstrap.result,alpha = 0.95)
for( cluster in r$cluster ){
  e <- which.edge( hcp,group = cluster)
  val[e] <- 1
  colors[e] <- labelColors[4]
  
  idx <- match( cluster, hcp$tip.label )
  groups[idx] <- max(groups)+1
}

plotBranchbyTrait(hcp,val,mode="edge",type="unrooted", legend=F,colors=colors, show.tip.label=FALSE)  
tiplabels(pch=21, bg=labelColors[clusters], cex=1.5)




################################################################################
#                                                                              #
#   Clustering Within Observed Groups                                          #
#                                                                              #
################################################################################

if( 0 ) {
  # remove constant within group 
  const <- c(4,11,59,120,301,310,336,366,376,463,
             485,533,555,660,663,667,797,810,870,
             879,917,918,927,1019,1099,1307,1411,
             1476,1479,1481,1529,1618,1641,1698,
             1737,1789,1907,1908,1911,1982,2050,2057,2104)
  const <- -1*rev(const)
  xp <- x[const]
  (fit.lgmethods.lda <- lda( x = xp, grouping = groups))
  lgmethods.pred <- predict( fit.lgmethods.lda,dimen=2)
  vals <- lgmethods.pred$x
  t <- table(lgmethods.pred$class, groups)
  cat("\nDiscriminant Separation of LG Methods Sections:")
  print(t)
  
  print("Missclassification Rate:")
  missclassification <- 1-sum(diag(t)) / sum(t)
  print(missclassification)
  
  
  lda.arrows <- function(x, myscale = 1, tex = 0.75, choices = c(1,2), largest=2, ...){
    ## adds `biplot` arrows to an lda using the discriminant function values
    heads <- coef(x)
    heads2 <- abs(heads)
    keepers <- c( rownames(heads2[ rev(order(heads2[,1]))[1:largest],]),
                  rownames(heads2[ rev(order(heads2[,2]))[1:largest],]))
    heads <- heads[ rownames(heads) %in% keepers, ]
    arrows(x0 = 0, y0 = 0, 
           x1 = myscale * heads[,choices[1]], 
           y1 = myscale * heads[,choices[2]], ...)
    text(myscale * heads[,choices], labels = row.names(heads), cex = tex)
  }
  plot(fit.lgmethods.lda,asp=1)
  lda.arrows(fit.lgmethods.lda,col=2,myscale=0.025, largest=3)
  
  
  ################################################################################
  #                                                                              #
  #   Stemmed Frequency Analyses  By Cluster Group                               #
  #                                                                              #
  ################################################################################
  
  group1 <- x[ groups==1,]
  group2 <- x[ groups==2,]
  group3 <- x[ groups==3,]
  
  freq_vec <- function( data ){
    ndata <- as.matrix(data)
    ndata[ ndata != 0 ] <- 1
    ret <- colSums(ndata) / nrow(ndata)
    names(ret) <- names(data)
    return(ret)
  } 
  
  
  freq.1 <- freq_vec( group1 )
  freq.2 <- freq_vec( group2 )
  freq.3 <- freq_vec( group3 )
  
  
  
}




