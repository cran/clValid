
## First put all methods/class definitions

setClassUnion("numeric or NULL", c("numeric", "NULL"))
setClassUnion("array or NULL", c("array", "NULL"))
setClassUnion("character or array or list", c("character","array","list"))
setClass("clValid",representation(clusterObjs="list",measures="array",measNames="character",clMethods="character",
                                  nClust="numeric",validation="character",
                                  metric="character",method="character",neighbSize="numeric",
                                  annotation="character or array or list",
                                  GOcategory="character",
                                  goTermFreq="numeric",
                                  call="call"))

## cluster methods accessor
setGeneric("clusterMethods", function(object, ...) standardGeneric("clusterMethods"))
setMethod("clusterMethods",signature(object="clValid"),
          function(object) return(object@clMethods))

## number of clusters accessor
setGeneric("nClusters", function(object, ...) standardGeneric("nClusters"))
setMethod("nClusters",signature(object="clValid"),
          function(object) return(object@nClust))
          
## measure names accessor
setGeneric("measNames", function(object, ...) standardGeneric("measNames"))
setMethod("measNames",signature(object="clValid"),
          function(object) return(object@measNames))

## clusters accessor
setGeneric("clusters", function(object, ...) standardGeneric("clusters"))
setMethod("clusters",signature(object="clValid"),
          function(object,method=clusterMethods(object)) {
            method <- match.arg(method,clusterMethods(object))
            return(object@clusterObjs[[method]])})

## measures accessor
setGeneric("measures", function(object, ...) standardGeneric("measures"))
setMethod("measures",signature(object="clValid"),
          function(object,measures=measNames(object)) {
            measures <- match.arg(measures,measNames(object),several.ok=TRUE)
            return(object@measures[measures,,,drop=FALSE])})


setMethod("print","clValid",
          function(x) {
            cat("\nCall:\n")
            print(x@call); cat("\n")
            cat("Clustering Methods:\n",clusterMethods(x),"\n\n")
            cat("Cluster sizes:\n",nClusters(x),"\n\n")
            cat("Validation measures:\n",measNames(x),"\n\n")
          })

setMethod("show","clValid",
          function(object) {
            cat("\nCall:\n")
            print(object@call); cat("\n")
            cat("Clustering Methods:\n",clusterMethods(object),"\n\n")
            cat("Cluster sizes:\n",nClusters(object),"\n\n")
            cat("Validation measures:\n",measNames(object),"\n\n")
          })

setMethod("summary","clValid",
          function(object, digits = max(3,getOption("digits")-3)) {
            cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
            cat("Cluster sizes:\n",nClusters(object),"\n\n")
            cat("Validation Measures:\n")
            print(ftable(round(measures(object),digits),row.vars=c(3,1)))
            cat("\n")
            ## Find best scores
            ## APN, AD, ADM, Connectivity, FOM minimized
            ## BHI, BSI, Dunn, Silhouette maximized
            measNames <- measNames(object)
            best <- numeric(length(measNames))
            bestMeth <- character(length(measNames))
            bestNc <- character(length(measNames))
            names(best) <- names(bestMeth) <- names(bestNc) <- measNames
            minmeas <- c("APN", "AD", "ADM", "FOM", "Connectivity")
            maxmeas <- c("BHI","BSI","Dunn","Silhouette")
                                        ## Measures to minimize
            if (any(a <- minmeas%in%measNames)) {
              best[minmeas[a]] <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,min)
              bestInd <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,function(x) which(x==min(x),arr.ind=TRUE)[1,])
              bestNc[minmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[minmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
                                        ## Measures to maximize
            if (any(a <- maxmeas%in%measNames)) {
              best[maxmeas[a]] <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,max)
              bestInd <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,function(x) which(x==max(x),arr.ind=TRUE)[1,])
              bestNc[maxmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[maxmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            
            cat("Optimal Scores:\n\n") 
            print(data.frame("Score"=round(best,digits),"Method"=bestMeth,"Clusters"=bestNc), right=FALSE)
            cat("\n")
          })




## Add 95% confidence level of random clusters?
## Will depend on whether goal is to minimize or maximize ...
## Can skip for now
setMethod("plot",c("clValid","missing"),
          function(x,y,measures=measNames(x), legend=TRUE, legendLoc="topright", main=NULL,
                   pch=NULL, type="b", ask=prod(par("mfcol")) < length(measures) && dev.interactive(), ...) {
            measures <- match.arg(measures,measNames(x),several.ok=TRUE)
            methods <- clusterMethods(x)
            nclust <- nClusters(x)
            k <- length(methods)
            if (ask) {
              op <- par(ask = TRUE)
              on.exit(par(op))
            }
#            if(is.null(main))
#              main <- paste("Validation Measures for ", deparse(substitute(x, sys.frame(-1))))
            if (is.null(pch)) 
              pch <- c(paste(c(1:9, 0)), letters)[1:k]            
            for(i in 1:length(measures)) {
              if (is.null(main)) {
                main <- switch(measures[i],
                               APN="Stability validation",
                               AD="Stability validation",
                               ADM="Stability validation",
                               FOM="Stability validation",
                               Connectivity="Internal validation",
                               Dunn="Internal validation",
                               Silhouette="Internal validation",
                               BHI="Biological validation",
                               BSI="Biological validation")
              }
              matplot(measures(x)[measures[i],,],type=type,ylab=measures[i],
                      xlab="Number of Clusters",col=1:k,
                      lty=1:k,main=main,xaxt="n", pch=pch, ...)
              axis(1,at=1:length(nclust),labels=nclust)
              if(legend) legend(x=legendLoc,methods,lty=1:k,col=1:k, pch=pch,...)
            }
          })


####################################################################################################
## Some other utility functions
## optimalScores
####################################################################################################



          
## optimalScores method
setGeneric("optimalScores",function(object, ...) standardGeneric("optimalScores"))
setMethod("optimalScores", signature(object="clValid"),
          function(object,measures=measNames(object)) {
            measNames <- match.arg(measures, measNames(object), several.ok=TRUE)
            best <- numeric(length(measNames))
            bestMeth <- character(length(measNames))
            bestNc <- character(length(measNames))
            names(best) <- names(bestMeth) <- names(bestNc) <- measNames
            minmeas <- c("APN", "AD", "ADM", "FOM", "Connectivity")
            maxmeas <- c("BHI","BSI","Dunn","Silhouette")
                                        ## Measures to minimize
            if (any(a <- minmeas%in%measNames)) {
              best[minmeas[a]] <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,min)
              bestInd <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,function(x) which(x==min(x),arr.ind=TRUE)[1,])
              bestNc[minmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[minmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
                                        ## Measures to maximize
            if (any(a <- maxmeas%in%measNames)) {
              best[maxmeas[a]] <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,max)
              bestInd <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,function(x) which(x==max(x),arr.ind=TRUE)[1,])
              bestNc[maxmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[maxmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            return(data.frame("Score"=best,"Method"=bestMeth,"Clusters"=bestNc))
          })






#####################################################################################
## Next put all the functions for cluster validation
#####################################################################################

clValid <- function(obj, nClust, clMethods="hierarchical", validation="stability", maxitems=600, 
                             metric="euclidean", method="average", neighbSize=10, annotation="entrezgene", GOcategory="all", 
                             goTermFreq=0.05,  ...) {

  clMethods <- tolower(clMethods)  
  clMethods <- match.arg(clMethods,c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"),
                              several.ok=TRUE)
  validation <- match.arg(validation,c("stability","internal","biological"),several.ok=TRUE)
  metric <- match.arg(metric,c("euclidean", "correlation", "manhattan")) ## used for hierarchical, diana, fanny, agnes, pam
  method <- match.arg(method,c("ward", "single", "complete", "average")) ## for hclust, agnes
  GOcategory <- match.arg(GOcategory, c("all","BP","CC","MF"))
  
  switch(class(obj),
         matrix = mat <- obj,                                        
         ExpressionSet = mat <- Biobase::exprs(obj),
         data.frame = {
           if(any(!sapply(obj,class)%in%c("numeric","integer")))
             stop("data frame 'obj' contains non-numeric data")
           mat <- as.matrix(obj)
         },
##         dist = Dist <- obj,
         stop("argument 'obj' must be a matrix, data frame, ExpressionSet, or dist object"))

  if (nrow(mat)>maxitems) {
    if (interactive()) {
      cat("\nThe number of items to be clustered is larger than 'maxitems'\n")
      cat("The memory and time required may be excessive, do you wish to continue?\n")
      cat("(y to continue, any other character to exit) ")
      ans <- tolower(substr(readLines(n=1),1,1))
      if (ans!="y") {
        stop("Exiting clValid, number of items exceeds 'maxitems'")
      }
    } else {
      stop("The number of items to be clustered is larger than 'maxitems'\n  Either decrease the number of rows (items) or increase 'maxitems'\n")
    }
  }
    
  
  if ("clara"%in%clMethods & metric=="correlation")
    warning("'clara' currently only works with 'euclidean' or 'manhattan' metrics - metric will be changed to 'euclidean'  ")

  if ("biological"%in%validation & is.character(annotation)) {
    if(!require(Biobase) | !require(GO) | !require(annotate)) {
      stop("packages 'Biobase', 'GO', and 'annotate' required for 2nd type of biological validation \n
these can be downloaded from Bioconductor (www.bioconductor.org)")
    }
  }
  
  if (!is.matrix(mat) | !is.numeric(mat))
    stop("argument 'mat' must be a numeric matrix")
  
  nClust <- floor(nClust)
  if (any(nClust<1))
    stop("argument 'nClust' must be a positive integer vector")

  if(metric=="correlation")
    Dist <- as.dist(1-cor(t(mat)))  else
  Dist <- dist(mat,method=metric)

  clusterObjs <- vector("list",length(clMethods))
  names(clusterObjs) <- clMethods

  measures <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
              if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
              if("biological"%in%validation) c("BHI","BSI"))
  validMeasures <- array(dim=c(length(measures),length(nClust),length(clMethods)))
  dimnames(validMeasures) <- list(measures,nClust,clMethods)

  for (i in 1:length(clMethods)) {
    
    cvalid <- vClusters(mat,clMethods[i],nClust, validation=validation,
                        Dist=Dist, method=method, metric=metric, annotation=annotation,
                        GOcategory=GOcategory, goTermFreq=goTermFreq, neighbSize=neighbSize, ...)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[,,i] <- cvalid$measures
  }

  new("clValid", clusterObjs=clusterObjs, measures=validMeasures, measNames=measures, clMethods=clMethods, nClust=nClust, validation=validation,
      metric=metric,method=method, neighbSize=neighbSize,  GOcategory=GOcategory, goTermFreq=goTermFreq, annotation=annotation, 
      call=match.call())
}

  
  

#########################################################################################
## vClusters() function
####################################################################################################
## Arguments
#########################################################################################
# mat - data matrix
# clMethod - clustering method
## nClust - minimun number of clusters to evaluate
## nclustMax - maximun number of clusters to evaluate
## ... - arguments to pass to clustering functions
#########################################################################################
## Value
#########################################################################################
## List with:
## clusterObj - clustering results
## measures - validation measures
#########################################################################################


vClusters <- function(mat,clMethod,nClust,nclustMax, validation,
                             Dist, method, metric, annotation, GOcategory, goTermFreq, neighbSize, ... ) {
  
  measNames <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
              if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
              if("biological"%in%validation) c("BHI","BSI"))
  measures <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  rownames(measures) <- measNames
  colnames(measures) <- nClust

  switch(clMethod,
         hierarchical = {
           clusterObj <- hclust(Dist,method)
         },
         diana = {
           clusterObj <- diana(Dist, ...)
         },
         kmeans = {
           clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust
           clusterObjInit <- hclust(Dist,method)
         },
         agnes = {
           clusterObj <- agnes(Dist, method=method, ...)
         },
         ## otherwise - sota, fanny, som, model, pam, clara
         { clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust })

  ind <- 1
  for (nc in nClust) {
    switch(clMethod,
           kmeans = {
             initial <- tapply(mat, list(rep(cutree(clusterObjInit,nc),ncol(mat)),col(mat)),mean)
             if(length(dup <- which(duplicated(initial)))>0) {
               for(dupi in dup) 
                 initial[dupi,] <- initial[dupi,] + jitter(initial[dupi,])
             }
             dimnames(initial) <- list(NULL,dimnames(mat)[[2]])
             clusterObj[[ind]] <- kmeans(mat,initial,...)
             cluster <- clusterObj[[ind]]$cluster
           },
           fanny = {
             clusterObj[[ind]] <- fanny(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           model = {
             clusterObj[[ind]] <- Mclust(mat,nc, ...)
             cluster <- clusterObj[[ind]]$classification
           },
           som = {
             clusterObj[[ind]] <- som(mat, grid=somgrid(1,nc), ...)
             cluster <- clusterObj[[ind]]$unit.classif
           },
           pam = {
             clusterObj[[ind]] <- pam(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           clara = {
             clusterObj[[ind]] <- clara(mat, nc, metric=ifelse(metric=="correlation","euclidean",metric), ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           sota = {
             clusterObj[[ind]] <- sota(mat,nc-1)
             cluster <- clusterObj[[ind]]$clust
           },
           ## otherwise - hierarchical, diana, agnes
           {cluster <- cutree(clusterObj,nc)})

    ## internal validation measures
    if ("internal"%in%validation) {
      measures["Dunn",ind] <- dunn(Dist ,cluster)
      measures["Silhouette",ind] <- mean(silhouette(cluster, dmatrix=as.matrix(Dist))[,3])
      measures["Connectivity",ind] <- connectivity(Dist ,cluster, neighbSize=neighbSize)
    }
    
    if("biological"%in%validation) {
      measures["BHI",ind] <- BHI(cluster,annotation=annotation, names=rownames(mat),
                                 category=GOcategory)
    }

    ## stability validation measures
    if ("stability"%in%validation | "biological"%in%validation) {
      for (del in 1:ncol(mat)) {
        matDel <- mat[,-del]               ## matDel <- as.matrix(matDel)
        if(metric=="correlation") DistDel <- as.dist(1-cor(t(matDel))) else DistDel <- dist(matDel,method=metric)
        switch(clMethod,
               hierarchical = clusterObjDel <- hclust(DistDel,method),
               kmeans = clusterObjInitDel <- hclust(DistDel,method),
               diana = clusterObjDel <- diana(DistDel, ...),
               agnes = clusterObjDel <- agnes(DistDel, method=method, ...),
               clara = clusterObjDel <- clara(matDel,nc,metric=ifelse(metric=="correlation","euclidean",metric), ...))


        switch(clMethod,
               kmeans = {
                 initialDel <- tapply(matDel, list(rep(cutree(clusterObjInitDel,nc),ncol(matDel)),col(matDel)),mean)
                 if(length(dup <- which(duplicated(initialDel)))>0) {
                   for(dupi in dup) 
                     initialDel[dupi,] <- initialDel[dupi,] + jitter(initialDel[dupi,])
                 }
                 dimnames(initialDel) <- list(NULL,dimnames(matDel)[[2]])
                 kmdel <- kmeans(matDel,initialDel, ...)
                 clusterDel <- kmdel$cluster
               },
               fanny = {
                 hfdel <- fanny(DistDel, nc, ...)
                 clusterDel <- hfdel$clustering
               },
               model = {
                 clusterDel <- Mclust(matDel,nc, ...)$classification
               },
               som = {
                 hsdel <- som(matDel, grid=somgrid(1,nc), ...)
                 clusterDel <- hsdel$unit.classif
               },
               pam = {
                 clusterDel <- pam(DistDel, nc, cluster.only=TRUE, ...)
               },
               clara = {
                 clusterDel <- clusterObjDel$clustering
               },
               sota = {
                 clusterDel <- sota(matDel,nc-1)$clust
               },
               ## otherwise - hierarchical, diana, agnes
               {clusterDel <- cutree(clusterObjDel,nc)})

        if("stability"%in%validation) {
          stabmeas <- stability(mat, Dist, del, cluster, clusterDel)
          measures["APN",ind] <- measures["APN",ind] + stabmeas["APN"]
          measures["AD",ind]  <- measures["AD",ind]  + stabmeas["AD"]
          measures["ADM",ind] <- measures["ADM",ind] + stabmeas["ADM"]
          measures["FOM",ind] <- measures["FOM",ind] + stabmeas["FOM"]
        }
        if("biological"%in%validation) {
          measures["BSI",ind] <- measures["BSI",ind] + BSI(cluster,clusterDel,annotation=annotation, names=rownames(mat),
                                                           category=GOcategory, goTermFreq=goTermFreq)
        }
        
      } #END OF del LOOP
    } #END of STABILITY measures
    ind <- ind+1  #ind tracks number clusters
  } #END OF NC LOOP

  if ("stability"%in%validation) {
    measures["APN",] <- measures["APN",]/ncol(mat)
    measures["AD",] <-  measures["AD",]/ncol(mat)
    measures["ADM",] <- measures["ADM",]/ncol(mat)
    measures["FOM",] <- measures["FOM",]/ncol(mat)  ## little different from Yeung paper (doesn't do this)
  }
  if ("biological"%in%validation) {
    measures["BSI",] <- measures["BSI",]/ncol(mat)
  }
    
  list(clusterObj=clusterObj, measures=measures)
}



#####################################################################################
## Functions for Validation Measures
## Stability, Internal, and Biological
#####################################################################################

#####################################################################################
## Stability Measures
## APN, AD, ADM, and FOM
## All measures in [0,infty] and should be minimized
#####################################################################################


stability <- function(mat, Dist=NULL, del, cluster, clusterDel, method="euclidean") {

  obsNum <- 1:nrow(mat)
  nc1 <- length(table(cluster))
  nc2 <- length(table(clusterDel))
  stabmeas <- numeric(4)
  names(stabmeas) <- c("APN","AD","ADM","FOM")

  ## measure APN
  ## calculate a ncxnc matrix of proportion of non-overlaps in the two collection of nc clusters
  overlap <- xtabs(~cluster + clusterDel)
  ## measure AD
  ## calculate a ncxnc matrix of average-distance in the two collection of nc clusters
  dij <- matrix(rep(NA,nc1*nc2),nc1,nc2)

  if (is.null(Dist)) matDist <- as.matrix(dist(mat, method=method))
  if (class(Dist)=="dist") matDist <- as.matrix(Dist)
  if (class(Dist)=="matrix") matDist <- Dist
  
  ## measure ADM
  ## calculate a ncxnc matrix of distance-average in the two collection of nc clusters
  dij2 <- matrix(rep(NA,nc1*nc2),nc1,nc2)
  ii <- 1
  for (i in sort(unique(cluster))) {
    jj <- 1
    xbari <- if(is.null(dim(mat[cluster==i,]))) mat[cluster==i,] else
    apply(mat[cluster==i,],2,mean)
    for (j in sort(unique(clusterDel))) {
      ## measure AD
      clusi <- obsNum[cluster==i]
      clusdelj <- obsNum[clusterDel==j]
      cl <- length(clusi)*length(clusdelj)
      if (cl>0) dij[ii,jj] <- mean(matDist[clusi,clusdelj])
##      if (cl>0) dij[ii,jj] <- mean(as.matrix(Dist)[clusi,clusdelj])
      ## measure ADM
      xbarj <- if(is.null(dim(mat[clusterDel==j,]))) mat[clusterDel==j,] else
      apply(mat[clusterDel==j,],2,mean)
      diff <- xbari-xbarj
      dij2[ii,jj] <- ifelse(length(diff)==0, 0, as.numeric(sqrt(diff%*%diff)))
      jj <- jj+1
    }
    ii <- ii+1
  }
  rs <- matrix(rowSums(overlap),nrow=nrow(overlap),ncol=ncol(overlap),byrow=FALSE)
  cs <- matrix(colSums(overlap),nrow=nrow(overlap),ncol=ncol(overlap),byrow=TRUE)
  stabmeas["APN"] <- 1-sum(overlap^2/rs)/sum(overlap)
  stabmeas["AD"] <- sum(overlap*dij)/nrow(mat)
  stabmeas["ADM"] <- sum(overlap*dij2)/nrow(mat)
  xbar <- tapply(mat[,del],clusterDel,mean)
  stabmeas["FOM"] <- sqrt(mean((mat[,del]-xbar[as.character(clusterDel)])^2))/sqrt((nrow(mat)-nc1)/nrow(mat))
  return(stabmeas)
}



####################################################################
## Internal Validation Functions
## Silhouette
## Connectivity
## Dunn index
####################################################################



########################################################################################
## Silhouette
## Value in [-1,1] to be maximized
## NOTE! cluster library also has 'silhouette' function!
## Probably use that instead
#####################################################################################


mysilhouette <- function(distance=NULL, clusters, Data=NULL, method="euclidean"){

  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  dista <- apply(distance,2,function(x) tapply(x, clusters, mean))
  nc <- ncol(dista); nr <- nrow(dista);
  a <- dista[matrix(c(clusters,1:nc),ncol=2,nrow=nc)]
  distb <- matrix(dista[-(clusters+(0:(nc-1))*nr)],ncol=nc,nrow=(nr-1))
  b <- apply(distb,2,min)
  s <- (b-a)/pmax(a,b)
  return(mean(s))
}



##########################################################################################
## Connectivity
## Value in [0,infty] to be *minimized*
#####################################################################################

connectivity <- function(distance=NULL, clusters, Data=NULL, neighbSize=10, method="euclidean"){
  
  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  nearest <- apply(distance,2,function(x) sort(x,ind=TRUE)$ix[2:(neighbSize+1)])
  nr <- nrow(nearest);nc <- ncol(nearest)
  same <- matrix(clusters,nrow=nr,ncol=nc,byrow=TRUE)!=matrix(clusters[nearest],nrow=nr,ncol=nc)
  conn <- sum(same*matrix(1/1:neighbSize,nrow=nr,ncol=nc))
  return(conn)
}



##########################################################################################
## Dunn index
## Value in [0,infty], to be maximized
#####################################################################################

dunn <- function(distance=NULL, clusters, Data=NULL, method="euclidean"){

  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  nc <- max(clusters)
  interClust <- matrix(NA, nc, nc)
  intraClust <- rep(NA, nc)

  for (i in 1:nc) {
    c1 <- which(clusters==i)
    for (j in i:nc) {
      if (j==i) intraClust[i] <- max(distance[c1,c1])
      if (j>i) {
        c2 <- which(clusters==j)
        interClust[i,j] <- min(distance[c1,c2])
      }
    }
  }
  dunn <- min(interClust,na.rm=TRUE)/max(intraClust)
  return(dunn)
}
       
#####################################################################################
## Biological validation functions
## BHI (homogeneity index)
## BSI (stability index)
## Requires annotate and GO packages
#####################################################################################


BHI <- function(statClust,annotation,names=NULL,category="all") {

  
  ## Case 1
  ## Biological clusters provided by user
  if(is.list(annotation)) {
    FF <- length(annotation)
    obs <- unique(unlist(annotation))
    if(is.null(names(statClust)))
      names(statClust) <- names
    stat <- statClust[obs]
    sc <- unique(stat)
    bhi <- 0
    for (k in sc) {
      Ck <- obs[stat==k]
      nk <- length(Ck)
      if (nk >1) {
        count <- 0
        for (i in 1:(nk-1)) {
          for (j in (i+1):nk) {
            f <- 1
            repeat{
              if (is.element(Ck[i],annotation[[f]]) &&
                  is.element(Ck[j],annotation[[f]])) {
                count <- count+1
                break
              }
              f <- f+1
              if (f > FF) break
            }
          }
        }
        bhi <- bhi + 2*count/(nk*(nk-1))
      }
    }
    return(bhi/length(unique(statClust)))
  }

  
  ## Case 2
  ## Name of annotation package provided by user
  ## Gene names assumed to correspond with rownames
  ## Requires Biobase, annotation
  ## Gene names and type of id provided by user
  ##  category <- match.arg(category,c("all","BP","CC","MF"))
  switch(annotation,
         entrezgene = {
           goTerms <- GOENTREZID2GO[names]
         },
         {
           if(!require(annotation,character.only=TRUE)) {
             cat(paste("package",annotation,"not found, attempting download from Bioconductor\n",
                       sep=" "))
             source("http://bioconductor.org/biocLite.R")
             try(biocLite(annotation))
           }
##           if(!require(annotation,character.only=TRUE)) {
##             stop(paste("package",annotation,"not found",sep=" "))
##           }
           goTerms <- getGO(names,annotation)

         })
  
  bhi <- tapply(goTerms,statClust,function(x)  matchGO(x,category))
  return(mean(bhi[bhi!=-9]))
} ## End BHI function


matchGO <- function(gg,category) {
  ## x[1] is x$GOID, x[3] is x$Ontology
  goIDs <- lapply(gg, function(a) sapply(a, function(x) x[1]))
  ont   <- lapply(gg, function(a) sapply(a, function(x) x[3]))
  switch(category,
         BP = {
           goBP <- sapply(ont, function(x) any(x%in%"BP"))
           goIDs <- goIDs[goBP]
         },
         CC = {
           goCC <- sapply(ont, function(x) any(x%in%"CC"))
           goIDs <- goIDs[goCC]
         },
         MF = {
           goMF <- sapply(ont, function(x) any(x%in%"MF"))
           goIDs <- goIDs[goMF]
         },
         all = {
           goAll <- sapply(goIDs, function(a) all(sapply(a,function(x) is.null(x))))
           goIDs <- goIDs[!goAll]
         })
  n <- length(goIDs)
  if (n<2) return(-9)
  sum <- 0
  for (i in 1:(length(goIDs)-1)) {
    for (j in (i+1):length(goIDs)) {
      switch(category,
             all =  sum <- sum + any(goIDs[[i]]%in%goIDs[[j]]),
             BP  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="BP"]%in%goIDs[[j]][ont[[j]]=="BP"]),
             CC  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="CC"]%in%goIDs[[j]][ont[[j]]=="CC"]),
             MF  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="MF"]%in%goIDs[[j]][ont[[j]]=="MF"]))
    }
  }
  return(sum/(n*(n-1)))
} ## End matchGO function


BSI <- function(statClust,statClustDel,annotation,names=NULL,category="all", goTermFreq=0.05) {

  ## Case 1
  ## Biological clusters provided by user

  ## fixme
  ## browser()
  if(is.list(annotation)) {
    nClasses <- sapply(annotation, length)
    FF <- length(annotation)
    if(is.null(names(statClust))) {
      names(statClust) <- names
      names(statClustDel) <- names
    }
    s <- 0
    overlap <- xtabs(~statClust + statClustDel)
    rsums <- rowSums(overlap)
    for (f in 1:FF)
      {
        osum <- 0
        for (gx in annotation[[f]] ) {
          for (gy in annotation[[f]] ) {
            if (gx != gy) {
              i <- statClust[gx]
              j <- statClustDel[gy]
              osum <- osum + overlap[i,j]/rsums[i]
            }
          }
        }
        s <- s + osum/(nClasses[f]*(max(nClasses[f]-1,1)))
      }
    return(s/FF)
  }

  
  ## Case 2
  ## Gene names and type of id provided by user
  ## For option 2 of BSI need to determine how many GO terms to use
  ## Cutoff determined by goTermFreq (cuts all GO terms w/ < 5% freq in Data set)
#  category <- match.arg(category,c("all","BP","CC","MF"))

  tab <- xtabs(~statClust + statClustDel)
  rs <- rowSums(tab)
  n <- length(statClust)
  
  switch(annotation,
         entrezgene = {
           goTerms <- GOENTREZID2GO[names]
         },
         {
         if(!require(annotation,character.only=TRUE)) {
           cat(paste("package",annotation,"not found, attempting download from Bioconductor\n",
                     sep=" "))
           source("http://bioconductor.org/biocLite.R")
           try(biocLite(annotation))
         }
##           if(!require(annotation,character.only=TRUE)) {
##             stop(paste("package",annotation,"not found",sep=" "))
##           }
           goTerms <- getGO(names,annotation)
         })


  ## Things to do
  ## 1. extract all relevant GO terms for each gene (bp,cc, etc)
  ## 2. get freq table of terms - keep most frequent
  ## 3. create indicator matrix for each term indicating which genes have that term

  ## fixme
  ## browser()
  switch(category,
         ## x[1] is x$GOID, x[3] is x$Ontology
         BP = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="BP"])),
         CC = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="CC"])),
         MF = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="MF"])),
         all = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1])))

  goTab <- table(unlist(goIDs))
  keepTerms <- names(goTab)[goTab>floor(n*goTermFreq)]
  termMat <- matrix(0,ncol=length(keepTerms),nrow=n)
  for (i in 1:length(keepTerms)) {
    termMat[,i] <- sapply(goIDs, function(x) keepTerms[i]%in%x)
  }

  bsi <- apply(termMat,2, function(a) {
    out <- outer(statClust[as.logical(a)],statClustDel[as.logical(a)], function(x,y) tab[cbind(x,y)]/rs[x])
    (sum(out) - sum(diag(out)))/(sum(a)*(ifelse(sum(a)>1,sum(a)-1,1)))
  })
  
  return(mean(bsi))
  
}


#####################################################################################
## Now functions for clustering
## SOTA 
#####################################################################################


sota.init <- function(data){
	nodes <- matrix(0, nrow(data)*2, 3+ncol(data))
	if(is.null(colnames(data)))
		colnames(data) <- paste("V", 1:ncol(data))
	colnames(nodes) <- c("ID", "anc", "cell", colnames(data))
	nodes[,"ID"]=1:(nrow(data)*2)	

	nodes[1,] <- c(1, 0, 0, apply(data,2,mean))
	nodes[2,] <- c(2, 1, 1, nodes[1,][-c(1,2,3)])
	nodes[3,] <- c(3, 1, 1, nodes[1,][-c(1,2,3)])
	return(nodes)
}

dist.fn <- function(input, profile, distance){
	if(distance=="correlation")
		return(1-cor(input,profile))
	else
		return(sqrt(sum((input-profile)^2)))
}

cl.ID <- function(clust, old.cl, new.cl){
	for(i in 1:length(clust))
		clust[i] <- new.cl[which(old.cl==clust[i])]
	clust
}

getResource <- function(data, tree, clust, distance, pr){
	dist <- rep(0, length(clust))
	resource <- rep(0, max(clust))
	
	for(i in unique(clust)){
		temp <- data[clust==i,]
		if(is.vector(temp))
			temp <- matrix(temp, nrow=1, ncol=ncol(data))
		if(distance=="correlation")
			resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr], 
				distance=distance))
		else
			resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr], 
				distance=distance))}
	resource
}

getCells <- function(tree, neighb.level, n){
	or.n <- n
	cells <- c(n-1,n)
	for(i in 1:(neighb.level+1)){
		n  <- tree[n, "anc"]
		if(n==1)
			break
	}
	for(j in 2:(or.n-2)){
		z <- j
		if(tree[j,"cell"]!=1)
			next
		while(z > 0){
			z <- tree[z, "anc"]
			if(z==n){
				cells <- c(cells, j)
				break}
		}
	}
	return(tree[cells,])
}


sota <- function(data, maxCycles, maxEpochs=1000, distance="euclidean",
			wcell=.01, pcell=.005, scell=.001, delta=.0001, neighb.level=0,
			maxDiversity = .9, unrest.growth=TRUE, ...){
	tree <- sota.init(data)
	pr <- 4:ncol(tree)
	n <- 3
	genes<- 1:nrow(data)
	clust <- rep(1, length(genes))
	Node.Split <- 1
	
	for(k in 1:maxCycles){                                #loop for the Cycles
		trainNode <- Node.Split
		trainSamp <- genes[clust==trainNode]
		curr.err <- 1e10
		ep <- 1
		while(ep <= maxEpochs){	      			#loop for the Epochs
			last.err <- 0
			left.ctr <- right.ctr <- 0
			left.d <- right.d <- 0
			for(i in trainSamp){
				cells <- tree[c(n-1,n),]
				dist <- rep(0, nrow(cells))
				for(j in 1:2)
					dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)
											
				or <- which.min(dist)
				if(or==1)
					left.ctr <- left.ctr + 1
				else
					right.ctr<- right.ctr + 1
				
				closest <- cells[or,1]				
				sis <- ifelse(closest%%2==0,closest+1,closest-1)
				sis.is.cell <- ifelse(tree[sis,"cell"]== 1, 1, 0)
				
				##   updating the cell and its neighbourhood
				if(sis.is.cell==1){
					parent <- tree[closest, "anc"]
					tree[closest, pr] <- tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
					tree[sis, pr] <- tree[sis, pr]+scell*(data[i,]-tree[sis,pr])
					tree[parent, pr] <- tree[parent, pr]+pcell*(data[i,]-tree[parent, pr])
				}
				else
				{
					tree[closest, pr] <- tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
				}
			}
			cells <- tree[c(n-1,n),]
			for(i in trainSamp){
				for(j in 1:2)
					dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)
				last.err <- last.err+min(dist)}
			last.err <- last.err/length(trainSamp)
			
			if(ifelse(last.err==0, 0, abs((curr.err-last.err)/last.err)) < delta
			   && left.ctr !=0 && right.ctr !=0)
				break
			ep <- ep + 1
			curr.err <- last.err
		}
		clust <- assignGenes(data, trainSamp, clust, tree, n, distance, pr, neighb.level)
		Res.V <- getResource(data, tree, clust, distance, pr)
		if(k==maxCycles)
			break  ## do not split the cell 
		newCells <- splitNode(Res.V, tree, n)
		tree <- newCells$tree
		n <- newCells$n
		Node.Split <- newCells$toSplit
		if(max(Res.V) < maxDiversity & unrest.growth==FALSE)
			break
	}

	tree <- trainLeaves(data, tree, clust, pr, wcell, distance, n, delta)
	Res.V <- getResource(data, tree, clust, distance, pr)
	Res.V <- Res.V[Res.V!=0]
	if(distance=="correlation")
		Res.V <- 1-Res.V

	treel <- tree[tree[,"cell"]==1,]
	old.cl <- treel[,1]
	treel[,1] <- 1:nrow(treel)
	old.clust <- clust
		
	clust <- cl.ID(old.clust, old.cl, 1:nrow(treel))
	
	totals <- table(clust)
	out <- list(data=data, c.tree=tree[1:n,], tree=treel, clust=clust, 
			totals=totals, dist=distance, diversity=Res.V)
	
	class(out) <- "sota"
	return(out)
}

trainLeaves <- function(data, tree, clust, pr, wcell, distance, n, delta){
	nc <- ncol(data)
	for(i in 1:n){
		if(!is.element(i, clust))
			next
		temp <- matrix(data[clust==i,], ncol=nc)
		converged <- FALSE
		init.err <- getCellResource(temp, tree[i,pr], distance)
		while(!converged){
			for(j in 1:nrow(temp))
				tree[i, pr] <- tree[i, pr]+wcell*(temp[j,]-tree[i, pr])	
		
			last.err <- getCellResource(temp, tree[i,pr], distance)
			converged <- ifelse(abs((last.err-init.err)/last.err) < delta, TRUE, FALSE)
			init.err <- last.err
		}
	}
	return(tree)
}


assignGenes <- function(data, Sample, clust, tree, n, distance, pr, neighb.level){
	if(neighb.level==0)
		cells <- tree[c(n-1,n),]	
	else	
		cells <- getCells(tree, neighb.level, n)

	for(i in Sample){
		dist <- rep(0, nrow(cells))
		for(j in 1:nrow(cells))
			dist[j] <- dist.fn(data[i,], cells[j,pr], distance)
		or <- which.min(dist)
		closest <- cells[or,1]	
		clust[i] <- closest
	}
	clust
}

splitNode <- function(Res.V, tree, n){
	maxheter <- which.max(Res.V)
	cl.to.split <- tree[maxheter,1]
	tree[n<-n+1,-1] <- tree[cl.to.split,-1]
	tree[n, "anc"] <- cl.to.split
	tree[n<-n+1,-1] <- tree[cl.to.split,-1]
	tree[n, "anc"] <- cl.to.split
	tree[cl.to.split, "cell"] <- 0
	return(list(tree=tree, n=n, toSplit=cl.to.split))
}


	
getCellResource <- function(temp, profile, distance){
	if(distance=="correlation")
		resource <- mean(apply(temp, 1, dist.fn, profile,
			distance=distance))
	else(distance=="euclidean")
		resource <- mean(apply(temp, 1, dist.fn, profile,
			distance=distance))
	resource
}

print.sota <- function(x, ...){
  results <- as.matrix(cbind(x$tree[,1], as.numeric(x$totals),
                             x$diversity))
  colnames(results) <- c("ID","Size", "Diversity")
  rownames(results) <- rep("", nrow(results))
  cat("\nClusters:\n")
  print(results, ...)
  cat("\nCentroids:\n")
  
  print(format(data.frame(x$tree[,-c(1:3)])), ...)
  cat("\n")
  cat(c("Distance: ", x$dist, "\n"))
  invisible(x)
}


plot.sota <- function(x, cl=0, ...){

  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if(cl!=0)
    par(mfrow=c(1,1)) else
  {	
    pdim <- c(0,0)
    for(i in 1:100){
      j <- i
      if(length(x$totals) > i*j)
        j <- j+1 
      else{
        pdim <- c(i,j)
        break}
      if(length(x$totals) > i*j)
        i <- i+1 
      else{
        pdim <- c(i,j)
        break}
    }
    par(mfrow=pdim)
  }
  
  ylim = c(min(x$data), max(x$data))
  pr <- 4:ncol(x$tree)
  if(cl==0)
    cl.to.print <- 1:max(x$clust) else
    cl.to.print <- cl

  for(i in cl.to.print){
    plot(1:ncol(x$data), x$tree[i, pr], col="red", type="l",
         ylim=ylim, xlab=paste("Cluster ",i), ylab="Expr. Level", ...)
    legend("topleft", legend=paste(x$totals[i], " Genes"), cex=.7, 
           text.col="navy", bty="n")
    cl <- x$data[x$clust==i,]
    if(is.vector(cl))
      cl <- matrix(cl, nrow=1)
    for(j in 1:x$totals[i])
      lines(1:ncol(x$data), cl[j,], col="grey")
    lines(1:ncol(x$data), x$tree[i, pr], col="red", ...)

  }
}
