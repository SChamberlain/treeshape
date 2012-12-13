#' Written by Liam J. Revell 2011
#' From the phytools package
#' citation: Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things) Methods Ecol. Evol. 3, 217-223 doi:10.1111/j.2041-210X.2011.00169.x
#' @keywords internal
#' @export 
nodeHeights <- function (tree) 
{
	if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
																												 "order"))) 
		t <- reorder(tree)
	else t <- tree
	root <- length(t$tip) + 1
	X <- matrix(NA, nrow(t$edge), 2)
	for (i in 1:nrow(t$edge)) {
		if (t$edge[i, 1] == root) {
			X[i, 1] <- 0
			X[i, 2] <- t$edge.length[i]
		}
		else {
			X[i, 1] <- X[match(t$edge[i, 1], t$edge[, 2]), 2]
			X[i, 2] <- X[i, 1] + t$edge.length[i]
		}
	}
	if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
																												 "order"))) 
		o <- apply(matrix(tree$edge[, 2]), 1, function(x, y) which(x == 
																															 	y), y = t$edge[, 2])
	else o <- 1:nrow(t$edge)
	return(X[o, ])
}

#' Written by Liam J. Revell 2011
#' From the phytools package
#' citation: Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things) Methods Ecol. Evol. 3, 217-223 doi:10.1111/j.2041-210X.2011.00169.x
#' @keywords internal
#' @export 
pbtree <- function (b = 1, n = NULL, t = NULL, scale = NULL, nsim = 1, 
										ape = TRUE) 
{
	if (nsim > 1) {
		trees <- replicate(nsim, pbtree(b, n, t, scale), simplify = FALSE)
		class(trees) <- "multiPhylo"
		return(trees)
	}
	else {
		if (!is.null(t)) 
			stop("time stop not yet implemented")
		else {
			node <- n + 1
			edge <- matrix(c(node, NA, node, NA), 2, 2, byrow = T)
			edge.length <- c(0, 0)
			node <- node + 1
			tip <- 0
			while (nrow(edge) < (2 * n - 2)) {
				o <- is.na(edge[, 2])
				p <- which(o)
				q <- sample(p)[1]
				edge[q, 2] <- node
				edge <- rbind(edge, matrix(c(node, NA, node, 
																		 NA), 2, 2, byrow = T))
				node <- node + 1
				l <- rexp(n = 1, sum(o) * b)
				edge.length[p] <- edge.length[p] + l
				edge.length <- c(edge.length, rep(0, 2))
			}
			o <- is.na(edge[, 2])
			p <- which(o)
			l <- rexp(n = 1, sum(o) * b)
			edge.length[p] <- edge.length[p] + l
			edge[is.na(edge[, 2]), 2] <- tip + 1:sum(is.na(edge[, 
																													2]))
			tip.label <- paste("t", 1:n, sep = "")
			tree <- list(edge = edge, edge.length = edge.length, 
									 tip.label = tip.label, Nnode = n - 1)
			class(tree) <- "phylo"
			if (!is.null(scale)) {
				h <- max(nodeHeights(tree))
				tree$edge.length <- scale * tree$edge.length/h
			}
			if (ape) 
				tree <- read.tree(text = write.tree(tree))
			return(tree)
		}
	}
}

#' Written by Liam J. Revell 2011
#' From the phytools package
#' citation: Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things) Methods Ecol. Evol. 3, 217-223 doi:10.1111/j.2041-210X.2011.00169.x
#' @keywords internal
#' @export 
fastBM<-function(tree,a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),internal=FALSE,nsim=1){
	
	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(bounds[2]<bounds[1]){
		warning("bounds[2] must be > bounds[1]. Simulating without bounds.")
		bounds<-c(-Inf,Inf)
	}
	if(bounds[1]==-Inf&&bounds[2]==Inf) no.bounds=TRUE
	else no.bounds=FALSE
	if(a<bounds[1]||a>bounds[2]){
		warning("a must be bounds[1]<a<bounds[2]. Setting a to midpoint of bounds.")
		a<-bounds[1]+(bounds[2]-bounds[1])/2
	}
	if(sig2<0){
		warning("sig2 must be > 0.  Setting sig2 to 1.0.")
		sig2=1.0
	}
	
	# function for reflection off bounds
	reflect<-function(yy,bounds){
		while(yy<bounds[1]||yy>bounds[2]){
			if(yy<bounds[1]) yy<-2*bounds[1]-yy
			if(yy>bounds[2]) yy<-2*bounds[2]-yy
		}
		return(yy)
	}
	
	# how many species?
	n<-length(tree$tip)
	
	# first simulate changes along each branch
	x<-matrix(data=rnorm(n=length(tree$edge.length)*nsim,mean=rep(mu*tree$edge.length,nsim),sd=rep(sqrt(sig2*tree$edge.length),nsim)),length(tree$edge.length),nsim)
	
	# now add them up
	y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),nsim))
	for(i in 1:nrow(x)){
		if(tree$edge[i,1]==(n+1))
			y[i,1,]<-a
		else
			y[i,1,]<-y[match(tree$edge[i,1],tree$edge[,2]),2,]
		
		y[i,2,]<-y[i,1,]+x[i,]
		if(!no.bounds) y[i,2,]<-apply(as.matrix(y[i,2,]),1,function(yy) reflect(yy,bounds))
	}
	
	rm(x); x<-matrix(data=rbind(y[1,1,],as.matrix(y[,2,])),length(tree$edge.length)+1,nsim)
	rownames(x)<-c(n+1,tree$edge[,2])
	x<-as.matrix(x[as.character(1:(n+tree$Nnode)),])
	rownames(x)[1:n]<-tree$tip.label
	
	# return simulated data
	if(internal==TRUE)
		return(x[1:nrow(x),]) # include internal nodes
	else
		return(x[1:length(tree$tip.label),]) # tip nodes only
	
}