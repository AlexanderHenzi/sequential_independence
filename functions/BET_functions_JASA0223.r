subsets=function(n,r,v=1:n){
  if(r<=0) NULL else
  if(r>=n) matrix(v[1:n],nrow=1) else
  rbind(cbind(v[1],subsets(n-1,r-1,v[-1])),subsets(n-1,r,v[-1]))
}

frac.2 <- function(a,d){				# converting a fraction 0<a<1 to binary expansions up to depth d 
	a0 <- a
	b <- rep(NA,d)
	for (i in 1:d){
		temp <- 2^i*a0
		if (temp>1){
			b[i] <- 1
			a0 <- a0-1/2^i
		}else{
			b[i] <- 0
		}
	}
	b
}

frac2 <- function(a,d){					# vectorizing frac.2
	if (d==1){
	t(t(sapply(a,function(x){frac.2(x,d)})))
	}else{
	t(sapply(a,function(x){frac.2(x,d)}))
	}
}

bex.centers <- function(depth=3){			# depth >=1
	cbind(rep((1:2^depth)/2^depth,2^depth),rep( (1:2^depth)/2^depth, rep(2^depth,2^depth)   ))-1/2^(depth+1)
}

plot.bet <- function(depth=3,be.ind1,be.ind2){
	xyc <- bex.centers(depth)
	BEx <- frac2(xyc[,1],depth)
	BEy <- frac2(xyc[,2],depth)
	
	RDx <- 2*BEx-1
	RDy <- 2*BEy-1
	
	be.ind1.num <- as.numeric(unlist(strsplit(be.ind1,":")))
	x.prod <- apply(RDx[,be.ind1.num,drop=F],1,prod)
	
	be.ind2.num <- as.numeric(unlist(strsplit(be.ind2,":")))
	y.prod <- apply(RDy[,be.ind2.num,drop=F],1,prod)
	
	col.ind <- x.prod*y.prod
	
	for (i.col in 1: nrow(xyc)){
		if (col.ind[i.col]<0){
			xycc <- xyc[i.col,]
			xp <- c(xycc[1]-1/2^(depth+1),xycc[1]+1/2^(depth+1),xycc[1]+1/2^(depth+1),xycc[1]-1/2^(depth+1))
			yp <- c(xycc[2]-1/2^(depth+1),xycc[2]-1/2^(depth+1),xycc[2]+1/2^(depth+1),xycc[2]+1/2^(depth+1))
			polygon(xp,yp,border=NA,col=rgb(0,0,1,1/4))
		}
	}
	
}

BET <- function(x,y,depth=3,plot=F,cex=0.5){
    n <- length(x)
	
	BEx <- frac2(x,depth)
	BEy <- frac2(y,depth)
	
	RDx <- 2*BEx-1
	RDy <- 2*BEy-1
	
	BEx.complete <- NULL
	BEx.complete.colnames <- NULL

	BEy.complete <- NULL
	BEy.complete.colnames <- NULL

	for (id in 1:depth){
        sets <- subsets(depth,id)
		if (is.null(dim(sets))){
			nint <- 1
			sets <- matrix(sets,nrow=1)
		}else{
			nint <- nrow(sets)
		}
        for(iint in 1:nint){
            temp1 <- apply(RDx[,sets[iint,],drop=F],1,prod)
			BEx.complete.colnames <- c(BEx.complete.colnames,paste(sets[iint,],collapse=":"))
            BEx.complete <- cbind(BEx.complete,temp1)
			
            temp2 <- apply(RDy[,sets[iint,],drop=F],1,prod)
			BEy.complete.colnames <- c(BEy.complete.colnames,paste(sets[iint,],collapse=":"))
            BEy.complete <- cbind(BEy.complete,temp2)			
        }
	}	
	
	colnames(BEx.complete) <- BEx.complete.colnames
	colnames(BEy.complete) <- BEy.complete.colnames
	
	count.interaction <- t(apply(BEx.complete,2,function(x){colSums(x*BEy.complete)})) # transpose so that x interactions are rows

	if (nrow(count.interaction)==1){
		rownames(count.interaction) <- "1"
		colnames(count.interaction) <- "1"
	}
	
	max.abs.count.interaction <- max(abs(count.interaction))
	max.ind <- which(abs(count.interaction)==max.abs.count.interaction,arr.ind=T)		# max.ind is a matrix
	max.index <- c(colnames(count.interaction)[max.ind[1,1]],colnames(count.interaction)[max.ind[1,2]])		
	p.value <- min((2^depth-1)^2*2*(1-pbinom((max.abs.count.interaction+n)/2-1,n,prob=1/2)),1) # Bonferroni p-value
	strongest.asymmetry <- (count.interaction[max.ind])[1]
	
	table22 <- matrix(c(max.abs.count.interaction/4+n/4, -max.abs.count.interaction/4+n/4,-max.abs.count.interaction/4+n/4,max.abs.count.interaction/4+n/4  ),2,2)
	FE22 <- fisher.test(table22,conf.int=FALSE)$p.value- dhyper(table22[1,1],n/2,n/2,n/2)/2
	FE.pvalue <- min((2^depth-1)^2*FE22,1)
	
	chisq.stat <- sum(count.interaction^2)/n
	chisq.pvalue <- 1-pchisq(chisq.stat, df=(2^depth-1)^2)
	
	cell.BEx <- apply(BEx,1,function(a){sum(a*rev(cumprod(c(1,rep(2,depth-1)))))})+1
	cell.BEy <- apply(BEy,1,function(a){sum(a*rev(cumprod(c(1,rep(2,depth-1)))))})+1
	cell.store <- matrix(0,2^depth,2^depth)
	
	for(ic in 1:n){
		cell.store[cell.BEx[ic],cell.BEy[ic]] <- cell.store[cell.BEx[ic],cell.BEy[ic]]+1
	}
	
	LRT <- sum(cell.store[cell.store>0]*log(cell.store[cell.store>0]))- sum((rowSums(cell.store)[rowSums(cell.store)>0])*log((rowSums(cell.store)[rowSums(cell.store)>0]))) - sum((colSums(cell.store)[colSums(cell.store)>0])*log((colSums(cell.store)[colSums(cell.store)>0])))+ n*log(n)
	
	LRT.pvalue <- 1-pchisq(2*LRT,df=(2^depth-1)^2)

	if (plot==T){
			par(mgp = c(1.8, 0.5, 0),mar=c(3,3,3,1))
			plot(c(0,1), c(0,1), xlab=expression(U[x]),ylab=expression(U[y]),type = "n")
			points(x,y,mgp = c(1.8, 0.5, 0),xlim=c(0,1),ylim=c(0,1),cex=cex,col=2, pch=16)
			plot.bet(depth,max.index[1],max.index[2])
	}
	return(list(p.value=p.value,FE.pvalue=FE.pvalue,count.interaction=count.interaction,strongest.asymmetry=strongest.asymmetry,max.index=max.index,chisq.stat=chisq.stat,chisq.pvalue=chisq.pvalue,cell.store=cell.store,LRT=LRT,LRT.pvalue=LRT.pvalue))
}

BETd <- function(u,v,d1=2,d2=2,plot.bid=F,cex=0.5){
    n <- length(u)
	
	BEx <- frac2(u,d1)
	BEy <- frac2(v,d2)
	
	RDx <- 2*BEx-1
	RDy <- 2*BEy-1

	count.interaction <- matrix(NA,2^(d1-1),2^(d2-1))		# x interactions are rows
	
	if (d1==1){
		rownames(count.interaction) <- index.names1  <- "1"	
	}else{
		index.matrix1 <- as.matrix(expand.grid( rep( list( 0:1), d1-1)))
		index.names1 <- apply(index.matrix1,1,function(a){paste(which(c(rev(a),1)>0),collapse=":")})
		rownames(count.interaction) <- index.names1
	}
	
	if (d2==1){
		colnames(count.interaction) <- index.names2  <- "1"		
	}else{
		index.matrix2 <- as.matrix(expand.grid( rep( list( 0:1), d2-1)))
		index.names2 <- apply(index.matrix2,1,function(a){paste(which(c(rev(a),1)>0),collapse=":")})
		colnames(count.interaction) <- index.names2
	}	
	

	for (int1 in 1:length(index.names1)){
		for (int2 in 1:length(index.names2)){
			count.interaction[int1,int2] <- sum(apply(cbind(RDx[,as.numeric(unlist(strsplit(index.names1[int1], "[:]")))],RDy[,as.numeric(unlist(strsplit(index.names2[int2], "[:]")))]),1,prod))
		}
	}
	
	
	max.abs.count.interaction <- max(abs(count.interaction))
	max.ind <- which(abs(count.interaction)==max.abs.count.interaction,arr.ind=T)		# max.ind is a matrix
	max.index <- c(rownames(count.interaction)[max.ind[1,1]],colnames(count.interaction)[max.ind[1,2]])		
	p.value <- min(2^(d1+d2-2)*2*(1-pbinom((max.abs.count.interaction+n)/2-1,n,prob=1/2)),1) # Bonferroni p-value
	strongest.asymmetry <- (count.interaction[max.ind])[1]
	
	table22 <- matrix(c(max.abs.count.interaction/4+n/4, -max.abs.count.interaction/4+n/4,-max.abs.count.interaction/4+n/4,max.abs.count.interaction/4+n/4  ),2,2)
	FE22 <- fisher.test(table22,conf.int=FALSE)$p.value- dhyper(table22[1,1],n/2,n/2,n/2)/2
	FE.pvalue <- min(2^(d1+d2-2)*FE22,1)
	
	if (plot.bid){
			par(mgp = c(1.8, 0.5, 0),mar=c(3,3,3,1))
			plot(c(0,1), c(0,1), xlab=expression(U[x]),ylab=expression(U[y]),type = "n")
			points(u,v,mgp = c(1.8, 0.5, 0),xlim=c(0,1),ylim=c(0,1),cex=cex,col=2, pch=16)
			plot.bet(max(d1,d2),max.index[1],max.index[2])
	}
	return(list(p.value=p.value,FE.pvalue=FE.pvalue,count.interaction=count.interaction,strongest.asymmetry=strongest.asymmetry,max.index=max.index))
}

BETs <- function(u,v,d.max=4){
    n <- length(u)
	temp <- BET(u,v,1)
	bet.adj.pvalues <- rep(NA,d.max)
	bet.max.ind <- matrix(NA,d.max,2)	
	bet.adj.pvalues[1] <- temp$FE.pvalue
	bet.max.ind[1,] <- temp$max.index
		temp.colnames <- colnames(temp$count.interaction)
		temp.rownames <- rownames(temp$count.interaction)

	if (d.max==1){
	return(list(bet.s.pvalue=bet.adj.pvalues[d.max],bet.s.index=bet.max.ind[d.max,]))
	}else{	
	for (id in 2:d.max){
		tempa <- BET(u,v,id)
		tempa.colnames <- colnames(tempa$count.interaction)
		tempa.rownames <- rownames(tempa$count.interaction)
		
		tempa.count.interaction <- tempa$count.interaction
		
		tempa.count.interaction[temp.rownames,temp.colnames] <- 0
		
	max.abs.count.interaction <- max(abs(tempa.count.interaction))
	max.ind <- which(abs(tempa.count.interaction)==max.abs.count.interaction,arr.ind=T)		# max.ind is a matrix
	max.index <- c(colnames(tempa.count.interaction)[max.ind[1,1]],colnames(tempa.count.interaction)[max.ind[1,2]])		
	p.value <- min(((2^id-1)^2-(2^(id-1)-1)^2)*2*(1-pbinom((max.abs.count.interaction+n)/2-1,n,prob=1/2)),1) # Bonferroni p-value
	strongest.asymmetry <- (tempa.count.interaction[max.ind])[1]
	
	table22 <- matrix(c(max.abs.count.interaction/4+n/4, -max.abs.count.interaction/4+n/4,-max.abs.count.interaction/4+n/4,max.abs.count.interaction/4+n/4  ),2,2)
	FE22 <- fisher.test(table22,conf.int=FALSE)$p.value- dhyper(table22[1,1],n/2,n/2,n/2)/2
	FE.pvalue <- min(FE22*((2^id-1)^2-(2^(id-1)-1)^2),1)
		
		bet.adj.pvalues[id] <- FE.pvalue
		bet.max.ind[id,] <- max.index
		
		temp.colnames <- tempa.colnames
		temp.rownames <- tempa.rownames
		
	}
	bet.s.pvalue <- min(min(bet.adj.pvalues)*d.max,1)
	bet.s.index <- bet.max.ind[which(bet.adj.pvalues==min(bet.adj.pvalues)),]
	return(list(bet.s.pvalue=bet.s.pvalue,bet.s.index=bet.s.index))
	}
}

