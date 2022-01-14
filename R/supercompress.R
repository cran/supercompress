# auxiliary function
L2=function(D,x,y){
  N=length(y)
  cluster=drop(knnx.index(D,x,k=1))
  ss=function(u) sum((u-mean(u))^2)
  val=sum(sapply(split(y,cluster),ss))
  return(val/N)
}
L2t=function(D,z){
  N=dim(D)[1]
  p=dim(D)[2]
  cluster=drop(knnx.index(D,matrix(z[,1:p],ncol=p),k=1))
  val=sum(as.matrix(stats::aggregate(z,by=list(cluster),sumofsq)[,-1]))
  return(val/N)
}
sumofsq=function(v) sum((v-mean(v))^2)


#' @title 
#' Supervised Data Compression via Clustering.
#'
#' @description
#' \code{supercompress} is the supervised data compression method proposed in Joseph 
#' and Mak (2021). It is a nonparametric compression method that incorporates
#' information of the response.
#' 
#' @details 
#' The \code{supercompress} algorithm finds the \code{n} compressed points by 
#' sequentially splitting the space into \code{n} Voronoi regions with centers
#' being the \code{n} compressed points. The splitting is done to minimize the 
#' total within-cluster sum of squares. The parameter \code{lam} 
#' controls the robustness of the splitting, with value 0 being fully 
#' supervised (objective based on response \code{y} only) and value 1 being fully 
#' unsupervised (objective based on feature \code{x} only), where the latter case 
#' reduces to the kmeans clustering. The Vornoi regions are identified 
#' by the fast nearest neighbor search implemented in the R package \code{FNN}.
#' Only continuous response and features are supported at this time. 
#' Default is to standardize the big data to have zero mean and unit variance
#' before processing. Please see Joseph and Mak (2021) for details. 
#' 
#' @author Chaofan Huang and V. Roshan Joseph
#' 
#' @import FNN
#' @importFrom stats aggregate kmeans sd
#' 
#' @param n number of compressed data points
#' @param x features of the input big data
#' @param y responses of the input big data
#' @param lam robustness parameter takes value between 0 (fully supervised) 
#' and 1 (fully unsupervised)
#' @param standardize should the big data be normalized to 
#' have zero mean unit variance
#' 
#' @return
#' \item{D}{features of compressed data points}
#' \item{ybar}{responses of compressed data points}
#' \item{cluster}{a vector of integers indicating assignment of each point 
#' to its nearest compressed data point}
#' \item{l2}{the total sum of squares}
#' 
#' @references 
#' Joseph, V. R. and Mak, S. (2021). Supervised compression of big data. 
#' \emph{Statistical Analysis and Data Mining: The ASA Data Science Journal}, 
#' 14(3), 217-229.
#' 
#' Beygelzimer, A., Kakadet, S., Langford, J., Arya, S., Mount, D., and Li, S. (2019). 
#' FNN: Fast nearest neighbor search algorithms and applications, R 1.1.3.
#' 
#' @export
#' 
#' @examples
#' #########################################################################
#' # One dimensional example
#' #########################################################################
#' # generate big data
#' set.seed(1)
#' N <- 3000
#' x <- seq(0,1,length=N)
#' f <- function(x) dnorm(x, mean = 0.4, sd = 0.01)
#' y <- f(x) + 0.1 * rnorm(N)
#' x <- matrix(x, ncol=1)
#' # visualize big data
#' plot(x,y,cex=.5,main="Big Data",cex.main=3,xlab="x",ylab="y",cex.lab=2, cex.axis=2)
#' # big data reduction via supercompress
#' n <- 30
#' sc <- supercompress(n,x,y,lam=0)
#' D <- sc$D # reduced data point input features
#' ybar <- sc$ybar # reduced data point response
#' points(cbind(D, ybar), pch=4,col=4,lwd=4, cex=1.5)
#' 
#' #########################################################################
#' # Two dimensional Michaelwicz function
#' #########################################################################
#' f=function(x) {
#'   p=length(x)
#'   x=pi*x
#'   val=-sum(sin(x)*(sin((1:p)*x^2/pi))^(20))
#'   return(val)
#' }
#' # generate big data
#' p=2
#' N=10000*p
#' set.seed(1)
#' x=NULL
#' for(i in 1:p) x=cbind(x,runif(N))
#' y=apply(x,1,f)+.0001*rnorm(N)
#' true=apply(x,1,f)
#' # groundtruth
#' N.plot=250
#' p1=seq(0,1,length=N.plot)
#' p2=seq(0,1,length=N.plot)
#' fc=matrix(apply(expand.grid(p1,p2),1,f),nrow = N.plot, ncol= N.plot)
#' # big data reduction via supercompress
#' n <- 100
#' sc <- supercompress(n,x,y,lam=1/(1+p))
#' D <- sc$D # reduced data point input features
#' ybar <- sc$ybar # reduced data point response
#' image(p1,p2,fc,col=cm.colors(5),xlab=expression(x[1]),ylab=expression(x[2]),
#' main="robust-supervised",cex.main=3,cex.lab=2, cex.axis=2)
#' points(D,pch=16,col=4,cex=2)
#' 
supercompress=function(n,x,y,lam=0,standardize=TRUE){
  if(!is.matrix(x)) x=matrix(x,ncol=1)
  p=dim(x)[2]
  N=length(y)
  
  # condition check
  if (n > N) stop("Number of compressed points must be smaller than the size of big data!")
  if (nrow(x) != N) stop("Number of points must be same for x and y!")
  if (lam < 0 | lam > 1) stop("Robustness parameter lam must be between 0 and 1!")
  
  if (standardize){
    mux <- apply(x,2,mean)
    sx <- apply(x,2,sd)
    muy <- mean(y)
    sy <- sd(y)
    x <- (x-rep(1,N)%*%t(mux))/(rep(1,N)%*%t(sx))
    y <- (y-muy)/sy
  }
  
  if(lam==0){
    n0=2
    km=stats::kmeans(x,n0)
    Dx=as.matrix(stats::aggregate(x,by=list(km$cluster),mean)[,-1])
    cluster=drop(knnx.index(Dx,x,k=1))
    l2=numeric(n)
    l2[1:(n0-1)]=L2(Dx,x,y)
    while(dim(Dx)[1]<n){
      i=dim(Dx)[1]
      ss=sapply(split(y,cluster),sumofsq)
      l2[i]=sum(ss)/N
      Nsub=sapply(split(rep(1,N),cluster),sum)
      ss[Nsub<3]=0
      ind=which.max(ss)
      u=matrix(x[cluster==ind,],ncol=p)
      km1=stats::kmeans(y[cluster==ind],2)
      Dx.new1=as.matrix(stats::aggregate(u,by=list(km1$cluster),mean)[,-1])
      km2=stats::kmeans(x[cluster==ind,],2)
      Dx.new2=km2$centers
      opt1=L2(Dx.new1,u,y[cluster==ind])
      #opt2=sum(sapply(split(y[cluster==ind],km2$cluster),sumofsq))/length(y[cluster==ind])
      opt2=L2(Dx.new2,u,y[cluster==ind])
      if(opt1<opt2) Dx.new=Dx.new1 else Dx.new=Dx.new2
      Dx=rbind(matrix(Dx[-ind,],ncol=p),Dx.new)
      cluster=drop(knnx.index(Dx,x,k=1))
      m=length(unique(cluster))
      if(m<dim(Dx)[1]){
        miss=(1:dim(Dx)[1])[!(1:dim(Dx)[1])%in%cluster]
        miss.near=knnx.index(x,matrix(Dx[miss,],ncol=p),k=1)
        Dx[miss,]=matrix(x[miss.near,],ncol=p)
        cluster=drop(knnx.index(Dx,x,k=1))
        m=length(unique(cluster))
        if(m<dim(Dx)[1]){
          m=length(unique(cluster))
          miss=(1:dim(Dx)[1])[!(1:dim(Dx)[1])%in%cluster]
          Dx=matrix(Dx[-miss,],ncol=p)
          cluster=drop(knnx.index(Dx,x,k=1))
        }
      }
    }
    ss=sapply(split(y,cluster),sumofsq)
    l2[n]=sum(ss)/N
  }else{
    #alpha=max(lam*p/(lam*p+1-lam),(p+1)/n)
    
    alpha=lam
    
    n0=ceiling(n*alpha)
    km=stats::kmeans(x,n0, iter.max = max(10,n0))
    Dx=km$centers
    cluster=km$cluster
    l2=numeric(n)
    y1=sqrt((1-lam)/lam)*y
    z=cbind(x,y1)
    l2[1:(n0-1)]=sum(as.matrix(stats::aggregate(z,by=list(cluster),sumofsq)[,-1]))
    while(dim(Dx)[1]<n){
      i=dim(Dx)[1]
      ss=apply(as.matrix(stats::aggregate(z,by=list(cluster),sumofsq)[,-1]),1,sum)
      l2[i]=sum(ss)/N
      Nsub=sapply(split(rep(1,N),cluster),sum)
      ss[Nsub<3]=0
      ind=which.max(ss)
      km=stats::kmeans(y[cluster==ind],2)
      Dx.new1=as.matrix(stats::aggregate(matrix(x[cluster==ind,],ncol=p),by=list(km$cluster),mean)[,-1])
      Dx.new2=stats::kmeans(x[cluster==ind,],2)$centers
      opt1=L2t(D=Dx.new1,z=z[cluster==ind,])
      opt2=L2t(D=Dx.new2,z=z[cluster==ind,])
      if(opt1<opt2) Dx.new=Dx.new1 else Dx.new=Dx.new2
      Dx=rbind(matrix(Dx[-ind,],ncol=p),Dx.new)
      cluster=drop(knnx.index(Dx,x,k=1))
      m=length(unique(cluster))
      if(m<dim(Dx)[1]){
        miss=(1:dim(Dx)[1])[!(1:dim(Dx)[1])%in%cluster]
        miss.near=knnx.index(x,matrix(Dx[miss,],ncol=p),k=1)
        Dx[miss,]=matrix(x[miss.near,],ncol=p)
        cluster=drop(knnx.index(Dx,x,k=1))
        m=length(unique(cluster))
        if(m<dim(Dx)[1]){
          m=length(unique(cluster))
          miss=(1:dim(Dx)[1])[!(1:dim(Dx)[1])%in%cluster]
          Dx=matrix(Dx[-miss,],ncol=p)
          cluster=drop(knnx.index(Dx,x,k=1))
        }
      }
    }
    ss=sum(as.matrix(stats::aggregate(z,by=list(cluster),sumofsq)[,-1]))
    l2[n]=sum(ss)/N
  }
  ybar=sapply(split(y,cluster),mean)
  
  if (standardize){
    ybar <- ybar*sy+muy
    Dx <- Dx*(rep(1,n)%*%t(sx))+rep(1,n)%*%t(mux)
  }
  
  return(list(D=Dx,ybar=ybar,cluster=cluster,l2=l2))
}
