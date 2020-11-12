MA_components <- function (AA,maxlag){
# AA must be a list
k<-nrow(AA[[1]])
p<-length(AA)
AA_n<-as.list(1:maxlag)
for(i in 1:maxlag){AA_n[[i]]<-matrix(0,nrow=k,ncol=k)}
for(i in 1:p){AA_n[[i]]<-AA[[i]]}
Fi<-as.list(1:(maxlag+1))
for(i in 1:(maxlag+1)){Fi[[i]]<-matrix(0,nrow=k,ncol=k)}
Fi[[1]]<-diag(k)
Fi[[2]]<-AA_n[[1]]
for (i in 2:maxlag){
for (j in 1:i){
FF<- Fi[[i+1-j]]%*%AA_n[[j]]
Fi[[i+1]]<-Fi[[i+1]]+FF
}
}
#Fi<-Fi[-1]
Fi
}


