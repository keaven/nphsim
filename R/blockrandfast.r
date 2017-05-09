blockrandfast <- function(n, num.lebvels=2, levels=c("Control","Experimental"), 
                          block.size=1)
{ block.num<-ceiling(n/block.size/length(levels)) #number of blocks
  block.length<-length(levels)*block.size #Length of each block
  block.content<-rep(levels,block.size) # block content trt, control etc
  n.id<-block.num*block.length
  tmp<-data.frame(id=1:n.id,rand=runif(n.id),trt=rep(block.content,block.num))#genderate random number
  tmp$block<-floor((tmp$id-1)/block.length)+1 #generate block 
  treatment<-tmp[order(tmp$block,tmp$rand),]$trt #generate block randomization 
  out<-data.frame(treatment=treatment) #output
  out$block.size<-block.size
  return(out)
}