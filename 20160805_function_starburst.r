################################################################################
### Function - Starburst
################################################################################
##Author: Mattias Aine (mattias.aine@med.lu.se)
##Date: 2016-08-05

##Get function
#source("C:/path/to/function/20160805_function_starburst.r")

##Usage:
#Starburst(p.a=NULL,     #p-values for diff(a-b) for data type A
#p.b=NULL,               #p-values for diff(a-b) for data type B
#fc.a=NULL,              #fold change for diff(a-b) for data type A (only sign of change used)
#fc.b=NULL,              #fold change for diff(a-b) for data type B (only sign of change used)
#sigLine=.01,            #Gene names, both data sets A and B need to match 100% in by row..
#gnames=NULL,            #plot names of significantly differing genes -> can be nessy with many siggenes
#plotNames=TRUE,         #your chosen significance cutoff, will be log10-converted
#MainText="MyTestplot",  #parameter for plot title
#PlotName="Testplot",    #name of created figure
#Vector=TRUE             #output as pdf (TRUE) or tif (FALSE)
#)

##The function:
Starburst<-function(p.a=NULL,p.b=NULL,fc.a=NULL,fc.b=NULL,sigLine=.01,gnames=NULL,plotNames=TRUE,MainText="MyTestplot",PlotName="Testplot",Vector=TRUE) {

  ##convert to log10(p)
  logp1<- -log10(p.a) #make sure there are no "zero" p-values or this will go to infinity..
  logp2<- -log10(p.b)
  ##get fold change direction
  fc1<-sign(fc.a)
  fc1[fc.a==0]<-1
  fc2<-sign(fc.b)
  fc2[fc.b==0]<-1
  ##modify log(p) based on FC
  logp1[fc1 ==1] <- logp1[fc1 ==1] * -1
  logp2[fc2 ==1] <- logp2[fc2 ==1] * -1

  xrange<-range(logp1)+c(-.5,.5)
  yrange<-range(logp2)+c(-.5,.5)
  xrange<-rep(max(abs(xrange)),2)*c(-1,1)
  yrange<-rep(max(abs(yrange)),2)*c(-1,1)
  
  siggenes<- abs(logp1)> -log10(sigLine) & abs(logp2)> -log10(sigLine)

  if(Vector) { #use vector graphics?
    pdf(paste(PlotName,".pdf",sep=""),height=12,width=12,useDingbats=F) ##adjust aspects by tweaking width/height
  } else {
    tiff(paste(PlotName,".tif",sep=""),height=12,width=12,units="in",res=600,compression="lzw") ##adjust aspects by tweaking width/height+res
  }
  par(fig=c(0,1,0,1),mar=c(6,5,3,1)+.1,font=2,font.lab=2,font.axis=2,font.sub=2,lwd=6,cex.main=2,cex.axis=2,cex.lab=2,cex.sub=1.5)
  plot(1,type="n",xlim=xrange,ylim=yrange,axes=F,
  xlab="Delta(Meth)>=0 -log10(p,fdr) or Delta(Meth)<0 log10(p,fdr)",  #
  ylab="Delta(GEX)>=0 -log10(p,fdr) or Delta(GEX)<0 log10(p,fdr)",    #
  sub="--- significance threshold",main=MainText
  )
  axis(1,lwd=6,las=1)
  axis(2,lwd=6,las=1)
  lines(x=rep(-log10(sigLine),2),y=c(-log10(sigLine),yrange[2]+1),lty=2)
  lines(x=rep(log10(sigLine),2),y=c(-log10(sigLine),yrange[2]+1),lty=2)
  lines(y=rep(-log10(sigLine),2),x=c(-log10(sigLine),xrange[2]+1),lty=2)
  lines(y=rep(log10(sigLine),2),x=c(-log10(sigLine),xrange[2]+1),lty=2)
  lines(y=rep(-log10(sigLine),2),x=c(log10(sigLine),xrange[1]-1),lty=2)
  lines(y=rep(log10(sigLine),2),x=c(log10(sigLine),xrange[1]-1),lty=2)
  lines(x=rep(-log10(sigLine),2),y=c(log10(sigLine),yrange[1]-1),lty=2)
  lines(x=rep(log10(sigLine),2),y=c(log10(sigLine),yrange[1]-1),lty=2)
  box()

  points(logp1,logp2,pch=16,cex=.5)
  points(logp1[siggenes],logp2[siggenes],pch=16,cex=1,col=2)
  if(plotNames) {
  text(logp1[siggenes & fc1 ==1],logp2[siggenes & fc1 ==1],,labels=gnames[siggenes & fc1 ==1],col=2,pos=2)
  text(logp1[siggenes & fc1 ==-1],logp2[siggenes & fc1 ==-1],,labels=gnames[siggenes & fc1 ==-1],col=2,pos=4)
  }
  text(x=xrange[1],y=yrange[2],labels="Methylated/Repressed",pos=4,cex=1.5)
  text(x=xrange[1],y=yrange[1],labels="Methylated/Expressed",pos=4,cex=1.5)
  text(x=xrange[2],y=yrange[2],labels="Demethylated/Repressed",pos=2,cex=1.5)
  text(x=xrange[2],y=yrange[1],labels="Demethylated/Expressed",pos=2,cex=1.5)
  dev.off()
}
###END