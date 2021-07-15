################################################################################
### Test - Starburst function
################################################################################
##Author: Mattias Aine (mattias.aine@med.lu.se)
##Date: 2016-08-05

##Get function
source("C:/path/to/function/20160805_function_starburst.r")

##Exaggerated Mock data -> you will need to figure out which data type is "a" and which is "b" for everything to be right..
set.seed(20160805)

a<-matrix(rnorm(200,0,2),ncol=10) ##mock methylation

b<-matrix(rnorm(200,0,2),ncol=10) ##mock gex

a<-rbind(
  cbind(matrix(rnorm(40,5,2),ncol=5),matrix(rnorm(40,-5,2),ncol=5)),
  a,
  cbind(matrix(rnorm(40,-5,2),ncol=5),matrix(rnorm(40,5,2),ncol=5)),
  cbind(matrix(rnorm(20,-5,2),ncol=5),matrix(rnorm(20,5,2),ncol=5)),
  cbind(matrix(rnorm(20,5,2),ncol=5),matrix(rnorm(20,-5,2),ncol=5))
  )
colnames(a)<-paste(rep(c("x","y"),each=5),1:5,sep="")
rownames(a)<-c(paste("upLeft",1:8,sep=""),paste("noDiff",1:20,sep=""),paste("dnRight",1:8,sep=""),paste("upRight",1:4,sep=""),paste("dnLeft",1:4,sep=""))

b<-rbind(
  cbind(matrix(rnorm(40,-5,2),ncol=5),matrix(rnorm(40,5,2),ncol=5)),
  b,
  cbind(matrix(rnorm(40,5,2),ncol=5),matrix(rnorm(40,-5,2),ncol=5)),
  cbind(matrix(rnorm(20,-5,2),ncol=5),matrix(rnorm(20,5,2),ncol=5)),
  cbind(matrix(rnorm(20,5,2),ncol=5),matrix(rnorm(20,-5,2),ncol=5))
  )


colnames(b)<-paste(rep(c("x","y"),each=5),1:5,sep="")
rownames(b)<-c(paste("upLeft",1:8,sep=""),paste("noDiff",1:20,sep=""),paste("dnRight",1:8,sep=""),paste("upRight",1:4,sep=""),paste("dnLeft",1:4,sep=""))

a[c(1:2,10:11,35:36),c(1:2,6:7)]
#                x1         x2        y1         y2
#upLeft1   3.968939  5.1015170 -4.083808 -2.3383394
#upLeft2   4.924658  4.9003711 -3.667048 -1.6516695
#noDiff2  -4.163530 -0.6336259 -2.984229  0.6770169
#noDiff3   2.162231  2.6772234 -3.082869  0.3929817
#dnRight7 -6.774964 -7.7278000  3.316252  5.6230534
#dnRight8 -6.533722 -2.3497151  2.951716  5.2231327

b[c(1:2,10:11,35:36),c(1:2,6:7)]
#                 x1        x2        y1        y2
#upLeft1  -5.2998658 -6.386704  7.363402  4.540020
#upLeft2  -5.2233814 -4.466345  7.363536  5.134730
#noDiff2   0.9822422  2.039219 -1.960021  1.845007
#noDiff3  -2.0256301 -2.410593 -3.857819 -1.790918
#dnRight7  5.6918021  4.384387 -3.100640 -3.453930
#dnRight8  4.7609186  4.645154 -6.823446 -2.366537

##Do sigtests (not pretty but works, get p-value and fold change however you wish)
testClass<-relevel(factor(sub("\\d","",colnames(a))),"y")  #define your test classes so that the delta becomes "right"

#testClass<-factor(sub("\\d","",colnames(a))) #flips the plot..

testClass
# [1] x x x x x y y y y y
#Levels: y x

sigTestA<-t(apply(a,1,function(x) {
  as.vector(summary(lm(x~testClass))$coeff[2,c(1,4)]) #gets folddiff and p(t) from linear model
  }))

sigTestA[c(1:2,10:11,35:36),]
#                [,1]         [,2]
#upLeft1    8.8055414 7.360001e-06
#upLeft2    9.0868214 1.044334e-04
#noDiff2   -0.8869128 5.254326e-01
#noDiff3    2.0318152 1.465067e-01
#dnRight7 -10.9214107 2.309634e-05
#dnRight8  -9.1026476 1.538005e-04

sigTestB<-t(apply(b,1,function(x) {
  as.vector(summary(lm(x~testClass))$coeff[2,c(1,4)]) #gets folddiff and p(t) from linear model
  }))

sigTestB[c(1:2,10:11,35:36),]
#                [,1]         [,2]
#upLeft1  -11.1707352 2.055820e-05
#upLeft2  -11.0017631 4.014330e-06
#noDiff2    0.3489725 8.094891e-01
#noDiff3    0.3929895 7.844901e-01
#dnRight7   9.7910012 1.024124e-06
#dnRight8   9.4344757 4.931708e-06

##Will use raw p-values but should use adjusted from e.g. "p.adjust"-function

##Running the function..
Starburst(p.a=sigTestA[,2],  #p-values for diff(a-b) for data type A
p.b=sigTestB[,2],            #p-values for diff(a-b) for data type B
fc.a=sigTestA[,1],           #fold change for diff(a-b) for data type A (only sign of change used)
fc.b=sigTestB[,1],           #fold change for diff(a-b) for data type B (only sign of change used)
gnames=rownames(sigTestA),   #Gene names, both data sets A and B need to match 100% in by row..
plotNames=TRUE,              #plot names of significantly differing genes -> can be nessy with many siggenes
sigLine=.01,                 #your chosen significance cutoff, will be log10-converted
MainText="MyTestplot",      #parameter for plot title
PlotName="Testplot",        #name of created figure
Vector=TRUE                 #output as pdf (TRUE) or tif (FALSE)
)

Starburst(p.a=sigTestA[,2],  #p-values for diff(a-b) for data type A
p.b=sigTestB[,2],            #p-values for diff(a-b) for data type B
fc.a=sigTestA[,1],           #fold change for diff(a-b) for data type A (only sign of change used)
fc.b=sigTestB[,1],           #fold change for diff(a-b) for data type B (only sign of change used)
gnames=rownames(sigTestA),   #Gene names, both data sets A and B need to match 100% in by row..
plotNames=TRUE,              #plot names of significantly differing genes -> can be nessy with many siggenes
sigLine=.01,                 #your chosen significance cutoff, will be log10-converted
MainText="MyTestplot",      #parameter for plot title
PlotName="Testplot",        #name of created figure
Vector=FALSE                 #output as pdf (TRUE) or tif (FALSE)
)

###END















###END
