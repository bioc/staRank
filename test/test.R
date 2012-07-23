# simple test case data sets to run stability on
# 
# Author: julianes
###############################################################################


# generate dataset: all stable
d<-matrix(c(1:10,1:10,1:10,1:10),ncol=4)
rownames(d)<-letters[1:10]
# get stable Ranking, stable setsizes and the Pi matrix for default settings
s<-getStability(d,0.9)
# get samples manually
sr<-replicate(100,getSampleScores(d,'mean',decreasing=TRUE,bootstrap=TRUE))
sr<-apply(sr,2,function(x){rank(x,na.last='keep')})

# generate dataset: all stable reverse
d<-matrix(c(10:1,10:1,10:1,10:1),ncol=4)
rownames(d)<-letters[1:10]
# get stable Ranking, stable setsizes and the Pi matrix for default settings
s<-getStability(d,0.9)
# get samples manually
sr<-replicate(100,getSampleScores(d,'mean',decreasing=TRUE,bootstrap=TRUE))
sr<-apply(sr,2,rank)

# generate dataset: a bit unstable
d<-matrix(c(1,1,2,2,2,2,1,1,3,3,4,4,4,4,3,3,6,6,5,5,5,5,6,6,
				rep(7:10,each=4)),ncol=4,byrow=TRUE)
rownames(d)<-letters[1:10]
# get stable Ranking, stable setsizes and the Pi matrix for default settings
s<-getStability(d,0.9)
# get samples manually
sr<-replicate(100,getSampleScores(d,'mean',decreasing=TRUE,bootstrap=TRUE))
sr<-apply(sr,2,rank)

# generate dataset: random
d<-replicate(4,sample(1:10,10,replace=FALSE))
rownames(d)<-letters[1:10]
# get stable Ranking, stable setsizes and the Pi matrix for default settings
s<-getStability(d,0.9)
# get samples manually
sr<-replicate(100,getSampleScores(d,'mean',decreasing=TRUE,bootstrap=TRUE))
sr<-apply(sr,2,rank)

# generate dataset: missing data
d<-replicate(4,sample(1:10,10,replace=FALSE))
rownames(d)<-letters[1:10]
d[4,]<-rep(NA,4)
# get stable Ranking, stable setsizes and the Pi matrix for default settings
s<-getStability(d,0.9)
# get samples manually
sr<-replicate(100,getSampleScores(d,'mean',decreasing=TRUE,bootstrap=TRUE))
sr<-apply(sr,2,rank)