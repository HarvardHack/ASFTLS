ADVANCED STATISTICS FOR THE LIFE SCIENCES: GSE5859Subset

WEEK 1

# q1: 5

sum(sampleInfo$date=="2005-06-27")


# q2: 21

sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

#q3: 8.233599

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

# q4: 0.041

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)
mean(pvals < 0.05)

# q5: 0.008

mean(pvals < 0.01)  #  using code from previous question

# q6: 1

set.seed(100)
pvals<- replicate(20,{
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.value
})
sum(pvals<=0.05)

# q7: 1.007

set.seed(100)
B = 1000
plessthan = replicate(B,{
  pvals = replicate(20,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.value
    })
  sum(pvals<=0.05)
})
table(plessthan) ##just for illustration
mean(plessthan)

# q8: 0.646

mean(plessthan>0)

# q9: "the identity line"

# q10: 1

B<-1000
minpval <- replicate(B, min(runif(8793,0,1))<0.05)
mean(minpval>=1)

# q11: 0.000005833407

B=10000
cutoffs = 10^seq(-7,-4,0.1) 
prob = sapply(cutoffs,function(cutoff){
    minpval =replicate(B, min(runif(8793,0,1))<=cutoff)
    mean(minpval>=1)
    })
cutoffs[which.min(abs(prob-0.05))]

# q12: Bonferroni's

alphas <- seq(0,0.25,0.01)
par(mfrow=c(2,2))
for(m in c(2,10,100,1000)){
  plot(alphas,alphas/m - (1-(1-alphas)^(1/m)),type="l")
  abline(h=0,col=2,lty=2)
}


# q13: 0.0467

set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- alpha/m
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

# q14: 0.0475

set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- (1-(1-alpha)^(1/m))
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

# q15: 1383

g <- factor(sampleInfo$group)
pvals = rowttests(geneExpression,g)$p.value
sum(pvals<0.05)

# q16: 10

k = 0.05/length(pvals)
sum(pvals<k)

# q17: 13

g = factor(sampleInfo$group)
pvals = rowttests(geneExpression,g)$p.value
fdr = p.adjust(pvals,method="fdr")
sum(fdr<0.05)

# q18: 17

library(qvalue)
res = qvalue(pvals)
qvals = res$qvalues
sum(qvals<0.05)

# q19: 0.6855371

qvalue(pvals)$pi0

# q20: "The qvalue function estimates the proportion of genes for which the null hypothesis is true and provides a less conservative estimate"

plot(qvalue(pvals)$qvalue/p.adjust(pvals,method="fdr"))
abline(h=qvalue(pvals)$pi0,col=2)
hist(pvals,breaks=seq(0,1,len=21))
expectedfreq <- length(pvals)/20 
abline(h=expectedfreq*qvalue(pvals)$pi0,col=2,lty=2)

# q21: 0.000005305679

set.seed(1)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)
  FN1 <- sum(pvals[1:positives]>0.05/m)
  qvals1 = p.adjust(pvals,method="fdr")
  FP2 <- sum(qvals1[-(1:positives)]<=0.05)
  FN2 <- sum(qvals1[1:positives]>0.05)
  qvals2 = qvalue(pvals)$qvalue
  FP3 <- sum(qvals2[-(1:positives)]<=0.05)
  FN3 <- sum(qvals2[1:positives]>0.05)  
  c(FP1,FN1,FP2,FN2,FP3,FN3)
  })

# q22: 0.763252

mean(result[2,]/(positives))

# q23: 0.002737851

mean(result[3,]/(m-positives))

# q24: 0.081418

mean(result[4,]/(positives))

# q25: 0.002926082

mean(result[5,]/(m-positives))

# q26: 0.07737

mean(result[6,]/(positives))
 
# q27: 0.9475336

mean(e[,1]<k & e[,2]<k)

# q28: "The tails do not dominate the plot"

# q29: 0.7767887

sd(e[,2]-e[,1]) 

# q30: 3057

sum(abs(e[,2]-e[,1])>1)
