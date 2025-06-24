#rm(list=ls())

#where are we working?
setwd("C:/Users/nicholls.NICHOLLS2389/Documents/collab - Kate/Oslo")

#so we can reset the plot window
gpars<-par()
i=match(c("cin","cra","csi","cxy","din","page"),names(gpars))
gpars<-gpars[-i]

#load libraries
library(MASS)
library(mnem)       #needed for transitive.X
library(igraph) 
library(Rgraphviz)
library(graph) 
library(coda)       #for MCMC output analysis
library(lecount)    #count linear extensions of a PO
library(partitions) #used by dimension functions

source(file="pofun.R")
source(file="pofun_aux.R")
source(file="modelfun.R")
source(file="outputfun.R")
source(file="dimfun.R")

################################################################
#PO basics

#start by generating a random PO with n items 
#we just need an example to look at and discuss
set.seed(10); 
n=8
h<-rZPO(n)   #we will come back to what is going on here
h

#how many relations in h?
sum(h)
#how many pairs of items have no relation?
choose(n,2)-sum(h)

#display PO h as a directed graph
showDAG(h,edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)
#this is the transitive closure as it includes edges implied by transitivity
title('transitive closure of h')

#find the list of top nodes - no "in" edges
(top=which(apply(h,2,sum)==0))

#find the list of bottom nodes - no "out" edges
(bot=which(apply(h,1,sum)==0))

#drop edges implied by transitivity for easier viewing
hr<-transitive.reduction(h)
showDAG(hr,edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)
title('transitive reduction of h')

#make a plot highlighting the redundant edges and marking top, bottom
vcols=rep(NA,n); 
vcols[top]<-'lightgreen'; vcols[bot]<-'pink'

#igraph orders the edges by first vertex
i=which(t(h)==1)
j=which(t(h)!=t(hr))
j.in.i=match(j,i)
ecols=rep('black',length(i))
ecols[j.in.i]<-'lightblue'

showDAG(h,vertex.color=vcols,edge.color=ecols,vertex.size=25,
        edge.arrow.size=0.5)
legend('topright',cex=0.9,lty=c(1,1,NA,NA),pch=c(NA,NA,16,16),col=c('black','lightblue','lightgreen','pink'),lwd=c(2,2,NA,NA),
       legend=c('reduction','closure','max set','min set'),bty="n")

################################################################
# POs and LEs
par(mfrow=c(1,2))
n=4
b<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),n,n,byrow=TRUE)
row.names(b)<-colnames(b)<-1:n
showDAG(b,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)
title('example partial order')
#get all the LEs of hc and check their intersection gives us back the original PO hc
(le<-le.expand(b)); 
#number of linear extensions
dim(le)[2]; 
#direct count of nle's without (explicit) enumeration
nle(b); 
#intersect the LEs and check you get back the PO
bc<-intersect.TO(le)
showDAG(bc,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)
title('intersection of its LEs')
check.dec(b,le) #this checks what we already see - bc, b obviously equal
par(gpars)

##################################################
#a little bit about PO dimension

#the dimension of a PO is the smallest number of its LEs 
#that intersect to give the PO

#calculating dimension is "hard" - our code calculates dim(PO)
#in reasonable time for any PO up to n=7 then dies - it simply 
#checks every subset of the LEs of size less or equal floor(n/2).
#It can only handle POs on n above 7 if they dont have many LEs.

le
(b.realiser<-decompose(b))
all(intersect.TO(b.realiser)==b) #intersect realising LEs gives back PO

dimension(b)

#the "crown" PO on n items has dim floor(n/2) - max possible
cpo=crown.PO(6) 
showDAG(cpo,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)
le<-le.expand(cpo)
le #48 LEs
(cpo.realiser<-decompose(cpo)) #takes a couple of seconds
all(intersect.TO(b.realiser)==b)

#knowing max dim PO is floor(n/2) is useful when we make a
#prior for POs - if it generates POs by taking K random total
#orders and intesecting them then we can represent any PO
#on n items using K=floor(n/2) total orders

##############################################################

# random linear extensions as realisation of a queue

#suppose actors 1,2,...,5 are in a queue contrained by a simple PO
n=5
h=matrix(c(0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0),n,n,byrow=TRUE)
colnames(h)<-rownames(h)<-1:n
showDAG(transitive.reduction(h),edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)

#set.seed(10); n=8; h<-rZPO(n) #try this - explore random PO

(le=le.expand(h))

#simulate the evolving queue - nbrs swap positions if OK by PO
q=le[,1] #start state must be a LE
T=10000  #simulate T steps
X<-matrix(NA,n,T) #store the state of the queue
for (t in 1:T) {
  i=sample(1:(n-1),1)         #propose to swap q[i] and q[i+1]
  if (h[q[i],q[i+1]]!=1) {    #if q[i+1] isnt below q[i] in the PO
    q[c(i,i+1)]=q[c(i+1,i)]   #swap q[i] and q[i+1]
  }
  X[,t]=q                     #save the current queue state
}
#go through the sampled queue realisations and label them by their LE
ord=apply(X,2,function(x){which(apply(le,2,function(y){all(x==y)}))})
#the sampled queue realisaitons should be uniformly distributed over LEs
ft<-table(ord)/T

le.str<-apply(le,2,function(x){paste(x,collapse = ' ')})
barplot(ft~le.str,las=2,cex.names=0.5,
        ylab='probability',xlab='',
        main='distribution of random orders - noise free case',
        cex.main=0.8,cex.axis=0.7);
abline(h=1/nle(h),col=2,lwd=2,lty=2)

##############################################################
#the MLE PO in the noise free case is the intersection order
#if we start with an unknown true PO "h.true" and have (noise
#free) observaitons y_i~U(LE(h.true)) then the interection order
#of y=(y_1,...,y_N) is the PO that maximises the probability
#to see the lists y we observed : h.hat = argmax_h p(y|h)

set.seed(10)
n=8
h.true<-rZPO(n)       #make a random PO with n items
par(mfrow=c(1,2)); 
showDAG(transitive.reduction(h.true),edge.color='black',vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)

N=10 #experiment with smaller N 
Y<-matrix(NA,n,N)
for (i in 1:N) {
  Y[,i]<-PunifLE(h.true)
}
Y
h.hat<-intersect.TO(Y)
showDAG(transitive.reduction(h.hat),edge.color='black',vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)

#notice the example had 48 LEs and we only have 10 rank-lists
nle(h.true)
#the dimension was just 2
decompose(h.true) 

#Would N=2 (or in general dim(h.true)) observed lists be enough?
#Repeat the exercise with less data - say N=2 or N=5 - you should find
#the MLE is sometimes wrong - not every subset of LEs of size
#dim(h) or more is a realiser for h

##############################################################
#suborders

#suppose our assessor has an unknown true PO h=([n],>) 
#expressing their preferences over n items. We present
#the assessor with subsets S_i i=1,...,N of the n items
#(called choice sets) and ask them to order the items 
#in each subset/choice set

#they will use the suborder for the items in the choice 
#set to form the order so in the noise free case h_i=(S_i,>)
#and y_i~U(LE(h_i))

set.seed(11)
N=10 #try increasing N - how big does it have to be?

#create N suborders of h.true from the last example
#first create N random-sized subsets of [n]
size=sample(2:n,N,replace=TRUE) 
S<-lapply(size,function(x){sample(1:n,x,replace=FALSE)})

#now pull out the suborders and display them
h.true.sub<-lapply(S,function(x){suborder(h.true,x)})
showDAGs(B=1,T=N,PO=h.true.sub)
par(gpars)

#now the observed lists have varying length
Y<-vector('list',N)
for (i in 1:N) {
  Y[[i]]<-PunifLE(h.true.sub[[i]])
}
Y

#it no longer holds that the intersection order is the MLE the 
#problem is that intersection isnt defined for POs on different
#ground sets of items 

#we can take all relations displayed and not contradicted by the data
#this will give a PO that admits all the lists as LEs. It isnt the
#MLE as it doesnt include relations not displayed and not contradicted 
#this will converge to the true PO in the limit that every pair of 
#items appears in an infite number of choice sets

h.est<-intersect.SO(Y)

par(mfrow=c(1,2)); 
showDAG(transitive.reduction(h.true),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
showDAG(transitive.reduction(h.est),edge.color='black',vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)

#try increasing N - now that we only have short preference orders
#in our data it takes a larger N to get a decent estimate of h.true

##############################################################

#the bad news: calculating the number of LEs of a PO is #P
#eg - calculate |LE(h_n)| where h_n is the crown PO on n items

nv<-seq(35,47,2); #unwise to try above 50 (memory)
nvl=length(nv); rt<-numeric(nvl)
pb<-txtProgressBar(min=1, max=nvl, style=3)
for (i in 1:nvl) {
  h<-crown.PO(nv[i])
  rt[i]<-system.time({nle(h)})[3] #how long to calculate #LEs?
  setTxtProgressBar(pb, i)
  flush.console()
}
par(gpars)
plot(nv,rt,log='y',type='l',ylab='run time',xlab='number of items'); 
points(nv,rt,pch=4)

#with careful coding exact analysis (ie converging MCMC estimates)
#is possible for up to around 20-30 items

#analysis of POs scales linearly with the number of lists but
#exponentially with the number of items in the longest list


##############################################################
#sample a PO UAR with n nodes and note informative depth dbn

set.seed(1)
n<-20

# N=number of samples, S=subsample step, T=number of steps
N=100
S=n^2
T=N*S

#run the MCMC targeting uniform dbn on POs with n nodes
POu<-rupo(n,N,S,T,DRAW=TRUE)

#if you try this with n=50 it will take about an hour to run
#but the convergence to depth 3 is more dramatic

#show some pictures of 10 sampled states
s=10; step=floor(N/s); 
showDAGs(B=1,T=N,poi=seq(step,N,step),PO=POu)

Du=sapply(POu,dagdepth) #calcuate the depth for each sample from MCMC
effectiveSize(Du)       #if interested in MCMC mixing and convergence 
#the start state was a total order - you see the depth coming down to small values
par(gpars); plot(Du,type='l')

hist(Du,breaks=seq(0.5,n+0.5,1),xlab='PO depth',
     main='Depth dbn, uniform PO prior',freq=FALSE) 
#concentrates on low depth, asymptotically in n the dpeth is always 3
table(Du)/length(Du)


#######################################################
#exploring the latent variable prior
# make the example 5 element PO used as a simple example

set.seed(6)
n=6 #number of items/nodes in PO
K=3 #number of features (columns of U/Z)

#just the basic PO from U without covariates
(rho=rRprior())  #the correlation within of rows of U
Sig=matrix(rho,K,K); diag(Sig)=1; #correlation matrix
U=mvrnorm(n,rep(0,K),Sig) #the feature matrix U, one row for each item
vp<-latent2order(U) #convert each column k+1:K to a total order, ranking by U[,k]
mu<-order2partial(vp,n); #intersect the column orders 
muc<-my.transitive.closure(mu) #take the transitive closure
mur<-transitive.reduction(mu)  #and transitive reduction

#plot the rows of U as paths (left) and the resulting PO (right) 
par(mfrow=c(2,2),mai=c(0,0,0,0))
plot(0,0,xlim=c(1,K+1),ylim=c(-4,4),type='n',axes=FALSE,ann=FALSE)
cl=0; apply(U,1,function(x){lines(x,lwd=2,col=(cl<<-cl+1))})
title(expression('U'),line=-3)
legend('topright',legend=1:n,lty=1,col=1:n,lwd=2)
showDAG(mur,vertex.color=NA,vertex.size=25)
title(expression('h(U)'),line=-6)

#now add covariates (just simulating the effects X*beta for the sake of example)
beta=rnorm(n)
Z=U+beta #features with covariate additive effect offset
vp<-latent2order(Z)
mz<-order2partial(vp,n); 
mzc<-my.transitive.closure(mz)
mzr<-transitive.reduction(mz)

plot(0,0,xlim=c(1,K+1),ylim=c(-4,4),type='n',axes=FALSE,ann=FALSE)
cl=0; apply(Z,1,function(x){lines(x,lwd=2,col=(cl<<-cl+1))})
title(expression(paste("Z=U+X",beta)),line=-3)
showDAG(mzr,vertex.color=NA,vertex.size=25)
title(expression('h(Z)'),line=-6)
par(gpars)

#how does the depth distribution look for the latent variable prior
n=20
K=10

N=1000
POl<-vector('list',N)

pb<-txtProgressBar(min=0, max=N, style=3)
for (t in 1:N) {
  # rho=rRprior(fac=1/16)
  # Sig=matrix(rho,K,K); diag(Sig)=1;
  # U=mvrnorm(n,rep(0,K),Sig)
  beta=rnorm(n)
  # Z=U+beta
  # vp<-latent2order(Z)
  # mu<-order2partial(vp,n); 
  POl[[t]]<-my.transitive.closure(rZPO(n,K,b=beta))
  setTxtProgressBar(pb, t/N); flush.console()
}

Dl=sapply(POl,dagdepth)
hist(Dl,breaks=seq(0.5,n+0.5,1),xlab='PO depth',
     main='Depth dbn, latent feature PO prior',freq=FALSE)
table(Dl)/length(Dl)
#depth is not unform a priori but we have removed some of the 
#biasing effect wrt depth

#compare prior depth distributions
plot(density(Du,bw=0.5,from=1,to=n),xlim=c(0,n+1),
     main='prior depth dbns',xlab='depth')
lines(density(Dl,bw=0.5,from=1,to=n),col=2)
abline(v=c(1,n),col=3)
legend('topright',legend=c('Uniform','Latent'),col=c(1,2),lty=c(1,1),lwd=2)

################################################

#sampling the queue jumping observation model

#it is likely the lists in the data dont perfectly respect the PO
#there might be recording errors or maybe occasionaly the actors
#or assessor disregard the constraint from the PO

#go back to simple example
n=4
h<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),n,n,byrow=TRUE)
row.names(h)<-colnames(h)<-1:n
showDAG(h,edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)

#generate T lists using queue-jumping noise
T=10000;
X=matrix(NA,n,T)

p.val=0.5 #the probability the next item in the list is chosen randomly
for (t in 1:T) {
  X[,t]=PunifLE(h,p=p.val)
}

#now we add noise the observed list can be any permutation of 1:n
perm=t(enum.seq(1:n)) #list all perms of 1:n
perm.str=apply(perm,2,function(x){paste(x,collapse = ' ')}) #turn them into strings

#which of the permutations correspond to LEs of h - there will
#have higher probability of being chosen when p is close to zero
le=le.expand(h)
is.le=1+apply(perm,2,function(x){any(apply(le,2,function(y){all(x==y)}))})
cols=c('lightgrey','lightblue')[is.le]
#color permutations that are LEs in blue, the rest in grey

#go through all samples - for each one identify it as a permutation
ord=apply(X,2,function(x){which(apply(perm,2,function(y){all(x==y)}))})
ft=table(ord)/T
barplot(ft~perm.str,las=2,cex.names=0.7,col=cols,
        ylab='probability',xlab='permutation',
        main=paste('distribution of random orders with noise probability = ',p.val),
        cex.main=0.8)
legend('topright',legend=c('LE of h','not LE'),
       col=c('lightblue','lightgrey'),pch=c(15,15))

#as p.val gets bigger you will see the probability mass moving away
#from the blue bars (the LEs of h) and getting smeared out across
#all the permutations - still the signal from the PO is still quite 
#clear even when the probability for noise is quite large, like 0.5

