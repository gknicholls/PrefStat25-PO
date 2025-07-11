###
#Practical for PrefStat25 
#Thursday July 3rd 2025
#Lecturers: Geoff Nicholls & Kate Lee
#Bayesian Inference for Partial Orders from Rank-Order Data (part I)
#15:30 – 17:00 LAB: Explore basic Partial Order properties in R

#https://github.com/gknicholls/PrefStat25-PO

#if you are unable to install lecount then there are a few
#adjustments needed below - this version of the prac doesnt use lecount

#######################################################################
### 1/12 Setup

# save default plot params so we can reset the plot window

gpars<-par()
i=match(c("cin","cra","csi","cxy","din","page"),names(gpars))
gpars<-gpars[-i]

# load libraries

library(MASS)
library(mnem)       # needed for transitive.reduction/closure
library(igraph) 
library(Rgraphviz)
library(graph) 
library(coda)       # for MCMC output analysis
#library(lecount)    # count linear extensions of a PO
library(partitions) # used by dimension functions

# where are we working?

#update this to suit yourself
wd<-"C:/Users/nicholls.NICHOLLS2389/OneDrive - Nexus365/Documents/GitHub/PrefStat25-PO/Code for Practical - Without lecount/"
setwd(wd)

# load function files

source(file="pofun.R")
source(file="pofun_aux.R")
source(file="modelfun.R")
source(file="dimfun.R")


#######################################################################
### 2/12 Partial Order basics

# Start by generating a random PO with M items - we just need an 
# example to look at and discuss

set.seed(10); 
M=8
h<-rZPO(M)   # for now this is just a random PO

# Partial orders are one to one with transitively closed directed 
# acyclic graphs. As a PO is a kind of DAG we represent it on the 
# computer using its incidence matrix.

#h[i,j]=1 if i>j
h

# How many relations/edges in h?
  
sum(h)

# How many pairs of items have no relation?
  
choose(M,2)-sum(h)

# We can visualise a PO by plotting it as a graph.

# display PO h as a directed graph

showDAG(h,edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)
title('transitive closure of h')

# This is the transitive closure as it includes edges implied 
# by transitivity.

# One basic operation is to list the items at the top and bottom 
# of the PO.

# find the list of top nodes - no "in" edges
(top=which(apply(h,2,sum)==0))

# list the bottom nodes - no "out" edges
(bot=which(apply(h,1,sum)==0))

# We also make use of the depth of a partial order. This is the length 
# of the longest chain. What is the depth of h?
  
dagdepth(h)

# The empty order has depth one and a complete order on n items has 
# depth n.

h.empty=0*h
h.total=h.empty; h.total[upper.tri(h.total)]<-1

par(mfrow=c(1,2))
showDAG(h.empty,vertex.color=NA,vertex.size=25)
title('Empty order')
showDAG(h.total,edge.color='black',vertex.color=NA,vertex.size=15,edge.arrow.size=0.25)
title('Complete order')
par(gpars)

dagdepth(h.empty) #depth is the length of longest chain, counting vertices in chain
dagdepth(h.total); dagdepth(h.total)==M

# For viewing we drop edges implied by transitivity.

hr<-transitive.reduction(h)
showDAG(hr,edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)
title('transitive reduction of h')

# Here is a plot highlighting the redundant edges and marking top and 
# bottom nodes (the max set and the min set).

vcols=rep(NA,M); 
vcols[top]<-'lightgreen'; vcols[bot]<-'pink'

# igraph orders the edges by first vertex
i=which(t(h)==1)
j=which(t(h)!=t(hr))
j.in.i=match(j,i)
ecols=rep('black',length(i))
ecols[j.in.i]<-'lightblue'

showDAG(h,vertex.color=vcols,edge.color=ecols,vertex.size=25,
        edge.arrow.size=0.5)
legend('topright',cex=0.9,lty=c(1,1,NA,NA),pch=c(NA,NA,16,16),col=c('black','lightblue','lightgreen','pink'),lwd=c(2,2,NA,NA),
       legend=c('reduction','closure','max set','min set'),bty="n")


#######################################################################
### 3/12 Intersecting orders

# We can intersect partial orders on the same choice set. The intersection 
# has all the relations that shared by the two orders.

# For example, if our two partial orders are (1\>2\>3) and (1\>3\>2), 
# ie two total orders, then we get a PO with relations 1\>2, 1\>3 
# and no other relations.

# Since total orders are just ordered lists we have a simpler notation. 
# For example the total order 1\>3\>2 is the list (1,3,2). Put the two 
# ordered lists together as columns of an nx2 matrix.

v=matrix(c(1,2,3,1,3,2),3,2)
v

# Now compute the intersection - the intersection order is 
# a partial order

hi<-intersect.TO(v)
hi

# and plot it.

showDAG(hi,edge.color='black',vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)


#######################################################################
### 4/12 Linear Extensions of a partial order

# The partial order will be the thing we dont know and want to estimate. 
# Our data will be linear extensions. These are total orders that respect 
# the partial order. The set of all LEs of a PO h is LE(h) with \|LE(h)\| 
# elements/LEs in total.

# Take a simple PO on 4 items.

par(mfrow=c(1,2))
M=4
b<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),M,M,byrow=TRUE)
row.names(b)<-colnames(b)<-1:M
showDAG(b,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)
title('example partial order')

# Can you list the LEs by hand? There are 5. We now enumerate all the LEs 
# of b and check their intersection gives us back the original PO.

(le<-le.expand(b)); 

# The number of LEs of b is the number of columns of le. A function nle() 
# directly counts the number of LEs.

# number of linear extensions

dim(le)[2]; 

# direct count of nle's without (explicit) enumeration

nle(b); 

# Intersect the LEs and check you get back the PO

bc<-intersect.TO(le) 
# the right function when the LEs have the same length - takes matrix
# like le[] as input - when the LEs have different length we put them
# in a list and intersect the using the function intersect.SO() - 
# the TO stands for total orders (ie, complete orders) and the SO stands
# for suborders

# poset bc should be the same as b
showDAG(bc,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)
title('intersection of its LEs')
par(gpars)


#######################################################################
### 5/12 PO dimension

# The dimension of a PO is the smallest number of its LEs that intersect 
# to give the PO.

# Calculating dimension is "hard" - our code calculates dim(PO) in 
# reasonable time for any PO up to M=7 then dies - it simply checks 
# every subset of the LEs of size less or equal floor(M/2).

# It can only handle POs on M above 7 if they don't have many LEs.

# Recall the LEs of the PO b above

le

# the smallest set of LEs intersecting to give back b is a "realiser"

(b.realiser<-decompose(b))       

all(intersect.TO(b.realiser)==b) # intersect LEs in realiser gives back PO

dimension(b)

# The "crown" PO on M items looks like a crown.

cpo=crown.PO(6) 
showDAG(cpo,edge.color='black',vertex.color=NA,vertex.size=45,edge.arrow.size=0.5)

# It has dimension equal floor(M/2) - the maximum possible. The crown 
# with M=6 has 48 linear extensions but you only need 3 to represent 
# this PO.

le<-le.expand(cpo)
le 
(cpo.realiser<-decompose(cpo))     # takes a couple of seconds
all(intersect.TO(cpo.realiser)==cpo)

# Knowing the maximum dimension of a PO is floor(n/2) is useful 
# when we make a prior for POs - if it generates POs by taking 
# K random total orders and intersecting them then we can represent 
# any PO on n items using K=floor(n/2) total orders.


#######################################################################
### 6/12 Data as a realisation of a queue

# Suppose actors 1,2,...,5 are in a queue contrained by a simple PO

M=5
h=matrix(c(0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0),M,M,byrow=TRUE)
colnames(h)<-rownames(h)<-1:M
showDAG(transitive.reduction(h),edge.color='black',vertex.color=NA,vertex.size=25,edge.arrow.size=0.5)

# set.seed(10); M=8; h<-rZPO(M) # try this - explore random PO

# Simulate the evolving queue - neighbors swap positions if not 
# constrained by PO

(le=le.expand(h))
q=le[,1] # start state must be a LE

T=10000  # simulate T steps
X<-matrix(NA,M,T) # store the state of the queue
for (t in 1:T) {
  i=sample(1:(M-1),1)         # propose to swap q[i] and q[i+1]
  if (h[q[i],q[i+1]]!=1) {    # if q[i+1] isnt below q[i] in the PO
    q[c(i,i+1)]=q[c(i+1,i)]   # swap q[i] and q[i+1]
  }
  X[,t]=q                     # save the current queue state
}

# the columns of X contain 10000 random orders

# The sampled queue realisations should be uniformly distributed over LEs

# Go through the sampled queue realisations (columns of X) and label them 
# by which LE they are

ord=apply(X,2,function(x){
  which(apply(le,2,function(y){all(x==y)}))
})

# how frequently did each LE appear

(ft<-table(ord)/T)

# Make a histogram. On the x-axis are all the distinct LEs and 
# the y-axis shows how frequently they appeared.

le.str<-apply(le,2,function(x){paste(x,collapse = ' ')})
barplot(ft~le.str,las=2,cex.names=0.5,
        ylab='probability',xlab='',
        main='distribution of random orders - noise free case',
        cex.main=0.8,cex.axis=0.7);
abline(h=1/nle(h),col=2,lwd=2,lty=2)


#######################################################################
### 7/12 Intersection order is MLE (sometimes)

# The MLE PO in the noise free case is the intersection order: 
# if we start with an unknown true PO "h.true" and have (noise free) 
# observations y_i~U(LE(h.true)) then the intersection order of 
# y=(y_1,...,y_N) is the PO that maximises the likelihood: 
# h.hat = argmax_h p(y\|h)

# Take a simple example PO

set.seed(10)
M=8
h.true<-rZPO(M)       # make a random PO with n items
showDAG(transitive.reduction(h.true),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
title("true PO")

# Simulate synthetic data which is N preference orders.

N=10                 # try experimenting with smaller N 
Y<-matrix(NA,M,N)
for (i in 1:N) {
  Y[,i]<-PunifLE(h.true)
}
Y

# h.hat is the MLE PO

h.hat<-intersect.TO(Y)

# compare the MLE with the truth

par(mfrow=c(1,2))
showDAG(transitive.reduction(h.hat),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
title("MLE PO")
showDAG(transitive.reduction(h.true),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
title("true PO")

# In this example the true PO had 48 LEs and dimension 2 and we 
# only observed 10 rank-order lists.

nle(h.true)
decompose(h.true) # smallest number of lists that identify the PO

# Would N=2 (or in general dim(h.true)) observed lists be enough? 
# Repeat the exercise with less data - say N=2 or N=5 - you should 
# find the MLE is sometimes wrong - not every subset of LEs of size 
# dim(h) or more is a realiser for h

#######################################################################
### 8/12 Suborders

# Suppose our assessor has an unknown true PO h=([M],\>) expressing 
# their preferences over M items. We present the assessor with 
# subsets S_i i=1,...,N of the M items (called choice sets) and ask 
# them to order the items in each subset/choice set. They return N 
# lists y_1,...,y_N with y_i a ranking of S_i.

# The assessor uses the suborder h_i=(S_i,\>) for the items in the 
# choice set to form the order so in the noise free case y_i~U(LE(h_i)).

set.seed(11)
N=10 

# Create N suborders of h.true from the last example. First create 
# N random-sized subsets of [n].

size=sample(2:M,N,replace=TRUE) 
S<-lapply(size,function(x){sample(1:M,x,replace=FALSE)})

# Now pull out the suborders and display them

h.true.sub<-lapply(S,function(x){suborder(h.true,x)})
showDAGs(B=1,T=N,PO=h.true.sub)
par(gpars)

# We will now simulate synthetic data - one linear extension for 
# each of the N suborders above.

Y<-vector('list',N)
for (i in 1:N) {
  Y[[i]]<-PunifLE(h.true.sub[[i]])
}
Y

# Now the observed lists have varying length. It no longer holds that 
# the intersection order is the MLE. we can take all relations displayed 
# and not contradicted by the data - this will give a PO that admits all 
# the lists as LEs. It isnt the MLE as it doesnt include relations not 
# displayed and not contradicted.

# This PO estimate will converge to the true PO in the limit that every 
# pair of items appears in an infinite number of choice sets.

h.est<-intersect.SO(Y)

par(mfrow=c(1,2)); 
showDAG(transitive.reduction(h.true),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
title("true PO")
showDAG(transitive.reduction(h.est),edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)
title("estimated PO")

# Try increasing N - now that we only have short preference orders 
# in our data it takes a larger N to get a decent estimate of h.true.


#######################################################################
### 9/12 Counting linear extensions for a general PO is hard

# Bad news for statistical inference: calculating the number of LEs 
# (ie, the likelihood) of a general PO is # P.

# For example, let's calculate the number of LEs nle=\|LE(h_m)\| 
# of h_m where is the crown PO on m items.

#you may need to stop the following and adjust the numbers down
#depending on the hardware you are using

mv<-seq(8,16,2);   # if no lecount()
#mv<-seq(35,47,2); # if lecount() loaded - unwise to try above 50 (memory)
mvl=length(mv); rt<-numeric(mvl)
pb<-txtProgressBar(min=1, max=mvl, style=3)
for (i in 1:mvl) {
  h<-crown.PO(mv[i])
  rt[i]<-system.time({nle(h)})[3] # how long to calculate # LEs?
  setTxtProgressBar(pb, i)
  flush.console()
}

# Time consuming! The run time grows exponentially fast.

# The runtime measured for m=8 may be zero - ignore warning
par(gpars) 
plot(mv,rt,log='y',type='l',ylab='run time',xlab='number of items'); 
points(mv,rt,pch=4)

# With careful coding, exact analysis (ie converging MCMC estimates) 
# is possible for up to around 20-30 items.

# Analysis of POs scales linearly with the number of lists but 
# exponentially with the number of items in the longest list.


#######################################################################
### 10/12 The uniform prior over partial orders

# Sample a PO UAR with M nodes and note informative depth distribution. 
# When M is large the uniform distribution concentrates on POs of 
# depth 3. In this example we only take M=20 (so we dont have to wait 
# too long for the sampler to converge). With M=20 the uniform prior 
# distribution is already quite concentrated on depth 4-6.

# The following uses MCMC to sample POs uniformly at random. 
# It initialises the MCMC state with a total order. Once in equilibrium 
# the PO depth stays around 4-6.

set.seed(1)
M<-20

#  N=number of samples, S=subsample step, T=number of steps
N=100
S=M^2
T=N*S

# run the MCMC targeting uniform dbn on POs with M nodes

POu<-rupo(M,N,S,T,DRAW=TRUE)

# If you try this with M=50 it will take about an hour to run but 
# the convergence to depth 3 is more dramatic

# Here are 10 sampled states taken at intervals along the MCMC run.

s=10; step=floor(N/s); 
showDAGs(B=1,T=N,poi=seq(step,N,step),PO=POu)

# Calculate the depth for each sample from MCMC.

Du=sapply(POu,dagdepth) 
effectiveSize(Du) # if interested in MCMC mixing and convergence 

# Plot the depth of the evolving PO state. The start state was 
# a total order - you see the depth coming down to small values.

par(gpars); plot(Du,type='l',ylab='Depth of random PO',xlab='MCMC step',
                 main='evolving depth of PO sampled U.A.R.')

hist(Du,breaks=seq(0.5,M+0.5,1),xlab='PO depth',
     main='Depth dbn, uniform PO prior',freq=FALSE) 

# It concentrates on low depth, asymptotically in n the depth is always 3.

table(Du)/length(Du)


#######################################################################
### 11/12 Exploring the latent variable prior

# The latent variable prior has three nice properties: it is 
# "marginally consistent" (AKA "projective"); it is easy to incorporate 
# covariates on the ordered items; it has a prior distribution 
# over depths that is less informative.

# The prior has a parameter, K, the embedding dimension. 
# Each ranked item i=1,...,M is assigned a real vector U_i of dimension K. 
# If i and j are two items then i->j in the PO iff 
# U_{i,k}>U_{j,k} for each k=1,...,K.

set.seed(6)
M=6 # number of items/nodes in PO
K=3 # number of features (columns of U/Z)

# Here is a PO simulated using latent variables U.

(rho=rRprior())  # the correlation within of rows of U
Sig=matrix(rho,K,K); diag(Sig)=1; # correlation matrix
U=mvrnorm(M,rep(0,K),Sig) # the feature matrix U, one row for each item
vp<-latent2order(U) # convert each column k+1:K to a total order, ranking by U[,k]
mu<-order2partial(vp,M); # intersect the column orders 
muc<-my.transitive.closure(mu) # take the transitive closure
mur<-transitive.reduction(mu)  # and transitive reduction

# Plot the rows of U as paths (left) and the resulting PO (right)

par(mfrow=c(2,2),mai=c(0,0,0,0))
plot(0,0,xlim=c(1,K+1),ylim=c(-4,4),type='n',axes=FALSE,ann=FALSE)
cl=0; apply(U,1,function(x){lines(x,lwd=2,col=(cl<<-cl+1))})
title(expression('U'),line=-3)
legend('topright',legend=1:M,lty=1,col=1:M,lwd=2)
showDAG(mur,vertex.color=NA,vertex.size=25)
title(expression('h(U)'),line=-6)

# Now add covariates (just simulating the effects X\*beta for 
# the sake of example). Each item gets a covariate. If the latent 
# feature vector for item i is U_i=(U\_{i,1},...,U\_{i,K}) and its 
# linear predictor value is alpha_i=x_i\^T beta then the latent 
# feature vector is updated to Z_i=U_i+alpha_i 1_K - a constant 
# alpha_i gets added to every entry in U_i so Z\_{i,k}=U\_{i,k}+alpha_i 
# for each k=1,...,K.

alpha=rnorm(M)
Z=U+alpha # features with covariate additive effect offset
vp<-latent2order(Z)
mz<-order2partial(vp,M); 
mzc<-my.transitive.closure(mz)
mzr<-transitive.reduction(mz)

# Now plot the PO simulated from the latent variable prior but this time 
# with covariate offsets

plot(0,0,xlim=c(1,K+1),ylim=c(-4,4),type='n',axes=FALSE,ann=FALSE)
cl=0; apply(Z,1,function(x){lines(x,lwd=2,col=(cl<<-cl+1))})
title(expression(paste("Z=U+X",beta)),line=-3)
showDAG(mzr,vertex.color=NA,vertex.size=25)
title(expression('h(Z)'),line=-10)
par(gpars)

# How does the depth distribution look for the latent variable prior? 
# For the uniform prior we saw it concentrated on a small range of (small) 
# depth values. We took M=20 in our example there, so take M=20 here.

# For the latent variable prior we use K=10 because that is n/2 
# - theory tells us that dim(h) is at most floor(n/2) so if we form h 
# by intersecting K total orders with K=floor(n/2) then our prior can 
# generate any PO on n items.

M=20
K=10

N=1000
POl<-vector('list',N)

pb<-txtProgressBar(min=0, max=N, style=3)
for (t in 1:N) {
  beta=rnorm(M)
  POl[[t]]<-my.transitive.closure(rZPO(M,K,b=beta))
  setTxtProgressBar(pb, t/N); flush.console()
}

Dl=sapply(POl,dagdepth)
hist(Dl,breaks=seq(0.5,M+0.5,1),xlab='PO depth',
     main='Depth dbn, latent feature PO prior',freq=FALSE)
table(Dl)/length(Dl)

# Depth is still not uniform a priori but we have removed some of 
# the biasing effect wrt depth we saw in the uniform PO prior.

# Compare prior depth distributions

plot(density(Du,bw=0.5,from=1,to=M),xlim=c(0,M+1),
     main='prior depth dbns',xlab='depth')
lines(density(Dl,bw=0.5,from=1,to=M),col=2)
abline(v=c(1,M),col=3)
legend('topright',legend=c('Uniform','Latent'),col=c(1,2),lty=c(1,1),lwd=2)


#######################################################################
### 12/12 The queue jumping observation model

# It is likely the lists in the data don't perfectly respect the PO 
# there might be recording errors or maybe occasionally the actors 
# or assessor disregard the constraint from the PO

# Go back to simple example

M=4
h<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),M,M,byrow=TRUE)
row.names(h)<-colnames(h)<-1:M
showDAG(h,edge.color='black',
        vertex.color=NA,vertex.size=35,edge.arrow.size=0.5)

# generate T lists using queue-jumping noise

T=10000;
X=matrix(NA,M,T)

p.val=0.5 # the probability the next item in the list is chosen randomly
for (t in 1:T) {
  X[,t]=PunifLE(h,p=p.val)
}

# As we added noise the observed list can be any permutation of 1:M. 
# We will make a histogram showing the probability for each possible 
# permutation to be the output of a simulation of the queue-jumping 
# noise process. For that we need to enumerate all permutations of [n] 
# and turn them into strings for plotting.

# list all perms of 1:M

perm=t(enum.seq(1:M)) 

# turn them into strings

perm.str=apply(perm,2,function(x){paste(x,collapse = ' ')}) 

# Which of the permutations correspond to LEs of h? There will be 
# a higher probability of simulating these when p is close to zero 
# as we get pure LEs when p=0.

# Color permutations that are LEs in blue, the rest in grey.

le=le.expand(h)
is.le=1+apply(perm,2,function(x){any(apply(le,2,function(y){all(x==y)}))})
cols=c('lightgrey','lightblue')[is.le]

# We have a list of all the permutations in perm. We have a list 
# of sampled permutations in X. Go through all the samples in X 
# - for each one identify it as one of the permutations in perm.

ord=apply(X,2,function(x){
    which(apply(perm,2,function(y){all(x==y)}))
  })
ft=table(ord)/T
#if p.val is small this may create an error as some orders will be missing
#either use a bigger T values or take a bigger p.val (or fix it properly yourself!)

barplot(ft~perm.str,las=2,cex.names=0.7,col=cols,
        ylab='probability',xlab='permutation',
        main=paste('distribution of random orders with noise probability = ',p.val),
        cex.main=0.8)
legend('topright',legend=c('LE of h','not LE'),
       col=c('lightblue','lightgrey'),pch=c(15,15))

# if you make the value of p.val small the distribution will
# concentrate on the LEs (the blue bars), if you make p.val
# close to one the distribution will converge to uniform and
# there will be no signal from the PO in the data
