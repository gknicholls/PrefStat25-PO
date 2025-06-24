dag.conc2<-function(mP,display.threshold=0.5,label.threshold=0.5,fontsize=24,show.closure=TRUE) {
  
  #compute reduction of concensus DAG from the n x n x T array
  #of adjacncy matrices X starting from mP, the T-wise average of X
  #mP<-meanPO(X)
  #include edges with support above display.threshold
  #label edges with support above label.threshold
  
  #n<-dim(X)[1]
  #who=colnames(X[,,1])
  #dag.con<-apply(X,3,my.transitive.closure) #ed 5-4-22 for mnem 
  ##each sample matrix is turned into a vector - we hope columnwise
  #dag.con<-matrix(apply(dag.con,1,mean),n,n,dimnames=list(who,who))
  #now written back into matrix columnwise - yuck, but seems OK 
  #dag.con<-meanPO(X)
  if (show.closure) {
    dag.con2<-my.transitive.closure(mP>display.threshold)
  } else {
    dag.con2<-transitive.reduction(mP>display.threshold)
  }
  dag.con<-round(mP,2)*dag.con2
  
  dag.con<-new("graphAM", adjMat=dag.con, edgemode="directed", values=list(weight=1))
  dag.con<-as(dag.con,"graphNEL") #graph_from_adjacency_matrix(dag.con,mode ="directed")#
  
  ewv<-unlist(edgeWeights(dag.con))
  #ew <- as.character(ewv); ew[ewv<label.threshold]<-""
  #ew <- ew[setdiff(seq(along = ew), removedEdges(dag.con))]
  #names(ew) <- edgeNames(dag.con)
  eAttrs <- list()
  #eAttrs$label <- ew
  #ft<-as.numeric(ew)*0+fontsize; names(ft)<-names(ew)
  #eAttrs$fontsize <- ft
  n.edges=length(edgeNames(dag.con))
  #e.style <- rep("dashed",n.edges); names(e.style)<-edgeNames(dag.con)
  #eAttrs$style<-e.style
  e.color <- rep("red",n.edges); names(e.color)<-edgeNames(dag.con)
  e.color[ewv<0.9]<-"black"
  eAttrs$color<-e.color
  
  #labeldistance <- rep(10,n.edges); names(labeldistance)<-edgeNames(dag.con)
  #eAttrs$labeldistance<-labeldistance
  
  #labelfloat<-rep("true",n.edges); names(labelfloat)<-edgeNames(dag.con)
  #eAttrs$labelfloat<-labelfloat
  
  nAttrs<-list()
  n.nodes=length(nodes(dag.con))
  nAttrs$height <- nAttrs$width <- rep("6", n.nodes)
  nAttrs$style <- rep("invis",n.nodes)
  nAttrs$fontsize <- rep(fontsize, n.nodes)
  nAttrs <- lapply(nAttrs, function(x) {names(x) <- nodes(dag.con); x })
  
  attrs <- list(node = list(shape = "circle", fixedsize = "false"), edge=list(headclip="true",
                                                      tailclip="true",style="dashed",arrowsize=1), graph=list(splines="true"))
  #plot(dag.con,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs,"dot") #how to plot
  
  list(dag.con=dag.con,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs)
}

meanPO<-function(X) {
  
  #compute "average" PO h over mcmc samples from the n x n x T array
  #of adjacncy matrices X 
  
  n<-dim(X)[1]
  who=colnames(X[,,1])
  a<-apply(X,3,my.transitive.closure) #ed 5-4-22 for mnem 
  #each sample matrix is turned into a vector - we hope column-wise
  H.hat<-matrix(apply(a,1,mean),n,n,dimnames=list(who,who))
  #and written back into matrix columnwise - yuck, but seems OK 
  
  return(H.hat)
}

dag2ig<-function(h,fontsize=8) {
  
  dag.ig<-new("graphAM", adjMat=h, edgemode="directed", values=list(weight=1))
  dag.ig<-as(dag.ig,"graphNEL") #graph_from_adjacency_matrix(dag.con,mode ="directed")#
  
  eAttrs <- list()
  n.edges=length(edgeNames(dag.ig))
  e.color <- rep("red",n.edges); names(e.color)<-edgeNames(dag.ig)
  eAttrs$color<-e.color
  
  nAttrs<-list()
  n.nodes=length(nodes(dag.ig))
  nAttrs$height <- nAttrs$width <- rep("20", n.nodes)
  nAttrs$style <- rep("invis",n.nodes)
  nAttrs$fontsize <- rep(fontsize, n.nodes)
  nAttrs <- lapply(nAttrs, function(x) {names(x) <- nodes(dag.ig); x })
  
  attrs <- list(node = list(shape = "plaintext", fixedsize = "false"), 
                edge=list(headclip="false",tailclip="true",style="dashed",arrowsize=2.0), 
                graph=list(splines="false"))
  #plot(dag.ig,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs,"dot") #how to plot
  
  list(dag.con=dag.ig,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs)
}


showDAGcon<-function(B,T,PO.con,years=1:T) {
  #display the concensus DAGs in PO.con arranged in a lattice and labeled by year
  n.years=length(years)
  n.row=min(4,floor(sqrt(n.years)))
  par(mfrow=c(n.row,ceiling(length(years)/n.row)),mai=c(0,0,0,0),xpd=NA)
  for (t in years) {
    h=PO.con[[t]]
    numEdges=length(unlist(edgeL(PO.con[[t]]$dag.con)))
    if (numEdges==0) {
      a=plot(h$dag.con,attrs=h$attrs,"dot")
    } else {
      a=plot(h$dag.con,nodeAttrs=h$nodeAttrs,edgeAttrs=h$edgeAttrs,attrs=h$attrs,"dot")
    }
    text(a@boundBox@upRight@x/2,1.04*a@boundBox@upRight@y,paste('[',as.character(B+t-1),']',sep=''),col=2,cex=1.5)
  }
}


outputanalysis<-function(out.dir,out.file,burn,doi.plot=NA,yoi=1:T,pdf.file=NA,P.samples=NA,full.analysis=TRUE,db=NA,dr=NA,U.con=FALSE) {
  
  #out.dir=note$RUNDIR;burn=2;pdf.file=NA; doi.plot=doi; yoi=1:T #yoi=c(1,2,3,23,24,25,35,36,37,38,39,40,41,52,53,54)
  #out.file="dbtpHB1a.RData"; full.analysis=TRUE; P.samples=NA
  #out.file="9506pcb.RData"
  #out.file='eightdaybigtrialrun.RData'
  
  cdir=getwd()
  setwd(out.dir)
  
  ################################################################################################
  #the file of mcmc output we are loading and the burnin
  
  loadfile=out.file     #this will be note$savefile for the run we want
  
  ################################################################################################
  #get mcmc output and cut any emply cells (if run hasnt finished)
  
  output=my.load(loadfile)
  
  setwd(cdir) #changed dirs at start turns out this has side effect in calling environment
  
  note=output$note; B=output$B;E=output$E;T=output$T;cil=output$cil;cla=output$cla;
  NB=output$NB;NF=output$NF;NL=output$NL;DB=output$DB;rank=output$rank;active=output$active;
  J=output$J;start=output$start;SS=output$SS;PO=output$PO;Uout=output$Uout;P=output$P;
  step.count=output$step.count; doi=output$doi
  
  if (exists("DOTAU",where=note)) {wide.lists=output$wide.lists; Nwl=length(wide.lists)}
  if (exists("accept",where=output)) {accept=output$accept;propose=output$propose}
  if (exists("st",where=output)) st=output$st
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) Sstate=output$Sstate
  if (is.na(doi.plot[1])) doi.plot=unique(sapply(cil,function(x) {x$diocese}))
  #attach(output)
  
  J.done=step.count #last complete step	
  (N.sample.done=J.done/SS+1) #starts at 0
  PO=PO[1:N.sample.done]
  Uout=Uout[1:N.sample.done]
  
  P=P[1:N.sample.done,] #cut any remaining NA's
  
  ############################################################
  
  #if (!is.na(pdf.file)) pdf(file=pdf.file,paper='a4')
  
  ############################################################
  #plot obvious stuff
  
  if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-mcmc.pdf',sep=''),6,4)}
  #print(colnames(P)) #to see what we have available to plot
  if (exists("DOQ",where=note) && note$DOQ) {par(mfrow=c(3,3),mai=0.1*c(8,8,1,1))} else {par(mfrow=c(2,3),mai=0.1*c(8,8,1,1))}
  plot(burn:N.sample.done,P[burn:N.sample.done,"oll"],xlab='MCMC-sample',ylab='log-lkd',type='l'); points(N.sample.done,P[N.sample.done,"oll"],pch=16,col=2)
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$oll.tot,col=2)
  plot(burn:N.sample.done,P[burn:N.sample.done,"lpU"],xlab='MCMC-sample',ylab='log U-prior',type='l'); points(N.sample.done,P[N.sample.done,"lpU"],pch=16,col=2)
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$oup$lpU.tot,col=2)
  plot(burn:N.sample.done,P[burn:N.sample.done,"obp"],xlab='MCMC-sample',ylab='log beta-prior',type='l'); points(N.sample.done,P[N.sample.done,"obp"],pch=16,col=2)
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$obp,col=2)
  plot(burn:N.sample.done,P[burn:N.sample.done,"rho"],xlab='MCMC-sample',ylab='rho',type='l'); points(N.sample.done,P[N.sample.done,"rho"],pch=16,col=2)
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$rho,col=2)
  plot(burn:N.sample.done,P[burn:N.sample.done,"theta"],xlab='MCMC-sample',ylab='theta',type='l'); points(N.sample.done,P[N.sample.done,"theta"],pch=16,col=2)
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$theta,col=2)
  if (exists("DOP",where=note)) {
    plot(burn:N.sample.done,P[burn:N.sample.done,"p"],xlab='MCMC-sample',ylab='p',type='l'); points(N.sample.done,P[N.sample.done,"p"],pch=16,col=2)
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$p,col=2)
  }
  if (!is.na(pdf.file)) dev.off()
  
  if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-rho-th-p.pdf',sep=''),6,2)}
  if (exists("DOQ",where=note) && note$DOQ) {par(mfrow=c(2,3),mai=0.1*c(8,8,1,1))} else {par(mfrow=c(1,3),mai=0.1*c(8,8,1,1))}
  plot(density(P[burn:N.sample.done,"rho"],from=0,to=1),xlab='rho',main='',ann=TRUE);
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(v=Sstate$rho,col=2)
  par(mai=0.1*c(8,4,1,1))
  plot(density(P[burn:N.sample.done,"theta"],from=0,to=1),xlab='theta',ylab='',main='',ann=TRUE);
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(v=Sstate$theta,col=2)
  if (exists("DOP",where=note)) {
    par(mai=0.1*c(8,4,1,1))
    plot(density(P[burn:N.sample.done,"p"],from=0,to=1),xlab='p',ylab='',main='',ann=TRUE); 
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(v=Sstate$p,col=2)
  }
  if (exists("DOQ",where=note) && note$DOQ) {
    plot(burn:N.sample.done,P[burn:N.sample.done,"q"]); points(N.sample.done,P[N.sample.done,"q"],pch=16,col=2)
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(h=Sstate$q,col=2)
    
  }
  if (exists("DOQ",where=note) && note$DOQ) {
    plot(density(P[burn:N.sample.done,"q"],from=0,to=1),ann=FALSE); 
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) abline(v=Sstate$q,col=2)
    plot(P[burn:N.sample.done,"p"],P[burn:N.sample.done,"q"]); points(N.sample.done,P[N.sample.done,"q"],pch=16,col=2)
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) points(Sstate$p,Sstate$q,col=2,pch=16,cex=3)
  }
  if (!is.na(pdf.file)) dev.off()
  
  beta.ind=grep("beta",colnames(P))  
  uncentred.beta=P[,beta.ind]   	#need to compute Z below
  M=matrix(-1/DB,DB,DB)+diag(DB)
  centred.beta=uncentred.beta%*%M
  
  if (TRUE) {
    Zout=Uout; 
    #centred.beta=uncentred.beta
    Hout=PO
    for (j in 1:N.sample.done) {
      um=mean(Uout[[j]],na.rm=TRUE)
      #us=sd(Uout[[j]],na.rm=TRUE)
      Uout[[j]]=(Uout[[j]]-um)#/us
      #centred.beta[j,]=(uncentred.beta[j,]+um)/us
      Zout[[j]]=ComputeZ(years=1:T,bishops=1:NB,Uout[[j]],rank,centred.beta[j,])
      #Hout[[j]]=ZtoPO(Zout[[j]],active,B,T)
    }
    #if (!identical(Hout,PO)) stop('centering messed up') 
  } else {
    #M=matrix(-1/DB,DB,DB)+diag(DB)
    #centred.beta=uncentred.beta%*%M
    for (j in 1:N.sample.done) Zout[[j]]=ComputeZ(years=1:T,bishops=1:NB,Uout[[j]],rank,centred.beta[j,])
  }
  
  #P[,beta.ind]=P[,beta.ind]%*%M #centering beta's for sensible plotting - but still work to centre beta, U jointly
  
  
  if (is.na(pdf.file)) {windows()} else {pdf(file=paste(pdf.file,'-beta.pdf',sep=''),paper='a4',4,6)} 
  boxplot(centred.beta[burn:N.sample.done,],xlab='seniority rank', ylab='effect, beta',cex.axis=0.8,boxwex=1)
  
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) points(1:length(beta.ind),Sstate$beta,col=2,pch=16)
  if (!is.na(pdf.file)) dev.off()
  
  tau.ind=grep("tau",colnames(P)) 
  if (exists("DOTAU",where=note) && note$DOTAU) {
    tau.pm=apply(P[burn:N.sample.done,tau.ind],2,mean)
    #############################################################################
    #make a nice graph of the HPD sets for tau
    
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-tau.pdf',sep=''),8,6)}
    P.tau=P[burn:N.sample.done,tau.ind[wide.lists]]
    mt=apply(P.tau,2,median); omt=order(mt)
    boxplot(B-1+P.tau[,omt],border="0",col=0,xaxt='n',xlab='',ylab='year range')
    title(xlab="list index", mgp=c(1,1,0))
    legend("topleft",legend=c('data constraint','HPD set for list time'),lty=c(1,NA),col=c("red","grey"),pch=c(NA,15))
    tl=sapply(output$cla,function(x) x$tl)[wide.lists[omt]]
    tu=sapply(output$cla,function(x) x$tu)[wide.lists[omt]]
    hpd=lapply(1:length(omt),function(i) {hpd.set(P.tau[,omt[i]],tl[i],tu[i],alpha=0.5,B=B)})
    dev.null<-lapply(1:length(omt),function(x){points(rep(x,length(hpd[[x]])),hpd[[x]],pch=15,cex=0.8,col="grey")})
    dev.null<-sapply(1:length(omt),function(x) {li=wide.lists[omt[x]]; lines(c(x,x),c(tl[x],tu[x]),col=2)})
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) {
      points(B-1+Sstate$tau[wide.lists[omt]],pch=16,col=1); 
    }
    if (!is.na(pdf.file)) dev.off()
  } else {
    tau.pm=NA
  }
  
  ############################################################
  #Can just do a light touch report with basic summary plot
  
  if (full.analysis) {
    
    ############################################################
    #reform output - entries in PO and Uout give all the years for one MCMC sample
    #reprocess Uout and PO into lists with each entry giving all the MCMC samples for one year 
    #(and compute same sort of thing for Z)
    
    #set concensus graph thresholds for black and red edges
    if (is.na(db)) db=0.5
    if (is.na(dr)) dr=0.9
    
    years=1:T
    N.years=length(years)
    N.po=length(PO)
    Z=U=PO.con=POU.con=PO.con.tr=Y=W=POU.hat=PO.hat=vector('list',N.years)
    Z.mean=Z.sd=U.sd=U.mean=matrix(NA,NB,T)
    POU=PO
    if (U.con) {
      for (n in burn:N.po) {
        POU[[n]]=ZtoPO(Uout[[n]],active,B,T)
      }
    }
    for (t in years) {
      n.active.t=sum(active[,t]) #number active in year t
      who.active=colnames(PO[[1]][[t]])
      if (length(who.active)!=n.active.t) stop('number active doesnt match who active')
      Y[[t]]=array(NA,c(n.active.t,n.active.t,N.po-burn+1),dimnames=list(who.active,who.active,NULL))
      if (U.con) W[[t]]=Y[[t]]
      Z[[t]]=U[[t]]=array(NA,c(NB,NF,N.po-burn+1),dimnames=list(1:NB,1:NF,NULL))
      count=0
      for (n in burn:N.po) {
        if (!all(dim(PO[[n]][[t]])==c(n.active.t,n.active.t))) stop('dimensions dont make sense')
        count=count+1
        Y[[t]][,,count]=PO[[n]][[t]]
        if (U.con) W[[t]][,,count]=POU[[n]][[t]]
        U[[t]][,,count]=Uout[[n]][,,t]
        Z[[t]][,,count]=Zout[[n]][,,t]  #ComputeZ(years=t,bishops=1:NB,Uout[[n]],rank,uncentred.beta[n,])
      }
      PO.hat[[t]]=meanPO(Y[[t]])   #real matrices with prob for order as entries
      if (U.con) POU.hat[[t]]=meanPO(W[[t]])
      PO.con[[t]]=dag.conc2(PO.hat[[t]],display.threshold=db,label.threshold=dr,fontsize=8,show.closure=TRUE) #a list with graph properties
      PO.con.tr[[t]]=dag.conc2(PO.hat[[t]],display.threshold=db,label.threshold=dr,fontsize=8,show.closure=FALSE) #a list with graph properties
      if (U.con) POU.con[[t]]=dag.conc2(POU.hat[[t]],display.threshold=db,label.threshold=dr,fontsize=8,show.closure=TRUE) #a list with graph properties
      
      U.mean[,t]=apply(U[[t]],1,mean,na.rm=TRUE)
      U.sd[,t]=apply(U[[t]],1,sd,na.rm=TRUE)
      Z.mean[,t]=apply(Z[[t]],1,mean,na.rm=TRUE)
      Z.sd[,t]=apply(Z[[t]],1,sd,na.rm=TRUE)
    }
    
    ############################################################
    #produce plot of U/Z values by diocese
    
    #doi #to see which dioceses we want
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) {
      sU=apply(Sstate$U,c(1,3),mean,na.rm=TRUE)
      sZ=apply(Sstate$Z,c(1,3),mean,na.rm=TRUE)
    } else {
      sU=sZ=NA
    }
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-Uav.pdf',sep=''),7.8,10)}
    showDOIs(doi.plot=doi.plot,V.mean=U.mean,V.sd=U.sd,cil.plot=cil,B,T,sV=sU,tau.pm.plot=tau.pm)
    if (!is.na(pdf.file)) dev.off()
    
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-Zav.pdf',sep=''),7.8,10)}
    showDOIs(doi.plot=doi.plot,V.mean=Z.mean,V.sd=Z.sd,cil.plot=cil,B,T,sV=sZ,tau.pm.plot=tau.pm)
    if (!is.na(pdf.file)) dev.off()
    
    #plot(1:T,rep(0,T),ylim=c(-2,2),type='n')
    #count=1; apply(U.mean,1,function(x){lines(1:T,x,col=count);count=<<-count+1})
    
    ############################################################
    #look at probs to come first or last as alternative to plotting U, Z etc
    
    if (!is.na(P.samples)) {
      subsub.sample=floor(seq(burn,N.sample.done,length.out=P.samples))
      if (FALSE){ #(note$DOPAR) {
        cl<-makeCluster(5) #note$nthread)
        clusterEvalQ(cl,{library(lecount); source(file="pofun.R")})
        P.tb<-clusterApplyLB(cl, PO[subsub.sample], p2, NB)
        stopCluster(cl); gc()
      } else {
        P.tb<-lapply(PO[subsub.sample], p2, NB)
      }
      P.diff.sd=P.bot.sd=P.top.sd=P.bot.hat=P.top.hat=matrix(0,NB,T)
      for (j in 1:P.samples) {
        P.top.hat=P.top.hat+P.tb[[j]]$P.top
        P.bot.hat=P.bot.hat+P.tb[[j]]$P.bot
      }
      P.top.hat=P.top.hat/P.samples
      P.bot.hat=P.bot.hat/P.samples
      P.diff.hat<-P.top.hat-P.bot.hat
      for (j in 1:P.samples) {
        P.top.sd=P.top.sd+(P.tb[[j]]$P.top-P.top.hat)^2
        P.bot.sd=P.bot.sd+(P.tb[[j]]$P.bot-P.bot.hat)^2
        P.diff.sd=P.diff.sd+(P.tb[[j]]$P.top-P.tb[[j]]$P.bot-P.diff.hat)^2
      }
      P.top.sd=sqrt(P.top.sd/(P.samples-1))
      P.bot.sd=sqrt(P.bot.sd/(P.samples-1))
      P.diff.sd=sqrt(P.diff.sd/(P.samples-1))
      
      P.top.hat[P.top.hat==0]<-NA
      P.bot.hat[P.bot.hat==0]<-NA
      P.diff.hat[P.diff.hat==0]<-NA
      P.top.sd[P.top.hat==0]<-NA
      P.bot.sd[P.bot.hat==0]<-NA
      P.diff.sd[P.diff.hat==0]<-NA
      
      #if (!is.na(pdf.file)) dev.off()
      #showDOIs(doi.plot=doi.plot,V.mean=P.top.hat,V.sd=P.top.sd,cil.plot=cil,B,T,sV=NA,plot.U=FALSE,tau.pm.plot=tau.pm)
      #if (is.na(pdf.file)) windows() 
      #showDOIs(doi.plot=doi.plot,V.mean=P.bot.hat,V.sd=P.bot.sd,cil.plot=cil,B,T,sV=NA,plot.U=FALSE,tau.pm.plot=tau.pm)
      if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-pFLdiff.pdf',sep=''),7.8,10)}
      showDOIs(doi.plot=doi.plot,V.mean=P.diff.hat,V.sd=P.diff.sd,cil.plot=cil,B,T,sV=NA,plot.U=FALSE,tau.pm.plot=tau.pm)
      if (!is.na(pdf.file)) dev.off()
    }
    
    ############################################################
    #Consensus PO 
    #TODO - function below needs work - not very pretty
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-POconB.pdf',sep=''),40,20)} #8,10
    showDAGcon(B,T,PO.con,years=yoi)
    if (!is.na(pdf.file)) dev.off()
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-POconB-tr.pdf',sep=''),40,20)} #8,10
    showDAGcon(B,T,PO.con.tr,years=yoi)
    if (!is.na(pdf.file)) dev.off()
    
    if (FALSE) {
      get.sgi.vsp<-function(tc) {
        tcg<-graph_from_adjacency_matrix(tc,mode ="directed")
        b<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),4,4,byrow=TRUE)
        bg<-graph_from_adjacency_matrix(b,mode ="directed")
        return(subgraph_isomorphisms(bg,tcg,induced=TRUE))
      }
      dag.con3<-lapply(PO.hat,function(x) my.transitive.closure(x>dr))
      sgi3=lapply(dag.con3,get.sgi.vsp)
      dag.con4<-lapply(PO.hat,function(x) my.transitive.closure(x>db))
      sgi4=lapply(dag.con4,get.sgi.vsp)
      
      for (t in 1:length(dag.con3)){
        if (length(sgi3[[t]])>0) {
          for (i in 1:length(sgi3[[t]])) {
            good=FALSE
            if (length(sgi4[[t]])>0) {
              for (j in 1:length(sgi4[[t]])) {
                if (all(sort(sgi3[[t]][[i]])==sort(sgi4[[t]][[j]]))) {
                  good=TRUE
                }
              }
            }
            if (good) print(sprintf('(%d,%d)',t,i))
          }
        }
      }
      t=58; i=103
      a=sgi3[[t]][[i]]
      par(mfrow=c(1,2))
      showDAG(dag.con3[[t]][a,a])
      showDAG(dag.con4[[t]][a,a])
      showDAGcon(B,T,PO.con,years=t)
      round(PO.hat[[t]][a,a],2)
      
      #>     round(PO.hat[[t]][a,a],2)
      #35   48 47   38
      #35 0.00 0.89  0 0.44
      #48 0.00 0.00  0 0.00
      #47 0.29 0.99  0 0.92
      #38 0.00 0.46  0 0.00
      
    }#end forbidden subgraph search
    
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-POmcmc.pdf',sep=''),6,3)}
    showDAGs(B,T,PO[[N.sample.done]],years=yoi)
    if (!is.na(pdf.file)) dev.off() 
    if (U.con) showDAGcon(B,T,POU.con,years=1:T)
    if (exists("DOSYNTH",where=note) && note$DOSYNTH) {
      POt=lapply(Sstate$h,function(x) {y<-transitive.reduction(x); dag2ig(y,fontsize=20)}) 
      if (is.na(pdf.file)) windows(); 
      showDAGcon(B,T,POt,years=yoi)
      Ht=Sstate$h
      #scoring=vector('list',T)
      #for (t in 1:T) scoring[[t]]<-sc(t)
      Scoring=lapply(1:T, 
                     function(t) { 
                       Htz=Ht[[t]]; diag(Htz)<-1
                       an=sum( (Htz==0 & t(Htz)==0 ) ) #num pairs of nodes with no reln
                       ap=sum( (Ht[[t]]==1 ) )                   #num of nodes on the > side of a reln
                       fn=sum( (Ht[[t]]==1) & ( PO.hat[[t]]<0.5 & t(PO.hat[[t]])<0.5 ) ) #reln in true not in est
                       fp=sum( (Htz==0 & t(Htz)==0) & ( PO.hat[[t]]>=0.5 | t(PO.hat[[t]])>=0.5 ) ) #no reln in true but yes in est
                       tp=sum( (Ht[[t]]==1) & (PO.hat[[t]]>=0.5 )) #reln in both - right way round
                       tn=sum( (Htz==0 & t(Htz)==0) & ( PO.hat[[t]]<0.5 & t(PO.hat[[t]])<0.5 ) ) #no reln in both
                       wr=sum( (Ht[[t]]==1) & (t(PO.hat[[t]])>=0.5) ) #reln in both but wrong way round
                       n.ht=dim(Ht[[t]])[1]
                       if ( (an+2*ap)!=(n.ht*(n.ht-1))) stop('ap+an issues')
                       if ( (wr+fn+tp)!= ap) stop('wwr+fn+tp issues')
                       if ( (fp+tn)!= an) stop('fp+tn issues')
                       #ap<-max(ap,1); an<-max(an,1)
                       return(list(fn=fn,fp=fp,tp=tp,tn=tn,wr=wr,ap=ap,an=an))
                     })
      if (is.na(pdf.file)) {windows()} else {pdf(file=paste(pdf.file,'-false-pos-neg.pdf',sep=''),6,6)} 
      dates=years+B-1
      Fn=sapply(Scoring,function(x) x$fn)
      Fp=sapply(Scoring,function(x) x$fp)
      Tp=sapply(Scoring,function(x) x$tp)
      Tn=sapply(Scoring,function(x) x$tn)
      Wr=sapply(Scoring,function(x) x$wr)
      y.max=max(c(Fn,Fp,Tp,Tn,Wr))
      #y.min=min(c(Fn,Fp,Tp,Tn,Wr))
      plot(dates,Fn,col=1,type='l',ylim=c(0,y.max),xlab='calendar year',ylab='positive and negative counts')
      lines(dates,Fp,col=2)
      lines(dates,Tp,col=3)
      lines(dates,Tn,col=4)
      lines(dates,Wr,col=5)
      legend('topleft',legend = c("false -ve order","false +ve order","true +ve order","true -ve order","order reversal"),lty=c(1,1,1,1,1),col=c(1,2,3,4,5))
      if (!is.na(pdf.file)) dev.off() 
    }
    
    b.names=sapply(cil,function(x) x$name); nrb=ceiling(NB/2)
    
    if (is.na(pdf.file)) {windows(10,10)} else {pdf(file=paste(pdf.file,'-names.pdf',sep=''),paper='a4')}
    par(mfrow=c(1,2),mai=0.1*c(1,1,1,1)); 
    plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',xaxt='n',yaxt='n',ann=FALSE) 
    for (k in 1:nrb) text(0,1-((k-1)/nrb),sprintf("%g) %s\n",k,b.names[k]),adj=c(0,1),cex=0.6)
    plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',xaxt='n',yaxt='n',ann=FALSE)  
    for (k in (nrb+1):NB) text(0,1-((k-nrb-1)/nrb),sprintf("%g) %s\n",k,b.names[k]),adj=c(0,1),cex=0.6)
    if (!is.na(pdf.file)) dev.off() 
    
    #################################################################################
    #traces for beta and U - also depth
    
    if (FALSE) {
      if (is.na(pdf.file)) windows() #beta traces are a bit of a mess
      plot(burn:N.sample.done,0*(burn:N.sample.done),type='n',ylim=c(-2,2),ylab="beta")
      count=1
      apply(centred.beta[burn:N.sample.done,],2,function(x){lines(burn:N.sample.done,x,col=count);count<<-count+1})
      
      #We have a problem with mixing of U in 1101 (year 1) and 1136 (year 36 in main runs)
      n.yoi=length(yoi)
      if (is.na(pdf.file)) windows(20,10) #U traces?
      par(mfrow=c(ceiling(n.yoi/5),min(n.yoi,5)),mai=c(0.2,0,0.2,0),oma=c(1,2,1,2))
      for (t in yoi) {
        plot(burn:N.sample.done,0*(burn:N.sample.done),type='n',ylim=c(-4,4),main=as.character(B+t-1))
        count=1
        #X=apply(Z[[t]],c(1,3),mean,na.rm=TRUE) #for Z #[,,burn:N.sample.done]
        X=apply(U[[t]],c(1,3),mean,na.rm=TRUE) #[,,burn:N.sample.done]
        apply(X,1,function(x){lines(burn:N.sample.done,x,col=count);count<<-count+1})
      }
      
      ############################################################
      
      #how does the depth vary with time
      PO.depth=t(sapply(PO[burn:N.sample.done],function(x) sapply(x,dagdepth)))
      PO.depth.mean=apply(PO.depth,2,mean)
      PO.depth.sd=apply(PO.depth,2,sd)
      n.active=apply(active,2,sum)
      
      if (is.na(pdf.file)) windows() 
      par(mfrow=c(2,2))
      plot(B-1+1:T,rep(NA,T),ylim=c(0,max(PO.depth)+1))
      #apply(PO.depth,1,function(x) lines(B-1+1:T,x))
      lines(B-1+1:T,PO.depth.mean,lwd=1,col=2); 
      lines(B-1+1:T,PO.depth.mean-PO.depth.sd,lty=2,lwd=1,col=2)
      lines(B-1+1:T,PO.depth.mean+PO.depth.sd,lty=2,lwd=1,col=2)
      if (exists("DOSYNTH",where=note) && note$DOSYNTH) {lines(B-1+1:T,sapply(Sstate$h,dagdepth),col=3,lwd=2)}
      
      plot(B-1+1:T,n.active,type='l')
      
      plot(B-1+1:T,rep(NA,T),ylim=c(0,1.1*max(PO.depth/n.active)))
      lines(B-1+1:T,PO.depth.mean/n.active,lwd=1,col=3); 
      lines(B-1+1:T,(PO.depth.mean-PO.depth.sd)/n.active,lty=2,lwd=1,col=3)
      lines(B-1+1:T,(PO.depth.mean+PO.depth.sd)/n.active,lty=2,lwd=1,col=3)
      abline(v=1135) #hmm... 1135-1154 King Stephen 
      
      plot(n.active,PO.depth.mean,type='n'); text(n.active,PO.depth.mean,as.character(B-1+1:T))
      
    } #skip some stuff
    
    ############################################################
    
  } # closing "if (full.analysis) {"
  
  ############################################################
  
  #if (!is.na(pdf.file)) dev.off()
  
  ############################################################
  
  #effective sample sizes
  ESS=effectiveSize(P[burn:N.sample.done,-tau.ind])
  
  #acceptance rates
  if (exists('accept') && exists('propose')) {
    acr=round(100*accept/propose)
  } else {acr=NA}
  
  if (exists('st')) {
    st.spu=sum(st)/step.count; 
    st.prop=round(100*st/sum(st))
  } else {
    st.spu=NA; st.prop=NA
  }
  
  #run setup summary
  if (exists('split.years',where=note)) {note$split.years<-NULL}
  if (exists('doi',where=note)) {cnote=note[-grep('doi',names(note))]} else {cnote=note}
  run.summary.all=as.data.frame(cnote); 
  run.summary.all$Usweeps2=run.summary.all$Usweeps[c(2,1)]
  fields=c("RANDSEED","B","E","NF","constrainbeta","model","MCMC.SWEEPS","SUBSAMPLE.INTERVAL","Usweeps","Usweeps2","UALLFRACUPD","STARTSTATE","init.p","init.theta","DOTHETA","DOSYNTH")
  show.entries=is.element(colnames(run.summary.all),fields)
  run.summary=run.summary.all[1,show.entries]
  tn=matrix(c(NF,NL,NB,DB),1,4)
  colnames(tn)<-c("NF","NL","NB","DB")
  if (exists("DOSYNTH",where=note) && note$DOSYNTH) {
    stn=matrix(c(output$NF,output$NL,output$NB,output$DB),1,4); colnames(stn)<-c("sNF","sNL","sNB","sDB")
    tn=cbind(tn,stn)
  }
  run.used=list(model.dims=tn)
  if (exists("PhPar",where=note)) {run.used=c(run.used,list(PhPar=unlist(note$PhPar)))}
  if (exists("DOQ",where=note)) {run.used=c(run.used,list(DOQ=note$DOQ))}
  
  print(list(ESS=ESS,acceptance.rate=acr,update.percent.time=st.prop,seconds.per.mcmc.iteration=st.spu,run.summary=run.summary,run.used))
  
  ############################################################
  if (full.analysis) {
    return(list(ESS=ESS,acceptance.rate=acr,update.percent.time=st.prop,seconds.per.mcmc.iteration=st.spu,run.summary=run.summary,run.used=run.used,PO.con=PO.con))
  } else {
    return(list(ESS=ESS,acceptance.rate=acr,update.percent.time=st.prop,seconds.per.mcmc.iteration=st.spu,run.summary=run.summary,run.used))
  }
}

p2<-function(x,NB) {
  #compute prob for each bishop in each year to be top or bottom
  X<-lapply(x,p1,'first')
  Y<-lapply(x, p1, 'last')
  P.top=sapply(X, function(x) {i=as.numeric(names(x)); a=rep(0,NB); a[i]<-x; return(a)})
  P.bot=sapply(Y, function(x) {i=as.numeric(names(x)); a=rep(0,NB); a[i]<-x; return(a)})
  return(list(P.top=P.top,P.bot=P.bot))
}

hpd.set<-function(x,l,u,alpha=0.66,B=1) {
  #x are integer samples l and u are upper and lower range (might exceed range in samples)
  #designed for tau - works for integers
  a=hist(B-1+x,breaks=(l-0.5):(u+0.5),plot=FALSE)$density
  a<-matrix(a,1,length(a)); colnames(a)<-l:u
  hpd=c(); tot=0
  while (tot<alpha) {i=which.max(a); tot<-tot+a[i]; hpd<-c(hpd,colnames(a)[i]); a<-a[1,-i,drop=FALSE]}
  return(sort(as.integer(hpd)))
}
