
library(deSolve)
library(igraph)

#Model function: Non-dimensional Rosenzweig MacArthur; Predator-Prey

Model_Equation = function(t,x,parms=NULL,patch,r,m.e,dh,dp,ei,d){
  
  time = t
  
  x[x<0]=0
  
  
  seq(from=1,to=patch*2,by=2)
  
  hi = matrix(x[seq(from=1,to=(patch*2-1),by=2)],nrow=patch,ncol=patch,byrow = TRUE)
  
  p.i = matrix(x[seq(from=2,to=patch*2,by=2)],nrow=patch,ncol=patch,byrow = TRUE)
  
  hi[hi<0]=0
  
  p.i[p.i<0]=0
  
  m.e = matrix(m.e, nrow = patch, ncol = patch)
  
  dH = hi[1,]*(1-r*hi[1,])-(p.i[1,]*hi[1,])/(1+hi[1,])+dh*rowSums(m.e*hi)
  
  dP = (ei*p.i[1,]*hi[1,])/(1+hi[1,])-d*p.i[1,]+dp*rowSums(m.e*p.i)
  
  
  res = c(rbind(dH,dP))
  list(res)
  
}

#Set equation parameters

    #predator conversion efficiency
    ei = 5
    
    #prey growth
    r = .3
    
    #predator mortality
    d = 1
    
    #prey dispersal
    dh = .018
    
    #predator dispersal
    dp = .018
    
    #number of patches
    patch = 10
    
    #calculate equilibrium values
    hhat = d/(ei-d)
    phat = (1+hhat)*(1-r*hhat)

#set simulation parameters

    #length of simulation
    sim.length = 5000
    
    #variability of initial conditions
    P = .05

#Generate random network
    
    #k.regular.game(#patches,#connections per patch) -> generate regular graphs
    NM = k.regular.game(patch,4)
    
    A = get.adjacency(NM,sparse=FALSE)
         
    dimnames(A)=NULL
    
    E = matrix(0,nrow = patch,ncol=patch)
    
    diag(E)=rowSums(A)
    
    m.e = A-E
    
    dimnames(m.e)=NULL
    
#Perturb starting densities
    
    hini = runif(patch,(hhat-P*hhat),(hhat+P*hhat))
    pini = runif(patch,(phat-P*phat),(phat+P*phat))

#Simulation

    t.ana = seq(from=1,to=sim.length,by=1)
  
    xstart = c(rbind(hini,pini))
  
    sim.ana = as.data.frame(lsoda(xstart,t.ana,parms = NULL, Model_Equation, rtol = 1e-6, atol = 1e-6,
                                jacfunc = NULL, jactype = "fullint", rootfunc = NULL,
                                verbose = FALSE, nroot = 0, tcrit = NULL,
                                hmin = 0, hmax = NULL, hini = 0, ynames = TRUE,
                                maxordn = 12, maxords = 5, bandup = NULL, banddown = NULL,
                                maxsteps = 10000, dllname = NULL, initfunc = NULL,
                                initpar = parms, rpar = NULL, ipar = NULL, nout = 0,
                                outnames = NULL, forcings = NULL, initforc = NULL,
                                fcontrol = NULL, events = NULL, lags = NULL,
                                patch,r,m.e,dh,dp,ei,d))[,-1]
    
    prey = sim.ana[,seq(1,ncol(sim.ana),by=2)]
    pred = sim.ana[,seq(2,ncol(sim.ana),by=2)]

#measure variability
    prey.var.total = sd(unlist(prey[(sim.length-1000:sim.length),]))/mean(unlist(prey[(sim.length-1000:sim.length),]))
    prey.var.patch = apply(prey[(sim.length-1000:sim.length),],2,sd)/apply(prey[(sim.length-1000:sim.length),],2,mean)
    
    pred.var.total = sd(unlist(pred[(sim.length-1000:sim.length),]))/mean(unlist(pred[(sim.length-1000:sim.length),]))
    pred.var.patch = apply(pred[(sim.length-1000:sim.length),],2,sd)/apply(pred[(sim.length-1000:sim.length),],2,mean)
    
#find time to equilibrium
    
    S.star = sim.ana[nrow(sim.ana),]
    
    s.find = 1
    s.test = 1
    
    phase.pos = NULL
    
    Last.period = as.numeric(seq(0,0,length.out=patch))
    
    N.period = as.numeric(seq(0,0,length.out=patch))
    
    Period.p = matrix(NA,nrow = patch,ncol=8000)
    
    Z.i = matrix(NA,nrow = patch,ncol=8000)
    
    p.trans = as.numeric(seq(0,0,length.out=patch))
    
    phase.pos = rbind(phase.pos,as.numeric(hhat*sim.ana[1,seq(from=2,to=(patch*2),by=2)]-phat*sim.ana[1,seq(from=1,to=(patch*2-1),by=2)]))
    phase.pos = rbind(phase.pos,as.numeric(hhat*sim.ana[2,seq(from=2,to=(patch*2),by=2)]-phat*sim.ana[2,seq(from=1,to=(patch*2-1),by=2)]))
    
    is.d = 3
    set = FALSE
    while(is.d<=nrow(sim.ana)){
      
      phase.pos = rbind(phase.pos,as.numeric(hhat*sim.ana[is.d,seq(from=2,to=(patch*2),by=2)]-phat*sim.ana[is.d,seq(from=1,to=(patch*2-1),by=2)]))
      
      ip =1
      while(ip<=patch){
        
        if(phase.pos[is.d,ip]>0&phase.pos[is.d-1,ip]<0){
          
          N.period[ip] = N.period[ip]+1
          
          Period.p[ip,N.period[ip]] = is.d-Last.period[ip]
          
          Last.period[ip] = is.d
          
        }
        
        ip=ip+1 
      }
      
      is.d = is.d+1
    }
    
    sinc = function(x,M=150,fc=.05,alpha=.167){
      
      h=NULL
      
      i.seq = seq(0,M)
      
      i=1
      while(i<=length(i.seq)){
        if(i.seq[i]==(M/2)){
          h=c(h,2*pi*fc)}else{
            h=c(h,((sin(2*pi*fc*(i.seq[i]-M/2)))/(i.seq[i]-(M/2)))*(.42-.5*cos((2*pi*i.seq[i])/M)+.08*cos((4*pi*i.seq[i])/M)))
          }
        i=i+1
      }
      
      h=h/sum(h)
      
      h.full = h
      
      h = NULL  
      y=NULL
      j=1
      while(j<=length(x)){
        
        
        h = h.full
        
        
        y=c(y,0)
        i=1
        while(i<=j&i<=M+1){
          
          y[j]=y[j]+x[j-i+1]*h[i]
          
          i=i+1
        }
        
        j=j+1 
      }
      
      y
    }
    
    ip=1
    while(ip<=patch){
      
      Obs.int = N.period[ip] - .167*N.period[ip]
      
      Z.raw = sinc(Period.p[ip,1:N.period[ip]])
      
      Z.alpha.mean = mean(Z.raw[Obs.int:N.period[ip]])
      
      Z.i[ip,1:N.period[ip]] = Z.raw - Z.alpha.mean
      
      ip=ip+1
    }
    
    "%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
    
    iz = 1
    while(iz<=patch){
      
      iz2 = 1
      while(iz2<=2000){
        
        if(FALSE%!in%(abs(na.omit(Z.i[iz,iz2:2000]))<.1)){
          p.trans[iz]=sum(na.omit(Period.p[iz,1:(iz2-1)]))
          break
        }
        
        iz2 = iz2+1
      }
      
      iz=iz+1
    }

    trans.length.patch = p.trans
    
    p.trans = sort(p.trans)
    
    trans.length.total = p.trans[8]
    
#find number of synchronized clusters    
    is.p = 1
    cluster.test = NULL
    while(is.p<=patch)
    {
      cluster.test = cbind(cluster.test,prey[(sim.length-50):sim.length,is.p])
      
      is.p=is.p+1
    }
    
    cluster.test = as.matrix(prey[(sim.length-50):sim.length,])
    
    cluster.cor = cor(cluster.test)
    
    cluster.id = seq(from=1,to=patch)
    patch.cluster = seq(from=1,to=patch)
    
    is.c.max = 1
    while(is.c.max<=length(cluster.id)){
      is.c = is.c.max+1
      n.cluster = is.c.max
      c.cache = cluster.id
      while(is.c<=patch){
        
        if(cluster.cor[c.cache[is.c.max],is.c]>.999){
          cluster.id = cluster.id[!cluster.id==is.c]
          patch.cluster[is.c]=c.cache[is.c.max]
        }
        
        is.c=is.c+1
      }
      is.c.max=is.c.max+1
    }
    
    unq.clust.ids = unique(patch.cluster)
    
    i.clust = 1
    while(i.clust<=length(unq.clust.ids))
    {
      
      patch.cluster[patch.cluster==unq.clust.ids[i.clust]] = i.clust
      
      i.clust=i.clust+1
    }
    
#Output
        
#set up plotting
    
    start.time = sim.length-1000
    
    cols = rainbow(ncol(prey))
    
    par(mfrow=c(2,2))
    
#plot network structure
    
    plot.igraph(NM,vertex.color=patch.cluster,vertex.size=40,edge.width=3)
    
#plot pred abundance vs prey abundance
    
    plot(prey[start.time:sim.length,1],pred[start.time:sim.length,1],pch=19,col=cols[1],xlim=c(min(prey[start.time:sim.length,]),max(prey[start.time:sim.length,])),ylim=c(min(pred[start.time:sim.length,]),max(pred[start.time:sim.length,])),ylab="Pred",xlab="Prey")
    
    i = 2
    while(i<=ncol(pred)){
      
      points(prey[start.time:sim.length,i],pred[start.time:sim.length,i],col=cols[i],pch=19)  
      
      i=i+1
    }
    
#plot prey abundance vs time
    
    plot(seq(start.time,sim.length),prey[start.time:sim.length,1],type="l",col=cols[1],ylim=c(min(prey[start.time:sim.length,]),max(prey[start.time:sim.length,])),ylab="Abundance",xlab="Time",main="Prey")
    
    i = 2
    while(i<=ncol(prey)){
      
      lines(seq(start.time,sim.length),prey[start.time:sim.length,i],col=cols[i])  
      
      i=i+1
    }
    
#plot pred abundance vs time
    
    plot(seq(start.time,sim.length),pred[start.time:sim.length,1],type="l",col=cols[1],ylim=c(min(pred[start.time:sim.length,]),max(pred[start.time:sim.length,])),ylab="Abundance",xlab="Time",main="Pred")
    
    i = 2
    while(i<=ncol(pred)){
      
      lines(seq(start.time,sim.length),pred[start.time:sim.length,i],col=cols[i])  
      
      i=i+1
    }
    
#Summary Data
    
    Total.Summary = c(prey.var.total,pred.var.total,trans.length.total,n.cluster)
    
    Patch.Summary = cbind(prey.var.patch,pred.var.patch,trans.length.patch,patch.cluster)