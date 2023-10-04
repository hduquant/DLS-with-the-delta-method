library(mvtnorm)
library(MASS)
library(matrixcalc)
library(doParallel)
library(CompQuadForm)
#library(Matrix)
library(expm)
library(combinat)
library(data.table)
#library(bigmemory)
count.dups <- function(DF){
  
  DT <- data.table(DF)
  DT[,.N, by = names(DT)]
}

RUNVARS <- commandArgs(T)
RUNMETHOD <- as.integer(RUNVARS[1])
id <- as.integer(RUNVARS[2]) #dim(grid)
ncores <- as.integer(RUNVARS[3])
FILEPATH <- RUNVARS[4]
#don't care (RUNVARS[5])

#----------------------------------------------------------------
# Computing transformed skewness and kurtosis, consistent with FML
#----------------------------------------------------------------
Mardia.tran1 <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x.c <- x-matrix(apply(x,2,mean),nrow=n,ncol=p,byrow=T)
  
  S <- cov(x)*(n-1)/n
  S_inv <- solve(S)
  D <- x.c%*%S_inv%*%t(x.c)
  
  b1 <- sum(D^3)
  b2 <- sum(diag(D^2))/n   #Mardia's kurtosis
  
  df <- p*(p+1)*(p+2)
  sk <- b1/(n*df)
  kt <- b2/(p*(p+2))
  
  return(list('skew'=sk,'kurt'=kt))
}

################

vech<-function(x) return(x[!upper.tri(x)]) ########################

################
# ----------------------------------- 
# normal data 
sim.norm <- function(seed,n,p,sig12){
  set.seed(seed)
  z <- rmvnorm(n=n, mean=rep(0,p))
  x <- z%*%sig12
  return(x)
}
# ----------------------------------- 
# elliptical data 
sim.enorm <- function(seed,n,p,sig12){
  set.seed(seed)
  z <- rmvnorm(n=n, mean=rep(0,p))
  x <- z%*%sig12              #\Lamb\Xi+\varepsilon
  # generating random \chi_5^2
  chis5 <- 2*rgamma(n=n,shape=2.5,scale=1)
  scal2 <- sqrt(chis5/3.0)                 #1/r
  scal2_mat <- matrix(rep(scal2,p),ncol=p)
  x <- x/scal2_mat  #elliptically distributed x;
  return(x)
}
# -----------------------------------
# skewed error data; #c3
sim.skew <- function(seed,n,p,m,lambphi_12, psi12){
  set.seed(seed)
  # one_p_mat <- matrix(rep(1,p*n), ncol=n)
  sqrt2 <- sqrt(2)
  f <- matrix(rnorm(m*n),ncol=n)
  f_s <- (f*f-1)/sqrt2;   #z_{\xi}, standardized chisquare;f~N(0,1),E(f*f)=1
  e <- matrix(rnorm(p*n),ncol=n)
  e_s <- (e*e-1)/sqrt2; 
  x_i <- lambphi_12%*%f_s + matrix(rep(psi12,n),ncol=n)*e_s
  # generating random \chi_5^2
  chis5 <- 2*rgamma(n=n,shape=2.5,scale=1)
  scal2 <- sqrt(chis5/3.0)
  scal2_mat <- matrix(rep(scal2,p),ncol=p)
  x <- t(x_i) / scal2_mat
  return(x)
}
###########


################################################################################
cond <- NULL
 cond[[1]] <- expand.grid(n_vec=c(40,60,100,200,300,500,1000),
                           m_vec=c(1),
                           p=5)



condition<-do.call(rbind, cond)
index.l<-c(1:nrow(condition))
################################################################################
methodlist<-c("gammaN.model")
dislist<-c("norm","enorm","skew")
it<-c(1:200)

#grid[(grid$method=="ridgeD"&grid$simdist=="enorm"),]

grid <- expand.grid(it,methodlist,dislist,index.l)
grid<-as.data.frame(grid)
colnames(grid)<-c("it","method","simdist","index")

args <- grid[id, ]

do_one <- function(it,method,simdist,index){
  
  m=condition[index,2]
  p=condition[index,3]
  n=condition[index,1]
  ms=m*(m-1)/2
  q=2*p+ms
  ps=p*(p+1)/2
  p2=p*2
  df=ps-q
  p_m=p/m
  
  p4<-p*(p+1)*(p+2)*(p+3)/24 ################# This is used many times later
  
  ################
  ################
  
  set.seed(which(simdist==dislist)*13^2+n*17)
  lamb_uni<-sort(sample(seq(70,95,5),p_m,replace=TRUE))/100
  lamb <- matrix(0,p,m)
  for (j in 1:m){
    lamb[(p_m*(j-1)+1):(p_m*j),j] <- lamb_uni
  }
  
  #lamb_t <- t(lamb) ###############
  phi <- 0.5*diag(m)+matrix(0.5,m,m) #desired phi
  
  eigen_phi <- eigen(phi)
  phi_val <- eigen_phi$values
  phi_vec <- eigen_phi$vectors
  phi_h <- phi_vec%*%diag(sqrt(phi_val))%*%t(phi_vec)
  lambphi_12 <- lamb%*%phi_h   # to simulate data from skewed distribution; phi^0.5
  rm(eigen_phi,phi_val,phi_vec,phi_h) ##################
  
  
  lamblamb <- lamb%*%phi%*%t(lamb) #################
  psi <- diag(p)-diag(diag(lamblamb)) # psi is specified to ensure sigma's diagonals are 1
  psi_vec <- diag(psi)
  psi12 <- sqrt(psi_vec) # to simulate data from skewed distribution
  
  sig <- lamblamb+psi    #BIG SIGMA
  eigen_sig <- eigen(sig)
  sig_vec <- eigen_sig$vectors
  sig_val <- eigen_sig$values
  sig12 <-  sig_vec%*%diag(sqrt(sig_val))%*%t(sig_vec)
  rm(lamblamb,eigen_sig,sig_val,sig_vec) ##################
  
  
  theta_p <- c(rep(lamb_uni,m),psi_vec,phi[lower.tri(phi, diag = FALSE)])
  Sigma0.t <- lamb%*%phi%*%t(lamb) + diag(psi_vec) 
  rm(lamb_uni)###############
  
  result.table<-matrix(NA,nrow=5*3,ncol=22+q*4)
  
  for (rep in 1:5){
    #  print(rep)
    ############################################# 
    ###############simulate x
    #### Below is the SAME for each method 
    seed=(it-1)*5+rep 
    
    if (simdist=='norm'){x <- sim.norm(seed,n,p,sig12)}
    if (simdist=='enorm'){x <- sim.enorm(seed,n,p,sig12)}
    if (simdist=='skew'){x <- sim.skew(seed,n,p,m,lambphi_12,psi12)}
    
    #get sample covariance
#    ps <- p*(p+1)/2
#    pss <- ps*(ps+1)/2 ################## This is never used.
    
    Scov <-cov(x)*(n-1)/n
    vscov <- vech(Scov)
    
    mean_xt<-colMeans(x)###################
    si=apply(x, 1, function(z) vech((z-mean_xt)%*%t(z-mean_xt)))##########
    
    # Note that the calculation of si is also redundant ################## 
    # since it is the same as sigmaele calculated later ##################
    
    #get model-implied covariance and first order derivative
    sigdsig <- function(p,m,theta0){
      ps <- p*(p+1)/2
      p_h <- p/m
      p2 <- p*2
      lamb <- matrix(0,p,m)
      for (j in 1:m){
        p_hj <- (j-1)*p_h+1
        p_hj1 <- j*p_h
        lamb[p_hj:p_hj1,j] <- theta0[p_hj:p_hj1]
      }
      
      psi_vec <- theta0[(p+1):p2] 
      phi <- matrix(1,m,m)
      k <- 0
      if (m>1){
        for (i in 1:(m-1)){
          for (j in (i+1):m){
            k <- k+1
            phi[j,i] <- theta0[p2+k] 
            phi[i,j] <- theta0[p2+k]
          }
        }
      }
      Sigma0 <- lamb%*%phi%*%t(lamb) + diag(psi_vec) #model-implied covariance
      
      # The following is to evaluate (dSigma0/dtheta) 
      # DL is the derivative with Lambda
      DL <- matrix(0,ps,p) 
      lambphi <- lamb%*%phi
      lambphi_t <- t(lambphi)
      for (j in 1:m){
        p_hj <- (j-1)*p_h+1
        p_hj1 <- j*p_h
        for (k in p_hj:p_hj1){
          tt <- matrix(0,p,m)
          tt[k,j] <- 1
          ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
          DL[,k] <- vech(ttt)
        }
      }
      #Dps is the derivative with Psi
      Dps <- matrix(0,ps,p)
      for (j in 1:p){
        tt <- matrix(0,p,p)
        tt[j,j] <- 1
        Dps[,j] <- vech(tt)
      }
      
      #Dphi is the derivative with phi
      if (m>1){
        ms <- m*(m-1)/2 
        Dphi <- matrix(0,ps,ms)  
        k <- 0 
        for (i in 1:(m-1)){
          for (j in (i+1):m){
            k <- k+1 
            tt <- matrix(0,m,m) 
            tt[j,i] <- 1 ; tt[i,j] <- 1 
            ttt <- lamb%*%tt%*%t(lamb)
            Dphi[,k] <- vech(ttt)
          }
        } 
        vdsig <- cbind(DL,Dps,Dphi)
      } else {  vdsig <- cbind(DL,Dps)}
      
      out <- list('Sigma0'=Sigma0, 'vdsig'=vdsig)
      return(out)
    }
    
    #######GLM
    ep <- 0.0001 ################# never used
#    ps <- p*(p+1)/2 ################## already calculated
#    q <- length(theta_p)################## already calculated
#    df <- ps-q################ already calculated
    
#    start_time <- Sys.time()
    ###### W matrix
    ## gammaadf 
    mean_xt <- apply(x,2,mean)
    x_c <-x-matrix(1,n,1)%*%mean_xt
    Scov <-cov(x_c)
    vscov <- vech(Scov)
   # x_c <- scale(x,center=mean_xt,scale = FALSE) ##################
   k=1
   sigmaele<-matrix(NA,nrow=ps,ncol=n) # second order ###############
   order2<-matrix(NA,nrow=ps,ncol=3)
   for(i in 1:p){
     for(j in i:p){
#print(c(i,j))
       sigmaele[k,]<- x_c[,i]*x_c[,j]
       order2[k,]=c(k,i,j)
       k=k+1
     }
   }
    #sigmaele<-si #################### You want to remove si or signmaele to minmize space
  #  order2<-cbind(rep(1:p,each=p),rep(1:p,p))[c(lower.tri(Scov,diag = TRUE)),]########
  #  order2<-cbind(1:nrow(order2),order2) #####################
    
    sigmaij=sigmakl=rowSums( sigmaele)/n #w_ij  ########### They are just vscov you calculated earlier
 #  sigmaij<-sigmakl<-vscov 
    gammaadf<-cov(t(sigmaele))*(n-1)/n ##################
#    sigmaijkl<- (sigmaele%*%t(sigmaele)/n)[lower.tri(gammaadf,diag = TRUE)] ########### This quantity is not used later
    
    k3=1
    sigmaele3<-matrix(NA,nrow=p*(p+1)*(p+2)/6,ncol=n) # 3th order
    order3<-matrix(NA,nrow=p*(p+1)*(p+2)/6,ncol=4)
    for(i in 1:p){
      for(j in i:p){
        for(r in j:p){
           # print(c(i,j,r))
            sigmaele3[k3,]<- x_c[,i]*x_c[,j]*x_c[,r]
            order3[k3,]=c(k3,i,j,r)
            k3=k3+1
          }
        }
      }
    sigmaijr=rowSums( sigmaele3)/n #w_ijr
    
    k=1
    sigmaele4<-matrix(NA,nrow=p*(p+1)*(p+2)*(p+3)/24,ncol=n) # 4th order
    order4<-matrix(NA,nrow=p*(p+1)*(p+2)*(p+3)/24,ncol=5)
    for(i in 1:p){
      for(j in i:p){
        for(r in j:p){
          for(s in r:p){
            #print(c(i,j,r,s))
        sigmaele4[k,]<- x_c[,i]*x_c[,j]*x_c[,r]*x_c[,s]
        order4[k,]=c(k,i,j,r,s)
        k=k+1
         }
        }
       }
    }
    
    kk=1
    index<-matrix(NA,nrow=p,ncol=p)
    for(i in 1:p){
      for(j in i:p){
        index[i,j]=index[j,i]=kk
        kk=kk+1
      }}
    
    D<-duplication.matrix(p) #############
    K<-scale(D,center=F,scale=colSums(D))##############
    gammaNS<-t(K)%*%kronecker(Scov,Scov)%*%K*2##############
    rm(D,K)
    
    ###################################
    
    matrix1<-cbind(x_c*n/(n-1),
                   t(sigmaele)*(n/(n-1))^2,
                   t(sigmaele4)*(n/(n-1))^4)
    
    #### matrix2
    k=1
   # matrix2 <- matrix(0,nrow=p+ps+p*(p+1)*(p+2)*(p+3)/24, ncol=p*(p+1)*(p+2)*(p+3)/24)
    matrix2 <- rep(list(matrix(0,nrow=p+ps+p*(p+1)*(p+2)*(p+3)/24, ncol=1)),p*(p+1)*(p+2)*(p+3)/24)
    #index2 <- matrix(NA,nrow=p*(p+1)*(p+2)*(p+3)/24,ncol=4) #############
    newx<-NULL
    pilist<-NULL

    for(i in 1:p){
      for(j in i:p){
        for(r in j:p){
          for(s in r:p){
          #  print(k)
 #        index2[k,]=c(i,j,r,s) ############### index2 is just order4[,2:5], no need to do again
            ############## for x's
            freq = combn(c(i,j,r,s), 1)
            freq = count.dups(t(freq))
            freq = as.matrix(freq)
            for (ind in 1:nrow(freq)){
              row.n = freq[ind,1]
              delete=which(c(i,j,r,s)==freq[ind,1])
              content = c(i,j,r,s)[-delete[1]]
              content=sort(content)
              matrix2[[k]][row.n,]=-freq[ind,2]*sigmaijr[which(order3[,2]==content[1]&order3[,3]==content[2]&order3[,4]==content[3])]
            }
            #######################
            ##### second order # starts at (p+1)th row
            freq2 = combn(c(i,j,r,s), 2)
            freq2 = count.dups(t(freq2))
            freq2 = as.matrix(freq2)
            for (ind2 in 1:nrow(freq2)){
              row.n = which(order2[,2]==min(freq2[ind2,1:2]) & order2[,3]==max(freq2[ind2,1:2]))
              delete1=which(c(i,j,r,s)==freq2[ind2,1])
              delete2=which(c(i,j,r,s)==freq2[ind2,2])
              if(freq2[ind2,1]!=freq2[ind2,2]){
              content = c(i,j,r,s)[-c(delete1[1],delete2[1])]} else {
                content = c(i,j,r,s)[-delete1[1:2]]
              }
            matrix2[[k]][(p+row.n),]=-freq2[ind2,3]*sigmaij[which(order2[,2]==min(content) &order2[,3]==max(content))]
              }
            #######################
            ##### fourth order # starts at p+ps+1th row
            row.n=which(order4[,2]==i&order4[,3]==j&order4[,4]==r&order4[,5]==s)
            matrix2[[k]][(p+ps+row.n),]=1
            newx[[k]]=matrix1%*%matrix2[[k]]
            k=k+1
        }
      }
     }
    }

   # newx=sapply(matrix2,function(matrix2) matrix1%*%matrix2)
    
    rm(matrix1,matrix2) ############# Remove objects not to be used later

################ compute v
    ######### recommended code #################
    
    index2= order4=order4[,2:5]    
    
    #### Now calculate the 24-term part and add to the existing matrix
    #  t0<-proc.time()
    V=matrix(0,p4,p4)
        term24index<-unlist(permn(c(1,2,3,4)))
        i<-term24index[seq(from=1,by=4,length.out=24)]
        j<-term24index[seq(from=2,by=4,length.out=24)]
        ii<-term24index[seq(from=3,by=4,length.out=24)]
        jj<-term24index[seq(from=4,by=4,length.out=24)]
    
        for (ele in 1:p4)
        {
          Scov.4rows<-Scov[order4[ele,],]
          for (ele2 in ele:p4)
          {
            Scov.block<-Scov.4rows[,order4[ele2,]]
            V[ele,ele2]<-V[ele2,ele]<-V[ele2,ele]+sum(
              Scov.block[1,i]*Scov.block[2,j]*Scov.block[3,ii]*Scov.block[4,jj])
          }
        }
     #   proc.time()-t0

############ estimate a
###################################
    pi.telda=cov(Reduce(cbind,newx)) ###############
#    v1=chol2inv(chol(v)) ##################
        
    V.upper<-chol(V) ###############
        
        newgamma <- (gammaadf-gammaNS)

        select.f= function(index2){
          newgamma[index[index2[1],index2[2]],index[index2[3],index2[4]]]
        }
        newnewgamma<-apply(index2,1,select.f)
   
             y<-backsolve(r=V.upper,x=newnewgamma,transpose = TRUE)##############
        tr<-sum(diag(backsolve(r=V.upper,x=t(backsolve(r=V.upper,x=pi.telda,transpose = TRUE)),
                       transpose = TRUE)))#############
             nu=t(y)%*%y/tr-1/(n-3)##############
             weight=1/nu/(1/nu+n-3)
             if (nu<0){ weight=1}
             
             a.n=1-weight
             alist=c(0,a.n,1)
             ###########################
             ###########################
             ###### W matrix
             ## gammaadf.raw 
             mean_xt <- apply(x,2,mean)
             x_c <- x-matrix(1,n,1)%*%mean_xt
             Scov <-cov(x)
             vscov <- vech(Scov)
             k=1
             sigmaele<-matrix(NA,nrow=ps,ncol=n) # second order
             order2<-matrix(NA,nrow=ps,ncol=3)
             for(i in 1:p){
               for(j in i:p){
                 # print(c(i,j))
                 sigmaele[k,]<- x_c[,i]*x_c[,j]
                 order2[k,]=c(k,i,j)
                 k=k+1
               }
             }
             sigmaij=sigmakl=rowSums( sigmaele)/n #w_ij
             
             sigmaijkl=c()
             gammaadf=matrix(NA,nrow=ps,ncol=ps)
             k=1
             for(i in 1:ps){
               for(j in i:ps){
                 sigmaijkl[k]<-sum(sigmaele[i,]*sigmaele[j,])/n
                 gammaadf[i,j]<-sigmaijkl[k]-sum(sigmaele[i,])*sum(sigmaele[j,])/n^2
                 gammaadf[j,i]<- gammaadf[i,j]
                 k=k+1
               }
             }
             
             
             ###############################
             
             set.seed(which(simdist==dislist)*13^2+n*17)
             
             for (aid in 1:3){
               a<-alist[aid]
               # print(a)
               
               theta0=theta_p-runif(q,0.1,0.2)
               if (sum(theta0<=0 | theta0>=1)!=0){
                 theta0[which(theta0<=0 | theta0>=1)]=theta_p[which(theta0<=0 | theta0>=1)]}
               
               sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
               vdsig<-sig$vdsig  #first order, sigma-dot
               Sigma0<-sig$Sigma0 # model_implied covariance matrix
               vsig0 <- vech(Sigma0)
               
               diconverg=0
               ###### 
               for(t in 1:200){
                 ###### ## gammaN.model; updated with theta0 and sigma0
                 if (method=="gammaN.model"){
                   
                   gammaNm<-matrix(NA,nrow=ps,ncol=ps)
                   for(i in 1:p){
                     for(j in i:p){
                       for (k in  1:p){
                         for (l in 1:p){
                           if ( index[k,l]>=index[i,j] & l>=k){
                             #  print(c(i,j,k,l))
                             gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                             gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                           }
                         }
                       }
                     }
                   }
                   
                   weight<-try(solve(a* gammaadf+(1-a)*gammaNm))
                   if( class( weight)=="try-error"){
                     diconverg <- 1   # not converge
                     break
                   }
                 }
                 
                 
                 stdi <- try(solve(t(vdsig) %*% weight%*%vdsig))
                 if( class(  stdi)=="try-error"){
                   diconverg <- 1   # not converge
                   break
                 }
                 eresdu <- vscov-vsig0 #s-sigma
                 dtheta <- t( eresdu) %*% weight %*%  vdsig%*%  stdi
                 theta0 <- theta0 + dtheta
                 delta <- max(abs(dtheta))
                 
                 sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
                 vdsig<-sig$vdsig 
                 Sigma0<-sig$Sigma0
                 vsig0 <- vech(Sigma0)
                 
                 if(delta<=ep) {
                   #      sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
                   #### ## gammaN.model; updated with theta0 and sigma0
                   if (method=="gammaN.model"){
                     gammaNm<-matrix(NA,nrow=ps,ncol=ps)
                     for(i in 1:p){
                       for(j in i:p){
                         for (k in  1:p){
                           for (l in 1:p){
                             if ( index[k,l]>=index[i,j] & l>=k){
                               #  print(c(i,j,k,l))
                               gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                               gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                             }
                           }
                         }
                       }
                     }
                     
                     weight<-try(solve(a* gammaadf+(1-a)*gammaNm))
                     if( class( weight)=="try-error"){
                       diconverg <- 1   # not converge
                       break
                     }
                   }#get final weight (gamma.N.modelbased)
                   #  vdsig<-sig$vdsig 
                   # Sigma0<-sig$Sigma0
                   # vsig0 <- vech(Sigma0)#first order
                   break};
                 
               }
               
               if(t<200){
                 if(diconverg==0){  #converge within 300 iteration 
                   # Test start here;
                   dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
                   if( class( dtd)=="try-error"){
                     break
                   }
                   r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
                   r.SE <- sqrt(diag(r.Var)/(n-1))
                   Var<-dtd
                   SE <- sqrt(diag(Var)/(n-1))
                   dwe<-t(vdsig)%*%weight
                   U<- weight-t(dwe)%*%dtd%*%dwe
                   ugamma <-  U %*% gammaadf
                   
                   rank<-qr(ugamma)$rank
                   
                   ###### calculate a2    
                   U5=sqrtm(U)
                   wi=apply(si, 2, function(si)  U5%*%si)
                   wi=Re(wi)
                   wbar=rowMeans(wi)
                   yi=wi-wbar
                   H= yi%*%t(yi)
                   
                   Dmatrix<-function(i){
                     t(yi[,i])%*%yi[,i]
                   } 
                   D<-lapply(1:n, Dmatrix)
                   D<-unlist(D)
                   a2c<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)*sum(diag(H%*%H))-
                                                   n*(n-1)*sum(D^2)+sum(diag(H))^2 )
                   ####################################
                   
                   Fstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))
                   
                   Tstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))*(n-1)
                   
                   c1 <- sum(diag(ugamma ))/df
                   rTstats1 <- Tstats/c1
                   
                   c2<-sum(diag(ugamma ))/qr(ugamma)$rank
                   rTstats2 <- Tstats/c2
                   
                   ugamma2 <- ugamma %*% ugamma
                   c3 <-sum(diag( ugamma2))/sum(diag(ugamma )) 
                   df2<-(sum(diag(ugamma )))^2/ sum(diag( ugamma2))
                   rTstats3 <- Tstats/c3
                   
                   c4 <-a2c/sum(diag(ugamma )) 
                   df2.2<-(sum(diag(ugamma )))^2/ a2c
                   rTstats4 <- Tstats/c4
                   
                   ms<-Mardia.tran1(x)$skew
                   mk<-Mardia.tran1(x)$kurt
                   
                   result.table[(rep-1)*3+aid,]<-c(rep,p,m,n,q,df,df2,df2.2,rank,ms,mk,simdist,method,a,
                                                   nu,a.n,
                                                   Fstats,
                                                   Tstats,rTstats1,rTstats2,rTstats3,rTstats4,
                                                   theta_p,c(theta0),
                                                   c(r.SE),c(SE))
                   
                 }  
                 ###################################
               }
               
               ###################################
               
             }
  }
  
  fn <- paste0("result/",simdist,"/p",p,"m",m,"n",n,"result",it,".csv")
  write.table(result.table, file = fn,  sep = ",", col.names = FALSE,append = TRUE)
}
########################
do.call(do_one, args)
