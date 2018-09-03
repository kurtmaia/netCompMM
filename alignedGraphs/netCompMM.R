library(BayesLogit)

lBETA <- function(x) sum(lgamma(x)) - lgamma(sum(x))
BETA <- function(x) exp(sum(lgamma(x)) - lgamma(sum(x)))


get_pi_group <- function(pi.h,vh1,vh2){
  pi1 <- vh1%*%pi.h
  pi2 <- vh2%*%pi.h
  return(list(pi1,pi2))
}

get_pi_diff <- function(pi.h,vh1,vh2, toMatrix = TRUE){
  pi.group <- get_pi_group(pi.h,vh1,vh2)
  if (toMatrix) {return(edge2matrix(pi.group[[1]]-pi.group[[2]]))}
  else {return(pi.group[[1]]-pi.group[[2]])}
}

get_local_hypt <- function(pi.h,vh1,vh2,phoy, toMatrix = TRUE,...) {
  pL_y <- get_pi_group(pi.h,vh1,vh2)
  pL <- phoy*(pL_y[[1]]) + (1-phoy)*(pL_y[[2]])
  
  y_1 <- ((pL_y[[1]] - pL)^2)/pL + (((1-pL_y[[1]]) - (1-pL))^2)/(1-pL)
  y_2 <- ((pL_y[[2]] - pL)^2)/pL + (((1-pL_y[[2]]) - (1-pL))^2)/(1-pL)
  
  pl2 <- phoy*y_1+(1-phoy)*y_2
  pl <- sqrt(pl2)
  if (toMatrix) {return(edge2matrix(pl,...))}
  else {return(pl)}
 }

get_lkhd_expressions <- function(x){
  if (x == "binom"){
    ans <- "dbinom(t$x,t$n,param,log = TRUE)"
  } else if (x == "poisson") {
      ans <-  "dpois(t$x,param,log = TRUE)"
  }
}

threshold_graphs <- function(x, threshold_level = .2){
ans <- x*0
for (i in 1:(nrow(x)-1)){
  for (j in (i+1):ncol(x)){
    p <- x[i,j]/min(x[i,i],x[j,j])
    p <- ifelse(is.na(p),0,p)
    ans[i,j] <- ifelse(p >= threshold_level,1,0)
    ans[j,i] <- ans[i,j]
  }
}
diag(ans) <- 1
return(ans)
}

binarize <- function(a,threshold = 0) {
  if (sum(a$n)>0) {
    p <- a$x/a$n
    an <- (a$n>0)*1
  }
  else{
    p <- a$x/max(a$x)
    an <- 1
  }
  p[is.nan(p)] <- 0 
  p <- (p > threshold)*1
  return(data.frame(n = an, x = p))
}



div <- function(dt) {
  ans <- ifelse(dt$n==0,0,dt$x/dt$n)
  ans[ans==1] <- 1 - 1e-3
  ans[ans==0] <- 1e-3
  return(ans)
}
  


netCompMM <- function(Ai,PopId, pH0 = .5, a1 = 2, a2 = 2, H = 5, R = 4, heterogenity = FALSE,threshold = NULL, family = "poisson", link = "exp", Tsamples = 500, burnIn = round(Tsamples*.05), alpha = 1, beta = 1,epsilon = .01,print.iter = FALSE, print.nys = FALSE,by.print = ceiling(.05*Tsamples),matrix.plot = FALSE,pp = TRUE, nCores = 4) {

if (!heterogenity){
  PopId <- data.frame(user = names(PopId), PopId = PopId, net.name = names(PopId),stringsAsFactors = FALSE)
} else {
  PopId <- PopId[,c("user","PopId","net.name")]
}


PopId <- dplyr::filter(PopId, net.name %in% names(Ai))

lkhd <- get_lkhd_expressions(family)

Words <- colnames(Ai[[1]])
n1 <- sum(PopId$PopId==1); n2 <-  sum(PopId$PopId==2)
n <- n1 + n2


pH1 <- 1 - pH0
a.H <- rep(1/H,H)
v.h <- list(rep(1/H,H),rep(1/H,H))

V <- ncol(Ai[[1]])
E <- V*(V-1)/2

muL <- rep(0,V*(V-1)/2)
sigmaL <- rep(1,V*(V-1)/2)
Z <- rnorm(E,muL,sqrt(sigmaL))


Ai <- llply(Ai, function(g) order_graphs(g,Words))
if (!is.null(threshold)) {
  Ai <- llply(Ai, function(x) threshold_graphs(x,threshold))
}
registerDoMC(nCores)
binom_input <- occ.matrix2dataframe(Ai,Words,.parallel = pp)


gen_start <- function(pitmp){
  dtmp0 <- matrix(0,V,V)
  Stmp <- log(pitmp/(1-pitmp))
  dtmp0[lower.tri(dtmp0)] <- Stmp - Z
  dtmp0 <- dtmp0 + t(dtmp0)
  dtmp <- eigen(dtmp0)
  Ld <- dtmp$values[1:R]
  X <- dtmp$vectors[,1:R]
  Ld[which(Ld<0)] <- 1
  Xbar <- X%*%diag(sqrt(Ld))
  vi <- 1/Ld
  Dd.h <- dtmp0[lower.tri(dtmp0)]
  S <- Z + Dd.h
  ans <- list(Xbar,Dd.h,Ld,vi,S)
  names(ans) <- c("Xbar","Dd.h","Ld","vi","S")
  return(ans)
}
pis <- llply(1:H, function(i) runif(E))

STARTS <- llply(pis, gen_start)
Xbar <- llply(STARTS,function(i) i$Xbar)
Ld <- llply(STARTS,function(i) i$Ld)
vi <- llply(STARTS,function(i) i$vi)
S <- llply(STARTS,function(i) i$S)
Dd.h <- laply(STARTS,function(i) i$Dd.h)

omega <- matrix(NA, H, E)
HypT <- NULL
pi <- list()
Y <- list()
V.H <- list()
Dzao <- list()
Piaglo <- list()
LD <- list()
Clusters <- list()
NYI <- list()
pho.y <- list()
local.hyp <- matrix(0,1,E)

t0 <- Sys.time()


Gi.sampler <- function(t,vhii,pi){
  if (family == "binom"){
  part1 <- laply(1:H,function(h) sum(dbinom(t$x, t$n, pi[[h]], log = TRUE))) 
} else if (family == "poisson") {
  part1 <- laply(1:H,function(h) sum(dpois(t$x, pi[[h]], log = TRUE)))
}
  prodsF <- part1 + log(vhii)
  prob <- softmax(prodsF)
  return(rcateg(1,prob))
}




PopId$dish <- factor(rcateg(n,a.H),levels = 1:H)
Gi <- PopId$dish %>% as.numeric
names(Gi) <- PopId$net.name

for (iter in 1:Tsamples){

  pho.y[[iter]] <- rbeta(1, alpha + n1, beta + n2)
  
    if (iter > 1) {
      Gi[names(Ai)] <- try(laply(names(Ai), function(x) Gi.sampler(binom_input[[x]],vhi[[x]],pi),.parallel = pp))
        if (str_split(Gi," ")[[1]][1] == "Error") {
          print('probleminha')
           return(list(binom_input,PopId,pi,vhi,S,Y,YT)); break}
      PopId$dish <- factor(Gi,levels = 1:H)
    }

  vals <- llply(1:H,function(i) names(Gi[Gi==i]))
  n.y.i <- xtabs(~PopId+dish+ user,PopId)
  m.y <- xtabs(~PopId+dish,PopId)
  m.h <- xtabs(~dish,PopId)

  if (print.nys) {print(m.y)}


 add_or_zero <- function(x) {
    if (length(x) > 0) {
      return(add(x))
    } else {
      return(data.frame(n = rep(0,E), x = rep(0,E)))
    }
  }


  get_mean_zero <- function(x) {
    if (length(x) > 0) {
      return(add(x)/length(x))
    } else {
      return(data.frame(n = rep(0,E), x = rep(0,E)))
    }
  }

  conjY <- llply(vals,function(i) add_or_zero(binom_input[i]))

  Hvals <- seq(1,H)[laply(vals, function(i) length(i) > 0)]
  if (iter == 1) {Hvals = 1:H}
  neg.Hvals <- setdiff(1:H,Hvals)

if (link == "exp"){
  Y <- laply(conjY,function(xx) xx$x)
  YT <- laply(vals,function(i) get_mean_zero(binom_input[i])$x)
} else if (link == "logit") {
  Y <- laply(conjY,function(xx) xx$n)
  YT <- laply(conjY,function(xx) xx$x)
}

  omega <- laply(1:H, function(h) laply(1:E,function(i) rpg.gamma(1,Y[h,i],S[[h]][i])),.parallel = pp)
  sigma.zl <- 1/(1/sigmaL+apply(omega,2,sum))
  sum.H <-  YT - omega*Dd.h
  mu.zl <- sigma.zl*(muL/sigmaL+apply(sum.H,2,sum))
  Z.l <- unlist(lapply(1:E,function(i){return(rnorm(1,mu.zl[i],sqrt(sigma.zl[i])))}))

  
  for (h in 1:H){
    Om <- edge2matrix(omega[h,],lower = FALSE)
    Z.ll <- edge2matrix(Z.l, lower = FALSE)
    Y.l <- edge2matrix(Y[h,],lower = FALSE)
    YT.l <- edge2matrix(YT[h,],lower = FALSE)

      for (v in 1:V){
      Om.v <- diag(Om[v,-v])
      Y.v <- Y.l[-v,v]
      YT.v <- YT.l[-v,v]
      Z.v <- Z.ll[-v,v]
      eta.v <- t(Xbar[[h]][-v,])%*%(YT.v  - diag(Om.v)*Z.v)
      Sigma0 <- t(Xbar[[h]][-v,])%*%Om.v%*%Xbar[[h]][-v,] 
      diag.inv <- diag(1/Ld[[h]])
      Sigma <- Sigma0 + diag.inv
      Sigma.inv <- try(chol2inv(chol(Sigma)))
      if (class(Sigma.inv) == "try-error") {
        print(c(h,v))
        return(list(Sigma0,Xbar[[h]],Om.v,Ld[[h]])); break}
      mean <- Sigma.inv%*%eta.v
        tmp <- mvrnorm(1,mean,Sigma.inv)
      Xbar[[h]][v,] <- tmp
      }

    theta_prod <- function(rmInd,vvi = vi[[h]]) {
      ans <- NULL
      for (i in 1:R){
        if (i < rmInd) {
          ans <- c(ans,prod(vvi[1:i]))
        } else if(i >= rmInd){
          ans <- c(ans,prod(vvi[1:i])/vvi[rmInd])
        }
      }
      return(ans)
    }
    
   S.xbarh <- as.matrix(apply(Xbar[[h]],2,function(i) sum(i^2)))
   thetas <- as.matrix(ldply(1:R, theta_prod))
   thetas[lower.tri(thetas)] <- 0 
   bb <- (thetas%*%S.xbarh)/2 + 1
   as <- c(a1,rep(a2,R-1))
  vi[[h]] <- unlist(lapply(1:R,function(i) rgamma(1,as[i]  + V*(R+1-i)/2,scale = bb)))
   Ld[[h]] <- cumprod(1/vi[[h]])
     Dd.tmp <- Xbar[[h]]%*%t(Xbar[[h]])
     Dd.tmp <- Dd.tmp[lower.tri(Dd.tmp)]
     Dd.h[h,] <- Dd.tmp
     S[[h]] <- Z.l + Dd.tmp
     if (link == "exp"){
         pi[[h]] <- exp(S[[h]])
    } else if (link == "logit") {
      pi[[h]] <- 1/(1+exp(-S[[h]]))
      pi[[h]][pi[[h]]>=1] <- 1 - 1e-5
      pi[[h]][pi[[h]]<=0] <- 1e-5   
    }
  }


  LD[[iter]] <- Ld
  Piaglo[[iter]] <- pi
  Clusters[[iter]] <- Gi


  if (iter > burnIn){
  local.tmp <- get_local_hypt(do.call('rbind',pi),v.h[[1]],v.h[[2]],pho.y[[iter]],toMatrix = FALSE)
  local.hyp <- local.hyp + (local.tmp>epsilon)
}
 

  H1 <- sum(apply(a.H + m.y,1,lBETA) - lBETA(a.H))
  H0 <- lBETA(a.H + m.h) - lBETA(a.H) 
  prob3.2 <- 1/(1+(pH0/pH1)*exp(H0-H1))

 
TT <- rbinom(1,1,prob3.2)
if (TT == 0){vcal <- rdirichlet(1, a.H + m.h) ; v.h <- list(vcal,vcal)}
if (TT == 1){v.h <- list(rdirichlet(1, a.H + m.y[1,]),rdirichlet(1, a.H + m.y[2,]))}

if (heterogenity){
  vhi <- dlply(PopId,.(net.name), function(x) rdirichlet(1,n.y.i[x$PopId,,x$user] + v.h[[x$PopId]]))
} else {
  vhi <- dlply(PopId,.(net.name), function(x) v.h[[x$PopId]])
}

Dzao[[iter]] <- Dd.h
V.H[[iter]] <- v.h
NYI[[iter]] <- n.y.i
# print(v.h)
HypT <- c(HypT,TT)

if (print.iter | iter%%by.print == 0) {cat(iter,"| p =",prob3.2,"| Pr(H1|-) =",sum(HypT)/iter,"\n")}

if (iter%%by.print == 0){
  t1 <- Sys.time()
  t.diff <- as.numeric(difftime(t1,t0,units = 'secs'))
  reman <- (Tsamples-iter)*(t.diff/iter)
  estimate <- as.POSIXct(as.numeric(t1,units = 'secs') + reman,origin = "1970-01-01")
  cat(paste("Iteration ",iter,' of ',Tsamples,'.',sep=""),
      '\n',
      paste("Expected to finish at ",estimate,'.',sep=""),'\n')

  if (matrix.plot) {
      mplot.dt <- laply(pi,function(p) edge2matrix(p)) %>% melt()
      mplot <- ggplot(mplot.dt,aes(x=Var2,y=Var3,fill=value)) + geom_tile()+facet_wrap(~Var1) 
      plot(mplot)
  }
}

}

get.cluster <- function(i,C) as.numeric(names(sort(table(laply(C, function(x) x[i])),decreasing = TRUE)[1]))
clus <- laply(1:n, function(i) get.cluster(i, Clusters[burnIn:Tsamples]))
names(clus) <- names(Ai)

PIans <- do.call('rbind',lapply(1:H, function(i) pi.hat(Piaglo[burnIn:Tsamples],i)))
Dans <- do.call('rbind',lapply(1:H, function(i) pi.hat(Dzao[burnIn:Tsamples],i)))
Local.HypT <- local.hyp/(Tsamples-burnIn)

Fans <- list(mean(HypT[burnIn:Tsamples]),V.H[burnIn:Tsamples],LD[burnIn:Tsamples],PIans,Dans,Words,mean(unlist(pho.y)[burnIn:Tsamples]),Local.HypT, clus,Clusters[burnIn:Tsamples],NYI[burnIn:Tsamples])
names(Fans) <- c("HypT","V.H","LD","Probs","Distances","Words","pho.y","Local.HypT","Clusters","PercClus","NYI")
return(Fans)
}

