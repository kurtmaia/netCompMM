packages <-  c("stringr","MASS" ,"reshape2" ,"plyr" ,"tm" ,"gtools" ,"ggplot2","gtools","doMC","VGAM")
# ins <- lapply(packages,install.packages)
# rm(ins)
req <- lapply(packages, require, character.only = TRUE)
rm(req)

# Defining time variables
def.tmcatt <- function(x,intervals){
  tmp.tm <- as.numeric(x)
  tm.catt <- matrix(NA,nrow = length(x), ncol = 1)
  for (i in 1:(length(intervals)-1)){
      tmp <- which(tmp.tm < intervals[i+1] & tmp.tm >= intervals[i])
      tm.catt[tmp] <- i
    }
  return(tm.catt)
}

# Count words per period
count.words <- function(xx,by){ 
  xx <- tolower(xx)
  b <-  unique(by)
  words <- NULL
  wrd <- NULL
    for (i in b){
      ln <- which(by==i); if (length(ln) == 0) {next}
      status <- xx[ln]
      if (nchar(as.character(status))==0){next}
      s_clean <- tm_map(Corpus(VectorSource(status)), removeWords, stopwords("pt"))
      dtm <- DocumentTermMatrix(s_clean,control = list(removePunctuation = TRUE,removeNumbers = TRUE))
      dtm2 <- as.matrix(dtm)
      wrd$freq <- colSums(dtm2)
      if (length(wrd$freq)==0){wrd$freq <- 0}
      wrd$by.cat <- i
      wrd$terms <- ifelse(wrd$freq==0,"NONAME",colnames(dtm2))
      words <- rbind(words,do.call(cbind.data.frame,wrd))
      print(i)
    }
  return(words)
}

ind2sub <- function(m, ind){
  r <- ((ind-1) %% m) + 1
  c <- floor((ind-1) / m) + 1
  return(cbind(r,c))
}

# Proximity matrix by var
prox.matrix <- function(a, by){
  bb <-  unique(by)
  P <- list()
  m <- 0
  for (b in bb){
    m <- m + 1
    idx <- which(by==b)
    aa <- a[idx,]
    v <- crossprod(aa, aa)
    vv <- v*upper.tri(v)
    index <- which(vv > 0,arr.ind = T, useNames = T)
  p <- v*0
  if (length(index)==0) {next}
  if (!is.matrix(index)) {index <- ind2sub(dim(vv)[1],index)}
    for (k in 1:nrow(index)){
      i <- index[k,1]; j <- index[k,2]
      u <- v[i,i] + v[j,j] - v[i,j]
      p[c(i,j),c(j,i)] <- v[i,j]/u
      # p[c(i,j),c(j,i)] <- 1
    }
    diag(p) <- 0
    P[[m]] <- p
    print(b)
  }
  names(P) <- paste("by.",bb,sep = "")
  return(P)
}

add <- function(x) Reduce("+", x)

print_graph <- function(m, order,title){
  m$terms <- factor(m$terms,order)
  m$terms.1 <- factor(m$terms.1,order)
  graph <- ggplot(m, aes(terms.1, terms),levels = order) + geom_tile(aes(fill = value),colour = "white") + scale_fill_gradient(low = "white",high = "steelblue", guide = FALSE) +  ggtitle(title)
  p <- graph + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + coord_fixed()

  return(p)
}

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...))
}


to.groups <- function(DD,w.names){
  ans <- NULL
  for (i in 1:nrow(DD)){
    dist <- matrix(0,V,V)
    dist[lower.tri(dist)] <- DD[i,]
    dist <- dist + t(dist)
    fit <- cmdscale(dist, k=2)
    points <- as.data.frame(fit)
    points$latent <- i
    points$names <- w.names
    ans <- rbind(ans, points)
  }
  ans$names <- factor(ans$names)
  ans$latent <- factor(ans$latent)
  return(ans)
}


pi.hat <- function(xx,i) apply(ldply(xx,function(x) x[[i]]),2,mean)

edge2matrix <-  function(x, lower = TRUE,dnames=NULL){
  E <- length(x)
  V <- (1 + sqrt(1 + 8*E))/2
  ans <- matrix(0,V,V)
  ans[lower.tri(ans)] <- x
  if (!is.null(dnames)) {dimnames(ans) <- list(dnames,dnames)}
  if (lower) {return(ans)}
  else {return(ans + t(ans))}
}

get_results <- function(model){
  con <- file(paste("./progs/python/results",model,".txt", sep = ""), "r", blocking = FALSE)
  lda_results <- readLines(con)
  close(con)
  # if (model > 1){
  # tt <- str_split(lda_results,"';'")
  # cod <- unlist(lapply(tt,function(i) i[2]))
  # lda_results <- unlist(lapply(tt,function(i) i[1]))
  # cod <- unlist(lapply(cod, function(i) as.numeric(str_extract_all(i,"[.0-9]+")[[1]])))
  # validate <- twdt[twdt$cod %in% cod,]
  # }
  # else {
  validate <- twdt #}
  extract_vals <- function(x){
    val <- matrix(str_extract_all(x,"[.0-9]+")[[1]],ncol = 2, byrow = T)
    ans <- as.data.frame(t(as.numeric(val[,2])))
    colnames(ans) <- val[,1]
    return(ans)
  }
  Lda <- llply(lda_results,extract_vals,.progress = "text")
  Lda <- rbind.fill(Lda,.progress = "text")
  Lda[is.na(Lda)] <- 0
  # sample_ind_c <- read.table("./progs/python/Validate_data_user.csv",header = F,skipNul = T)$V2
  # validate <- twdt[which(twdt$cod %in% sample_ind_c),]
  
  
  ans <- validate
  ans <- cbind(ans,Lda)
  ans$side <- as.factor(ans$side)
  ans$day <- as.Date(ans$date,  origin = "1960-01-01")
  ans$model <- model
  return(ans)
}

# Words
get_words <- function(file,nT = 3){
  con <- file(file, "r", blocking = FALSE) 
  words <- readLines(con)
  close(con)
  topics <- str_split(str_sub(words,2,str_length(words)),"[]]")[[1]]
  trim.leading <- function (x)  sub("^\\s+", "", x)
  get_probs <- function(x){
    xx <- str_split(x,"\\[")[[1]]
    name <- str_extract_all(xx[1],"[.0-9]+")[[1]]
    tmp <- data.frame(matrix(str_split(xx[2],",")[[1]],ncol = 2, byrow = T))
    colnames(tmp) <- c("Words", "Prob")
    tmp[,2] <-  as.numeric(str_extract_all(tmp[,2],"[.0-9e-]+"))
    tmp[,1] <- trim.leading(str_replace_all(tmp[,1],"[(]u'",""))
    tmp[,1] <- as.factor(str_replace_all(tmp[,1],"'",""))
    # tmp$topic <- as.factor(name)
    return(tmp)
  }
  topic <- NULL
  for (i in 1:nT){
    topic[[i]] <- get_probs(topics[i])
    if (i == 1) {tt <- topic[[1]]}
    else if(i > 1) {
      tt <- merge(x = tt, y = topic[[i]], by = "Words", all = TRUE)
      colnames(tt) <- c("Words",(1:i)-1)
      }
  }
  # topic1 <- get_probs(topics[1])
  # topic2 <- get_probs(topics[2])
  # topic3 <- get_probs(topics[3])
  # t1 <- merge(x = topic1, y = topic2, by = "Words", all = TRUE)
  # t2 <- merge(t1,topic3,by = "Words", all=T)
  ans <- melt(tt,id.vars = "Words")
  colnames(ans) <- c("Words","topic","Prob")
  return(ans)
}

prep_input <- function(cat1, cat2 = NULL, gg = graphs,s = side,time = FALSE,...){
  # cat1 = "golpe"; cat2 = "forapt" ; gg = graphs; s = side[side$user%in%g1,]
  if (time) {return(prep_input_time(cat1,cat2,gg=gg,s=s))}
uu <- s[s$side%in%c(cat1,cat2),]
Ai <- gg[uu$user]
if (is.null(cat2)) {
  PopId <- as.numeric(rbinom(nrow(uu),1,.5))+1
  } else {
  PopId <- NULL
  PopId[uu$side %in% cat1] <- 1
  PopId[uu$side %in% cat2] <- 2
  }
names(PopId) <- uu$user
ans <- list(PopId,Ai)
names(ans) <- c("PopId","Ai")
return(ans)
}

prep_input_time <- function(cat1, cat2 = NULL, gg = graphs,s = side,...){
  # cat1 = "golpe";cat2 <- "fora";gg = graphs; s = side[side$user%in%g1,]
  uu <- s[s$side%in%c(cat1,cat2),] 
  Ai <- gg[uu$net.name%>% as.character()]
  if (is.null(cat2)) {
    nuser = length(unique(uu$user))
    n.2 = round(nuser/2)
    PopId <- data.frame(user = unique(uu$user),
                        PopId = c(rep(1,n.2),rep(2,nuser-n.2)),stringsAsFactors = FALSE)
    PopId <- merge(uu,PopId,by = 'user')
  } else {
    PopId <- uu
    PopId$PopId[uu$side %in% cat1] <- 1
    PopId$PopId[uu$side %in% cat2] <- 2
  }
  ans <- list(PopId,Ai)
  names(ans) <- c("PopId","Ai")
  return(ans)
}

sparsity <- function(x){
  xx <- x[lower.tri(x)]
  return( sum(xx == 0)/length(xx))
}

sample.user <- function(group, nProb = .1, s = side,time = FALSE,rand.time=FALSE){
  # group = "golpe"; nProb = .1; s = side
  if (time){
    if (rand.time) {
      s <- s[,c("net.name","side")]
      colnames(s)[1] <- "user"  
    } else {
    s <- unique(s[,c("user","side")])
  }
  }
  a <- s$user[which(s$side %in% group)]
  ans <- sample(a,length(a)*nProb)
  return(ans)
}
get.HypT <- function(...,graphs,PopId){
  polya.nets(graphs,PopId, R = 20, H = 7,Tsamples = 3000,burnIn = 1000)$HypT
}

runParallel <- function(nServers, nTimes, nProb, prefix = "Hypt", group = "forapt", ...){
  for (i in 1:nTimes){
    print(i)
    aa <- llply(1:nServers, function(i) prep_input(group,gg = graphs[sample.user(group,s.size = nProb/100)]))
    registerDoMC(nServers)
    HypT <- llply(aa,function(x) get.HypT(graphs = x$Ai,PopId = x$PopId),.parallel = TRUE)
    registerDoMC(1)
    saveRDS(unlist(HypT),paste(paste(prefix,nProb,as.numeric(Sys.time()),sep = "_"),".RDS", sep = ""))
  }
}

runParallelSpark <- function(nServers, nTimes, nProb, prefix = "Hypt", group = "forapt", ...){
  for (i in 1:nTimes){
    aa <- llply(1:nServers, function(i) prep_input(group,gg = graphs[sample.user(group,s.size = nProb/100)]))
    HypT <- spark.lapply(aa,function(x) get.HypT(graphs = x$Ai,PopId = x$PopId),.parallel = TRUE)
    saveRDS(unlist(HypT),paste(paste(prefix,nProb,as.numeric(Sys.time()),sep = "_"),".RDS", sep = ""))
  }
}

mergeParallel <- function(pattern, rm = FALSE, ...){
  fn <- list.files(pattern = pattern, ...)
  ans <- unlist(lapply(fn,readRDS))
  if (rm == TRUE) {file.remove(fn)}
  return(ans)
}

genHC <- function(nProb){
  paste(" hc_add_series_density(den[[\"",nProb,"\"]],area = TRUE, name = ", nProb,") %>%",sep="")
}


gen_harmonic <- function(k,s) sum((1:k)^(-s)) 
dzipf <- function(x,s,N = NULL){
  if (is.null(N)) {x^(-s)/zeta(s)
    } else {
  x^(-s)/gen_harmonic(N,s) 
    }
}
pzipf <- function(x,s=NULL,N = NULL){
  if (is.null(S)) {get_s(x)}
  if (is.null(N)) {gen_harmonic(x,s)/zeta(s)
  } else {
  gen_harmonic(x,s)/gen_harmonic(N,s)
  }
}
rzipf <- function(n,s,N){
  # Using Rejection method: x ~ U(1,100)
  ans <- NULL
  while(length(ans) < n){
    M <- N/gen_harmonic(N,s)
    x <- sample(1:N,1)
    r <- dzipf(x,N,s)/(M*(1/N))
    if (rbinom(1,1,r)==1){
      ans <- c(ans,x)
    }
  }
  return(ans)
}


occ.matrix2dataframe <- function(Ai,Words,...){
  V <- length(Words)
  combs <- llply(1:(V-1), function(j) list(Words[j],Words[(j+1):V]))
  get_nij <- function(x, v1,v2) min(x[v1,v1],x[v2,v2])
  get_all_nij <- function(x){
    ans <- ldply(combs,function(xx) ldply(xx[[2]],function(j) c(get_nij(x,xx[[1]],j),x[xx[[1]],j])))
    colnames(ans) <- c("n","x")
    return(ans)
  }
  return(llply(Ai,get_all_nij,...))
}

get_s <- function(x) {
  m <- mean(x)
  N <- max(x)
  mm <- function(s) {
    obj <- gen_harmonic(N,s-1)/gen_harmonic(N,s) - m
    return(t(obj)%*%obj)
  }
  s <- optim(1,mm,method = 'Brent',lower = 1,upper = 30)$par
  return(s)
}

get_pvalue_chisq <- function(x, s = NULL){
  x <- x[x>0]
  t <- table(x)
  if (length(t)<2) {
    return(NA)
  } else {
    if (is.null(s)) {s <- get_s(x)}
    tt <- dzipf(as.numeric(names(t)),s,max(x))
    # tt <- dpois(as.numeric(names(t)),s)
    return(suppressWarnings(chisq.test(t,p = tt,rescale.p = TRUE)$p.value))
  }
}

get_pvalue_chisq_pois <- function(x){
 #  x <- x[x>0]
  t <- table(x)
  if (length(t)<2) {
    return(NA)
  } else {
    tt <- dpois(as.numeric(names(t)),mean(x))
    return(suppressWarnings(chisq.test(t,p = tt,rescale.p = TRUE)$p.value))
  }
}

get_zipfs_pvals <- function(AA, ... ){
  pvals <- na.omit(laply(AA, function(xx) get_pvalue_chisq(xx[lower.tri(xx)],...)))
  return(pvals)
}

get_zipfs_measure <- function(AA, ... ){
  pvals <- na.omit(laply(AA, function(xx) get_pvalue_chisq(xx[lower.tri(xx)],...)))
  ptest <- prop.test(sum(pvals > .15),length(pvals),alternative = "greater")
  return(unlist(ptest[c("estimate","conf.int","p.value")]))
}

thm1 <- function(x){
  alpha <- get_s(x)
  A = (alpha^2)*(alpha+2)/4
  N  = sum(x)
  (A*N/log(N))^(1/(alpha+2))
}

get_thm1 <- function(x) {
  xx <- x[lower.tri(x)]
  xx <- xx[xx>0]
  return(thm1(xx))
}

get_pois_measure <- function(l,BB = binom_input){
  edgeden <- laply(BB,function(x) x[l,]$n)
  tt <- table(edgeden)
  tt <- tt[tt > 1]
  ttt <- as.numeric(names(tt[order(tt,decreasing = TRUE)]))
  ttt <- ttt[ttt>0]
  get_gf_pois <- function(NN){
  ans <- laply(BB,function(xx) ifelse(xx[l,]$n==NN,xx[l,]$x,NA)) %>%
        na.omit() %>%
        get_pvalue_chisq_pois()
  return(ans)
  }
  m <- laply(ttt, get_gf_pois) %>%
        na.omit() 
  return(sum(m>.15)/length(ttt))
}


order_graphs <- function(g,words.order) g[words.order,words.order]
standardize <- function(u) u/max(u)
lowerize <- function(u) u[lower.tri(u)]
laplacian <- function(u) {
  diag_ans <- u%*%ones(nrow(u))
  ans <- -u
  diag(ans) <- diag_ans
  return(ans)
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  # alpha = 1
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
b.est <- function(x) estBetaParams(mean(x),var(x))

d.gamma <- function(x) exp(lgamma(sum(x)) - lgamma(x[2]))
diff.gamma <- function(x) {
  dd <- d.gamma(x)
  return(x[2] - dd^(1/x[1]))
}

 get.diff.gamma <- function(NN){ 
     if (is.matrix(NN)) {
      nn <- apply(NN,2,mean)
      } else {
      nn <- mean(NN)     
      }
      ppp <- 1/(1+nn)
      vals <- b.est(ppp) %>% unlist()
       return(diff.gamma(vals))
 }

z.gof <- function(N,S = 1000,s = nrow(N)*.8){
  # V <- N %>% ncol()
  nds <- llply(1:S,function(...) sample(1:length(N),s))
  zM <- laply(nds,function(args) get.diff.gamma(N[args])) %>% na.omit()
  qq <- quantile(zM,c(.025,.975))
  return(c(mean(zM),qq))
 } 

zM.gof <- function(N,S = 1000,s = nrow(N)*.8, pvalue = FALSE){
  # V <- N %>% ncol()
  nds <- llply(1:S,function(...) sample(1:nrow(N),s))
  zM <- laply(nds,function(args) get.diff.gamma(N[args,])) %>% na.omit()
  zMn <- na.omit(zM) 
  zMn <- zMn[!is.infinite(zMn)]
  qq <- quantile(zMn,c(.025,.975))
  if (pvalue) {
    # f <- ecdf(zM)
    # ans <- ifelse(mean(zM)>0,f(0),1 -f(0))
      ans <- 1-pnorm(abs(mean(zMn)/sqrt(var(zMn))))
    # ans <- list(p.value = ans,zM = zM)
    } else {ans <- c(mean(zMn),qq)}
  return(ans)
 } 

ones <- function(r,c = 1) matrix(rep(1,r*c),nrow = r,ncol = c)


logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}


gen.N <- function(u) {
  m <- length(u)
  ans <- llply(1:(m-1), function(i) llply((i+1):m, function(j) min(u[i],u[j])))
  return(unlist(ans))
}

binarize <- function(g,p){
  ans <- g*0
  lg <- g[lower.tri(g)]/gen.N(diag(g))
  lg[is.nan(lg)] <- 0
  ans[lower.tri(ans)] <- ifelse(lg>p,1,0)
  return(ans+t(ans))
}

toDesikan <- function(u,regions = c(2:4,6:35,1)) {
  rows <- c(regions,regions + 35)
  ans <- u[rows,rows]
  return(ans)
}

toBrainNetViewer <- function(u,file,regions = c(2:4,6:35,1)) {
  rows <- c(regions,regions + 35)
  ans <- u[rows,rows]
  dimnames(ans) <- NULL
  write.table(ans,file,col.names = FALSE,row.names = FALSE)
}


rcateg <- function(n,prob) rmultinom(n,1,prob) %>% apply(2,which.max)
toMultinomial <- function(x,H) (x==(1:H))*1
as.Rexpression <- function(x) eval(parse(text=x))

reorder_factor <- function(x,order) x = factor(x,levels(x)[order])

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

m_matrix <- function(x) add(x)/length(x)
weighted_m_matrix <- function(x) {
  wm <- llply(x, function(xx) xx*(1/nrow(xx)))
  return(add(wm))
}