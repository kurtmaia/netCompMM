source("https://raw.githubusercontent.com/kurtmaia/sna/master/sna_utils.R")

pzipfIV <- function(x, k0,s,g,a) 1- (1 + ((x-k0)/s)^(1/g))^(-a)

rzipfIV <- function(n,k0,s,g,a, sort = FALSE,...){
  u <- runif(n)
  k <- round(s*((1-u)^(-1/a)-1)^g + k0)
  if (sort) {k <- k %>% sort(decreasing = TRUE,...)}
return(k)
}

n.edges <- function(V) V*(V-1)/2
gen.M <- function(u, L = u[length(u)]) {
  ans <- laply(1:(length(u)-1), function(j) min(u[j],L))
  return(ans)
}
gen.N <- function(u) {
  m <- length(u)
  ans <- llply(1:(m-1), function(i) llply((i+1):m, function(j) min(u[i],u[j])))
  return(unlist(ans))
}
get.N <- function(U,m1,m2) c(gen.N(U[1:m1]),gen.N(U[(m1+1):(m1+m2)]))
gen.cooc <- function(n,p) rbinom(length(n),n,p)
gen.adj <- function(u,la,V,lower = FALSE,inv = FALSE) {
  ans <- matrix(0,V,V)
  if (lower) {
    ans[lower.tri(ans)] <- la
    } else {
      ans[upper.tri(ans)] <- la
    }
  ans <- ans + t(ans)
  diag(ans) <- u
  seq <- 1:nrow(ans)
  if (inv) {seq <- nrow(ans):1}
  names <- paste("w",seq,sep = '')
  dimnames(ans) <- list(names,names)
  return(ans)
}

rMzipfIV <- function(r,V,a,lambda){ 
   p <- rbeta(V,a,lambda)%>% sort()
   N <- llply(1:r, function(...) {
           return(rgeom(V,p))
    }) 
  return(N)
}

rZipfGraph <- function(r,V, a=1,lambda=5, p= NULL, m1 = V/2,m2 = V-m1,p1 = .95,p2 = .2,p3 = .01,...) {
  N <- rMzipfIV(r,V,a,lambda)
  minN <- llply(N,function(x) gen.N(x))
  a <-  add(minN) %>% edge2matrix() %>% melt()
  a <- a[melt(lower.tri(matrix(NA,V,V)))$value,]
  names(a) <- c("row","col","N")
  if (is.null(p)) {
    a$p <- p1
   a$p[a$row > m1 & a$col > m1] <- runif(sum(a$row > m1 & a$col > m1),0,p3)
   a$p[a$row > m1 & a$col <= m1] <- runif(sum(a$row > m1 & a$col <= m1),0,p2)
  } else {
    a$p <- p
  }
  X <- llply(minN, function(x) rbinom(V*(V-1)/2,x,a$p))
  gr <- llply(1:r, function(i) gen.adj(N[[i]],X[[i]],V,lower = TRUE,...))
 return(list(graphs = gr,struc = a))
}


gen_words_seq <- function(words,seq){
  ans <- c(words[seq],sample(words[-seq],length(words)-length(seq)))
  return(ans)
}

single.rearrange <- function(net,words) {
  dimnames(net) <- list(words,words)
  return(net)
}
rearrange <- function(graphs,w){
  ans <- llply(graphs, function(x) single.rearrange(x,w))
  return(ans)
}

# General synthetic data setup
V <- 100  # number of nodes
T <- 4    # number of time points
n <- 50   # number of graphs
nG <- 2   # number of groups

# Sample graphs
## Sample a random matrix where cells are indep
rMat <- matrix(rbeta(V*V,1,15),V,V)
## Include dependencies
rMat[1:20,1:20] <- runif(400,.7,1)
## Generate the edge probability matrix 
lowerTri <- melt(lower.tri(rMat))$value
p0 <- melt(rMat)[lowerTri,3]
## Sample nG*T*n graphs of size V
gr <- rZipfGraph(nG*T*n,V,a  = 1, lambda = 10,p = p0)
ggplot(gr$struc,aes(x = row, y = col, fill=p)) + geom_tile()  + scale_fill_gradient2(limit=c(0,1),low="white", high = "blue")
words <- rownames(gr$graphs[[1]])

# Setting up group differences for 2 groups
seqs <- list(1:10,c(11:20,1:10),c(21:30,c(11:20,1:10)),c(31:40,c(21:30,c(11:20,1:10))))
group1 <- llply(1:T, function(i) rearrange(gr$graphs[(i-1)*n + 1:n],gen_words_seq(words,seqs[[i]]))) %>% unlist(recursive = FALSE)

seqs2 <- llply(seqs, function(x) x+30)
group2 <- llply(1:T, function(i) rearrange(gr$graphs[(i-1)*n + 1:n],gen_words_seq(words,seqs2[[i]]))) %>% unlist(recursive = FALSE)

graphs <- c(group1,group2)


side <-  data.frame(
           user = c(rep(paste("u",1:n,sep = ""),T),rep(paste("u",n+1:n,sep = ""),T)),
           time = rep(1:T,each = n),
           PopId = rep(1:2,each = T*n),
           stringsAsFactors = FALSE)
side$net.name <- paste(side$user,side$time,sep=".t")
names(graphs) <- side$net.name
graphs <- llply(graphs,function(x) x[words,words])
