### function pffrSim from package refund with small changes
### see ?pffrSim
pffrSimVar <- function (scenario = "all", n = 100, nxgrid = 40, nygrid = 60, 
                        SNR = 10, propmissing = 0, limits = NULL, bd = 7,
                        simFun=function(s, t){s * cos(pi * abs(s - t)) - 0.19},
                        irrTP = FALSE) 
{
  mc <- match.call()
  for (i in 2:length(mc)) if (is.symbol(mc[[i]])) 
    mc[[i]] <- get(deparse(mc[[i]]), envir = parent.frame())
  rf <- function(x = seq(0, 1, length = 100), bs.dim = 7, center = FALSE) {
    nk <- bs.dim - 2
    xu <- max(x)
    xl <- min(x)
    xr <- xu - xl
    xl <- xl - xr * 0.001
    xu <- xu + xr * 0.001
    dx <- (xu - xl)/(nk - 1)
    kn <- seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 
                2)
    X <- splines::spline.des(kn, x, 4, x * 0)$design
    drop(X %*% rnorm(bs.dim))
  }
  test1 <- simFun
  if(irrTP){
    
    s <- (1:nxgrid-1)^2/(nxgrid-1)^2
    t <- (1:nygrid-1)^2/(nygrid-1)^2
    
  }else{
    
    s <- seq(0, 1, length = nxgrid)
    t <- seq(0, 1, length = nygrid)
    
  }
  mu.t <- matrix(1 + dbeta(t, 2, 7), nrow = n, ncol = nygrid, 
                 byrow = TRUE)
  data <- list()
  etaTerms <- list()
  etaTerms$int <- mu.t
  data$X1 <- I(t(replicate(n, rf(s, bs.dim = bd))))
  L <- if(irrTP) FDboost:::integrationWeightsLeft(X1=data$X1, xind=s) else 
    matrix(1/nxgrid, ncol = nxgrid, nrow = n)
  LX1 <- L * data$X1
  beta1.st <- outer(s, t, test1)
  if (!is.null(limits)) {
    range <- outer(s, t, limits)
    beta1.st <- beta1.st * range
  }
  etaTerms$X1 <- LX1 %*% beta1.st
  data$xlin <- I(rnorm(n))
  beta.t <- matrix(scale(-dnorm(4 * (t - 0.2))), nrow = n, 
                   ncol = nygrid, byrow = T)
  etaTerms$xlin <- data$xlin * beta.t
  data$xsmoo <- I(rnorm(n))
  etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t - 
                                                           0.5), "*")
  data$xte1 <- I(rnorm(n))
  data$xte2 <- I(rnorm(n))
  etaTerms$xte <- matrix(drop(scale(-data$xte1 * data$xte2^2)), 
                         ncol = nygrid, nrow = n)
  data$xconst <- I(rnorm(n))
  etaTerms$xconst <- matrix(2 * data$xconst, ncol = nygrid, 
                            nrow = n)
  if (length(scenario) == 1) {
    eta <- mu.t + switch(scenario, int = 0, all = Reduce("+", 
                                                         etaTerms), ff = etaTerms$X1, lin = etaTerms$xlin, 
                         smoo = etaTerms$xsmoo, te = etaTerms$xte, const = etaTerms$xconst)
  }
  else {
    stopifnot(all(scenario %in% c("int", "ff", "lin", "smoo", 
                                  "te", "const")))
    eta <- 0 * mu.t
    if ("int" %in% scenario) 
      eta <- eta + mu.t
    if ("ff" %in% scenario) 
      eta <- eta + etaTerms$X1
    if ("lin" %in% scenario) 
      eta <- eta + etaTerms$xlin
    if ("smoo" %in% scenario) 
      eta <- eta + etaTerms$xsmoo
    if ("te" %in% scenario) 
      eta <- eta + etaTerms$xte
    if ("const" %in% scenario) 
      eta <- eta + etaTerms$xconst
  }
  eps <- sd(as.vector(eta))/sqrt(SNR) * matrix(scale(rnorm(n * 
                                                             nygrid)), nrow = n)
  data$Y <- I(eta + eps)
  if (propmissing == 0) {
    return(structure(as.data.frame(data, rownames = 1:n), 
                     xindex = s, yindex = t, truth = list(eta = eta, etaTerms = etaTerms), 
                     call = mc))
  }
  else {
    missing <- sample(c(rep(T, propmissing * n * nygrid), 
                        rep(F, n * nygrid - propmissing * n * nygrid)))
    data <- as.data.frame(data, rownames = 1:n)
    ydata <- data.frame(.obs = rep(1:n, each = nygrid)[!missing], 
                        .index = rep(t, times = n)[!missing], .value = as.vector(t(data$Y))[!missing])
    return(structure(list(data = data, ydata = ydata), xindex = s, 
                     yindex = t, truth = list(eta = eta, etaTerms = etaTerms), 
                     call = mc))
  }
}