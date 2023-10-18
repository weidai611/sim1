simrec <- function(N,
                   fu.min,
                   fu.max,
                   cens.prob = 0,
                   dist.x = "binomial",
                   par.x = 0,
                   beta.x = 0,
                   dist.z = "gamma",
                   par.z = 0,
                   dist.rec,
                   par.rec,
                   pfree = 0,
                   dfree = 0) {
  if ((cens.prob > 0) & (fu.min != fu.max)) {
    warning(
      paste0(
        "The censoring scheme for parameters cens.prob greater than 0 and fu.min != fu.max is undefined. \n",
        "The package properly implements two censoring schemes depending on parameters cens.prob, fu.min, and fu.max: \n",
        "a) cens.prob>0 and fu.min=fu.max: follow-up ends at time fu.max with a probability of 1-cens.prob and follow-up ends uniformly distributed in [0, fu.max] with a probability of cens.prob. \n",
        "b) cens.prob=0 and fu.min<fu.max: follow-up ends uniformly distributed in [0,fu.max] for each subject."
      )
    )
  }
  
  ID <- c(1:N)
  
  # generating the follow-up
  # follow-up uniformly distributed in [fu.min, fu.max] if not censored
  # or uniformly distributed in [0, fu.max] if censored
  if (cens.prob < 0 || cens.prob > 1) {
    stop("cens.prob must be a probability between 0 and 1")
  }
  if (fu.min > fu.max || fu.min < 0) {
    stop("fu.min must be a non-negative value smaller or equal fu.max")
  }
  fu <- rbinom(N, 1, cens.prob) # 1 = censored
  nr.cens <- sum(fu)
  if (nr.cens == 0) { # nobody censored
    fu <- runif(N, min = fu.min, max = fu.max)
  } else {
    index.cens <- which(fu == 1)
    fu[-index.cens] <- runif((N - nr.cens), min = fu.min, max = fu.max)
    fu[index.cens] <- runif(nr.cens, min = 0, max = fu.max)
  }
  
  if (length(beta.x) != length(dist.x)) {
    stop("dimensions of beta.x and dist.x differ")
  }
  if (length(beta.x) != length(par.x)) {
    stop("dimensions of beta.x and par.x differ")
  }
  
  # generating the covariate-matrix x
  nr.cov <- length(beta.x) # number of covariates
  x <- matrix(0, N, nr.cov) # matrix with N lines and one column for each covariate
  for (i in 1:nr.cov) {
    dist.x[i] <- match.arg(dist.x[i], choices = c("binomial", "normal"))
    if (dist.x[i] == "binomial") {
      if (length(par.x[[i]]) != 1) {
        stop("par.x has wrong dimension")
      }
      if (par.x[[i]] < 0 || par.x[[i]] > 1) {
        stop("par.x must be a probability between 0 and 1 for the binomial distributed covariate")
      }
      x[, i] <- c(rbinom(N, 1, par.x[[i]]))
    } else { # normally distributed covariate
      if (length(par.x[[i]]) != 2) {
        stop("par.x has wrong dimension")
      }
      mu.x <- par.x[[i]][1]
      sigma.x <- par.x[[i]][2]
      x[, i] <- c(rnorm(N, mean = mu.x, sd = sigma.x))
    }
  }
  
  # generating the frailty variable z
  z <- rep(1, N)
  dist.z <- match.arg(dist.z, choices = c("gamma", "lognormal"))
  if (length(par.z) != 1) {
    stop("par.z has wrong dimension")
  }
  if (par.z < 0) {
    stop("par.z must be non-negative")
  }
  if (par.z != 0) { # if par.z=0 then frailty=1 for all
    if (dist.z == "gamma") { # gamma-frailty
      aGamma <- 1 / par.z
      z <- rgamma(N, shape = aGamma, scale = 1 / aGamma)
    } else { # lognormal frailty
      mu.z <- log(1 / sqrt(par.z + 1))
      sigma.z <- sqrt(log(par.z + 1))
      z <- exp(rnorm(N, mean = mu.z, sd = sigma.z))
    }
  }
  
  # derivation of the distributional parameters for the recurrent event data
  dist.rec <- match.arg(dist.rec, choices = c("weibull", "lognormal", "gompertz", "step"))
  if (dist.rec == "lognormal") { # lognormal
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    }
    mu <- par.rec[1]
    sigma <- par.rec[2]
    if (any(beta.x != 0)) {
      warning("lognormal together with covariates specified: this does not define the usual lognormal model! see help for details")
    }
  } else if (dist.rec == "weibull") { # weibull
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    }
    lambda <- par.rec[1]
    nu <- par.rec[2]
  } else if (dist.rec == "gompertz") {
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    } # gompertz
    lambdag <- par.rec[1]
    alpha <- par.rec[2]
  } else if (dist.rec == "step") { # step
    if (length(par.rec) != 3) {
      stop("par.rec has wrong dimensions")
    }
    fc <- par.rec[1]
    sc <- par.rec[2]
    jump <- par.rec[3]
    jumpinv <- jump * fc
  }
  
  if (length(pfree) != 1) {
    stop("pfree has wrong dimension")
  }
  if (length(dfree) != 1) {
    stop("dfree has wrong dimension")
  }
  if (pfree < 0 || pfree > 1) {
    stop("pfree must be a probability between 0 and 1")
  }
  
  # initial step: simulation of N first event times
  U <- runif(N)
  Y <- (-1) * log(U) * exp((-1) * x %*% beta.x) * 1 / z
  if (dist.rec == "lognormal") { # lognormal
    t <- exp(qnorm(1 - exp((-1) * Y)) * sigma + mu)
  } else if (dist.rec == "weibull") { # weibull
    t <- ((lambda)^(-1) * Y)^(1 / nu)
  } else if (dist.rec == "gompertz") { # gompertz
    t <- (1 / alpha) * log((alpha / lambdag) * Y + 1)
  } else if (dist.rec == "step") { # step
    t <- rep(NA, N)
    indexTr1 <- which(Y <= jumpinv)
    if (length(indexTr1 > 0)) {
      t[indexTr1] <- Y[indexTr1] / fc
    }
    indexTr2 <- which(Y > jumpinv)
    if (length(indexTr2 > 0)) {
      t[indexTr2] <- (Y[indexTr2] - (fc - sc) * jump) / sc
    }
  }
  T <- matrix(t, N, 1)
  dirty <- rep(TRUE, N)
  T1 <- NULL
  
  # recursive step: simulation of N subsequent event times
  while (any(dirty)) {
    pd <- rbinom(N, 1, pfree)
    U <- runif(N)
    Y <- (-1) * log(U) * exp((-1) * x %*% beta.x) * 1 / z
    t1 <- t + pd * dfree
    if (dist.rec == "lognormal") { # lognormal
      t <- (t1 + exp(qnorm(1 - exp(log(1 - pnorm((log(t1) - mu) / sigma)) - Y)) * sigma + mu) - (t1))
    } else if (dist.rec == "weibull") { # weibul.l
      t <- (t1 + ((Y + lambda * (t1)^(nu)) / lambda)^(1 / nu) - (t1))
    } else if (dist.rec == "gompertz") { # gompertz
      t <- (t1 + ((1 / alpha) * log((alpha / lambdag) * Y + exp(alpha * t1))) - (t1))
    } else if (dist.rec == "step") { # step
      indexTr3 <- which((t1 <= jump) & (Y <= (jump - t1) * fc))
      if (length(indexTr3 > 0)) {
        t[indexTr3] <- t1[indexTr3] + Y[indexTr3] / fc
      }
      indexTr4 <- which((t1 <= jump) & (Y > (jump - t1) * fc))
      if (length(indexTr4 > 0)) {
        t[indexTr4] <- t1[indexTr4] + (Y[indexTr4] + (fc - sc) * (t1[indexTr4] - jump)) / sc
      }
      indexTr5 <- which(t1 > jump)
      if (length(indexTr5 > 0)) {
        t[indexTr5] <- t1[indexTr5] + Y[indexTr5] / sc
      }
    }
    T1 <- cbind(T1, ifelse(dirty, t1, NA))
    dirty <- ifelse(dirty, (t(t) < fu) & (t(t1) < fu), dirty)
    if (!any(dirty)) break
    T <- cbind(T, ifelse(dirty, t, NA))
  }
  
  # start times
  start.t <- cbind(0, T1)
  start.t <- as.vector(t(start.t))
  tab.start.t <- start.t[!is.na(start.t)]
  
  # stop times
  stop.t <- cbind(T, NA)
  d <- apply(!is.na(T), 1, sum) # number of events per individual
  f <- d + 1
  for (i in 1:N) {
    stop.t[i, f[i]] <- fu[i]
  }
  stop.t <- as.vector(t(stop.t))
  tab.stop.t <- stop.t[!is.na(stop.t)]
  
  # deriving the censoring indicator variable and truncating stop times that are larger than FU
  e <- NULL
  for (i in 1:N) {
    e <- cbind(e, t(rep(1, d[i])), 0)
  }
  
  tab.ID <- rep(ID, f)
  tab.X <- x[rep(1:nrow(x), f), ]
  tab.Z <- rep(z, f)
  tab.Fu <- rep(fu, f)
  
  w <- tab.start.t > tab.stop.t
  v <- rep(0, length(w))
  for (i in 1:length(w)) {
    if (w[i]) {
      v[i - 1] <- 1
    }
  }
  
  l <- tab.stop.t > tab.Fu
  for (i in 1:length(l)) {
    if (l[i]) {
      tab.stop.t[i] <- tab.Fu[i]
      e[i] <- 0
    }
  }
  
  tab <- cbind(tab.ID, tab.X, tab.Z, tab.start.t, tab.stop.t, t(e), tab.Fu)
  for (i in 1:length(w)) {
    if (w[i]) {
      tab[i, ] <- rep(NA, nr.cov + 6)
    }
  }
  tab <- data.frame(id = tab[, 1], x = tab[, 2:(nr.cov + 1)], z = tab[, (nr.cov + 2)], start = tab[, (nr.cov + 3)], stop = tab[, (nr.cov + 4)], status = tab[, (nr.cov + 5)], fu = tab[, (nr.cov + 6)])
  tab <- na.omit(tab)
  return(tab)
}

####

simrecPlot <- function(data,
                       id = "id",
                       start = "start",
                       stop = "stop",
                       status = "status") {
  if (!(id %in% colnames(data))) {
    stop("Please give the name of the id-column")
  }
  if (!(start %in% colnames(data))) {
    stop("Please give the name of the start-column")
  }
  if (!(stop %in% colnames(data))) {
    stop("Please give the name of the stop-column")
  }
  if (!(status %in% colnames(data))) {
    stop("Please give the name of the status-column")
  }
  
  colnames(data)[colnames(data) == id] <- "id"
  colnames(data)[colnames(data) == start] <- "start"
  colnames(data)[colnames(data) == stop] <- "stop"
  colnames(data)[colnames(data) == status] <- "status"
  
  data <- data[order(data$id), ] # data ordered by id
  t <- table(data$id) # the table entries will also be ordered by id
  idvec <- names(t) # all occuring IDs just one time
  idnum <- seq(along = idvec) # number the patients consecutively (corresponding to the ordering by id)
  idtable <- data.frame(idvec, idnum)
  data$idnum <- NULL
  for (i in idvec) {
    data$idnum[data$id == i] <- idtable$idnum[idtable$idvec == i] # new column with a numerical id for each patient
  }
  
  events <- data$stop[data$status == 1] # all event times in one vector
  cens <- data$stop[data$status == 0] # all censoring times in one vector
  
  idevents <- data$idnum[data$status == 1] # all numerical patient-IDs corresponding to the event times in one vector
  idcens <- data$idnum[data$status == 0] # all numerical patient-IDs corresponding to the censoring times in one vector
  
  nevents <- length(events) # how many event times    ...  if no events: events = numeric(0) and length(events)=0
  ncens <- length(cens) # how many censoring times
  
  opar <- par(las = 1, cex.axis = 0.5)
  on.exit(par(opar))
  
  plot(c(events, cens), c(idevents, idcens), # plot time points vs. corresponding patient-id
       pch = c(rep(20, nevents), rep(1, ncens)), # symbol for event: filled circle / symbol for censoring: circle
       xlab = "time", ylab = "machines", yaxt = "n", # axes=FALSE,
       main = "event history of machines in standard model "
  )
  for (i in idvec) {
    datai <- subset(data, id == i)
    for (j in seq_along(datai$id)) {
      lines(c(datai$start[j], datai$stop[j]), c(datai$idnum[j], datai$idnum[j]), lty = "solid") # add lines for each start-stop intervall
    }
  }
  # axis(1, labels=TRUE, at=0:max(data$fu), tick=TRUE)
  axis(2, labels = idvec, at = idnum, tick = TRUE)
}


###

simrecPlot(simdata_ga_we)





N <- 500
#'
#' ### with a binomially distributed covariate with a regression coefficient
#' ### of beta=0.3, and a standard normally distributed covariate with a
#' ### regression coefficient of beta=0.2,
#'
dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
par.x <- list(0.5, 0.5,0.5,0.5,0.5)
beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
#'
#' ### a gamma distributed frailty variable with variance 0.25
#'
 dist.z <- "gamma"
 par.z <- 0.2
#'
#' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
#' ### and scale parameter nu=2.
#'
 dist.rec <- "weibull"
 par.rec <- c(1, 2)
#'
#' ### Subjects are to be followed for two years with 20% of the subjects
#' ### being censored according to a uniformly distributed censoring time
#' ### within [0,2] (in years).
#'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
#'
#' ### After each event a subject is not at risk for experiencing further events
#' ### for a period of 30 days with a probability of 50%.
#'
 dfree <- 30 / 365
pfree <- 0.5
#'
 simdata_ga_we <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_ga_we)  # only run for small N!


 
 
 ##
 N <- 500
 #'
 #' ### with a binomially distributed covariate with a regression coefficient
 #' ### of beta=0.3, and a standard normally distributed covariate with a
 #' ### regression coefficient of beta=0.2,
 #'
 dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
 par.x <- list(0.5, 0.5,0.5,0.5,0.5)
 beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
 #'
 #' ### a gamma distributed frailty variable with variance 0.25
 #'
 dist.z <- "lognormal"
 par.z <- 0.2
 #'
 #' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
 #' ### and scale parameter nu=2.
 #'
 dist.rec <- "weibull"
 par.rec <- c(1, 2)
 #'
 #' ### Subjects are to be followed for two years with 20% of the subjects
 #' ### being censored according to a uniformly distributed censoring time
 #' ### within [0,2] (in years).
 #'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
 #'
 #' ### After each event a subject is not at risk for experiencing further events
 #' ### for a period of 30 days with a probability of 50%.
 #'
 dfree <- 30 / 365
 pfree <- 0.5
 #'
 simdata_log_we <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_log_we)  # only run for small N!
 
 
 
 
 
 ##
 N <- 500
 #'
 #' ### with a binomially distributed covariate with a regression coefficient
 #' ### of beta=0.3, and a standard normally distributed covariate with a
 #' ### regression coefficient of beta=0.2,
 #'
 dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
 par.x <- list(0.5, 0.5,0.5,0.5,0.5)
 beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
 #'
 #' ### a gamma distributed frailty variable with variance 0.25
 #'
 dist.z <- "lognormal"
 par.z <-0.2
 #'
 #' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
 #' ### and scale parameter nu=2.
 #'
 dist.rec <- "lognormal"
 par.rec <- c(0, 1)
 #'
 #' ### Subjects are to be followed for two years with 20% of the subjects
 #' ### being censored according to a uniformly distributed censoring time
 #' ### within [0,2] (in years).
 #'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
 #'
 #' ### After each event a subject is not at risk for experiencing further events
 #' ### for a period of 30 days with a probability of 50%.
 #'
 dfree <- 30 / 365
 pfree <- 0.5
 #'
 simdata_log_log <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_log_log)  # only run for small N!
 
 
 
 ##
 N <- 500
 #'
 #' ### with a binomially distributed covariate with a regression coefficient
 #' ### of beta=0.3, and a standard normally distributed covariate with a
 #' ### regression coefficient of beta=0.2,
 #'
 dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
 par.x <- list(0.5, 0.5,0.5,0.5,0.5)
 beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
 #'
 #' ### a gamma distributed frailty variable with variance 0.25
 #'
 dist.z <- "lognormal"
 par.z <- 0.2
 #'
 #' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
 #' ### and scale parameter nu=2.
 #'
 dist.rec <- "gompertz"
 par.rec <- c(1, 2)
 #'
 #' ### Subjects are to be followed for two years with 20% of the subjects
 #' ### being censored according to a uniformly distributed censoring time
 #' ### within [0,2] (in years).
 #'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
 #'
 #' ### After each event a subject is not at risk for experiencing further events
 #' ### for a period of 30 days with a probability of 50%.
 #'
 dfree <- 30 / 365
 pfree <- 0.5
 #'
 simdata_log_go <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_log_go)  # only run for small N!
 
 
 
 N <- 500
 #'
 #' ### with a binomially distributed covariate with a regression coefficient
 #' ### of beta=0.3, and a standard normally distributed covariate with a
 #' ### regression coefficient of beta=0.2,
 #'
 dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
 par.x <- list(0.5, 0.5,0.5,0.5,0.5)
 beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
 
 #' ### a gamma distributed frailty variable with variance 0.25
 #'
 dist.z <- "gamma"
 par.z <- 0.2
 #'
 #' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
 #' ### and scale parameter nu=2.
 #'
 dist.rec <- "gompertz"
 par.rec <- c(1, 2)
 #'
 #' ### Subjects are to be followed for two years with 20% of the subjects
 #' ### being censored according to a uniformly distributed censoring time
 #' ### within [0,2] (in years).
 #'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
 #'
 #' ### After each event a subject is not at risk for experiencing further events
 #' ### for a period of 30 days with a probability of 50%.
 #'
 dfree <- 30 / 365
 pfree <- 0.5
 #'
 simdata_ga_go <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_ga_go)  # only run for small N!
 
 
 
 N <- 500
 #'
 #' ### with a binomially distributed covariate with a regression coefficient
 #' ### of beta=0.3, and a standard normally distributed covariate with a
 #' ### regression coefficient of beta=0.2,
 #'
 dist.x <- c("binomial", "binomial","binomial", "binomial","binomial")
 par.x <- list(0.5, 0.5,0.5,0.5,0.5)
 beta.x <- c(0.3, 0.2,0.1,0.4,0.25)
 #'
 #' ### a gamma distributed frailty variable with variance 0.25
 #'
 dist.z <- "gamma"
 par.z <- 0.2
 #'
 #' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
 #' ### and scale parameter nu=2.
 #'
 dist.rec <- "lognormal"
 par.rec <- c(0, 1)
 #'
 #' ### Subjects are to be followed for two years with 20% of the subjects
 #' ### being censored according to a uniformly distributed censoring time
 #' ### within [0,2] (in years).
 #'
 fu.min <- 2
 fu.max <- 2
 cens.prob <- 0.2
 #'
 #' ### After each event a subject is not at risk for experiencing further events
 #' ### for a period of 30 days with a probability of 50%.
 #'
 dfree <- 30 / 365
 pfree <- 0.5
 #'
 simdata_ga_log <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
   dist.rec, par.rec, pfree, dfree )
 print(simdata_ga_log)  # only run for small N!
 
 
 
 
 
 
 
 
 
 
 
 


 




  N <- 1000
  dist.x <- c("binomial", "normal")
  par.x <- list(0.5, c(0, 1))
 beta.x <- c(0.3, 0.2)
 dist.z <- "gamma"
  par.z <- 0.25
 dist.rec <- "weibull"
  par.rec <- c(1, 2)
 fu.min <- 2
 fu.max <- 2
  cens.prob <- 0.2
  dfree <- 30 / 365
 pfree <- 0.5
 simdata_ga_we <- simrec(
   N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x,
  dist.z, par.z, dist.rec, par.rec, pfree, dfree
  )
 print(simdata_ga_we)

 
 ####
simdata_ga_we<-as.data.frame(simdata_ga_we)
ag_ga_we<- aggregate(simdata_ga_we$status,by=list(id=simdata_ga_we$id),FUN = sum)
f_ga_we<-aggregate(simdata_ga_we$z,by=list(id=simdata_ga_we$id),FUN = mean)

plot(ag_ga_we)
plot(f_ga_we)
mean(f_ga_we$x)#0.9985773
var(f_ga_we$x)#0.2336942
max(ag_ga_we$x)  #16
min(ag_ga_we$x)#0
mean(ag_ga_we$x) #3.546
sum(ag_ga_we$x>1)#733   intotal 1000n, 73.3%recurrent



#
simdata_ga_go<-as.data.frame(simdata_ga_go)
ag_ga_go<- aggregate(simdata_ga_go$status,by=list(id=simdata_ga_go$id),FUN = sum)
f_ga_go<-aggregate(simdata_ga_go$z,by=list(id=simdata_ga_go$id),FUN = mean)

plot(ag_ga_go)
plot(f_ga_go)
mean(f_ga_go$x)#0.9940622
var(f_ga_go$x)#0.2537558
max(ag_ga_go$x)  #34
min(ag_ga_go$x)#0
mean(ag_ga_go$x) #12.494
sum(ag_ga_go$x>1)#918   intotal 1000n, 91.8%recurrent


#
simdata_ga_log<-as.data.frame(simdata_ga_log)
ag_ga_log<- aggregate(simdata_ga_log$status,by=list(id=simdata_ga_log$id),FUN = sum)
f_ga_log<-aggregate(simdata_ga_log$z,by=list(id=simdata_ga_log$id),FUN = mean)

plot(ag_ga_log)
plot(f_ga_log)
mean(f_ga_log$x)#0.990703
var(f_ga_log$x)#0.2558677
max(ag_ga_log$x)  #6
min(ag_ga_log$x)#0
mean(ag_ga_log$x) #0.586
sum(ag_ga_log$x>1)#127   intotal 1000n, 12.7%recurrent


#
simdata_log_go<-as.data.frame(simdata_log_go)
ag_log_go<- aggregate(simdata_log_go$status,by=list(id=simdata_log_go$id),FUN = sum)
f_log_go<-aggregate(simdata_log_go$z,by=list(id=simdata_log_go$id),FUN = mean)

plot(ag_log_go)
plot(f_log_go)
mean(f_log_go$x)#1.007694
var(f_log_go$x)#0.2785436
max(ag_log_go$x)  #36
min(ag_log_go$x)#0
mean(ag_log_go$x) #12.525
sum(ag_log_go$x>1)#919   intotal 1000n, 91.9%recurrent


#
simdata_log_we<-as.data.frame(simdata_log_we)
ag_log_we<- aggregate(simdata_log_we$status,by=list(id=simdata_log_we$id),FUN = sum)
f_log_we<-aggregate(simdata_log_we$z,by=list(id=simdata_log_we$id),FUN = mean)

plot(ag_log_we)
plot(f_log_we)
mean(f_log_we$x)#1.006624
var(f_log_we$x)#0.2549697
max(ag_log_we$x)  #23
min(ag_log_we$x)#0
mean(ag_log_we$x) #3.621
sum(ag_log_we$x>1)#764 intotal 1000n, 12.7%recurrent



#
simdata_log_log<-as.data.frame(simdata_log_log)
ag_log_log<- aggregate(simdata_log_log$status,by=list(id=simdata_log_log$id),FUN = sum)
f_log_log<-aggregate(simdata_log_log$z,by=list(id=simdata_log_log$id),FUN = mean)

plot(ag_log_log)
plot(f_log_log)
mean(f_log_log$x)#1.001378
var(f_log_log$x)#0.2463812
max(ag_log_log$x)  #5
min(ag_log_log$x)#0
mean(ag_log_log$x) #0.643
sum(ag_log_log$x>1)#147 intotal 1000n, 14.7%recurrent






write.csv(x=simdata_ga_go,file="simdata_ga_go.csv")




































x<-c(1,2,3)
y<-c(0.0037558,0.0058677,0.163058)
plot(x,y,col="blue",main="Standard frailty model with lognormal frailty VS gamma frailty")


x2<-c(1,2,3)
y2<-c(0.0285436,0.006188,0.0049697)
lines(x2,y2,col="red")





x1<-c("Gompertz","Lognormal","Weibull")
y1<-c(0.0159953,0.0044328,0.0038424)
plot(x,y,xlab="Baseline function",ylB="
     distance",col="blue")



