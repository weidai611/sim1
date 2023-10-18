
simreccompPlot <- function(data,
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
  compev <- data$stop[data$status == 2] # all competing event times in one vector
  
  idevents <- data$idnum[data$status == 1] # all numerical patient-IDs corresponding to the event times in one vector
  idcens <- data$idnum[data$status == 0] # all numerical patient-IDs corresponding to the censoring times in one vector
  idcompev <- data$idnum[data$status == 2] # all numerical patient-IDs corresponding to the competing event times in one vector
  
  nevents <- length(events) # how many event times                if no events: events = numeric(0) and length(events)=0
  ncens <- length(cens) # how many censoring times
  ncompev <- length(compev) # how many competing events
  
  opar <- par(las = 1, cex.axis = 0.5)
  on.exit(par(opar))
  
  plot(c(events, cens, compev), c(idevents, idcens, idcompev), # plot time points vs. corresponding patient-id
       pch = c(rep(20, nevents), rep(1, ncens), rep(4, ncompev)), # symbol for event: filled circle / symbol for censoring: circle / symbol for comp. event: x
       xlab = "time", ylab = "machines", yaxt = "n", # axes=FALSE,
       main = "event history of machines in joint model"
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

simreccompPlot(simdata1)
simreccompPlot(simdata2)



simreccomp <- function(N,
                       fu.min,
                       fu.max,
                       cens.prob = 0,
                       dist.x = "binomial",
                       par.x = 0,
                       beta.xr = 0,
                       beta.xc = 0,
                       dist.zr = "gamma",
                       par.zr = 0,
                       a = NULL,
                       dist.zc = NULL,
                       par.zc = NULL,
                       dist.rec,
                       par.rec,
                       dist.comp,
                       par.comp,
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
  # generating the follow-up  *****************************************************************
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
  if (length(beta.xr) != length(dist.x)) {
    stop("dimensions of beta.xr and dist.x differ")
  }
  if (length(beta.xr) != length(par.x)) {
    stop("dimensions of beta.xr and par.x differ")
  }
  if (length(beta.xr) != length(beta.xc)) {
    stop("dimensions of beta.xr and beta.xc differ")
  }
  
  # generating the covariate-matrix x   *****************************************************
  nr.cov <- length(beta.xr) # number of covariates
  x <- matrix(0, N, nr.cov) # matrix with N lines and one column for each covariate
  
  for (i in 1:nr.cov) {
    dist.x[i] <- match.arg(dist.x[i], choices = c("binomial", "normal"))
    if (dist.x[i] == "binomial") {
      if (length(par.x[[i]]) != 1) {
        stop("par.x has wrong dimension")
      }
      if (par.x[[i]] < 0 || par.x[[i]] > 1) {
        stop("par.x must be a probability between 0 and 1 for the binomially distributed covariate")
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
  
  # generating the frailty variables zr and zc  ************************************************************
  if (length(a) != 0 & (length(dist.zc) != 0 | length(par.zc) != 0)) {
    stop("enter either a or dist.zc and par.zc")
  }
  zr <- rep(1, N)
  dist.zr <- match.arg(dist.zr, choices = c("gamma", "lognormal"))
  if (length(par.zr) != 1) {
    stop("par.zr has wrong dimension")
  }
  if (par.zr < 0) {
    stop("par.zr must be non-negative")
  }
  if (par.zr != 0) { # if par.zr=0 then frailty=1 for all
    if (dist.zr == "gamma") { # gamma-frailty
      aGamma.r <- 1 / par.zr
      zr <- rgamma(N, shape = aGamma.r, scale = 1 / aGamma.r)
    } else { # lognormal frailty
      mu.zr <- log(1 / sqrt(par.zr + 1))
      sigma.zr <- sqrt(log(par.zr + 1))
      zr <- exp(rnorm(N, mean = mu.zr, sd = sigma.zr))
    }
  }
  
  if (length(a) == 0) {
    if (length(dist.zc) == 0 | length(par.zc) == 0) {
      stop("enter either a or dist.zc and par.zc")
    }
    zc <- rep(1, N)
    dist.zc <- match.arg(dist.zc, choices = c("gamma", "lognormal"))
    if (length(par.zc) != 1) {
      stop("par.zc has wrong dimension")
    }
    if (par.zc < 0) {
      stop("par.zc must be non-negative")
    }
    if (par.zc != 0) { # if par.zc=0 then frailty=1 for all
      if (dist.zc == "gamma") { # gamma-frailty
        aGamma.c <- 1 / par.zc
        zc <- rgamma(N, shape = aGamma.c, scale = 1 / aGamma.c)
      } else { # lognormal frailty
        mu.zc <- log(1 / sqrt(par.zc + 1))
        sigma.zc <- sqrt(log(par.zc + 1))
        zc <- exp(rnorm(N, mean = mu.zc, sd = sigma.zc))
      }
    }
  } else {
    zc <- zr**a
  }
  # generating the recurrent event times *************************************************************
  # derivation of the distributional parameters for the recurrent event data
  dist.rec <- match.arg(dist.rec, choices = c("weibull", "lognormal", "gompertz", "step"))
  if (dist.rec == "lognormal") { # lognormal
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    }
    mu <- par.rec[1]
    sigma <- par.rec[2]
    if (any(beta.xr != 0)) {
      warning("lognormal together with covariates specified: this does not define the usual lognormal model! see help for details")
    }
  } else if (dist.rec == "weibull") { # weibull
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    }
    lambda <- par.rec[1]
    nu <- par.rec[2]
  } else if (dist.rec == "gompertz") { # gompertz
    if (length(par.rec) != 2) {
      stop("par.rec has wrong dimension")
    }
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
  Y <- (-1) * log(U) * exp((-1) * x %*% beta.xr) * 1 / zr
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
    Y <- (-1) * log(U) * exp((-1) * x %*% beta.xr) * 1 / zr
    t1 <- t + pd * dfree
    if (dist.rec == "lognormal") { # lognormal
      t <- (t1 + exp(qnorm(1 - exp(log(1 - pnorm((log(t1) - mu) / sigma)) - Y)) * sigma + mu) - (t1))
    } else if (dist.rec == "weibull") { # weibull
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
  
  # comp. events simulation ***********************************************************
  # derivation of the distributional parameters for the comp. events
  dist.comp <- match.arg(dist.comp, choices = c("weibull", "lognormal", "gompertz", "step"))
  
  if (dist.comp == "lognormal") { # lognormal
    if (length(par.comp) != 2) {
      stop("par.comp has wrong dimension")
    }
    mu2 <- par.comp[1]
    sigma2 <- par.comp[2]
  } else if (dist.comp == "weibull") { # weibull
    if (length(par.comp) != 2) {
      stop("par.comp has wrong dimension")
    }
    lambda2 <- par.comp[1]
    nu2 <- par.comp[2]
  } else if (dist.comp == "gompertz") { # gompertz
    if (length(par.comp) != 2) {
      stop("par.comp has wrong dimension")
    }
    lambdag2 <- par.comp[1]
    alpha2 <- par.comp[2]
  } else if (dist.comp == "step") { # step
    if (length(par.comp) != 3) {
      stop("par.comp has wrong dimensions")
    }
    fc2 <- par.comp[1]
    sc2 <- par.comp[2]
    jump2 <- par.comp[3]
    jumpinv2 <- jump2 * fc2
  }
  
  # simulation of N comp. events
  U2 <- runif(N)
  Y2 <- (-1) * log(U2) * exp((-1) * x %*% beta.xc) * 1 / zc
  if (dist.comp == "lognormal") { # lognormal
    t2 <- exp(qnorm(1 - exp((-1) * Y2)) * sigma2 + mu2)
  } else if (dist.comp == "weibull") { # weibull
    t2 <- ((lambda2)^(-1) * Y2)^(1 / nu2)
  } else if (dist.comp == "gompertz") { # gompertz
    t2 <- (1 / alpha2) * log((alpha2 / lambdag2) * Y2 + 1)
  } else if (dist.comp == "step") { # step
    t2 <- rep(NA, N)
    indexTr12 <- which(Y2 <= jumpinv2)
    if (length(indexTr12 > 0)) {
      t2[indexTr12] <- Y2[indexTr12] / fc2
    }
    indexTr22 <- which(Y2 > jumpinv2)
    if (length(indexTr22 > 0)) {
      t2[indexTr22] <- (Y2[indexTr22] - (fc2 - sc2) * jump2) / sc2
    }
  }
  T2 <- matrix(t2, N, 1)
  comp.event <- as.vector(t(T2))
  comp.event <- comp.event[!is.na(comp.event)]
  
  # **************************************************************************************
  
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
  tab.zr <- rep(zr, f)
  tab.zc <- rep(zc, f)
  tab.Fu <- rep(fu, f)
  tab.comp.event <- rep(comp.event, f)
  
  w <- tab.start.t > tab.stop.t
  # v <- rep(0,length(w))
  # for (i in 1:length(w)){
  #  if (w[i]) {v[i-1] <- 1}
  # }
  
  l <- tab.stop.t > tab.Fu
  for (i in 1:length(l)) {
    if (l[i]) {
      tab.stop.t[i] <- tab.Fu[i]
      e[i] <- 0
    }
  }
  
  s <- tab.start.t > tab.comp.event # create vector which remembers the times which are after the comp.event and therefore do not exist.
  # r<-rep(0,length(s))
  # for (i in 1:length(s)){
  #  if (s[i]) {r[i-1]<-1}
  # }
  
  m <- tab.stop.t > tab.comp.event # modify status vector, whenever the stop.time is greater than the comp. event it leads to status 2
  for (i in 1:length(m)) {
    if (m[i]) {
      tab.stop.t[i] <- tab.comp.event[i]
      e[i] <- 2
    }
  }
  
  tab <- cbind(
    tab.ID, tab.X, tab.zr, tab.zc, tab.start.t, tab.stop.t,
    t(e),
    pmin(tab.Fu, tab.comp.event),
    tab.comp.event
  )
  
  
  for (i in 1:length(w)) {
    if (w[i]) {
      tab[i, ] <- rep(NA, nr.cov + 8)
    } # delete times, which are after the FU and therefore don't exist
  }
  for (i in 1:length(w)) {
    if (s[i]) {
      tab[i, ] <- rep(NA, nr.cov + 8)
    } # delete times, which are after the comp. event and therefore don't exist
  }
  
  tab <- data.frame(id = tab[, 1], x = tab[, 2:(nr.cov + 1)], zr = tab[, (nr.cov + 2)], zc = tab[, (nr.cov + 3)], start = tab[, (nr.cov + 4)], stop = tab[, (nr.cov + 5)], status = tab[, (nr.cov + 6)], fu = tab[, (nr.cov + 7)])
  tab <- na.omit(tab)
  
  return(tab)
}











































N <- 300
#'
#' ### with a binomially distributed covariate and a standard normally distributed covariate
#' ### with regression coefficients of beta.xr=0.3 and beta.xr=0.2, respectively,
#' ### for the recurrent events,
#' ### as well as regression coefficients of beta.xc=0.5 and beta.xc=0.25, respectively,
#' ### for the competing event.
#'
 dist.x <- c("binomial", "binomial","binomial","binomial","binomial")
 par.x <- list(0.5, 0.5, 0.5,0.5,0.5)
beta.xr <- c(0.3, 0.2,0.1,0.4,0.25)
beta.xc <- c(0.5, 0.25, 0.15,0.65,0.4)
#'
#' ### a gamma distributed frailty variable for the recurrent event with variance 0.25
#' ### and for the competing event with variance 0.3,
#'
 dist.zr <- "gamma"
par.zr <- 0.2
#'
 dist.zc <- "gamma"
 par.zc <- 0.3
#'
#' ### alternatively the frailty variable for the competing event can be computed via a:
a <- 0.5
#'
#' ### Furthermore a Weibull-shaped baseline hazard for the recurrent event with shape parameter
#' ### lambda=1 and scale parameter nu=2,
#'
dist.rec <- "gompertz"
 par.rec <- c(1, 2)
#'
#' ### and a Weibull-shaped baseline hazard for the competing event with shape parameter lambda=1
#' ### and scale parameter nu=2
#
dist.comp <- "gompertz"
 par.comp <- c(1,2)
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
pfree <- 0.2
#'
simdata1 <- simreccomp(
   N = N, fu.min = fu.min, fu.max = fu.max, cens.prob = cens.prob,
   dist.x = dist.x, par.x = par.x, beta.xr = beta.xr, beta.xc = beta.xc,
   dist.zr = dist.zr, par.zr = par.zr, a = a,
   dist.rec = dist.rec, par.rec = par.rec, dist.comp = dist.comp, par.comp = par.comp,
   pfree = pfree, dfree = dfree
 )
#'
 simdata2 <- simreccomp(
   N = N, fu.min = fu.min, fu.max = fu.max, cens.prob = cens.prob,
   dist.x = dist.x, par.x = par.x, beta.xr = beta.xr, beta.xc = beta.xc,
   dist.zr = dist.zr, par.zr = par.zr, dist.zc = dist.zc, par.zc = par.zc,
   dist.rec = dist.rec, par.rec = par.rec, dist.comp = dist.comp, par.comp = par.comp,
   pfree = pfree, dfree = dfree
 )
#'
 simdata1
 simdata2
 
 #
 
 
d1<- read.csv("simdata1.csv",header = TRUE)
d2<- read.csv("simdata2.csv",header = TRUE)
 write.csv(x=simdata1,file="simdata1.csv")
 
 write.csv(x=simdata2,file="simdata2.csv")
 
 write.table(d1,"C:/Users/freedom1106/d1.txt")
 
 write.table(d2,"C:/Users/freedom1106/d2.txt")
 
data1<- read.table("C:/Users/freedom1106/data_failure_1.txt",header = TRUE)
 
write.csv(x=data1,file="data1.csv")
 
 
 
 
 #
 simdata1<-as.data.frame(simdata1)
 ag1<- aggregate(simdata1$status,by=list(id=simdata1$id),FUN = sum)
 fr<-aggregate(simdata1$zr,by=list(id=simdata1$id),FUN = mean)
 fc<-aggregate(simdata1$zc,by=list(id=simdata1$id),FUN = mean)
 
 mean(fr$x)#1.013466
 var(fr$x)# 0.2461576
 mean(fc$x)#0.9757493
 var(fc$x)#0.06144049

 mean(simdata1$status==1) #1.466546
 mean(simdata2$status==2)
 sum(simdata1$status==0)#99  intotal 1000n, 73.3%recurrent
 sum(simdata1$status==1) #687
 sum(simdata1$status==2)#873    1659 in total
 
 
 
 
 
 #
 simdata2<-as.data.frame(simdata2)
 ag2<- aggregate(simdata2$status,by=list(id=simdata2$id),FUN = sum)
 fr2<-aggregate(simdata2$zr,by=list(id=simdata2$id),FUN = mean)
 fc2<-aggregate(simdata2$zc,by=list(id=simdata2$id),FUN = mean)

 mean(fr2$x)# 1.020242
 var(fr2$x)#  0.2678889
 mean(fc2$x)#0.9657091
 var(fc2$x)#0.2888465
 
 
 
 
 mean(simdata2$status) #1.351452
 sum(simdata2$status==0)#133  
 sum(simdata2$status==1) #983
 sum(simdata2$status==2)#840
 
 
 
 
 
 
 
 
 
 