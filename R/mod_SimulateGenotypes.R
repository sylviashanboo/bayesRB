sg <- function (or, f, p, n = 1e+05, nog = 400, varyEffects = FALSE, 
                seed = NULL, silent = FALSE) 
{
  error = ""
  if (!varyEffects & length(or) > 1) {
    error = paste(error, "length(or) must ==1 if not varying effects\n", 
                  sep = "")
  }
  if (!varyEffects & length(f) > 1) {
    error = paste(error, "length(f) must ==1 if not varying effects\n", 
                  sep = "")
  }
  if (length(p) > 1) {
    error = paste(error, "length(p) must ==1\n", sep = "")
  }
  if (error != "") {
    cat(error)
    return
  }
  nold <- n
  if (n > 10000) 
    n = ceiling(n/1000) * 1000
  if (nold != n) 
    cat(paste("n adjusted to ", n, sep = ""))
  lNog = length(nog)
  nogV = nog
  if (lNog > 1) {
    nog = max(nog)
  }
  makeVV <- function(x) {
    xout <- x
    lx <- length(x)
    if (lx == 1) {
      xout <- rep(x, nog)
    }
    else if (lx < nog) {
      xout <- c(x, rep(x[lx], nog - lx))
    }
    else {
      xout <- x[1:nog]
    }
    xout
  }
  if (!varyEffects | length(or) == 1 & length(f) == 1) {
    varyEffects <- FALSE
  }
  else {
    or <- makeVV(or)
    f <- makeVV(f)
  }
  risk = dis = matrix(NA, ncol = lNog, nrow = n)
  sampleDivide <- c(1, 5, 10, 20, 100)
  iSS <- 1
  done <- FALSE
  loopCount <- 0
  while (!done) {
    loopCount <- loopCount + 1
    sampleSucceed <- FALSE
    while (!sampleSucceed) {
      if (loopCount == 1 & !is.null(seed)) {
        set.seed(seed)
      }
      subSize <- n/sampleDivide[iSS]
      if (!varyEffects | length(unique(f)) == 1) {
        genePool <- c(rep(1, n * f^2), rep(2, 2 * n * 
                                             f * (1 - f)), rep(3, n * (1 - f)^2))
        nid <- try(sample(genePool, nog * subSize, replace = TRUE), 
                   silent = silent)
        whichEE <- try(matrix(nid == 1, ncol = nog), 
                       silent = silent)
        whichEe <- try(matrix(nid == 2, ncol = nog), 
                       silent = silent)
        whichee <- try(matrix(nid == 3, ncol = nog), 
                       silent = silent)
      }
      else {
        sampleGenes <- function(ff, nn) {
          sapply(ff, function(x) sample(1:3, nn, replace = TRUE, 
                                        prob = c(x^2, 2 * x * (1 - x), (1 - x)^2)))
        }
        nid <- try(sampleGenes(f, subSize), silent = silent)
        whichEE <- NULL
        whichEe <- NULL
        whichee <- NULL
      }
      if (is(nid)[1] == "try-error" | is(whichEE)[1] == 
          "try-error" | is(whichEe)[1] == "try-error" | 
          is(whichee)[1] == "try-error") {
        if (loopCount > 1) {
          cat("Try again using a different memory setting\n")
          break
        }
        iSS <- iSS + 1
      }
      else {
        sampleSucceed <- TRUE
      }
    }
    if (loopCount == 1) {
      baselineLOdds <- list()
      if (!varyEffects | length(or) == 1 & length(f) == 
          1) {
        for (i in 1:lNog) {
          baselineLOdds[[i]] <- calcBaselineLOdds(or, 
                                                  f, nogV[i], p)
        }
      }
      else {
        nogForFunc <- rep(1, nog)
        for (i in 1:lNog) {
          baselineLOdds[[i]] <- calcBaselineLOddsFromSample(or[1:nogV[i]], 
                                                            f[1:nogV[i]], nogForFunc[1:nogV[i]], p, nid[, 
                                                                                                        1:nogV[i]])
        }
      }
    }
    offset <- (loopCount - 1) * subSize
    subGroup <- offset + 1:subSize
    if (varyEffects) {
      if (n > 10000) {
        subsubSize <- ceiling(subSize/10)
        for (i in 1:10) {
          subsub <- (i - 1) * subsubSize + (1:subsubSize)
          logORi = log(or)
          logOR <- t(apply(3 - nid[subsub, ], 1, function(x) (x - 
                                                                2 * f) * logORi))
          bLO = unlist(lapply(baselineLOdds, function(x) x$bLO))
          logodds <- t(bLO + apply(logOR, 1, function(x) cumsum(x)[nogV]))
          risk[offset + subsub, ] <- exp(logodds)/(1 + 
                                                     exp(logodds))
        }
      }
      else {
        logORi = log(or)
        logOR <- t(apply(3 - nid, 1, function(x) (x - 
                                                    2 * f) * logORi))
        bLO = unlist(lapply(baselineLOdds, function(x) x$bLO))
        logodds <- t(bLO + apply(logOR, 1, function(x) cumsum(x)[nogV]))
        risk[subGroup, ] <- exp(logodds)/(1 + exp(logodds))
      }
    }
    else {
      for (j in 1:lNog) {
        nogJ = nogV[j]
        numEE = rowSums(whichEE[, 1:nogJ])
        numEe = rowSums(whichEe[, 1:nogJ])
        numee = rowSums(whichee[, 1:nogJ])
        numalleles = 2 * numEE + numEe
        logodds = baselineLOdds[[j]]$bLO + (numalleles - 
                                              2 * f * nogJ) * log(or)
        risk[subGroup, j] = exp(logodds)/(1 + exp(logodds))
      }
    }
    dis[subGroup, ] = 1 * (risk[subGroup, ] > matrix(runif(subSize * 
                                                             lNog, 0, 1), nrow = subSize))
    if (rev(subGroup)[1] == n) 
      done = TRUE
  }
  parms <- baselineLOdds
  parms$nog <- nogV
  parms$seed <- seed
  list(parms = parms, risk = risk, disease = dis,g=nid, or=or, bLO=bLO)
}