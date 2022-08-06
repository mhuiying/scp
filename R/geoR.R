# These are function taken from geoR version 1.8-1
#
# The code is also here https://rdrr.io/cran/geoR/src/R/variogram.R
#

#' Matern function
#'
#' @param u
#' @param phi
#' @param kappa
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
matern <- function (u, phi, kappa){
    if (is.vector(u))
        names(u) <- NULL
    if (is.matrix(u))
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf,
        gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)),
        1)
    uphi[u > 600 * phi] <- 0
return(uphi)}


#' variog
#'
#' @param geodata
#' @param coords
#' @param data
#' @param uvec
#' @param breaks
#' @param trend
#' @param lambda
#' @param option
#' @param estimator.type
#' @param nugget.tolerance
#' @param max.dist
#' @param pairs.min
#' @param bin.cloud
#' @param direction
#' @param tolerance
#' @param unit.angle
#' @param angles
#' @param messages
#' @param ...
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
variog <- function (geodata, coords = geodata$coords, data = geodata$data,
    uvec = "default", breaks = "default", trend = "cte", lambda = 1,
    option = c("bin", "cloud", "smooth"), estimator.type = c("classical",
        "modulus"), nugget.tolerance, max.dist, pairs.min = 2,
    bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,
    unit.angle = c("radians", "degrees"), angles = FALSE, messages,
    ...){

    if (missing(geodata))
        geodata <- list(coords = coords, data = data)
    call.fc <- match.call()
    if (missing(messages))
        messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")),
            TRUE, getOption("geoR.messages")))
    else messages.screen <- messages
    keep <- list(...)
    if (is.null(keep$keep.NA))
        keep.NA <- FALSE
    else keep.NA <- keep$keep.NA
    unit.angle <- match.arg(unit.angle)
    if (mode(direction) == "numeric") {
        if (length(direction) > 1)
            stop("only one direction is allowed")
        if (length(tolerance) > 1)
            stop("only one tolerance value is allowed")
        if (unit.angle == "degrees") {
            ang.deg <- direction
            ang.rad <- (ang.deg * pi)/180
            tol.deg <- tolerance
            tol.rad <- (tol.deg * pi)/180
        }
        else {
            ang.rad <- direction
            ang.deg <- (ang.rad * 180)/pi
            tol.rad <- tolerance
            tol.deg <- (tol.rad * 180)/pi
        }
        if (ang.rad > pi | ang.rad < 0)
            stop("direction must be an angle in the interval [0,pi[ radians")
        if (tol.rad > pi/2 | tol.rad < 0)
            stop("tolerance must be an angle in the interval [0,pi/2] radians")
        if (tol.deg >= 90) {
            direction <- "omnidirectional"
            cat("variog: computing omnidirectional variogram\n")
        }
        else {
            if (messages.screen) {
                cat(paste("variog: computing variogram for direction = ",
                  round(ang.deg, digits = 3), " degrees (", round(ang.rad,
                    digits = 3), " radians)\n", sep = ""))
                cat(paste("        tolerance angle = ", round(tol.deg,
                  digits = 3), " degrees (", round(tol.rad, digits = 3),
                  " radians)\n", sep = ""))
            }
        }
    }
    else if (messages.screen)
        cat("variog: computing omnidirectional variogram\n")
    coords <- as.matrix(coords)
    data <- as.matrix(data)
    if (nrow(coords) != nrow(data))
        stop("coords and data have incompatible dimensions")
    data.var <- apply(data, 2, var)
    n.data <- nrow(coords)
    n.datasets <- ncol(data)
    data <- drop(data)
    option <- match.arg(option)
    estimator.type <- match.arg(estimator.type)
    if (abs(lambda - 1) > 1e-04) {
        if (abs(lambda) < 1e-04)
            data <- log(data)
        else data <- ((data^lambda) - 1)/lambda
    }
    xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
    if (nrow(xmat) != n.data)
        stop("coords and trend have incompatible sizes")
    if (trend != "cte") {
        if (is.vector(data)) {
            temp.fit <- lm(data ~ xmat + 0)
            beta.ols <- temp.fit$coeff
            data <- temp.fit$residuals
            temp.fit <- NULL
            names(data) <- NULL
        }
        else {
            only.res <- function(y, x) lm(y ~ xmat + 0)$residuals
            data <- apply(data, 2, only.res, x = xmat)
            only.beta <- function(y, x) lm(y ~ xmat + 0)$coef
            beta.ols <- apply(data, 2, only.beta, x = xmat)
        }
    }
    else beta.ols <- colMeans(as.matrix(data))
    u <- as.vector(dist(as.matrix(coords)))
    if (missing(nugget.tolerance) || nugget.tolerance < 1e-11) {
        nugget.tolerance <- 1e-12
        nt.ind <- FALSE
    }
    else {
        if (mode(nugget.tolerance) != "numeric")
            stop("nugget.tolerance must be numeric")
        nt.ind <- TRUE
    }
    min.dist <- min(u)
    if (min.dist < nugget.tolerance)
        nt.ind <- TRUE
    if (direction != "omnidirectional" | angles) {
        u.ang <- .C("tgangle", as.double(as.vector(coords[, 1])),
            as.double(as.vector(coords[, 2])), as.integer(dim(coords)[1]),
            res = as.double(rep(0, length(u))), PACKAGE = "geoR")$res
        if (any(is.na(u.ang)))
            stop("NA returned in angle calculations maybe due to co-located data")
        u.ang <- atan(u.ang)
        u.ang[u.ang < 0] <- u.ang[u.ang < 0] + pi
    }
    if (option == "bin" && bin.cloud == FALSE && direction ==
        "omnidirectional") {
        if (missing(max.dist))
            umax <- max(u)
        else umax <- max(u[u < max.dist])
        dbins <- .define.bins(max.dist = umax, uvec = uvec, breaks = breaks,
            nugget.tolerance = nugget.tolerance)
        uvec <- dbins$uvec
        bins.lim <- dbins$bins.lim
        nbins <- length(bins.lim) - 1
        if (missing(max.dist))
            max.dist <- max(bins.lim)
        if (bins.lim[1] < 1e-16)
            bins.lim[1] <- -1
        bin.f <- function(data) {
            cbin <- vbin <- sdbin <- rep(0, nbins)
            .C("binit", as.integer(n.data), as.double(as.vector(coords[,
                1])), as.double(as.vector(coords[, 2])), as.double(as.vector(data)),
                as.integer(nbins), as.double(as.vector(bins.lim)),
                as.integer(estimator.type == "modulus"), as.double(max.dist),
                cbin = as.integer(cbin), vbin = as.double(vbin),
                as.integer(TRUE), sdbin = as.double(sdbin), PACKAGE = "geoR")[c("vbin",
                "cbin", "sdbin")]
        }
        result <- array(unlist(lapply(as.data.frame(data), bin.f)),
            dim = c(nbins, 3, n.datasets))
        indp <- (result[, 2, 1] >= pairs.min)
        result[!indp, 1, ] <- NA
        if (bins.lim[1] < 0)
            bins.lim[1] <- 0
        if (!nt.ind) {
            uvec <- uvec[-1]
            indp <- indp[-1]
            bins.lim <- bins.lim[-1]
            result <- result[-1, , , drop = FALSE]
        }
        if (keep.NA)
            result <- list(u = uvec, v = result[, 1, ], n = result[,
                2, 1], sd = result[, 3, ], bins.lim = bins.lim,
                ind.bin = indp)
        else result <- list(u = uvec[indp], v = result[indp,
            1, ], n = result[indp, 2, 1], sd = result[indp, 3,
            ], bins.lim = bins.lim, ind.bin = indp)
    }
    else {
        data <- as.matrix(data)
        v <- matrix(0, nrow = length(u), ncol = n.datasets)
        for (i in 1:n.datasets) {
            v[, i] <- as.vector(dist(data[, i]))
            if (estimator.type == "modulus")
                v[, i] <- v[, i, drop = FALSE]^(0.5)
            else v[, i] <- (v[, i, drop = FALSE]^2)/2
        }
        if (!missing(max.dist)) {
            v <- v[u <= max.dist, , drop = FALSE]
            if (direction != "omnidirectional")
                u.ang <- u.ang[u <= max.dist]
            u <- u[u <= max.dist]
        }
        if (direction != "omnidirectional") {
            ang.lower <- ang.rad - tol.rad
            ang.upper <- ang.rad + tol.rad
            if (ang.lower >= 0 & ang.upper < pi)
                ang.ind <- (!is.na(u.ang) & ((u.ang >= ang.lower) &
                  (u.ang <= ang.upper)))
            if (ang.lower < 0)
                ang.ind <- (!is.na(u.ang) & ((u.ang < ang.upper) |
                  (u.ang > (pi + ang.lower))))
            if (ang.upper >= pi)
                ang.ind <- (!is.na(u.ang) & ((u.ang > ang.lower) |
                  (u.ang < (ang.upper - pi))))
            v <- v[ang.ind, , drop = FALSE]
            u <- u[ang.ind]
        }
        data <- drop(data)
        v <- drop(v)
        if (option == "cloud") {
            result <- list(u = u, v = v)
            if (angles)
                result$angles <- u.ang
        }
        if (option == "bin") {
            if (missing(max.dist))
                umax <- max(u)
            else umax <- max(u[u < max.dist])
            if (bin.cloud == "diff")
                dd <- diffpairs(coords, data)$diff
            else dd <- 0
            result <- .rfm.bin(cloud = list(u = u, v = v, d = dd),
                estimator.type = estimator.type, uvec = uvec,
                breaks = breaks, nugget.tolerance = nugget.tolerance,
                bin.cloud = bin.cloud, max.dist = umax, keep.NA = keep.NA)
            if (keep.NA) {
                if (pairs.min > 0) {
                  indp <- (result$n < pairs.min)
                  if (!nt.ind) {
                    for (i in 1:5) result[[i]] <- result[[i]][-1]
                    indp <- indp[-1]
                  }
                  if (is.matrix(result$v)) {
                    result$v[indp, ] <- result$sd[indp, ] <- NA
                  }
                  else {
                    result$v[indp] <- result$sd[indp] <- NA
                  }
                }
                result$ind.bin <- indp
            }
            else {
                if (pairs.min > 0) {
                  if (!nt.ind) {
                    for (i in 1:5) result[[i]] <- result[[i]][-1]
                  }
                  indp <- (result$n >= pairs.min)
                  if (is.matrix(result$v)) {
                    result$v <- result$v[indp, ]
                    result$sd <- result$sd[indp, ]
                  }
                  else {
                    result$v <- result$v[indp]
                    result$sd <- result$sd[indp]
                  }
                  result$u <- result$u[indp]
                  result$n <- result$n[indp]
                }
                result$ind.bin <- indp
            }
        }
        if (option == "smooth") {
            if (is.matrix(v))
                stop("smooth not yet available for more than one data-set")
            temp <- ksmooth(u, v, ...)
            result <- list(u = temp[[1]], v = temp[[2]])
        }
        if (missing(max.dist))
            max.dist <- max(u)
    }
    if (nt.ind) {
        if (!exists(".variog4.nomessage"))
            cat("variog: co-locatted data found, adding one bin at the origin\n")
        if (all(result$u[1:2] < 1e-11))
            result$u[2] <- sum(result$bins.lim[2:3])/2
    }
    result <- c(result, list(var.mark = data.var, beta.ols = beta.ols,
        output.type = option, max.dist = max.dist, estimator.type = estimator.type,
        n.data = n.data, lambda = lambda, trend = trend, pairs.min = pairs.min))
    result$nugget.tolerance <- nugget.tolerance
    if (direction != "omnidirectional")
        result$direction <- ang.rad
    else result$direction <- "omnidirectional"
    if (direction != "omnidirectional")
        result$tolerance <- tol.rad
    else result$tolerance <- "none"
    result$uvec <- uvec
    result$call <- call.fc
    oldClass(result) <- "variogram"
return(result)}





#' Trend Spatial
#'
#' @param trend
#' @param geodata
#' @param add.to.trend
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
trend.spatial <- function (trend, geodata, add.to.trend){
    if (!missing(geodata)) {
        gd <- geodata
        if (any(class(geodata) %in% c("SpatialGridDataFrame",
            "SpatialLinesDataFrame", "SpatialPolygonsDataFrame",
            "SpatialPointsDataFrame", "SpatialPixelsDataFrame"))) {
            geodata <- as.data.frame(gd)
            geodata$coords <- coordinates(gd)
        }
        if (any(class(geodata) == "geodata")) {
            geodata <- as.data.frame(gd)
            geodata$coords <- gd$coords
        }
    }
    else geodata <- data.frame()
    if (inherits(trend, "formula")) {
        trend.mat <- try(model.matrix(trend, data = geodata),
            silent = TRUE)
        if (inherits(trend.mat, "try-error"))
            stop("\ntrend elements/variables not found")
    }
    else {
        if (mode(trend) == "numeric")
            trend.mat <- unclass(trend)
        else if (trend == "cte") {
            if (missing(geodata))
                stop("argument geodata must be provided with trend=\"cte\"")
            trend.mat <- as.matrix(rep(1, nrow(geodata$coords)))
        }
        else if (trend == "1st") {
            if (missing(geodata))
                stop("argument geodata must be provided with trend=\"1st\"")
            trend.mat <- cbind(1, geodata$coords)
        }
        else if (trend == "2nd") {
            if (missing(geodata))
                stop("argument geodata must be provided with trend=\"2nd\"")
            trend.mat <- cbind(1, geodata$coords, geodata$coords[,
                1]^2, geodata$coords[, 2]^2, geodata$coords[,
                1] * geodata$coords[, 2])
        }
        else stop("external trend must be provided for data locations to be estimated using the arguments trend.d and trend.l. Allowed values are the strings \"cte\", \"1st\", \"2nd\" or  a model formula")
    }
    trend.mat <- as.matrix(trend.mat)
    if (!missing(add.to.trend)) {
        if (missing(geodata))
            trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend)[,
                -1])
        else trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend,
            geodata = geodata)[, -1])
    }
    dimnames(trend.mat) <- list(NULL, NULL)
    oldClass(trend.mat) <- "trend.spatial"
return(trend.mat)}


#' Define Bins
#'
#' @param max.dist
#' @param uvec
#' @param breaks
#' @param nugget.tolerance
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
.define.bins <- function(max.dist, uvec = "default", breaks = "default", nugget.tolerance){
  if(all(breaks ==  "default")){
    if (all(uvec == "default")) uvec <- 13
    if (mode(uvec) == "numeric"){
      if(length(uvec) == 1){
        bins.lim <- seq(0, max.dist, l = uvec+1)
        bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
        uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
      }
      else{
        uvec <- c(0, uvec)
        nvec <- length(uvec)
        d <- 0.5 * diff(uvec[2:nvec])
        bins.lim <- c(0, (uvec[2:(nvec - 1)] + d), (d[nvec - 2] + uvec[nvec]))
        bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
      }
    }
    else stop("argument uvec can only take a numeric vector")
  }
  else{
    if(mode(breaks) != "numeric") stop("argument breaks can only take a numeric vector")
    else bins.lim <- breaks
    bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
    uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
  }
  return(list(uvec = uvec, bins.lim = bins.lim))
}
