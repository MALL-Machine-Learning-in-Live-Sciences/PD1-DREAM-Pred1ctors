# Utils functions
standarize = function(df, log2 = T, scale = T){
  
  # Details:
  #   Standarize expression matrix (from counts to scale)
  # Arguments:
  #   df    : expression data.frame with patients in rows and genes in columns
  #   log2  : logical value. If TRUE, log2 normalization will be calculated. If FALSE, only scale will be carried out
  #   scale : logical value. If TRUE, scale function will be applied
  # Value:
  #   Data.frame with standarize genes
  
  if (log2 == T){
    print('Applying log2 transformation...')
    df = apply(df, 2, function(x) log2(x + 1))
  }
  
  if (scale == T){
    print('Scaling dataset...')
    df = scale(df)
  }
  
  res = df
  return(res)
}

dataTrans <- function(data, x, y, z, tt, std.x, std.i, std.tt, inter, trace = TRUE){
  
  #########################################################
  ### Control and checking
  
  if(length(y) != 2)
    stop("\nThe outcome must contains 2 columns: 'time' and 'status'.")
  if(min(data[, y[1]]) < 0)
    stop("\nTimes must be non-negative.")
  if(sum(data[, y[2]] %in% c(0, 1)) != nrow(data))
    stop("\nStatus must be either 0 (censor) or 1 (event).")
  
  if(!is.null(tt)){
    TT <- as.data.frame(data[, tt])
    if(length(tt) > 1)
      stop("\nThe treatment must be a single variable.")
    if(length(tt) == 1){
      if(length(unique(as.vector(t(data[, tt])))) > 2)
        stop(paste0("\nThe treatment variable must consider only two groups."))
      # if(length(unique(as.vector(t(data[, tt])))) == 1){
      #   warning(paste0("\nAll patients are in the same treatment group. The analysis is then switch to a prognostic setting."))
      #   tt <- NULL
      #   inter <- FALSE
      # }
      if((sum(unique(as.vector(t(data[, tt]))) %in% c(-0.5, +0.5)) != 2) & std.tt == TRUE)
        data[, tt] <- as.numeric(factor(as.vector(t(data[, tt])))) - 1.5
    }
  }
  
  #########################################################
  
  itt <- tt; iz <- z; ix <- x; iy <- y
  iptt <- which(colnames(data) %in% tt)
  ipz <- which(colnames(data) %in% z)
  ipx <- which(colnames(data) %in% x)
  ipy <- which(colnames(data) %in% y)
  isSim <- (!is.null(attributes(data)$isSim))
  
  if(std.x == TRUE)
    data[, x] <- scale(data[, x], center = T, scale = T)
  
  if(inter == TRUE){
    XT <- as.matrix(data[, x]) * matrix(data[, tt], nrow = nrow(data), ncol = length(x))
    if(std.i == TRUE)
      XT <- scale(XT, center = T, scale = T)
    colnames(XT) <- paste0("bi", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
    xt <- colnames(XT)
  }else{
    xt <- NULL
  }
  
  data <- cbind(data[, c(tt, z, x, y)])
  if(inter == TRUE){
    data <- cbind(data, XT)
    colnames(data)[1] <- "treat"
  }
  
  tnames <- c(rep("tt", length(tt)), rep("z", length(z)), rep("x", length(x)), rep("y", length(y)), rep("xt", length(xt)))
  if(!is.null(z))
    colnames(data)[tnames == "z"] <- paste0("cl", gsub(" ", "0", format(c(length(z), 1:length(z))))[-1])
  colnames(data)[tnames == "x"] <- paste0("bm", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
  colnames(data)[tnames == "y"] <- c("time", "status")
  
  data <- na.omit(data)
  if(!is.null(attributes(data)$na.action) & trace == TRUE){
    nmiss <- length(attributes(data)$na.action)
    message(paste0(
      "\rData management: ", nmiss, " observation", ifelse(nmiss > 1, "s were", " was"), " excluded due to missing data."))
  }
  
  if(!(class(unlist(data)) %in% c("numeric", "integer")))
    stop("\nAll variables must be numerical.")
  
  attributes(data) <- append(
    x = attributes(data),
    values = list(
      inter = inter,
      inames = list(
        tt = itt,
        z = iz,
        x = ix,
        y = iy),
      ipos = list(
        tt = iptt,
        z = ipz,
        x = ipx,
        y = ipy
      ),
      tnames = tnames,
      pos = list(
        z = grep("cl", colnames(data)),
        x = grep("bm", colnames(data)),
        xt = grep("bi", colnames(data)),
        X = (1:ncol(data))[-which(colnames(data) %in% c("time", "status"))],
        y = which(colnames(data) %in% c("time", "status"))
      ),
      weights = c(rep(0, length(c(tt, z))), rep(1, length(x) * ((inter == TRUE) + 1))),
      std.x = std.x,
      std.tt = std.tt)
  )
  if(inter == TRUE){
    attributes(data)$pos$tt <- grep("treat", colnames(data))
    attributes(data)$std.i <- std.i
    attributes(data)$inames$xt <- paste0(ix, ":", itt)
  }
  
  return(data)
}

predBiospear= function(res, method, newdata, randomPFS){
  
  # Arguments
  #   model: resBMsel object
  #   method: selected method to make the predictions
  #   newdata: validation data
  
  if (randomPFS == T){
    newdata$PFS = sample(c(1:1825), size = nrow(newdata), replace = T)
    newdata$PFS.Event = sample(c(0,1), size = nrow(newdata), replace = T)
  }
  
  tt = attributes(res)$inames[which(attributes(res)$tnames ==  "tt")]
  x = attributes(res)$inames[which(attributes(res)$tnames ==  "x")]
  z = attributes(res)$inames[which(attributes(res)$tnames ==  "z")]
  y = attributes(res)$inames[which(attributes(res)$tnames ==  "y")]
  isRidge <- (unique(!is.na(attributes(res)$ridge)))
  
  Res <- data.frame(summary(res, show = FALSE, add.ridge = isRidge))
  Res <- Res[, gsub("-", ".", method), drop = FALSE]
  
  if (attributes(res)$inter == TRUE) {
    Res.i <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
    Res.i <- Res.i[, gsub("-", ".", method), drop = FALSE]
  }
  nmeth <- ncol(Res)
  
  newdata <- dataTrans(data = newdata, x = x, y = y, z = z, 
                       tt = tt, std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, 
                       std.tt = attributes(res)$std.tt, inter = attributes(res)$inter, 
                       trace = TRUE)
  
  colnames(newdata) <- attributes(res)$inames
  surv.new <- Surv(newdata[, attributes(res)$inames[which(attributes(res)$tnames == "y")][1]],
                   newdata[, attributes(res)$inames[which(attributes(res)$tnames == "y")][2]])
  
  lp.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
  if (nrow(Res) > 0) 
    lp.new <- as.matrix(newdata[, rownames(Res)]) %*% 
    as.matrix(Res)
  if (attributes(res)$inter == TRUE) {
    lpint.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
    if (nrow(Res.i) > 0 & sum(Res.i) != 0) 
      lpint.new <- as.matrix(newdata[, gsub(paste0(":", attributes(res)$inames[1]), "", rownames(Res.i))]) %*% as.matrix(Res.i)
  }
  
  return(lp.new)
  
}