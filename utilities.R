# Run the main loop for the complete data with or without zeros

lagfit = function(data, ydata, zeros=TRUE, plotlag=FALSE, plotfreq=FALSE){ 
  
  
  # plotlag - TRUE if you want plot.lagphase for each species
  # plotfreq - TRUE if you want frequency plot for each species
  
  outspecies<-NULL
  outscene <- NULL
  outlag <- NULL
  outlaglength <- NULL
  outfirstyear <- NULL
  outendyear <- NULL
  
  
    Species = unique(as.character(data$Species))
    
    # Run the main loop
    
    for( i in 1: length(Species)){
      
      sdata = get.species(data,ydata, species = Species[i], zeros = zeros)
      #print(sdata$Species)
      #sdata$Island=island
      if(length(sdata$Year) >= 2) { #Insufficient data check to fit models
      
        #outisland = c(outisland,island)
        outspecies = c(outspecies, Species[i])
        fit0 = lagphase(sdata, zeros=zeros)
        # Check for the scenario
        outlag = c(outlag, fit0$lagphase)
        outlaglength = c(outlaglength,fit0$lengthlag)
        outfirstyear = c(outfirstyear,fit0$Year[1])
        if(!fit0$lagphase){
          outendyear = c(outendyear, "NA")
          if(fit0$scene=="constant"){
            scene = 1
          }
          if(fit0$scene=="linear"){
            if(fit0$coef[2] >0) scene = 2
            if(fit0$coef[2] <0) scene = 3
          }
        }
        else{
          if(fit0$coef[2] < 0) scene = 4
          if(fit0$coef[2] > 0) scene = 5
          outendyear = c(outendyear, fit0$knots[1])
        
        }
        outscene = c(outscene, scene)
        
        if(plotlag) plot.lagphase(fit0)
        if(plotfreq) freqplot(fit0)
      }
    }
  
  out = data.frame(Species = outspecies, Scene = outscene, Lag =outlag, Laglength = outlaglength, FirstYear=outfirstyear, EndYear=outendyear)
  out
}

# Extract Frequency data for given island and Species from data
# If zeros=TRUE, include zeros in returned data

get.species <- function(x, y, species, zeros=TRUE)
{
  out <- subset(x, x$Species==species)
  out <- out[,c("Year","Frequency", "Specimens")]
  
  # Sort the data by Year
  yorder = order(out$Year)
  out = out[yorder,]
  
  indx = which(diff(out$Year)==0)
  if(length(indx)>0){
    out$Frequency[indx]=out$Frequency[indx]+out$Frequency[indx+1]
    out = out[-(indx+1),]
  }
  
  
  
  
  if(zeros)
  {
    yrs <- min(out$Year):max(out$Year)
    zeros <- as.data.frame(matrix(0,nrow=length(yrs),ncol=3))
    colnames(zeros) <- colnames(out)
    zeros[,"Year"] <- yrs
    j <- is.element(zeros[,"Year"],out[,"Year"])
    zeros[j,"Frequency"] <- out[,"Frequency"]
    j <- is.element(y[,"Year"],zeros[,"Year"])
    zeros[,"Specimens"] <- y[j,"Specimens"]
    out <- zeros
  }
  # Either way, Frequency is missing if Specimens=0
  #out$Frequency[out$Specimens==0] <- NA
  out <- as.list(out)
  out$Species <- species
  #out$Island <- island
  return(out)
}

# Main function. Give it a set of data where the columns include
#  Year
#  Frequency
#  Specimens
# It will find appropriate knots if not specified
# It will choose an appropriate order if not specified
# Just don't give it knots but no order
# If gam=TRUE, it will return a gam model instead.

lagphase <- function(data, knots=NULL, order=1, gam=FALSE,zeros=TRUE)
{
  # Set zeros to missing
  if(!zeros)
    data$Frequency[data$Frequency==0] <- NA
  else
    data$Frequency[data$Specimens==0] <- NA
  
  # Fit gam
  if(gam)
  {
    gamfit <- gam(Frequency ~ s(Year), offset=log(Specimens), data=data, family=poisson)
    gamfit$Year <- data$Year
    gamfit$name <- data$Species
    return(gamfit)
  }
  
  # Otherwise fit a glm
  # Check if knots==0 meaning no knots to be included
  if(!is.null(knots))
  {
    if(length(knots)==1)
    {
      if(knots==0) # i.e., no knots to be included
      {
        # Fit model with no knots
        suppressWarnings(fit <- glm(Frequency ~ 1, offset=log(Specimens), data=data, family=poisson, na.action=na.omit))
        fit$Year <- data$Year
        
        fit$name <- data$Species
        
        fit$data <- data
        fit$lengthlag <- NA
        fit$lagphase <- FALSE        
        class(fit) <- c("lagphase","glm","lm")
        return(fit)
      }
    }
  }
  
  # Otherwise proceed 
  # Choose order if not provided
  if(is.null(order))
  {
    if(!is.null(knots))
      stop("Not implemented. If you specify the knots, you need to specify the order.")
    fit1 <- lagphase(data, order=1)
    fit3 <- lagphase(data, order=3)
    bestfit <- fit1
    if(AICc(fit3) < AICc(bestfit))
      bestfit <- fit3
    return(bestfit)    
  }
  # Otherwise proceed with specified order
  if(!is.null(knots))
  {
    return(lagphase.knots(knots, data, order))
  }
  # Otherwise order specified but knots unspecified
  
  # Choose some initial knots
  knots <- quantile(data$Year,prob=c(0.2,0.4,0.6,0.8))
  names(knots) <- NULL
  
  # Fit best 4, 3, 2, 1 and 0 knot models
  
  fit4 <- optim(knots, tryknots, data=data, order=order)
  fit3 <- optim(knots[2:4], tryknots, data=data, order=order)
  fit2 <- optim(knots[c(2,4)], tryknots, data=data, order=order)
  fit1 <- optim(knots[2], tryknots, data=data, order=order, method="Brent", 
                lower=min(data$Year), upper=max(data$Year))
  suppressWarnings(fit0 <- glm(Frequency ~ 1, offset=log(Specimens), family=poisson, data=data, na.action=na.omit))
  fitl <- glm(Frequency ~ Year, offset=log(Specimens), family=poisson, data=data, na.action=na.omit)
  
  # Find best of these models:
  bestfit <- fit4
  if(fit3$value < bestfit$value)
    bestfit <- fit3
  if(fit2$value < bestfit$value)
    bestfit <- fit2
  if(fit1$value < bestfit$value)
    bestfit <- fit1
  if(AICc(fit0) < bestfit$value)
  {
    bestfit <- fit0
    bestfit$Year <- data$Year
    bestfit$name <- data$Species
    
    bestfit$data <- data
    bestfit$scene <- "constant"
  }
  else if(AICc(fitl) < bestfit$value)
  {
    bestfit <- fitl
    bestfit$Year <- data$Year
    bestfit$name <- data$Species
    bestfit$data <- data
    bestfit$scene <- "linear"
    
  }
  else  # Refit best model
    bestfit <- lagphase.knots(bestfit$par, data=data, order=order)
  
  if(is.null(bestfit$knots))
  {
    bestfit$lagphase <- FALSE
    bestfit$lengthlag <- NA
  }
  return(bestfit)
}


# Fit model with lag phase followed by growth
# where knots and order are specified
lagphase.knots <- function(knots, data, order)
{
  x <- matrix(NA,ncol=length(knots),nrow=length(data$Year))
  #x[,1] <- as.numeric(data$Year < knots[1])
  for(i in 1:length(knots))
    x[,i] <- pmax((data$Year-knots[i])^order,0)
  
  suppressWarnings(fit <- glm(Frequency ~ x, offset=log(Specimens), family=poisson, data=data, na.action=na.omit))
  fit$knots <- knots
  names(fit$knots) <- paste("K",1:length(knots),sep="")
  fit$Year <- data$Year
  fit$order <- order
  fit$name <- data$Species
  fit$data <- data
  
  # Check if there is a lag phase and record it
  fit$lengthlag <- NA
  if(length(fit$knots) > 0)
  {
    if(fit$coef[2] != 0)
      fit$lengthlag <- fit$knots[1] - min(data$Year)
  }
  fit$lagphase <- !is.na(fit$lengthlag)
  
  class(fit) <- c("lagphase","glm","lm")
  return(fit)
}

# Check that the specified knots make sense.
# Then use lagphase.knots to fit the model
# Returns AIC of fitted model
tryknots <- function(knots, data, order)
{
  # Knots must be interior to the data
  if(min(knots) < min(data$Year))
    return(1e50)
  if(max(knots) > max(data$Year))
    return(1e50)
  # Knots must be five Years apart and ordered
  if(length(knots) > 1)
  {
    dk <- diff(knots)
    if(min(diff(knots)) < 5)
      return(1e50)
  }
  
  # OK. Now fit the model
  fit <- lagphase.knots(knots, data, order)
  # Return the AICc of the fitted model
  return(AICc(fit))
}

# Function to return corrected AIC from a fitted object
AICc <- function(object)
{
  aic <- object$aic
  k <- length(object$coefficients)
  n <- object$df.residual+k
  aicc <- aic + 2*k*(k+1)/(n-k-1)
  if(aicc == Inf) aicc=1e50
  return(aicc)  
}

# Produces plot of the fitted spline function after adjusting for 
# number of Specimens

plot.lagphase <- function(fit,ylim=NULL,xlab="Year", ylab="Adjusted Frequency", main=fit$name,...)
{
  fits <- predict(fit, se.fit=TRUE)
  
  #Specimens <- model.matrix(fit)[,"Specimens"]
  #adjfits <- fits$fit - coef(fit)["Specimens"]*(Specimens - mean(Specimens))
  ladjfits <- fits$fit - fit$offset + mean(fit$offset,na.rm=TRUE)
  adjfits <- exp(ladjfits)
  up <- exp(ladjfits + 2*fits$se.fit)
  lo <- exp(ladjfits - 2*fits$se.fit)
  if(is.null(ylim))
    ylim <- range(lo,pmin(up,3*adjfits))
    
  j <- (fit$data$Specimens > 0)
 
  plot(fit$data$Year[j],adjfits, ylim=ylim,  type="n", xlab=xlab,ylab=ylab,
       main=main,...)
  polygon(c((fit$data$Year[j]),rev(fit$data$Year[j])),c(lo,rev(up)),border=FALSE,col="gray")
  lines(fit$data$Year[j],adjfits)
  if(!is.null(fit$knots))
    abline(v=fit$knots[1],col="gray")
  rug(fit$Year[fit$data$Frequency > 0])
}


freqplot <- function(fit1, fit2=NULL, fit3=NULL, fit4=NULL, 
                     xlab="Year", ylab="Frequency", main=fit1$name, cols=2:5, ...)
{
  if(is.element("data",names(fit1)))
    data <- fit1$data
  else
  {
    data <- fit1
    fit1 <- NULL
  }
  
  plot(Frequency ~ Year, data=data, xlab=xlab, ylab=ylab, main=main, 
       ylim=range(0,data$Frequency,na.rm=TRUE),...)
  j <- (data$Specimens > 0)
  if(!is.null(fit1))
    lines(data$Year[j],fitted(fit1),col=cols[1])
  if(!is.null(fit2))
    lines(data$Year[j],fitted(fit2),col=cols[2])
  if(!is.null(fit3))
    lines(data$Year[j],fitted(fit3),col=cols[3])
  if(!is.null(fit4))
    lines(data$Year[j],fitted(fit4),col=cols[4])
}

# Fit piecewise linear model
pwlm <- function(x,y)
{
  # choose knot
  minmse <- Inf
  xgrid <- seq(min(x),max(x),l=102)[-c(1,102)]
  for(k in xgrid)
  {
    x2 <- pmax(x-k,0)
    fit <- lm(y ~ x + x2)
    res <- residuals(fit)
    mse <- mean(res^2)
    if(mse < minmse)
    {
      bestfit <- fit
      bestfit$k <- k
      minmse <- mse
    }
  }
  return(bestfit)
}

print.lagphase <- function(x,...)
{
  if(x$lagphase)
  {
    cat(paste("Lag phase: ",x$data$Year[1],"-",round(x$data$Year[1]+x$lengthlag),
              " (",round(x$lengthlag)," Years)\n",sep=""))
  }
  else
  {  
    cat("No lag phase identified\n")
  }
  if(length(x$coef) > 1)
  {
    cat("  Knots:",x$knots,"\n")
  }
  stats:::print.lm(x,...) 
}
