

## ------------------------------------------------------------------------
#OBS! this installs the pacman package, and any of the others listed!
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, zoo, Kendall, grid, boot, RColorBrewer, NADA)

## ------------------------------------------------------------------------
#ggplot theme used
taCont_theme <- function() {

  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[1]
  color.grid.major = palette[2]
  color.facet = palette[4]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]

  # Begin construction of chart
  theme_bw(base_size=9) +

    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +

    # Format the grid
    theme(panel.grid.major=element_line(colour = color.grid.major)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +

    # Format the legend, but hide by default
    theme(legend.position="right") +
    theme(legend.title = element_text(size =7 , colour = color.axis.title))+
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=8,color=color.axis.title)) +

    # Format facets
    theme(strip.background = element_rect(fill = color.facet, color = color.facet)) +
    theme(strip.text = element_text(face = "bold", colour = color.background)) +

    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
    theme(axis.text.x=element_text(size=8,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=8,color=color.axis.text)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_text(size=9,color=color.axis.title, vjust=1.25)) +
    theme(axis.line = element_blank())+

    # Plot margins
    theme(plot.margin = unit(c(0.2, 0.3, 0.3, 0.2), "cm"))
}

#Additional funtions needed for the temporal trend analysis

#Function to calculate delta
## ------------------------------------------------------------------------
delta <- function(phi, slope, yrs) {
  (log(1+slope/100))^2*((yrs-1)*yrs*(yrs+1)/(12*phi^2))
}

#f-distribution at alpha 0.05 for yrs amount of years
## ------------------------------------------------------------------------
fdist <-function(yrs) {
  qf(0.95,1,yrs-2,ncp=0)
}

#Function with known phi and slope to calculate how many years are needed
#for a power of 80%
## ------------------------------------------------------------------------
poweryrs<-function(phi, slope) {
  i<-seq(2.2,50, by=0.1)
  pow<-numeric(length(i))
  # no. of sampling years with particular slope (in %)
  delta<-log(1+slope/100)^2*(((i-1)*i*(i+1))/(12*phi^2))
  f<-qf(0.95,1,i-2,ncp=0)
  pow<-1-pf(f,1,i-2,ncp=delta)
  data<-data.frame(i,pow)
  yrs80 <- subset(data, pow>=0.8)
  yrs80$i[1]
}

#Function with known phi and years to calculate the slope at a power of 80%
## ------------------------------------------------------------------------
powerslope<-function(phi, yrs) {
  slope<-seq(0.1,1000, by=0.1)
  pow<-numeric(length(slope))

  delta<-log(1+slope/100)^2*(((yrs-1)*yrs*(yrs+1))/(12*phi^2))
  f<-qf(0.95,1,yrs-2,ncp=0)
  pow<-1-pf(f,1,yrs-2,ncp=delta)
  data<-data.frame(slope,pow)
  yrs80 <- subset(data, pow>=0.8)
  yrs80$slope[1]
}

#Function with known phi and slope to get powercurve

#this results in a data frame which has the powercurve for a particular phi and slope that can easily be plottet
#i = number of years (0-50)
#pow =  the calculated power (0-1)
## ------------------------------------------------------------------------
powerdf<-function(phi, slope) {
  i<-seq(2.2,50, by=0.1)
  pow<-numeric(length(i))
  # no. of sampling years with particular slope (in %)
  delta<-log(1+slope/100)^2*(((i-1)*i*(i+1))/(12*phi^2))
  f<-qf(0.95,1,i-2,ncp=0)
  pow<-1-pf(f,1,i-2,ncp=delta)
  data<-data.frame(i,pow)
}

#Bootsrapping used to calculate the confidence interval (CI) of the median and the ROS model

#The number of replicas is either R1 = 1000 or R2 = 9999
#depending on the computer power you have, these can be increased here
##R1 used for the ROS models
##R2 used everywhere else

## ------------------------------------------------------------------------
median.boot <- function(x,i){
  median(x[i], na.rm = TRUE)
}

R1 <- 1000
R2 <- 9999

## ------------------------------------------------------------------------

#The temporal analysis of contaminants function taCont()
taCont <- function(date, y, censored = NULL,
                   plot=TRUE, pub=TRUE, palmost=TRUE,
                   cenPerc = NULL, onlyRes=FALSE) {
  cen <- if(is.null(censored)) {rep(FALSE, length(y))} else {censored}
  ylab <- deparse(substitute(y))
  Perclimit <- if(is.null(cenPerc)) {80} else {cenPerc}
  if(is.logical(cen) == FALSE) {warning("censored has to be a logical vector")
  }

  #if the y data is negative, transformations are used for the calculations, this is only valid if all of the y data are negative!!
  ynegative <- ifelse(sum(y > 0, na.rm = TRUE) == 0, TRUE, FALSE)
  if(ynegative == FALSE & sum(y < 0, na.rm = TRUE) > 0) {
    warning("your y data consists of both negative and positive numbers, you have to transform the data by moving it up or down the y-axis")
  }
  y <- if(ynegative == TRUE){abs(y)} else {y}

  #if the y data consists of small numbers <1, transformations are used to handle error resulting from log()
  ytransformed <- ifelse(min(y, na.rm = TRUE) < 1, TRUE, FALSE)
  y <- if(ytransformed == TRUE){y*1000} else{y}

  #date or year
  x <- if(lubridate::is.POSIXct(date)){as.numeric(lubridate::year(date))} else{date}

  #set up a dataframe
  datafr<-data.frame(date,y,cen,x)
  datafr <- datafr %>%
    filter(!is.na(y), date > 0) %>%
    mutate(yplot = ifelse(cen == TRUE, y/sqrt(2), y), #this is used for plotting <LOQ
           lny = log(y))
  datafr$yplot <- if(ytransformed == TRUE) {datafr$yplot/1000} else {datafr$yplot}
  datafr$yplotplot <- if(ynegative == TRUE){datafr$yplot * -1} else {datafr$yplot}

  yrs_start <- min(datafr$x)
  yrs_end <- max(datafr$x)

  #check for <LOQ above yearly mean, because this drops values in cenros() for bootstrapping
  datafr <- datafr %>%
    group_by(x) %>%
    arrange(lny) %>%
    mutate(mean = mean(lny, na.rm = TRUE),
           LOQLrgMean = ifelse(cen == TRUE & lny > mean, TRUE, FALSE))

  #summarize the data pr. year
  datamaAll <- datafr %>%
    group_by(x) %>%
    summarise(N = sum(!is.na(lny)),
              N_cen = sum(cen, na.rm = TRUE),
              PercCen = sum(cen, na.rm = TRUE)/sum(!is.na(lny))*100,
              N_LOQexcl = sum(LOQLrgMean, na.rm = TRUE),
              mean = mean(lny, na.rm = TRUE),
              median_unCen = median(lny, na.rm = TRUE),
              sum = sum(lny, na.rm = TRUE),
              min = min(lny, na.rm = TRUE),
              max = max(lny, na.rm = TRUE))

  #xcen is the years where censored data is above a certain percentage
  xcen <- (datamaAll %>% filter(PercCen > Perclimit))$x
  if(length(xcen)) {
    datafr <- datafr %>%
      mutate(YearExcluded = ifelse(x %in% xcen, TRUE, FALSE))
  } else {
    datafr$YearExcluded <- rep(FALSE, length(datafr$x))
  }

  #cenros model results of filtered data
  datafr <- datafr %>%
    group_by(x) %>%
    arrange(lny) %>%
    mutate(ROSmodel = ifelse(YearExcluded == TRUE, NA,
                             ifelse(LOQLrgMean == TRUE, NA,
                                    tryCatch(cenros(lny,cen)$modeled, error=function(e)NA))))

  dataCI <- datafr %>%
    group_by(x) %>%
    summarise(median = median(ROSmodel, na.rm = TRUE),
              med_5CI = tryCatch(quantile(boot::boot(ROSmodel,median.boot,R=R1)$t, 0.05),
                                 error=function(e)NA)[[1]],
              med_95CI = tryCatch(quantile(boot::boot(ROSmodel,median.boot,R=R1)$t, 0.95),
                                  error=function(e)NA)[[1]]
    )

  datamaAll <- left_join(datamaAll, dataCI, by = "x")

  yrs <- paste(yrs_start, "-", yrs_end, sep = "")
  n_yrs = c(paste(length(datamaAll$x)," (", sum(!is.na(datamaAll$median)), ")", sep = ""), NA)
  n <- c("all", "last 10")
  censored <- if(is.null(cen)){c(NA,NA)
  } else {
    c(sum(datafr$cen, na.rm = TRUE)/sum(!is.na(datafr$lny))*100, NA)
  }

  results <- as.data.frame((matrix(NA, nrow = 2, ncol = 20)))
  colnames(results) <- c("n", "censored", "yrs", "n_yrs", "n_tot",
                         "median", "median_CI_5", "median_CI_95",
                         "annual_change", "ac_CI_low", "ac_CI_upp",
                         "CV_linear", "r2", "p_linear", "lowest_linear_det_change", "pow", "adequacy",
                         "lowest_nonlinear_det_change", "CV_nonlin", "p_nonlin")
  results$n_yrs <- n_yrs
  results$n <- n
  results$censored <- censored
  datamaPerc <- datamaAll %>% filter(PercCen >0)
  datamaFirstPlot <- datamaAll %>%
    group_by(x) %>%
    filter(!is.na(median)) %>%
    mutate(median = ifelse(ytransformed == TRUE, exp(median)/1000, exp(median)),
           med_5CI = ifelse(ytransformed == TRUE, exp(med_5CI)/1000, exp(med_5CI)),
           med_95CI = ifelse(ytransformed == TRUE, exp(med_95CI)/1000, exp(med_95CI)))

  # rules for the plot
  lines <- 0.6
  coldot <- "black"
  colbar <- "grey10"
  colline <- "violetred1"
  collinesm <- "springgreen"
  colline10 <- "royalblue"
  censdot <- "violetred4"
  miny <- if(ynegative == TRUE){min(datafr$yplotplot)-1} else {0}
  maxy <- if(ynegative == TRUE){max(datafr$yplotplot)-1} else {max(datafr$yplot)}
  maxypub <- if(ynegative == TRUE){abs(maxy)*2} else {maxy*2.5}
  minx <-if((min(datafr$x)%%2)==0){min(datafr$x)}else{min(datafr$x)-1}
  by <- if((max(datafr$x)-minx)<10){1}else{4}
  lowhorizontal <- if(ynegative == TRUE){min(datafr$yplotplot)-1} else {0}


  gcen<-ggplot(datamaFirstPlot, aes(x=x, y=median)) +
    geom_jitter(data=datafr %>% filter(YearExcluded == FALSE, LOQLrgMean == FALSE),
                aes(x=x, y=yplotplot),
                position = position_jitter(0.1), size = 1,
                colour= coldot,
                alpha = 0.25)+
    geom_jitter(data=datafr %>% filter(YearExcluded == TRUE),
                aes(x=x, y=yplotplot),
                position = position_jitter(0.1), size = 1,
                colour= censdot,
                alpha = 0.25)+
    geom_jitter(data=datafr %>% filter(LOQLrgMean == TRUE),
                aes(x=x, y=yplotplot, colour = cen),
                position = position_jitter(0.1), size = 1,
                colour= censdot,
                shape = 8,
                alpha = 0.5)+
    geom_errorbar(data=datamaFirstPlot, aes(ymin=med_5CI, ymax=med_95CI),
                  width = 0.1, size = 0.8, colour = coldot, alpha = 0.6) +
    geom_point(data = datamaFirstPlot, aes(x= x, y=median),
               size = 2.5,
               #shape = 21,
               fill = coldot,
               colour = coldot,
               na.rm = TRUE)+
    geom_hline(yintercept = lowhorizontal, colour = "grey50")+
    scale_y_continuous(limits = c(miny,maxy))+
    scale_x_continuous(breaks = seq(minx,max(x), by=by))+
    taCont_theme()+
    ylab(ylab)+
    theme(axis.title.x = element_blank(),
          plot.title = element_text(face = "bold"))
  gcen <- if(length(datamaPerc$PercCen) < 1) {gcen} else {
    gcen+geom_text(data = subset(datamaFirstPlot, PercCen >0),
                   aes(label = paste(round(PercCen,0),"%", sep = ""), y = 0),
                   vjust = 1, size = 2.5, colour = censdot)}

  if(sum(!is.na(datamaAll$median))<=4 & plot == TRUE & length(datamaAll$x)>4) {
    return(list(gcen,results))
  } else {
    if(sum(!is.na(datamaAll$median))<=4) {
      return(results)
    } else {

      # moving average of k years - smoother, k set to 3
      datama <- datamaAll %>% filter(!is.na(median))
      tem.zoo<-zoo(datama$median, datama$x)
      m.av<- rollmean(tem.zoo, 3, fill = c(NA, NULL, NA))
      m.av[1]<-datama$median[1]*1/2+datama$median[2]*1/2
      m.av[length(m.av)]<-datama$median[length(m.av)-1]*1/2+datama$median[length(m.av)]*1/2
      datama$m.av<-coredata(m.av)
      mav <- datama$m.av
      y<-datama$median
      x<-datama$x

      #mean or median
      if(ytransformed == TRUE) {
        c(model1 <- median(datafr$y, na.rm=TRUE)/1000,
          model1mean <- exp(mean(datama$mean, na.rm=TRUE))/1000,
          model1median <- exp(median(datama$median, na.rm=TRUE))/1000,
          model1med_5CI <- exp(quantile(boot::boot(datama$median,median.boot,R=R2)$t, 0.05)[[1]])/1000,
          model1med_95CI <- exp(quantile(boot::boot(datama$median,median.boot,R=R2)$t, 0.95)[[1]])/1000)
      } else {
        c(model1 <- median(datafr$y, na.rm=TRUE),
          model1mean <- exp(mean(datama$mean, na.rm=TRUE)),
          model1median <- exp(median(datama$median, na.rm=TRUE)),
          model1med_5CI <- exp(quantile(boot::boot(datama$median,median.boot,R=R2)$t, 0.05)[[1]]),
          model1med_95CI <- exp(quantile(boot::boot(datama$median,median.boot,R=R2)$t, 0.95)[[1]]))
      }
      model1median.log <- median(datama$median, na.rm = TRUE)

      #log-linear regression
      model2 <- lm(y~x)
      m2 <- summary(model2)
      p2 <- predict(model2)
      cip2 <- predict(model2, interval = "confidence") # 95% confidence interval of fit, lwr and upr
      slope <-as.vector(coef(model2)[2])
      CIslope <- confint(model2, "x")*100
      if(ytransformed == TRUE) {
        c(lastyr <- tail(exp(cip2), n=1)[1]/1000,
          lastyrlwr <- tail(exp(cip2), n=1)[2]/1000,
          lastyrupr <- tail(exp(cip2), n=1)[3]/1000)
      } else {
        c(lastyr <- tail(exp(cip2), n=1)[1],
          lastyrlwr <- tail(exp(cip2), n=1)[2],
          lastyrupr <- tail(exp(cip2), n=1)[3])
      }
      plin_m2 <- m2$coefficients[8]

      #log-linear regression of the latest 10 years of data, skips years if no data recorded!
      y10 <- tail(y, n=10); x10 <- tail(x, n= 10)
      if(ytransformed == TRUE) {
        c(model10median <- exp(median(y10, na.rm = TRUE))/1000,
          model10med_5CI <- exp(quantile(boot::boot(y10,median.boot,R=R2)$t, 0.05)[[1]])/1000,
          model10med_95CI <- exp(quantile(boot::boot(y10,median.boot,R=R2)$t, 0.95)[[1]])/1000)
      } else {
        c(model10median <- exp(median(y10, na.rm = TRUE)),
          model10med_5CI <- exp(quantile(boot::boot(y10,median.boot,R=R2)$t, 0.05)[[1]]),
          model10med_95CI <- exp(quantile(boot::boot(y10,median.boot,R=R2)$t, 0.95)[[1]]))
      }
      yrs10_start <- min(x10)
      yrs10_end <- max(x10)
      n_tot_10 <- sum(tail(datama$N, n = 10))
      model10 <- lm(y10~x10)
      m10 <- summary(model10)
      p10 <- predict(model10)
      cip10 <- predict(model10, interval = "confidence") # 95% confidence interval of fit, lwr and upr
      slope10 <-as.vector(coef(model10)[2])
      CIslope10 <- confint(model10, "x10")*100
      if(ytransformed == TRUE) {
        c(lastyr10 <- tail(exp(cip10), n=1)[1]/1000,
          lastyrlwr10 <- tail(exp(cip10), n=1)[2]/1000,
          lastyrupr10 <- tail(exp(cip10), n=1)[3]/1000,
          data10 <- data.frame(x10, "exp.p10" = exp(p10)/1000))
      } else {
        c(lastyr10 <- tail(exp(cip10), n=1)[1],
          lastyrlwr10 <- tail(exp(cip10), n=1)[2],
          lastyrupr10 <- tail(exp(cip10), n=1)[3],
          data10 <- data.frame(x10, "exp.p10" = exp(p10)))
      }
      plast <- if(length(x10)<10){1} else {m10$coefficients[8]}

      #residual sum of squares
      rss1 <- sum((y - model1median.log)^2)
      rss2 <- sum((y - p2)^2)
      #rss3 <- sum((y - p3)^2) # this is for loess smoother
      rss3 <- sum((y - mav)^2) # this is for a k-year moving average smoother

      #degrees of freedom
      df1 <- length(y) - 1
      df2 <- length(y) - 2
      df3 <- (2*length(y) - 1)/3

      #analysis of variance is given by:
      #systematic year effect
      dfsys <- df1-df3
      sssys <- rss1 - rss3
      fsys <- (sssys * df3)/(rss3*dfsys)
      psys <- 1 - pf(fsys, dfsys, df3)
      #non-linearity
      dfnlin <- df2 - df3
      ssnlin <- rss2 - rss3
      fnlin <- (ssnlin * df3)/(rss3*dfnlin)
      pnlin <- 1 - if(is.nan(pf(fnlin, dfnlin, df3))) {0} else {pf(fnlin, dfnlin, df3)}
      #linearity
      dflin <- df1 - df2
      sslin <- rss1 - rss2
      flin <- (sslin * df2)/(rss2 * dflin)
      plin <- 1 - pf(flin, dflin, df2)
      #Error
      dferror <- df3
      sserror <- rss3
      s2error <- rss3/df3

      if(ytransformed == TRUE){
        c(datama$exp.p2 <- exp(p2)/1000,
          datama$exp.mav <- exp(mav)/1000,
          datama$exp.med <- exp(datama$median)/1000,
          datama$med_5CIplot <- exp(datama$med_5CI)/1000,
          datama$med_95CIplot <- exp(datama$med_95CI)/1000)
      } else {
        c(datama$exp.p2 <- exp(p2),
          datama$exp.mav <- exp(mav),
          datama$exp.med <- exp(datama$median),
          datama$med_5CIplot <- exp(datama$med_5CI),
          datama$med_95CIplot <- exp(datama$med_95CI))
      }

      if(ynegative == TRUE){
        c(datama$exp.p2 <- datama$exp.p2*-1,
          datama$exp.mav <- datama$exp.mav*-1,
          datama$exp.med <- datama$exp.med*-1,
          datama$med_5CIplot <- datama$med_5CIplot*-1,
          datama$med_95CIplot <- datama$med_95CIplot*-1)
      } else {
        c(datama$exp.p2 <- datama$exp.p2,
          datama$exp.mav <- datama$exp.mav,
          datama$exp.med <- datama$exp.med,
          datama$med_5CIplot <- datama$med_5CIplot,
          datama$med_95CIplot <- datama$med_95CIplot)
      }

      datama$date<-as.Date.yearmon(datama$x,6)

      q <- slope*100
      q10 <- slope10*100
      philin <- sqrt(rss2/df2) # residuals of regression line
      phinlin <- sqrt(s2error) # residuals of smoothed line
      phi10 <- m10$sigma
      #different power statistics calculated for the time series
      #power for the current time series
      powercurrent <- 1-pf(fdist(length(y)),1,length(y)-2,ncp=delta(philin,q,length(y)))
      #power for the timeseries to detect 5% and 10% slopes
      powerreg5 <- 1-pf(fdist(length(y)),1,length(y)-2,ncp=delta(philin,5,length(y)))
      powerreg10 <- 1-pf(fdist(length(y)),1,length(y)-2,ncp=delta(philin,10,length(y)))
      #power for the current time series phi to detect 5% slope in 10 years
      powerreg510yrs <- 1-pf(fdist(10),1,10-2,ncp=delta(philin,5,10))
      #power for the last 10 yrs
      powercurrentlast <- 1-pf(fdist(length(y10)),1,length(y10)-2,ncp=delta(phi10,q10,length(y10)))
      powerreg5last <- 1-pf(fdist(length(y10)),1,length(y10)-2,ncp=delta(phi10,5,length(y10)))
      powerreg5yrslast <- 1-pf(fdist(10),1,10-2,ncp=delta(phi10,5,10))

      datafr10 <- datafr %>% filter(x > yrs10_start)
      datamaPerc <- datama %>% filter(PercCen >0)

      yrs <- c(paste(yrs_start, "-", yrs_end, sep = ""), paste(yrs10_start, "-", yrs10_end, sep = ""))
      censored <- if(is.null(cen)){c(NA,NA)
      } else {
        c(NADA::censummary(datafr$y, datafr$cen)$all[[3]],
          NADA::censummary(datafr10$y, datafr10$cen)$all[[3]])
      }
      n_yrs <- c(paste(length(datamaAll$x), " (", sum(!is.na(datamaAll$median)), ")", sep = ""), 10)
      n_tot <- c(length(datafr$y), n_tot_10)
      median <- if(ynegative == TRUE){
        c(model1median*-1, model10median*-1)
      } else {
        c(model1median, model10median)
      }
      median_CI_5 <- if(ynegative == TRUE){
        c(model1med_5CI*-1, model10med_5CI*-1)
      } else {
        c(model1med_5CI, model10med_5CI)
      }
      median_CI_95 <- if(ynegative == TRUE){
        c(model1med_95CI*-1, model10med_95CI*-1)
      } else {
        c(model1med_95CI, model10med_95CI)
      }
      annual_change <- if(ynegative == TRUE){c(q*-1, q10*-1)
      } else {
        c(q, q10)
      }
      ac_CI_low <- if(ynegative == TRUE){
        c(CIslope[1]*-1, CIslope10[1]*-1)
      } else {
        c(CIslope[1], CIslope10[1])
      }
      ac_CI_upp <- if(ynegative == TRUE){
        c(CIslope[2]*-1, CIslope10[2]*-1)
      } else {
        c(CIslope[2], CIslope10[2])
      }
      CV_linear <- c(m2$sigma, m10$sigma)
      r2 <- if(ynegative == TRUE){
        c(m2$r.squared*-1, m10$r.squared*-1)
      } else {
        c(m2$r.squared, m10$r.squared)
      }
      p_linear <- c(m2$coefficients[8],m10$coefficients[8])
      lowest_linear_det_change <- c(powerslope(philin,length(y)), powerslope(phi10,length(y10)))
      pow <- c(powercurrent, powercurrentlast)
      adequacy <- c(length(y)/poweryrs(philin, 5), length(y10)/poweryrs(phi10, 5))
      lowest_nonlinear_det_change <- c(powerslope(phinlin,length(y)), NA)
      CV_nonlin <- c(phinlin, NA)
      p_nonlin <- c(pnlin, NA)

      results<-data.frame(n, censored, yrs, n_yrs, n_tot, median, median_CI_5, median_CI_95,
                          annual_change, ac_CI_low, ac_CI_upp,
                          CV_linear, r2, p_linear, lowest_linear_det_change, pow, adequacy,
                          lowest_nonlinear_det_change, CV_nonlin, p_nonlin)
      results[2,]<-if(length(datama$x)<10){NA} else {results[2,]}

      # rules for the plot
      lines <- 0.6
      coldot <- "black"
      colbar <- "grey10"
      colline <- "violetred1"
      collinesm <- "springgreen"
      colline10 <- "royalblue"
      censdot <- "violetred4"
      miny <- if(ynegative == TRUE){min(datafr$yplotplot)-1} else {0}
      maxy <- if(ynegative == TRUE){max(datafr$yplotplot)+1} else {max(datafr$yplot)}
      maxypub <- if(ynegative == TRUE){abs(maxy)} else {maxy*2.5}
      minx <-if((min(datafr$x)%%2)==0){min(datafr$x)}else{min(datafr$x)-1}
      by <- if((max(datafr$x)-minx)<10){1}else{4}
      lowhorizontal <- if(ynegative == TRUE){min(datafr$yplotplot)-1} else {0}

      ken <- Kendall(x,y)
      MK <- MannKendall(y)
      kentau <- if(ynegative == TRUE){ken$tau[1]*-1} else {ken$tau[1]}
      kenp <- ken$sl[1]

      model1median <- if(ynegative == TRUE){model1median*-1} else {model1median}

      message("psys=", psys, "\n")
      message("plin=", plin, "\n")
      message("pnlin=", pnlin, "\n")
      message("plin(m2)=", plin_m2, "\n")

      a1 <- paste("n(tot) =",length(datafr$y),
                  ", n(yrs) =",length(datamaAll$x),
                  " (", sum(!is.na(datamaAll$median)), ")" , sep = "")
      a2 <- if(ynegative == TRUE){
        paste("y(",tail(datama$x, n=1),") = ", format(lastyr*-1, digits = 2),
              " (", format(lastyrlwr*-1, digits = 2),
              ", ", format(lastyrupr*-1, digits = 2), ")", sep="")
      } else {
        paste("y(",tail(datama$x, n=1),") = ", format(lastyr, digits = 2),
              " (", format(lastyrlwr, digits = 2),
              ", ", format(lastyrupr, digits = 2), ")", sep="")
      }
      a3 <- paste("median = ",format(model1median, digits = 2),
                  " (", format(median_CI_5[1], digits = 2),
                  ", ", format(median_CI_95[1], digits = 2),")", sep = "")
      a4 <- paste("slope = ", format(annual_change[1], digits = 2),"%",
                  " (", format(ac_CI_low[1], digits = 2),
                  ", ", format(ac_CI_upp[1], digits = 2),")", sep = "")
      a5 <- paste("SD(lr) =",format(m2$sigma, digits = 2), ", ",
                  powerslope(philin,length(y)),"%, ",
                  format(poweryrs(philin, 5), digits = 1)," yrs", sep="")
      a61 <- paste("Power (5%, 10 yrs) = ", format(powerreg5, digits = 2),
                   ", ", format(powerreg510yrs, digits = 2),
                   ", ", powerslope(philin,10),"%", sep = "")
      a6 <- paste("Adequacy = ", format(length(y)/poweryrs(philin, 5), digits = 2),
                  ", Power = ", format(powercurrent, digits = 2), sep = "")
      a7 <- paste("r2 = ", format(r2[1], digits = 2),
                  ", p < ", format(m2$coefficients[8], digits = 2), sep ="")
      a8 <- paste("tau = ", format(kentau, digits = 2),
                  ", p < ", format(kenp, digits = 2), sep = "")
      a9 <- paste("SD(sm) = ", format(phinlin, digits = 2),
                  ", p < ", format(pnlin, digits = 2),
                  ", ",powerslope(phinlin,10),"%", sep = "")
      a10 <- "- - -"
      a11 <- paste("slope = ", format(annual_change[2], digits = 2),"%",
                   " (", format(ac_CI_low[2], digits = 2),
                   ", ", format(ac_CI_upp[2], digits = 2),")", sep = "")
      a12 <- paste("SD(lr) =",format(m10$sigma, digits = 2),
                   ", ", powerslope(phi10,length(y10)),"%, ",
                   format(poweryrs(phi10, 5), digits = 1)," yrs", sep="")
      a121 <- paste("Adequacy = ", format(length(y10)/poweryrs(phi10, 5), digits = 2),
                    ", Power = ", format(powercurrentlast, digits = 2), sep = "")
      a13 <- paste("Power (5%, 10 yrs) = ", format(powerreg5last, digits = 2),
                   ", ", format(powerreg5yrslast, digits = 2),
                   ", ", powerslope(phi10,10),"%", sep = "")
      a14 <- paste("r2 = ", format(r2[2], digits = 2),
                   ",p < ", format(m10$coefficients[8], digits = 2), sep ="")
      z <- if(length(x10)==10 & length(x) > 10) {
        paste(a1,a2,a3,a4,a5,a6,a61,a7,a8,a9,a10,a11,a12,a121,a13,a14,sep = "\n")
      } else {
        paste(a1,a2,a3,a4,a5,a6,a61,a7,a8,a9,sep = "\n")
      }
      my_grob <- grobTree(textGrob(z,x=0.05, y=0.95, hjust=0, vjust = 1,
                                   gp=gpar(col="grey10", fontsize = 8)))

      g<-ggplot(datama, aes(x=x, y=exp.med)) +
        geom_jitter(data=datafr %>% filter(YearExcluded == FALSE, LOQLrgMean == FALSE),
                    aes(x=x, y=yplotplot),
                    position = position_jitter(0.1), size = 1,
                    colour= coldot,
                    alpha = 0.25)+
        geom_jitter(data=datafr %>% filter(YearExcluded == TRUE),
                    aes(x=x, y=yplotplot),
                    position = position_jitter(0.1), size = 1,
                    colour= censdot,
                    alpha = 0.25)+
        geom_jitter(data=datafr %>% filter(LOQLrgMean == TRUE),
                    aes(x=x, y=yplotplot, colour = cen),
                    position = position_jitter(0.1), size = 1,
                    colour= censdot,
                    shape = 8,
                    alpha = 0.5)+
        geom_errorbar(data=datama, aes(ymin=med_5CIplot, ymax=med_95CIplot),
                      width = 0.1, size = 0.8, colour = coldot, alpha = 0.6) +
        geom_point(data = datama, aes(x= x, y=exp.med),
                   size = 2.5,
                   #shape = 21,
                   fill = coldot,
                   colour = coldot,
                   na.rm = TRUE)+
        geom_hline(yintercept = model1median, linetype = 2, colour= "gray10") +
        geom_hline(yintercept = lowhorizontal, colour = "grey50")+
        scale_y_continuous(limits = c(miny,maxy))+
        scale_x_continuous(breaks = seq(minx,max(x), by=by))+
        taCont_theme()+
        ylab(ylab)+
        theme(axis.title.x = element_blank(),
              plot.title = element_text(face = "bold"))
      g <- if(plin_m2<0.05){
        g+geom_line(data= datama, aes(x=x, y=exp.p2), colour=colline, size =lines)
      } else {g}
      g <- if(plin_m2>0.05 & plin_m2<0.1 & palmost == TRUE){
        g+geom_line(data= datama, aes(x=x, y=exp.p2), colour=colline, size =lines, linetype = 2)
      } else {g}
      g <- if(pnlin<0.05){
        g+geom_line(data= datama, aes(x=x, y=exp.mav), colour=collinesm, size =lines)
      } else {g}
      g <- if(pnlin>0.05 & pnlin<0.1 & palmost == TRUE){
        g+geom_line(data= datama, aes(x=x, y=exp.mav), colour=collinesm, linetype = 5, size =lines)
      } else {g}
      g <- if(plast<0.05 & length(x) > 10){
        g+geom_line(data = data10, aes(x=x10, y=exp.p10), colour=colline10, size =lines)
      } else {g}
      g <- if(plast>0.05 & plast<0.1 & length(x) > 10 & palmost == TRUE){
        g+geom_line(data = data10, aes(x=x10, y=exp.p10), colour=colline10, linetype = 2, size =lines)
      } else {g}
      g <- if(length(datamaPerc$PercCen) < 1) {g
      } else {
        g+geom_text(data = subset(datamaAll, PercCen >0),
                    aes(label = paste(round(PercCen,0),"%", sep = ""), y = 0),
                    vjust = 1, size = 2.5, colour = censdot)
      }
      g <- if(pub == TRUE){
        g+scale_y_continuous(limits = c(miny,maxypub))+annotation_custom(my_grob)
      } else {g}

      #g <- if(pnlin<2){g+geom_line(data= datama, aes(x=x, y=exp.p3.), colour=colline, linetype = 2)} else {g}

      if(plot == TRUE & onlyRes == FALSE) {
        message(z)
        plot(g)
        return(list(g, results, model2, model10))
      } else if (plot == FALSE & onlyRes == FALSE) {
        return(list(g, results, model2, model10))
      } else {
        return(results)
      }
    }
  }
}



