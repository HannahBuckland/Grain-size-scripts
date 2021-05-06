# R code for fitting distributions to binned grain size data

# Script to fit normal distributions to phi binned data
# Author: Hannah Buckland, School of Earth Sciences University of Bristol

# Input data must be in tab separated text file with three columns with following headers:

#    | mm        | phi       | PDF      |
#    | 0.00024   | 12        | 0        |
#    | ...       | ...       | ...      |
#    | 32         | -5       | 0        |

# The first and final row must contain no wt% in the PDF column (expand range of bins if necessary)
# The data should be ordered low to high in mm
# The PDF column must sum to 100%
# An example input file is found in the examples folder


library(fitdistrplus)
library(mixR)
library(tidyverse)
library(gridExtra)
library(patchwork)


##############################################
# USER INPUTS

GSDfile <- "examples/AP1_xjet_xcminhalfphi.txt" # path to grain size distribution file (see format guidelines above)
GSDname <- "AP1_xcmin" # name for script to name outputs
bimodal <- "FALSE" # TRUE or FALSE

# Shouldn't need to change below this line (but feel free to personalise if you would like)
######################################################
# 

SIMDATAFUN <- function(GSDfile) {
  
  # Read in and manipulate data
  
  gsd <- read.table(GSDfile, sep="", header = TRUE)
  gsd$freq <- gsd$PDF*1000
  
  phi_bins <- gsd$phi
  phi_bins <- append(phi_bins, values = tail(phi_bins,1)-(phi_bins[1]-phi_bins[2]), after = length(phi_bins))
  
  mm_bins <- gsd$mm
  mm_bins <- append(mm_bins, values = 2^-(tail(phi_bins,1)), after = length(mm_bins))
  
  # Set midpoint of grain size bins
  
  gsd$midphi <- gsd$phi - (abs(gsd$phi[2]-gsd$phi[1])/2)
  
  gsd$midmm <- 2^-gsd$midphi
  
  # Manipulate data into artificial distribution in phi scale using wt% PDF
  # Simulates ~100000 grain size measurements that produce histogram that matches measured PDF
  # The frequency in each bin is a factor of the the wt or vol% in each bin 
  # The range of simulated values is uniformly distributed between the max and min of each grain size bin
  output <- list()
  for (i in 1:length(gsd$phi)) {
    output[[i]] <- runif((as.numeric(gsd$freq[i])), min=phi_bins[i + 1], max=phi_bins[i])
  }
  
  fout <- list(gsd=gsd,phi_bins = phi_bins, mm_bins = mm_bins, output=output)
  return(fout)
  
} # function to simulate GSD data

simulate <- SIMDATAFUN(GSDfile) # applying function to GSD file

simdata <- unlist(simulate$output) # getting into long format
simdatared <- sample(simdata,size = 10000,replace=FALSE) # reducing sample size for 'mixfit' functionality

phi_bins <- simulate$phi_bins # add to environment

mm_bins <- simulate$mm_bins # add to environment

simdf <- data.frame(size=simdata)

gsd <- simulate$gsd # add to environment

# Some global plotting parameters
micron_plotting <- data.frame(breaks= -log2(c(10e-3,100e-3,1000e-3)), labels = c(10,100,1000))
phi_plotting <- data.frame(breaks = seq(-3,12))

NORMFITFUN <- function(gsd){
  
  gsd <- gsd
  
  # Calculate Method of Moments mean and standard deviation
  
  log_mean <- sum(gsd$PDF*gsd$midphi)/100
  log_sd <- sqrt(sum(gsd$PDF*(gsd$midphi-log_mean)^2)/100)
  
  # Graphical Folk and Ward method using linear interpolation of cumulative plot
  
  gsd <- gsd[order(gsd$phi),]
  gsd$cum <- cumsum(gsd$PDF)
  gsd <- gsd[order(-gsd$phi),]
  
  phi5 <- approx(x=gsd$cum, y=gsd$phi, xout=5, method = "linear", ties = mean)
  phi16 <- approx(gsd$cum, gsd$phi, xout=16, method = "linear", ties = mean)
  phi50 <- approx(gsd$cum, gsd$phi, xout=50, method = "linear", ties = mean)
  phi84 <- approx(gsd$cum, gsd$phi, xout=84, method = "linear", ties = mean)
  phi95 <- approx(gsd$cum, gsd$phi, xout=95, method = "linear", ties = mean)
  
  fwmean <- (phi16[[2]] + phi50[[2]] + phi84[[2]])/3
  fwsd <- ((phi84[[2]] - phi16[[2]])/4) + ((phi95[[2]] - phi5[[2]])/6.6)
  
  inmean <- 0.5*(phi16[[2]]+phi84[[2]])
  insd <- 0.5*(phi84[[2]]-phi16[[2]])
  
  dens <- hist(simdata,breaks = phi_bins,plot=FALSE)
  dens <- dens$density
  
  gsd$dens <- rev(dens)
  
  simuplt <-
    ggplot(data = simdf) +
    geom_histogram(aes(x=size, fill = "Simulated Data"), breaks = phi_bins, colour = "grey50") +
    geom_line(data = gsd,
              aes(x=midphi, y = freq)) +
    geom_point(data = gsd,
               aes(x=midphi, y = freq, colour = "GSD")) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "grey60") +
    ylab("Frequency/Counts") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~./1000,
                                           name = "Measured (%)")) +
    scale_x_reverse(breaks = phi_plotting$breaks, 
                    labels = phi_plotting$breaks,
                    name = "Grain size (phi)",
                    sec.axis = sec_axis(trans = ~.*1, 
                                        name = "Grain Size (µm)", 
                                        breaks = micron_plotting$breaks,
                                        labels = micron_plotting$labels)) +
    ggtitle("Simulated vs measured") +
    theme_bw() +
    theme(legend.position = c(0.15,0.7),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black"),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
  
  norm_fit <- fitdistrplus::fitdist(simdata,"norm") # fits normal distribution, setups up list with parameters
  
  x <-seq(min(phi_bins), max(phi_bins),length.out = length(simdata)) 
  
  nfit <- dnorm(x,
                mean = (norm_fit$estimate['mean']),
                sd = (norm_fit$estimate['sd']))
  ncheck <- 2*(pnorm(norm_fit$estimate['mean'], 
                     mean = norm_fit$estimate['mean'], 
                     sd = norm_fit$estimate['sd']))
  
  n_mean = norm_fit$estimate['mean']
  n_sd = norm_fit$estimate['sd']
  
  phifit <- dnorm(x,
                  log_mean,
                  log_sd) # use phi mean and sd logarithmic Moment Method
  
  fwfit <- dnorm(x,
                 fwmean,
                 fwsd) # use Folk and Ward mean and sd 
  
  infit <- dnorm(x,
                 inmean,
                 insd) # use Inman mean and sd 
  
  nfit_param <- c(n_mean,n_sd) # for adding to table of parameters
  logM <- c(log_mean,log_sd)
  FW <- c(fwmean,fwsd)
  IN <- c(inmean,insd)
  
  normal_fits <- data.frame(x=x,nfit=nfit,phifit = phifit, fwfit= fwfit,infit=infit)
  normal_pivot <- pivot_longer(normal_fits, cols = c(nfit,phifit,fwfit,infit))
  
  mean_fits <- data.frame(nfit = n_mean, phifit = log_mean, fwfit = fwmean, infit = inmean)
  mean_pivot <- pivot_longer(mean_fits, cols = c(nfit,phifit,fwfit,infit))
  
  fit_labels <- c("Folk and Ward","Inman","fitdist","Method of Moments")
  fit_colours <- c("#264653","#2a9d8f","#e9c46a","#e63946")
  
  nfitplt <-
    ggplot(data = simdf) +
    geom_histogram(aes(x=size, y =..density..,fill = "Simulated Data"), breaks = phi_bins, colour = "grey50") +
    geom_line(data = gsd,
              aes(x=midphi, y = PDF/50)) +
    geom_point(data = gsd,
               aes(x=midphi, y = PDF/50), colour = "grey30") +
    geom_line(data = normal_pivot,
              aes(x=x,y=value,colour=name)) +
    geom_vline(data = mean_pivot,
               aes(xintercept = value, colour = name), linetype = 2) +
    scale_fill_manual(values = "grey80") +
    scale_colour_manual(values = fit_colours,
                        labels = fit_labels) +
    ylab("Density") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*50,
                                           name = "Measured (%)")) +
    scale_x_reverse(breaks = phi_plotting$breaks, 
                    labels = phi_plotting$breaks,
                    name = "Grain size (phi)",
                    sec.axis = sec_axis(trans = ~.*1, 
                                        name = "Median Grain Size (µm)", 
                                        breaks = micron_plotting$breaks,
                                        labels = micron_plotting$labels)) +
    ggtitle("Fitted normal distributions") +
    theme_bw() +
    theme(legend.position = c(0.15,0.65),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black"),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
  
  # Setup dataframe of parameters to print
  
  params <- as.data.frame(cbind(logM,FW,IN,nfit_param))
  colnames(params) <- c("Method of Moments (phi)", "Folk and Ward (phi)","Inman (phi)","Normal Fit (phi)")
  rownames(params) <- c("Mean","SD")
  
  params <- round(params,3)
  fout <- list(x=x,nfitparams=params, simuplot = simuplt, nfitplot = nfitplt)
  return(fout)
  
} # function to fit unimodal normal distributions 

normfit <- NORMFITFUN(gsd) # applying normal fitting functions
x <- normfit$x # add x to environment


if (bimodal == TRUE) {
  print("Results in bimodal_results and bimodal_plots")
  
  BINORMFITFUN <- function(simdatared) {
    
    mixmdl <- mixfit(simdatared, ncomp = 2, family = "norm")
    
    nfit1 <- mixmdl$pi[1]*dnorm(x,mixmdl$mu[1],mixmdl$sd[1])
    nfit2 <- mixmdl$pi[2]*dnorm(x,mixmdl$mu[2],mixmdl$sd[2])
    nsumfit <- nfit1+nfit2
    
    bimodal_fit <- data.frame(x=x,nfit1=nfit1,nfit2=nfit2,nsumfit=nsumfit)
    bimodal_pivot <- pivot_longer(bimodal_fit,cols =c(nfit1,nfit2,nsumfit))
    
    bimodalfitplt <-
      ggplot(data = simdf) +
      geom_histogram(aes(x=size, y =..density..,fill = "Simulated Data"), breaks = phi_bins, colour = "grey50") +
      geom_line(data = gsd,
                aes(x=midphi, y = PDF/50)) +
      geom_point(data = gsd,
                 aes(x=midphi, y = PDF/50), colour = "grey30") +
      geom_line(data = bimodal_pivot,
                aes(x=x,y=value,colour=name, size = name)) +
      scale_fill_manual(values = "grey80") +
      scale_colour_manual(values = c("darkcyan","coral","grey10")) +
      scale_size_manual(values = c(1,1,2)) + 
      ylab("Density") +
      scale_y_continuous(sec.axis = sec_axis(trans = ~.*50,
                                             name = "Measured (%)")) +
      scale_x_reverse(breaks = phi_plotting$breaks, 
                      labels = phi_plotting$breaks,
                      name = "Grain size (phi)",
                      sec.axis = sec_axis(trans = ~.*1, 
                                          name = "Median Grain Size (µm)", 
                                          breaks = micron_plotting$breaks,
                                          labels = micron_plotting$labels)) +
      ggtitle("Bimodal normal distributions") +
      theme_bw() +
      theme(legend.position = c(0.15,0.65),
            legend.title = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black"),
            plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
            plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
    
    nmixcheck <- auc(x,nsumfit)
    
    fit1param <- c(mixmdl$pi[1],mixmdl$mu[1],mixmdl$sd[1])
    fit2param <- c(mixmdl$pi[2],mixmdl$mu[2],mixmdl$sd[2])
    
    bimod <- data.frame(cbind(fit1param,fit2param))
    colnames(bimod) <- c("nfit1","nfit2")
    rownames(bimod) <- c("Prop", "Mean","SD")
    
    bimod[1:2] <- round(bimod[1:2], 3)
    
    fout <- list(bimodalfitplt = bimodalfitplt,bimodalparams=bimod)
    
  } # function to fit bimodal normal distribution
  
  bimodalfit <- BINORMFITFUN(simdatared)
  
  bimodal_results <- list(unimodalresults = normfit$nfitparams,bimodalresults = bimodalfit$bimodalparams)
  bimodal_table <- bind_rows(bimodal_results)
  rownames(bimodal_table) <- c("UniMean","UniSd","Prop","BimodMean","BimodSd")
  bimodal_plots <- list(unimodal_simulated = normfit$simuplot, unimodal_plot = normfit$nfitplot, bimodal_plot = bimodalfit$bimodalfitplt)
  
  # Build output plot
  unimodtab <- tableGrob(bimodal_results$unimodalresults)
  bimodtab <- tableGrob(bimodal_results$bimodalresults)
  
  finalplot <- (bimodal_plots$unimodal_simulated + bimodal_plots$unimodal_plot + bimodal_plots$bimodal_plot) +
    plot_layout(ncol = 3, nrow = 2, widths = c(1,1,1)) +
    plot_annotation(title = paste("Distribution plots and data for bimodal", GSDname),
                    theme = theme(plot.title = element_text(hjust =0.5,
                                                            size = 14,
                                                            margin = margin(1,0,0.5,0,"cm")))) +
    inset_element(bimodtab,1,-1,0,0) +
    inset_element(unimodtab,-1.5,-1,0,0)
  
  ggsave(filename = sprintf("plots/%s.png",GSDname), finalplot, width = 20, height = 9)
  write.table(bimodal_table,file =sprintf("output/%s.txt",GSDname),sep = "\t")
  
} else {
  print("Results in unimodal_results and unimodal_plots")
  
  # Built output plot
  unimodal_results <- list(unimodalresults = normfit$nfitparams)
  unimodal_table <- bind_rows(unimodal_results)
  unimodal_plots <- list(unimodal_simulated = normfit$simuplot, unimodal_plot = normfit$nfitplot)
  
  unimodtab <- tableGrob(unimodal_results$unimodalresults)
  
  finalplot <- unimodal_plots$unimodal_simulated + unimodal_plots$unimodal_plot + 
    plot_layout(ncol = 2, nrow = 2, widths = c(1,1)) +
    plot_annotation(title = paste("Distribution plots for unimodal", GSDname),
                    theme = theme(plot.title = element_text(hjust =0.5,
                                                            size = 12,
                                                            margin = margin(0,0,0.5,0,"cm")))) + 
    inset_element(unimodtab,1,-1,0,0)
  
  ggsave(filename = sprintf("plots/%s.png",GSDname), finalplot, width = 15, height = 9)
  write.table(unimodal_table,file =sprintf("output/%s.txt",GSDname),sep = "\t")
}

######################################################
# OUTPUTS ARE SAVED INTO "outputs/" and "plots/"
