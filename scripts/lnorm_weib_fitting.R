# R code for fitting distributions to binned grain size data

# Script to fit log-normal and weibull distributions to mm binned data
# Author: Hannah Buckland, School of Earth Sciences University of Bristol

# Input data must be in tab separated text file with three columns with following headers:

#    | mm       | phi         | PDF     |
#    | 0        | 6.085       | 0.805   |
#    | ...      | ...         | ...     |
#    | 0.255    | 1.971       | 0       |

# The final row must contain no wt% in the PDF column (expand range of bins if necessary)
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

GSDfile <- "examples/Starbuck_xjet_xcmin5um.txt" # path to grain size distribution file (see format guidelines above)
GSDname <- "Starbuck_xjet_xcmin5um" # name for script to name outputs
bimodal <- "TRUE" # TRUE or FALSE
method_min <- 0.0008 # set the minimum limit of the grain size method in mm

# Shouldn't need to change below this line (but feel free to personalise if you would like)
######################################################
# 

SIMDATAFUN <- function(GSDfile) {
  
  # Read in and manipulate data
  
  gsd <- read.table(GSDfile, sep="", header = TRUE)
  gsd$freq <- round(gsd$PDF*100,0)
  
  mm_bins <- gsd$mm
  mm_bins[1] <- method_min # set resolution of measurement 
  mm_bins <- append(mm_bins, values = tail(mm_bins,1)-(mm_bins[1]-mm_bins[2]), after = length(mm_bins))
  
  
  # Set midpoint of grain size bins
  
  gsd$midmm<- gsd$mm + (abs(gsd$mm[3]-gsd$mm[2])/2)
  
  gsd$midphi <- -log2((gsd$midmm))
  
  # Manipulate data into artificial distribution in phi scale using wt% PDF
  # Simulates ~100000 grain size measurements that produce histogram that matches measured PDF
  # The frequency in each bin is a factor of the the wt or vol% in each bin 
  # The range of simulated values is uniformly distributed between the max and min of each grain size bin
  output <- list()
  for (i in 1:length(gsd$phi)) {
    output[[i]] <- runif((as.numeric(gsd$freq[i])), 
                         min=mm_bins[i], 
                         max=mm_bins[i+1])
  }
  
  fout <- list(gsd=gsd, mm_bins = mm_bins, output=output)
  return(fout)
  
} # function to simulate GSD data

simulate <- SIMDATAFUN(GSDfile) # applying function to GSD file

simdata <- unlist(simulate$output) # getting into long format

mm_bins <- simulate$mm_bins # add to environment

mm_bin_fact <- 100*(mm_bins[3]-mm_bins[2])

simdf <- data.frame(size=simdata)

gsd <- simulate$gsd # add to environment

# Some global plotting parameters
micron_plotting <- data.frame(breaks= seq(0,1,0.05), labels = seq(0,1000,50))
phi_plotting <- data.frame(breaks = c(2^-1,2^-2,2^-3,2^-4, 2^-5, 2^-6), labels = c(1,2,3,4,5,6))

DISTFITFUN <- function(gsd){
  
  gsd <- gsd
  
  simuplt <-
    ggplot(data = simdf) +
    geom_histogram(aes(x=size, fill = "Simulated Data"), breaks = mm_bins, colour = "grey50") +
    geom_line(data = gsd,
              aes(x=midmm, y = freq)) +
    geom_point(data = gsd,
               aes(x=midmm, y = freq, colour = "GSD")) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "grey60") +
    ylab("Frequency/Counts") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~./100,
                                           name = "Measured (%)")) +
    scale_x_continuous(breaks = micron_plotting$breaks, 
                       labels = micron_plotting$labels,
                       name = "Grain size (µm)",
                       sec.axis = sec_axis(trans = ~.*1, 
                                           name = "Grain Size (phi)", 
                                           breaks = phi_plotting$breaks,
                                           labels = phi_plotting$labels)) +
    ggtitle("Simulated vs measured") +
    theme_bw() +
    theme(legend.position = c(0.85,0.7),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black"),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
  
  lnorm_fit <- fitdistrplus::fitdist(simdata,"lnorm") # fits log-normal distribution, setups up list with parameters
  weib_fit <- fitdistrplus::fitdist(simdata,"weibull") # fits Weibull distribution, setups up list with parameters
  
  x <-seq(0, max(mm_bins),length.out = length(simdata))
  
  # Log-normal distribution fitting and parameter calculations
  
  lnfit <- dlnorm(x,
                  mean = (lnorm_fit$estimate['meanlog']),
                  sd = (lnorm_fit$estimate['sdlog']))
  ncheck <- 2*(plnorm(lnorm_fit$estimate['meanlog'],
                      mean = lnorm_fit$estimate['meanlog'],
                      sd = lnorm_fit$estimate['sdlog']))
  
  ln_mu_dash <- as.numeric(lnorm_fit$estimate['meanlog'])
  ln_sig_dash <- as.numeric(lnorm_fit$estimate['sdlog'])
  
  ln_median <- as.numeric(exp(lnorm_fit$estimate['meanlog']))
  ln_mean <- as.numeric(exp(lnorm_fit$estimate['meanlog'] + 0.5*(lnorm_fit$estimate['sdlog'])^2))
  ln_mod <- as.numeric(exp(lnorm_fit$estimate['meanlog'] - (lnorm_fit$estimate['sdlog'])^2))
  
  # Weibull distribution fitting and parameter calculations
  
  wbfit <- dweibull(x,
                    shape = weib_fit$estimate['shape'],
                    scale = weib_fit$estimate['scale'])
  
  wb_shape <- as.numeric(weib_fit$estimate['shape'])
  wb_scale <- as.numeric(weib_fit$estimate['scale'])
  
  wb_median <- wb_scale*(log(2))^(1/wb_shape)
  wb_mean <- wb_scale*gamma(1+1/wb_shape)
  wb_mod <- wb_scale*((wb_shape-1)/wb_shape)^(1/wb_shape)
  
  # Assign parameters
  
  lnfit_param <- c(ln_mu_dash,ln_sig_dash,ln_median,ln_mean,ln_mod) # for adding to table of parameters
  wbfit_param <- c(wb_shape,wb_scale,wb_median,wb_mean,wb_mod) # for adding to table of parameters
  
  dist_fits <- data.frame(x=x,lnfit=lnfit,wbfit=wbfit)
  dist_pivot <- pivot_longer(dist_fits, cols = c(lnfit,wbfit))
  
  stat_fits <- data.frame(id = c("median","mean","mode"),lnfit = lnfit_param[3:5],wbfit = wbfit_param[3:5])
  stat_pivot <- pivot_longer(stat_fits, cols = c(lnfit,wbfit))
  
  
  fitdistplt <-
    ggplot(data = simdf) +
    geom_histogram(aes(x=size, y =..density..,fill = "Simulated Data"), breaks = mm_bins, colour = "grey50") +
    geom_line(data = gsd,
              aes(x=midmm, y = PDF/mm_bin_fact)) +
    geom_point(data = gsd,
               aes(x=midmm, y = PDF/mm_bin_fact), colour = "grey30") +
    geom_line(data = dist_pivot,
              aes(x=x,y=value,colour=name)) +
    geom_vline(data = stat_pivot,
               aes(xintercept = value, colour =name, linetype = id)) +
    scale_fill_manual(values = "grey80") +
    scale_colour_manual(values = c("#264653","#e63946")) +
    ylab("Density") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*mm_bin_fact,
                                           name = "Measured (%)")) +
    scale_x_continuous(breaks = micron_plotting$breaks, 
                       labels = micron_plotting$labels,
                       name = "Grain size (µm)",
                       sec.axis = sec_axis(trans = ~.*1, 
                                           name = "Grain Size (phi)", 
                                           breaks = phi_plotting$breaks,
                                           labels = phi_plotting$labels)) +
    ggtitle("Fitted log-normal distributions") +
    theme_bw() +
    theme(legend.position = c(0.85,0.60),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black"),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
  
  # Setup dataframe of parameters to print
  
  params <- data.frame(lnfit_param,wbfit_param)
  colnames(params) <- c("Log-normal", "Weibull")
  rownames(params) <- c("mu' or shape","sigma' or scale","median (mm)","mean (mm)","mode (mm)")
  
  params <- round(params,3)
  fout <- list(x=x,distparams=params, simuplot = simuplt, fitdistplot = fitdistplt)
  return(fout)
  
} # function to fit unimodal normal distributions 

fitdist_res <- DISTFITFUN(gsd) # applying normal fitting functions
x <- fitdist_res$x # add x to environment


if (bimodal == TRUE) {
  print("Results in bimodal_results and bimodal_plots")
  
  BINORMFITFUN <- function(simdatared) {
    
    lnmix <- mixfit(simdata, ncomp = 2, family = "lnorm")
    weibmix <- mixfit(simdata, ncomp = 2, family = "weibull")
    
    lnfit1 <- lnmix$pi[1]*dlnorm(x,lnmix$mulog[1],lnmix$sdlog[1])
    lnfit2 <- lnmix$pi[2]*dlnorm(x,lnmix$mulog[2],lnmix$sdlog[2])
    lnsumfit <- lnfit1+lnfit2
    
    weibfit1 <- weibmix$pi[1]*dweibull(x,shape=weibmix$k[1],scale=weibmix$lambda[1])
    weibfit2 <- weibmix$pi[2]*dweibull(x,shape=weibmix$k[2],scale=weibmix$lambda[2])
    weibsumfit <- weibfit1+weibfit2
    
    bimodal_fit <- data.frame(x=x,lnfit1=lnfit1,lnfit2=lnfit2,lnsumfit=lnsumfit,weibfit1=weibfit1,weibfit2=weibfit2,weibsumfit=weibsumfit)
    bimodal_pivot <- pivot_longer(bimodal_fit,cols =c(lnfit1,lnfit2,lnsumfit,weibfit1,weibfit2,weibsumfit))
    
    bimodalfitplt <-
      ggplot(data = simdf) +
      geom_histogram(aes(x=size, y =..density..,fill = "Simulated Data"), breaks = mm_bins, colour = "grey50") +
      geom_line(data = gsd,
                aes(x=midmm, y = PDF/mm_bin_fact)) +
      geom_point(data = gsd,
                 aes(x=midmm, y = PDF/mm_bin_fact), colour = "grey30") +
      geom_line(data = bimodal_pivot,
                aes(x=x,y=value,colour=name, size = name), alpha = 0.6) +
      scale_fill_manual(values = "grey80") +
      scale_colour_manual(values = c("#264653","#2a9d8f","grey10","#e9c46a","#e63946","grey20")) +
      scale_size_manual(values = c(1,1,2,1,1,2)) + 
      ylab("Density") +
      scale_y_continuous(sec.axis = sec_axis(trans = ~.*mm_bin_fact,
                                             name = "Measured (%)")) +
      scale_x_continuous(breaks = micron_plotting$breaks, 
                         labels = micron_plotting$labels,
                         name = "Grain size (µm)",
                         sec.axis = sec_axis(trans = ~.*1, 
                                             name = "Grain Size (phi)", 
                                             breaks = phi_plotting$breaks,
                                             labels = phi_plotting$labels)) +
      ggtitle("Bimodal log-normal distributions") +
      theme_bw() +
      theme(legend.position = c(0.85,0.65),
            legend.title = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black"),
            plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
            plot.title = element_text(margin = margin(0,0,0.5,0,"cm")))
    
    
    lnfit1_mod <- exp(lnmix$mulog[1] - lnmix$sdlog[1]^2)
    lnfit2_mod <- exp(lnmix$mulog[2] - lnmix$sdlog[2]^2)
    
    lnfit1_mean <- exp(lnmix$mulog[1] + 0.5*(lnmix$sdlog[1]^2))
    lnfit2_mean <- exp(lnmix$mulog[2] + 0.5*(lnmix$sdlog[2]^2))
    
    wbfit1median <- weibmix$lambda[1]*(log(2))^(1/weibmix$k[1])
    wbfit1mod <- weibmix$lambda[1]*((weibmix$k[1]-1)/weibmix$k[1])^(1/weibmix$k[1])
    
    wbfit2median <- weibmix$lambda[2]*(log(2))^(1/weibmix$k[2])
    wbfit2mod <- weibmix$lambda[2]*((weibmix$k[2]-1)/weibmix$k[2])^(1/weibmix$k[2])
    
    lnfit1param <- c(lnmix$pi[1],lnmix$mulog[1],lnmix$sdlog[1],exp(lnmix$mulog[1]),lnfit1_mean,lnfit1_mod)
    lnfit2param <- c(lnmix$pi[2],lnmix$mulog[2],lnmix$sdlog[2],exp(lnmix$mulog[2]),lnfit2_mean,lnfit2_mod)
    
    wbfit1param <- c(weibmix$pi[1],weibmix$k[1],weibmix$lambda[1],wbfit1median,weibmix$mu[1],wbfit1mod)
    wbfit2param <- c(weibmix$pi[2],weibmix$k[2],weibmix$lambda[2],wbfit2median,weibmix$mu[2],wbfit2mod)
    
    
    bimod <- data.frame(cbind(lnfit1param,lnfit2param,wbfit1param,wbfit2param))
    colnames(bimod) <- c("lnfit1","lnfit2","weibfit1","weibfit2")
    rownames(bimod) <- c("Prop","mu' or shape","sigma' or scale", "median (mm)","mean (mm)","mode (mm)")
    
    bimod <- round(bimod, 3)
    
    fout <- list(bimodalfitplt = bimodalfitplt,bimodalparams=bimod)
    
  } # function to fit bimodal normal distribution
  
  bimodalfit <- BINORMFITFUN(simdata)
  
  bimodal_results <- list(unimodalresults = fitdist_res$distparams,bimodalresults = bimodalfit$bimodalparams)
  bimodal_table <- bind_rows(bimodal_results)
  rownames(bimodal_table) <- c("Uni_mu'_or_shape","Uni_sigma'_or_scaale","Uni_median (mm)","Uni_mean (mm)","Uni_mode (mm)","Bimod_Prop","Bimod_mu'_or_shape","Bimod_sigma'_or_scale", "Bimod_median (mm)","Bimod_mean (mm)","Bimod_mode (mm)")
  bimodal_plots <- list(unimodal_simulated = fitdist_res$simuplot, unimodal_plot = fitdist_res$fitdistplot, bimodal_plot = bimodalfit$bimodalfitplt)
  
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
  write.table(bimodal_table,file =sprintf("output/%s.txt",GSDname))
  
} else {
  print("Results in unimodal_results and unimodal_plots")
  
  # Built output plot
  unimodal_results <- list(unimodalresults = fitdist_res$distparams)
  unimodal_table <- bind_rows(unimodal_results)
  unimodal_plots <- list(unimodal_simulated = fitdist_res$simuplot, unimodal_plot = fitdist_res$fitdistplot)
  
  unimodtab <- tableGrob(bimodal_results$unimodalresults)
  
  finalplot <- unimodal_plots$unimodal_simulated + unimodal_plots$unimodal_plot + 
    plot_layout(ncol = 2, nrow = 2, widths = c(1,1)) +
    plot_annotation(title = paste("Distribution plots for unimodal", GSDname),
                    theme = theme(plot.title = element_text(hjust =0.5,
                                                            size = 12,
                                                            margin = margin(0,0,0.5,0,"cm")))) + 
    inset_element(unimodtab,-1.5,-1,0,0)
  
  ggsave(filename = sprintf("plots/%s.png",GSDname), finalplot, width = 15, height = 9)
  write.table(unimodal_table,file =sprintf("output/%s.txt",GSDname))
}




######################################################
# OUTPUTS ARE SAVED INTO "outputs/" and "plots/"