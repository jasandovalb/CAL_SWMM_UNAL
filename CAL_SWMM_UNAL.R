
# Created on Monday Nov 6 11:11:00 2023

# @author: John Alexander Sandoval Barrera
# email: jasandovalb@unal.edu.co

# This script allows to calibrate a SWMM model applying the GLUE, SCE and DDS optimization Methods and compare its performance
# Beven, K. (2012). Rainfall-Runoff Modelling: The Primer: Second Edition. In Rainfall-Runoff Modelling: The Primer: Second Edition. 
# https://doi.org/10.1002/9781119951001
# Duan, Q. Y., Gupta, V. K., & Sorooshian, S. (1993). Shuffled complex evolution approach for effective and efficient global minimization. Journal of Optimization Theory and Applications, 76(3). 
# https://doi.org/10.1007/BF00939380
# Tolson, B. A., & Shoemaker, C. A. (2007). Dynamically dimensioned search algorithm for computationally efficient watershed model calibration. Water Resources Research, 43(1). 
# https://doi.org/10.1029/2005WR004723

# Install the required libraries
# install.packages("devtools")
# library(devtools)
# install_github("dleutnant/swmmr")
# install_github("dkneis/mcu")
# install.packages("xts")
# install.packages("tcltk")
# install.packages("Hmisc")
# install.packages("ggplot2")
# install.packages("SoilHyP")
# install.packages("DEoptim")
# install.packages("pso")
# install.packages("optimization")
# install.packages("GenSA")

# Import required libraries
library(swmmr)
library(xts)
library(tcltk)
library(Hmisc)
library(ggplot2)
library(SoilHyP)
library(mcu)
library(DEoptim)
library(pso)
library(optimization)
library(GenSA)
library(parallelly)

# Set the working directory
setwd("G:/My Drive/Tesis/4_EjecutablesSWMMBase")

# Define a function to calculate the Kling-Gupta efficiency (KGE)
# An additional function to obtain the Pearson correlation coefficient is provided
# Input x is a two column xts object, col1: observed data, col2: simulated data
r <- function(x) {
  sum((x[,1] - mean(x[,1], na.rm = T)) * (x[,2] - mean(x[,2], na.rm = T)), na.rm = T) / 
    sqrt(sum((x[,1] - mean(x[,1], na.rm = T)) ^ 2, na.rm = T) * sum((x[,2] - mean(x[,2], na.rm = T)) ^ 2, na.rm = T))
}

KGE <- function(x) {
  1 - sqrt((mean(x[,2], na.rm = T) / mean(x[,1], na.rm = T) - 1) ^ 2 + (sd(x[,2], na.rm = T) / sd(x[,1], na.rm = T) - 1) ^ 2 + (r(x) - 1) ^ 2)
}

# Define a function to run the SWMM model, get the simulated data time series and evaluate the objective function
#     x: Numeric vector with the values for the calibrated parameters
#     inp: list created with the read_inp function of the swmmr library
#     inp_warmup: list created with the read_inp function of the swmmr library for the warmup period
#     obs: xts object containing the observed data
#     CN_Cal: A number which can be 0 for independent calibration of CN values or 1 for utilization of a scaling factor
#     n_subcatch: number of subcatchments
#     OptimMethod: A number which can be 0 for GLUE optimization and 1 for the rest of methods
#     InpPath: Path where inp files are located
Cal_SWMM <- function(x, inp, inp_warmup, obs, CN_Cal, n_subcatch, OptimMethod, InpPath, TraceSim = FALSE) {
  
  if (CN_Cal == 0) {
    
    # Update parameters in the main project file
    inp$conduits$Roughness <- x[1]
    inp$infiltration$CurveNum <- x[2:(n_subcatch + 1)]
    inp$subareas$`N-Imperv` <- x[n_subcatch + 2]
    inp$subareas$`N-Perv` <- x[n_subcatch + 3]
    inp$subareas$`S-Imperv` <- x[n_subcatch + 4] 
    inp$subareas$`S-Perv` <- x[n_subcatch + 5]
    
    # Update parameters in the warmup project file
    inp_warmup$conduits$Roughness <- x[1]
    inp_warmup$infiltration$CurveNum <- x[2:(n_subcatch + 1)]
    inp_warmup$subareas$`N-Imperv` <- x[n_subcatch + 2]
    inp_warmup$subareas$`N-Perv` <- x[n_subcatch + 3]
    inp_warmup$subareas$`S-Imperv` <- x[n_subcatch + 4] 
    inp_warmup$subareas$`S-Perv` <- x[n_subcatch + 5]
    
  } else if (CN_Cal == 1) {
    
    # Update parameters in the main project file
    inp$conduits$Roughness <- x[1]
    inp$infiltration$CurveNum <- inp$infiltration$CurveNum * x[2] / 100
    inp$infiltration$CurveNum[which(inp$infiltration$CurveNum < 60)] <- 60
    inp$infiltration$CurveNum[which(inp$infiltration$CurveNum > 100)] <- 100
    inp$subareas$`N-Imperv` <- x[3]
    inp$subareas$`N-Perv` <- x[4]
    inp$subareas$`S-Imperv` <- x[5] 
    inp$subareas$`S-Perv` <- x[6]
    
    # Update parameters in the warmup project file
    inp_warmup$conduits$Roughness <- x[1]
    inp_warmup$infiltration$CurveNum <- inp_warmup$infiltration$CurveNum * x[2] / 100
    inp_warmup$infiltration$CurveNum[which(inp_warmup$infiltration$CurveNum < 60)] <- 60
    inp_warmup$infiltration$CurveNum[which(inp_warmup$infiltration$CurveNum > 100)] <- 100
    inp_warmup$subareas$`N-Imperv` <- x[3]
    inp_warmup$subareas$`N-Perv` <- x[4]
    inp_warmup$subareas$`S-Imperv` <- x[5] 
    inp_warmup$subareas$`S-Perv` <- x[6]
    
  }
  
  # write project files (main and warmup) with new parameters in a temporary files
  tmp_inp <- tempfile()
  write_inp(inp, tmp_inp)
  tmp_inp_warmup <- tempfile()
  write_inp(inp_warmup, tmp_inp_warmup)
  
  # run warmup and prepare hotstart file
  run_swmm(tmp_inp_warmup, stdout = NULL)
  file.rename(paste0(InpPath, list.files(InpPath, "Hotstart")), paste0(InpPath ,gsub("_TMP", "", list.files(InpPath, "Hotstart"))))
  
  # run swmm with the new parameter values. stdout = NULL discards output
  swmm_files <- suppressMessages(run_swmm(tmp_inp, stdout = NULL))
  
  # remove files when function exits to avoid heavy disk usage
  on.exit(file.remove(unlist(swmm_files)))
  
  # read and save the simulated series
  sim <- read_out(
    file = swmm_files$out, # path to out file
    iType = 2, # type: link
    object_name = "L8", # name of link
    vIndex = 0 # parameter: flow rate
  )[["L8"]]$flow_rate # directly access to xts object
  
  # Update vector with F.O values for each call
  if (TraceSim) {
    SIM.FO.Values <<- c(SIM.FO.Values, 1 - KGE(merge(obs, sim)))
    print(length(SIM.FO.Values))
  }
  
  # calculate goodness-of-fit and return its value
  if (OptimMethod == 0) {
    return(list(Qsim = sim, FO = 1- KGE(merge(obs, sim))))
  } else if (OptimMethod == 1) {
    return(setNames(1 - KGE(merge(obs, sim)), "1-KGE"))
  }

}

#---------------------------------------------------------------------------------------------------------------------------------
#Initiation of simulations

# Path to the Inp files
InpPath <- "./Nivel_3/08Oct2007/"
# Number of subcatchments for the level of discretization simulated (L1:3, L2:16, L3:84)
n_subcatch <- 84

# Set a seed to get reproducible results
# Set.seed(84)

# Set path to the inp files
inp_file <- paste0(InpPath,"UNAL_SWMM.inp")
inp_file_warmup <- paste0(InpPath,"UNAL_SWMM_WarmUp.inp")

# Read the model structure from the inp files
inp <- read_inp(inp_file)
inp_warmup <- read_inp(inp_file_warmup)

# Explore the inp file
#summary(inp)

# Read the observed data into a xts object
obs_08Oct <- as.data.frame(read.csv("./TimeSeries_Qout_08102007.csv", sep = "", header = F))
obs_08Oct <- cbind(c(rep("2007-10-08", dim(obs_08Oct)[1] - 1), "2007-10-09"), obs_08Oct)
obs_08Oct <- cbind(paste(obs_08Oct[,1], obs_08Oct[,2]), obs_08Oct[,3])
obs_08Oct <- xts(x = as.numeric(obs_08Oct[,2]), order.by = as.POSIXct(obs_08Oct[,1], format = "%Y-%m-%d %H:%M", origin="1970-01-01", tz="GMT"))

obs_13Oct <- as.data.frame(read.csv("./TimeSeries_Qout_13102007.csv", sep = "", header = F))
obs_13Oct <- cbind(c(rep("2007-10-13", dim(obs_13Oct)[1] - 1), "2007-10-14"), obs_13Oct)
obs_13Oct <- cbind(paste(obs_13Oct[,1], obs_13Oct[,2]), obs_13Oct[,3])
obs_13Oct <- xts(x = as.numeric(obs_13Oct[,2]), order.by = as.POSIXct(obs_13Oct[,1], format = "%Y-%m-%d %H:%M", origin="1970-01-01", tz="GMT"))

# Define the lower-limit and upper-limit for the range of variation of each calibrable parameter
#     n: Manning coefficient for the conduits of the sewer system [dimensionless]
#     CN: Curve number for soil conservation service (SCS) loss method [dimensionless]
#     CN_Factor: Define a percentage of change above and below the initial values
#                 For instance, 75 and 125 to consider values 25% above and below the initial value of each subcatchment
#     n_imperv: Manning coefficient for the overland flow through the impervious surfaces [dimensionless]
#     n_perv: Manning coefficient for the overland flow through the pervious surfaces [dimensionless]
#     S_imperv: Depth of depression storage on impervious area [mm]
#     S_perv: Depth of depression storage on pervious area [mm]
# For CN calibration with a scale factor
par_names <- c("n", "CN_Factor", "n_imperv", "n_perv", "s_imperv", "s_perv")
lower_par <- c(0.008, 80, 0.01, 0.01, 1, 1)
upper_par <- c(0.017, 120, 0.03, 0.15, 15, 15)
ini_par <- c(0.015, 100, 0.015, 0.05, 2, 2)
# For individual CN calibation
par_names <- c("n", paste0(rep("CN_", n_subcatch), (1:n_subcatch)), "n_imperv", "n_perv", "s_imperv", "s_perv")
lower_par <- c(0.008, rep(40, n_subcatch), 0.01, 0.01, 1, 1)
upper_par <- c(0.017, rep(100, n_subcatch), 0.03, 0.15, 15, 15)
ini_par <- c(0.015, inp$infiltration$CurveNum, 0.015, 0.05, 2, 2)
# ini_par <- c(runif(1,lower_par[1],upper_par[1]), 
#              runif(n_subcatch,40,100), 
#              runif(1,lower_par[n_subcatch + 2],upper_par[n_subcatch + 2]), 
#              runif(1,lower_par[n_subcatch + 3],upper_par[n_subcatch + 3]),
#              runif(1,lower_par[n_subcatch + 4],upper_par[n_subcatch + 4]),
#              runif(1,lower_par[n_subcatch + 5],upper_par[n_subcatch + 5]))

#---------------------------------------------------------------------------------------------------------------------------------
# Generalized Likelihood Uncertainty Estimation (GLUE)

# Number of Montecarlo simulations
n <- 1000

# Initialize matrix of parameters for MonteCarlo simulations
parameters_set <- matrix(ncol = length(par_names), nrow = n)
colnames(parameters_set) <- par_names

# Define random values for the parameters following the uniform distribution
for(i in 1:length(par_names)){
  parameters_set[,i] <- runif(n, min=lower_par[i], max=upper_par[i])
}

# Create a bar to track calibration progress
pb <- tkProgressBar(title = "Calibrating SWMM", min = 0, max = n, width = 300)

# Create empty lists to save calibration results
Qsim <- list()
FO <- list()

# Perform MonteCalo simulations
system.time(
for(j in 1:n){
  
  SIM <- Cal_SWMM(x = parameters_set[j,], 
                  inp = inp, 
                  inp_warmup = inp_warmup, 
                  obs = obs_08Oct,
                  CN_Cal = 1,
                  n_subcatch = n_subcatch,
                  OptimMethod = 0,
                  InpPath = InpPath)

  Qsim[[j]] <- SIM$Qsim
  FO[[j]] <- SIM$FO

  #Modify progress bar
  setTkProgressBar(pb, j, label=paste(round(j/n*100, 2),"% done"))
  
}
)

#Close progress bar
close(pb)

#Prepare results to post-processing
Qsim <- as.data.frame(do.call(cbind, Qsim))
colnames(Qsim) <- NULL
FO <- do.call(rbind, FO)

#Define behavioral simulations
threshold <- 0.3
parameters_beh <- parameters_set[FO > threshold,]
NSE_beh <- FO[FO > threshold]
Qsim_cal_beh <- Qsim[,FO > threshold]
weights <- NSE_beh - threshold
weights <- weights / sum(weights)
#Calculate confidence interval for 90%
limits <- apply(Qsim_cal_beh, 1, "wtd.quantile", weights = weights,
                probs = c(0.05,0.95), normwt=T)
limits <- t(limits)

#---------------------------------------------------------------------------------------------------------------------------------
# Differential Evolution (DE)

SIM.FO.Values <<- vector(mode = "numeric", length = 0)
CN_Cal <- 1

if (CN_Cal == 1) {
  
  initialpop <- vector(mode = "numeric", length = 0)
  
  for (i in 1:(length(ini_par)*10 - 1)) {
    
    initialpop <- c(initialpop,
                    runif(1,lower_par[1],upper_par[1]),
                    runif(1,lower_par[2],upper_par[2]),
                    runif(1,lower_par[3],upper_par[3]),
                    runif(1,lower_par[4],upper_par[4]),
                    runif(1,lower_par[5],upper_par[5]),
                    runif(1,lower_par[6],upper_par[6]))
    
  }
  
  initialpop <- matrix(data = c(ini_par, initialpop),
                       nrow = length(ini_par)*10, 
                       ncol = length(ini_par),
                       byrow = TRUE)
  
} else if (CN_Cal == 0){
  
  initialpop <- vector(mode = "numeric", length = 0)
  
  for (i in 1:(length(ini_par)*10 - 1)) {
    
    initialpop <- c(initialpop,
                    runif(1,lower_par[1],upper_par[1]),
                    runif(n_subcatch,40,100),
                    runif(1,lower_par[n_subcatch + 2],upper_par[n_subcatch + 2]),
                    runif(1,lower_par[n_subcatch + 3],upper_par[n_subcatch + 3]),
                    runif(1,lower_par[n_subcatch + 4],upper_par[n_subcatch + 4]),
                    runif(1,lower_par[n_subcatch + 5],upper_par[n_subcatch + 5]))
    
  }
  
  initialpop <- matrix(data = c(ini_par, initialpop),
                       nrow = length(ini_par)*10, 
                       ncol = length(ini_par),
                       byrow = TRUE)
  
}

system.time(
SIM <- DEoptim(Cal_SWMM, 
               lower = lower_par, 
               upper = upper_par, 
               control = DEoptim.control(itermax = 16, 
                                         storepopfrom = 1,
                                         initialpop = initialpop), 
               inp = inp, 
               inp_warmup = inp_warmup, 
               obs = obs_08Oct,
               CN_Cal = CN_Cal,
               n_subcatch = n_subcatch,
               OptimMethod = 1,
               InpPath = InpPath,
               TraceSim = TRUE)
)

#---------------------------------------------------------------------------------------------------------------------------------
# Shuffled complex evolution (SCE-UA)

SIM.FO.Values <<- vector(mode = "numeric", length = 0)

system.time(
SIM <- SCEoptim(Cal_SWMM, 
                par = ini_par, 
                lower = lower_par, 
                upper = upper_par, 
                control = list(maxeval = 1000, trace = 1, returnpop = TRUE, ncomplex = 2), #CAMBIAR n!!
                inp = inp, 
                inp_warmup = inp_warmup, 
                obs = obs_08Oct,
                CN_Cal = 0,
                n_subcatch = n_subcatch,
                OptimMethod = 1,
                InpPath = InpPath,
                TraceSim = TRUE)
)

#---------------------------------------------------------------------------------------------------------------------------------
# Dynamically dimensioned search (DDS)

system.time(
SIM <- dds(Cal_SWMM, 
           p = data.frame(name = par_names, initial = ini_par, min = lower_par, max = upper_par), 
           m = 1000, 
           r = 0.2, 
           plot = TRUE,
           inp = inp, 
           inp_warmup = inp_warmup, 
           obs = obs_08Oct,
           CN_Cal = 1,
           n_subcatch = n_subcatch,
           OptimMethod = 1,
           InpPath = InpPath)
)

#---------------------------------------------------------------------------------------------------------------------------------
# Particle Swarm Optimization (PSO)

SIM <- psoptim(par = ini_par, 
               fn = Cal_SWMM, 
               inp = inp, 
               inp_warmup = inp_warmup, 
               obs = obs_08Oct,
               CN_Cal = 1,
               n_subcatch = n_subcatch,
               OptimMethod = 1,
               InpPath = InpPath, 
               lower = lower_par, 
               upper = upper_par, 
               control = list(trace = 1, maxit = 2, REPORT = 1, trace.stars = TRUE))

#---------------------------------------------------------------------------------------------------------------------------------
# Simulated Annealing (SA)

SIM.FO.Values <<- vector(mode = "numeric", length = 0)

system.time(
SIM <- GenSA(par = ini_par, 
             fn = Cal_SWMM, 
             lower = lower_par, 
             upper = upper_par, 
             control = list(verbose = TRUE, max.call = 1000, smooth = TRUE),
             inp = inp, 
             inp_warmup = inp_warmup, 
             obs = obs_08Oct,
             CN_Cal = 1,
             n_subcatch = n_subcatch,
             OptimMethod = 1,
             InpPath = InpPath,
             TraceSim = TRUE)
)

#---------------------------------------------------------------------------------------------------------------------------------
# Graphical outputs

# Graph of confidence interval
data.graph <- data.frame(time = index(obs)[-1], Lower = limits[,1], Upper = limits[,2], Obs = as.numeric(obs)[-1])

ConfidenceInterval <- ggplot(data.graph, aes(x=time)) + geom_ribbon(aes(ymin=Lower, ymax=Upper, fill = "Simulated flow"), color = "black", alpha = 0.5) + 
                        geom_point(aes(y=Obs, color = "Observed flow"), size = 1, alpha = 0.5) +
                        theme_light() + labs(y = "Flow (LPS)") +
                        scale_fill_manual(values = c("Simulated flow" = "black")) +
                        scale_color_manual(values = c("Observed flow" = "red")) +
                        theme(aspect.ratio = 0.50,
                              legend.title = element_blank(),
                              legend.spacing.y = unit(0.005, 'cm'),
                              legend.box.background = element_rect(color = "black", fill = NA ,size = 0.5, linetype = "solid"),
                              legend.box.margin = unit(c(0.05,0.05,0.05,0.05), "cm"),
                              legend.position = c(0.85, 0.8),
                              axis.title = element_text(size = 18, family = "Times New Roman"),
                              axis.text = element_text(size = 14, family = "Times New Roman"),
                              legend.text = element_text(size = 14, family = "Times New Roman"),
                              axis.title.x = element_blank(),
                              axis.ticks.length = unit(.20, "cm"),
                              axis.ticks = element_line(color = "black", linetype = "solid"),
                              panel.border = element_rect(color = "black", size = 0.5, linetype = "solid"))

ggsave(filename = "SWMM_UNAL_ConfidenceInterval.png", plot = ConfidenceInterval, path = ".")

# Graph of best simulation
NSE_max <- max(FO)
NSE_max_pars <- parameters_set[which(FO==NSE_max),]
Qsim_best <- Qsim[,which(FO==NSE_max)]

data.graph_2 <- data.frame(time = index(obs)[-1], Sim = Qsim_best, Obs = as.numeric(obs)[-1])

BestSimulation <- ggplot(data.graph_2, aes(x=time)) + geom_point(aes(y=Obs, fill = "Observed flow"), color = "red", size = 1, alpha = 0.5) + 
                    geom_line(aes(y=Sim, color = "Simulated flow"), size = 0.75) +
                    scale_color_manual(values = c("Simulated flow" = "black")) +
                    theme_light() + labs(y = "Flow (LPS)") +
                    theme(aspect.ratio = 0.50,
                          legend.title = element_blank(),
                          legend.spacing.y = unit(0.005, 'cm'),
                          legend.box.background = element_rect(color = "black", fill = NA ,size = 0.5, linetype = "solid"),
                          legend.box.margin = unit(c(0.05,0.05,0.05,0.05), "cm"),
                          legend.position = c(0.85, 0.8),
                          axis.title = element_text(size = 18, family = "Times New Roman"),
                          axis.text = element_text(size = 14, family = "Times New Roman"),
                          legend.text = element_text(size = 14, family = "Times New Roman"),
                          axis.title.x = element_blank(),
                          axis.ticks.length = unit(.20, "cm"),
                          axis.ticks = element_line(color = "black", linetype = "solid"),
                          panel.border = element_rect(color = "black", size = 0.5, linetype = "solid"))

ggsave(filename = "SWMM_UNAL_BestSimulation.png", plot = BestSimulation, path = ".")

# Dotty plots
data.graph_3 <- as.data.frame(cbind(parameters_set, 1 - FO))
colnames(data.graph_3)[length(par_names) + 1] <- "Obj_Funct" 

for (j in 1:length(par_names)) {
  
  DottyPlot <- ggplot(data.graph_3, aes(x=data.graph_3[,j])) + geom_point(aes(y=Obj_Funct), color = "black", size = 1, alpha = 0.5) + 
                  geom_point(aes(x=data.graph_3[which.max(FO),j], y=1-max(FO)), colour="red", size = 3) +
                  scale_color_manual(values = c("Simulated flow" = "black")) +
                  theme_light() + labs(y = "1 - NSE", x = par_names[j]) +
                  theme(aspect.ratio = 1.00,
                        legend.title = element_blank(),
                        legend.spacing.y = unit(0.005, 'cm'),
                        legend.box.background = element_rect(color = "black", fill = NA ,size = 0.5, linetype = "solid"),
                        legend.box.margin = unit(c(0.05,0.05,0.05,0.05), "cm"),
                        legend.position = c(0.85, 0.8),
                        axis.title = element_text(size = 18, family = "Times New Roman"),
                        axis.text = element_text(size = 14, family = "Times New Roman"),
                        legend.text = element_text(size = 14, family = "Times New Roman"),
                        axis.ticks.length = unit(.20, "cm"),
                        axis.ticks = element_line(color = "black", linetype = "solid"),
                        panel.border = element_rect(color = "black", size = 0.5, linetype = "solid"))
  
  ggsave(filename = paste0("DottyPlot-",par_names[j],".png"), plot = DottyPlot, path = ".")
  
}






















