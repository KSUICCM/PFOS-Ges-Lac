#===========================================================================================================================
# Figure 3: Sensitivity analiysis
# - Author : Wei-Chun Chou
# - Date   : March, 2020
#===========================================================================================================================

##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for ggplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(gridExtra)   ## R-package for function "grid.arrange"   

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")
source (file = "HMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_R <- mcode ("PreG_RatPBPK.code", PreG_RatPBPK.code)
Gmod_R    <- mcode ("GRatPBPK.code", GRatPBPK.code)
Lmod_R    <- mcode ("LRatPBPK.code", LRatPBPK.code)
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code) 
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Loading datasets
RDat <- read.csv(file = "Data_R.csv")
HDat <- read.csv(file = "Data_H.csv")

## Loading Fit results from rds files
GFit_R <- readRDS(file = "GFit_R.rds")
LFit_R <- readRDS(file = "LFit_R.rds")
GFit_H <- readRDS(file = "GFit_H.rds")
LFit_H <- readRDS(file = "LFit_H.rds")


## All parameters for rats and humans
Rat.theta.G <- log(c(
# Physiological parameters
    BW                             = 0.185,
    QLC                            = 0.183,
    QKC                            = 0.141,
    QMC                            = 0.002,
    QFC                            = 0.07,
    VLC                            = 0.035,
    VKC                            = 0.0084,
    VMC                            = 0.01,
    VFC                            = 0.07,
    VPlasC                         = 0.0466,
    VFilC                          = 8.4E-4,
    VPTCC                          = 1.35E-4,
    GFRC                           = 41.04,
    QLC_Fet                        = 0.061,
    VPlasC_Fet                     = 0.047,
# Chmical-specific parameters
    Vmax_baso_invitro              = 221, # fitting parameters                      
    Km_baso                        = 19.9, # fitting parameters                       
    Vmax_apical_invitro            = 1808,                        
    Km_apical                      = 278,                        
    RAFbaso                        = 4.15,                       
    RAFapi                         = 1.90,                         
    KeffluxC                       = 2.09,                        
    KbileC                         = 0.007, # fitting parameters                       
    KurineC                        = 1.6,                        
    Free                           = 0.019, #fitting parameters                        
    PL                             = 3.17,  #fitting parameters                      
    PK                             = 0.80,                         
    PM                             = 0.16,
    PF                             = 0.13,
    PRest                          = 0.22,                         
    PPla                           = 0.41,
    K0C                            = 1,                           
    Kabsc                          = 2.12,                        
    Kdif                           = 1e-3,                       
    KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 1.27, # fitting parameters
    Ktrans2C                       = 1,
    Ktrans3C                       = 0.23, # fitting parameters
    Ktrans4C                       = 0.001,
    Free_Fet                       = 0.022, 
    PL_Fet                         = 1.30, # fitting parmaeters
    PRest_Fet                      = 0.11  # fitting parameters
))

Rat.theta.L <- log(c(
# Physiological parameters
    BW0                            = 0.247,
    QLC0                           = 0.196,
    QKC0                           = 0.1161,
    QMC0                           = 0.0887,
    QFC                            = 0.07,
    VLC0                           = 0.046,
    VKC0                           = 0.011,
    VMC0                           = 0.049,
    VFC0                           = 0.1245,
    VPlasC                         = 0.0466,
    VFilC                          = 8.4E-4,
    VPTCC                          = 1.35E-4,
    GFRC                           = 41.04,
# Chemical-specific parameters
    Vmax_baso_invitro              = 393.45,
    Km_baso                        = 27.2,  
    Vmax_apical_invitro            = 4141,# fitting parameter
    Km_apical                      = 278,
    RAFbaso                        = 1.90,
    RAFapi                         = 4.15,
    KeffluxC                       = 2.09,
    KbileC                         = 0.0026,
    KurineC                        = 1.60,
    Free                           = 0.0037, # fitting parameter
    PL                             = 2.72,   # fitting parameter
    PK                             = 0.8,
    PM                             = 0.16,
    PF                             = 0.13,
    PRest                          = 0.033,# fitting parameter
    PMilkM                         = 1.9,  
    PMilkP                         = 0.11,
    PAMilkC                        = 0.5,
    KMilk0                         = 0.28, # fitting parameter 
    K0C                            = 1,
    KabsC                          = 2.12,
    KunabsC                        = 7.05e-5,
    Kdif                           = 1e-3,
    Free_pup                       = 0.022,
    PL_pup                         = 2.55,  # fitting parameter
    PK_pup                         = 0.80,
    PRest_pup                      = 0.16,
    KabsC_pup                      = 2.12,
    KbileC_pup                     = 0.0026,
    KurineC_pup                    = 1.6,
    KeffluxC_p                     = 2.09,
    Kdif_pup_0                     = 0.001,
    Vmax_baso_invitro_p            = 393.45,
    Km_baso_p                      = 27.2,
    Vmax_apical_invitro_p          = 1808,
    Km_apical_p                    = 278
))

Human.theta.G <- log(c(
# Physiological parameters
    BW                             = 67.7,
    QCC                            = 16.4,
    QLC                            = 0.25, 
    QKC                            = 0.141,
    QMC                            = 0.027,
    QFC                            = 0.052,
    VLC                            = 0.026,
    VKC                            = 0.004,
    VFC                            = 0.214,
    VMC                            = 0.0062,
    VFileC                         = 8.4E-4,  
    VPTCC                          = 1.35E-4, 
    GFRC                           = 24.19,
# Chemical-specific parameters
    Vmax_baso_invitro              = 479,                      
    Km_baso                        = 20.1,                        
    Vmax_apical_invitro            = 51803,                        
    Km_apical                      = 248,   # fitting parameter                        
    RAFapi                         = 0.001,                       
    RAFbaso                        = 1,                        
    KeffluxC                       = 0.15,                        
    KbileC                         = 1.3e-4,                       
    KurineC                        = 0.17,                        
    Free                           = 0.0027, # fitting parameter                        
    PL                             = 2.03,                        
    PK                             = 1.26,                         
    PM                             = 0.13,
    PF                             = 0.16,
    PRest                          = 0.20,                         
    PPla                           = 0.13, # fitting parameter
    K0C                            = 1,                           
    Kabsc                          = 2.12,                        
    Kdif                           = 0.001,                       
    KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 0.79, # fitting parameter
    Ktrans2C                       = 1.12, # fitting parameter
    Ktrans3C                       = 0.006,
    Ktrans4C                       = 0.001,
    Free_Fet                       = 0.0038, # fitting parameter
    PL_Fet                         = 0.58,  # fitting parameter
    PRest_Fet                      = 2.3   # fitting parameter
))


Human.theta.L <- log(c(
# Physiological parameters
    BW0                            = 68.05,
    QCC                            = 16.4,
    QLC                            = 0.25,
    QKC                            = 0.141,
    QMC                            = 0.027,
    QFC                            = 0.052,
    VLC                            = 0.026,
    VKC                            = 0.004,
    VMC                            = 0.0062,
    VMilk                          = 0.008,
    VFilC                          = 8.4E-4,
    VPTCC                          = 1.35E-4,
    GFRC                           = 24.19,
# Chemical-specific parameters
    Vmax_baso_invitro              = 479,
    Km_baso                        = 20.1,
    Vmax_apical_invitro            = 51803,
    Km_apical                      = 64.4,
    RAFbaso                        = 1,
    RAFapi                         = 0.525, # fitting parameter
    KeffluxC                       = 0.15, 
    KbileC                         = 1.3e-4,
    KurineC                        = 0.096,
    Free                           = 0.037, # fitting parameter
    PL                             = 2.03,
    PK                             = 1.26,
    PM                             = 0.13,
    PF                             = 0.16,
    PRest                          = 0.20, 
    PMilkM                         = 1.9,
    PMilkP                         = 0.11,
    PAMilkC                        = 0.0028, # fitting parameter
    K0C                            = 1,
    KabsC                          = 2.12,
    KunabsC                        = 7.05e-5,
    Kdif                           = 0.001,
    Free_neo                       = 0.014,
    PL_neo                         = 2.03,
    PK_neo                         = 1.26,
    PRest_neo                      = 0.20,
    Vmax_baso_invitro_neo          = 479,
    Km_baso_neo                    = 20.1,
    Vmax_apical_invitro_neo        = 51803,
    Km_apical_neo                  = 64.4,
    KeffluxC_neo                   = 0.150,
    Kdif_neo                       = 0.001,
    KbileC_neo                     = 1.3e-4,
    KurineC_neo                    = 0.001
))


## Defined the function

pred.Rat <- function(Gpars, Lpars, DOSE) {
    
    ## Get out of log domain
    Gpars <- lapply(Gpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    Lpars <- lapply(Lpars, exp)           ## Return a list of exp (parametrs for lactational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.267                  ## Body weight before gestation
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose 
    GDOSE        = DOSE                   ## Repeat oral dose 
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a exposure scenario
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time 
    
    ## Simulation of exposure scenaior A and B (oral single dose to 15 mg/kg)
    Gout <- 
        Gmod_R %>%
        param (Gpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Goutdf = cbind.data.frame(Time   = Gout$time, 
                              CPlas  = Gout$AUC_CPlas, CPlas_pup = Gout$Plasma_Fet,
                              AUC_CPlas = Gout$AUC_CPlas,AUC_CPlas_pup = Gout$AUC_CPlas_Fet)
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.267                  ## Rat body weight on PND0 (GD21); from Thibodeaux et al., 2003
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from PND0 - PND21
    LDOSE        = DOSE                   ## Repeat oral dose of GDOSE mg/kg/day;  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Lex.oral <- ev (ID   = 1,            ## One individual
                    amt  = LDOSEoral,    ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = LTDOSE - 1,   ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 1) ## Simulation time from PND0 to PND74 with dosing of LTDOSE - 1+ 2 days (No dosing)
    
    
    Lout <- Lmod_R %>% 
        param (Lpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    Loutdf = cbind.data.frame(Time      = Lout$time, 
                              CPlas     = Lout$Plasma,CPlas_pup = Lout$Plasma_pup,
                              AUC_CPlas = Lout$AUC_CPlas,AUC_CPlas_pup = Lout$AUC_CPlas_pup)
    
    Goutdf <- Goutdf %>% filter (Time == 21*24) ## The model output on GD21
    Loutdf <- Loutdf %>% filter (Time == 21*24) ## The model output on PND21
    
    return (list("G" = Goutdf, "L" = Loutdf))
}



pred.Human <- function(Gpars, Lpars, DOSE) {
    
    ## Get out of log domain
    Gpars <- lapply(Gpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    Lpars <- lapply(Lpars, exp)           ## Return a list of exp (parametrs for lactational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 67                     ## Body weight before gestation 
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; 
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create exposure scenario
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 0,             ## Dossed strat time (GD0)
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time 
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>%
        param (Gpars) %>%
        update(atol = 1E-3, maxsteps = 500000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Goutdf = cbind.data.frame(Time   = Gout$time, 
                              CPlas  = Gout$AUC_CPlas, CPlas_pup = Gout$CordB,
                              AUC_CPlas = Gout$AUC_CPlas,AUC_CPlas_pup = Gout$AUC_CPlas_Fet)
    
    
    ## Exposure scenario for Lactational exposure;  (to the end of lactation)
    LBW          = 67                     ## Rat body weight on PND0 
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 41*7                   ## Total dosing/Dose times; R
    LDOSE        = DOSE                   ## Repeat oral dose  (mg/kg/day);  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose (mg)
    
    
    # To create a data set of 1 subject receiving DOSE every 24 hours
    Lex.oral <- ev (ID   = 1,             ## One individual
                    amt  = LDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = LTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 1) ## Simulation time
    
    
    Lout <- Lmod_H %>% # Lactational model
        param (Lpars) %>% # Update the parameter list with pars
        update(atol = 1E-3, maxsteps = 500000) %>% # Atol: absolute tolerance parameter; hmax: The maximum step size; maxstep: maximum number of steps the solver will take when advancing from one time to the next         
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    ## Extract the concentration 
    Loutdf = cbind.data.frame(Time      = Lout$time, 
                              CPlas     = Lout$Plasma, CPlas_neo = Lout$CPneo,
                              AUC_CPlas = Lout$AUC_CPlas,AUC_CPlas_pup = Lout$AUC_CPlas_neo)
    
    Goutdf <- Goutdf %>% filter (Time == 39*24*7) ## Gestational model output on GA39
    Loutdf <- Loutdf %>% filter (Time == 25*24*7) ## Lactational model output on 25 weeks after birth
    return (list("G" = Goutdf, "L" = Loutdf))
    
}

## Extract the best parameters and instaed the inital parmaeters
R.par.G <- GFit_R$par
R.par.L <- LFit_R$par
H.par.G <- GFit_H$par
H.par.L <- LFit_H$par

Rat.theta.G[names(R.par.G)]  <- as.numeric(R.par.G) 
Rat.theta.L[names(R.par.L)]  <- as.numeric(R.par.L) 
Human.theta.G[names(H.par.G)]<- as.numeric(H.par.G) 
Human.theta.L[names(H.par.L)]<- as.numeric(H.par.L) 


## Define the sensitivity function
NSC_func <- function (Gpars, Lpars, Pred, DOSE) {
    nG <- length(Gpars)
    nL <- length(Lpars)
    NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 2)
    NSC_LCA     = matrix(NA, nrow = length(Lpars) , ncol = 2)

    for (i in 1:nG) {
        Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
        Rnew.G         <- Pred(Gpars.new, Lpars, DOSE)
        R.G            <- Pred(Gpars, Lpars, DOSE)
        delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
        
        ## Estimated the AUC
        AUC.GCA.new       =  Rnew.G$G %>% select (AUC_CPlas)
        AUC.GCA.ori       =  R.G$G    %>% select (AUC_CPlas)
        AUC.GCA_pup.new   =  Rnew.G$G %>% select (AUC_CPlas_pup)
        AUC.GCA_pup.ori   =  R.G$G    %>% select (AUC_CPlas_pup)
        
        delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
        delta.AUC.GCA_pup =  AUC.GCA_pup.new - AUC.GCA_pup.ori
        
        NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
        NSC_GCA     [i, 2]   <- as.numeric((delta.AUC.GCA_pup /AUC.GCA_pup.ori) * delta.Gpars)
        
    }
    for (j in 1:nL) {
        Lpars.new      <- Lpars %>% replace(j, log(exp(as.numeric(Lpars[j]))*1.01))
        Rnew.L         <- Pred(Gpars, Lpars.new, DOSE)
        R.L            <- Pred(Gpars, Lpars, DOSE)
        delta.Lpars    <- exp(Lpars[j])/(exp(Lpars[j])*0.01)
        
        ## Estimated the AUC
        AUC.LCA.new       =  Rnew.L$L %>% select (AUC_CPlas)
        AUC.LCA.ori       =  R.L$L    %>% select (AUC_CPlas)
        AUC.LCA_pup.new   =  Rnew.L$L %>% select (AUC_CPlas_pup)
        AUC.LCA_pup.ori   =  R.L$L    %>% select (AUC_CPlas_pup)
        
        delta.AUC.LCA    =  AUC.LCA.new - AUC.LCA.ori
        delta.AUC.LCA_pup    =  AUC.LCA_pup.new - AUC.LCA_pup.ori
        
        NSC_LCA     [j, 1]   <- as.numeric((delta.AUC.LCA/AUC.LCA.ori) * delta.Lpars)
        NSC_LCA     [j, 2]   <- as.numeric((delta.AUC.LCA_pup /AUC.LCA_pup.ori) * delta.Lpars)
        
    }
    return (list(NSC_GCA = NSC_GCA, NSC_LCA = NSC_LCA))
}

# Model results
A <- NSC_func (Rat.theta.G, Rat.theta.L, pred.Rat, 0.1)
B <- NSC_func (Human.theta.G, Human.theta.L, pred.Human, 1.9E-7)

rownames (A$NSC_GCA)  = names(Rat.theta.G)
colnames (A$NSC_GCA)  = c("NSC_CPlas", "NSC_CPlas_pup")
rownames (A$NSC_LCA)  = names(Rat.theta.L)
colnames (A$NSC_LCA)  = c("NSC_CPlas", "NSC_CPlas_pup")

rownames (B$NSC_GCA)  = names(Human.theta.G)
colnames (B$NSC_GCA)  = c("NSC_CPlas", "NSC_CPlas_pup")
rownames (B$NSC_LCA)  = names(Human.theta.L)
colnames (B$NSC_LCA)  = c("NSC_CPlas", "NSC_CPlas_pup")

NSC_GCA_R <- data.frame(A$NSC_GCA)
NSC_LCA_R <- data.frame(A$NSC_LCA)
NSC_GCA_H <- data.frame(B$NSC_GCA)
NSC_LCA_H <- data.frame(B$NSC_LCA)


# Manipulate the data for the use of plot
NSC_GCA_R<- NSC_GCA_R %>% gather(key = "matrix", value = "value")  %>% mutate (Pars = rep (names(Rat.theta.G),2)) %>% mutate (Pregnacy = "Gestation_R")
NSC_LCA_R<- NSC_LCA_R %>% gather(key = "matrix", value = "value")  %>% mutate (Pars = rep (names(Rat.theta.L),2)) %>% mutate (Pregnacy = "Lactation_R")
NSC_GCA_H<- NSC_GCA_H %>% gather(key = "matrix", value = "value")  %>% mutate (Pars = rep (names(Human.theta.G),2)) %>% mutate (Pregnacy = "Gestation_H")
NSC_LCA_H<- NSC_LCA_H %>% gather(key = "matrix", value = "value")  %>% mutate (Pars = rep (names(Human.theta.L),2)) %>% mutate (Pregnacy = "Lactation_H")

NSC_total <- rbind (NSC_GCA_H, NSC_LCA_H,NSC_GCA_R, NSC_LCA_R) 

label <- c("Human gestational model", "Human lactational model", "Rat gestational model", "Rat lactational model")
names(label) <- c("Gestation_H", "Lactation_H","Gestation_R", "Lactation_R")

# Find sensitivity parameters
G_rat_senpars <- unique (NSC_GCA_R %>% filter (abs(value) > 0.2) %>% select(Pars))
L_rat_senpars <- unique (NSC_LCA_R %>% filter (abs(value) > 0.2) %>% select(Pars))
G_H_senpars   <- unique (NSC_GCA_H %>% filter (abs(value) > 0.2) %>% select(Pars))
L_H_senpars   <- unique (NSC_LCA_H %>% filter (abs(value) > 0.2) %>% select(Pars))


# Plot 
p1<- ggplot(data = NSC_total%>%filter (abs(value)>0.01), aes(x = Pars, y = value)) +
    geom_bar (aes(fill = matrix), position = "dodge", stat="identity", show.legend = TRUE) +  
    facet_wrap(~ Pregnacy,  scales = "free_y", labeller = labeller(Pregnacy = label)) +
    scale_y_continuous(expand=c(0,0),limits = c(-1, 1)) +
    coord_flip()



## Add the font to font database
windowsFonts("Times" = windowsFont("Times New Roman"))
Plot_theme <-
    theme_bw () +  
    theme (
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        axis.text.x             = element_text (size = 10, colour = "black", face = "bold"),    # tick labels along axes 
        axis.text.y             = element_text (size = 10, colour = "black", face = "bold"),    # tick labels along axes 
        axis.title              = element_blank(),   
        strip.text.x            = element_text(size = 12, face = "bold", colour = "black"),
        legend.position='none') 
   



p1 <- p1 + Plot_theme+ xlab("") + ylab("") 

# Save the figure
ggsave("Fig.3.tiff",scale = 1,
       plot = p1 ,
       path = "C:/Users/weichunc/Desktop",
       width = 20, height = 25, units = "cm", dpi = 320)




