## Bioenergetics model for spring Chinook (R code was created by Rich Zabel and
## was based on Fish Bioenegetics 3.0 Wisconsin model. R code was 
## modiefied by Paul Chittaro, Nov 5, 2013).

library(dplyr)
library(zoo)

# How many objects are available in current directory?
objects()
# Remove objects
rm (list = ls(all = TRUE))
# get working directory
getwd()
#set working directory to...
setwd("C:/Users/lkbob/Documents/NOAA Hollings/Internship/data/")

#Consumption parameters
G1 <- 0.5274
G2 <- 1.4724
CQ <- 5
CTM <- 18
CTL <- 24
CK1 <- 0.36
CK4 <- 0.01
CA <- 0.303
CB <- -0.275
#Respiration parameters
RA <- 0.00264
RB <- -0.217
RQ <- 0.06818
RTO <- 0.0234
RTM <- 0
RTL <- 25
RK1 <- 1
RK4 <- 0.13
ACT <- 9.7
BACT <- 0.0405
SDA <- 0.172
OXY <- 3240
TC <- 0.239
#Egestion/Excretion Parameters
FA <- 0.212
FB <- -0.222
FG <- 0.631
UA <- 0.0314
UB <- 0.58
UG <- -0.299
#Predator energy density Parameters
alpha2 <- 7602
beta2 <- 0.5266

write(c("FishID", "row.name", "doy", "Stream","Year","Gg.est", "T", "p", "ED", "C", "Cj"), file = "Bioenergetics output.csv", sep = ",", append = FALSE, ncolumns = 11)

xFull <- read.csv("chinook data for Chittaro bioenergetics 3000 ED.csv") 

# xFull <- filter(xFull, fl >= 40)
xFull <- filter(xFull, fl >= 45)

xFull <- filter(xFull, T != "NA")

#Runs all the way through with all of these steps and temps<10 changed to 10
#xFull$T[xFull$T < 5] <- 5
#xFull$T[xFull$T < 6] <- 6
#xFull$T[xFull$T < 7] <- 7
#xFull$T[xFull$T < 8] <- 8
#xFull$T[xFull$T < 10] <- 10

#Need to remove this fish (only has one line)
#xFull <- filter(xFull, id != "MAR.2006.9.C9" & id != "SFS.2009.8.C8" & id != "SFS.2011.7.C3") #when FL >= 45 mm
#xFull <- filter(xFull, id != "LAK.2005.7.C7") #when FL >= 50 mm
xFull <- xFull %>% 
  add_count(id) %>%
  filter(n > 8)

n <- 1

for(i in unique(xFull$id)){
  x <- subset(xFull, id == i)
  print(i)
  print(paste0(as.character(n),"/",as.character(length(unique(xFull$id)))))
  n <- n + 1
  x$mass.rollmean <- rollmean(x$mass, 7, fill = NA, align = "right")
  x$T.rollmean <- rollmean(x$T, 7, fill = NA, align = "right")
  x <- filter(x, T.rollmean != "NA")
  o <- 1
  for (j in 1:(length(x$mass)-1)){
    print(j)
    print(paste0(as.character(o),"/",as.character(length(unique(x$mass)))))
    o <- o + 1
    
    T <- x$T.rollmean[j]
    Gg <- (x$mass.rollmean[j+1] - x$mass.rollmean[j])/x$mass.rollmean[j]
    M <- x$mass.rollmean[j]
    ED <- x$ED[j]
    
    Gg.est <- 0
    p <- 0
    while (Gg.est < Gg){
      p.old <- p
      Gg.est.old <- Gg.est
      p <- p + 0.0001
      #Consumption
      L1 <- exp(G1*(T-CQ))
      L2 <- exp(G2*(CTL-T))
      KA <- (CK1*L1)/(1+CK1*(L1-1))
      KB <- (CK4*L2)/(1+CK4*(L2-1))
      fc <- KA*KB
      Cmax <- CA*M^CB
      C <- Cmax*p*fc
      Cj <- C*ED
      
      #Respiration
      VEL <- ACT*(M^RK4)*exp(BACT*T)
      ACTIVITY <- exp(RTO*VEL)
      ft <- exp(RQ*T)
      R <- RA*(M^RB)*ft*ACTIVITY*(OXY/TC)
      
      #Egestion
      FP <- FA*(T^FB)*exp(FG*p)
      F <- FP*Cj
      
      #Excretion
      UP <- UA*(T^UB)*exp(UG*p)
      U <- UP*(Cj-F)
      
      #Specific Dynamic Action
      S <- SDA*(Cj-F)
      
      #Growth in joules
      G <- Cj - R - F - U - S
      
      #Growth in grams
      PED <- alpha2 + beta2*M
      Gg.est <- G/PED
      #cat(Gg.est, " ")
    }
    #print(c(p, Gg , Gg.est))
    prop <- (Gg - Gg.est.old)/(Gg.est - Gg.est.old)
    p <- p.old + prop*(p-p.old)
    #print(p)
    
    #Consumption
    L1 <- exp(G1*(T-CQ))
    L2 <- exp(G2*(CTL-T))
    KA <- (CK1*L1)/(1+CK1*(L1-1))
    KB <- (CK4*L2)/(1+CK4*(L2-1))
    fc <- KA*KB
    Cmax <- CA*M^CB
    C <- Cmax*p*fc
    Cj <- C*ED
    
    #Respiration
    VEL <- ACT*(M^RK4)*exp(BACT*T)
    ACTIVITY <- exp(RTO*VEL)
    ft <- exp(RQ*T)
    R <- RA*(M^RB)*ft*ACTIVITY*(OXY/TC)
    
    #Egestion
    FP <- FA*(T^FB)*exp(FG*p)
    F <- FP*Cj
    
    #Excretion
    UP <- UA*(T^UB)*exp(UG*p)
    U <- UP*(Cj-F)
    
    #Specific Dynamic Action
    S <- SDA*(Cj-F)
    
    #Growth in joules
    G <- Cj - R - F - U - S
    
    #Growth in grams
    PED <- alpha2 + beta2*M
    Gg.est <- G/PED
    
    write(c(i,j,x$doy[j], as.character(x$stream[1]),x$year[1], Gg.est, T, p, ED, C, Cj), file = "Bioenergetics output.csv", sep = ",", append = TRUE, ncolumns = 11)
  }
}
