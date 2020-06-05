## Bioenergetics model for spring Chinook (R code was created by Rich Zabel and
## was based on Fish Bioenegetics 3.0 Wisconsin model. R code was 
## modiefied by Paul Chittaro, Nov 5, 2013).

# How many objects are available in current directory?
objects()
# Remove objects
rm (list = ls(all = TRUE))
# get working directory
getwd()
#set working directory to...
setwd("C:/Users/lkbob/Documents/NOAA Hollings/Internship/")

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

elk <- c("ELK_CH1_04.csv", "ELK_CH4_04.csv", "ELK_CH5_04.csv")

sink("elk.2004.out")
for(i in 1:3){
        x <- read.table(paste("C:/Users/lkbob/Documents/NOAA Hollings/Internship/",elk[i], sep=""),
		row.names = NULL, header = T)
        for (j in 1:(length(x$mass)-1)){

T <- x$T[j]
Gg <- (x$mass[j+1] - x$mass[j])/x$mass[j]
M <- x$mass[j]
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

cat("elk", " ", elk[i], " ", x$date[j], " ", M, " ", Gg.est, " ", T, " ", p, " ", ED, " ", C, " ", Cj, "\n") 
}
}

