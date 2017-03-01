library(openair)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# Data input and output locations
input  <- "C:/Users/Mike/Dropbox/Paper Cordbaia/Revision/Rev2/"
output <- "C:/Users/Mike/Desktop/Cordbaia_graphs/"

## Input data sources -------------------------------------------------------------------------------------------

input.path <- paste(input, "Lathuilliere_etal_SM_Data_Measurements.csv", sep = "")

# output for graphs and tables
graph1       <- paste(output, "diffusivities.pdf", sep = "")
graph2       <- paste(output, "effluxes.pdf", sep = "")
graph3       <- paste(output, "soilcurve.pdf", sep = "")
graph4       <- paste(output, "figure3.tif", sep = "")
graph5       <- paste(output, "boxplot.pdf", sep = "")
graph6       <- paste(output, "figureS10.tif", sep = "")
graph7       <- paste(output, "figureS7.tif", sep ="")
graph8       <- paste(output, "figureS8.tif", sep ="")
graph9       <- paste(output, "figure4.tif", sep = "")
graph10      <- paste(output, "figureS5.tif", sep ="")
graph11      <- paste(output, "figureS6.tif", sep ="")
graph12      <- paste(output, "F0vsT.pdf", sep ="")
graph13      <- paste(output, "CO2vsT.pdf", sep ="")

output.table <- paste(output, "Station.csv", sep = "")
output.table2 <- paste(output, "validation30min.csv", sep = "")
output.table3 <- paste(output, "validationday.csv", sep = "")

Station <- read.table(input.path, header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

## Apply data formats and corrections ---------------------------------------------------------------------------

Station$timestamp <- as.POSIXct(Station$timestamp, "%m/%d/%Y %H:%M", tz="GMT")      
Station$date <- as.Date(Station$timestamp, "%Y-%m-%d %H:%M:%S", tz="GMT")  

# Adjustment of CO2.1 with 1500 ppm according to Text S1 in the Supporting Information
# to be used only for efflux calculations
Station$CO2.1.corr <- Station$CO2.1+1500

#----------------------------------------------------------------------------------------------------------------
## soil CO2 efflux calculation

Pair <- Station$Pa*100        # convert barometric pressure into Pascals (Pa)
Tair <- Station$Ta + 273.15   # convert air temperature into Kelvin temperature (K)

CO2.air.ppm <- 400            # assumption for atmospheric CO2 concentration
R  <- 8.3145                  # Ideal gas constant (J mol-1 K-1)

# Error of Vaisala sensors +/-[1.5% range + 2% reading] = +/-[3000 + 2% reading] ppm
Station$e.CO2.1 <- 3000 + 0.02*Station$CO2.1    # at 10 cm depth
Station$e.CO2.2 <- 3000 + 0.02*Station$CO2.2    # at 30 cm depth

# Convert CO2 concentration in ppm to umol/m3
# Temperature measurement error of +/- 1 oC from the MPS2 (Decagon Devices Inc., Pullman, WA, USA)
CO2.air <- (CO2.air.ppm*Pair)/(R*Tair)     # atmospheric concentration under the canopy (umol/m3)
Station$CO2.1.vol <- (Station$CO2.1.corr*Pair)/(R*(Station$MPST.1 + 273.15))   # at 10 cm depth
Station$CO2.2.vol <- (Station$CO2.2*Pair)/(R*(Station$MPST.2 + 273.15))   # at 30 cm depth

# error in CO2.air (umol/m3)
num    <- CO2.air.ppm*Pair
e.num  <- num*sqrt( (0.5/Pair)^2 )
den    <- R*Tair
e.den  <- R*0.3
e.CO2.air <- CO2.air*sqrt( (e.num/num)^2 + (e.den/den)^2 ) 

# error in CO2.1.vol and CO2.2.vol
Station$e.CO2.1.vol <- Station$CO2.1.vol*sqrt( (Station$e.CO2.1/Station$CO2.1)^2 + (0.5*100/Pair)^2 + (0.3/Tair)^2 )   # at 10 cm depth
Station$e.CO2.2.vol <- Station$CO2.2.vol*sqrt( (Station$e.CO2.2/Station$CO2.2)^2 + (0.5*100/Pair)^2 + (0.3/Tair)^2 )   # at 30 cm depth

# Determine CO2 diffusion in soil (Ds) using 3 separate models:

#---- Model 1: Millington & Quirk (1961) - see Text S1 in the Supporting Information ----

fi.10 <- 0.565   #soil porosity at 10 cm depth (g/cm3)
fi.30 <- 0.443   #soil porosity at 30 cm depth (g/cm3)

# Calculate air-filled porosity (m3/m3), following equation (7)
Station$alpha.10 <- fi.10 - Station$VWC.1   #at 10 cm depth
Station$alpha.30 <- fi.30 - Station$VWC.2   #at 30 cm depth
Station$e.alpha.10 <- 0.03   #at 10 cm depth, error assumed similar to the GS3 sensor (Decagon Devices Inc., Pullman, WA, USA)
Station$e.alpha.30 <- 0.03   #at 30 cm depth

# Calculation tortuosity (E) from air-filled porosity
Station$E.10 <- signif((Station$alpha.10^(10/3))/((0.440^2)), digits = 3)    # Tortuosity at 10 cm depth
Station$E.30 <- signif((Station$alpha.30^(10/3))/((0.349^2)), digits = 3)   # Tortuosity at 30 cm depth

# Error in tortuosity E
# error in the total porosity^2 is 0.001 (g/cm3)^2
Station$e.E.10 <- signif(Station$E.10*sqrt( (((10/3)*Station$e.alpha.10)/Station$alpha.10)^2 + ((-2*0.001)/0.440)^2 ), digits = 2)
Station$e.E.30 <- signif(Station$E.30*sqrt( (((10/3)*Station$e.alpha.30)/Station$alpha.30)^2 + ((-2*0.001)/0.349)^2 ), digits = 2)

# Determine the diffusion coefficient of CO2 in air (Da, m2/s) following equation (6)
Station$Da <- signif(14.7*10^-6*((Tair/293.15)^1.75)*(Station$Pa/1013), digits = 3) 

# Calculate the error in Da, as e.Da (m2/s)
A.1   <- (Tair/293.15)
e.A.1 <- (0.3/293.15)
A.2   <- A.1^1.75 
e.A.2 <- A.1*sqrt( ((1.75*e.A.1)/A.1)^2 )
A.3   <- 14.7*10^-6*(Station$Pa/1013)
e.A.3 <- 14.7*10^-6*(0.5/1013)
Station$e.Da <- signif(Station$Da*sqrt( (e.A.2/A.2)^2 + (e.A.3/A.3)^2 ), digits = 1)

# Calculate soil CO2 diffusion (D, m2/s) at 10 and 30 cm following the Millington and Quirk (1961) model
# Remove any values of D when soil volumetric water content > 0.378 m3/m3 and replace by the 
# literature value of 1.92 10^-9 m2/s
Station$D.10 <- signif(ifelse(Station$VWC.1 > 0.378, NA, Station$Da*Station$E.10), digits = 3)  # at 10 cm depth 
Station$D.30 <- signif(ifelse(Station$VWC.2 > 0.378, NA, Station$Da*Station$E.30), digits = 3)  # at 30 cm depth 

D.20 <- vector()
for (i in 1:length(Station$timestamp)){
  d20 <- mean(c(Station$D.10[i], Station$D.30[i])) #, na.rm = TRUE)
  D.20 <- append(D.20, d20)    
}  
Station$D.20 <- D.20

# Calculate the error estimate in D, as e.D.10 (at 10 cm) and e.D.30 (at 30 cm)
Station$e.D.10 <- signif(ifelse(Station$VWC.1 > 0.378, NA, Station$D.10*sqrt( (Station$e.E.10/Station$E.10)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)
Station$e.D.30 <- signif(ifelse(Station$VWC.2 > 0.378, NA, Station$D.30*sqrt( (Station$e.E.30/Station$E.30)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)


#---- Model 2: Moldrup et al. (1999) - see Text S1 in the Supporting Information ---- 

# Calculation tortuosity (E) from air-filled porosity alpha - see Text S1 in the Supporting Information
Station$E.10.2 <- signif( ( (Station$alpha.10)^(2+(3/7.75)) )/(0.44^(3/7.75)), digits = 3)    # Tortuosity at 10 cm depth
Station$E.30.2 <- signif( ( (Station$alpha.30)^(2+(3/7.75)) )/(0.349^(3/7.75)), digits = 3)   # Tortuosity at 30 cm depth

# Error in tortuosity E.2
# error in the total porosity^2 is 0.001
# at 10 cm depth
A.4 <- 2+(3/7.75)
A.5 <- 3/7.75
A.6 <- Station$alpha.10^A.4*sqrt( ((A.4*Station$e.alpha.10)/Station$alpha.10)^2 )
A.7 <- 0.440^A.5*sqrt( ((A.5*0.001)/0.440)^2 )
Station$e.E.10.2 <- signif(Station$E.10.2*sqrt( (A.6/(Station$alpha.10^A.4))^2 + (A.7/(0.440^A.5))^2 ), digits = 2)  #at 10 cm depth

# at 30 cm depth
A.8 <- Station$alpha.30^A.4*sqrt( ((A.4*Station$e.alpha.30)/Station$alpha.30)^2 )
A.9 <- 0.349^A.5*sqrt( ((A.5*0.001)/0.349)^2 )
Station$e.E.30.2 <- signif(Station$E.30.2*sqrt( (A.8/(Station$alpha.30^A.4))^2 + (A.9/(0.349^A.5))^2 ), digits = 2)  #at 10 cm depth

# Calculate soil CO2 diffusion (D, m2/s) at 10 and 30 cm following the Moldrup (1999) model
# Remove any values of D when soil volumetric water content > 0.378 m3/m3 and replace by the 
# literature value of 1.92 10^-9 m2/s
Station$D.10.2 <- signif(ifelse(Station$VWC.1 > 0.378, NA, Station$Da*Station$E.10.2), digits = 3)  # at 10 cm depth
Station$D.30.2 <- signif(ifelse(Station$VWC.2 > 0.378, NA, Station$Da*Station$E.30.2), digits = 3)  # at 30 cm depth


D.20.2 <- vector()
for (i in 1:length(Station$timestamp)){
  d20.2 <- mean(c(Station$D.10.2[i], Station$D.30.2[i])) #, na.rm = TRUE)
  D.20.2 <- append(D.20.2, d20.2)    
}  
Station$D.20.2 <- D.20.2

# Calculate the error estimate in D, as e.D.10 (at 10 cm) and e.D.30 (at 30 cm)
Station$e.D.10.2 <- signif(ifelse(Station$VWC.1 > 0.378, NA, Station$D.10.2*sqrt( (Station$e.E.10.2/Station$E.10.2)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)   #at 10 cm depth
Station$e.D.30.2 <- signif(ifelse(Station$VWC.2 > 0.378, NA, Station$D.30.2*sqrt( (Station$e.E.30.2/Station$E.30.2)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)   #at 30 cm depth


#---- Model 3: from best fit of data - see Text S1 in the Supporting Information ----

# Coefficients of best fit at R2=0.63 - see Text S1 in the Supporting Information
a.1 <- 0.794
b.1 <- 1.97

# Calculation tortuosity (E) from air-filled porosity 
Station$E.10.3 <- signif(a.1*Station$alpha.10^b.1, digits = 3)   #at 10 cm depth
Station$E.30.3 <- signif(a.1*Station$alpha.30^b.1, digits = 3)   #at 30 cm depth

# Error in tortuosity E
Station$e.E.10.3 <- signif(a.1*Station$E.10.3*sqrt( ((b.1*Station$e.alpha.10)/Station$alpha.10)^2 ), digits = 3)   #at 10 cm depth
Station$e.E.30.3 <- signif(a.1*Station$E.30.3*sqrt( ((b.1*Station$e.alpha.30)/Station$alpha.30)^2 ), digits = 3)   #at 30 cm depth

# Calculate soil CO2 diffusion (D, m2/s) at 10/30 cm following the Best fit model
Station$D.10.3 <- ifelse(Station$VWC.1 > 0.378, NA, Station$Da*Station$E.10.3)  # at 10 cm depth
Station$D.30.3 <- ifelse(Station$VWC.2 > 0.378, NA, Station$Da*Station$E.30.3)  # at 30 cm depth

D.20.3 <- vector()
for (i in 1:length(Station$timestamp)){
  d20.3 <- mean(c(Station$D.10.3[i], Station$D.30.3[i])) #, na.rm = TRUE)
  D.20.3 <- append(D.20.3, d20.3)    
}  
Station$D.20.3 <- D.20.3

# Calculate the error estimate in D, as e.D.10 (at 10 cm) and e.D.30 (at 30 cm)
Station$e.D.10.3 <- signif(ifelse(Station$VWC.1 > 0.378, NA, Station$D.10.3*sqrt( (Station$e.E.10.3/Station$E.10.3)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)   # at 10 cm depth
Station$e.D.30.3 <- signif(ifelse(Station$VWC.2 > 0.378, NA, Station$D.30.3*sqrt( (Station$e.E.30.3/Station$E.30.3)^2 + (Station$e.Da/Station$Da)^2 )), digits = 3)   # at 30 cm depth


#---- Calculate CO2 efflux from diffusion coefficients using the gradient method (Tang et al., 2005), from equations (8) and (9) ----

# Tortuosity model from Millington and Quirk (1961) at 5 cm and 20 cm for soil efflux at surface (umol(m2s))
Station$F5  <- signif(-Station$D.10*((CO2.air - Station$CO2.1.vol)/0.10), digits = 3)
Station$F20 <- signif(-Station$D.20*((Station$CO2.1.vol - Station$CO2.2.vol)/0.20), digits = 3)

#Error in F5 and F20 assuming 0.01 m measurement error in depth
B.1 <- (CO2.air - Station$CO2.1.vol)/0.10
B.2 <- B.1*sqrt( ((Station$e.CO2.1.vol + e.CO2.air)/(CO2.air - Station$CO2.1.vol))^2 + (0.01/0.10)^2 )  # error in the ratio for F5
B.3 <- Station$CO2.1.vol - Station$CO2.2.vol
B.4 <- Station$e.CO2.1.vol + Station$e.CO2.2.vol                             
B.5 <- (Station$CO2.1.vol - Station$CO2.2.vol)/0.20
B.6 <- B.5*sqrt( (B.4/B.3)^2 + (0.01/0.20)^2 )          # error in the ratio for F20

Station$e.F5  <- signif(Station$F5*sqrt( (Station$e.D.10/Station$D.10)^2 + (B.2/B.1)^2 ), digits = 3)
Station$e.F20 <- signif(Station$F20*sqrt( (Station$e.D.30/Station$D.30)^2 + (B.6/B.5)^2 ), digits = 3)

# Tortuosity model from Moldrup (1999) at 5 and 20 cm for soil efflux at surface (umol/(m2s))
Station$F5.2  <- signif(-Station$D.10.2*((CO2.air - Station$CO2.1.vol)/0.10), digits = 3)
Station$F20.2 <- signif(-Station$D.20.2*((Station$CO2.1.vol - Station$CO2.2.vol)/0.20), digits = 3)

#Error in efflux assuming 0.01 m measurement error in depth
Station$e.F5.2  <- signif(Station$F5.2*sqrt( (Station$e.D.10.2/Station$D.10.2)^2 + (B.2/B.1)^2 ), digits = 3)
Station$e.F20.2 <- signif(Station$F20.2*sqrt( (Station$e.D.30.2/Station$D.30.2)^2 + (B.6/B.5)^2 ), digits = 3)

# Tortuosity model (our best fit) model for efflux at 5 cm and 20 cm for soil efflux at surface (umol/m2s))
Station$F5.3  <- signif(-Station$D.10.3*((CO2.air - Station$CO2.1.vol)/0.10), digits = 3)
Station$F20.3 <- signif(-Station$D.20.3*((Station$CO2.1.vol - Station$CO2.2.vol)/0.20), digits = 3)

#Error in efflux assuming 0.01 m measurement error in depth
Station$e.F5.3  <- signif(Station$F5.3*sqrt( (Station$e.D.10.3/Station$D.10.3)^2 + (B.2/B.1)^2 ), digits = 3)
Station$e.F20.3 <- signif(Station$F20.3*sqrt( (Station$e.D.30.3/Station$D.30.3)^2 + (B.6/B.5)^2 ), digits = 3)

# Estimate of soil efflux (F0, umol/(m2s))
Station$F0   <- signif((0.20*Station$F5 - 0.05*Station$F20)/0.15, digits = 3)      # from Milington and Quirk (1961)
Station$F0.2 <- signif((0.20*Station$F5.2 - 0.05*Station$F20.2)/0.15, digits = 3)  # from Moldrup (1999)
Station$F0.3 <- signif((0.20*Station$F5.3 - 0.05*Station$F20.3)/0.15, digits = 3)  # from our best fit model

#Remove negative effluxes and convert to g CO2-C/(m2s)
Station$F0   <- ifelse(Station$F0 < 0, NA, Station$F0*43.999*0.273*10^-6)
Station$F0.2 <- ifelse(Station$F0.2 < 0, NA, Station$F0.2*43.999*0.273*10^-6)
Station$F0.3 <- ifelse(Station$F0.3 < 0, NA, Station$F0.3*43.999*0.273*10^-6)

# Error in efflux estimate assuming 0.01 m error in depth
# F0 - from Millington and Quirk (1961) model
C.1 <- 0.20*Station$F5    
C.2 <- C.1*sqrt( (0.01/0.20)^2 + (Station$e.F5/Station$F5)^2 )  # error in first term
C.3 <- 0.05*Station$F20
C.4 <- C.3*sqrt( (0.01/0.05)^2 + (Station$e.F20/Station$F20)^2 ) #error in second term
C.5 <- 0.20*Station$F5 - 0.05*Station$F20
C.6 <- C.5*sqrt( (C.2/C.1)^2 + (C.4/C.3)^2 )  #error in difference (numerator)
Station$e.F0 <- signif((abs(Station$F0))*sqrt( (C.6/C.5)^2 + (0.02/0.15)^2 ), digits = 3)

# F0.2 - from Moldrup et al. (1999) model
G.1 <- 0.20*Station$F5.2    
G.2 <- G.1*sqrt( (0.01/0.20)^2 + (Station$e.F5.2/Station$F5.2)^2 )  # error in first term
G.3 <- 0.05*Station$F20.2
G.4 <- G.3*sqrt( (0.01/0.05)^2 + (Station$e.F20.2/Station$F20.2)^2 ) #error in second term
G.5 <- 0.20*Station$F5.2 - 0.05*Station$F20.2
G.6 <- G.5*sqrt( (G.2/G.1)^2 + (G.4/G.3)^2 )  #error in difference (numerator)
Station$e.F0.2 <- signif((abs(Station$F0.2))*sqrt( (G.6/G.5)^2 + (0.02/0.15)^2 ), digits = 3)

# F0.3 - from best fit model
H.1 <- 0.20*Station$F5.3    
H.2 <- H.1*sqrt( (0.01/0.20)^2 + (Station$e.F5.3/Station$F5.3)^2 )  # error in first term
H.3 <- 0.05*Station$F20.3
H.4 <- H.3*sqrt( (0.01/0.05)^2 + (Station$e.F20.3/Station$F20.3)^2 ) #error in second term
H.5 <- 0.20*Station$F5.3 - 0.05*Station$F20.3
H.6 <- H.5*sqrt( (H.2/H.1)^2 + (H.4/H.3)^2 )  #error in difference (numerator)
Station$e.F0.3 <- signif((abs(Station$F0.3))*sqrt( (H.6/H.5)^2 + (0.02/0.15)^2 ), digits = 3)

# Create a graph with all diffussivity estimates based on diffusion model used
attach(Station)
par(mfrow = c(4,1), mar = c(2,4,2,1))
plot(date, CO2.1/10000, type = "l", yaxs = "i", ylab="")
lines(date, CO2.2/10000, type = "l", col = "blue")
mtext(expression(paste("[CO"[2], "] (10,000 ppm)", sep="")), side=2, line=2, cex=0.75)
legend(date[length(date)/2], 28, legend = c("10 cm", "30 cm"), pch = c(0.5, 0.5),
       bty = "n", lty = 1, ncol = 2, col = c("black", "blue")) 
plot(date, D.10, type = "l", main = " Millington and Quirk [1961]", ylim = c(0, 7E-06), yaxs = "i", ylab="")
lines(date, D.30, type = "l", col = "blue")
mtext(expression(paste(paste(italic("D"[s]), sep=""), " (m"^{2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
legend(date[length(date)/2], 9.25*10^-6, legend = c("10 cm", "30 cm"), pch = c(0.5, 0.5),
       bty = "n", lty = 1, ncol = 2, col = c("black", "blue")) 
plot(date, D.10.2, type = "l", main = "Moldrup et al. [1999]", ylim = c(0, 7E-06), xlab = "", yaxs = "i", ylab="")
lines(date, D.30.2, type = "l", col = "blue")
mtext(expression(paste(paste(italic("D"[s]), sep=""), " (m"^{2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
legend(date[length(date)/2], 8*10^-6, legend = c("10 cm", "30 cm"), pch = c(0.5, 0.5),
       bty = "n", lty = 1, ncol = 2, col = c("black", "blue")) 
plot(date, D.10.3, type = "l", main = "Best fit (this study)", ylim = c(0, 7E-06), yaxs = "i", ylab="")
lines(date, D.30.3, type = "l", col = "blue")
mtext(expression(paste(paste(italic("D"[s]), sep=""), " (m"^{2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
legend(date[length(date)/2], 8*10^-6, legend = c("10 cm", "30 cm"), pch = c(0.5, 0.5),
       bty = "n", lty = 1, ncol = 2, col = c("black", "blue")) 
par(mfrow = c(1,1))
detach(Station)

dev.print(pdf, file=graph1, width=8, height=8, pointsize=10)

# Create a graph with all Efflux estimates based on diffusion model used
attach(Station)
par(mfrow = c(4,1), mar = c(2,4,2,1))
plot(date, CO2.1/10000, type = "l", yaxs = "i", ylab="")
lines(date, CO2.2/10000, type = "l", col = "blue")
mtext(expression(paste("[CO"[2], "] (10,000 ppm)", sep="")), side=2, line=2, cex=0.75)
legend(date[length(date)/2], 28, legend = c("10 cm", "30 cm"), pch = c(0.5, 0.5),
       bty = "n", lty = 1, ncol = 2, col = c("black", "blue")) 
plot(date, F0*1000, type = "l", main = " Millington and Quirk [1961]", ylim = c(0, max(F0*1000, na.rm = TRUE)), yaxs = "i", ylab="")
mtext(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
plot(date, F0.2*1000, type = "l", main = "Moldrup et al. [1999]", ylim = c(0, max(F0*1000, na.rm = TRUE)), xlab = "", yaxs = "i", ylab="")
mtext(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
plot(date, F0.3*1000, type = "l", main = "Best fit (this study)", ylim = c(0, max(F0*1000, na.rm = TRUE)), yaxs = "i", ylab="")
mtext(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
par(mfrow = c(1,1))
detach(Station)

dev.print(pdf, file=graph2, width=8, height=8, pointsize=10)

#----------------------------------------------------------------------------------------------------------------
## Water flow in the soil using hydraulic conductivity

#Calculate the water potential gradient between depths
Station$dWP.dz <- (0.30*Station$WP.2*0.1 - 0.10*Station$WP.1*0.1)/0.2   

# Define the hydraulic conductivity as a function of water potential

### Empirical estimate based on HYPROP(R) data 
#(pF = logWP with WP in cm)
upper1 <- -2.930*log10(-10.197*Station$WP.1) + 5.247 
upper2 <- -2.930*log10(-10.197*Station$WP.2) + 5.247

Station$K.10 <-  10^(upper1)   #in cm/d
Station$K.30 <-  10^(upper2)   #in cm/d

# Derive the mean hydraulic conductivity (cm/d) for 20-25 cm depth
K <- vector()
for (i in 1:length(Station$timestamp)){
  k1 <- mean(c(Station$K.10[i], Station$K.30[i]), na.rm = TRUE)
  K <- append(K, k1)    
}  
Station <- cbind(Station, K)
                                                                            
# Vertical water flow between 10 and 30 cm depths in the soil 

# using data from the HYPROP(R) (Decagon Devices Inc., Pullman, WA, USA) K we estimate the soil water flux J (cm/d)
Station$J <-  - K*Station$dWP.dz 

#----------------------------------------------------------------------------------------------------------------
## Data separation and averaging of results

# half-hourly data for 2014 inundation
Station.2014 <- selectByDate(Station, start = "2013-10-18", end = "2014-10-17")     

# half-hourly data for 2015 inundation
Station.2015 <- selectByDate(Station, start = "2014-10-18", end = "2015-10-17")     

# daily average for entire time series 2013-2015
Station.full.daily <- timeAverage(Station, avg.time = "day")                        

# daily average for 2014 inundation
Station.2014.full.daily <- selectByDate(Station.full.daily, start = "2013-10-18", end = "2014-10-17")

# daily average for 2015 inundation
Station.2015.full.daily <- selectByDate(Station.full.daily, start = "2014-10-18", end = "2015-10-17")

# Daily accumulated precipitation 

PPT.2014 <- data.frame(Station.2014.full.daily$date, Station.2014.full.daily$PPT * 48) 
colnames(PPT.2014) <- c("Date", "PPT")

PPT.2015 <- data.frame(Station.2015.full.daily$date, Station.2015.full.daily$PPT * 48)
colnames(PPT.2015) <- c("Date", "PPT")

PPT.all <- data.frame(Station.full.daily$date, Station.full.daily$PPT * 48)
colnames(PPT.all) <- c("Date", "PPT")


#----------------------------------------------------------------------------------------------------------------
# Soil conditions obtained from the HYPROP(R) (Decagon Devices Inc., Pullman, WA, USA) (m3/m3)

PWP.1 <- min(Station$VWC.1, na.rm = TRUE)   #permanent wilting point based on min volumetric water content (m3/m3)
PWP.2 <- min(Station$VWC.2, na.rm = TRUE)   #permanent wilting point based on min volumetric water content (m3/m3)
FC.1  <- max(Station$VWC.1, na.rm = TRUE)   #field capacity based on max volumetric water content (m3/m3) 
FC.2  <- max(Station$VWC.2, na.rm = TRUE)   #field capacity based on max volumetric water content (m3/m3)

# Calculate available water fraction (Awf between 0-1)
Station.2014$Awf.1 <- signif((Station.2014$VWC.1 - PWP.1)/(FC.1 - PWP.1), digits = 1)
Station.2014$Awf.2 <- signif((Station.2014$VWC.2 - PWP.2)/(FC.2 - PWP.2), digits = 1)

Station.2015$Awf.1 <- signif((Station.2015$VWC.1 - PWP.1)/(FC.1 - PWP.1), digits = 1)
Station.2015$Awf.2 <- signif((Station.2015$VWC.2 - PWP.2)/(FC.2 - PWP.2), digits = 1)

Station$Awf.1 <- signif((Station$VWC.1 - PWP.1)/(FC.1 - PWP.1), digits = 1)
Station$Awf.2 <- signif((Station$VWC.2 - PWP.2)/(FC.2 - PWP.2), digits = 1)

#----------------------------------------------------------------------------------------------------------------
## Soil retention curves and available soil water

# Provide soil water retention curve (VWC vs. WP) curve at 0.10 m and 0.30 m
soil10 <- ggplot() +
  geom_point(data=Station.2014.full.daily, aes(x=WP.1, y=VWC.1, color="2014")) +
  geom_point(data=Station.2015.full.daily, aes(x=WP.1, y=VWC.1, color="2015")) +
  xlab(expression(paste(italic(psi), " (kPa) - at 10 cm depth"), sep="")) +
  ylab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ") - at 10 cm depth"), sep="")) +
  theme_bw() + 
  scale_y_continuous(limits=c(0,0.6)) + scale_x_continuous(limits=c(-33,-10)) +
  labs(color="") + theme(legend.position=c(0.1,0.9)) +
  scale_colour_manual(values=c("blue","black")) +
  theme(legend.title=element_blank(), legend.key = element_blank(), legend.key.size = unit(0.5,"cm"), 
        legend.text = element_text(size = 12))

soil30 <- ggplot() +
  geom_point(data=Station.2014.full.daily, aes(x=WP.2, y=VWC.2), colour="blue") +
  geom_point(data=Station.2015.full.daily, aes(x=WP.2, y=VWC.2)) +
  xlab(expression(paste(italic(psi), " (kPa) - at 30 cm depth"), sep="")) +
  ylab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ") - at 30 cm depth"), sep="")) +
  theme_bw() + 
  scale_y_continuous(limits=c(0,0.6)) + scale_x_continuous(limits=c(-33,-10)) 
  
grid.draw(rbind(ggplotGrob(soil10), ggplotGrob(soil30), size = "last"))

dev.print(pdf, file=graph3, width=8, height=8, pointsize=10)

#----------------------------------------------------------------------------------------------------------------
## Anomalous Accumulation applied to precipitation to define seasons following equation (10)

# Determine average daily precipitation
daily.Precip.2014 <- signif(tail(cumsum(PPT.2014$PPT), 1)/365, 3)
daily.Precip.2015 <- signif(tail(cumsum(PPT.2015$PPT), 1)/365, 3)

#Assess the anomalous accumulation of precipitation following equation (10)
PPT.2014$AA <- cumsum(PPT.2014$PPT - daily.Precip.2014)
PPT.2015$AA <- cumsum(PPT.2015$PPT - daily.Precip.2015)

test.df <- data.frame(1:length(PPT.2014$Date), PPT.2014$Date, PPT.2014$AA)
colnames(test.df) <- c("day", "Date", "AA")
minAA.2014 <- subset(test.df, AA == min(AA))
maxAA.2014 <- subset(test.df, AA == max(AA))

test.df2 <- data.frame(1:length(PPT.2015$Date), PPT.2015$Date, PPT.2015$AA)
colnames(test.df2) <- c("day", "Date", "AA")
minAA.2015 <- subset(test.df2, AA == min(AA))
maxAA.2015 <- subset(test.df2, AA == max(AA))

#--------------------------------------------------------------------------------------------------
## Separation of 2014 and 2015 according to wetting, wet and dry periods in the time series

wet.2014.start      <- minAA.2014$Date
wet.2014.end        <- maxAA.2014$Date - 1
wetting.2014.start  <- Station.2014$date[1]
wetting.2014.end    <- wet.2014.start - 1
dry.2014.start      <- maxAA.2014$Date
dry.2014.end        <- tail(Station.2014$date,1)
wet.2015.start      <- minAA.2015$Date
wet.2015.end        <- maxAA.2015$Date - 1
wetting.2015.start  <- Station.2015$date[1]
wetting.2015.end    <- wet.2015.start - 1
dry.2015.start      <- maxAA.2015$Date
dry.2015.end        <- tail(Station.2015$date,1)

# All 2014-2015 data considered
wet.2014.all        <- selectByDate(Station.2014, start = wet.2014.start, end = wet.2014.end)
wetting.2014.all    <- selectByDate(Station.2014, start = wetting.2014.start, end = wetting.2014.end)
dry.2014.all        <- selectByDate(Station.2014, start = dry.2014.start, end = dry.2014.end)
wet.2015.all        <- selectByDate(Station.2015, start = wet.2015.start, end = wet.2015.end)
wetting.2015.all    <- selectByDate(Station.2015, start = wetting.2015.start, end = wetting.2015.end)
dry.2015.all        <- selectByDate(Station.2015, start = dry.2015.start, end = dry.2015.end)

# Daily average considered for each period
wet.2014.full.daily      <- timeAverage(wet.2014.all, avg.time = "day")
wetting.2014.full.daily  <- timeAverage(wetting.2014.all, avg.time="day")
dry.2014.full.daily      <- timeAverage(dry.2014.all, avg.time="day")
wet.2015.full.daily      <- timeAverage(wet.2015.all, avg.time="day")
wetting.2015.full.daily  <- timeAverage(wetting.2015.all, avg.time="day")
dry.2015.full.daily      <- timeAverage(dry.2015.all, avg.time="day")

#-------------------------------------------------------------------------------------------------
## Plot graphs described in the manuscript and Supporting Information 

# Figure 3 in the main document (daily means represented)

# Precipitation 2013-2015 (P, mm/d)
r1 <- ggplot(PPT.all, aes(x=Date, y=PPT)) +
      geom_bar(stat = "identity") +
      xlab("") + theme_bw() +
      ylab(expression(paste(italic(P), " (mm ", "d"^{-1}, ")"), sep="")) +
      scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2015-10-17")))+
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# CO2 concentration in the soil 2013-2015 ([CO2], 10,000 pppm)
r2 <- ggplot()+
      geom_line(data=Station.full.daily, aes(x=date, y=CO2.1/10000, color="0.10 m")) +
      geom_line(data=Station.full.daily, aes(x=date, y=CO2.2/10000, color="0.30 m")) +
      xlab("") + ylab(expression(paste("[CO"[2], "] (10,000 ppm)"), sep="")) +
      scale_colour_manual(values = c("black", "blue")) +
      theme_bw() +  labs("") +  
      scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2015-10-17"))) +
      ggtitle(expression(bold("B"))) + 
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
      legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

# Modeled efflux using Moldrup et al. (1999) 2013-2015 (F0, umol/(m2s))
r3 <- ggplot() +
      geom_line(data=Station.full.daily, aes(x=date, y=F0.2*1000)) +
      xlab("") + ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep="")) +
      theme_bw() +    
      scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2015-10-17")))+    
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# Oxidation Reduction Potential 2013-2015 (Eh, mV)
r4 <- ggplot() +
      geom_line(data=Station.full.daily, aes(x=date, y=ORP.1, colour="0.10 m")) +
      geom_line(data=Station.full.daily, aes(x=date, y=ORP.2, color="0.30 m")) +
      xlab("") + ylab(expression(paste(italic("E"[h]), " (mV)"), sep="")) +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs("") + 
      scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2015-10-17")))+
      ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
      legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1),
      axis.text.x = element_blank())

# Soil volumetric water content 2013-2015 (VWC, m3/m3)
r5 <- ggplot() +
      geom_line(data=Station.full.daily, aes(x=date, y=VWC.1, colour="0.10 m")) +
      geom_line(data=Station.full.daily, aes(x=date, y=VWC.2, colour="0.30 m")) +
      xlab("") + ylab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ")"), sep="")) +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs(color="") + 
      scale_x_date(breaks = date_breaks("3 months"), labels = date_format("%b-%y"), limits=c(as.Date("2013-10-18"), as.Date("2015-10-17")))+
      ggtitle(expression(bold("E"))) + 
      theme(plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
      legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

grid.draw(rbind(ggplotGrob(r1), ggplotGrob(r2), ggplotGrob(r3), ggplotGrob(r4), ggplotGrob(r5), size = "last"))

dev.print(tiff, file=graph4, width=1500, height=1500, res=150) 


# Boxplot showing differences in CO2 concentrations and F0 by period

df3 <- data.frame(c(wetting.2014.all$CO2.1, wetting.2015.all$CO2.1), 
                  c(wetting.2014.all$CO2.2, wetting.2015.all$CO2.2), as.factor("wetting-up"))
colnames(df3) <- c("CO2.1", "CO2.2", "Season")
df4 <- data.frame(c(wet.2014.all$CO2.1, wet.2015.all$CO2.1),
                  c(wet.2014.all$CO2.2, wet.2015.all$CO2.2), as.factor("wet"))
colnames(df4) <- c("CO2.1", "CO2.2", "Season")
df5 <- data.frame(c(dry.2014.all$CO2.1, dry.2015.all$CO2.1),
                  c(dry.2014.all$CO2.2, dry.2015.all$CO2.2), as.factor("dry"))
colnames(df5) <- c("CO2.1", "CO2.2", "Season")

df6 <- data.frame(c(wetting.2014.all$F0.2, wetting.2014.all$F0.2),
                  as.factor("wetting-up"))
colnames(df6) <- c("F0", "Season")
df7 <- data.frame(c(wet.2014.all$F0.2, wet.2014.all$F0.2),
                  as.factor("wet"))
colnames(df7) <- c("F0", "Season")
df8 <- data.frame(c(dry.2014.all$F0.2, dry.2014.all$F0.2),
                  as.factor("dry"))
colnames(df8) <- c("F0", "Season")

df9 <- rbind(df3, df4, df5)
df10 <- rbind(df6, df7, df8)

# Bloxplot showing seasonal differences in CO2 concentrations at 10 cm depth ([CO2], 10,000 ppm)
u1 <- ggplot(df9, aes(x = Season, y = CO2.1/10000))+
      geom_boxplot()+
      xlab("") + ylab(expression(paste("[CO"[2], "] - 0.10 m (10,000 ppm)"), sep=""))+
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# Bloxplot showing seasonal differences in CO2 concentrations at 30 cm depth ([CO2], 10,000 ppm)
u2 <- ggplot(df9, aes(x = Season, y = CO2.2/10000))+
      geom_boxplot()+
      xlab("") + ylab(expression(paste("[CO"[2], "] - 0.30 m (10,000 ppm)"), sep=""))+
      ggtitle(expression(bold("B"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# Bloxplot showing seasonal differences in soil CO2 efflux (F0, umol/(m2s))
u3 <- ggplot(df10, aes(x = Season, y = F0*1000))+
      geom_boxplot()+
      xlab("") + ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep="")) +
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) 

grid.arrange(u1, u2, u3, ncol=1)

dev.print(pdf, file=graph5, width=8, height=8, pointsize=10)


# Figure S10 in the supporting information (daily means represented)

# Removed values of daily mean efflux when volumetric water content > 0.378 m3/m for the graph
# (affects < 10 points)

# plot modeled soil CO2 efflux (F0, mg CO2-C/(m2s)) as a function soil water potential measured at 10 cm depth (psi, kPa) 
t1 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=WP.1, y=F0.2*1000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=WP.1, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=WP.1, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=WP.1, y=F0.2*1000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=WP.1, y=F0.2*1000), fill=I("red")) +
      theme_bw() +  
      xlab("") + ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")", " at 0.10 m depth"), sep="")) +
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank()) 

# plot modeled soil CO2 efflux (F0, mg CO2-C/(m2s)) as a function soil volumetric water content measured at 10 cm depth (fi, m3/m3) 
t2 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=VWC.1, y=F0.2*1000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=VWC.1, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=VWC.1, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=VWC.1, y=F0.2*1000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=VWC.1, y=F0.2*1000), fill=I("red")) +
      theme_bw() +
      xlab("") + ylab("") +
      ggtitle(expression(bold("B"))) + theme(plot.title=element_text(hjust=0)) +
      scale_x_continuous(limits=c(0.1,0.5)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_blank())

# plot modeled soil CO2 efflux (F0, mg CO2-C/(m2s)) as a function soil water potential measured at 30 cm depth (psi, kPa) 
t3 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=WP.2, y=F0.2*1000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=WP.2, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=WP.2, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=WP.2, y=F0.2*1000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=WP.2, y=F0.2*1000), fill=I("red")) +
      theme_bw() +
      xlab(expression(paste(psi, " (kPa)"), sep="")) +
      ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")", " at 0.30 m depth"), sep="")) +
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) 

# plot modeled soil CO2 efflux (F0, mg CO2-C/(m2s)) as a function soil volumetric water potential measured at 30 cm depth (fi, m3/m3) 
t4 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=VWC.2, y=F0.2*1000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=VWC.2, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=VWC.2, y=F0.2*1000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=VWC.2, y=F0.2*1000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=VWC.2, y=F0.2*1000), fill=I("red")) +
      theme_bw() +
      xlab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ")"), sep="")) + ylab("") +
      ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.y = element_blank()) +
      scale_x_continuous(limits=c(0.1,0.5),breaks=c(0.1,0.2,0.3,0.4,0.5))

grid.arrange(t1,t2, t3, t4, ncol=2)

dev.print(tiff, file=graph6, width=1500, height=1500, res=150)


# Figure S7 in the Supporting Information, date range: 18-Oct-2013 to 17-Oct-2014

# Plot Precipitation (P, mm/d)
p1 <- ggplot(Station.2014.full.daily, aes(x=date, y=PPT.2014$PPT))+
      geom_bar(stat = "identity") +
      xlab("") + ylab(expression(paste(italic(P), " (mm ", "d"^{-1}, ")"), sep="")) +
      theme_bw() +
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(limits = c(0, 80)) +
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

## Anomalous accumulation of precipitation
p2 <- ggplot() +
      geom_line(data=PPT.2014, aes(x=Date, y=AA))+
      xlab("") + ylab(expression(paste(italic(AA(t)), " (mm)", sep=""))) +
      theme_bw() +
      scale_x_date(breaks = date_breaks("3 months"))+
      ggtitle(expression(bold("B"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# Available water fraction
p3 <- ggplot() +
      geom_line(data=Station.2014, aes(x=date, y=Awf.1, colour="at 0.10 m")) +
      geom_line(data=Station.2014, aes(x=date, y=Awf.2, colour="at 0.30 m")) +
      xlab("") + ylab("AWF") +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs(color="") + 
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(labels = function(x){as.character(round(x,1))})+
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(), 
      legend.key = element_blank(), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 7), legend.justification=c(0.5,1))

# Oxidation Reduction Potential
p4 <- ggplot() +
      geom_line(data=Station.2014, aes(x=date, y=ORP.1, colour="at 0.10 m")) +
      geom_line(data=Station.2014, aes(x=date, y=ORP.2, colour="at 0.30 m")) +
      xlab("") + ylab(expression(paste(italic("E"[h]), " (mV)"), sep="")) +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs(color="") + 
      ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(limits=c(-400, 600)) +
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(), 
      legend.key = element_blank(), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 7), legend.justification=c(0.5,1))

# Soil water flux
p5 <- ggplot() +
      geom_line(data=Station.2014.full.daily, aes(x=date, y=J)) +
      ggtitle(expression(bold("E"))) + theme(plot.title=element_text(hjust=0)) +
      theme_bw() +
      xlab("") + ylab(expression(paste(italic(J), " (cm", " d"^{-1}, ")"), sep="")) +
      scale_x_date(breaks = date_breaks("3 months"), labels = date_format("%b-%Y")) +
      scale_y_continuous(limits=c(-1.5, 0.5)) +
      theme(plot.title=element_text(hjust=0))

grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), ggplotGrob(p5), size = "last"))

dev.print(tiff, file=graph7, width=1500, height=1500, res=150)


# Figure S8 in the Supporting Information, date range: 18-Oct-2014 to 17-Oct-2015

## Precipitation (P, mm/d)
q1 <- ggplot(Station.2015.full.daily, aes(x=date, y=PPT.2015$PPT))+
      geom_bar(stat = "identity") +
      xlab("") + ylab(expression(paste(italic(P), " (mm ", "d"^{-1}, ")"), sep="")) +
      theme_bw() +
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(limits=c(0, 80)) +
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

## Anomalous accumulation of precipitation
q2 <- ggplot() +
      geom_line(data=PPT.2015, aes(x=Date, y=AA))+
      xlab("") + ylab(expression(paste(italic(AA(t)), " (mm)", sep=""))) +
      theme_bw() +
      scale_x_date(breaks = date_breaks("3 months"))+
      ggtitle(expression(bold("B"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# Available water fraction
q3 <- ggplot() +
      geom_line(data=Station.2015, aes(x=date, y=Awf.1, colour="at 0.10 m")) +
      geom_line(data=Station.2015, aes(x=date, y=Awf.2, colour="at 0.30 m")) +
      xlab("") + ylab("AWF") +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs(color="") + 
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(labels = function(x){as.character(round(x,1))})+
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(), 
      legend.key = element_blank(), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 7), legend.justification=c(0.5,1))

# Oxidation Reduction Potential
q4 <- ggplot() +
      geom_line(data=Station.2015, aes(x=date, y=ORP.1, colour="at 0.10 m")) +
      geom_line(data=Station.2015, aes(x=date, y=ORP.2, colour="at 0.30 m")) +
      xlab("") + ylab(expression(paste(italic("E"[h]), " (mV)"), sep="")) +
      scale_colour_manual(values=c("black","blue")) + 
      theme_bw() +  labs(color="") + 
      ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
      scale_x_date(breaks = date_breaks("3 months"))+
      scale_y_continuous(limits=c(-400, 600)) +
      theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(), 
      legend.key = element_blank(), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 7), legend.justification=c(0.5,1))

# Soil water flux
q5 <- ggplot() +
      geom_line(data=Station.2015.full.daily, aes(x=date, y=J)) +
      ggtitle(expression(bold("E"))) + theme(plot.title=element_text(hjust=0)) +
      theme_bw() +
      xlab("") + ylab(expression(paste(italic(J), " (cm", " d"^{-1}, ")"), sep="")) +
      scale_x_date(breaks = date_breaks("3 months"), labels = date_format("%b-%Y")) +
      scale_y_continuous(limits=c(-1.5, 0.5)) +
      theme(plot.title=element_text(hjust=0))

grid.draw(rbind(ggplotGrob(q1), ggplotGrob(q2), ggplotGrob(q3), ggplotGrob(q4), ggplotGrob(q5), size = "last"))

dev.print(tiff, file=graph8, width=1500, height=1500, res=150)


# Figure 4 in the main document

# plot soil CO2 concentration ([CO2], 10,000 ppm) as a function soil water potential measured at 10 cm depth (psi, kPa) 
s1 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=WP.1, y=CO2.1/10000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=WP.1, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=WP.1, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=WP.1, y=CO2.1/10000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=WP.1, y=CO2.1/10000), fill=I("red")) +
      theme_bw() +
      xlab("") + ylab(expression(paste("[CO"[2], "] (10,000 ppm) at 0.10 m depth"), sep="")) +
      ggtitle(expression(bold("A"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.x = element_blank())

# plot soil CO2 concentration ([CO2], 10,000 ppm) as a function soil volumetric water content measured at 10 cm depth (fi, kPa) 
s2 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=VWC.1, y=CO2.1/10000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=VWC.1, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=VWC.1, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=VWC.1, y=CO2.1/10000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=VWC.1, y=CO2.1/10000), fill=I("red")) +
      theme_bw() +
      xlab("") + ylab("") +
      ggtitle(expression(bold("B"))) + theme(plot.title=element_text(hjust=0)) +
      scale_x_continuous(limits=c(0.1,0.5)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_blank())

# plot soil CO2 concentration ([CO2], 10,000 ppm) as a function soil water potential measured at 30 cm depth (psi, kPa) 
s3 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=WP.2, y=CO2.1/10000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=WP.2, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=WP.2, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=WP.2, y=CO2.1/10000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=WP.2, y=CO2.1/10000), fill=I("red")) +
      theme_bw() +
      xlab(expression(paste(psi, " (kPa)"), sep="")) +
      ylab(expression(paste("[CO"[2], "] (10,000 ppm) at 0.30 m depth"), sep="")) +
      ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) 

# plot soil CO2 concentration ([CO2], 10,000 ppm) as a function soil volumetric water content measured at 30 cm depth (fi, kPa)   
s4 <- ggplot()+
      geom_point(data=Station.full.daily, shape=1, cex=3, aes(x=VWC.2, y=CO2.1/10000)) +
      geom_point(data=wet.2014.full.daily, shape=21, cex=3, aes(x=VWC.2, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=wet.2015.full.daily, shape=21, cex=3, aes(x=VWC.2, y=CO2.1/10000), fill=I("blue")) +
      geom_point(data=dry.2014.full.daily, shape=21, cex=3, aes(x=VWC.2, y=CO2.1/10000), fill=I("red")) +
      geom_point(data=dry.2015.full.daily, shape=21, cex=3, aes(x=VWC.2, y=CO2.1/10000), fill=I("red")) +
      theme_bw()+
      xlab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ")"), sep="")) + ylab("") +
      ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
      theme(axis.text.y = element_blank()) +
      scale_x_continuous(limits=c(0.1,0.5),breaks=c(0.1,0.2,0.3,0.4,0.5))

grid.arrange(s1,s2, s3, s4, ncol=2)

dev.print(tiff, file=graph9, width=1500, height=1500, res=150)


# Figure S5 in the Supporting Information

# CO2 concentration in the soil 2013-2014 ([CO2], 10,000 pppm)
f2 <- ggplot()+
  geom_line(data=Station.2014.full.daily, aes(x=date, y=CO2.1/10000, color="0.10 m")) +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=CO2.2/10000, color="0.30 m")) +
  xlab("") + ylab(expression(paste("[CO"[2], "] (10,000 ppm)"), sep="")) +
  scale_colour_manual(values = c("black", "blue")) +
  theme_bw() +  labs("") +  
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2014-10-17"))) +
  scale_y_continuous(limits = c(0, 23)) +
  ggtitle(expression(bold("B"))) + 
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

# soil CO2 diffussivities 2013-2014 as determined by Moldrup et al. (1999) (Ds, m2/s)
f3 <- ggplot() +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=D.10.2*1000000, colour="0.10 m")) +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=D.30.2*1000000, color="0.20 m")) +
  xlab("") + ylab(expression(paste(italic("D"[s]), " (10"^{-6}, " m"^{2}, " s"^{-1}, ")"), sep="")) +
  scale_colour_manual(values=c("black","blue")) + 
  theme_bw() +  labs("") + 
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2014-10-17")))+
  ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1),
        axis.text.x = element_blank())

# Modeled efflux using Moldrup et al. (1999) 2013-2014 (F0, mg CO2-C/(m2s))
f4 <- ggplot() +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=F0.2*1000)) +
  xlab("") + ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep="")) +
  theme_bw() +    
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2013-10-18"), as.Date("2014-10-17")))+    
  ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) +
  theme(axis.text.x = element_blank())

# Soil volumetric water content 2013-2015 (VWC, m3/m3)
f5 <- ggplot() +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=VWC.1, colour="0.10 m")) +
  geom_line(data=Station.2014.full.daily, aes(x=date, y=VWC.2, colour="0.20 m")) +
  xlab("") + ylab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ")"), sep="")) +
  scale_colour_manual(values=c("black","blue")) + 
  theme_bw() +  labs(color="") + 
  scale_x_date(breaks = date_breaks("3 months"), labels = date_format("%b-%y"), limits=c(as.Date("2013-10-18"), as.Date("2014-10-17")))+
  scale_y_continuous(limits=c(0, 0.55)) +
  ggtitle(expression(bold("E"))) + 
  theme(plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

grid.draw(rbind(ggplotGrob(p1), ggplotGrob(f2), ggplotGrob(f4), ggplotGrob(f3), ggplotGrob(f5), size = "last"))

dev.print(tiff, file=graph10, width=1500, height=1500, res=150)


# Figure S6 in the Supporting Information

# CO2 concentration in the soil 2014-2015 ([CO2], 10,000 pppm)
f6 <- ggplot()+
  geom_line(data=Station.2015.full.daily, aes(x=date, y=CO2.1/10000, color="0.10 m")) +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=CO2.2/10000, color="0.30 m")) +
  xlab("") + ylab(expression(paste("[CO"[2], "] (10,000 ppm)"), sep="")) +
  scale_colour_manual(values = c("black", "blue")) +
  theme_bw() +  labs("") +  
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2014-10-18"), as.Date("2015-10-17"))) +
  scale_y_continuous(limits = c(0, 23)) +
  ggtitle(expression(bold("B"))) + 
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

# soil CO2 diffussivities 2013-2014 as determined by Moldrup et al.(1999) (Ds, m2/s)
f7 <- ggplot() +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=D.10.2*1000000, colour="0.10 m")) +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=D.30.2*1000000, color="0.20 m")) +
  xlab("") + ylab(expression(paste(italic("D"[s]), " (10"^{-6}, " m"^{2}, " s"^{-1}, ")"), sep="")) +
  scale_colour_manual(values=c("black","blue")) + 
  theme_bw() +  labs("") + 
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2014-10-18"), as.Date("2015-10-17")))+
  ggtitle(expression(bold("D"))) + theme(plot.title=element_text(hjust=0)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1),
        axis.text.x = element_blank())

# Modeled efflux using Moldrup et al. (1999) 2013-2014 (F0, mg CO2-C/(m2s))
f8 <- ggplot() +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=F0.2*1000)) +
  xlab("") + ylab(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep="")) +
  theme_bw() + 
  scale_x_date(breaks = date_breaks("3 months"), limits=c(as.Date("2014-10-18"), as.Date("2015-10-17")))+    
  scale_y_continuous(breaks=round(seq(0, 0.6, by = 0.2), digits = 1), limits = c(0, 0.6)) +
  ggtitle(expression(bold("C"))) + theme(plot.title=element_text(hjust=0)) +
  theme(axis.text.x = element_blank())

# Soil volumetric water content 2013-2015 (VWC, m3/m3)
f9 <- ggplot() +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=VWC.1, colour="0.10 m")) +
  geom_line(data=Station.2015.full.daily, aes(x=date, y=VWC.2, colour="0.30 m")) +
  xlab("") + ylab(expression(paste(italic(theta), " (m"^{3}, " m"^{-3}, ")"), sep="")) +
  scale_colour_manual(values=c("black","blue")) + 
  theme_bw() +  labs(color="") + 
  scale_x_date(breaks = date_breaks("3 months"), labels = date_format("%b-%y"), limits=c(as.Date("2014-10-18"), as.Date("2015-10-17")))+
  scale_y_continuous(limits=c(0, 0.55)) +
  ggtitle(expression(bold("E"))) + 
  theme(plot.title=element_text(hjust=0), legend.position = c(0.95,0.75), legend.title=element_blank(),
        legend.key = element_rect(colour="white"), legend.key.size = unit(0.4,"cm"), legend.text = element_text(size = 9), legend.justification=c(0.5,1))

grid.draw(rbind(ggplotGrob(q1), ggplotGrob(f6), ggplotGrob(f8), ggplotGrob(f7), ggplotGrob(f9), size = "last"))

dev.print(tiff, file=graph11, width=1500, height=1500, res=150)

#---------------------------------------------------------------
# Check correlation of soil CO2 efflux F0.2 with soil temperatures

# soil temperature at 0.10 m depth
a1<-round(coef(lm(Station$MPST.1~Station$F0.2))[2],digits=2) #coef(lm(y~x))[2] is the slope of the regression
b1<-round(coef(lm(Station$MPST.1~Station$F0.2))[1],digits=2)   
r2.1<-round(summary(lm(Station$MPST.1~Station$F0.2))$ r.squared,digits=2)
lm_eq1<-paste0("y=",a1,"x",ifelse(b1>0,"+",""),b1)
R2.1<-bquote(R^2 == .(r2.1)) 

results.T1 <- paste(lm_eq1, " R2=", as.character(R2.1[3]), sep = "")

# soil temperature at 0.30 m depth
a2<-round(coef(lm(Station$MPST.2~Station$F0.2))[2],digits=2) #coef(lm(y~x))[2] is the slope of the regression
b2<-round(coef(lm(Station$MPST.2~Station$F0.2))[1],digits=2)   
r2.2<-round(summary(lm(Station$MPST.2~Station$F0.2))$ r.squared,digits=2)
lm_eq2<-paste0("y=",a2,"x",ifelse(b2>0,"+",""),b2)
R2.2<-bquote(R^2 == .(r2.2)) 

results.T2 <- paste(lm_eq2, " R2=", as.character(R2.2[3]), sep = "")

# plot both graphs for 0.10 m and 0.30 m depths
par(mfrow = c(2,1), mar = c(2,3,2,1), oma = c(1,1,1,1))
plot(Station$MPST.1, Station$F0.2*1000, pch = 20, ylab="", xlab="", xlim = c(18, 34), cex.axis = 0.75)
abline(lm(Station$MPST.1~Station$F0.2))
mtext(results.T1, side=3,line=0.15, cex=0.9)
mtext(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
mtext("Soil temperature, 0.10 m (oC)", side=1, line=2, cex=0.75)
plot(Station$MPST.2, Station$F0.2*1000, pch = 20, ylab = "", xlab = "", xlim = c(18, 34), cex.axis = 0.75)
abline(lm(Station$MPST.2~Station$F0.2))
mtext(results.T2, side=3,line=0.15, cex=0.9)
mtext(expression(paste(paste(italic("F"[0]), sep=""), " (mg CO"[2], "-C m"^{-2}, " s"^{-1}, ")" ), sep=""), side=2, line=2, cex=0.75)
mtext("Soil temperature, 0.30 m (oC)", side=1, line=2, cex=0.75)
par(mfrow = c(1,1))

dev.print(pdf, file=graph12, width=8, height=8, pointsize=10)

#---------------------------------------------------------------
# Check correlation of soil CO2 concentration with soil temperatures

# soil temperature at 0.10 m depth
a3<-round(coef(lm(Station$MPST.1~Station$CO2.1))[2],digits=2) #coef(lm(y~x))[2] is the slope of the regression
b3<-round(coef(lm(Station$MPST.1~Station$CO2.1))[1],digits=2)   
r2.3<-round(summary(lm(Station$MPST.1~Station$CO2.1))$ r.squared,digits=2)
lm_eq3<-paste0("y=",a3,"x",ifelse(b3>0,"+",""),b3)
R2.3<-bquote(R^2 == .(r2.3)) 

results.T3 <- paste(lm_eq3, " R2=", as.character(R2.3[3]), sep = "")

# soil temperature at 0.30 m depth
a4<-round(coef(lm(Station$MPST.2~Station$CO2.2))[2],digits=2) #coef(lm(y~x))[2] is the slope of the regression
b4<-round(coef(lm(Station$MPST.2~Station$CO2.2))[1],digits=2)   
r2.4<-round(summary(lm(Station$MPST.2~Station$CO2.2))$ r.squared,digits=2)
lm_eq4<-paste0("y=",a4,"x",ifelse(b4>0,"+",""),b4)
R2.4<-bquote(R^2 == .(r2.2)) 

results.T4 <- paste(lm_eq4, " R2=", as.character(R2.4[3]), sep = "")

# plot both graphs for 0.10 m and 0.30 m depths
par(mfrow = c(2,1), mar = c(2,3,2,1), oma = c(1,1,1,1))
plot(Station$MPST.1, Station$CO2.1, pch = 20, ylab="", xlab="", xlim = c(18, 34), cex.axis = 0.75)
abline(lm(Station$MPST.1~Station$CO2.1))
mtext(results.T3, side=3,line=0.15, cex=0.9)
mtext(expression(paste("[CO"[2], "] (10,000 ppm)"), sep=""), side=2, line=2, cex=0.75)
mtext("Soil temperature, 0.10 m (oC)", side=1, line=2, cex=0.75)
plot(Station$MPST.2, Station$CO2.2, pch = 20, ylab = "", xlab = "", xlim = c(18, 34), cex.axis = 0.75)
abline(lm(Station$MPST.2~Station$CO2.2))
mtext(results.T4, side=3,line=0.15, cex=0.9)
mtext(expression(paste("[CO"[2], "] (10,000 ppm)"), sep=""), side=2, line=2, cex=0.75)
mtext("Soil temperature, 0.30 m (oC)", side=1, line=2, cex=0.75)
par(mfrow = c(1,1))

dev.print(pdf, file=graph13, width=8, height=8, pointsize=10)

#---------------------------------------------------------------
# Create the data frame to compare field measurements for validation

# select dates from field measurements

date.validation <- c("2014-04-08 08:00:00 GMT",
                     "2014-06-11 08:00:00 GMT",
                     "2014-07-10 08:00:00 GMT",
                     "2014-08-06 08:00:00 GMT",
                     "2014-09-05 08:00:00 GMT",
                     "2014-10-07 08:00:00 GMT",
                     "2014-11-05 08:00:00 GMT",
                     "2014-12-04 08:00:00 GMT",
                     "2015-06-12 08:00:00 GMT",
                     "2015-02-06 08:00:00 GMT",
                     "2015-03-10 08:00:00 GMT",
                     "2015-05-20 08:00:00 GMT",
                     "2015-06-02 08:00:00 GMT",
                     "2015-07-16 08:00:00 GMT",
                     "2015-08-24 08:00:00 GMT",
                     "2015-09-25 08:00:00 GMT")

# Field measurements taken with the EGM4 infrared gas analyzer (PP Systems, UK)
measurements <- c(7.73E-2, 5.90E-2, 2.43E-2, 4.46E-2, 9.06E-2, 3.85E-2, 9.27E-2, 1.01E-1, 
                  8.03E-2, 7.66E-2, 6.45E-2, 5.99E-2, 1.02E-1, 2.05E-2, 2.43E-2, 3.34E-2)

validation <- data.frame()
for (i in date.validation[1:16]) {
  ti <- subset(Station, Station$timestamp == i)
  validation <- rbind(validation, ti)
  }

date.validation.day <- as.Date(date.validation, format = "%Y-%m-%d")

validation2 <- data.frame()
for (i in date.validation.day[1:16]) {
  ti <- subset(Station.full.daily, Station.full.daily$date == i)
  validation2 <- rbind(validation2, ti)
}

validation.30min <- validation[, c("timestamp", "F0", "e.F0", "F0.2", "e.F0.2", "F0.3", "e.F0.3")]
validation.30min$measurements <- measurements

validation.day              <- validation2[,c("date", "F0", "e.F0", "F0.2", "e.F0.2", "F0.3", "e.F0.3")]
validation.day$measurements <- measurements

# Write data

write.table(Station, file= output.table, sep = ",", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(validation.30min, file= output.table2, sep = ",", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(validation.day, file= output.table3, sep = ",", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )

### END ###########################################################