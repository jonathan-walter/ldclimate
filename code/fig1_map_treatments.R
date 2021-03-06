# make a figure showing a map of population locations and temperature accumulation in different treatmetns

rm(list=ls())

library(rgdal)
library(disturbhf)
library(viridis)
library(RColorBrewer)

#Load needed data files

states<-readOGR("~/Box Sync/Grayson Lab - Box/GM - Dispar Moth/Jeffress Spread Forecasting/Data/Raw/statesp020.shp")

#temp.VAB<-read.csv("/Users/jonathanwalter/Box Sync/Grayson Lab - Box/Gypsy Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/VA_Baseline.csv")
#temp.VA4.5<-read.csv("/Users/jonathanwalter/Box Sync/Grayson Lab - Box/Gypsy Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/VA4.5_Merlin_15min_5-1-19_8-7-19.csv")
#temp.VA8.5<-read.csv("/Users/jonathanwalter/Box Sync/Grayson Lab - Box/Gypsy Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/VA8.5_Guinevere_15min_5-1-19_8-7-19.csv")
#temp.WIB<-read.csv("/Users/jonathanwalter/Box Sync/Grayson Lab - Box/Gypsy Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/WI_Baseline.csv")
#temp.WI4.5
#temp.WI8.5<-read.csv("/Users/jonathanwalter/Box Sync/Grayson Lab - Box/Gypsy Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/WI8.5_Arthur_15min_5-1-19_8-7-19.csv")
climDat<-read.csv("~/Box Sync/Grayson Lab - Box/GM - Dispar Moth/NSF Macrosystems 2018-2019/Climate Simulation Summer 2019/Data/all_ECCC_truncate.csv")

popInfo<-read.csv("~/Box Sync/Grayson Lab - Box/GM - Dispar Moth/Jeffress Spread Forecasting/ClimateSuitabilitySynthesis/source_population_info.csv")

popInfo<-popInfo[popInfo$Population %in% c("AL","IR","MN","WI1","WV1","NC1","NC2"),]

#Select states to include in figure

#mich<-states[states$STATE=="Michigan",]
#mich<-mich[!mich$STATESP020 %in% c("1295"),]
#plot(mich, col=rainbow(length(mich)))

states<-states[states$STATE %in% c("Minnesota","Wisconsin","Michigan","Illinois","Indiana",
                                   "Ohio","West Virginia","Virginia","North Carolina","Maryland",
                                   "Pennsylvania","Delaware","New Jersey","New York","Vermont",
                                   "New Hampshire","Rhode Island","Massachusetts","Maine",
                                   "Connecticut","Iowa","Kentucky"),]

drop<-c( "1718" #STATESP020
        ,"1735"
        ,"1311"
        ,"1524"
        ,"1720"
        ,"1716"
        ,"1323"
        ,"1500"
        ,"1295") 

states<-states[!states$STATESP020 %in% drop,]
#TODO: order populations by hotness
popInfo<-popInfo[match(c("NC2","NC1","WV1","WI1","MN","AL","IR"),popInfo$Population),]


popPal<-brewer.pal(8, "RdYlBu")[-4]


# #fix dates
# temp.VAB$DateTime<-as.POSIXct(temp.VAB$DateTime, format="%m/%d/%Y %H:%M")
# temp.VA4.5$DateTime<-as.POSIXct(temp.VA4.5$DateTime, format="%m/%d/%Y %H:%M")
# temp.VA8.5$DateTime<-as.POSIXct(temp.VA8.5$DateTime, format="%m/%d/%Y %H:%M")
# temp.WIB$DateTime<-as.POSIXct(temp.WIB$DateTime, format="%m/%d/%Y %H:%M")
# #temp.WI4.5$DateTime<-as.POSIXct(temp.WI4.5$DateTime, format="%m/%d/%Y %H:%M")
# temp.WI8.5$DateTime<-as.POSIXct(temp.WI8.5$DateTime, format="%m/%d/%Y %H:%M")
# 
# #fill NAs with medians
# temp.VAB$Temp[is.na(temp.VAB$Temp)]<-median(temp.VAB$Temp, na.rm=T)
# temp.VA4.5$Temp[is.na(temp.VA4.5$Temp)]<-median(temp.VA4.5$Temp, na.rm=T)
# temp.VA8.5$Temp[is.na(temp.VA8.5$Temp)]<-median(temp.VA8.5$Temp, na.rm=T)
# temp.WIB$Temp[is.na(temp.WIB$Temp)]<-median(temp.WIB$Temp, na.rm=T)
# #temp.WI4.5$Temp[is.na(temp.WI4.5$Temp)]<-median(temp.WI4.5$Temp, na.rm=T)
# temp.WI8.5$Temp[is.na(temp.WI8.5$Temp)]<-median(temp.WI8.5$Temp, na.rm=T)
# 
# ## want to truncate lines according to last 5th instar date for that treatment
# # VAB = 178
# # WIB = 183
# # VA4.5 = 171
# # WI4.5 = 169
# # VA8.5 = 161
# # WI8.5 = 163
# 
# temp.VAB<-temp.VAB[decimal_doy(temp.VAB$DateTime)<179,]
# temp.VA4.5<-temp.VA4.5[decimal_doy(temp.VA4.5$DateTime)<172,]
# temp.VA8.5<-temp.VA8.5[decimal_doy(temp.VA8.5$DateTime)<162,]
# temp.WIB<-temp.WIB[decimal_doy(temp.WIB$DateTime)<184,]
# #temp.WI4.5<-temp.WI4.5[decimal_doy(temp.WI4.5$DateTime)<170,]
# temp.WI8.5<-temp.WI8.5[decimal_doy(temp.WI8.5$DateTime)<164,]

TrtLocs <- data.frame(lat = c(37.1122, 45.7992),
                      lon = c(-77.2017, -90.9947))


trtPal<-magma(9)[3:8]



#  Treatment FifthTime  simDOY
#1       WIB     54.25  199
#2     WI4.5     42.00  177
#3     WI8.5     35.65  168
#4       VAB     49.00  148
#5     VA4.5     42.00  136
#6     VA8.5     34.00  130


png("~/Box Sync/Grayson Lab - Box/GM - Dispar Moth/NSF Macrosystems 2018-2019/Climate Simulation Summer 2019/Manuscript files/fig_map_treatments.png",
    res=300, units="in", width=6.5, height=3.25)

layout(matrix(1:2,ncol=2,byrow=T), widths=c(0.6,0.4))

par(mar=rep(0,4))

plot(states, border="grey90", col="grey", lwd=0.5)
points(popInfo$Longitude, popInfo$Latitude, pch=21, col="black", bg=popPal)
points(TrtLocs$lon, TrtLocs$lat, pch=8, cex=0.7)
legend("bottomleft", bty="n", legend=c("NC2","NC1","WV","WI","MN","AL","IR"),
       pch=21, pt.bg=popPal, cex=0.9, inset=c(0.03,0))
legend("bottom", bty="n", legend="Climate change simulation location", pch=8, cex=0.7)
text(TrtLocs$lon, TrtLocs$lat, c("VA","WI"), pos=c(3,4))

text(par("usr")[1] + 0.05*abs(diff(par("usr")[1:2])),
     par("usr")[4] - 0.07*abs(diff(par("usr")[3:4])),"a)")

par(mar=c(3.5,3.5,3.5,1.1), mgp=c(2,0.6,0), tcl=-0.3, xpd=F)



plot(NA, NA, xlim=c(90,240), ylim=c(20,36), xlab="Day of year", ylab="Daily max temperature")
abline(h=29, lty=2, col="darkgrey",lwd=1.5)
text(par("usr")[1] + 0.07*abs(diff(par("usr")[1:2])),
     par("usr")[4] - 0.07*abs(diff(par("usr")[3:4])),"b)")

par(xpd=T)
lines(climDat$DOY[climDat$Simulation=="WI.Baseline"], climDat$MaximumAirTemperature[climDat$Simulation=="WI.Baseline"]
      , col=trtPal[1], lwd=2)
points(199,climDat$MaximumAirTemperature[climDat$Simulation=="WI.Baseline" & climDat$DOY==199], pch=5, cex=0.5)
lines(climDat$DOY[climDat$Simulation=="VA.Baseline"], climDat$MaximumAirTemperature[climDat$Simulation=="VA.Baseline"]
      , col=trtPal[4], lwd=2)
points(148,climDat$MaximumAirTemperature[climDat$Simulation=="VA.Baseline" & climDat$DOY==148], pch=5, cex=0.5)
lines(climDat$DOY[climDat$Simulation=="WI.ECCC4.5"], climDat$MaximumAirTemperature[climDat$Simulation=="WI.ECCC4.5"]
      , col=trtPal[2], lwd=2)
points(177,climDat$MaximumAirTemperature[climDat$Simulation=="WI.ECCC4.5" & climDat$DOY==177], pch=5, cex=0.5)
lines(climDat$DOY[climDat$Simulation=="WI.ECCC8.5"], climDat$MaximumAirTemperature[climDat$Simulation=="WI.ECCC8.5"]
      , col=trtPal[3], lwd=2)
points(168,climDat$MaximumAirTemperature[climDat$Simulation=="WI.ECCC8.5" & climDat$DOY==168], pch=5, cex=0.5)
lines(climDat$DOY[climDat$Simulation=="VA.ECCC4.5"], climDat$MaximumAirTemperature[climDat$Simulation=="VA.ECCC4.5"]
      , col=trtPal[5], lwd=2)
points(136,climDat$MaximumAirTemperature[climDat$Simulation=="VA.ECCC4.5" & climDat$DOY==136], pch=5, cex=0.5)
lines(climDat$DOY[climDat$Simulation=="VA.ECCC8.5"], climDat$MaximumAirTemperature[climDat$Simulation=="VA.ECCC8.5"]
      , col=trtPal[6], lwd=2)
points(130,climDat$MaximumAirTemperature[climDat$Simulation=="VA.ECCC8.5" & climDat$DOY==130], pch=5, cex=0.5)
legend("top", ncol=2, inset=-0.35,bty="n",lwd=1,col=trtPal,
        legend=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"), cex=0.75)


dev.off()