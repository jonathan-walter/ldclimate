## test and evaluate statistical models of gm physiological performance
## in growth chamber climate simulation experiments

## This version adds interactions between sex and treatment

rm(list=ls())

library(here)
library(lsmeans)
library(RColorBrewer)
library(lmerTest)
library(car)
library(viridis)

dat.raw<-read.csv(here("./data/Climate Simulation Exp-Master Data.csv"),
                  stringsAsFactors = FALSE)

neonate.masses<-read.csv(here("./data/neonate_masses.csv"))

str(dat.raw)
dat<-dat.raw

#correct some data entry anomalies
dat$AdultDOY[dat$AdultDOY < 150]<-NA
dat$Sex[dat$Sex=="FP"]<-"F"
dat$Sex[dat$Sex=="MP"]<-"M"

#create some new time-to-event variables
dat$ThirdTime = dat$ThirdDOY - dat$StartDOY
dat$FifthTime = dat$FifthDOY - dat$StartDOY
dat$PupaTime = dat$PupaDOY - dat$StartDOY
dat$AdultTime = dat$AdultDOY - dat$StartDOY

#add population-specific average neonate mass
dat$NeonateMass <- rep(NA, nrow(dat))
for(ii in 1:nrow(dat)){
  dat$NeonateMass[ii] <- neonate.masses$mass_g[neonate.masses$Population==dat$Population[ii]]
  
}

#compute mass gain/time at 3rd and 5th instar
dat$ThirdGain <- (dat$ThirdMass-dat$NeonateMass)/dat$ThirdTime
dat$FifthGain <- (dat$FifthMass-dat$NeonateMass)/dat$FifthTime

#binary adult survival variable
dat$SurvBin<-ifelse(dat$Survival=="A",1,0)

#dat<-dat[dat$Sex!="U",]

dat$Population[dat$Population=="WI1"]<-"WI"
dat$Population[dat$Population=="WV1"]<-"WV"

dat$Treatment<-factor(dat$Treatment, levels=c("WIB","WI4.5","WI8.5","VAB","VA4.5","VA8.5"), ordered=T)
dat$Population<-factor(dat$Population, levels=c("NC2","NC1","WV","WI","MN","AL","IR"), ordered=T)


##analyses on full dataset -------------------------------------------------------------------------------

popPal<-brewer.pal(8, "RdYlBu")[-4]


ThirdGain.lmer<-lmer(ThirdGain ~ Treatment + Population + Treatment:Population + (1|Sex), data=dat,
                 contrasts=list(Treatment ='contr.sum', Population='contr.sum'))
Anova(ThirdGain.lmer, type="III")
summary(ThirdGain.lmer)


ThirdGain.trt.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Treatment))
ThirdGain.pop.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Population))
ThirdGain.trtPop.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Treatment:Population))

trtPal<-magma(9)[3:8]


# png(here("./figures/FigS1_ThirdGain_OLD.png"), units="in", res=300, width=8.5, height=3)
# 
# par(mar=c(6.1,3.6,1.1,1.1), mfrow=c(1,3), oma=c(0,2,0,0), mgp=c(4.8,1,0), cex.lab=1.2)
# boxplot(ThirdGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(0,0.0025),
#         names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
#         col=trtPal)
# text(x=c(1,4,5),y=0.00015,"a")
# text(x=c(2,3,6),y=0.00015,"b")
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"a)")
# 
# boxplot(ThirdGain ~ Population, las=2, data=dat, ylab="", ylim=c(0,0.0025), col=popPal)
# text(x=c(4,5),y=0.00015,"a")
# text(x=c(2,6,7),y=0.00015,"b")
# text(x=c(1,3),y=0.00015,"c")
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"b)")
# 
# plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0,0.0025))
# axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
# lines(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC2"], type="b", col=popPal[1], lwd=2) #NC2
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC1"], type="b", col=popPal[2], lwd=2) #NC1
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WV"], type="b", col=popPal[3], lwd=2) #WV1
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WI"], type="b", col=popPal[4], lwd=2) #WI1
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="MN"], type="b", col=popPal[5], lwd=2) #MN
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="AL"], type="b", col=popPal[6], lwd=2) #AL
# points(1:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="IR"], type="b", col=popPal[7], lwd=2) #IR
# legend("topright", lty=1, col=popPal, legend=c("NC2","NC1","WV","WI","MN","AL","IR"), cex=1, bty="n", ncol=2)
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")
# 
# mtext("Growth to 3rd instar (g/day)", 2, outer=T, at=0.6, cex=0.85, line=0.7)
# 
# dev.off()
# 



png(here("./figures/FigS1_ThirdGain.png"), units="in", res=300, width=6.5, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T))
par(mar=c(5.2,2.6,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(ThirdGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(0,0.0027),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=c(1,4,5),y=0.00015,"a")
text(x=c(2,3,6),y=0.00015,"b")
#text(x=c(4,5,6),y=0.007,"c")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(ThirdGain ~ Population, las=2, data=dat, ylab="", ylim=c(0,0.0027), col=popPal)
text(x=c(4,5),y=0.00015,"a")
text(x=c(2,6,7),y=0.00015,"b")
text(x=c(1,3),y=0.00015,"c")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0.0005,0.002))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
lines(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC2"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC1"][1:3], type="b", col=popPal[2], lwd=2) #NC1
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WV"][1:3], type="b", col=popPal[3], lwd=2) #WV1
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WI"][1:3], type="b", col=popPal[4], lwd=2) #WI1
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="MN"][1:3], type="b", col=popPal[5], lwd=2) #MN
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="AL"][1:3], type="b", col=popPal[6], lwd=2) #AL
points(1:3,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="IR"][1:3], type="b", col=popPal[7], lwd=2) #IR
lines(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC2"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="NC1"][4:6], type="b", col=popPal[2], lwd=2) #NC1
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WV"][4:6], type="b", col=popPal[3], lwd=2) #WV1
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="WI"][4:6], type="b", col=popPal[4], lwd=2) #WI1
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="MN"][4:6], type="b", col=popPal[5], lwd=2) #MN
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="AL"][4:6], type="b", col=popPal[6], lwd=2) #AL
points(4:6,ThirdGain.trtPop.emm$emmean[ThirdGain.trtPop.emm$Population=="IR"][4:6], type="b", col=popPal[7], lwd=2) #IR
legend("top", lty=1, col=popPal, legend=c("NC2","NC1","WV","WI","MN","AL","IR"), cex=1, bty="n", ncol=7)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Growth to 3rd instar (g/day)", 2, outer=T, at=0.6, cex=0.9, line=0.7)

dev.off()



### Fifth Instar ----------------------------------------------------------------------------------

FifthGain.lmer<-lmer(FifthGain ~ Treatment + Population + Treatment:Population + (1|Sex), data=dat, 
                 contrasts=list(Treatment ='contr.sum', Population='contr.sum'))
hist(residuals(FifthGain.lmer))
Anova(FifthGain.lmer, type="III")


FifthGain.trt.emm<-print(emmeans(FifthGain.lmer, specs= ~ Treatment))
FifthGain.pop.emm<-print(emmeans(FifthGain.lmer, specs= ~ Population))
FifthGain.trtPop.emm<-print(emmeans(FifthGain.lmer, specs= ~ Treatment:Population))

# png(here("./figures/Fig2_FifthGain_OLD.png"), units="in", res=300, width=8.5, height=3)
# 
# par(mar=c(6.1,2.6,1.1,1.1), mfrow=c(1,3), oma=c(0,2,0,0), mgp=c(4.8,1,0), cex.lab=1.2)
# boxplot(FifthGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(0,0.02),
#         names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
#         col=trtPal)
# text(x=1,y=0,"a")
# text(x=4,y=0,"b")
# text(x=5,y=0,"c")
# text(x=c(2,3,6),y=0,"d")
# #text(x=c(4,5,6),y=0.007,"c")
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"a)")
# 
# boxplot(FifthGain ~ Population, las=2, data=dat, ylab="", ylim=c(0,0.02), col=popPal)
# text(x=5,y=0,"a")
# text(x=4,y=0,"b")
# text(x=c(1,3,7),y=0,"c")
# text(x=c(2,6),y=0,"bc")
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"b)")
# 
# plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0,0.02))
# axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
# lines(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC2"][1:3], type="b", col=popPal[1], lwd=2) #NC2
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC1"][1:3], type="b", col=popPal[2], lwd=2) #NC1
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WV"][1:3], type="b", col=popPal[3], lwd=2) #WV1
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WI"][1:3], type="b", col=popPal[4], lwd=2) #WI1
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="MN"][1:3], type="b", col=popPal[5], lwd=2) #MN
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="AL"][1:3], type="b", col=popPal[6], lwd=2) #AL
# points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="IR"][1:3], type="b", col=popPal[7], lwd=2) #IR
# lines(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC2"][4:6], type="b", col=popPal[1], lwd=2) #NC2
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC1"][4:6], type="b", col=popPal[2], lwd=2) #NC1
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WV"][4:6], type="b", col=popPal[3], lwd=2) #WV1
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WI"][4:6], type="b", col=popPal[4], lwd=2) #WI1
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="MN"][4:6], type="b", col=popPal[5], lwd=2) #MN
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="AL"][4:6], type="b", col=popPal[6], lwd=2) #AL
# points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="IR"][4:6], type="b", col=popPal[7], lwd=2) #IR
# legend("topright", lty=1, col=popPal, legend=c("NC2","NC1","WV","WI","MN","AL","IR"), cex=1, bty="n", ncol=2)
# text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")
# 
# mtext("Growth to 5th instar (g/day)", 2, outer=T, at=0.6, cex=0.85, line=0.7)
# 
# dev.off()



png(here("./figures/Fig2_FifthGain.png"), units="in", res=300, width=6.5, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T))
par(mar=c(5.2,2.6,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(FifthGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(-0.001,0.021),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=1,y=-0.0005,"a")
text(x=4,y=-0.0005,"b")
text(x=5,y=-0.0005,"c")
text(x=c(2,3,6),y=-0.0005,"d")
#text(x=c(4,5,6),y=0.007,"c")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(FifthGain ~ Population, las=2, data=dat, ylab="", ylim=c(-0.001,0.021), col=popPal)
text(x=5,y=-0.0005,"a")
text(x=4,y=-0.0005,"b")
text(x=c(1,3,7),y=-0.0005,"c")
text(x=c(2,6),y=-0.0005,"bc")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0.003,0.014))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
lines(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC2"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC1"][1:3], type="b", col=popPal[2], lwd=2) #NC1
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WV"][1:3], type="b", col=popPal[3], lwd=2) #WV1
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WI"][1:3], type="b", col=popPal[4], lwd=2) #WI1
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="MN"][1:3], type="b", col=popPal[5], lwd=2) #MN
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="AL"][1:3], type="b", col=popPal[6], lwd=2) #AL
points(1:3,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="IR"][1:3], type="b", col=popPal[7], lwd=2) #IR
lines(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC2"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="NC1"][4:6], type="b", col=popPal[2], lwd=2) #NC1
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WV"][4:6], type="b", col=popPal[3], lwd=2) #WV1
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="WI"][4:6], type="b", col=popPal[4], lwd=2) #WI1
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="MN"][4:6], type="b", col=popPal[5], lwd=2) #MN
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="AL"][4:6], type="b", col=popPal[6], lwd=2) #AL
points(4:6,FifthGain.trtPop.emm$emmean[FifthGain.trtPop.emm$Population=="IR"][4:6], type="b", col=popPal[7], lwd=2) #IR
legend("top", lty=1, col=popPal, legend=c("NC2","NC1","WV","WI","MN","AL","IR"), cex=1, bty="n", ncol=7)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Growth to 5th instar (g/day)", 2, outer=T, at=0.6, cex=0.9, line=0.7)

dev.off()


aggregate(FifthTime ~ Treatment, data=dat, FUN=quantile, probs=0.95, na.rm=TRUE)



## look at survival by school
dat.ur<-dat[dat$School=="UR",]
dat.vcu<-dat[dat$School=="VCU",]

png(here("./figures/survival_by_lifestage_school.png"), units="in", res=300, width=3.5, height=3.5)
par(mar=c(5.1,4.1,1.1,1.1))
plot(NA,NA,xlim=c(1,5), ylim=c(0,525), xaxt="n", ylab="Number of survivors", xlab="")
axis(1, at=1:5, labels=c("Neonate","3rd instar","5th instar","Pupa","Adult"), las=2)
points(x=1:5,
       y=c(525
           ,sum(dat.ur$Survival=="T", na.rm=T) + sum(dat.ur$Survival=="F", na.rm=T) + sum(dat.ur$Survival=="P", na.rm=T) + sum(dat.ur$Survival=="A", na.rm=T)
           ,sum(dat.ur$Survival=="F", na.rm=T) + sum(dat.ur$Survival=="P", na.rm=T) + sum(dat.ur$Survival=="A", na.rm=T)
           ,sum(dat.ur$Survival=="P", na.rm=T) + sum(dat.ur$Survival=="A", na.rm=T)
           ,sum(dat.ur$Survival=="A", na.rm=T)
       )
       ,pch=19, col="red")
points(x=1:5,
       y=c(525
           ,sum(dat.vcu$Survival=="T", na.rm=T) + sum(dat.vcu$Survival=="F", na.rm=T) + sum(dat.vcu$Survival=="P", na.rm=T) + sum(dat.vcu$Survival=="A", na.rm=T)
           ,sum(dat.vcu$Survival=="F", na.rm=T) + sum(dat.vcu$Survival=="P", na.rm=T) + sum(dat.vcu$Survival=="A", na.rm=T)
           ,sum(dat.vcu$Survival=="P", na.rm=T) + sum(dat.vcu$Survival=="A", na.rm=T)
           ,sum(dat.vcu$Survival=="A", na.rm=T)
       )
       ,pch=19, col="black")
legend("bottomleft",pch=19,col=c("red","black"),legend=c("UR","VCU"), bty="n")
dev.off()




