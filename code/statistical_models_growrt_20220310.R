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
library(MuMIn)
library(effectsize)

dat.raw<-read.csv(here("../Data/Climate Simulation Exp/Climate Simulation Exp-Master Data.csv"),
                  stringsAsFactors = FALSE)

neonate.masses<-read.csv(here("../Data/neonate_masses.csv"))

str(dat.raw)
dat<-dat.raw

#correct some data entry anomalies
dat$AdultDOY[dat$AdultDOY < 150]<-NA
dat$Sex[dat$Sex=="FP"]<-"F"
dat$Sex[dat$Sex=="MP"]<-"M"

dat$SexDummy <- NA
dat$SexDummy[dat$Sex=="M"] <- -1
dat$SexDummy[dat$Sex=="U"] <- 0
dat$SexDummy[dat$Sex=="F"] <- 1

dat$Region <- NA
dat$Region[dat$Population %in% c("AL","IR","MN","WI1")] <- "North"
dat$Region[dat$Population %in% c("NC2","NC1","WV1")] <- "South"

any(is.na(dat$Region))

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


dat$Population[dat$Population=="WI1"]<-"WI"
dat$Population[dat$Population=="WV1"]<-"WV"

dat$Treatment<-factor(dat$Treatment, levels=c("WIB","WI4.5","WI8.5","VAB","VA4.5","VA8.5"), ordered=T)
dat$Population<-factor(dat$Population, levels=c("NC2","NC1","WV","WI","MN","AL","IR"), ordered=T)
dat$Region<- factor(dat$Region, levels=c("South","North"), ordered=T)

any(is.na(dat$Region))

dat2<-dat[dat$Sex!="U",]


trtPal<-magma(9)[3:8] #some parameters for plotting
popPal<-brewer.pal(3, "RdYlBu")[-2]


dat3<-dat[!is.na(dat$FifthTime),]

table(dat3$Treatment, dat3$Population)#/25


966/1050

# dat$FifthSurv <- ifelse(!is.na(dat$FifthGain),1,0)
# 
# Anova(glm(FifthSurv ~ Population + Treatment, data=dat, family="binomial"))
# summary(glm(FifthSurv ~ Population + Treatment, data=dat, family="binomial"))


##analyses on full dataset -------------------------------------------------------------------------------



cor.test(dat$ThirdGain, dat$ThirdTime)

ThirdGain.lmer<-lmer(ThirdGain ~ Treatment + Region + Treatment:Region + SexDummy + (1|Population), data=dat,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum"))
Anova(ThirdGain.lmer, type="III")
#summary(ThirdGain.lmer)
eta_squared(ThirdGain.lmer)

ThirdTime.lmer<-lmer(ThirdTime ~ Treatment + Region + Treatment:Region + SexDummy + (1|Population), data=dat,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum"))
Anova(ThirdTime.lmer, type="III")
eta_squared(ThirdTime.lmer)


ThirdGain.lmer2<-lmer(ThirdGain ~ Treatment + Region + Treatment:Region + Sex + (1|Population), data=dat2,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum", Sex="contr.sum"))
Anova(ThirdGain.lmer2, type="III")

ThirdTime.lmer2<-lmer(ThirdTime ~ Treatment + Region + Treatment:Region + Sex + (1|Population), data=dat2,
                      contrasts=list(Treatment='contr.sum', Region="contr.sum", Sex="contr.sum"))
Anova(ThirdTime.lmer2, type="III")


ThirdGain.trt.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Treatment))
ThirdGain.reg.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Region))
ThirdGain.trtReg.emm<-print(emmeans(ThirdGain.lmer, specs= ~ Treatment:Region))

ThirdTime.trt.emm<-print(emmeans(ThirdTime.lmer, specs= ~ Treatment))
ThirdTime.reg.emm<-print(emmeans(ThirdTime.lmer, specs= ~ Region))
ThirdTime.trtReg.emm<-print(emmeans(ThirdTime.lmer, specs= ~ Treatment:Region))



png(here("../Manuscript files/FigSX_ThirdGain_R1.png"), units="in", res=300, width=6, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T), widths=c(0.65,0.45))
par(mar=c(5.3,2.6,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(ThirdGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(0,0.0027),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=c(1,4),y=0.0001,"a")
text(x=c(2,3,5),y=0.0001,"b")
text(x=c(6),y=0.0001,"c")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(ThirdGain ~ Region, las=2, data=dat, ylab="", ylim=c(0,0.0027), col=popPal)
text(x=1,y=0.0001,"a")
text(x=2,y=0.0001,"b")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0.0005,0.002))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
arrows(x0=1:3, x1=1:3, y0=ThirdGain.trtReg.emm$lower.CL[ThirdGain.trtReg.emm$Region=="South"][1:3],
       y1=ThirdGain.trtReg.emm$upper.CL[ThirdGain.trtReg.emm$Region=="South"][1:3], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=1:3, x1=1:3, y0=ThirdGain.trtReg.emm$lower.CL[ThirdGain.trtReg.emm$Region=="North"][1:3],
       y1=ThirdGain.trtReg.emm$upper.CL[ThirdGain.trtReg.emm$Region=="North"][1:3], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(1:3,ThirdGain.trtReg.emm$emmean[ThirdGain.trtReg.emm$Region=="South"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,ThirdGain.trtReg.emm$emmean[ThirdGain.trtReg.emm$Region=="North"][1:3], type="b", col=popPal[2], lwd=2) #NC1

arrows(x0=4:6, x1=4:6, y0=ThirdGain.trtReg.emm$lower.CL[ThirdGain.trtReg.emm$Region=="South"][4:6],
       y1=ThirdGain.trtReg.emm$upper.CL[ThirdGain.trtReg.emm$Region=="South"][4:6], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=4:6, x1=4:6, y0=ThirdGain.trtReg.emm$lower.CL[ThirdGain.trtReg.emm$Region=="North"][4:6],
       y1=ThirdGain.trtReg.emm$upper.CL[ThirdGain.trtReg.emm$Region=="North"][4:6], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(4:6,ThirdGain.trtReg.emm$emmean[ThirdGain.trtReg.emm$Region=="South"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,ThirdGain.trtReg.emm$emmean[ThirdGain.trtReg.emm$Region=="North"][4:6], type="b", col=popPal[2], lwd=2) #NC1

legend("top", lty=1, col=popPal, legend=c("South","North"), cex=1, bty="n", ncol=2)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Growth to 3rd instar (g/day)", 2, outer=T, at=0.6, cex=0.9, line=0.9)

dev.off()



png(here("../Manuscript files/FigSX_ThirdTime_R1.png"), units="in", res=300, width=6, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T), widths=c(0.65,0.45))
par(mar=c(5.3,2.1,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(ThirdTime ~ Treatment, las=2, data=dat, ylab="", ylim=c(10,47),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=c(1,4),y=11,"a")
text(x=c(2,5),y=11,"b")
text(x=c(3),y=11,"c")
text(x=c(6),y=11,"d")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(ThirdTime ~ Region, las=2, data=dat, ylab="", col=popPal, ylim=c(10,47))
text(x=1,y=11,"a")
text(x=2,y=11,"b")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(12,35))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
arrows(x0=1:3, x1=1:3, y0=ThirdTime.trtReg.emm$lower.CL[ThirdTime.trtReg.emm$Region=="South"][1:3],
       y1=ThirdTime.trtReg.emm$upper.CL[ThirdTime.trtReg.emm$Region=="South"][1:3], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=1:3, x1=1:3, y0=ThirdTime.trtReg.emm$lower.CL[ThirdTime.trtReg.emm$Region=="North"][1:3],
       y1=ThirdTime.trtReg.emm$upper.CL[ThirdTime.trtReg.emm$Region=="North"][1:3], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(1:3,ThirdTime.trtReg.emm$emmean[ThirdTime.trtReg.emm$Region=="South"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,ThirdTime.trtReg.emm$emmean[ThirdTime.trtReg.emm$Region=="North"][1:3], type="b", col=popPal[2], lwd=2) #NC1

arrows(x0=4:6, x1=4:6, y0=ThirdTime.trtReg.emm$lower.CL[ThirdTime.trtReg.emm$Region=="South"][4:6],
       y1=ThirdTime.trtReg.emm$upper.CL[ThirdTime.trtReg.emm$Region=="South"][4:6], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=4:6, x1=4:6, y0=ThirdTime.trtReg.emm$lower.CL[ThirdTime.trtReg.emm$Region=="North"][4:6],
       y1=ThirdTime.trtReg.emm$upper.CL[ThirdTime.trtReg.emm$Region=="North"][4:6], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(4:6,ThirdTime.trtReg.emm$emmean[ThirdTime.trtReg.emm$Region=="South"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,ThirdTime.trtReg.emm$emmean[ThirdTime.trtReg.emm$Region=="North"][4:6], type="b", col=popPal[2], lwd=2) #NC1

legend("top", lty=1, col=popPal, legend=c("South","North"), cex=1, bty="n", ncol=2)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Development time to 3rd instar (days)", 2, outer=T, at=0.6, cex=0.9, line=0.7)

dev.off()


### Fifth Instar ----------------------------------------------------------------------------------

cor.test(dat$FifthGain, dat$FifthTime)

FifthGain.lmer<-lmer(FifthGain ~ Treatment + Region + SexDummy + Treatment:Region +
                     (1|Population), data=dat,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum"))
hist(residuals(FifthGain.lmer))
Anova(FifthGain.lmer, type="III")
#summary(FifthGain.lmer)
eta_squared(FifthGain.lmer)

FifthTime.lmer<-lmer(FifthTime ~ Treatment + Region + SexDummy + Treatment:Region
                     + (1|Population), data=dat,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum"))
hist(residuals(FifthTime.lmer))
Anova(FifthTime.lmer, type="III")
#summary(FifthTime.lmer)
eta_squared(FifthTime.lmer)



FifthGain.lmer2<-lmer(FifthGain ~ Treatment + Region + Sex + Treatment:Region +
                       (1|Population), data=dat2,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum", Sex="contr.sum"))
hist(residuals(FifthGain.lmer2))
Anova(FifthGain.lmer2, type="III")
#summary(FifthGain.lmer)


FifthTime.lmer2<-lmer(FifthTime ~ Treatment + Region + Sex + Treatment:Region
                     + (1|Population), data=dat2,
                     contrasts=list(Treatment='contr.sum', Region="contr.sum", Sex="contr.sum"))
hist(residuals(FifthTime.lmer2))
Anova(FifthTime.lmer2, type="III")
#summary(FifthTime.lmer)



# FifthGain.lm<-lm(FifthGain ~ Treatment + Region/Population + SexDummy + Treatment:Region, data=dat,
#                  contrasts=list(Treatment ='contr.sum', Region='contr.sum'))
# #plot(FifthGain.lm)
# hist(residuals(FifthGain.lm))
# anova(FifthGain.lm, test="LRT")
# Anova(FifthGain.lm, type="III")
# summary(FifthGain.lm)

FifthGain.trt.emm<-print(emmeans(FifthGain.lmer, specs= ~ Treatment))
FifthGain.reg.emm<-print(emmeans(FifthGain.lmer, specs= ~ Region))
FifthGain.trtReg.emm<-print(emmeans(FifthGain.lmer, specs= ~ Treatment:Region))

FifthTime.trt.emm<-print(emmeans(FifthTime.lmer, specs= ~ Treatment))
FifthTime.reg.emm<-print(emmeans(FifthTime.lmer, specs= ~ Region))
FifthTime.trtReg.emm<-print(emmeans(FifthTime.lmer, specs= ~ Treatment:Region))



png(here("../Manuscript files/Fig2_FifthGain_R1.png"), units="in", res=300, width=6, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T), widths=c(0.65,0.45))
par(mar=c(5.3,2.6,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(FifthGain ~ Treatment, las=2, data=dat, ylab="", ylim=c(-0.002, 0.021),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=c(1),y=-0.001,"a")
text(x=c(4),y=-0.001,"b")
text(x=c(5),y=-0.001,"c")
text(x=c(2,3,6),y=-0.001,"d")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(FifthGain ~ Region, las=2, data=dat, ylab="", col=popPal, ylim=c(-0.002, 0.021))
text(x=1,y=-0.001,"a")
text(x=2,y=-0.001,"b")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0.002,0.015))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
arrows(x0=1:3, x1=1:3, y0=FifthGain.trtReg.emm$lower.CL[FifthGain.trtReg.emm$Region=="South"][1:3],
       y1=FifthGain.trtReg.emm$upper.CL[FifthGain.trtReg.emm$Region=="South"][1:3], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=1:3, x1=1:3, y0=FifthGain.trtReg.emm$lower.CL[FifthGain.trtReg.emm$Region=="North"][1:3],
       y1=FifthGain.trtReg.emm$upper.CL[FifthGain.trtReg.emm$Region=="North"][1:3], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(1:3,FifthGain.trtReg.emm$emmean[FifthGain.trtReg.emm$Region=="South"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,FifthGain.trtReg.emm$emmean[FifthGain.trtReg.emm$Region=="North"][1:3], type="b", col=popPal[2], lwd=2) #NC1

arrows(x0=4:6, x1=4:6, y0=FifthGain.trtReg.emm$lower.CL[FifthGain.trtReg.emm$Region=="South"][4:6],
       y1=FifthGain.trtReg.emm$upper.CL[FifthGain.trtReg.emm$Region=="South"][4:6], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=4:6, x1=4:6, y0=FifthGain.trtReg.emm$lower.CL[FifthGain.trtReg.emm$Region=="North"][4:6],
       y1=FifthGain.trtReg.emm$upper.CL[FifthGain.trtReg.emm$Region=="North"][4:6], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(4:6,FifthGain.trtReg.emm$emmean[FifthGain.trtReg.emm$Region=="South"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,FifthGain.trtReg.emm$emmean[FifthGain.trtReg.emm$Region=="North"][4:6], type="b", col=popPal[2], lwd=2) #NC1

legend("top", lty=1, col=popPal, legend=c("South","North"), cex=1, bty="n", ncol=2)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Growth to 5th instar (g/day)", 2, outer=T, at=0.6, cex=0.9, line=0.7)

dev.off()



png(here("../Manuscript files/Fig3_FifthTime_R1.png"), units="in", res=300, width=6, height=4.5)

layout(matrix(c(1,2,3,3), 2, 2, byrow=T), widths=c(0.65,0.45))
par(mar=c(5.3,2.1,0.6,1.1), oma=c(0.1,2,0,0), mgp=c(4.3,0.6,0), cex.lab=1, tcl=-0.4)
boxplot(FifthTime ~ Treatment, las=2, data=dat, ylab="", ylim=c(16,60),
        names=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),
        col=trtPal)
text(x=c(1),y=18,"a")
text(x=c(4),y=18,"b")
text(x=c(2),y=18,"c")
text(x=5,y=18, "d")
text(x=3,y=18, "e")
text(x=6,y=18, "f")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"a)")

boxplot(FifthTime ~ Region, las=2, data=dat, ylab="", col=popPal, ylim=c(16,60))
text(x=1,y=18,"a")
text(x=2,y=18,"b")
text(par("usr")[1]+0.05*(par("usr")[2]-par("usr")[1]),0.93*par("usr")[4],"b)")

plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(25,55))
axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
arrows(x0=1:3, x1=1:3, y0=FifthTime.trtReg.emm$lower.CL[FifthTime.trtReg.emm$Region=="South"][1:3],
       y1=FifthTime.trtReg.emm$upper.CL[FifthTime.trtReg.emm$Region=="South"][1:3], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=1:3, x1=1:3, y0=FifthTime.trtReg.emm$lower.CL[FifthTime.trtReg.emm$Region=="North"][1:3],
       y1=FifthTime.trtReg.emm$upper.CL[FifthTime.trtReg.emm$Region=="North"][1:3], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(1:3,FifthTime.trtReg.emm$emmean[FifthTime.trtReg.emm$Region=="South"][1:3], type="b", col=popPal[1], lwd=2) #NC2
points(1:3,FifthTime.trtReg.emm$emmean[FifthTime.trtReg.emm$Region=="North"][1:3], type="b", col=popPal[2], lwd=2) #NC1

arrows(x0=4:6, x1=4:6, y0=FifthTime.trtReg.emm$lower.CL[FifthTime.trtReg.emm$Region=="South"][4:6],
       y1=FifthTime.trtReg.emm$upper.CL[FifthTime.trtReg.emm$Region=="South"][4:6], 
       length=0.07, angle=90, col=popPal[1], code=3)
arrows(x0=4:6, x1=4:6, y0=FifthTime.trtReg.emm$lower.CL[FifthTime.trtReg.emm$Region=="North"][4:6],
       y1=FifthTime.trtReg.emm$upper.CL[FifthTime.trtReg.emm$Region=="North"][4:6], 
       length=0.07, angle=90, col=popPal[2], code=3)
points(4:6,FifthTime.trtReg.emm$emmean[FifthTime.trtReg.emm$Region=="South"][4:6], type="b", col=popPal[1], lwd=2) #NC2
points(4:6,FifthTime.trtReg.emm$emmean[FifthTime.trtReg.emm$Region=="North"][4:6], type="b", col=popPal[2], lwd=2) #NC1

legend("top", lty=1, col=popPal, legend=c("South","North"), cex=1, bty="n", ncol=2)
text(par("usr")[1]+0.025*(par("usr")[2]-par("usr")[1]),0.95*par("usr")[4],"c)")

mtext("Development time to 5th instar (days)", 2, outer=T, at=0.6, cex=0.9, line=0.7)

dev.off()


#aggregate(FifthTime ~ Treatment, data=dat, FUN=quantile, probs=0.95, na.rm=TRUE)

## try looking at growth relative to baseline
# 
# FifthGain.trtPop.rel <- FifthGain.trtPop.emm[,1:3]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC2"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC2"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC2" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC1"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC1"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC1" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WV"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WV"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WV" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="AL"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="AL"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="AL" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="IR"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="IR"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="IR" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WI"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WI"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WI" & FifthGain.trtPop.rel$Treatment=="WIB"]
# FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="MN"] <- FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="MN"]/FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="MN" & FifthGain.trtPop.rel$Treatment=="WIB"]
# 
# plot(NA, NA, xlim=c(1,6), xlab="Treatment", ylab="", xaxt="n", las=2, ylim=c(0.5,2.5))
# axis(1,at=1:6,labels=c("Base WI","CC4.5WI","CC8.5WI","Base VA","CC4.5VA","CC8.5VA"),las=2)
# lines(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC2"], type="b", col=popPal[1], lwd=2) #NC2
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="NC1"], type="b", col=popPal[2], lwd=2) #NC1
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WV"], type="b", col=popPal[3], lwd=2) #WV1
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="WI"], type="b", col=popPal[4], lwd=2) #WI1
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="MN"], type="b", col=popPal[5], lwd=2) #MN
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="AL"], type="b", col=popPal[6], lwd=2) #AL
# points(1:6,FifthGain.trtPop.rel$emmean[FifthGain.trtPop.rel$Population=="IR"], type="b", col=popPal[7], lwd=2) #IR
# legend("top", lty=1, col=popPal, legend=c("NC2","NC1","WV","WI","MN","AL","IR"), cex=1, bty="n", ncol=7)

## look at survival by school
dat.ur<-dat[dat$School=="UR",]
dat.vcu<-dat[dat$School=="VCU",]

png(here("../Manuscript files/survival_by_lifestage_school.png"), units="in", res=300, width=3.5, height=3.5)
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




## Analyses on split datasets -------------------------------------------------------------------------------

dat.ur<-dat[dat$School=="UR",] #Richmond first



ThirdGain.lm.UR<-glm(ThirdGain ~ Sex + Treatment + Population + Treatment:Population, data=dat.ur, family=Gamma)
#plot(ThirdGain.lm)
hist(residuals(ThirdGain.lm.UR))
anova(ThirdGain.lm.UR, test="LRT") #interaction is significant; this is consistent with full model
summary(ThirdGain.lm.UR)



FifthGain.lm.UR<-glm(FifthGain ~ Sex + Treatment + Population + Treatment:Population, data=dat.ur, family=Gamma)
#plot(FifthGain.lm)
hist(residuals(FifthGain.lm.UR))
anova(FifthGain.lm.UR, test="LRT") #interaction is significant; this is consistent with full model
summary(FifthGain.lm.UR)




#now for VCU ----------------------------
dat.vcu<-dat[dat$School=="VCU",] #VCU



ThirdGain.lm.VCU<-glm(ThirdGain ~ Sex + Treatment + Population + Treatment:Population, data=dat.vcu, family=Gamma)
#plot(ThirdGain.lm)
hist(residuals(ThirdGain.lm.VCU))
anova(ThirdGain.lm.VCU, test="LRT") #interaction is not significant; this is not consistent with full model
#summary(ThirdGain.lm.VCU)


FifthGain.lm.VCU<-glm(FifthGain ~ Sex + Treatment + Population + Treatment:Population, data=dat.vcu, family=Gamma)
#plot(FifthGain.lm)
hist(residuals(FifthGain.lm.VCU))
anova(FifthGain.lm.VCU, test="LRT") #interaction is not significant; this is not consistent with full model
#summary(FifthGain.lm.VCU)
