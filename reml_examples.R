library(lme4)

#response = otu count
#fixed effects = T1Dstatus
#random effects = DaysG + Nulliparous + Age_LMP + BMI_conception + HLA.6DRML

mixed.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange), data = dragons)
summary(mixed.lmer)

#for endia:

mixed.lmer <- lmer(otucount ~ T1Dstatus*DaysG + T1Dstatus*Age_LMP + T1Dstatus*BMI_conception + (1|Nulliparous) + (1|HLA.6DRML))

or??????????

mixed.lmer <- lmer(otucount ~ T1Dstatus + DaysG + Age_LMP + BMI_conception + (1|Nulliparous) + (1|HLA.6DRML), reml=True)