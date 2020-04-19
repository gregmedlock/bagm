library(zoo) # AUC calculation


# read in the data from Biggs et al., ISME 2016
relativeIntensities = read.table('/home/glm5uh/asf_characterization/Metabolomics/NMR_Integrals/Integrals_UVA.csv',
                                 sep=',',header=T,row.names=NULL)
relativeIntensities = relativeIntensities[,-c(2,3,4,5,6)]
relativeIntensities[,1] = as.character(relativeIntensities[,1])

peakInfo = read.table('/home/glm5uh/asf_characterization/Metabolomics/NMR_Integrals/Peak_name_associations.txt',
                      sep='\t',header=T,row.names=NULL)

sampleInfo = read.table('/home/glm5uh/asf_characterization/Metabolomics/NMR_Integrals/sample_info_7Oct15.tsv',
                        sep='\t',header=F,row.names=NULL)
names(sampleInfo) = c("Grower","Medium","Rep","Round")
sampleInfo = sampleInfo[-253,]  # Last sample in "sample_info_7Oct15.tsv" is a blank buffer, so no need to keep it

### center all samples by the mean of the media in which they were grown
RI_media_centered = data.frame(matrix(ncol = 86, nrow = 0))
media = c(0,356,360,361,492,500,502,519)
species = as.character(c(356,360,361,492,500,502,519))
relativeIntensities.centered = relativeIntensities
for (s in species)
{
  # get the metabolite means for media produced by this species
  grower_means = apply(spentMedia[spentMediaInfo$Grower==s & spentMediaInfo$Medium=="0",],2,mean)
  
  # center samples from the normalizing condition itself
  #-c(1) to get rid of sample label
  relativeIntensities.centered[sampleInfo$Grower==s & sampleInfo$Medium=='0',-c(1)] = 
          relativeIntensities.centered[sampleInfo$Grower==s & sampleInfo$Medium=='0',-c(1)] - grower_means 
  
  # center the intensities for any sample that was grown in this media
  relativeIntensities.centered[sampleInfo$Medium==s,-c(1)] = relativeIntensities.centered[sampleInfo$Medium==s,-c(1)] - grower_means
}
  
RI_media_centered = relativeIntensities.centered
colnames(RI_media_centered) = c('Name',as.character(peakInfo$Name))
rowNames_RI_media_centered = RI_media_centered[,1]

#rename RI sample labels to match the growth curves
RI_media_centered$Name = paste0(sampleInfo$Grower,'in',sampleInfo$Medium,"_",sampleInfo$Rep)

### filter out the samples that didn't contain ASF519
RI_media_centered_519 = RI_media_centered[startsWith(RI_media_centered$Name,'519'),]

### filter out the samples that didn't contain ASF356
RI_media_centered_356 = RI_media_centered[startsWith(RI_media_centered$Name,'356'),]

### Get the outcome data. We are going to try the AUC for each sample's growth curve.
# This function interpolates all growth curves to a common reference.
correctTime <- function(od,t,timeRef)
{
  if (length(t) != length(timeRef))
  {
    fn = splinefun(t, y = od, method = "natural", ties = mean)
    interp_od = fn(timeRef)
    interp_od[interp_od<0.001] = 0.001
  }
  else
    interp_od = od
  return(interp_od)
}

# This function collects the 4 growth curves collected for each condition,
# calling correctTime to interpolate to a common reference time.
collectGrowthCurves <- function(fileName,timeRef)
{ 
  # initiate a matrix with one column and length(timeRef) rows
  remove_path = basename(fileName)
  condition = substr(remove_path,1,nchar(remove_path)-4)
  gc = read.table(fileName,sep='\t',header=F,row.names=NULL)
  gcdf = data.frame(gc)
  names(gcdf) = c("time","OD","group")
  t1 = gcdf$time[gcdf$group == 1]
  od1 = gcdf$OD[gcdf$group == 1]
  t2 = gcdf$time[gcdf$group == 2]
  od2 = gcdf$OD[gcdf$group == 2]
  t3 = gcdf$time[gcdf$group == 3]
  od3 = gcdf$OD[gcdf$group == 3]
  t4 = gcdf$time[gcdf$group == 4]
  od4 = gcdf$OD[gcdf$group == 4]
  # set zero/negative values as an arbitrarily low positive value.
  od1 = correctTime(od1,t1,timeRef)
  od2 = correctTime(od2,t2,timeRef)
  od3 = correctTime(od3,t3,timeRef)
  od4 = correctTime(od4,t4,timeRef)
  ODs = matrix(0,ncol=4,nrow=length(timeRef))
  ODs[,1] = od1
  ODs[,2] = od2
  ODs[,3] = od3
  ODs[,4] = od4
  ODs = data.frame(ODs)
  names(ODs) = paste0(condition,"_",c("1","2","3","4"))
  return(ODs)
}

# Process the non-spent media samples first, using the first replicate as the time reference.
fileName = "/home/glm5uh/asf_characterization/GrowthCurves/SpentMediaGrowthCurveData/ASF519in0.tsv"
gc = read.table(fileName,sep='\t',header=F,row.names=NULL)
gcdf = data.frame(gc)
names(gcdf) = c("time","OD","group")
t1 = gcdf$time[gcdf$group == 1]
ODs_519 = collectGrowthCurves(fileName,t1)
# od1 = gcdf$OD[gcdf$group == 1]
# t2 = gcdf$time[gcdf$group == 2]
# od2 = gcdf$OD[gcdf$group == 2]
# t3 = gcdf$time[gcdf$group == 3]
# od3 = gcdf$OD[gcdf$group == 3]
# t4 = gcdf$time[gcdf$group == 4]
# od4 = gcdf$OD[gcdf$group == 4]
# # set zero/negative values as an arbitrarily low positive value.
# od1 = correctTime(od1,t1,t1)
# od2 = correctTime(od2,t2,t1)
# od3 = correctTime(od3,t3,t1)
# od4 = correctTime(od4,t4,t1)
# ODs = matrix(0,ncol=4,nrow=length(t1))
# ODs[,1] = od1
# ODs[,2] = od2
# ODs[,3] = od3
# ODs[,4] = od4
# ODs = data.frame(ODs)
# names(ODs) = paste0("ASF519in0_",c("1","2","3","4"))

timeRef = t1
for (condition in paste0("ASF519in",c("356","360","361","492","500","502","519")))
{
  str(condition)
  fileName = paste0("/home/glm5uh/asf_characterization/GrowthCurves/SpentMediaGrowthCurveData/",condition,".tsv")
  new_ODs = collectGrowthCurves(fileName,timeRef=timeRef)
  ODs_519 = cbind(ODs_519,new_ODs)
}

# calculate AUCs
AUCs_519 = apply(ODs_519,2,function(x) sum(diff(timeRef)*rollmean(abs(x-0.001),2)))
# calculate max
maxes_519 = apply(ODs_519,2,function(x) max(x))


# Repeat growth curve processing for ASF356
fileName = "/home/glm5uh/asf_characterization/GrowthCurves/SpentMediaGrowthCurveData/ASF356in0.tsv"
gc = read.table(fileName,sep='\t',header=F,row.names=NULL)
gcdf = data.frame(gc)
names(gcdf) = c("time","OD","group")
t1 = gcdf$time[gcdf$group == 1]
ODs_356 = collectGrowthCurves(fileName,t1)

timeRef = t1
for (condition in paste0("ASF356in",c("356","360","361","492","500","502","519")))
{
  str(condition)
  fileName = paste0("/home/glm5uh/asf_characterization/GrowthCurves/SpentMediaGrowthCurveData/",condition,".tsv")
  new_ODs = collectGrowthCurves(fileName,timeRef=timeRef)
  ODs_356 = cbind(ODs_356,new_ODs)
}


#calculate AUC
AUCs_356 = apply(ODs_356,2,function(x) sum(diff(timeRef)*rollmean(abs(x-0.001),2)))
# calculate maxes
maxes_356 = apply(ODs_356,2,function(x) max(x))


# Rename rownames in metabolite DF to match AUCs
for (condition in unique(RI_media_centered_519$Name))
{
  RI_media_centered_519[RI_media_centered_519$Name == condition,]$Name = paste0("ASF",condition) 
}

for (condition in unique(RI_media_centered_356$Name))
{
  RI_media_centered_356[RI_media_centered_356$Name == condition,]$Name = paste0("ASF",condition) 
}


preprocess <- function(integrals,outcomes,outcome_name)
{
  #cbind with reordering of labels to match each other
  data = cbind(integrals,outcomes[integrals$Name])
  # rename the new column
  colnames(data)[length(colnames(data))] = outcome_name 
  # Scale the data. Min/max scaling seems to do the trick, since it preserves information about metabolite consumption/production.
  data[,!colnames(data) %in% c("Name",outcome_name)] = apply(data[,!colnames(data) %in% c("Name",outcome_name)],2,function(x) (x)/max(abs(x)))
  #scale the outcome as a z-score
  data[,outcome_name] = (data[,outcome_name] - mean(data[,outcome_name]))/sd(data[,outcome_name])
  return(data)
}
data_519 = data.frame(preprocess(RI_media_centered_519,AUCs_519,"AUC"))
data_356 = data.frame(preprocess(RI_media_centered_356,AUCs_356,"AUC"))
datamax_519 = data.frame(preprocess(RI_media_centered_519,maxes_519,"maxOD"))
datamax_356 = data.frame(preprocess(RI_media_centered_356,maxes_356,"maxOD"))
# remove the sample label column
#data = data[,-1]

# #top10 = names(tail(sort(abs(apply(RI_media_centered_519[,-1],2,sd))),15))
# data = RI_media_centered_519
# # Merge to add the AUC value
# data = cbind(data,AUCs_519[data$Name])
# colnames(data)[length(colnames(data))] = "AUC"
# # remove the sample label column
# data = data[,-1]
# # Scale the data. Min/max scaling seems to do the trick, since it preserves information about metabolite consumption/production.
# data = apply(data,2,function(x) (x)/max(abs(x)))

# get the top 10 most variable metabolites
top10_519 = names(tail(sort(apply(data_519[,!(colnames(data_519) %in% c("AUC","Name",colnames(data_519)[startsWith(colnames(data_519),"Unknown")]))],2,sd)),20))
top10_356 = names(tail(sort(apply(data_356[,!(colnames(data_356) %in% c("AUC","Name"))],2,sd)),20))

top10_519_max = names(tail(sort(apply(datamax_519[,!(colnames(datamax_519) %in% c("maxOD","Name",colnames(datamax_519)[startsWith(colnames(datamax_519),"Unknown")]))],2,sd)),20))
top10_356_max = names(tail(sort(apply(datamax_356[,!(colnames(datamax_356) %in% c("maxOD","Name",colnames(datamax_356)[startsWith(colnames(datamax_356),"Unknown")]))],2,sd)),20))

# select only the top 10 metabolites
data_519 = data_519[,c(top10_519,"AUC","Name")]
data_356 = data_356[,c(top10_356,"AUC","Name")]
datamax_519 = datamax_519[,c(top10_519_max,"maxOD","Name")]
datamax_356 = datamax_356[,c(top10_356_max,"maxOD","Name")]


# shortcut: lets use lm()
mod_519 = lm(AUC ~ . - Name, data=data_519)
# inspect the fit
plot(data_519$AUC,predict.lm(mod_519,data_519))

mod_519_max = lm(maxOD ~ . - Name, data=datamax_519)
# inspect the fit
plot(datamax_519$maxOD,predict.lm(mod_519_max,datamax_519))

# for the model, remove the sample with an erroneous spike in OD (ASF356in519_4)
data_356_prep = data_356[!(data_356$Name == "ASF356in519_4"),]
datamax_356_prep = datamax_356[!(datamax_356$Name == "ASF356in519_4"),]
mod_356 = lm(AUC ~ . - Name, data=data_356_prep)
plot(data_356_prep$AUC,predict.lm(mod_356,data_356_prep))

mod_356_max = lm(maxOD ~ . - Name, data=datamax_356_prep)
plot(datamax_356_prep$maxOD,predict.lm(mod_356_max,datamax_356_prep))


# actually use a bayesian model now
library(rethinking)


# declare the model
predictors <- reformulate(termlabels = c("a",paste(top10_519,"*",paste0("r",top10_519))), response = mu)[[3]] # [[3]] to get the RHS only
# construct the prior for each coefficient in the predictor formula
priors = vector('list',length(top10_519)) #preallocate the list for speed
for (i in 1:length(top10_519))
{
  prior = reformulate(termlabels = "dnorm(-1,1)", response = paste0("r",top10_519[i]))
  priors[[i]] <- prior
}

# Give up on handling formulas and just copy/paste the output from predictors into the assignment of mu.
coef_priors = lapply(paste0("r",top10_519),function(x) reformulate(termlabels = "dnorm(-1,1)", response = x))
metGauss <- map(
  append(alist(
    AUC ~ dnorm(mu, sigma),
    mu <- a + Leucine * rLeucine + Ethanol * rEthanol + Tryptophan * rTryptophan + Isovalerate * rIsovalerate + Propionate * rPropionate + Valine * rValine + Proline * rProline + Alanine * rAlanine + Methionine * rMethionine + Isoleucine * rIsoleucine + Acetate * rAcetate + Lactate * rLactate + Succinate * rSuccinate + Phenylalanine * rPhenylalanine + Lysine * rLysine + Tyrosine * rTyrosine + Butyrate * rButyrate + Betaine * rBetaine + Glycine * rGlycine + X3.hydroxybutyrate * rX3.hydroxybutyrate,
    a ~ dnorm(-10,10),
    sigma ~ dunif(0,10)
  ), priors) , 
  data = data_519, debug=TRUE, verbose = TRUE
)

mu = link(metGauss)
mu.mean = apply(mu,2,mean)
mu.PI = apply(mu,2,PI)
auc.sim = sim(metGauss,n=1e4)
auc.PI = apply(auc.sim, 2 , PI)

plot(mu.mean ~ data_519$AUC, col=rangi2, ylim=range(mu.PI),
     xlab='Observed AUC', ylab='Predicted AUC')
abline(a=0,b=1,lty=2)
for (i in 1:nrow(data_519))
{
  lines( rep(data_519$AUC[i],2), c(mu.PI[1,i],mu.PI[2,i]),col=rangi2)
}

# generate counterfactual plots for all metabolites from the 519 model
single_counterfact_plot <- function(metname, outcome, d, lenpred) {
  # get the mean metabolite values for all metabolites other than the target metabolite
  met.avgs = apply(d[,!colnames(d) %in% c(outcome,"Name",metname)],2,mean)
  # generate a range of values for the target metabolite, which we'd like to predict the outcome over
  targetmet = seq(from=-1, to=1, length.out=lenpred)
  met.avgs = do.call("rbind", replicate(lenpred, met.avgs, simplify = FALSE)) # match dimensions of the target variable
  pred.data = data.frame(targetmet.simulated=targetmet, nontarget.mets=met.avgs)
  colnames(pred.data) = c(metname,colnames(d)[!colnames(d) %in% c(outcome,"Name",metname)])
  
  # predict the mean and PI (credible interval, default = 0.89) of the coefficient
  mu = link(metGauss, data=pred.data)
  mu.mean = apply(mu,2,mean)
  mu.PI = apply(mu,2,PI)
  
  # sample from the posterior distribution to generate a distribution of real samples at each value of the target met
  R.sim = sim(metGauss, data=pred.data, n=1e4)
  R.PI = apply(R.sim, 2, PI)
  
  plot(d[,metname],d[,outcome] , pch = 16, xlab="", ylab="", 
                xlim=c(min(d[,metname]),max(d[,metname])),
                ylim=c(min(d[,outcome]),max(d[,outcome])))
  lines(targetmet, mu.mean)
  shade(mu.PI, targetmet)
  shade(R.PI, targetmet)
  #axis(2,las=2) # rotate y axis tick labels
  mtext(side=1, line=2, outcome, cex=2)
  mtext(side=2, line=3, paste0("Change in ",metname), cex=2)
}

for (metname in top10_519) {
  png(paste0('/home/glm5uh/bagm/results/counterfactual_AUC_',metname,".png"))
  single_counterfact_plot(metname,"AUC",data_519,50)
  dev.off()
}


