
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

### filter out the samples that didn't contain ASF519, the only strain we'll model for now
RI_media_centered_519 = RI_media_centered[startsWith(RI_media_centered$Name,'519'),]


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
od1 = gcdf$OD[gcdf$group == 1]
t2 = gcdf$time[gcdf$group == 2]
od2 = gcdf$OD[gcdf$group == 2]
t3 = gcdf$time[gcdf$group == 3]
od3 = gcdf$OD[gcdf$group == 3]
t4 = gcdf$time[gcdf$group == 4]
od4 = gcdf$OD[gcdf$group == 4]
# set zero/negative values as an arbitrarily low positive value.
od1 = correctTime(od1,t1,t1)
od2 = correctTime(od2,t2,t1)
od3 = correctTime(od3,t3,t1)
od4 = correctTime(od4,t4,t1)
ODs = matrix(0,ncol=4,nrow=length(t1))
ODs[,1] = od1
ODs[,2] = od2
ODs[,3] = od3
ODs[,4] = od4
ODs = data.frame(ODs)
names(ODs) = paste0("ASF519in0_",c("1","2","3","4"))

timeRef = t1
for (condition in paste0("ASF519in",c("356","360","361","492","500","502","519")))
{
  str(condition)
  fileName = paste0("/home/glm5uh/asf_characterization/GrowthCurves/SpentMediaGrowthCurveData/",condition,".tsv")
  new_ODs = collectGrowthCurves(fileName,timeRef=timeRef)
  ODs = cbind(ODs,new_ODs)
  str(names(ODs))
}


library(zoo)

#calculate AUCs
AUCs = apply(ODs,2,function(x) sum(diff(timeRef)*rollmean(abs(x-0.001),2)))

# Rename rownames in metabolite DF to match AUCs
for (condition in unique(RI_media_centered_519$Name))
{
  RI_media_centered_519[RI_media_centered_519$Name == condition,]$Name = paste0("ASF",condition) 
}

# get the top 10 most abundant metabolites based on average change from all conditions

top10 = names(tail(sort(abs(apply(RI_media_centered_519[,-1],2,sd))),15))
data = RI_media_centered_519
# Merge to add the AUC value
data = cbind(data,AUCs[data$Name])
colnames(data)[length(colnames(data))] = "AUC"
# remove the sample label column
data = data[,-1]
# Scale the data. Min/max scaling seems to do the trick, since it preserves information about metabolite consumption/production.
data = apply(data,2,function(x) (x)/max(abs(x)))

# get the top 10 most variable metabolites
top10 = names(tail(sort(apply(data[,!(colnames(data) %in% c("AUC",colnames(data)[startsWith(colnames(data),"Unknown")]))],2,sd)),20))

# select only the top 10 metabolites
data = data[,c(top10,"AUC")]


# shortcut: lets use lm()
data = data.frame(data)
mod = lm(AUC ~ ., data=data)

# inspect the fit
plot(data$AUC,predict.lm(mod,data))

# actually use a bayesian model now
library(rethinking)


# declare the model
metGauss <- map(
  alist(
    AUC ~ dnorm(mu, sigma),
    mu <- a+bR*Leucine+bA*Alanine+bG*Unknown39,
    a ~ dnorm(-10,10),
    bR ~ dnorm(-1,1),
    bA ~ dnorm(-1,1),
    bG ~ dnorm(-1,1),
    sigma ~ dunif(0,10)
  ) , 
  data = data
)
