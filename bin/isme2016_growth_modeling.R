
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

### center all samples by the media mean
RI_media_centered = data.frame(matrix(ncol = 86, nrow = 0))
media = c(0,356,360,361,492,500,502,519)
species = c(356,360,361,492,500,502,519)
blankMedia = relativeIntensities[(sampleInfo$Grower==0 & sampleInfo$Medium==0),-1]
spentMedia = relativeIntensities[(sampleInfo$Grower>0 | sampleInfo$Medium>0),-1]
spentMediaInfo = sampleInfo[(sampleInfo$Grower>0 | sampleInfo$Medium>0),]
FMmeans = colMeans(blankMedia)

# center each sample by the media average
for(r in 1:nrow(spentMedia))
{
  # Rename 
  tmpName = paste0(spentMediaInfo$Grower[r],'in',spentMediaInfo$Medium[r])
  
  # Calculate centered value
  tmp_center = spentMedia[r,] - FMmeans
  tmpdf = data.frame(matrix(ncol = 86, nrow = 1))
  tmpdf[1,2:86] = tmp_center
  tmpdf[1,1] = tmpName
  RI_media_centered = rbind(RI_media_centered,tmpdf)
}
names(RI_media_centered) = c('Name',as.character(peakInfo$Name))
rowNames_RI_media_centered = RI_media_centered[,1]

### filter out the samples that didn't contain ASF519, the only strain we'll model for now
RI_media_centered_519 = RI_media_centered[startsWith(RI_media_centered$Name,'519'),]


### Get the outcome data. We are going to try the AUC for each sample's growth curve.
# This function interpolates all growth curves to a common reference.
correctTime <- function(od,t,timeRef)
{
  if (length(t) != length(timeRef))
  {
    fn = splinefun(t, y = od, method = "natural", ties = mean)
    od2 = fn(timeRef)
    od2[od2<0.001] = 0.001
  }
  return(od2)
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

auc_inhibition_mat = matrix(0,ncol=7,nrow=7)
for(i in asfIDs)
{
  fileName = paste0("/home/glm5uh/asf_characterization/SpentMediaGrowthCurveData/ASF",i,"in",0,".tsv")
  curve0_df = averageGrowthCurve(fileName)
  auc0 = sum(diff(curve0_df$time)*rollmean(abs(curve0_df$ODave-0.001),2))
  
  for(j in asfIDs)
  {
    fileName = paste0("SpentMediaGrowthCurveData\\ASF",i,"in",j,".tsv")
    curve_df= averageGrowthCurve(fileName)
    aucj = sum(diff(curve_df$time)*rollmean(abs(curve_df$ODave-0.001),2))
    inhibition = -(auc0-aucj)/auc0
    auc_inhibition_mat[match(i,asfIDs),match(j,asfIDs)] = inhibition
    
    # Visualize, just for QC
    print(paste(fileName,": ",inhibition))
    plot(curve0_df$time,curve0_df$ODave,col='red',type='l')
    lines(curve_df$time,curve_df$ODave,col='blue',type='l')
    
  }
}
colnames(auc_inhibition_mat) = paste0("Spent",asfIDs)
rownames(auc_inhibition_mat) = paste0("ASF",asfIDs)
write.table(auc_inhibition_mat, file = "auc_inhibitions_all_pairs.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

