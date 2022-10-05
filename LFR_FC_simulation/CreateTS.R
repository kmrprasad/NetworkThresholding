
require(neuRosim)
library("oro.nifti")
library("R.matlab")


args <- commandArgs(trailingOnly = TRUE)
args[2]
TemplateNifti=readNIfTI('TemplateEpi.nii') #Brain Template with used Crossley (2013) and it 638 areas in our case we resized it to an EPI scale and used only the 600 first areas

dim <- dim(TemplateNifti)

##To change if necessary
nscan <- 250 #nb of point for the TS
TR <- 2
total.time <- nscan*TR


####onsets

##correlation Matrix to get the time series from
data<-readMat("PowerLawNetMat.mat")



#create resting state time serie for each areas/nodes
#simTSrestingstate from the package neuRosim
restAll<-NULL
for (iRegion in 1:ncol(data$corMat)){
restAll<-rbind(restAll,t(simTSrestingstate(nscan, base=0, TR, SNR=NULL, noise = "none",verbose = TRUE)))
}


# Specify the value of each ROI of the templates .... could be loaded from a file
indVal<-seq(1,ncol(data$corMat),1)



##Adjust time series for them to fit our connectivity matrix
L =data$corMatchol
r = t(L) %*% restAll
r = t(r)
 
rdata = as.data.frame(r)
#cor(rdata)
#data$corMat
####



####Design
baseline<-apply(TemplateNifti,1:3,mean)
baseline.bin<-ifelse(baseline>0,1,0)
ix<-which(baseline==1)
baseline[-ix]<-0


# Prepare for the different noise
resp <- runif(1, 12, 20)  #repiratory
heart.beat <- runif(1, 60, 90)
valSNR<-as.numeric(args[1])
wnoise<-c(0.5,0.5,0,0,0,0)

#Create the simulated volumes with only the noise
# simVOLfmri from the package neuRosim
if(valSNR>0){  
print('before out')
ptm <- proc.time()
out <- simVOLfmri(design=list(),  dim=dim(TemplateNifti),nscan=nscan, base=100, SNR=valSNR, noise="mixture", type="rician", TR=2,rho.temp = 0.2, rho.spat = 0.4, w=wnoise, template=baseline.bin, freq.heart = heart.beat/60,freq.resp = resp/60)
proc.time() - ptm
print('after out')
}else{   # empty volumes
tmpdim<-dim(baseline)
out<-array(0, dim=c(tmpdim[1],tmpdim[2],tmpdim[3],nscan))
}

########Inject timeseries in the volumes

#Change here if you want just one time series by area or if you want to simulate really all the voxels
nbVoxByAreas<-1  #dim(which(TemplateNifti>0, arr.ind=TRUE))[1]

timeSeries<- matrix(data=NA,nrow=ncol(data$corMat)*(nbVoxByAreas),ncol=nscan)
ind<-1;
for (iRegion in 1:ncol(data$corMat)){

vox<-which(TemplateNifti==indVal[iRegion], arr.ind=TRUE)
#for (iVox in 1:dim(vox)[1]){
for (iVox in 1:nbVoxByAreas){
	timeSeries[ind,]<-rdata[,iRegion]+out[vox[iVox,1],vox[iVox,2],vox[iVox,3],]
	ind<-ind+1
}
}


filename <- paste('./',args[2],'/Subject',args[3],".mat", sep="")
print(filename)
writeMat(filename,timeSeries=timeSeries,Areas=indVal,nbVoxByArea=nbVoxByAreas,infoHeart=heart.beat/60,infoResp=resp/60,RealTS=rdata,ExpectedCorMat=data$corMat,RealCorMat=cor(rdata),NoiseParams=wnoise)
print('Finished createTS')

