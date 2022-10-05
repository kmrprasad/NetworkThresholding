function CreatePowerLawRingNetwork(N,k,maxk,mut,muw,minc,maxc)

lfr_path='/data1/lfrwmx-master/build';
addpath(lfr_path);

if mut>1
	mut=mut/10;
end
if muw>1
	muw=mut;
end
lfrInfo.N=N;
lfrInfo.k=k;
lfrInfo.maxk=maxk;
lfrInfo.mut=mut;
lfrInfo.muw=muw;
lfrInfo.minc=minc;
lfrInfo.maxc=maxc;

% You need the toolbox from https://github.com/CarloNicolini
[corMatOrig,m]=lfrw_mx('N',N,'k',k,'maxk',maxk,'mut',mut,'muw',muw,'minc',minc,'maxc',maxc);
corMatOrigw=corMatOrig;
corMatOrig=double(corMatOrig>0);
corMat=nearestSPD(corMatOrig);
corMatchol=chol(corMat);

save( 'PowerLawNetMat.mat','corMat','m','corMatOrigw','corMatOrig','corMatchol','lfrInfo')

