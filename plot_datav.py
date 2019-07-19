from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#add location of data vector file for plotting
plotfile = "plots/WFIRST_version_comparison.png"
datavfile1 = "datav/WFIRST_3x2pt_clusterN_clusterWL_fid_opti_test"
datavfile2 = "datav/WFIRST_3x2pt_clusterN_clusterWL_fid_opti"


d1_s = np.genfromtxt(datavfile1)[:,1]
d2_s = np.genfromtxt(datavfile2)[:,1]
d1=np.asarray([i for i in d1_s if i>0.0])
d2=np.asarray([i for i in d2_s if i>0.0])
#use this covariance for redmagic lens sample
covfile = "./cov/cov_WFIRST_only_final"



# data = np.genfromtxt(covfile)
# cov =np.zeros((ndata,ndata))
# for i in range(0,data.shape[0]):
# 	cov[int(data[i,0]),int(data[i,1])] = data[i,8]+data[i,9]
# 	cov[int(data[i,1]),int(data[i,0])] = data[i,8]+data[i,9]
# 	if (int(data[i,0])-int(data[i,1])):
# 		cov[int(data[i,0]),int(data[i,1])]*= m[int(data[i,0])]*m[int(data[i,1])]  	
# 		cov[int(data[i,1]),int(data[i,0])]*= m[int(data[i,0])]*m[int(data[i,1])]  	
# #s now contains sqrt(cov[i,i])
# s = np.sqrt(np.diag(cov))
ind = np.arange(0,len(d1))
ndata=len(d1)
# chi =0.0
# inv = LA.inv(cov)
# for i in range(0,ndata):
# 	for j in range(0,ndata):
# 		chi +=(d1[i]-d2[i])*inv[i,j]*(d1[j]-d2[j])*m[i]*m[j]
# print("3x2pt: Delta chi2 = %f" %(chi))

# nshear = 1375
# nggl = 800
# nclustering = 250
# ncluster = ndata - nshear - nggl -nclustering
#s = np.sqrt(np.diag(cov))

plt.figure(figsize=(4,4), dpi=400)
fs = 18
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(0,ndata-1)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
#plt.errorbar(d1[],d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind,(d2-d1)/d2, marker='x', color='k',linestyle = '',markersize = 1.0)



plt.savefig(plotfile,dpi=400)
