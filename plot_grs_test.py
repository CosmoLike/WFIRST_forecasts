from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#add location of data vector file for plotting
datavfile1 = "GRS_data_vector"
datavfile2 = "GRS_pred_vector"

d1 = np.genfromtxt(datavfile1)[:,3]
d2 = np.genfromtxt(datavfile2)[:,3]

zbin= np.genfromtxt(datavfile2)[:,0]
kvalues= np.genfromtxt(datavfile2)[:,1]
muvalues=np.genfromtxt(datavfile2)[:,2]

#define plot ranges
Nk=100
Nmu=10
Nz=7
plotmax=Nk*Nmu

plotfile = "plots/GRS_test.png"

#variance of the GRS data points
varfile = "./GRS_variance"
var = np.sqrt(np.genfromtxt(varfile)[:,3])
ndata = d1.shape[0]
#file with Y1 scale cuts
	
print "chi2 calculated"
chi =0.0
for i in range(0,ndata):
		chi +=(d1[i]-d2[i])*(d1[i]-d2[i])/var[i]
print "GRS: Delta chi2 = %f" %(chi)

print kvalues[0:Nk*Nmu:Nmu]

plt.figure(figsize=(6,6), dpi=1000)
fs = 18

plt.subplot(2,2,1)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1.0e+01,2.5e+03)
plt.xlim(0.002,0.33)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$P(k)$', fontsize = fs)
plt.errorbar(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],var[0:Nk*Nmu:Nmu],marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1)
plt.plot(kvalues[1:Nk*Nmu:Nmu],d1[1:Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1)
plt.plot(kvalues[2:Nk*Nmu:Nmu],d1[2:Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)
plt.plot(kvalues[3:Nk*Nmu:Nmu],d1[3:Nk*Nmu:Nmu],marker='o', color='orange',linestyle = '-',markersize = 1)
plt.plot(kvalues[4:Nk*Nmu:Nmu],d1[4:Nk*Nmu:Nmu],marker='o', color='brown',linestyle = '-',markersize = 1)
plt.plot(kvalues[5:Nk*Nmu:Nmu],d1[5:Nk*Nmu:Nmu],marker='o', color='cyan',linestyle = '-',markersize = 1)
plt.plot(kvalues[6:Nk*Nmu:Nmu],d1[6:Nk*Nmu:Nmu],marker='o', color='yellow',linestyle = '-',markersize = 1)
plt.plot(kvalues[7:Nk*Nmu:Nmu],d1[7:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1)
plt.plot(kvalues[8:Nk*Nmu:Nmu],d1[8:Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1)
plt.plot(kvalues[9:Nk*Nmu:Nmu],d1[9:Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)

plt.subplot(2,2,2)
#plt.yscale('log')
plt.ylim(3.0e+01,2.5e+03)
plt.xlim(0.002,0.33)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$P(k)$', fontsize = fs)
plt.errorbar(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],var[0:Nk*Nmu:Nmu],marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1.)
plt.plot(kvalues[1:Nk*Nmu:Nmu],d1[1:Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1.)
plt.plot(kvalues[2:Nk*Nmu:Nmu],d1[2:Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)
plt.plot(kvalues[3:Nk*Nmu:Nmu],d1[3:Nk*Nmu:Nmu],marker='o', color='orange',linestyle = '-',markersize = 1)
plt.plot(kvalues[4:Nk*Nmu:Nmu],d1[4:Nk*Nmu:Nmu],marker='o', color='brown',linestyle = '-',markersize = 1)
plt.plot(kvalues[5:Nk*Nmu:Nmu],d1[5:Nk*Nmu:Nmu],marker='o', color='cyan',linestyle = '-',markersize = 1)
plt.plot(kvalues[6:Nk*Nmu:Nmu],d1[6:Nk*Nmu:Nmu],marker='o', color='yellow',linestyle = '-',markersize = 1)
plt.plot(kvalues[7:Nk*Nmu:Nmu],d1[7:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1)
plt.plot(kvalues[8:Nk*Nmu:Nmu],d1[8:Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1)
plt.plot(kvalues[9:Nk*Nmu:Nmu],d1[9:Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)

plt.subplot(2,2,3)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1.0e+01,2.5e+03)
plt.xlim(0.002,0.33)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$P(k)$', fontsize = fs)
plt.errorbar(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],var[0:Nk*Nmu:Nmu],marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1)
plt.plot(kvalues[1:Nk*Nmu:Nmu],d1[1:Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1)
plt.plot(kvalues[2:Nk*Nmu:Nmu],d1[2:Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)
plt.plot(kvalues[3:Nk*Nmu:Nmu],d1[3:Nk*Nmu:Nmu],marker='o', color='orange',linestyle = '-',markersize = 1)
plt.plot(kvalues[4:Nk*Nmu:Nmu],d1[4:Nk*Nmu:Nmu],marker='o', color='brown',linestyle = '-',markersize = 1)
plt.plot(kvalues[5:Nk*Nmu:Nmu],d1[5:Nk*Nmu:Nmu],marker='o', color='cyan',linestyle = '-',markersize = 1)
plt.plot(kvalues[6:Nk*Nmu:Nmu],d1[6:Nk*Nmu:Nmu],marker='o', color='yellow',linestyle = '-',markersize = 1)

plt.subplot(2,2,4)
#plt.yscale('log')
plt.ylim(3.0e+01,2.5e+03)
plt.xlim(0.002,0.33)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$P(k)$', fontsize = fs)
plt.errorbar(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],var[0:Nk*Nmu:Nmu],marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(kvalues[0:Nk*Nmu:Nmu],d1[0:Nk*Nmu:Nmu],marker='o', color='r',linestyle = '-',markersize = 1.)
plt.plot(kvalues[1*Nk*Nmu:2*Nk*Nmu:Nmu],d1[1*Nk*Nmu:2*Nk*Nmu:Nmu],marker='o', color='b',linestyle = '-',markersize = 1.)
plt.plot(kvalues[2*Nk*Nmu:3*Nk*Nmu:Nmu],d1[2*Nk*Nmu:3*Nk*Nmu:Nmu],marker='o', color='g',linestyle = '-',markersize = 1)
plt.plot(kvalues[3*Nk*Nmu:4*Nk*Nmu:Nmu],d1[3*Nk*Nmu:4*Nk*Nmu:Nmu],marker='o', color='orange',linestyle = '-',markersize = 1)
plt.plot(kvalues[4*Nk*Nmu:5*Nk*Nmu:Nmu],d1[4*Nk*Nmu:5*Nk*Nmu:Nmu],marker='o', color='brown',linestyle = '-',markersize = 1)
plt.plot(kvalues[5*Nk*Nmu:6*Nk*Nmu:Nmu],d1[5*Nk*Nmu:6*Nk*Nmu:Nmu],marker='o', color='cyan',linestyle = '-',markersize = 1)
plt.plot(kvalues[6*Nk*Nmu:7*Nk*Nmu:Nmu],d1[6*Nk*Nmu:7*Nk*Nmu:Nmu],marker='o', color='yellow',linestyle = '-',markersize = 1)

plt.savefig(plotfile,dpi=1000)

# plt.figure(figsize=(8,8), dpi=400)
# fs = 18
# plt.subplot(4,2,1)
# plt.yscale('log')
# plt.ylim(2.e-7,1.2e-4)
# plt.xlim(0.001,0.3)
# #plt.title(r'$\xi_+$')
# plt.ylabel(r'$\xi_+$', fontsize = fs)
# plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
# plt.plot(kvalues,d1,marker='o', color='r',linestyle = '',markersize = 1.5)


# plt.subplot(4,2,2)
# plt.ylim(-0.25,0.25)
# plt.plot([0,1000],[0,0],linestyle ='--',color='k')
# plt.xlim(0,nxip-1)
# plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
# plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
# plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
# plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

# plt.subplot(4,2,3)
# plt.yscale('log')
# plt.ylim(2.e-7,6.e-5)
# plt.xlim(nxip,nxip+nxim-1)
# #plt.title(r'$\xi_-$')
# plt.ylabel(r'$\xi_-$', fontsize = fs)
# plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
# plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)

# plt.subplot(4,2,4)
# plt.ylim(-0.25,0.25)
# plt.plot([0,1000],[0,0],linestyle ='--',color='k')
# plt.xlim(nxip,nxip+nxim-1)
# plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
# plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
# plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
# plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

# plt.subplot(4,2,5)

# plt.yscale('log')
# plt.ylim(2.e-6,2.5e-3)
# plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
# #plt.title(r'$\gamma_t$')
# plt.ylabel(r'$\gamma_t$', fontsize = fs)
# plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.2)
# plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
# #plt.plot(ind,d3,linestyle = '-')

# plt.subplot(4,2,6)
# plt.ylim(-0.25,0.25)
# plt.plot([0,1000],[0,0],linestyle ='--',color='k')
# plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
# #plt.title(r'$\gamma_t$')
# plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
# plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
# plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
# plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)


# plt.subplot(4,2,7)
# plt.yscale('log')
# plt.ylim(1.e-4,0.6)
# plt.xlim(nxip+nxim+nggl,ndata)
# #plt.title(r'$w$')
# plt.ylabel(r'$w$', fontsize = fs)
# plt.xlabel(r'bin number', fontsize = fs)
# plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.4)
# plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
# #plt.plot(ind,d3,linestyle = '-')


# plt.subplot(4,2,8)
# plt.ylim(-0.04,0.04)
# plt.plot([0,1000],[0,0],linestyle ='--',color='k')
# plt.xlim(nxip+nxim+nggl,ndata)
# #plt.title(r'$w$')
# plt.xlabel(r'bin number', fontsize = 18)
# plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
# plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
# plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
# plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

# plt.tight_layout()
# plt.savefig(plotfile,dpi=400)
