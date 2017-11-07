#!/usr/bin/python
import sys
sys.path.insert(0, '/home/teifler/anaconda2/')
import math, numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np

# infile =['cov/cov_3x2pt_3.300000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_3.600000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_3.900000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.200000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.800000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_5.100000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_5.400000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_1.500000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_2.500000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_3.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_3.500000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_4.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_1.000000e+04_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_1.800000e+04_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20','cov/cov_3x2pt_4.500000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20']

infile =['cov/cov_3x2pt_cluster_4.500000e+01_2.000000e+03_WFIRST']

nggl = 45
ngcl = 26
nlens = 10
nlenscl=4 # 3 if set_cluster_DES, 4 if set_cluster_LSST 
nshear = 55
ncl=15
nclgcl=2
nrich=7 # 7 for WFIRST

ndata = (nshear+nggl+nlens)*ncl+nlenscl*nrich+nrich*ngcl*nclgcl #??? for WFIRST
n2pt = (nshear+nggl+nlens)*ncl #(55+?+4)*20 for WFIRST
ncluster = nlenscl*nrich #4*7 fir WFIRST
n2ptcl=n2pt+ncluster
nclusterN_WL=ncluster+nrich*ngcl*nclgcl
# ncluster =0
# n2ptcl=0
# nclusterN_WL=0

for k in range(0,15):
  	covfile = np.genfromtxt(infile[k])
	cov = np.zeros((ndata,ndata))

	print ndata,n2pt,int(np.max(covfile[:,0])+1)

	for i in range(0,covfile.shape[0]):
	  	cov[int(covfile[i,0]),int(covfile[i,1])] = covfile[i,8]+covfile[i,9]
	  	cov[int(covfile[i,1]),int(covfile[i,0])] = covfile[i,8]+covfile[i,9]
	 
	# for i in range(n2ptcl,ndata):
	#   for j in range(0,nshear*ncl):  
	#     cov[i,j]=0.0
	#     cov[j,i]=0.0

	# for i in range(0,covfile.shape[0]):
	#   cov[int(covfile[i,0]),int(covfile[i,1])] += covfile[i,9]
	#   cov[int(covfile[i,1]),int(covfile[i,0])] += covfile[i,9]

	cor = np.zeros((ndata,ndata))
	for i in range(0,ndata):
	    for j in range(0,ndata):
	       if (cov[i,i]*cov[j,j] >0):
	         cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])


	a = np.sort(LA.eigvals(cor))
	print "min+max eigenvalues full cor:"
	print np.min(a), np.max(a)
	print "neg eigenvalues full cor:"
	for i in range(0,ndata):
	  if (a[i]< 0.0): print a[i]


	# ############### invert shear covariance #################
	inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
	a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
	print "min+max eigenvalues shear cov:"
	print np.min(a), np.max(a)
	outfile = infile[k]+"_shear_inv"
	f = open(outfile, "w")
	for i in range(0,nshear*ncl):
	  for j in range(0,nshear*ncl):
	    f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()

	# ############### invert ggl covariance #################
	# inv = LA.inv(cov[nshear*ncl:(nshear+nggl)*ncl,nshear*ncl:(nshear+nggl)*ncl])
	# a = np.sort(LA.eigvals(cov[nshear*ncl:(nshear+nggl)*ncl,nshear*ncl:(nshear+nggl)*ncl]))
	# print "min+max eigenvalues ggl cov:"
	# print np.min(a), np.max(a)
	# outfile = infile[k]+"_ggl_inv"
	# f = open(outfile, "w")
	# for i in range(0,nggl*ncl):
	#   for j in range(0,nggl*ncl):
	#     f.write("%d %d %e\n" %(i,j, inv[i,j]))
	# f.close()

	# # ############### invert clustering covariance #################
	# inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
	# a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
	# print "min+max eigenvalues clustering cov:"
	# print np.min(a), np.max(a)
	# outfile = infile[k]+"_clustering_inv"
	# f = open(outfile, "w")
	# for i in range(0,nlens*ncl):
	#   for j in range(0,nlens*ncl):
	#     f.write("%d %d %e\n" %(i,j, inv[i,j]))
	# f.close()

	# ############### invert 2pt covariance #################
	a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
	print "min+max eigenvalues 2pt cov:"
	print np.min(a), np.max(a)
	inv = LA.inv(cov[0:n2pt,0:n2pt])
	outfile = infile[k]+"_2pt_inv" 
	f = open(outfile, "w")
	for i in range(0,n2pt):
	  for j in range(0,n2pt):
	    f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()


	# # ############### invert clusterN covariance #################
	# inv = LA.inv(cov[n2pt:n2pt+ncluster,n2pt:n2pt+ncluster])
	# a = np.sort(LA.eigvals(cov[n2pt:n2pt+ncluster,n2pt:n2pt+ncluster]))
	# print "min+max eigenvalues of clusterN matrix:"
	# print np.min(a), np.max(a)
	# if (np.min(a)<0):
	#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

	# outfile = infile[k]+"_clusterN_inv"
	# f = open(outfile, "w")
	# for i in range(0,ncluster):
	#   for j in range(0,ncluster):
	#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
	# f.close()


	# # ############### invert full2pt+clusterN+clusterWL covariance #################
	precond = 1.e-7
	for i in range(0,ncluster):
	  cov[n2pt+i,:]*= precond
	  cov[:,n2pt+i]*= precond
	inv = LA.inv(cov)
	a = np.sort(LA.eigvals(cov))
	print "min+max eigenvalues of full 2ptclusterN+clusterWL pre-conditioned matrix:"
	print np.min(a), np.max(a)
	if (np.min(a)<0):
	  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])
	for i in range(0,ncluster):
	  inv[n2pt+i,:]*= precond
	  inv[:,n2pt+i]*= precond

	outfile = infile[k]+"_2pt_clusterN_clusterWL_inv"
	f = open(outfile, "w")
	for i in range(0,ndata):
	  for j in range(0,ndata):
	    f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()



	# # ############### invert clusterN+clusterWL covariance #################
	inv = LA.inv(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL])
	a = np.sort(LA.eigvals(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL]))
	print "min+max eigenvalues of clusterN_WL pre-conditioned matrix:"
	print np.min(a), np.max(a)
	if (np.min(a)<0):
	  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

	for i in range(0,ncluster):
	  inv[i,:]*= precond
	  inv[:,i]*= precond

	outfile = infile[k]+"_clusterN_clusterWL_inv"
	f = open(outfile, "w")
	for i in range(0,nclusterN_WL):
	  for j in range(0,nclusterN_WL):
	    f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()

	# ############### invert 2pt+clusterN covariance #################
	# inv = LA.inv(cov[0:n2pt+ncluster,0:n2pt+ncluster])
	# a = np.sort(LA.eigvals(cov[0:n2pt+ncluster,0:n2pt+ncluster]))
	# print "min+max eigenvalues of 2pt+N matrix:"
	# print np.min(a), np.max(a)
	# if (np.min(a)<0):
	#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

	# for i in range(0,ncluster):
	#   inv[n2pt+i,:]*= precond
	#   inv[:,n2pt+i]*= precond
	# outfile = infile[k]+"_2pt_clusterN_inv"
	# f = open(outfile, "w")
	# for i in range(0,n2pt+ncluster):
	#   for j in range(0,n2pt+ncluster):
	#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
	# f.close()


	# cor = np.zeros((ndata,ndata))
	# for i in range(0,ndata):
	#     for j in range(0,ndata):
	#       if (cov[i,i]*cov[j,j] >0):
	#         cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

	# plt.figure()
	# #plt.imshow(np.log10(cov[0:1500,2000:]), interpolation="nearest",vmin=-25, vmax=-10)
	# plt.imshow(np.log10(cov), interpolation="nearest",vmin=-25, vmax=-10)
	# #plt.imshow(cor[n2ptcl:n2ptcl+200,300:nshear*ncl], interpolation="nearest",vmax=0.5)
	# plt.colorbar()
	# plt.show()
