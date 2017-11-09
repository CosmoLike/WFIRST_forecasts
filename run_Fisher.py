import sys
#sys.path.append('/home/teifler/CosmoLike/WFIRST_forecasts/')

from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from cosmolike_libs import * 

write_datav = lib.write_vector_wrapper
write_datav.argtypes = [ctypes.c_char_p,InputCosmologyParams, InputNuisanceParams]

def init_WFIRST(file_source_z,file_lens_z):
	initcosmo()
	initbins(15,20.0,5000.0,5000.0,10.0,7,10)
	initsurvey("WFIRST")
	initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","source")
	initclusters()
	initia("none","DEEP2")
	initpriors("none","none","PhotoBAO","none")
	initprobes("all_2pt")

def read_invcov(cov_filename):
	covfile = np.genfromtxt(cov_filename)
	ndata = int(np.max(covfile[:,0]))+1
	invcov = np.zeros((ndata,ndata))
	for i in range(0,covfile.shape[0]):
		invcov[int(covfile[i,0]),int(covfile[i,1])]=covfile[i,2]
	return invcov

def FM(FM_params,invcov,step_width =1.0):
	print "\n\n--------------------------------------------"
	print "Running Cosmolike inFisher Matrix Mode"
	print "--------------------------------------------\n"
	print "step width = %.2f*fiducial_sigma" % (step_width)
	ndata = invcov.shape[0]
	npar = len(FM_params)
	derivs = np.zeros((npar,ndata))
	FM = np.zeros((npar,npar))
	tmpfile ="datav_FM"
	cosmo_fid = InputCosmologyParams().fiducial()
	cosmo_sigma = InputCosmologyParams().fiducial_sigma()
	nuisance_fid = InputNuisanceParams().fiducial()
	file1 = "FM_datav"
	n = 0
	for p in FM_params:
		print "FM: evaluting derivative for parameter %s"%(p)
		cosmo_var = InputCosmologyParams().fiducial()
		cosmo_max = InputCosmologyParams().fiducial()
		p0 = getattr(cosmo_fid,p)
		dp = getattr(cosmo_sigma,p)*step_width

		setattr(cosmo_var, p, p0-2.*dp)
		write_datav(file1,cosmo_var,nuisance_fid)
		dv_mm = np.genfromtxt(file1)[:,1]

		setattr(cosmo_var, p, p0-dp)
		write_datav(file1,cosmo_var,nuisance_fid)
		dv_m = np.genfromtxt(file1)[:,1]

		setattr(cosmo_var, p, p0+dp)
		write_datav(file1,cosmo_var,nuisance_fid)
		dv_p = np.genfromtxt(file1)[:,1]

		setattr(cosmo_var, p, p0+2.*dp)
		write_datav(file1,cosmo_var,nuisance_fid)
		dv_pp = np.genfromtxt(file1)[:,1]

		#five point method for the first derivative
		derivs[n,:] = (-dv_pp +8.*dv_p -8.*dv_m+dv_mm)/(12.*dp)
		if (np.sum(np.abs(derivs[n,:]))==0):
			print "derivate is zero\nEXIT!\n"
			exit(1)
		n += 1
	for i in range(0,npar):
		for j in range(0,npar):
			FM[i,j] = np.dot(np.dot(invcov,derivs[i,:]),derivs[j,:])
	FMinv= LA.inv(FM)
	for n in range(0,npar):
		print "sigma(%s) = %e (umarg: %e)" % (FM_params[n],np.sqrt(FMinv[n,n]),1./np.sqrt(FM[n,n]))
	ind_w0 = FM_params.index("w0") if "w0" in FM_params else -1	
	ind_wa = FM_params.index("wa") if "wa" in FM_params else -1	
	if ((ind_w0 > -1) & (ind_wa > -1)):
		FOM = 1./np.sqrt(FMinv[ind_w0,ind_w0]*FMinv[ind_wa,ind_wa] -FMinv[ind_w0,ind_wa]*FMinv[ind_wa,ind_w0])
		print "FoM = %e" %(FOM)
	return FM
file_source_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
file_lens_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
cov_file = os.path.join(dirname, "cov/cov_3x2pt_4.500000e+01_1.800000e+04_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20_2pt_inv")


init_WFIRST(file_source_z,file_lens_z)

invcov = read_invcov(cov_file)
FM_params= sample_cosmology_only(MG=False)
FM1 = FM(FM_params,invcov,step_width =0.03)
#do it once more with different step width to check stability of derivatives
FM2 = FM(FM_params,invcov,step_width =0.04)

#cor1 = np.zeros((npar,npar))
#cor2 = np.zeros((npar,npar))
#for i in range(0,npar):
#	for j in range(0,npar):
#		cor1[i,j] = FM1[i,j]/np.sqrt(FM1[i,i]*FM1[j,j])
#		cor2[i,j] = FM2[i,j]/np.sqrt(FM2[i,i]*FM2[j,j])
#plt.figure()
#plt.subplot(2,1,1)
#plt.imshow(cor1, interpolation="nearest")
#plt.subplot(2,1,2)
#plt.imshow(cor2, interpolation="nearest")
#plt.colorbar()
#plt.show()


#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
