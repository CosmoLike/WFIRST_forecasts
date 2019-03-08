import numpy as np
from chainconsumer import ChainConsumer

def DES_WFIRST_single_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])

	d2[:,0]=d2[:,0]-0.0486
	d3[:,0]=d3[:,0]-0.0486
	d2[:,1]=d2[:,1]
	d3[:,1]=d3[:,1]
	
	Abbott = d1[:,(0,26)]
	DES = d2[:,(0,1)]
	WFIRST = d3[:,(0,1)]
	
	Abbott[:,1]=Abbott[:,1]*(Abbott[:,0]/0.3)**0.5
	Abbott_weights=d1[:,28]

	DES[:,1]=((DES[:,1])*(DES[:,0]/0.3)**0.5)-0.0109
	WFIRST[:,1]=((WFIRST[:,1])*(WFIRST[:,0]/0.3)**0.5)-0.0109
	
	c = ChainConsumer()
	c.add_chain(Abbott[:,(0,1)],weights=Abbott_weights,parameters=paranames,name =chainnames[0])
	c.add_chain(DES[20000:,(0,1)], name =chainnames[1])
	c.add_chain(WFIRST[20000:,(0,1)], name =chainnames[2])	
	c.configure(shade=True, kde=1.5, shade_alpha=0.2, colors=["g", "r","b"], linewidths=1.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(truth=[0.267,0.773],figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)

def DES_WFIRST_Planck_single_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	d4 = np.genfromtxt(filename[3])

	d2[:,0]=d2[:,0]-0.0486
	d3[:,0]=d3[:,0]-0.0486
	d2[:,1]=d2[:,1]
	d3[:,1]=d3[:,1]
	
	Abbott = d1[:,(0,26)]
	DES = d2[:,(0,1)]
	WFIRST = d3[:,(0,1)]
	Planck = d4[:,(0,8)]
	
	Abbott[:,1]=Abbott[:,1]*(Abbott[:,0]/0.3)**0.5
	Abbott_weights=d1[:,28]

	Planck[:,1]=Planck[:,1]*(Planck[:,0]/0.3)**0.5
	Planck_weights=d4[:,10]

	DES[:,1]=((DES[:,1])*(DES[:,0]/0.3)**0.5)-0.0109
	WFIRST[:,1]=((WFIRST[:,1])*(WFIRST[:,0]/0.3)**0.5)-0.0109
	
	c = ChainConsumer()
	c.add_chain(Abbott[:,(0,1)],weights=Abbott_weights,parameters=paranames,name =chainnames[0])
	c.add_chain(Planck[:,(0,1)],weights=Planck_weights, name =chainnames[3])
	c.add_chain(DES[20000:,(0,1)], name =chainnames[1])
	c.add_chain(WFIRST[20000:,(0,1)], name =chainnames[2])
		
	c.configure(shade=True, kde=1.5, shade_alpha=0.2, colors=["g","o","r","b"], linewidths=1.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(extents=[[0.2, 0.48], [0.7, 0.98]],figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)



def twochain_single_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[:,(0,1)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[:,(0,1)], name =chainnames[1])
	c.configure(shade=True,kde=2.0,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)
	
def threechain_single_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	zpivot=0.0
	d1[:,3]= d1[:,3]+(zpivot/(1+zpivot))*d1[:,4]
	d2[:,3]= d2[:,3]+(zpivot/(1+zpivot))*d2[:,4]
	d3[:,3]= d3[:,3]+(zpivot/(1+zpivot))*d3[:,4]

	c = ChainConsumer()
	c.add_chain(d1[20000:,(3,4)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[20000:,(3,4)], name =chainnames[1])
	c.add_chain(d3[20000:,(3,4)], name =chainnames[2])

	c.configure(shade=True,kde=2.0,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)


def twochain_multi_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[100000:,(0,1,6,3,4)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[100000:,(0,1,6,3,4)], name =chainnames[1])
	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	
	fig = c.plotter.plot(filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)



# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_DES_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_LCDM"]
# chainnames=[r"DES Y1 3x2 Abbott et al",r"DES Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier (WFIRST sys control)"]
# paranames=[r"$\Omega_m$", r"$S_8$"]
# DES_WFIRST_single_plot(filename,"WF_DES2.pdf",chainnames,paranames)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_DES_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_LCDM_sys"]
# chainnames=[r"DES Y1 3x2 Abbott et al",r"DES Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier no sys"]
# paranames=[r"$\Omega_m$", r"$S_8$"]
# DES_WFIRST_single_plot(filename,"WF_DES3.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_DES_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LCDM_3x2pt_sys"]
# chainnames=[r"DES Y1 3x2 Abbott et al",r"DES Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier (WFIRST sys control)"]
# paranames=[r"$\Omega_m$", r"$S_8$"]
# DES_WFIRST_single_plot(filename,"WF_DES3.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_LCDM"]
# chainnames=[r"DES Y1 3x2 Abbott et al",r"WFIRST Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier (WFIRST no sys)"]
# paranames=[r"$\Omega_m$", r"$S_8$"]
# DES_WFIRST_single_plot(filename,"WF_DES4.pdf",chainnames,paranames)


filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_DES_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_Planck_nolensing_LCDM"]
chainnames=[r"DES Y1 3x2 Abbott et al 2018",r"DES Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier (WFIRST sys control)", "Planck 2015 no lensing, as in Abbott et al 2018"]
paranames=[r"$\Omega_m$", r"$S_8$"]
DES_WFIRST_Planck_single_plot(filename,"WF_DES_Planck.pdf",chainnames,paranames)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_nosys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_LSST_3x2pt_nosys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LSST_3x2pt_nosys"]
# chainnames=[r"WFIRST HLS",r"LSST",r"WFIRST wide+LSST"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_LSST_orig.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/llike_LSST_3x2pt_photo_bias_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti"]
# chainnames=[r"WFIRST HLS sys",r"LSST sys",r"WFIRST wide+LSST sys"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_LSST_sys.pdf",chainnames,paranames)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_SN10_sys"]
# chainnames=[r"WFIRST all",r"WFIRST all sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$",r"$h_0$",r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"WF_ocelote_all_vs_allsys.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_3x2pt_SN10_MG","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_3x2pt_SN10"]
# chainnames=[r"WFIRST all",r"WFIRST all sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$",r"$h_0$",r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"WF_ocelote_MG_vsnonMG.pdf",chainnames,paranames)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote"]
# chainnames=[r"WFIRST multi-band",r"WFIRST single band+LSST"]
# paranames=[r"$w_0$", r"$w_a$"]
# twochain_single_plot(filename,"WF_ocelote.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_shear_shear"]
# chainnames=[r"WFIRST 3x2pt+clusterN+clusterWL",r"WFIRST WL"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$",r"$h_0$",r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"WF_ocelote_all_vs_shear.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_3x2pt","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_3x2pt_LSST_WFIRST"]
# chainnames=[r"WFIRST multi-band 3x2pt",r"WFIRST single band+LSST 3x2pt"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$",r"$h_0$",r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"WF_ocelote_LSST_WFIRST.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_ocelote_3x2pt"]
# chainnames=[r"WFIRST 3x2pt+clusterN+clusterWL",r"WFIRST 3x2pt"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$",r"$h_0$",r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"WF_ocelote_all_vs_3x2.pdf",chainnames,paranames)


# filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like_IFC/like_WFIRST_ALL_0cut_4.428000e+01_2.000000e+03_Rmin10_Ncl20_Ntomo10_no_sys","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like_IFC/like_WFIRST_nonifc_0cut_3.543000e+01_2.000000e+03_Rmin10_Ncl20_Ntomo10_no_sys"]
# chainnames=[r"WFIRST Cosmic Shear without IFC",r"WFIRST Cosmic Shear with IFC"]
# paranames=[r"$w_0$", r"$w_a$"]
# twochain_single_plot(filename,"WF_IFC_study.pdf",chainnames,paranames)






# filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_1.500000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_4.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG"]
# chainnames=[r"2000 deg$^2$", r"1500 deg$^2$", r"4000 deg$^2$"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"WF_area.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_1.000000e+04_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LSST_4.500000e+01_1.800000e+04_Rmin10_Ncl15_Ntomo10_no_sys_MG"]
# chainnames=["WFIRST std", r"WFIRST ext (10000 deg$^2$)", "WFIRST+LSST"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"WF_ext_LSST.pdf",chainnames,paranames)






