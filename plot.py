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



def twochain_single_plot(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[200000:,(0,1)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[200000:,(0,1)], name =chainnames[1])
	c.configure(shade=True,kde=2.0,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(truth=true,figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)


def GRS_single_plot(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	# zpivot=0.0
	# d1[:,3]= d1[:,3]+(zpivot/(1+zpivot))*d1[:,4]
	# d2[:,3]= d2[:,3]+(zpivot/(1+zpivot))*d2[:,4]
	# d3[:,3]= d3[:,3]+(zpivot/(1+zpivot))*d3[:,4]
	# d4[:,3]= d4[:,3]+(zpivot/(1+zpivot))*d4[:,4]

	c = ChainConsumer()
	c.add_chain(d1[400000:,], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[400000:,], name =chainnames[1])
	c.add_chain(d3[400000:,], name =chainnames[2])
	

	c.configure(shade=True,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(truth=true,figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)

def GRS_cosmology_plot(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	# zpivot=0.0
	# d1[:,3]= d1[:,3]+(zpivot/(1+zpivot))*d1[:,4]
	# d2[:,3]= d2[:,3]+(zpivot/(1+zpivot))*d2[:,4]
	# d3[:,3]= d3[:,3]+(zpivot/(1+zpivot))*d3[:,4]
	# d4[:,3]= d4[:,3]+(zpivot/(1+zpivot))*d4[:,4]

	c = ChainConsumer()
	c.add_chain(d1[400000:,0:7], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[400000:,0:7], name =chainnames[1])
	c.add_chain(d3[400000:,0:7], name =chainnames[2])
	

	c.configure(shade=True,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(truth=true,figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)


def threechain_single_plot(filename, out, chainnames, paranames, true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	zpivot=0.0
	d1[:,3]= d1[:,3]+(zpivot/(1+zpivot))*d1[:,4]
	d2[:,3]= d2[:,3]+(zpivot/(1+zpivot))*d2[:,4]
	d3[:,3]= d3[:,3]+(zpivot/(1+zpivot))*d3[:,4]

	c = ChainConsumer()
	c.add_chain(d1[700000:800000,(3,4)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[800000:900000,(3,4)], name =chainnames[1])
	c.add_chain(d3[4500000:,(3,4)], name =chainnames[2])

	c.configure(shade=True,kde=2.0,shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(truth=true,figsize=2.0,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)


def twochain_multi_plot(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[400000:,(0,1,6,3,4)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[400000:,(0,1,6,3,4)], name =chainnames[1])
	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	
	fig = c.plotter.plot(truth=true,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)

def twochain_multi_plot_clusters_only(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[400000:,0:34], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[400000:,0:34], name =chainnames[1])
	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	
	fig = c.plotter.plot(truth=true,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)

def twochain_multi_plot_3x2pt_only(filename, out, chainnames, paranames,true):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	c = ChainConsumer()
	c.add_chain(d1[400000:,0:49], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[400000:,0:49], name =chainnames[1])
	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	
	fig = c.plotter.plot(truth=true,filename="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_long_wideproposal","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_short_narrowproposal"]
# chainnames=[r"WF 3x2 sys long wide proposal",r"WF 3x2 sys short narrow proposal"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_S$", r"$w_0$",r"$w_a$",r"$\Omega_b$",r"$h_0$",r"$b_1$",r"$b_2$",r"$b_3$",r"$b_4$",r"$b_5$",r"$b_6$",r"$b_7$",r"$b_8$",r"$b_9$",r"$b_{10}$",r"$spz_1$",r"$spz_2$",r"$spz_3$",r"$spz_4$",r"$spz_5$",r"$spz_6$",r"$spz_7$",r"$spz_8$",r"$spz_9$",r"$spz_{10}$",r"$spz \sigma$",r"$lpz_1$",r"$lpz_2$",r"$lpz_3$",r"$lpz_4$",r"$lpz_5$",r"$lpz_6$",r"$lpz_7$",r"$lpz_8$",r"$lpz_9$",r"$lpz_{10}$",r"$lpz \sigma$",r"$m_1$",r"$m_2$",r"$m_3$",r"$m_4$",r"$m_5$",r"$m_6$",r"$m_7$",r"$m_8$",r"$m_9$",r"$m_{10}$"]
# true=[0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.01,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.01,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
# twochain_multi_plot_3x2pt_only(filename,"like_WFIRST_3x2pt_sys_opti_long_vs.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_long_wideproposal","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_long_wideproposal","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_long_wideproposal"]
# chainnames=[r"WF 3x2 sys 3M",r"WF 3x2 sys 4M",r"WF 3x2 sys 5M"]
# paranames=[r"$w_0$",r"$w_p$"]
# true=[-1.,0.]
# threechain_single_plot(filename,"like_WFIRST_3x2pt_sys_opti_long_convtest.pdf",chainnames,paranames,true)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade1_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade2_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade1_sys"]
# chainnames=[ r"WFIRST GRS sys std",r"WFIRST GRS sys deep",r"WFIRST GRS sys wide"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$",r"$w_a$", r"$\Omega_b$", r"$h_0$", r"$b_1$", r"$b_2$", r"$b_3$", r"$b_4$", r"$b_5$", r"$b_6$", r"$b_7$", r"$\sigma p_1$", r"$\sigma p_2$", r"$\sigma p_3$", r"$\sigma p_4$", r"$\sigma p_5$", r"$\sigma p_6$", r"$\sigma p_7$",r"$\sigma_z$", r"$KSTAR$"]
# true=[0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.2,1.3,1.45,1.6,1.70,1.8,1.9,290.,290.,290.,290.,290.,290.,290.,0.001,0.24]
# GRS_single_plot(filename,"WF_GRS_sys.pdf",chainnames,paranames,true)

filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade2_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade3_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade4_sys"]
chainnames=[ r"WFIRST GRS std",r"WFIRST GRS deep",r"WFIRST GRS wide"]
paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$",r"$w_a$",r"$\Omega_b$", r"$h_0$"]
true=[0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727]
GRS_cosmology_plot(filename,"WF_GRS_nosys.pdf",chainnames,paranames,true)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade2","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade3","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade4"]
# chainnames=[r"WFIRST GRS std",r"WFIRST GRS deep",r"WFIRST GRS wide"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$",r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# true=[0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727]
# GRS_single_plot(filename,"WF_GRS_nosys.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi"]
# chainnames=[r"WFIRST HLS 3x2 all sys opti",r"WFIRST HLS 3x2 all sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# twochain_single_plot(filename,"WF_3x2_allsys_opti_vs_pessi.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_shearsys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_shearsys_pessi"]
# chainnames=[r"WFIRST HLS 3x2 shear sys opti",r"WFIRST HLS 3x2 shear sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# twochain_single_plot(filename,"WF_3x2_shearsys_opti_vs_pessi.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_photosys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_photosys_pessi"]
# chainnames=[r"WFIRST HLS 3x2 photo sys opti",r"WFIRST HLS 3x2 photo sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# twochain_single_plot(filename,"WF_3x2_photosys_opti_vs_photopessi.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti_MG"]
# chainnames=[r"WFIRST HLS 3x2 sys opti",r"WFIRST HLS 3x2 sys opti MG"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# twochain_single_plot(filename,"WF_3x2_sys_opti_vs_opti_MG.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_flat_MORprior","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior"]
# chainnames=[r"WF clusterN+clusterWL sys self-calib MOR",r"WF clusterN+clusterWL sys std MOR"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_S$", r"$w_0$",r"$w_a$",r"$\Omega_b$",r"$h_0$",r"$spz_1$",r"$spz_2$",r"$spz_3$",r"$spz_4$",r"$spz_5$",r"$spz_6$",r"$spz_7$",r"$spz_8$",r"$spz_9$",r"$spz_{10}$",r"$spz \sigma$",r"$m_1$",r"$m_2$",r"$m_3$",r"$m_4$",r"$m_5$",r"$m_6$",r"$m_7$",r"$m_8$",r"$m_9$",r"$m_{10}$",r"$A$",r"$B$",r"$C$",r"$\sigma_0$",r"$q_M$",r"$q_z$"]
# true=[0.3156,0.831,0.6727,-1.,0.,0.0491685,0.6727,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.01,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3.207, 0.993, 0.0, 0.456, 0.0, 0.0]
# twochain_multi_plot_clusters_only(filename,"like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior_multipanel.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_flat_MORprior","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior"]
# chainnames=[r"WF clusterN+clusterWL sys self-calib MOR",r"WF clusterN+clusterWL sys std MOR"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$h_0$", r"$w_0$",r"$w_a$"]
# true=[0.3156,0.831,0.6727,-1.,0.]
# twochain_multi_plot(filename,"like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_flat_MORprior","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior"]
# chainnames=[r"WF clusterN+clusterWL sys self-calib MOR",r"WF clusterN+clusterWL sys std MOR"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# true=[0.3156,0.831]
# twochain_single_plot(filename,"like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior_single.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti_MG"]
# chainnames=[r"WFIRST HLS 3x2 sys opti",r"WFIRST HLS 3x2 sys opti MG"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$h_0$", r"$w_0$",r"$w_a$"]
# true=[0.3156,0.831,0.6727,-1.,0.]
# twochain_multi_plot(filename,"WF_3x2_sys_opti_vs_opti_MG_multi_panel.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti_MG"]
# chainnames=[r"WFIRST HLS 3x2 sys opti",r"WFIRST HLS 3x2 sys opti MG"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# twochain_single_plot(filename,"WF_3x2_sys_opti_vs_opti_MG.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_pessi_MG"]
# chainnames=[r"WFIRST HLS 3x2 sys pessi",r"WFIRST HLS 3x2 sys pessi MG"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$h_0$", r"$w_0$",r"$w_a$"]
# true=[0.3156,0.831,0.6727,-1.,0.]
# twochain_multi_plot(filename,"WF_3x2_sys_pessi_vs_pessi_MG_multi_panel.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_clusterN_clusterWL_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti"]
# chainnames=[r"WFIRST HLS 3x2+clusters sys opti",r"WFIRST HLS 3x2 sys opti",r"WFIRST HLS clusters sys opti"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# threechain_single_plot(filename,"WF_3x2clusters_vs_3x2_vsclusters_sys_opti.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_clusterN_clusterWL_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_pessi"]
# chainnames=[r"WFIRST HLS 3x2+clusters sys pessi",r"WFIRST HLS 3x2 sys pessi",r"WFIRST HLS clusters sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# threechain_single_plot(filename,"WF_3x2clusters_vs_3x2_vsclusters_sys_pessi.pdf",chainnames,paranames,true)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi"]
# chainnames=[r"WFIRST HLS 3x2+clusters sys pessi",r"WFIRST HLS 3x2 sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# true=[-1.,0.]
# threechain_single_plot(filename,"WF_3x2clusters_vs_3x2_sys_pessi.pdf",chainnames,paranames,true)

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


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_DES_3x2pt_LCDM","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_DES_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LCDM_3x2pt_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_Abbott_Planck_nolensing_LCDM"]
# chainnames=[r"DES Y1 3x2 Abbott et al 2018",r"DES Y1 3x2, Fourier (WFIRST sys control)", r"WFIRST 3x2, Fourier (WFIRST sys control)", "Planck 2015 no lensing, as in Abbott et al 2018"]
# paranames=[r"$\Omega_m$", r"$S_8$"]
# DES_WFIRST_Planck_single_plot(filename,"WF_DES_Planck.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_nosys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi"]
# chainnames=[r"WFIRST HLS no sys",r"WFIRST HLS opti",r"WFIRST HLS pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_3x2_sys_no_vs_opti.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_nosys_pessi"]
# chainnames=[r"WFIRST HLS sys pessi",r"WFIRST HLS sys opti",r"WFIRST HLS nosys opti"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_3x2_sys_optivspessi.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_clusterN_clusterWL_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_clusterN_clusterWL_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_clusterN_clusterWL_nosys_pessi"]
# chainnames=[r"WFIRST HLS sys pessi",r"WFIRST HLS sys opti",r"WFIRST HLS nosys opti"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_3x2_clusterN_clusterWL_sys_optivspessi.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_shear_shear_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_clustering_sys_pessi"]
# chainnames=[r"WF 3x2 sys pessi",r"WF weak lensing sys pessi",r"WF clustering sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_single_vs_multi_pessi.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade2_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade3_sys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/trade4_sys"]
# chainnames=[r"WFIRST GRS std",r"WFIRST GRS deep",r"WFIRST GRS wide"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_GRS_sys.pdf",chainnames,paranames)



# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_pessi"]
# chainnames=[r"WFIRST HLS sys opti",r"WFIRST HLS sys pessi"]
# paranames=[r"$w_0$", r"$w_a$"]
# twochain_single_plot(filename,"WF_3x2_sys.pdf",chainnames,paranames)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_nosys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_LSST_3x2pt_nosys","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LSST_3x2pt_nosys"]
# chainnames=[r"WFIRST HLS",r"LSST",r"WFIRST wide+LSST"]
# paranames=[r"$w_0$", r"$w_a$"]
# threechain_single_plot(filename,"WF_LSST_orig.pdf",chainnames,paranames)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_LSST_3x2pt_photo_bias_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LSST_3x2pt_obssys"]
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






