import numpy as np
from chainconsumer import ChainConsumer

def sevenchain_multi_plot(start,filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	d3read = np.genfromtxt(filename[2])
	d4read = np.genfromtxt(filename[3])
	d5read = np.genfromtxt(filename[4])	
	d6read = np.genfromtxt(filename[5])
	d7read = np.genfromtxt(filename[6])	
	
	d1=d1read[start:,(3,4)]-fid
	d2=d2read[start:,(3,4)]-fid
	d3=d3read[start:,(3,4)]-fid
	d4=d4read[start:,(1,2)]-fid
	d5=d5read[start:,(3,4)]-fid
	d6=d6read[start:,(3,4)]-fid
	d7=d7read[start:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	c.add_chain(d3,name =chainnames[2])
	c.add_chain(d4,name =chainnames[3])
	c.add_chain(d5,name =chainnames[4])
	c.add_chain(d6,name =chainnames[5])
	c.add_chain(d7,name =chainnames[6])
	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['c','y','r','g','b','k','brown'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])


def sixchain_multi_plot(start,filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	d3read = np.genfromtxt(filename[2])
	d4read = np.genfromtxt(filename[3])
	d5read = np.genfromtxt(filename[4])	
	d6read = np.genfromtxt(filename[5])	
	
	d1=d1read[start:,(3,4)]-fid
	d2=d2read[start:,(3,4)]-fid
	d3=d3read[start:,(3,4)]-fid
	d4=d4read[start:,(1,2)]-fid
	d5=d5read[start:,(3,4)]-fid
	d6=d6read[start:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	c.add_chain(d3,name =chainnames[2])
	c.add_chain(d4,name =chainnames[3])
	c.add_chain(d5,name =chainnames[4])
	c.add_chain(d6,name =chainnames[5])
	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['c','y','r','g','b','k'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])

def threechain_multi_plot(start,filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	d3read = np.genfromtxt(filename[2])
	
	d1=d1read[start:,(3,4)]-fid
	d2=d2read[start:,(3,4)]-fid
	d3=d3read[start:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	c.add_chain(d3,name =chainnames[2])

	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['b','g','r'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])
	
def twochain_multi_plot(start,filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	
	d1=d1read[start[0]:,(3,4)]-fid
	d2=d2read[start[1]:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	
	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['b','g'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])

filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_shear_shear_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_clustering_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_clusterN_clusterWL_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_SN","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/GRS_Chris_nuisance","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_3x2pt_clusterN_clusterWL_sys_opti"]

chainnames=[r"WL",r"LSS",r"Clusters",r"SN",r"GRS",r"3x2pt",r"3x2+clusters"]
paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
plotrange=[(-0.5,0.5),(-1.4,1.4)]
sevenchain_multi_plot(500000,filename,"WFIRST_individual_vs_multi.pdf",chainnames,paranames,plotrange)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_shear_shear_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clustering_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_clusterN_clusterWL_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_SN","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_3x2pt_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_3x2pt_clusterN_clusterWL_sys_opti"]

# chainnames=[r"WL",r"LSS",r"Clusters",r"SN",r"3x2pt",r"3x2+clusters"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.5,0.5),(-1.4,1.4)]
# sixchain_multi_plot(500000,filename,"WFIRST_individual_vs_multi.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/LSSTC_emu/like/like_LSST_3x2pt_sys_pessi","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_LSST_3x2pt_obssys"]
# chainnames=[r"WF HLS",r"LSST",r"WF wide + LSST"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.4,0.4),(-1.3,1.3)]
# threechain_multi_plot(filename,"WFIRST_LSST_incl_sys.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_shear_shear_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_shear_shear_sys_pessi"]
# chainnames=[r"shear shear opti",r"shear shear pessi"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.6,0.6),(-1.6,1.6)]
# start=500000
# twochain_multi_plot(start,filename,"WF_shear_shear_pessi_vs_opti.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_nosys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/like_WFIRST_3x2pt_sys_opti"]
# chainnames=[r"3x2 WF HLS no sys",r"3x2 WF HLS 42-dim-sys"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.6,0.6),(-1.6,1.6)]
# start=[50000,500000]
# twochain_multi_plot(start,filename,"WF_3x2_nosys_vs_optisys.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_photosys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_photosys_pessi"]
# chainnames=[r"3x2 $\sigma_z =0.01$, $\sigma (\Delta_z)=0.002$, $\sigma (\sigma_z) =0.002$", r"3x2 $\sigma_z =0.05$, $\sigma (\Delta_z)=0.05$, $\sigma (\sigma_z) =0.02$"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.6,0.6),(-1.6,1.6)]
# start=500000
# twochain_multi_plot(start,filename,"WF_3x2_photo_pessi_vs_opti.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_shearsys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_3x2pt_shearsys_pessi"]
# chainnames=[r"3x2 $\sigma (m)=0.002$",r"3x2 $\sigma (m)=0.01$"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.6,0.6),(-1.6,1.6)]
# start=500000
# twochain_multi_plot(start,filename,"WF_3x2_shear_pessi_vs_opti.pdf",chainnames,paranames,plotrange)


# filename=["/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_flat_MORprior","/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/like_WFIRST_clusterN_clusterWL_sys_opti_std_MORprior"]
# chainnames=[r"CL opti",r"CL opti flat prior",r"CL opti std prior"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.6,0.6),(-1.6,1.6)]
# start=200000
# threechain_multi_plot(start,filename,"WFIRST_CL_MOR_study.pdf",chainnames,paranames,plotrange)

