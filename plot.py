import numpy as np
from chainconsumer import ChainConsumer


def threechain_multi_plot(filename, out, chainnames, paranames):
	d1 = np.genfromtxt(filename[0])
	d2 = np.genfromtxt(filename[1])
	d3 = np.genfromtxt(filename[2])
	c = ChainConsumer()
	c.add_chain(d1[50000:,(0,8)], parameters=paranames, name =chainnames[0])
	c.add_chain(d2[50000:,(0,8)], name =chainnames[1])
	c.add_chain(d3[50000:,(0,8)], name =chainnames[2])
	c.configure(shade=True, shade_alpha=0.2, bar_shade=True)
	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
	fig = c.plotter.plot(figsize=2.0,filename="/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/plots/"+out)
	


filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3.300000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_5.400000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG"]
chainnames=[r"45 gal/arcmin$^2$", r"33 gal/arcmin$^2$", r"54 gal/arcmin$^2$"]
paranames=[r"$\Omega_m$", r"$\sigma_8$"]
threechain_multi_plot(filename,"WF_ngal.pdf",chainnames,paranames)

filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_1.500000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_4.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG"]
chainnames=[r"2000 deg$^2$", r"1500 deg$^2$", r"4000 deg$^2$"]
paranames=[r"$\Omega_m$", r"$\sigma_8$"]
threechain_multi_plot(filename,"WF_area.pdf",chainnames,paranames)

filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_1.000000e+04_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_LSST_4.500000e+01_1.800000e+04_Rmin10_Ncl15_Ntomo10_no_sys_MG"]
chainnames=["WFIRST std", r"WFIRST ext (10000 deg$^2$)", "WFIRST+LSST"]
paranames=[r"$\Omega_m$", r"$\sigma_8$"]
threechain_multi_plot(filename,"WF_ext_LSST.pdf",chainnames,paranames)






# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_12Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_16Mpc.txt"]
# chainnames=["8/8 Mpc", "8/12 Mpc", "8/16 Mpc"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"scale_cuts_info_loss.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_1H_conti.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_12Mpc_1H_conti.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_16Mpc_1H_conti.txt"]
# chainnames=["8/8 Mpc", "8/12 Mpc", "8/16 Mpc"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"scale_cuts_contaminated.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_AGN.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_CX.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_AD.txt"]
# chainnames=["AGN", "CX", "AD"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"baryons.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_AGN_geo.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_CX_geo.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_AD_geo.txt"]
# chainnames=["AGN", "CX", "AD"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"baryons_geo.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_lph_conti.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_sph_conti.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_shear_conti.txt"]
# chainnames=["lens photo-z wrong", "source photo-z wrong", "shear calib wrong"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"obs_sys.pdf",chainnames,paranames)

# def threechain_multi_plot(filename, out, chainnames, paranames):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	d3 = np.genfromtxt(filename[2])
# 	c = ChainConsumer()
# 	c.add_chain(d1[100000:,(0,17)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[100000:,(0,17)], name =chainnames[1])
# 	c.add_chain(d3[100000:,(0,17)], name =chainnames[2])
# 	c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
# 	fig = c.plot(filename="/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/plots/"+out,truth=[0.295,0.8344])

# filename=["/Users/teifler/Dropbox/1_eli_chain_exchange/1x2pt_LCDM_WL.txt","/Users/teifler/Dropbox/1_eli_chain_exchange/1x2pt_LCDM_WL_small_scales.txt","/Users/teifler/Dropbox/1_eli_chain_exchange/1x2pt_LCDM_WL_small_scales_AGN.txt"]
# chainnames=["WL", "WL incl small", "WL AGN incl small"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_multi_plot(filename,"WL.pdf",chainnames,paranames)


# def threechain_indivsmulti_plot(filename, out, chainnames, paranames):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	d3 = np.genfromtxt(filename[2])
# 	c = ChainConsumer()
# 	c.add_chain(d1[100000:,(0,17)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[100000:,(0,27)], name =chainnames[1])
# 	c.add_chain(d3[100000:,(0,27)], name =chainnames[2])
# 	c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
# 	fig = c.plot(filename="/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/plots/"+out,truth=[0.295,0.8344])

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/1x2pt_LCDM_WL.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/2x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc.txt"]
# chainnames=["Weak Lensing", "GGL+Clustering", "3x2pt"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_indivsmulti_plot(filename,"indivsmulti.pdf",chainnames,paranames)

# def twochain_multi_plot(filename, out, chainnames, paranames):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	c = ChainConsumer()
# 	c.add_chain(d1[100000:,(0,27)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[100000:,(0,27)], name =chainnames[1])
# 	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	
# 	fig = c.plot(filename="/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/plots/"+out,truth=[0.295,0.8344])


# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_flask_cov.txt"]
# chainnames=["CosmoLike cov", "FLASK cov"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_multi_plot(filename,"cov impact.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/2x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/2x2pt_LCDM_8Mpc_RSD_NonLimber.txt"]
# chainnames=["no RSD in data", "RSD in data"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_multi_plot(filename,"RSD2x2.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_RSD_NonLimber.txt"]
# chainnames=["no RSD in data", "RSD in data"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_multi_plot(filename,"RSD3x2.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_coyote.txt"]
# chainnames=["baseline", "Cosmic Emu in data"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_multi_plot(filename,"emu.pdf",chainnames,paranames)


# def twochain_multi_plot(filename, out, chainnames, paranames):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	c = ChainConsumer()
# 	c.add_chain(d1[100000:,(0,17)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[100000:,(0,27)], name =chainnames[1])
# 	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.])	
# 	fig = c.plot(filename="/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/plots/"+out,truth=[0.295,0.8344])

# filename=[,"/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/1x2pt_LCDM_WL.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/2x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_8Mpc_12Mpc.txt"]
# chainnames=["Weak Lensing", "GGL+Clustering", "3x2pt"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# threechain_indivsmulti_plot(filename,"indivsmulti.pdf",chainnames,paranames)



# def twochain_single_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])

# 	c = ChainConsumer()
# 	c.add_chain(d1[5000:,7:9], parameters=[r"$\Sigma_0$", r"$\mu_0$"], name ="")
# 	c.configure(shade=[True],shade_alpha=[0.2],linestyles=["--", "-"],linewidths=[0.5])	


# 	fig = c.plot(figsize=[7,6],filename='plots/'+out,truth=[0.0,0.0])
# 	fig.set_size_inches(8.0 + fig.get_size_inches())
# 	fig.show()

# filename=["/Users/teifler/Downloads/like_LSST_all_2pt_noredshifterr_nuisance.txt"]

# twochain_single_plot(filename,"chain_LSST_mu_Sigma_nosys.pdf")

# def twochain_multi_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,(0,1,2,3,4,27)], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$h_0$",r"$\sigma_8$"], name ="no contamination")
# 	c.add_chain(d2[20000:,(0,1,2,3,4,27)], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$h_0$",r"$\sigma_8$"], name ="1H-term contamination")
# 	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	


# 	fig = c.plot(filename='plots/'+out,truth=[0.295,2.260574e-09,0.9676,0.0468,0.6881,0.8344])
# 	fig.set_size_inches(8 + fig.get_size_inches())
	

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_1H_conti.txt"]
# twochain_multi_plot(filename,"scale_cuts8_8MPC.pdf")

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_12Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_12Mpc_1H_conti.txt"]
# twochain_multi_plot(filename,"scale_cuts8_12MPC.pdf")

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_16Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_16Mpc_1H_conti.txt"]
# twochain_multi_plot(filename,"scale_cuts8_16MPC.pdf")




# def twochain_multi_plot_new(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,(23,24,25,26)], parameters=[r"$A_1$", r"$A_2$", r"$A_3$",r"$A_4$"], name ="no contamination")
# 	c.add_chain(d2[20000:,(23,24,25,26)], parameters=[r"$A_1$", r"$A_2$", r"$A_3$",r"$A_4$"], name ="1H-term contamination")
# 	c.configure(shade=[True,False],shade_alpha=[0.2,0.2],linestyles=["--", "-"],linewidths=[0.5,1.])	


# 	fig = c.plot(filename='plots/'+out,truth=[0.0,0.0,0.0,0.0])
# 	fig.set_size_inches(8 + fig.get_size_inches())
	

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_1H_conti.txt"]
# twochain_multi_plot_new(filename,"scale_cuts8_8MPC_IA.pdf")

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_12Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_12Mpc_1H_conti.txt"]
# twochain_multi_plot_new(filename,"scale_cuts8_12MPC_IA.pdf")

# filename=["/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_16Mpc.txt","/Users/teifler/Dropbox/1_chain_exchange/cosmolike/3x2pt_LCDM_8Mpc_16Mpc_1H_conti.txt"]
# twochain_multi_plot_new(filename,"scale_cuts8_16MPC_IA.pdf")


# def twochain_walk_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,0:5], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$h_0$"])
# 	c.add_chain(d2[20000:,0:5], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$h_0$"])

# 	fig = c.plot_walks(filename='plots/'+out,truth=[0.295,2.260574e-09,0.9676,0.0468,0.6881], convolve=100)
# 	fig.set_size_inches(4.5 + fig.get_size_inches())
	
# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut3.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut3_contaminated_data.txt"]
# twochain_walk_plot(filename,"walk_scale_cuts8MPC.pdf")


# def threechain_multi_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	d3 = np.genfromtxt(filename[2])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,0:15], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$h_0$",r"$b1_1$",r"$b1_2$",r"$b1_3$",r"$b1_4$",r"$b1_5$",r"$b2_1$",r"$b2_2$",r"$b2_3$",r"$b2_4$",r"$b2_5$"], name ="4 Mpc")
# 	c.add_chain(d2[20000:,0:15], name ="4/8 Mpc")
# 	c.add_chain(d3[20000:,0:15], name ="8 Mpc")
# 	c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-"],linewidths=[0.5, 1., 1.])	


# 	fig = c.plot(filename='plots/'+out,truth=[0.295,2.260574e-09,0.9676,0.0468,0.6881,1.45,1.55,1.65,1.8,2.0,0.0,0.0,0.0,0.0,0.0])
# 	fig.set_size_inches(4.5 + fig.get_size_inches())
	
# def threechain2_multi_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	d3 = np.genfromtxt(filename[2])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,0:11], parameters=[r"$\Omega_m$", r"$A_s$", r"$n_s$",r"$\Omega_b$",r"$\omega_\nu$",r"$h_0$",r"$b1_1$",r"$b1_2$",r"$b1_3$",r"$b1_4$",r"$b1_5$"], name ="2 Mpc")
# 	c.add_chain(d2[20000:,0:11], name ="4 Mpc")
# 	c.add_chain(d3[20000:,0:11], name ="8 Mpc")
# 	c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-"],linewidths=[0.5, 1., 1.])	


# 	fig = c.plot(filename='plots/'+out,truth=[0.295,0.8344,0.9676,0.0468,0.0013,0.6881,1.45,1.55,1.65,1.8,2.0])
# 	fig.set_size_inches(4.5 + fig.get_size_inches())
	

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut2_nosys.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut3_nosys.txt"]

# threechain_multi_plot(filename,"chain_consumerY1_cuts_nosys.pdf")

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut2.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut3.txt"]

# threechain_multi_plot(filename,"chain_consumerY1_cuts.pdf")

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like_old/3x2pt_LCDM_cut2Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like_old/3x2pt_LCDM_cut4Mpc.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like_old/3x2pt_LCDM_cut8Mpc.txt"]

# threechain2_multi_plot(filename,"chain_consumerY1_cut_nob2.pdf")



# def threechain_single_plot(filename, out):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	d3 = np.genfromtxt(filename[2])

# 	c = ChainConsumer()
# 	c.add_chain(d1[20000:,0:2], parameters=[r"$\Omega_m$", r"$A_s$"], name ="base")
# 	c.add_chain(d2[20000:,0:2], parameters=[r"$\Omega_m$", r"$A_s$"], name ="b2 marg")
# 	c.add_chain(d3[20000:,0:2], parameters=[r"$\Omega_m$", r"$A_s$"], name ="b1b2 marg")

# 	c.configure(shade=[True,False,False],shade_alpha=[0.5,0.5,0.5],linestyles=["--", "-","-."],linewidths=[1.,1.,1.])	


# 	fig = c.plot(filename='plots/'+out,figsize="column",truth=[0.295,2.260574e-09])
# #	fig = c.plot(figsize=(12, 12))
# 	#fig.set_size_inches(4.5 + fig.get_size_inches())

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_AsOm.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_AsOmb2.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_AsOmb2b1.txt"]

# threechain_single_plot(filename,"testb2_1.pdf")

# filename=["/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_AsOm.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_b2bs_AsOmb2.txt","/Users/teifler/Dropbox/cosmolike/top-level/des_mpp/like/3x2pt_LCDM_scut1_nosys_v2_b2bs_AsOmb2b1.txt"]

# threechain_single_plot(filename,"testb2_2.pdf")



