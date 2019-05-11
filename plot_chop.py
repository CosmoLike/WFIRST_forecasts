import numpy as np
from chainconsumer import ChainConsumer

def chop_chain(filename, intervals):
	infile="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/"+filename[0]
	d1 = np.genfromtxt(infile)
	for i in range(0,len(intervals)-1):
		outfile="/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/"+filename[0]+"_"+str(intervals[i])+"_"+str(intervals[i+1])
		f = open(outfile,"w")
		for row in d1[intervals[i]:intervals[i+1],]:
			p_text = '  '.join(str(r) for r in row)
			f.write('%s\n' % (p_text))
		f.close() 
	

filename=["like_WFIRST_3x2pt_sys_pessi_long_wideproposal"]
intervals=[400000,600000,800000,1000000,1200000,1400000,1600000,1800000,2000000,2200000,2400000,2600000,2800000,3000000,3200000,3400000,3600000,3800000,4000000,4200000,4400000]
chop_chain(filename,intervals)





