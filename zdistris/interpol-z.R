library(stats)
mat<- read.table("zdistris/zdistri_WFIRST_LSST_lensing")
norm <- sum(0.2*mat[,4])
inter<-splinefun(mat[,2], mat[,4]/norm)

zmax <- 4.0
zmin <- 0.0

dz<-((zmax-zmin)/300)
col1<-formatC(seq(zmin, zmax-dz, dz),format="f",digits=5,zero.print=TRUE)
col2<-formatC(seq(zmin+0.5*dz, zmax-0.5*dz, dz),format="f",digits=5,zero.print=TRUE)
col3<-formatC(seq(zmin+dz, zmax, dz),format="f",digits=5,zero.print=TRUE)
new<-inter(col2)*dz
sum(new)
col4<-formatC(new,format="f",digits=5,zero.print=TRUE)

write(t(cbind(col1,col2,col3,col4)),"zdistris/zdistri_WFIRST_LSST_lensing_fine_bin",ncol=4)


#################################
## rebins zdistribution 
#################################
mat<- read.table("zdistribution_WFIRST_etc_wl_dn_dA_12_val_0")

#number of original bins
Nbin<-300
#number of bins merged
binmerge<-15

col1<-mat[,1]
col2<-mat[,2]
col3<-mat[,3]
col4<-mat[,4]

col1new <- col1[seq(1, length(col1), binmerge)]
col3new <- col3[seq(binmerge, length(col1), binmerge)]

col2new <- 0.5*(col3new-col1new)+col1new
col4new <- array(0, dim=length(col2new))
for (i in 1:length(col4new)) col4new[i]=sum(col4[((i-1)*binmerge+1):(i*binmerge)])

write(t(cbind(col1new,col2new,col3new,col4new)),"zdistribution_WFIRST_ETC_1",ncol=4)


barplot(col4new,space=0)

