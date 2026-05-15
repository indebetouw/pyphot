# must execfile phot.py before this!

r=Region()
from astropy.table import Table
import os
import numpy as np
import pylab as pl
pl.ion() #interactive mode


do_UL=False # turn low SNR detections into ULs


# region name for output files
regname='bica3'
tfile="bica08.table3.dat"
mytp="CN"

regname='bica5'
tfile="bica08.table5.dat"
mytp="NC"
#mytp="NA"

t=Table.read(tfile,format="ascii.fixed_width_no_header",col_starts=(0,32,35,38,42,46,49,52,56,62,69),col_ends=(31,34,37,41,45,48,51,55,61,68,72),names=("Name","rah","ram","ras","ded","dem","des","typ","bmaj","bmin","pa"))
#LMC
#t=t[638:]

# select type and LMC
cni= np.array([i for i in range(len(t)) if (mytp in t['typ'][i] and t['rah'][i]>3)])
t=t[cni]

ra=(t['rah'].data+t['ram'].data/60+t['ras'].data/3600)*15
de=(t['ded'].data-t['dem'].data/60-t['des'].data/3600)
rmaj=t['bmaj'].data*60/2 # arcsec
rmin=t['bmin'].data*60/2
radasec=np.sqrt(rmaj*rmin)
pa=t['pa'].data
sourcename=t['Name'].data

outdir=regname+".plots/"
if not os.path.exists(outdir):
    os.mkdir(outdir)

apphot_band=['u','Ha','K','I3.6','I8.0','M24']
apphot_wave=[0.35,0.656,2.16,3.6,8.0,24.]

mdir="../ir_mosaics/"
apphot_file=[mdir + s for s in [
    
    'LMC.ha.csub.fix2.fits',
    'lmc.2massK.fits',
    'lmc.irac3.6_2_mosaic.fits',
    'lmc.irac8.0_2_mosaic.fits',
    'SAGE_LMC_MIPS24_E12.fits'
]]


#=============================================================
# set up apphot - load ims 

try:
    print(len(apphot_ims))
except:
    print("loading images for aperture photometry")
    apphot_ims=loadims(apphot_file)
apphot_wave=pl.array(apphot_wave)


nimage=len(apphot_ims)
nsrc=len(ra)

dists=[1,5,10,20] # nominally Mpc, where 1Mpc is the native Spitzer at 50kpc -> JWST at 1Mpc
outcsv="bica."+mytp+".dists.apphot.csv"

dists=[1]
outcsv="bica."+mytp+".apphot.csv"

rpc=0.16*1e6/206265*np.array(dists)
ndist=len(dists)

app_f   =pl.zeros([nimage,nsrc,ndist])
app_df  =pl.zeros([nimage,nsrc,ndist])
app_raw =pl.zeros([nimage,nsrc,ndist])
app_flag=pl.zeros([nimage,nsrc,ndist])
app_bg  =pl.zeros([nimage,nsrc,ndist])

nsrc=10 # test

for isrc in range(nsrc):
    #print("starting",isrc)
    pl.clf() #clear figure
    pl.figure(1,figsize=[8,10])
    r.debug=False
    debugph=r.debug
    r.setbgfact([1.2,1.4]) # super narrow annulus

    # go through aperture radii
    for idist in range(ndist):

        f=[]
        df=[]
        raw=[]
        bg=[]
        
        minr=rpc[idist]*4 # -> to arcsec in LMC
        r0=rmaj[isrc]
        if r0<minr: r0=minr
        r1=rmin[isrc]
        if r1<minr: r1=minr
        r0/=3600
        r1/=3600
        
        if '-' in pa[isrc]:
            thispa=0
            r.setcircle([ra[isrc],de[isrc],r0])
        else:
            thispa=float(pa[isrc])
            r.setellipse([ra[isrc],de[isrc],r0,r1,thispa])
        print(r.imextents(apphot_ims[0]))
        xmin=np.min(r.imextents(apphot_ims[0]))
        if xmin<0 or np.isnan(xmin): 
            continue
        
        for k in range(nimage):
            f2,df2,raw2,bg2=phot1(r,[apphot_ims[k]],names=[apphot_band[k]],panel=[4,4,2*k+1],debug=debugph)
            f=pl.concatenate([f,f2])
            df=pl.concatenate([df,df2])
            raw=pl.concatenate([raw,raw2])
    
        # to mJy
        f  =1000*pl.array(f)
        df =1000*pl.array(df)
        raw=1000*pl.array(raw)
    
        # plot (also set flag for apphot)
        flag=pl.ones(nimage)
    
        if do_UL:
            z=pl.where((f>raw)*(apphot_wave>0))[0]
            if len(z)>0:
                flag[z]=3
            # set nondet flags 
            z=pl.where((f<=0)*(apphot_wave>0))[0]
            if len(z)>0:
                flag[z]=0        
        
        pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)
    
        # save apphot to big array
        app_f[:,isrc,idist]=f
        app_df[:,isrc,idist]=df
        app_flag[:,isrc,idist]=flag
        app_raw[:,isrc,idist]=raw
        # app_bg[:,isrc]=bg
        app_bg[:,isrc,idist]= raw-f # in mJy units
    
        
        pl.subplot(3,2,6)
        z=np.where(flag==1)[0]
        pl.errorbar(apphot_wave[z],f[z],yerr=df[z],fmt='.',label="apphot")
        #pl.plot(apphot_wave[z],raw[z],'b.')
        z=np.where(flag==3)[0]
        pl.plot(apphot_wave[z],f[z],'bv')
    
        pl.xscale("log")
        pl.yscale("log")
    
        pl.legend(loc="best",prop={"size":8})
        pl.xlabel("um")
        pl.ylabel("mJy")
    
            
    
        pl.subplot(6,3,16)
        pl.plot([0,1],[0,1],',')
        pl.text(0,0.5,"%10.6f %10.6f"%(ra[isrc],de[isrc]))
        ax=pl.gca()
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        pl.savefig(outdir+ ("/%s" % str(sourcename[isrc]))+".cutouts.png")


# MCELS in 1e-15 e/s/cm2

t=Table([ra,de,rmaj,rmin,pa,sourcename],names=["ra","de","rmaj","rmin","pa","name"])
for j in range(ndist):
    for i in range(nimage):
        t.add_columns([app_f[i,:,j],app_df[i,:,j],app_raw[i,:,j],app_flag[i,:,j],app_bg[i,:,j]],
                      names=["f_"+("d%i"%dists[j])+"_"+apphot_band[i],"df_"+("d%i"%dists[j])+"_"+apphot_band[i],
                             "raw_"+("d%i"%dists[j])+"_"+apphot_band[i],"flag_"+("d%i"%dists[j])+"_"+apphot_band[i],
                             "bg_"+("d%i"%dists[j])+"_"+apphot_band[i]])

t.write(outcsv,overwrite=True)
pl.subplots_adjust(left=0.13,bottom=0.1)

import pickle
pickle.dump( {"ra":ra,
              "dec":de,
              "name":sourcename,
              "app_f":app_f,
              "app_df":app_df,
              "app_raw":app_raw,
              "app_flag":app_flag,
              "app_bg":app_bg,
              "app_wave":apphot_wave,
              "app_files":apphot_file,
              }, open( outdir+"/"+regname+".pkl", "wb" ) )


    
    


