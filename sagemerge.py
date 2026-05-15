# must execfile phot.py before this!

r=Region()
from astropy.table import Table
import os
import numpy as np
import pylab as pl
pl.ion() #interactive mode


do_UL=False # turn low SNR detections into ULs


regname='mergesage5'
tfile="mergesage5.blend10.00.8.csv"
tfile="mergesage5p.blend10.00.85.csv"

t=Table.read(tfile,format="ascii")

ra=t['ira'].data
de=t['idec'].data
sourcename=t['iname'].data

outdir=regname+".plots/"
if not os.path.exists(outdir):
    os.mkdir(outdir)

apphot_band=['Ha','K','I3.6','I8.0'] #,'M24']
apphot_wave=[0.656,2.16,3.6,8.0    ] #,24.]

mdir="/Users/ri3e/cv/magellanic/ir_mosaics/"
apphot_file=[mdir + s for s in [
    # in the first attempt, this was fix1, but that removes YSOs!
    'LMC.ha.csub.fix1.fits',
    'lmc.2massK.fits',
    'lmc.irac3.6_2_mosaic.fits',
    'lmc.irac8.0_2_mosaic.fits'
   # 'SAGE_LMC_MIPS24_E12.fits'
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

#dists=[1,5,10,20]
dists=[0.6,1,5]
ndist=len(dists)

#nsrc=5 # test


restart=False
if not restart:
    print("allocating arrays")
    app_f   =pl.zeros([nimage,nsrc,ndist])
    app_df  =pl.zeros([nimage,nsrc,ndist])
    app_raw =pl.zeros([nimage,nsrc,ndist])
    app_flag=pl.zeros([nimage,nsrc,ndist])
    app_bg  =pl.zeros([nimage,nsrc,ndist])
    
#for isrc in np.arange(nsrc-5099)+5099:
for isrc in range(nsrc):
    #print("starting",isrc)
    pl.clf() #clear figure
    pl.figure(1,figsize=[8,10])
    r.debug=False
    debugph=r.debug
    r.setbgfact([2,3]) # annulus

    # go through aperture radii
    for idist in range(ndist):

        f=[]
        df=[]
        raw=[]
        bg=[]
        
        r0=0.17 *dists[idist]/.05/3600 # arcsec F1000 at distance -> LMC deg
        # that's r=3.4asec at 1Mpc
        # if we put in dists=0.6, we'll get r=2" which is ~native for mcels
        
        r.setcircle([ra[isrc],de[isrc],r0])
        #print(r.imextents(apphot_ims[0]))
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

t=Table([ra,de,sourcename],names=["ra","de","name"])
for j in range(ndist):
    for i in range(nimage):
        t.add_columns([app_f[i,:,j],app_df[i,:,j],app_raw[i,:,j],app_flag[i,:,j],app_bg[i,:,j]],
                      names=["f_"+("d%i"%dists[j])+"_"+apphot_band[i],"df_"+("d%i"%dists[j])+"_"+apphot_band[i],
                             "raw_"+("d%i"%dists[j])+"_"+apphot_band[i],"flag_"+("d%i"%dists[j])+"_"+apphot_band[i],
                             "bg_"+("d%i"%dists[j])+"_"+apphot_band[i]])

t.write("mergesage5.dists.apphot.csv",overwrite=True)
pl.subplots_adjust(left=0.13,bottom=0.1)



    
    


