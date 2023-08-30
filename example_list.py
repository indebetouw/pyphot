# must execfile phot.py before this!

r=Region()
from astropy.table import Table
import os
import pylab as pl
pl.ion() #interactive mode


# region name for output files
regname="example_test"


outdir=regname+".plots"
if not os.path.exists(outdir):
    os.mkdir(outdir)

apphot_band=['i1','i2','i3','i4']

apphot_wave=[3.55, 4.49, 5.73, 7.87] # microns

apphot_radius=[2.,3., 3., 4.] # arcsec

# what are the image files?
mdir="/Users/remy/dorado/mosaics/sage/"
apphot_file=[mdir + s for s in [
        'SAGE_LMC_IRAC3.6_1.2_mosaic.fits',
        'SAGE_LMC_IRAC4.5_1.2_mosaic.fits',
        'SAGE_LMC_IRAC5.8_1.2_mosaic.fits',
        'SAGE_LMC_IRAC8.0_1.2_mosaic.fits'
        ]]


# read in a list of positions
t = Table.read('example_list.csv', format='ascii')
ra=t['ra']
de=t['dec']
sourcename=t['sourcename']


#==============================================================
# set up apphot - load ims 

try:
    print len(apphot_ims)
except:
    print "loading images for aperture photometry"
    apphot_ims=loadims(apphot_file)
apphot_wave=pl.array(apphot_wave)


nimage=len(apphot_ims)
nsrc=len(ra)

app_f   =pl.zeros([nimage,nsrc])
app_df  =pl.zeros([nimage,nsrc])
app_raw =pl.zeros([nimage,nsrc])
app_flag=pl.zeros([nimage,nsrc])


for isrc in range(nsrc):
    print "starting",isrc
    pl.clf() #clear figure
    pl.figure(1,figsize=[8,10])
#    pl.show(block=False)
    r.debug=False
    debugph=r.debug


    # go through aperture radii
    k=0
    f=[]
    df=[]
    raw=[]
    panel=[6,6,0]  # number of subpanels
    while k<nimage:
        l=k
        while l<nimage and apphot_radius[l]==apphot_radius[k]: l=l+1
        r.setcircle([ra[isrc],de[isrc],apphot_radius[k]/3600.])
#        print "k=",k,l-1," r=",apphot_radius[k],apphot_radius[l-1]
#        if k>=12:
#            pdb.set_trace()
        panel[2]=2*k+1
        f2,df2,raw2=phot1(r,apphot_ims[k:l],names=apphot_band[k:l],panel=panel,debug=debugph)
        f=pl.concatenate([f,f2])
        df=pl.concatenate([df,df2])
        raw=pl.concatenate([raw,raw2])
        k=l
#    print "size of f,f2 = ",len(f),len(f2)

#    for k in range(nimage):
#        r.setcircle([ra[isrc],de[isrc],apphot_radius[k]/3600.])
#        f2,df2,raw2=phot1(r,apphot_ims[k],names=apphot_band[k],panel=[5,5,k+1],debug=debugph)
#        f=pl.concatenate([f,f2])
#        df=pl.concatenate([df,df2])
#        raw=pl.concatenate([raw,raw2])
        

# to mJy to match catalog phot
    f =1000*pl.array(f)
    df=1000*pl.array(df)
    raw=1000*pl.array(raw)


    # plot (also set flag for apphot)
    flag=pl.ones(nimage)

    z=pl.where(pl.isnan(f))[0]
    if len(z)>0:
        flag[z]=0
    # set Ul flags but not on the residual images w/wave=0
    z=pl.where((f<=df)*(apphot_wave>0))[0]
    if len(z)>0:
        flag[z]=3
    # set nondet flags 
    z=pl.where((f<=0)*(apphot_wave>0))[0]
    if len(z)>0:
        flag[z]=0

    pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)
    pl.savefig(outdir+ ("/%s" % str(sourcename[isrc]))+".cutouts.png")

    # save apphot to big array
    app_f[:,isrc]=f
    app_df[:,isrc]=df
    app_flag[:,isrc]=flag
    app_raw[:,isrc]=raw



    


import pickle
pickle.dump( {"ra":ra,
              "dec":de,
              "app_f":app_f,
              "app_df":app_df,
              "app_raw":app_raw,
              "app_flag":app_flag,
              "app_wave":apphot_wave,
              "app_files":apphot_file,
              }, open( outdir+regname+".pkl", "wb" ) )


    
    


