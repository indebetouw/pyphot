# must execfile phot.py before this!

r=Region()
from astropy.table import Table
import os
import numpy as np
import pylab as pl
pl.ion() #interactive mode


do_UL=False # turn low SNR detections into ULs

phangs_apcorr=np.array([1.94,2.31,2.15,2.00,2.17, 1.78,1.87,1.86,1.85, 2,2,2,2])


# region name for output files
regname='Phot_test_10_90_ul_off_ac_off_f_ext_off_annulus_2-2.5ap_mjyngc628c'
regname='Phot_test_10_90_ul_off_ac_off_f_ext_off_annulus_2-4ap_mjy_sclip_5iterngc628c'
regname='Phot_test_10_90_ul_off_ac_off_mjy_ngc628c'
regname='Phot_test_10_90_ul_off_ac_off_f_ext_off_annulus_2-2.5ap_mjy_sclip_5iterngc628c'
jtable=Table.read(regname+'.csv',format="ascii")
j_apcorr=np.ones(13)

#regname="ngc628.100Myr"
#jtable=Table.read('intermediate.csv',format="ascii")
#j_apcorr=phangs_apcorr

u=np.argsort(jtable['ID_phangs'])
jtable=jtable[u]

ra=jtable['raj2000']
de=jtable['dej2000']
sourcename=jtable['ID_phangs']

filter_name=np.array(['F275W','F336W','F435W','F555W','F814W','F200W','F300M','F335M','F360M','F770W','F1000W','F1130W','F2100W'])

jfcols=['flux_'+x for x in filter_name]
jdfcols=['er_flux_9010_'+x for x in filter_name]
# special case, no 9010:
jdfcols[7]='er_flux_F335M'
# this is in the whole annulus - would better as surface brightness etc
jfbgcols=['flux_bkg_'+x for x in filter_name]

jfbcols=['flux_'+x+'_b' for x in filter_name]


hsttable=Table.read("IR4_ngc628c_human_subsample.csv",format="ascii")
hstfcols= ['PHANGS_'+x+'_mJy_TOT' for x in filter_name[0:5]]
hstdfcols=['PHANGS_'+x+'_mJy_TOT_ERR' for x in filter_name[0:5]]



        

outdir=regname+".plots/"
if not os.path.exists(outdir):
    os.mkdir(outdir)


apphot_band=[x[:-1] for x in filter_name]
apphot_wave=[.275,.336,.435,.555,.814, 2,3,3.35,3.6, 7.7,10,11.3,21]
# 658: r=.076
# 50%EE radii
apphot_radius=np.array([.11,.067,.076,.068,.081, .042,.061,.070,.073, .17,.21,.24,.42])*2
# "standard" radii adopted by jimena
apphot_radius=np.array([.158,.158,.158,.158,.158, .124,.124,.124,.124, .168,.21,.24,.42])

# jr
# for all HST bands an  apertures of 0.158488" plus aperture corrections of
# ac_555= -0.72               *1.94
# ac_275 = ac_555-0.19 =-.91  *2.31
# ac_336 = ac_555-0.12 =-.83  *2.15
# ac_435 = ac_555-0.03 =-.75  *2.00
# ac_814 = ac_555-0.12 =-.84  *2.17
# For all  NIRCam bands I used an aperture of  0.124" plus aperture corrections of
# ac_200=-0.63
# ac_300=-0.68
# ac_360=-0.67


# what are the image files?
#mdir="/Users/ri3e/cv/jwst_clusters/phangs/"
mdir="ims_2023b/"

apphot_file=[mdir + s for s in [
    'ngc628c_uvis_f275w_err_drc_sci.Mjysr.fits',
    'ngc628c_uvis_f336w_err_drc_sci.Mjysr.fits',
     'ngc628c_acs_f435w_err_drc_sci.Mjysr.fits',
     'ngc628c_acs_f555w_err_drc_sci.Mjysr.fits',
     'ngc628c_acs_f814w_err_drc_sci.Mjysr.fits',
    'ngc0628_nircam_lv3_f200w_i2d_align.fits',
    'ngc0628_nircam_lv3_f300m_i2d_align.fits',
    'ngc0628_nircam_lv3_f335m_i2d_align.fits',
    'ngc0628_nircam_lv3_f360m_i2d_align.fits',
    'ngc0628_miri_lv3_f770w_i2d_align.fits',
    'ngc0628_miri_lv3_f1000w_i2d_align.fits',
    'ngc0628_miri_lv3_f1130w_i2d_align.fits',
    'ngc0628_miri_lv3_f2100w_i2d_align.fits'
]]




#==============================================================
# set up apphot - load ims 

try:
    print(len(apphot_ims))
except:
    print("loading images for aperture photometry")
    apphot_ims=loadims(apphot_file)
apphot_wave=pl.array(apphot_wave)


nimage=len(apphot_ims)
nsrc=len(ra)

app_f   =pl.zeros([nimage,nsrc])
app_df  =pl.zeros([nimage,nsrc])
app_raw =pl.zeros([nimage,nsrc])
app_flag=pl.zeros([nimage,nsrc])
app_bg  =pl.zeros([nimage,nsrc])

# "wiggle" 
wapp_f   =pl.zeros([nimage,nsrc])
wapp_df  =pl.zeros([nimage,nsrc])
wapp_raw =pl.zeros([nimage,nsrc])
wapp_flag=pl.zeros([nimage,nsrc])
wapp_bg  =pl.zeros([nimage,nsrc])

# jimena
jr_f   =pl.zeros([nimage,nsrc])
jr_df  =pl.zeros([nimage,nsrc])
jr_raw =pl.zeros([nimage,nsrc])
jr_flag=pl.zeros([nimage,nsrc])
jr_bg  =pl.zeros([nimage,nsrc])

# HST cat
cat_f   =pl.zeros([nimage,nsrc])
cat_df  =pl.zeros([nimage,nsrc])


#nsrc=40 # test
do_wiggle=False


for isrc in range(nsrc):
    print("starting",isrc)
    pl.clf() #clear figure
    pl.figure(1,figsize=[8,10])
#    pl.show(block=False)
    r.debug=False
    debugph=r.debug
    r.setbgfact([2,4])
    # larger but thinner bg than 2,2.5
    #r.setbgfact([3,3.2])
    # jimena 21um
    # r.setbgfact([3.40,3.67])
    

    # go through aperture radii
    k=0
    f=[]
    df=[]
    raw=[]
    bg=[]
    wf=[]
    wdf=[]
    wraw=[]
    wbg=[]

    panel=[6,6,0]  # number of subpanels
    while k<nimage:
        l=k
        while l<nimage and apphot_radius[l]==apphot_radius[k]: l=l+1
        r.setcircle([ra[isrc],de[isrc],apphot_radius[k]/3600.])

        
#        print "k=",k,l-1," r=",apphot_radius[k],apphot_radius[l-1]
#        if k>=12:
#            pdb.set_trace()
        panel[2]=2*k+1
        f2,df2,raw2,bg2=phot1(r,apphot_ims[k:l],names=apphot_band[k:l],panel=panel,debug=debugph,offsetfrac=0.0,buffr=2.)
        f=pl.concatenate([f,f2])
        df=pl.concatenate([df,df2])
        bg=pl.concatenate([bg,bg2])
        raw=pl.concatenate([raw,raw2])

        if do_wiggle:
            f2,df2,raw2,bg2=phot1(r,apphot_ims[k:l],names=apphot_band[k:l],panel=panel,debug=debugph,offsetfrac=0.1,buffr=2.,showmask="none",plot=False)
            wf=pl.concatenate([wf,f2])
            wdf=pl.concatenate([wdf,df2])
            wbg=pl.concatenate([wbg,bg2])
            wraw=pl.concatenate([wraw,raw2])


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
    wf =1000*pl.array(wf)
    wdf=1000*pl.array(wdf)
    wraw=1000*pl.array(wraw)


    # plot (also set flag for apphot)
    flag=pl.ones(nimage)
    wflag=pl.ones(nimage)

    z=pl.where(pl.isnan(wf))[0]
    if len(z)>0:
        wflag[z]=0
    if do_wiggle:
        z=pl.where(pl.isnan(wf))[0]
        if len(z)>0:
            wflag[z]=0

    if do_UL:
        z=pl.where((f>raw)*(apphot_wave>0))[0]
        if len(z)>0:
            flag[z]=3
        # set nondet flags 
        z=pl.where((f<=0)*(apphot_wave>0))[0]
        if len(z)>0:
            flag[z]=0        
    
        if do_wiggle:
            z=pl.where((wf>wraw)*(apphot_wave>0))[0]
            if len(z)>0:
                wflag[z]=3
            z=pl.where((wf<=0)*(apphot_wave>0))[0]
            if len(z)>0:
                wflag[z]=0 
       


    pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    # save apphot to big array
    app_f[:,isrc]=f
    app_df[:,isrc]=df
    app_flag[:,isrc]=flag
    app_raw[:,isrc]=raw
    # app_bg[:,isrc]=bg
    app_bg[:,isrc]= raw-f # in mJy units

    if do_wiggle:
        wapp_f[:,isrc]=wf
        wapp_df[:,isrc]=wdf
        wapp_flag[:,isrc]=wflag
        wapp_raw[:,isrc]=wraw
        #wapp_bg[:,isrc]=wbg
        wapp_bg[:,isrc]=wraw-wf

    
    pl.subplot(3,2,6)
    z=np.where(flag==1)[0]
    pl.errorbar(apphot_wave[z],f[z],yerr=df[z],fmt='.',label="apphot")
    #pl.plot(apphot_wave[z],raw[z],'b.')
    z=np.where(flag==3)[0]
    pl.plot(apphot_wave[z],f[z],'bv')

    if do_wiggle:
        z=np.where(wflag==1)[0]
        pl.errorbar(apphot_wave[z]*1.05,wf[z],yerr=wdf[z],fmt='m.',label="wiggle")
        #pl.plot(apphot_wave[z],oraw[z],'m.')
        z=np.where(wflag==3)[0]
        pl.plot(apphot_wave[z]*1.05,wf[z],'mv')

    pl.xscale("log")
    pl.yscale("log")

    jf =[jtable[isrc][cc] for cc in jfcols]
    djf=[jtable[isrc][cc] for cc in jdfcols]

    # jimena apperture correction
    jf =jf /j_apcorr 
    djf=djf/j_apcorr

    z=np.where(jf>0)[0]    
    pl.errorbar(apphot_wave[z],jf[z],yerr=djf[z],fmt='x',color='orange',label="Jimena")
    z=np.where(jf<0)[0]
    if len(z)>0:
        pl.plot(apphot_wave[z],np.absolute(jf[z]),'v',color='orange')

    jr_f[:,isrc]=jf
    jr_df[:,isrc]=djf
    jr_raw[:,isrc]= [jtable[isrc][cc] for cc in jfbcols]
    jr_bg[:,isrc] = [jtable[isrc][cc] for cc in jfbgcols]


    z=np.where(sourcename[isrc]==hsttable['ID_PHANGS_CLUSTERS'])[0]
    if len(z)>1:
        print("source name error")
    if len(z)>0:
        hf=[hsttable[z[0]][cc] for cc in hstfcols]
        hdf=[hsttable[z[0]][cc] for cc in hstdfcols]
        
        # jimena aperture correction?
        hf=hf/phangs_apcorr[0:5]
        hdf=np.absolute(hdf)/phangs_apcorr[0:5]

        pl.errorbar(apphot_wave[0:5],hf,yerr=hdf,fmt='o',mfc='none',color='green',label="HST")
        
        cat_f[0:5,isrc]=hf
        cat_df[0:5,isrc]=hdf


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
              "filter":filter_name,
              "jr_f":jr_f,
              "jr_raw":jr_raw,
              "jr_df":jr_df,
              "jr_bg":jr_bg,
              "cat_f":cat_f,
              "cat_df":cat_df,
              }, open( outdir+"/"+regname+".pkl", "wb" ) )


    
    


