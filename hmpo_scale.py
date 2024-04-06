# scale everything to full-region and 10kpc

import numpy as np
from astropy.table import Table
import pandas,pdb
import matplotlib.pyplot as pl
pl.ion()

pd=pandas.read_excel("hmpo24.xlsx") #,dtype=dtype_in)
hmpo=Table.from_pandas(pd)
hmpo.remove_column(hmpo.colnames[0])



swsfilts=[
    'irac1','irac2',
    'jwst.nircam.f300m',
    'jwst.nircam.f335m',
    'jwst.nircam.f360m',
    'jwst.nircam.f444w',
    'jwst.miri.f770w',
    'jwst.miri.f1000w',
    'jwst.miri.f1130w',
    'jwst.miri.f2100w']

akfilts=[
    'irac1','irac2',
    'jwst.nircam.f300m',
    'jwst.nircam.f335m',
    'jwst.nircam.f360m',
    'jwst.nircam.f444w']


nufnu=True

lwsfilts=['pacs1','pacs2','pacs3']
pfilts=np.array(['pg','pr','pi','pz','py','aJ','aH','aK','I1'])

outfilts=['pg','pr','pi','pz','py','aJ','aH','aK',
          'f300m','f335m','f360m','f444w','f770w','f1000w','f1130w','f2100w']
outwaves=[0.481,0.617,0.752,0.866,0.962, 1.235,1.663,2.159,
          3.,3.35,3.6,4.44,7.7,10.,11.3,21.0]
# add all lws filters - we may add another around 30um or so later.
outfilts=np.concatenate([outfilts,lwsfilts])
outwaves=np.concatenate([outwaves,[65,100,160]])


# preferentially use near distances 
hmpo_distance=hmpo['dist'].filled(0).data
hmpo_lum=hmpo['lum_e4_near'].filled(0).data


copycolnames=['viziername','galacticname','RA','Dec','morpho-','log_Q0']

# first scale everything that has "src":
z=np.where(hmpo['src_I1'].data>0)[0]

pl.clf()
pl.gcf().set_size_inches(7,5)
pl.subplots_adjust(left=0.13,bottom=0.1)
for iz in z:

    if iz==z[0]:
        hmpo_scl=Table({"scaled":[True],"lum_e4":[hmpo_lum[iz]]})
        for cc in copycolnames:
            hmpo_scl.add_column(hmpo[iz][cc],name=cc)
        for outfilt in outfilts:
            hmpo_scl.add_column(0.,name=outfilt)
            hmpo_scl.add_column(0.,name="d_"+outfilt)
    else:
        hmpo_scl.add_row({"scaled":True,"lum_e4":hmpo_lum[iz]})
        for cc in copycolnames:
            hmpo_scl[-1][cc]=hmpo[iz][cc]
        
            
    distscale=(hmpo_distance[iz]/10)**2 # to 10 kpc
    
    for pfilt in pfilts[:-1]: # don't add irac1
        hmpo_scl[-1][pfilt]     =hmpo[iz]["src_"+pfilt]  * distscale
        hmpo_scl[-1]["d_"+pfilt]=hmpo[iz]["d_src_"+pfilt]* distscale


    if hmpo['irac1_sws'][iz]>0:
        sfactor=hmpo['src_I1'][iz]*distscale /hmpo['irac1_sws'].filled(0).data[iz]
    else:
        sfactor=0

    if hmpo['irac1_akari'][iz]>0:
        afactor=hmpo['src_I1'][iz]*distscale /hmpo['irac1_akari'][iz].data
    else:
        afactor=0

    if afactor<=0 and sfactor<=0:
        print("no spectroscopy found for source "+hmpo[iz]['viziername'])
        pdb.set_trace()

    for ofilt in ['f300m','f335m','f360m','f444w']:
        if sfactor>0:
            if afactor>0:
                hmpo_scl[-1][ofilt]=0.5*(hmpo[iz][ofilt+"_akari"]*afactor+
                                         hmpo[iz][ofilt+"_sws"]*sfactor)

                hmpo_scl[-1]["d_"+ofilt]= \
                    np.sqrt( (hmpo[iz]["d_"+ofilt+"_akari"]*afactor)**2 + \
                             (hmpo[iz]["d_"+ofilt+"_sws"]*sfactor)**2 + \
                             (hmpo[iz][ofilt+"_akari"]*afactor - \
                              hmpo[iz][ofilt+"_sws"]*sfactor)**2/2 ) /np.sqrt(3)
            else:
                hmpo_scl[-1][ofilt]=hmpo[iz][ofilt+"_sws"]*sfactor
                hmpo_scl[-1]["d_"+ofilt]= hmpo[iz]["d_"+ofilt+"_sws"]*sfactor
        else:
            hmpo_scl[-1][ofilt]=hmpo[iz][ofilt+"_akari"]*afactor
            hmpo_scl[-1]["d_"+ofilt]= hmpo[iz]["d_"+ofilt+"_akari"]*afactor
            
                                         
        
    for ofilt in ['f770w','f1000w','f1130w','f2100w']:
        hmpo_scl[-1][ofilt]=hmpo[iz][ofilt+"_sws"]*sfactor
        hmpo_scl[-1]["d_"+ofilt]= hmpo[iz]["d_"+ofilt+"_sws"]*sfactor

    for ofilt in lwsfilts:
        hmpo_scl[-1][ofilt]=hmpo[iz][ofilt+"_lws"]*sfactor
        hmpo_scl[-1]["d_"+ofilt]= hmpo[iz]["d_"+ofilt+"_lws"]*sfactor


    outf=np.zeros(len(outfilts))
    outdf=np.zeros(len(outfilts))
    for k,ofilt in enumerate(outfilts):
        if nufnu:
            outf[k] =hmpo_scl[-1][ofilt]     *c/outwaves[k]*1e-23
            outdf[k]=hmpo_scl[-1]["d_"+ofilt]*c/outwaves[k]*1e-23
        else:
            outf[k] =hmpo_scl[-1][ofilt]
            outdf[k]=hmpo_scl[-1]["d_"+ofilt]
                  

    zz=np.where(outf>0)[0]
    myplot,=pl.plot(outwaves[zz],outf[zz],label=hmpo_scl[-1]['galacticname'])
    pl.errorbar(outwaves[zz],outf[zz],yerr=outdf[zz],fmt='.',color=myplot.get_color())
        

pl.xscale("log")
pl.yscale("log")
if nufnu:
    pl.ylabel("nu*Fnu [e/s/cm2] @10kpc")
else:
    pl.ylabel("Fnu [Jy] @10kpc")
pl.xlabel("wave [um]")
pl.legend(loc="best",prop={"size":8})

hmpo_scl.write("hmpo_scaled.tbl",format="ascii",overwrite=True)

