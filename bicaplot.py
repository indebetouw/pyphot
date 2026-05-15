from astropy.table import Table
import matplotlib.pyplot as pl
pl.ion()

typ="CN"
typ="NC"
t=Table.read("bica."+typ+".dists.apphot.csv")

fluxes=[x for x in t.colnames if 'f_d1'==x[0:2]]

for i in range(len(fluxes)):
    pl.plot(t["raw_"+fluxes[i][2:]],t[fluxes[i]],'.',label=fluxes[i][2:])

pl.xscale("log")
pl.yscale("log")
pl.legend(loc="best",prop={"size":8})

radpc=np.sqrt(t['rmaj'].data*t['rmin'].data)/4

pl.clf()
for i in range(len(fluxes)):
    pl.plot(radpc,t[fluxes[i]],'.',label=fluxes[i][2:])

pl.yscale("log")
pl.xscale("log")
pl.xlabel("radius [pc]")
pl.ylabel("flux [mJy @ 50kpc, 1e-15 e/s/cm2 for Ha]")
pl.legend(loc="best",prop={"size":8})


pl.clf()
for i in range(len(fluxes)):
    pl.plot(radpc,t[fluxes[i]]/np.pi/radpc**2,'.',label=fluxes[i][2:])

pl.yscale("log")
pl.xscale("log")
pl.xlabel("radius [pc]")
pl.ylabel("ave brightness [mJy/pc2, 1e-15 e/s/cm2/pc2 for Ha]")
pl.legend(loc="best",prop={"size":8})
pl.savefig("bica."+typ+".sbright.png")


# PHANGS aperture is 4 pix * 0.04" =

dists=[5,10,20]
yy=pl.ylim()
for i in range(len(dists)):
    rpc=dists[i]*0.16*1e6/206265
    pl.plot([rpc,rpc],[yy[1]*0.1,yy[1]],'k',alpha=0.5)
    pl.text(rpc,yy[1],"%iMpc"%dists[i],color='k',verticalalignment="top")
    if i==0:
        pl.text(rpc,yy[1],"HST ap@",color='k',horizontalalignment="right",verticalalignment="top")

pl.ylim(yy)
pl.savefig("bica."+typ+".sbright.png")



pl.clf()
pl.plot(radpc,t['f_d1_M24']/t['f_d1_I8.0']  ,'.')
pl.yscale("log")
pl.xscale("log")
pl.xlabel("radius [pc]")
pl.ylabel("Fnu(24)/Fnu(8)")
pl.xlim(1,30)
pl.plot([1.2,1.2],[1.8,11.2],label="Stephens")
pl.legend(loc="best",prop={"size":8})
pl.savefig("bica."+typ+".d1.sbright.png")


pl.clf()
for d in [1,5,10,20]:
    z=np.where(radpc<4)[0]
    pl.hist(np.log10(t['f_'+("d%i"%d)+'_M24']/t['f_'+("d%i"%d)+'_I8.0'])[z],bins=5,histtype="step",label="d=%iMpc"%d)
    #pl.plot(radpc,t['f_'+("d%i"%d)+'_M24']/t['f_'+("d%i"%d)+'_I8.0']  ,'.')
#pl.yscale("log")
#pl.xscale("log")
#pl.xlabel("radius [pc]")
pl.xlabel("log[ Fnu(24)/Fnu(8) ]")
#pl.xlim(1,30)
#pl.plot([1.2,1.2],[1.8,11.2],label="Stephens")
pl.legend(loc="best",prop={"size":8})
pl.savefig("bica."+typ+".brighthist.png")
