import matplotlib.pylab as pl
import astropy.io.fits as pyfits
import pdb

### TODO replace hextract with spectral-cube.subcube
### https://github.com/radio-astro-tools/spectral-cube/
### consider switching to that class as the general Image wrapper here

### TODO rewrite to be more OOP - have a base class with virtuals and then 
### inherit a Circle, Polygon, whatever

# TODO use ginsberg fits_utils instead of this:?
#def hextract(hin,x0,x1,y0,y1,outfile=None):
def hextract(hin,crds,outfile=None):
    """
    works in image coordinates; takes HDU, returns HDU
    """
    x0=crds[0]
    y0=crds[1]
    x1=crds[2]
    y1=crds[3]
    imshape=hin.shape
    # check and sanitize inputs
    if len(imshape)>2: raise Exception("can't extract subim from >2D image")
    if x0<0 or y0<0:
        x00=x0
        y00=y0
        if x0<0: x0=0
        if y0<0: y0=0
        print "blc (%f,%f) is out of bounds, correcting to (%f,%f)" % (x00,y00,x0,y0)
    if x1>(imshape[1]-1) or y1>(imshape[0]-1):
        x10=x1
        y10=y1
        if x1>(imshape[1]-1): x1=imshape[1]-1
        if y1>(imshape[0]-1): y1=imshape[0]-1
        print "trc (%f,%f) is out of bounds, correcting to (%f,%f)" % (x10,y10,x1,y1)
        
    # translate cdelt to center
    from astropy import wcs
    w=wcs.WCS(hin.header)
    c=[0.5*(x0+x1),0.5*(y0+y1)]
    # world coords of ctr of the cutout
    ctr=w.wcs_pix2world([c[0:2]],1)[0]

    hdrout=hin.header.copy(strip=True)
    # this was the old, just translate crpix:
#    hdrout['crpix1']=hdrout['crpix1']-x0
#    hdrout['crpix2']=hdrout['crpix2']-y0
    # better:
    hdrout['crpix1']=0.5*(x1-x0)
    hdrout['crpix2']=0.5*(y1-y0)
    hdrout['crval1']=ctr[0]
    hdrout['crval2']=ctr[1]
    # check
    w=wcs.WCS(hdrout)
#    print ctr,w.wcs_pix2world([[0.5*(x1-x0),0.5*(y1-y0)]],1)[0]
    
    datout=hin.data[y0:y1,x0:x1].copy()
    hduout=pyfits.PrimaryHDU(data=datout,header=hdrout)
    if outfile: hduout.writeto(outfile)
    return pyfits.HDUList([hduout])


#===========================================================================
# combine cat and app phot, also interactively edit fit phot if edit=True
def photcombine(a_wave,a_f,a_df,a_fl,c_wave,c_f,c_df,c_fl,f_wave,f_dwave,edit=False,preference=None):
    nfit=len(f_wave)
    f_f=pl.zeros(nfit)
    f_df=pl.zeros(nfit)
    f_fl=pl.zeros(nfit)

    for i in range(nfit):
        # TODO can't deal with flag=2,4 
        # are there any detections:
        a_det=pl.where((abs(a_wave-f_wave[i])<(0.5*f_dwave[i]))*(a_fl==1))[0]
        c_det=pl.where((abs(c_wave-f_wave[i])<(0.5*f_dwave[i]))*(c_fl==1))[0]
        # are there any UL:                                    
        a_ul =pl.where((abs(a_wave-f_wave[i])<(0.5*f_dwave[i]))*(a_fl==3))[0]
        c_ul =pl.where((abs(c_wave-f_wave[i])<(0.5*f_dwave[i]))*(c_fl==3))[0]

        # any cat UL?
        if len(c_ul)>0:
            # more than one?
            if len(c_ul)>1:
                d=abs(c_wave[c_ul]-f_wave[i])
                closest_c_ul=c_det[ pl.where(d==d.min())[0] ]
                print "ambiguous catalog upper limits, choosing %f for fitter %f" % (c_wave[closest_c_ul], f_wave[i])
                print "     set=",c_wave[c_ul]
            else:
                closest_c_ul=c_ul
        else:
            closest_c_ul=-1

        # any app UL?
        if len(a_ul)>0:
            # more than one?
            if len(a_ul)>1:
                d=abs(a_wave[a_ul]-f_wave[i])
                closest_a_ul=a_det[ pl.where(d==d.min())[0] ]
                print "ambiguous apphot upper limits, choosing %f for fitter %f" % (a_wave[closest_a_ul], f_wave[i])
                print "     set=",a_wave[a_ul]
            else:
                closest_a_ul=a_ul
        else:
            closest_a_ul=-1

        # any app detections?
        if len(a_det)>0:
            # more than one?
            if len(a_det)>1:
                d=abs(a_wave[a_det]-f_wave[i])
                closest_a_det=a_det[ pl.where(d==d.min())[0] ]
                print "ambiguous apphot photometry, choosing %f for fitter %f" % (a_wave[closest_a_det], f_wave[i])
                print "     set=",a_wave[a_det]
            else:
                closest_a_det=a_det
        else:
            closest_a_det=-1

        # any cat detections?
        if len(c_det)>0:
            # more than one?
            if len(c_det)>1:
                d=abs(c_wave[c_det]-f_wave[i])
                closest_c_det=c_det[ pl.where(d==d.min())[0] ]
                print "ambiguous catalog photometry, choosing %f for fitter %f" % (c_wave[closest_c_det], f_wave[i])
                print "     set=",c_wave[c_det]
            else:
                closest_c_det=c_det
        else:
            closest_c_det=-1

        # combine:
        if preference=="cat":
            # user wants cat, there's cat det, done.
            if closest_c_det>=0:                
                f_f[i]=c_f[closest_c_det]
                f_df[i]=c_df[closest_c_det]
                f_fl[i]=1
                # throw away apphot det silently here
                # TODO check if apphot UL is lower than cat_phot?
            elif closest_c_ul>=0:
                # there's no det, but a cat UL - is there app det?
                if closest_a_det>=0:
                    if a_f[closest_a_det]<=c_f[closest_c_ul]:
                        # there's an appdet below the cat UL:
                        f_f[i]=a_f[closest_a_det]
                        f_df[i]=a_df[closest_a_det]
                        f_fl[i]=1
                    else:
                        # there's an appdet _above_ the cat UL - WTF?
                        print "apphot detection brighter than catalog UL at ",f_wave[i]
                        # assume apphot is wrong
                        f_f[i]=c_f[closest_c_ul]
                        f_df[i]=c_df[closest_c_ul]
                        f_fl[i]=3
                else:
                    # start with that cat UL
                    f_f[i]=c_f[closest_c_ul]
                    f_df[i]=c_df[closest_c_ul]
                    f_fl[i]=3
                    # now if there's also an app UL
                    if closest_a_ul>=0:
                    # and its lower
                        if a_f[closest_a_ul]<=c_f[closest_c_ul]:
                            # use lower app UL instead of cat UL:
                            f_f[i]=a_f[closest_a_ul]
                            f_df[i]=a_df[closest_a_ul]
                            f_fl[i]=3
            else:
                # user wanted cat, but there's no cat.
                if closest_a_det>=0:
                    f_f[i]=a_f[closest_a_det]
                    f_df[i]=a_df[closest_a_det]
                    f_fl[i]=1
                elif closest_a_ul>=0:
                    f_f[i]=a_f[closest_a_ul]
                    f_df[i]=a_df[closest_a_ul]
                    f_fl[i]=3
                # otherwise they get nothing - f_fl stays=0
        
        elif preference=="app":
            # user wants app, there's app det, done.
            if closest_cadet>=0:                
                f_f[i]=a_f[closest_a_det]
                f_df[i]=a_df[closest_a_det]
                f_fl[i]=1
                # throw away catphot det silently here
                # TODO check if catphot UL is lower than appphot?
            elif closest_a_ul>=0:
                # there's no det, but an app UL - is there a cat det?
                if closest_c_det>=0:
                    if c_f[closest_c_det]<=a_f[closest_a_ul]:
                        # there's an catdet below the app UL:
                        f_f[i]=c_f[closest_c_det]
                        f_df[i]=c_df[closest_c_det]
                        f_fl[i]=1
                    else:
                        # there's an catdet _above_ the app UL - WTF?
                        print "catalog detection brighter than appphot UL at ",f_wave[i]
                        # assume apphot is wrong
                        f_f[i]=c_f[closest_c_det]
                        f_df[i]=c_df[closest_c_det]
                        f_fl[i]=1
                else:
                    # start with that app UL
                    f_f[i]=a_f[closest_a_ul]
                    f_df[i]=a_df[closest_a_ul]
                    f_fl[i]=3
                    # now if there's also a cat UL
                    if closest_c_ul>=0:
                    # and its lower
                        if c_f[closest_c_ul]<=a_f[closest_a_ul]:
                            # use lower app UL instead of cat UL:
                            f_f[i]=c_f[closest_c_ul]
                            f_df[i]=c_df[closest_c_ul]
                            f_fl[i]=3
            else:
                # user wanted app, but there's no app.
                if closest_c_det>=0:
                    f_f[i]=c_f[closest_c_det]
                    f_df[i]=c_df[closest_c_det]
                    f_fl[i]=1
                elif closest_c_ul>=0:
                    f_f[i]=c_f[closest_c_ul]
                    f_df[i]=c_df[closest_c_ul]
                    f_fl[i]=3
                # otherwise they get nothing - f_fl stays=0
        
        else: # preference is neither cat nor app:
            # implicit preference for cat but some averaging
            if closest_c_det>=0:
                if closest_a_det>=0:
                    # 2 dets -average
                    f_f[i] =0.5( c_f[closest_c_det]+a_f[closest_a_det] )
                    f_df[i]=max([ c_df[closest_c_det],
                                  a_df[closest_a_det],
                                  abs(c_f[closest_c_det]-a_f[closest_a_det]) ])
                    f_fl[i]=1
                else:
                    # cat det; is there an app UL?
                    if closest_a_ul>=0:
                        if a_f[closest_a_ul]<=c_f[closest_c_det]:
                            print "apphot UL below cat detection at ",f_wave[i]
                            # in case of discrepency, assum cat correct
                        f_f[i]=c_f[closest_c_det]
                        f_df[i]=c_df[closest_c_det]
                        f_fl[i]=1
            elif closest_c_ul>=0:
                # there's a catalog UL, but no det:
                # start by assuming cat right
                f_f[i]=c_f[closest_c_ul]
                f_df[i]=c_df[closest_c_ul]
                f_fl[i]=3                
                if closest_a_det>=0:
                    if a_f[closest_a_det]<=c_f[closest_c_ul]:
                        # apphot det below cat UL- replace with that
                        f_f[i]=a_f[closest_a_det]
                        f_df[i]=a_df[closest_a_det]
                        f_fl[i]=1
                elif closest_a_ul>=0:
                    if a_f[closest_a_ul]<=c_f[closest_c_ul]:
                        # apphot UL below cat UL- replace with that
                        f_f[i]=a_f[closest_a_ul]
                        f_df[i]=a_df[closest_a_ul]
                        f_fl[i]=3


    # next, set uncert minima to 10%
    z=pl.where(f_fl==1)[0]
    for zz in z:
        f_df[zz]=max([f_df[zz],0.1*f_f[zz]])

    # todo check for and set UL confidence levels?
    if edit:
        global whatx,fit_1,fit_3,startpos,endpos,fits_1,xs_1,x1,y1,x3,y3
    # for interactive editing
    #  plot fit phot and prepare to edit it
    z=pl.where(f_fl==1)[0]
    if len(z)>0:
        fit_1=pl.plot(f_wave[z],f_f[z],'r.',markersize=8,label="fitter")[0]
        fits_1=[]
        for j in range(len(z)):
            uncert=f_f[z[j]]+pl.array([-1,1])*f_df[z[j]]
            #if uncert[0]<pl.ylim()[0]: uncert[0]=pl.ylim()[0]
            fits_1.append(pl.plot(f_wave[z[j]]*pl.array([1,1]),uncert,'r')[0])
    else:
        fit_1=None            
    z=pl.where(f_fl==3)[0]
    if len(z)>0:
        fit_3=pl.plot(f_wave[z],f_f[z],'rv')[0]
    else:
        fit_3=None

    ndets=len(fits_1)
    xs_1=pl.zeros(ndets)# x locations of the error bars
    for k in range(ndets):
        xs_1[k]=fits_1[k].get_data()[0][0]
    
    pl.legend(loc=4,prop={'size':8},numpoints=1)

    if edit:
        def click(event):
            if not event.inaxes: return
            global whatx,fit_1,fit_3,startpos,endpos,fits_1,xs_1,x1,y1,x3,y3
            startpos=event.xdata, event.ydata       
        
            # find closest existing pt
            if fit_1==None:
#                print "no fit_1?!"
                x1=[]
                y1=[]
                d1=pl.array([1e10])
            else:
                x1,y1=fit_1.get_data()
                d1=abs(event.xdata-x1)
            if fit_3==None:
#                print "no fit_3?!"
                x3=[]
                y3=[]
                d3=pl.array([1e10])
            else:
                x3,y3=fit_3.get_data()
                d3=abs(event.xdata-x3)
    
            # todo: for deletions, make sure we have all avail wavelength pts
            # i suppose that the flux combination step that creates fit_wave 
            # will do that...
    
#            print "x1=",x1
#            print "x3=",x3

            if len(d1)<=0:
                d1=pl.array([1e10])
            if len(d3)<=0:
                d3=pl.array([1e10])

            if d1.min()<=d3.min():
                whatpoint=pl.where(d1==d1.min())[0][0]
                whatx=x1[whatpoint]
                print "deleting detection %d @ "%whatpoint,whatx
                fit_1.set_data(pl.delete(x1,whatpoint),pl.delete(y1,whatpoint))
                # delete the uncert error line too
#                ds_1=abs(event.xdata-xs_1)
#                k=pl.where(ds_1==ds_1.min())[0][0]
                k=whatpoint
                fits_1[k].remove()
                fits_1=pl.delete(fits_1,k)
                xs_1=pl.delete(xs_1,k)
            else:
                whatpoint=pl.where(d3==d3.min())[0][0]
                whatx=x3[whatpoint]
                print "deleting UL %d @ "%whatpoint,whatx
                x3=pl.delete(x3,whatpoint)
                y3=pl.delete(y3,whatpoint)
                fit_3.set_data(x3,y3)
    
            if event.button==3: #R-click
                x3=pl.append(x3,whatx)
                y3=pl.append(y3,startpos[1])
                if fit_3==None:
                    fit_3=pl.plot(x3,y3,'rv')[0]
                else:
                    fit_3.set_data(x3,y3)
    
            pl.draw()
#            print x3
#            print x1
#            print xs_1
                
        def unclick(event):
            if not event.inaxes: return
            global whatx,fit_1,fit_3,startpos,endpos,fits_1,xs_1,x1,y1,x3,y3
            endpos=event.xdata, event.ydata
            if event.button==1:
                if fit_1:
                    x1,y1=fit_1.get_data()
                    x1=pl.append(x1,whatx)
                    y1=pl.append(y1,0.5*(startpos[1]+endpos[1]))
                    fit_1.set_data(x1,y1)
                else:
                    fit_1=pl.plot(whatx,0.5*(startpos[1]+endpos[1]),'r.')[0]
                    fits_1=[]
                # add this to the list of uncert lines plots
                fits_1=pl.append(fits_1,pl.plot([whatx,whatx],[startpos[1],endpos[1]],'r')[0])
                xs_1=pl.append(xs_1,whatx)
#                print "xs_1 = ",xs_1
                # XXX TODO also set the uncert somewhere                    
            pl.draw()
#            print x3
#            print x1
#            print xs_1
        
        cid0=pl.connect('button_press_event',click)
        cid1=pl.connect('button_release_event',unclick)
    
        print "edit fitter points and then press enter in the terminal"
        x=raw_input()
    
        pl.disconnect(cid0)
        pl.disconnect(cid1)

    if not fit_3==None:
        x3,y3=fit_3.get_data()
        print "upper limits are now:",x3,y3
    
        for j in range(len(x3)):
            d=abs(f_wave-x3[j])
            z=pl.where(d==d.min())[0][0]
            f_f[z]=y3[j]
            f_fl[z]=3
            f_df[z]=0.999 # XXXX
    
    x1,y1=fit_1.get_data()
    print "detections are now:",x1,y1

    for j in range(len(x1)):
        d=abs(f_wave-x1[j])
        z=pl.where(d==d.min())[0][0]
        f_f[z]=y1[j]
        f_fl[z]=1
        f_df[z]= 0.1*f_f[z] # XXXX need real uncert from drawing!

    return f_f,f_df,f_fl
                


                



#===========================================================================
class Region:
    def __init__(self):
        self.bgfact=[2,2.5]
        self.type=None
        self.coords=[]
        self.bg0coords=[] # inner part of annulus
        self.bg1coords=[] # outer part of annulus
        self.mask=None
        self.debug=False

    #-------------------------------------------------------
    def setcircle(self,args):
        """
        args array = [ra,dec,rad in decimal degrees]
        """
        # TODO add radius units?
        self.type="circle"
        self.coords=pl.array(args)
        self.setbgcoords()

    #-------------------------------------------------------
    def setpoly(self,coordarr):
        """
        ra1,de1,ra2,de2... in decimal degrees
        """
        self.type="polygon"
        n=len(coordarr)/2
        c=pl.array(coordarr)
        xi=2*pl.array(range(n))
        yi=xi+1
        self.coords=pl.array([c[xi],c[yi]]).T
        # TODO checks on coordarr - at least check if its even
        # TODO do we need to close the poly to make other things work?
        self.setbgcoords()


# on hold - needs post-processing of coords to get them in the expected format
# for polygon i.e. array of [n,2] xy pairs
#    def setds9(self,str):
#        """
#        parse ds9 region string
#        """
#        # TODO raise exception for something other than circle and poly, 
#        # or leave that to the functions that deal with different types?
#        try:
#            import pyregion
#            r=pyregion.parse(str)
#            self.type=r[0].name
#            self.coords=r[0].coord_list
#            # TODO r.check_imagecoord() better be false - we want to store
#            # coordinates in radec not image coords
#        except:
#            raise Exception("could not import pyregion")
#        self.setbgcoords()

    #-------------------------------------------------------
    def setbgfact(self,bgfact):
        """
        set background scale factor to the 2-element array input parameter
        """
        if len(bgfact)!=2:
            raise Exception("input parameter should be 2-element array")
        self.bgfact=bgfact
      
    #-------------------------------------------------------
    def setbgcoords(self):
        if self.type==None: 
            raise Exception("region type=None - has it been set?")
        if self.type=="circle":
            if len(self.coords)!=3:
                raise Exception("region coords should be ctr_ra, ctr_dec, rad_arcsec - the coord array has unexpected length %d" % len(self.coords))
            self.bg0coords=pl.array(self.coords)
            self.bg1coords=pl.array(self.coords)
            # set larger radii for annulus
            self.bg0coords[2]=self.coords[2]*self.bgfact[0]
            self.bg1coords[2]=self.coords[2]*self.bgfact[1]
        elif self.type=="polygon":
            n=self.coords.shape[1]
            self.coords=pl.array(self.coords)
            ctr=[ self.coords[:,0].mean(), self.coords[:,1].mean() ]
            x=self.coords[:,0]-ctr[0]
            y=self.coords[:,1]-ctr[1]
            r=pl.sqrt(x**2+y**2)
            th=pl.arctan2(y,x)

            ct=pl.cos(th)
            st=pl.sin(th)
            # inner and outer background regions
            b=self.bgfact
            self.bg0coords=pl.array([r*b[0]*ct, r*b[0]*st]).T+ctr
            self.bg1coords=pl.array([r*b[1]*ct, r*b[1]*st]).T+ctr
            
        else: raise Exception("unknown region type %s" % self.type)

    #-------------------------------------------------------
    def plotradec(self):
        if self.type=="polygon":
            pl.plot(self.coords[:,0]   ,self.coords[:,1])
            pl.plot(self.bg0coords[:,0],self.bg0coords[:,1])
            pl.plot(self.bg1coords[:,0],self.bg1coords[:,1])
        elif self.type=="circle":
            n=23
            t=pl.arange(n)*pl.pi*2/n
            r=self.coords[2]
            ct=pl.cos(t)
            st=pl.sin(t)
            cdec=pl.cos(self.coords[1]*pl.pi/180)
            pl.plot(self.coords[0]+r*ct/cdec, self.coords[1]+r*st)
            r=self.bg0coords[2]
            pl.plot(self.bg0coords[0]+r*ct/cdec, self.bg0coords[1]+r*st)
            r=self.bg1coords[2]
            pl.plot(self.bg1coords[0]+r*ct/cdec, self.bg1coords[1]+r*st)

    #-------------------------------------------------------
    def plotimx(self,im):
        if self.type=="polygon":
            for reg in ("ap","bg0","bg1"):
                ci=self.imcoords(im,reg=reg)
                pl.plot(ci[:,0],ci[:,1])
        elif self.type=="circle":
            n=33
            t=pl.arange(n)*pl.pi*2/n            
            ct=pl.cos(t)
            st=pl.sin(t)
            for reg in ("ap","bg0","bg1"):
                ci=self.imcoords(im,reg=reg)
                #print "center of circle= %f %f" % ci[0:2]
                r=ci[2] # in pix
                pl.plot(ci[0]+r*ct, ci[1]+r*st)

            # point north:  XXX TODO make general
            from astropy import wcs
            w=wcs.WCS(im.header)
            # use origin=0 i.e. NOT FITS convention, but pl.imshow sets origin
            # to 0,0 so do that here so we can overplot on pl.imshow axes
            origin=0
            c=self.bg1coords # only works below for circle
            ctr=w.wcs_world2pix([c[0:2]],origin)[0]
            # north
            ctr2=w.wcs_world2pix([c[0:2]+pl.array([c[2],0])],origin)[0]
            pl.plot([ctr[0],ctr2[0]],[ctr[1],ctr2[1]])



    #-------------------------------------------------------
    def imcoords(self,im,reg="ap"):
        """
        given an image, return the region coords in image coordinates (pixels)
        can optionally specify reg=["ap","bg0","bg1"] (default=ap)
        for the aperture, the inner background and the outer background
        """
        if reg=="ap":
            c=self.coords
        elif reg=="bg0":
            c=self.bg0coords
        elif reg=="bg1":
            c=self.bg1coords
        else: raise Exception("unknown input reg=%s" % reg)

        from astropy import wcs
        w=wcs.WCS(im.header)
        # use origin=0 i.e. NOT FITS convention, but pl.imshow sets origin
        # to 0,0 so do that here so we can overplot on pl.imshow axes
        origin=0

        if self.type=="circle":
            ctr=w.wcs_world2pix([c[0:2]],origin)[0]
            # I couldn't find a simple way to convert from arcsec to pix
            # probably because it only is well-defined if the pixels are square
            # and non-distorted.  so assume that for now:
            ctr2=w.wcs_world2pix([c[0:2]+pl.array([0,c[2]])],origin)[0]
            rad=ctr2[1]-ctr[1] 
            return ctr[0],ctr[1],rad # should all be in pix now
        elif self.type=="polygon":
            return w.wcs_world2pix(c,origin)
        else: raise Exception("unknown region type %s" % self.type)
                

    #-------------------------------------------------------
    def imextents(self,im,buffer=1.):
        """
        extents in pixel space - buffer in units of inner ap radius
        """
        ca=self.imcoords(im,reg="ap")
        c1=self.imcoords(im,reg="bg1") # outer bg
        if self.type=="circle":
            dr=ca[2]
            xra=c1[0] + (dr*buffer+ c1[2])*pl.array([-1,1])
            yra=c1[1] + (dr*buffer+ c1[2])*pl.array([-1,1])
        elif self.type=="polygon":                
            dx=minmax(ca[0,:])
            dy=minmax(ca[1,:])
            ctr=[mean(dx),mean(dy)]
            dx=dx[1]-dx[0]
            dy=dy[1]-dy[0]
            dr=0.5*max(dx,dy)
            xra = minmax(c1[0,:]) +dr*buffer*pl.array([-1,1])
            yra = minmax(c1[1,:]) +dr*buffer*pl.array([-1,1])

        xra[0]=pl.floor(xra[0])+1
        yra[0]=pl.floor(yra[0])+1
        xra[1]=pl.ceil(xra[1]) +1
        yra[1]=pl.ceil(yra[1]) +1
            
        if xra[0]<0: xra[0]=0
        if yra[0]<0: yra[0]=0
        s=im.shape-pl.array([1,1]) # remember, transposed y,x
        if xra[1]>s[1]: xra[1]=s[1]
        if yra[1]>s[0]: yra[1]=s[0]

        return int(xra[0]),int(yra[0]),int(xra[1]),int(yra[1])
        

            


    #-------------------------------------------------------
    def setmask(self,im):
        """
        input an image (for now an HDU) and set self.mask to 
        an array the size of the image with the phot region =1
          and expanded background annulus =2
        for now we also create a mask the size of the image, so I recommend
          to extract a subimage and call this method with that input
        this method well trim the polyon to fit in the image
        """
        imshape=im.shape
        mask=pl.zeros(imshape)
        
        if self.type=="circle":
            x,y,r=self.imcoords(im) 
            x0=int(x); y0=int(y)
            dx=x-x0; dy=y-y0
            # grr pixel centers again - is this right?
#            dx=dx-0.5; dy=dy-0.5
            
            bg0_r=self.imcoords(im,reg="bg0")[2] #-0.2 # fudge
            bg1_r=self.imcoords(im,reg="bg1")[2] #+0.2 # fudge
            bg1_r0=int(pl.ceil(bg1_r))
            r2=r**2
            bg0_r2=bg0_r**2
            bg1_r2=bg1_r**2
            for i in     pl.array(range(2*bg1_r0+1))-bg1_r0:
                for j in pl.array(range(2*bg1_r0+1))-bg1_r0:
                    if y0+j>=0 and x0+i>=0 and y0+j<(imshape[0]-1) and x0+i<(imshape[1]-1):
                        d2=(1.*i-dx)**2+(1.*j-dy)**2
                        # d2 = (i-x)**2 + (j-y)**2 -> (i-x0-(x-x0))**2 + ...
                        if d2<=r2:
                            mask[y0+j,x0+i]=1 # remember indices inverted
                        if d2>=bg0_r2 and d2<=bg1_r2:
                            mask[y0+j,x0+i]=2 # remember indices inverted
#                        if x0+i==6: 
#                           print i,j,x0+i,y0+j,dx,dy,d2,bg0_r2,bg1_r2
                        
        elif self.type=="polygon":
            # turn annulus back into mask, will trim at edges of image
            from matplotlib.path import Path
            from matplotlib import __version__ as mpver
            v=mpver.split('.')
            if v[0]<1:
                raise Exception("need matplotlib >=1.3.1, or tell remy to add fallback nxutils option for Path.contains_points")
            elif v[1]<3:
                raise Exception("need matplotlib >=1.3.1, or tell remy to add fallback nxutils option for Path.contains_points")
            elif v[2]<1:
                raise Exception("need matplotlib >=1.3.1, or tell remy to add fallback nxutils option for Path.contains_points")

            # Create vertex coordinates for each grid cell
            x, y = pl.meshgrid(pl.arange(imshape[1]), pl.arange(imshape[0]))
            x, y = x.flatten(), y.flatten()
            points = pl.vstack((x,y)).T
            mask1 = Path(self.imcoords(im,reg="bg1")).contains_points(points)
            mask1 = mask1.reshape((imshape[0],imshape[1]))
            mask0 = Path(self.imcoords(im,reg="bg0")).contains_points(points)
            #,radius=1)        
            mask0 = mask0.reshape((imshape[0],imshape[1]))
            mask = Path(self.imcoords(im,reg="ap")).contains_points(points)
            mask = mask.reshape((imshape[0],imshape[1]))
                        
            mask = mask + (1*mask1-1*mask0)*2 
        else: raise Exception("unknown region type %s" % self.type)
        self.mask=mask
        return mask

    #-------------------------------------------------------
    def phot(self,im,showmask=True):
        # TODO if we switch to astropy.photometry then we can have that 
        # do the work with subpixels properly, but for now they don't 
        # do rms of the bg correctly so we can't use their stuff yet.

        mask=self.setmask(im)

        if showmask:
            cmap1=pl.matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',["black","blue"],2)
            cmap1._init()
            cmap1._lut[:,-1] = pl.array([0,0.5,0,0,0])
            pl.imshow(mask>0,origin="bottom",interpolation="nearest",cmap=cmap1)

        from scipy import ndimage
        from scipy.ndimage import measurements as m
        nin=len(pl.where(mask==1)[0])
        nout=len(pl.where(mask==2)[0])

        floor=pl.nanmin(im.data)
        if floor<0: floor=0
        raw=m.sum(im.data,mask,1)-floor*nin

        #bg=m.mean(im.data,mask,2)
        #bgsig=m.standard_deviation(im.data,mask,2)

#        from astropy.stats import sigma_clip
#        clipped = sigma_clip(im.data,sig=3,iters=2)
#        # http://astropy.readthedocs.org/en/latest/api/astropy.stats.sigma_clip.html#astropy.stats.sigma_clip
#        # TODO what we really want is to sigma-clip only the BG array/mask
#        # because including the source will probably just be domimated by the
#        # source...
#        bg   =m.mean(              clipped,mask,2)-floor
#        bgsig=m.standard_deviation(clipped,mask,2)

        # sigma_clip doesn't handle nans
        from scipy import stats
        def mymode(x):
            return stats.mode(x,axis=None)[0][0]

#        pdb.set_trace()
#        xx=stats.mode(im.data,axis=None)
#        print xx

        bg = ndimage.labeled_comprehension(im.data,mask,2,mymode,"float",0)-floor
#        bg = ndimage.labeled_comprehension(im.data,mask,2,pl.mean,"float",0)
        bgsig=m.standard_deviation(im.data,mask,2)
        

        # assume uncert dominated by BG level.
        # TODO add sqrt(cts in source) Poisson - need gain or explicit err/pix
        uncert = bgsig*nin/pl.sqrt(nout)
        

        results=raw, bg, raw-bg*nin, uncert

        
        f=self.photfactor(im)
        if f:
            if self.debug: print "phot factor = ",f
            results=pl.array(results)*f

        if self.debug:
#            print "max=", m.maximum(im.data,mask,1), m.maximum(im.data,mask,2)
#            print "nin,nout=",nin,nout 
            print "raw, bg, bgsubbed, uncert=", results
            pdb.set_trace()

        return results


    #-------------------------------------------------------
    def pixarea(self,im):
        """
        area of pixel in deg2; this doesn't belong in the Region ...       
        """
        from astropy import wcs
        mywcs=wcs.WCS(im.header)
        # TODO check that get_pc() returns unity matrix in CDELT header
        # this works for a CD matrix header AFAICT
        cd=mywcs.wcs.get_pc()*mywcs.wcs.get_cdelt()
        return -pl.det(cd)

    #-------------------------------------------------------
    def photfactor(self,im):
        """
        attempt to determine the multiplicative factor from pixel values to Jy
        this doesn't belong in the Region ...
        """
        h=im.header
        bunit=h.get("bunit")
        bunit=bunit.strip().upper()
        if bunit==None:
            bunit=h.get("qtty____")
            if bunit==None:
                return None
                # raise Exception("can't find BUNIT or QTTY____ in your image")
        # pixel area in sr
        pixsr=self.pixarea(im) * (pl.pi/180)**2

        if bunit=="MJY/SR":
            return 1.e6 * pixsr
        elif bunit=="JY/PIXEL":
            return 1.
        elif bunit=="JY/BEAM":
            # SPIRE only for the moment
            desc=h.get("desc")
            if desc=="PSW map":
                return 112.197 * 1.e6 * pixsr
            elif desc=="PMW map":
                return 61.415 * 1.e6 * pixsr
            elif desc=="PLW map":
                return 24.336 * 1.e6 * pixsr
            else:
                raise Exception("can't figure out beam area for Jy/beam bunit")
        elif bunit=="2MASS":
            return 7.e-06*4 # for 4x pixelization from 1s-2s 
        elif bunit=="IRSF":
            return 3.98e5 * pixsr
        else:
            #return None
            raise Exception("don't understand your bunit")


    #-------------------------------------------------------        
    def selftest(self,im,fig=True):        
        """
        assuming the region is set, plot transformations on the image
        """
        if fig:
            pl.clf()
            pl.ion()
            pl.subplot(222)
            self.plotradec()
            ax=pl.gca()
            ax.set_aspect("equal")
            pl.title("ra/dec")
        
            pl.subplot(221)
            self.plotimx(im)
            xlim=pl.xlim()
            ylim=pl.ylim()
            ax=pl.gca()
            ax.set_aspect("equal")
            pl.title("image coords")
        
            pl.subplot(223)
            pl.imshow(r.setmask(im),origin="bottom",interpolation="nearest")
            pl.xlim(xlim)
            pl.ylim(ylim)
            pl.title("image mask")
            pl.colorbar()
        
            pl.subplot(224)
            subim=hextract(im,[int(xlim[0]),int(ylim[0]),int(xlim[1]),int(ylim[1])])
            im=subim[0]
            s=im.shape
            # don't mess with "extent" it screws up where the pixel centers are
            pl.imshow(im.data,origin="bottom",interpolation="nearest")
            pl.xlim([0,s[1]])
            pl.ylim([0,s[0]])
            self.plotimx(im)

        self.phot(im)
        if fig:
            pl.subplot(222)
            pl.imshow(self.mask,origin="bottom",interpolation="nearest")
            self.plotimx(im)
            pl.xlim([0,s[1]])
            pl.ylim([0,s[0]])
            pl.title("subimage mask")        
            pl.show()
    
    #-------------------------------------------------------
    def selftest_dor(self,fig=True):
        """
        graphic test of transformations; needs a couple of test files
        """
        import os.path
        regfile="30dor.pacs100.test.reg"
        datfile="30dor.pacs160.fits"
        if not os.path.exists(regfile): raise Exception("Need test file "+regfile)
        if not os.path.exists(datfile): raise Exception("Need test file "+datfile)
        
        # need an image
        f=pyfits.open(datfile)
        # if hdu 0 has no data, go to hext hdu - if wcs is in hdu0 and data
        # in hdu1 then we'll be in trouble.
        i=0
        while len(f[i].data)<1: i=i+1

        import pyregion
        pyreg=pyregion.open(regfile)

        print ("a polygon region")
        self.setpoly(pyreg[1].coord_list)
        if fig: pl.figure()
        self.selftest(f[i],fig=fig)

        print "press enter"
        x=raw_input()

        print ("a circle region")
        self.setcircle(pyreg[0].coord_list)
        if fig: pl.figure()
        self.selftest(f[i],fig=fig)




#===========================================================================
# TODO some kind of ImList class that has band and wave info?

def loadims(imfiles):
    import os.path
    imlist=[]
    for f in imfiles:
        print f
        if not os.path.exists(f): raise Exception("Need file "+f)
        im=pyfits.open(f)
        # if hdu 0 has no data, go to hext hdu - if wcs is in hdu0 and data
        # in hdu1 then we'll be in trouble.
        i=0
        while len(im[i].data)<1: i=i+1
#        z=pl.where(pl.isnan(im[i].data))
#        pdb.set_trace()
#        im[i].data[z[0]]=0.
        imlist.append(im[i])
    return imlist
        
def phot1(r,imlist,plot=True,names=None,panel=None,debug=None,showmask="both"):
    """
    give me a region and an imlist and tell me whether to plot cutouts
    """
    nim=len(imlist)
    # TODO alg for how many subpanels.
    if panel==None: panel=[5,4,1] # for vertical page
    f=[]
    df=[]
    raw=[]
    for j in range(nim):  
        im=imlist[j]
        # if plotting, need to make image cutouts; even if not, this tells
        # us if the source is off the edge
        xtents=r.imextents(imlist[j])
        minsize=5
        if (xtents[3]-xtents[1])<minsize or (xtents[2]-xtents[0])<minsize:
            raw0,f0,df0=0,0,0
            print "phot region too small - %f,%f less than %d pixels" % ((xtents[3]-xtents[1]),(xtents[2]-xtents[0]),minsize)
            print xtents,r.imcoords(im,reg="bg1")
        else:
            if plot: 
                pl.subplot(panel[0],panel[1],panel[2])
                im=hextract(imlist[j],xtents)[0]
                ## ROTATE
                rotate=True
                if rotate:
                    from astropy import wcs
                    w=wcs.WCS(im.header)
                    from scipy.ndimage.interpolation import rotate
                    if w.wcs.has_crota():
                        t0=w.wcs.crota[1]
                    else:
                        t0=pl.arctan2(w.wcs.cd[0,1],-w.wcs.cd[0,0])*180/pl.pi                
                    theta=-1*t0
                    im.data=rotate(im.data,theta,reshape=False)
                    ct=pl.cos(pl.pi*theta/180)
                    st=pl.sin(pl.pi*theta/180)
                    if w.wcs.has_crota():
                        w.wcs.crota[1]=w.wcs.crota[1]+theta
                        im.header['CROTA2']=w.wcs.crota[1]
                        print "rotating crota by "+str(theta)
                    else:
                        w.wcs.cd=pl.matrix(w.wcs.cd)*pl.matrix([[ct,-st],[st,ct]])
                        im.header['CD1_1']=w.wcs.cd[0,0]
                        im.header['CD1_2']=w.wcs.cd[0,1]
                        im.header['CD2_1']=w.wcs.cd[1,0]
                        im.header['CD2_2']=w.wcs.cd[1,1]
                        print "rotating cd    by "+str(theta)
                    
                    #pdb.set_trace()

                # ugh need minmax of aperture region...
                # estimate as inner 1/2 for now
                s=im.shape
                z=im.data[int(s[0]*0.25):int(s[0]*0.75),int(s[1]*0.25):int(s[1]*0.75)]
                if len(z[0])<=0:
                    print z

                z=pl.where(pl.isnan(im.data))
                if len(z[0])>0:
                    z=pl.where(pl.isnan(im.data)==False)
                    std=im.data[z[0],z[1]].std()
                else:
                    std=im.data.std()
                rg=pl.median(im.data)+pl.array([-0.5,5])*std
                # marta wants them less saturated
                rg[1]=pl.nanmax(im.data)

                if rg[0]<0: rg[0]=0
                if rg[1]<=rg[0]:
                     rg=[pl.nanmin(z),pl.nanmax(z)]
                if showmask==False or showmask=="both": # show the jet one                    
                    pl.imshow(im.data,origin="bottom",interpolation="nearest",vmin=rg[0],vmax=rg[1])
                elif showmask==True: # only show the mask
                    pl.imshow(im.data,origin="bottom",interpolation="nearest",vmin=rg[0],vmax=rg[1],cmap="YlGn")
                ax=pl.gca()
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
                if names:
                    pl.text(0.05, 0.99, names[j], horizontalalignment='left',
                            verticalalignment='top',
                            transform=ax.transAxes,
                            bbox=dict(facecolor='white', alpha=0.5))
                pl.xlim([-0.5,s[1]-0.5])
                pl.ylim([-0.5,s[0]-0.5])
                if showmask==False:
                    r.plotimx(im) # overplot the apertures

                if showmask=="both":
                    panel[2]=panel[2]+1
                    pl.subplot(panel[0],panel[1],panel[2])
                    rg=pl.median(im.data)+pl.array([-0.5,5])*std
                    if rg[0]<0: rg[0]=0
                    if rg[1]<=rg[0]:
                        rg=[pl.nanmin(z),pl.nanmax(z)]
                    pl.imshow(im.data,origin="bottom",interpolation="nearest",vmin=rg[0],vmax=rg[1],cmap="YlGn")
                    ax=pl.gca()
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                    if names:
                        pl.text(0.05, 0.99, names[j], horizontalalignment='left',
                                verticalalignment='top',
                                transform=ax.transAxes,
                                bbox=dict(facecolor='white', alpha=0.5))
                    pl.xlim([-0.5,s[1]-0.5])
                    pl.ylim([-0.5,s[0]-0.5])

            if debug:
#                r.debug=True
                if names:
                    print names[j]
            
            raw0,bg,f0,df0=r.phot(im)
        raw.append(raw0)
        f.append(f0)
        df.append(df0)
        panel[2]=panel[2]+1
    if plot:
        pl.subplots_adjust(wspace=0.02,hspace=0.02,left=0.1,right=0.97,top=0.95,bottom=0.05)    
    return f,df,raw
            



#===========================================================================

# test region transformations (radec, im space, masking)
#r.selftest_dor(fig=False)




