REQUIRES:
pyregion
https://pyregion.readthedocs.io/en/latest/installation.html


photometry on lists of images,
either from a list of source positions and circular apertures
or from a (list of) ds9 regions containing polygons

execfile("phot.py")
r=Region()
r.selftest_dor()

and read phot.py::Region::selftest_dor as an example of what to do.