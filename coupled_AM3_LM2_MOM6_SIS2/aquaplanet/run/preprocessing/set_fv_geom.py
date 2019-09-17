def set_fv_geom(mlat,nlon):
# Sets the geometric factors used in FV core for the lat-lon grid
# This is based on fv_pack.f90 in the AM2 source code, initially
# edited by Andrew Shao (andrew.shao@noaa.gov)
   import numpy as np

   pi = 4.*np.arctan(1.)
   # Set aliases for trignometric functions
   lonb  = np.zeros(nlon+1)
   rlonb = np.zeros(nlon+1)
   latb  = np.zeros(mlat+1)

   for i in range(1,nlon+2):
      lonb[i-1] = (i-1)*360./(nlon)

   for j in range( 2, mlat+1 ):
#      latb[j-1] = -90. + (((j-1)-0.5) * 180.) / (mlat-1)
      latb[j-1] = -90. + (((j-1)-0.5) * 180.) / (mlat-1)

   latb[0] = -90
   latb[-1] = 90

# To be consistent with the coupler and the atmospheric grid, we must first convert the lat/lon to radians
# and then back to degrees using very slightly different conversion factors
   latb = latb*(pi/180.)
   latb = latb/(np.arctan(1.)/45.)
# Refine latb and lonb so that they resemble a grid with twice the resolution
   latb_super = np.zeros(2*len(latb)-1)
   lonb_super = np.zeros(2*len(lonb)-1)
   latb_super[::2] = latb
   lonb_super[::2] = lonb
   latb_super[1::2] = 0.5*(latb[0:-1]+latb[1:])
   lonb_super[1::2] = 0.5*(lonb[0:-1]+lonb[1:])

   longrid, latgrid = np.meshgrid(lonb_super, latb_super)
   return latgrid, longrid
