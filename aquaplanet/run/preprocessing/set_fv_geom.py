import numpy as np
def set_fv_geom(mlat,nlon):
# Sets the geometric factors used in FV core for the lat-lon grid
# This is based on fv_pack.f90 in the AM2 source code, initially
# converted into python by tools found at http://www.numericalexpert.com/blog/fotran90_to_python/
# and manually edited by Andrew Shao (andrew.shao@noaa.gov)

   # Set aliases for trignometric functions
   sin = np.sin
   cos = np.cos
   atan = np.arctan
   lon   = np.zeros(nlon)
   lonb  = np.zeros(nlon+1)
   rlonb = np.zeros(nlon+1)

   omega = 7.292e-5 # Earth's rotation rate in 1/s)
   radius = 6371.0e3 # Radius of the earth

   lat   = np.zeros(mlat)
   f_d   = np.zeros(mlat)
   latb  = np.zeros(mlat+1)
   rlatb = np.zeros(mlat+1)

   latgrid  = np.zeros( (mlat,nlon) )
   longrid  = np.zeros( (mlat,nlon) )
   area  = np.zeros( (mlat,nlon) )

   cose  = np.zeros( mlat )
   sine  = np.zeros( mlat )
   cosp  = np.zeros( mlat )
   sinp  = np.zeros( mlat )
   acosp = np.zeros( mlat )
   acosu = np.zeros( mlat )
   grid_weight = np.zeros( mlat )

   coslon = np.zeros( nlon )
   sinlon = np.zeros( nlon )
   cosl5  = np.zeros( nlon )
   sinl5  = np.zeros( nlon )

   dx = np.zeros( mlat )
   dy = np.zeros( nlon )
   one_by_dx = np.zeros( mlat )
   one_by_dy = np.zeros( nlon )

   beglat = 0
   endlat = mlat-1

   pi = 4.0 * atan(1.0)
   dl = (pi+pi) / nlon
   dp = pi/(mlat-1)

   sine[1] = -1.
   for j in range(1, mlat):
      ph5  = -0.50*pi + ((j-1)-0.50)*(pi/(mlat-1))
      sine[j] = sin(ph5)

   for j in range(0, mlat-1):
      grid_weight[j] = sine[j+1] - sine[j]
   grid_weight[-1] = 1. - sine[-1]

   for j in range(0, mlat):
      cosp[j] = grid_weight[j] / dp
      acosp[j] = dp / grid_weight[j]
      grid_weight[j] = grid_weight[j] / (2*nlon)

   cose[0] = 0.
   cose[1] = 0.5 * cosp[2]
   for j in range(2, mlat-1):
      cose[j] = 0.5 * (cosp[j-1] + cosp[j])
   cose[-1] = cose[1]

   for j in range(1, mlat-1):
      acosu[j] = 2. / (cose[j] + cose[j+1])

   for j in range(0, mlat-1):
      sinp[j] = 0.5 * (sine[j] + sine[j+1])
   sinp[-1] = 0.5 * (sine[-1] + 1.)

   acap = nlon*(1.+sine[1])/dp               # acap = nlon * cosp[1]
   rcap = 1.0 / acap

   for j  in range(0,mlat):
      f_d[j] = (omega+omega) * sinp[j]

   for j  in range(1, mlat):
      latb[j] = -90. + ((j-1)-0.5) * 180. / (mlat-1)

   latb[ 0 ] = -90.
   latb[ mlat ] =  90.

   for j in range(0, mlat):
      lat[j] = 0.5 * ( latb[j] + latb[j+1] )

   one_by_dy = 1. / ( radius*dp )

   for j in range(1, mlat-1):
      dx[j] = dl*radius*cosp[j]
      one_by_dx[j] = 1. / dx[j]

   dg2r = pi/180.

   for j in range(0, mlat+1):
      rlatb[j] = dg2r * latb[j]

   for i in range(0, nlon+1):
      lonb[i] = (i-1)*360./(nlon)

   for i in range(0, nlon):
      lon[i] = 0.5 * ( lonb[i] + lonb[i+1] )

   for i in range(0, nlon+1):
      rlonb[i] = dg2r * lonb[i]

   imh = round(nlon/2)
   if(nlon  !=  4*(imh/2)) :
      raise NameError('FV_INIT','nlon must be dividable by 4')
   # endif

   for i in range(0, imh):
      zam5          = (i-1) * dl
      zamda         = (0.50+(i-1))*dl
      cosl5[i]      =  cos(zam5)
      cosl5[i+imh-1]  = -cosl5[i]
      sinl5[i]      =  sin(zam5)
      sinl5[i+imh-1]  = -sinl5[i]
      coslon[i]     =  cos(zamda)
      coslon[i+imh-1] = -coslon[i]
      sinlon[i]     =  sin(zamda)
      sinlon[i+imh-1] = -sinlon[i]

   for j in range(0,mlat):
      for i in range(0, nlon):
         latgrid[j,i] = dg2r * lat[j]
         longrid[j,i] = dg2r * lon[i]          # Cubed sphere ready
         area[j,i] = (radius*dl*cosp[j]) * (dp*radius)

   asum = area.mean()
   print("Mean cell width (km)=", 1.e-3*np.sqrt(asum))

   # Return quantities as 2D grids
   dx = np.tile(dx,(1,nlon))
   dy = np.tile(dy,(mlat,1))
   angle_dx = np.zeros(mlat,nlon)

   return dx, dy, latgrid, longrid, area, angle_dx
