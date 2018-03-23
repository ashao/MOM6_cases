## Make ocean and atmospheric supergrids. Needs numpy, midas, and netCDF4

import midas
import netCDF4 as nc
import numpy as np
import set_fv_geom

## Ocean Grid
# Define the parameters of a cartesian grid
xrefine = 1
yrefine = 1
lat0 = -90 # Start at the south pole
lon0 = 0   # Start the grid at the prime meridian
lenlat = 180 #total latitude range (-90 to 90)
lenlon = 360 #total longitude range (0 to 360)
# Define parameters related to topography
ocean_depth = 5000. # Depth of ocean points in meters
land_depth = -0.    # "depth" of the land points
cap_extent = 80 # How far equatorward the land at poles extends

# Number of points in the super grid (twice as much as actual model grid)
# 'round' used here in case of roundoff error
nx_super = round(360*xrefine) # Nominally 2-degree zonal resolution
ny_super = round(180*yrefine) # Nominally 2-degree meridional resolution

# Make a equally spaced lat-lon grid
ocean_grid = midas.rectgrid.supergrid(nx_super,ny_super,'spherical','degrees',lat0,lenlat,lon0,lenlon,cyclic_x=True)
ocean_grid.grid_metrics()
ocean_grid.write_nc('ocean_supergrid.nc')

# Make a tile variable and a string dimension...needed to generate exchange grids
fout=nc.Dataset('ocean_supergrid.nc','r+',format='NETCDF3_CLASSIC')
string=fout.createDimension('string',255)
tile=fout.createVariable('tile','S1',('string'))
tile[0:5]=nc.stringtochar(np.array(['tile1'],'S5'))
fout.sync()
fout.close()

# Generate topography based on a two pole planet
# Generate the T, U, V grid from the supergrid
grid = midas.rectgrid.quadmesh(supergrid = ocean_grid, is_latlon = True, cyclic = True)
ny_grid, nx_grid = grid.x_T.shape
# Initialize topography to land to start, shape based on size of a T-grid
topo = np.zeros(grid.x_T.shape) + land_depth
# Set ocean depth for all points equatorward of the cap_extent
topo[ np.abs(grid.y_T) <= cap_extent ] = ocean_depth

# Write out the two-pole topography
fout = nc.Dataset('twopole_topography.nc', 'w', format='NETCDF3_CLASSIC')
# Set x,y dimensions
yax = fout.createDimension('ny',ny_grid)
xax = fout.createDimension('nx',nx_grid)
# Create the 'ntiles' dimension needed for mosaic generation
fout.createDimension('ntiles',1)
# Create variables storing mean and variance of topography
meantopo = fout.createVariable('depth','f8',('ny','nx'))
meantopo.units = 'meters'
meantopo.standard_name = 'topographic depth at T-cell centers'
meantopo[:] = topo
stdtopo = fout.createVariable('std','f8',('ny','nx'))
stdtopo.units = 'meters'
stdtopo.standard_name = 'subgrid variance of topography at T-cell centers'
stdtopo[:] = topo * 0. # No variance of topography within a cell
fout.sync()
fout.close()

## Atmospheric Grid
# Define the parameters.
xrefine = 0.8 # 2.5-degree resolution
yrefine = 1.0 # 2-degree resolution
lat0 = -90 # Start at the south pole
lon0 = 0   # Start the grid at the prime meridian
lenlat = 180 #total latitude range (-90 to 90)
lenlon = 360 #total longitude range (0 to 360)
lenlat = 180 #total latitude range (-90 to 90)
lenlon = 360 #total longitude range (0 to 360)
# Define parameters related to topography

# Number of points in the super grid (twice as much as actual model grid)
# 'round' used here in case of roundoff error
nx_super = round(360*xrefine)
ny_super = round(180*yrefine)

# Make a equally spaced lat-lon grid
atmos_grid = midas.rectgrid.supergrid(nx_super,ny_super,'spherical','degrees',lat0,lenlat,lon0,lenlon,cyclic_x=True)
atmos_grid.grid_metrics()
# Overwrite the metrics with those from set_fv_geoms
print(atmos_grid.area.shape)
y, x = set_fv_geom.set_fv_geom(int(ny_super), int(nx_super))
print( x.shape, y.shape)
atmos_grid.x[:,:] = x
atmos_grid.y[:,:] = y
atmos_grid.angle_dx[:,:] = 0.
atmos_grid.write_nc('atmos_supergrid.nc')

# Make a tile variable and a string dimension...needed to generate exchange grids
fout=nc.Dataset('atmos_supergrid.nc','r+',format='NETCDF3_CLASSIC')
string=fout.createDimension('string',255)
tile=fout.createVariable('tile','S1',('string'))
tile[0:5]=nc.stringtochar(np.array(['tile1'],'S5'))
fout.sync()
fout.close()
