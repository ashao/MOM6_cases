
subroutine set_fv_geom
! Sets the geometric factors used in FV core for the lat-lon grid
! This routine should only be called once per run

   integer i, j, imh
   real dp, dl
   real ph5, zamda, zam5
   real dg2r, asum, gmean

   allocate (  lon(nlon) )
   allocate (  lonb(nlon+1) )
   allocate ( rlonb(nlon+1) )

   allocate (   lat(mlat) )
   allocate (  latb(mlat+1) )
   allocate ( rlatb(mlat+1) )

   allocate ( rlat(nlon,beglat:endlat) )
   allocate ( rlon(nlon,beglat:endlat) )
   allocate ( area(nlon,beglat:endlat) )

   allocate ( cose(mlat) )
   allocate ( sine(mlat) )
   allocate ( cosp(mlat) )
   allocate ( sinp(mlat) )
   allocate (acosp(mlat) )
   allocate (acosu(mlat) )
   allocate ( grid_weight(mlat) )

   allocate ( coslon(nlon) )
   allocate ( sinlon(nlon) )

   allocate ( cosl5(nlon) )
   allocate ( sinl5(nlon) )

   allocate ( dx(mlat) )
   allocate ( one_by_dx(mlat) )


   pi = 4.d0 * atan(1.d0)
   dl = (pi+pi) / nlon
   dp = pi/(mlat-1)

   sine(1) = -1.
   do j=2,mlat
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(mlat-1))
      sine(j) = sin(ph5)
   enddo

   do j=1,mlat-1
      grid_weight(j) = sine(j+1) - sine(j)
   enddo
      grid_weight(mlat) = 1. - sine(mlat)

! Dne area-conservative cosine at "cell center"
   do j=1,mlat
             cosp(j) = grid_weight(j) / dp
            acosp(j) = dp / grid_weight(j)
      grid_weight(j) = grid_weight(j) / (2*nlon)
   enddo

! Dne cosine at edges..

      cose(1) = 0.
      cose(2) = 0.5 * cosp(2)
   do j=3,mlat-1
      cose(j) = 0.5 * (cosp(j-1) + cosp(j))
   enddo
      cose(mlat) = cose(2)

   do j=2,mlat-1
      acosu(j) = 2. / (cose(j) + cose(j+1))
   enddo

   do j=1,mlat-1
      sinp(j) = 0.5 * (sine(j) + sine(j+1))
   enddo
      sinp(mlat) = 0.5 * (sine(mlat) + 1.)

   acap = nlon*(1.+sine(2))/dp               ! acap = nlon * cosp(1)
   rcap = 1.d0 / acap

   do j = max(1,beglat-ng_d), min(mlat,endlat+ng_d)
      f_d(j) = (omega+omega) * sinp(j)
   enddo

   do j = 2, mlat
      latb(j) = -90. + (real(j-1)-0.5) * 180. / real(mlat-1)
   enddo

   latb( 1 ) = -90.
   latb( mlat+1 ) =  90.

   do j=1, mlat
      lat(j) = 0.5 * ( latb(j) + latb(j+1) )
   enddo

   one_by_dy = 1. / ( radius*dp )

   do j=2, mlat-1
             dx(j) = dl*radius*cosp(j)
      one_by_dx(j) = 1. / dx(j)
   enddo

         dg2r = pi/180.

   do j=1,mlat+1
      rlatb(j) = dg2r * latb(j)
   enddo

   do i=1,nlon+1
      lonb(i) = real(i-1)*360./real(nlon)
   enddo

   do i=1,nlon
      lon(i) = 0.5 * ( lonb(i) + lonb(i+1) )
   enddo

   do i=1,nlon+1
      rlonb(i) = dg2r * lonb(i)
   enddo


     imh = nlon/2
   if(nlon /= 4*(imh/2)) then
      call error_mesg('FV_INIT','nlon must be dividable by 4', FATAL)
   endif

   do i=1,imh
      zam5          = (i-1) * dl
      zamda         = (0.5d0+(i-1))*dl
      cosl5(i)      =  cos(zam5)
      cosl5(i+imh)  = -cosl5(i)
      sinl5(i)      =  sin(zam5)
      sinl5(i+imh)  = -sinl5(i)
      coslon(i)     =  cos(zamda)
      coslon(i+imh) = -coslon(i)
      sinlon(i)     =  sin(zamda)
      sinlon(i+imh) = -sinlon(i)
   enddo

   do j=beglat, endlat
      do i=1,nlon
         rlat(i,j) = dg2r * lat(j)
         rlon(i,j) = dg2r * lon(i)          ! Cubed sphere ready
         area(i,j) = (radius*dl*cosp(j)) * (dp*radius)
      enddo
   enddo

      asum = gmean(nlon, mlat, beglat, endlat, area)
   if(master) then
      write(6,*) 'Mean cell width (km)=', 1.e-3*sqrt(asum)
   endif

end subroutine set_fv_geom
