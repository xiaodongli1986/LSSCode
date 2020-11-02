

module LSS_2pcf
use LSS_chisq
implicit none
contains

subroutine xyz2sigpi(v1,v2,sig,pi)
real(dl), dimension(3) :: v1,v2
real(dl):: sig,pi, dist2
 dist2=sqrt((0.5*(v1(1)+v2(1)))**2 + (0.5*(v1(2)+v2(2)))**2 + (0.5*(v1(3)+v2(3)))**2)
 sig=0.5*((v1(1)**2 - v2(1)**2)+(v1(2)**2 - v2(2)**2)+ (v1(3)**2 - v2(3)**2))
 sig=abs(sig/dist2)
!! theta=1.0-1.0/(4.0*dist2*dist2)
!! dist2=1.0+1.0/(4.0*dist2*dist2)
 pi= v1(1)**2 + v2(1)**2 - 2.*v1(1)*v2(1) &
    +v1(2)**2 + v2(2)**2 - 2.*v1(2)*v2(2) &
    +v1(3)**2 + v2(3)**2 - 2.*v1(3)*v2(3)
 pi=sqrt(abs(pi-sig*sig))
return
end subroutine xyz2sigpi


        subroutine Tpcf(inputfilename, rmin,rmax,numrbin, numtbin, counts, decomp, printinfo)
                !<< Dummy >>        
                character(len=char_len), intent(in) :: inputfilename
                real(dl), intent(in) :: rmin,rmax
                integer, intent(in)  :: numrbin,numtbin
                integer, intent(in)  :: decomp ! by default it is smu; decomp = 0; if decom ==1 then sigpi, numtbin 
                real(dl), intent(out)::  counts(numrbin, numtbin)
                logical, intent(in)  :: printinfo
                !<< Local >>
                real(16) ::  datanorm
                real(dl) ::  omegam,w, x,y,z, x2,y2,z2, sep,dist, theta, v1(3), v2(3), v3(3),v4(3),&
                         deltar,odeltar, deltar2,odeltar2, dist1,dist1sq, dist2,dist2sq, sig,pi,mind, rmaxfact
                integer :: i1, i2, i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, di1, irbin,itbin, sigbin,pibin, dijk

                gb_usenumdensity = .false.
                gb_dodensitynorm = .false.

                counts = 0.0_dl

                deltar = (rmax-rmin)/dble(numrbin)
                deltar2 = (rmax-rmin)/dble(numtbin)
                odeltar = 1.0_dl/deltar
                odeltar2 = 1.0_dl/deltar2

                print *, '  (Tpcf) Begins.'
		print *, '    numrbin / numtbin = ', numrbin, numtbin
		print *, '    deltar  / deltar2 = ', deltar, deltar2
                if (decomp .eq. 0) then
                    print *, '    (Ttcf) decompsition = smu'
                elseif (decomp .eq. 1) then
                    print *, '    (Ttcf) decompsition = sigpi'
		else
		    print *, '    (Tpcf ERROR) Wrong decomp! must be 0 (smu) or 1 (sigpi)'
		    stop
		endif
                gb_datafile = inputfilename
                gb_minimalr_cut = -1.0e30
                gb_maximalr_cut =  1.0e30
                omegam=om_dft; w=w_dft
                if(.true.) then
                        call cosmo_funs_init(printinfo)
                        call readin_dataran(printinfo)
                        call init_cosmo(omegam,w,h_dft,printinfo)
                        call init_mult_lists(printinfo)
                        !call do_cell_init((real(gb_numdata)**0.33_dl*2), printinfo)               
                        call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)               
                endif

                !do i = 1, 100
                !        print *, i, gb_mass_list(i)
                !enddo
                !stop

                if(printinfo) then
                        !print *, '  (Tpcf) contours = ', counts
                endif
                di1 = gb_numdata / 20
                do i1 = 1, gb_numdata
                       ! avoid redudant calculation
                       if(mod(i1,di1).eq.1) then
                               write(*,'(f6.3,A,$)') (i1/dble(gb_numdata))*100, '%==>'
                       endif
                       x=gb_xyz_list(1,i1)
                       y=gb_xyz_list(2,i1)
                       z=gb_xyz_list(3,i1)


                       rmaxfact = 1.42
                       imin = floor((x-rmax*rmaxfact-gbgridxmin)/gbdeltax +1.0_dl)
                       imax = floor((x+rmax*rmaxfact-gbgridxmin)/gbdeltax +1.0_dl)
                       jmin = floor((y-rmax*rmaxfact-gbgridymin)/gbdeltay +1.0_dl)
                       jmax = floor((y+rmax*rmaxfact-gbgridymin)/gbdeltay +1.0_dl)
                       kmin = floor((z-rmax*rmaxfact-gbgridzmin)/gbdeltaz +1.0_dl)
                       kmax = floor((z+rmax*rmaxfact-gbgridzmin)/gbdeltaz +1.0_dl)

                       dijk = 1
                       do i = max(1,imin-dijk), min(gb_n_cellx,imax+dijk)
                       do j = max(1,jmin-dijk), min(gb_n_celly,jmax+dijk)
                       do k = max(1,kmin-dijk), min(gb_n_cellz,kmax+dijk)
                          if (decomp .eq. 0) then
                                do l = 1, gb_cell_mat(i,j,k)%numdata
                                        i2 = gb_cell_mat(i,j,k)%idatalist(l)
                                        if(i2 .eq. i1) cycle

                                        x2 = gb_xyz_list(1,i2)
                                        y2 = gb_xyz_list(2,i2)
                                        z2 = gb_xyz_list(3,i2)
                                        
                                        v4 = [(x-x2), (y-y2), (z-z2)]
                                        sep  = dsqrt(v4(1)**2.0_dl + v4(2)**2.0_dl + v4(3)**2.0_dl)
                                        if(sep .ge. rmax .or. sep .le. rmin) cycle

                                        v3 = [(x+x2)*0.5_dl, (y+y2)*0.5_dl, (z+z2)*0.5_dl]
                                        !v3 = [x,y,z]
                                        dist = dsqrt(v3(1)**2.0_dl + v3(2)**2.0_dl + v3(3)**2.0_dl)
                                        theta = abs(v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)) / dist / sep
                                        if (theta > 1.0_dl) cycle
                                        !if (theta > 1.0_dl) theta = abs(mod(theta,1.0))

                                        irbin = ceiling((sep-rmin)*odeltar)
                                        itbin = ceiling(theta*numtbin)
                                !        print*, i1,i2, dist, theta, irbin, itbin, gb_mass_list(i1), gb_mass_list(i2)
                                        counts(irbin,itbin) = counts(irbin,itbin) + gb_mass_list(i1)*gb_mass_list(i2)
                                enddo
                                !if(i1>10) stop
                          elseif (decomp .eq. 1) then
                                dist1sq = x*x + y*y + z*z 
                                dist1 = dsqrt( dist1sq )
                                do l = 1, gb_cell_mat(i,j,k)%numdata
                                        i2 = gb_cell_mat(i,j,k)%idatalist(l)
                                        if(i2 .eq. i1) cycle

                                        x2 = gb_xyz_list(1,i2)
                                        y2 = gb_xyz_list(2,i2)
                                        z2 = gb_xyz_list(3,i2)
                                        dist2sq =  x2*x2 + y2*y2 + z2*z2 
                                        dist2   = dsqrt( dist2sq )

                                        !v4 = [(x-x2), (y-y2), (z-z2)]
                                        !sep  = dsqrt(v4(1)**2.0_dl + v4(2)**2.0_dl + v4(3)**2.0_dl)
                                        !if(sep .ge. rmax .or. sep .le. rmin) cycle
                                        !if(sep .ge. rmax .or. sep .le. rmin) cycle
                                        if(.true.) then

                                         v4 = [(x-x2), (y-y2), (z-z2)]
                                         sep  = dsqrt(v4(1)**2.0_dl + v4(2)**2.0_dl + v4(3)**2.0_dl)
                                         if(sep .ge. 1.42_dl*rmax .or. sep .le. rmin) cycle

                                         v3 = [(x+x2)*0.5_dl, (y+y2)*0.5_dl, (z+z2)*0.5_dl]
                                         dist = dsqrt(v3(1)**2.0_dl + v3(2)**2.0_dl + v3(3)**2.0_dl)
                                         theta = abs(v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)) / dist / sep
                                         !if (theta > 1.0_dl) cycle
                                         sig = sep*theta; pi = dsqrt(sep*sep - sig*sig)
					 !print *, sig, pi
					 !v1(1)=x;v1(2)=y;v1(3)=z;
					 !v2(1)=x2;v2(2)=y2;v2(3)=z2;
					 !call xyz2sigpi(v1,v2,sig,pi)
					 !print *, sig, pi

                                        else
                                         theta = abs((x*x2 + y*y2 + z*z2) / dist1 / dist2 )
                                         theta = mod(theta, 1.0)
                                         mind=min(dist1,dist2)
                                         sig=abs(dist1-dist2)
                                         pi=dsqrt(2.d0*(mind*mind - mind*mind*theta))
                                        endif
 
                                        
                                        if(sig > rmax .or. pi > rmax) cycle
                                        sigbin = ceiling((sig-rmin)*odeltar)
                                        pibin  = ceiling((pi-rmin)*odeltar2)
                                        !itbin = ceiling(theta*numtbin)
                                        !sigbin=floor(sig*odeltar)+1
                                        !pibin=floor(pi*odeltar2)+1

                                        counts(sigbin,pibin) = counts(sigbin,pibin) + gb_mass_list(i1)*gb_mass_list(i2)
                                enddo
                          endif

                        enddo
                        enddo
                        enddo
                enddo
                datanorm = sum(gb_mass_list(1:gb_numdata))**2.0 * (gb_numdata-1)/dble(gb_numdata)

                counts = counts/datanorm

                print*, 'result: '
                print*, counts                          
        end subroutine Tpcf
end module LSS_2pcf
