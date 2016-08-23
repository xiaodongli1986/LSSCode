
module LSS_2pcf
use LSS_chisq
implicit none
contains
        subroutine Tpcf(inputfilename, rmin,rmax,numrbin, numtbin, counts, printinfo)
                !<< Dummy >>        
                character(len=char_len), intent(in) :: inputfilename
                real(dl), intent(in) :: rmin,rmax
                integer, intent(in)  :: numrbin,numtbin
                real(16), intent(out)::  counts(numrbin, numtbin)
                logical, intent(in)  :: printinfo
                !<< Local >>
                real(16) ::  datanorm
                real(dl) ::  omegam,w, x,y,z, x2,y2,z2, sep,dist, theta, v3(3),v4(3),&
                         deltar,odeltar,deltat,odeltat
                integer :: i1, i2, i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, di1, irbin,itbin

                gb_usenumdensity = .false.
                gb_dodensitynorm = .false.

                counts = 0.0_dl

                deltar = (rmax-rmin)/dble(numrbin)
                deltat = const_pi / dble(numtbin)
                odeltar = 1.0_dl/deltar
                odeltat = 1.0_dl/deltat

                print *, '  (Tpcf) Begins.'
                gb_datafile = inputfilename
                gb_minimalr_cut = -1.0e30
                gb_maximalr_cut =  1.0e30
                omegam=om_dft; w=w_dft
                if(.true.) then
                        call cosmo_funs_init(printinfo)
                        call readin_dataran(printinfo)
                        call init_cosmo(omegam,w,h_dft,printinfo)
                        call init_mult_lists(printinfo)
                        call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)               
                endif

                !do i = 1, 100
                !        print *, i, gb_mass_list(i)
                !enddo
                !stop

                if(printinfo) then
                        print *, '  (Tpcf) contours = ', counts
                endif
                di1 = 10000
                do i1 = 1, gb_numdata
                       ! avoid redudant calculation
                       if(mod(i1,di1).eq.1) then
                               write(*,'(f6.3,A,$)') (i1/dble(gb_numdata))*100, '%==>'
                       endif
                       x=gb_xyz_list(1,i1)
                       y=gb_xyz_list(2,i1)
                       z=gb_xyz_list(3,i1)

                       imin = int((x-rmax-gbgridxmin)/gbdeltax +1.0_dl)
                       imax = int((x+rmax-gbgridxmin)/gbdeltax +1.0_dl)
                       jmin = int((y-rmax-gbgridymin)/gbdeltay +1.0_dl)
                       jmax = int((y+rmax-gbgridymin)/gbdeltay +1.0_dl)
                       kmin = int((z-rmax-gbgridzmin)/gbdeltaz +1.0_dl)
                       kmax = int((z+rmax-gbgridzmin)/gbdeltaz +1.0_dl)

                       do i = max(1,imin), min(gb_n_cellx,imax)
                       do j = max(1,jmin), min(gb_n_celly,jmax)
                       do k = max(1,kmin), min(gb_n_cellz,kmax)
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
