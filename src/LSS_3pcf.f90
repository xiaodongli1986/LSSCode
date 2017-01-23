

module LSS_3pcf
use LSS_chisq
implicit none
contains
        subroutine Threepcf(inputfilename, rmin,rmax,numrbin, numtbin, counts, decomp, printinfo)
                !<< Dummy >>        
                character(len=char_len), intent(in) :: inputfilename
                real(dl), intent(in) :: rmin,rmax
                integer, intent(in)  :: numrbin,numtbin
                integer, intent(in)  :: decomp ! by default it is smu; decomp = 0; if decom ==1 then sigpi, numtbin 
                real(dl), intent(out)::  counts(numrbin, numtbin)
                logical, intent(in)  :: printinfo
                !<< Local >>
                real(16) ::  datanorm
                real(dl) ::  omegam,w, x,y,z, x2,y2,z2, x3,y3,z3, r1,r2,r3, sep,dist, theta, v3(3),v4(3),vr1(3),vr2(3),vr3(3),& !v3 is distance
                         deltar,odeltar, deltar2,odeltar2, dist1,dist1sq, dist2,dist2sq, sig,pi,mind
                integer :: i1, i2, i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, di1, irbin,itbin, sigbin,pibin,xyzextra,&
			ii,jj,kk,l2,l3,i3, i2min,i2max,j2min,j2max,k2min,k2max, ir1bin,ir2bin,ir3bin

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
                di1 = 50000
                do i1 = 1, gb_numdata
                       ! avoid redudant calculation
                       if(mod(i1,di1).eq.1) then
                               write(*,'(f6.3,A,$)') (i1/dble(gb_numdata))*100, '%==>'
                       endif
                       x=gb_xyz_list(1,i1)
                       y=gb_xyz_list(2,i1)
                       z=gb_xyz_list(3,i1)

                       xyzextra = 0
                       imin = int((x-rmax-gbgridxmin)/gbdeltax +1.0_dl-xyzextra)
                       imax = int((x+rmax-gbgridxmin)/gbdeltax +1.0_dl+xyzextra)
                       jmin = int((y-rmax-gbgridymin)/gbdeltay +1.0_dl-xyzextra)
                       jmax = int((y+rmax-gbgridymin)/gbdeltay +1.0_dl+xyzextra)
                       kmin = int((z-rmax-gbgridzmin)/gbdeltaz +1.0_dl-xyzextra)
                       kmax = int((z+rmax-gbgridzmin)/gbdeltaz +1.0_dl+xyzextra)
                       !print *, xyzextra; stop

                       do i = max(1,imin), min(gb_n_cellx,imax)
                       do j = max(1,jmin), min(gb_n_celly,jmax)
                       do k = max(1,kmin), min(gb_n_cellz,kmax)
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
                                        !theta = abs(v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)) / dist / sep
                                        theta = abs(x*v4(1) + y*v4(2) + z*v4(3)) / sqrt(x*x+y*y+z*z) / sep
                                        if (theta > 1.0_dl) cycle
                                        !if (theta > 1.0_dl) theta = abs(mod(theta,1.0))

                                        irbin = ceiling((sep-rmin)*odeltar)
                                        itbin = ceiling(theta*numtbin)
                                !        print*, i1,i2, dist, theta, irbin, itbin, gb_mass_list(i1), gb_mass_list(i2)
		                        do ii = max(1,imin), min(gb_n_cellx,imax)
                		        do jj = max(1,jmin), min(gb_n_celly,jmax)
		                        do kk = max(1,kmin), min(gb_n_cellz,kmax)
					do l3 = 1, gb_cell_mat(ii,jj,kk)%numdata
                                        	i3 = gb_cell_mat(ii,jj,kk)%idatalist(l3)
                                                !if(i3 <= i2) cycle
	                                        if(i3.eq.i1 .or. i3.eq.i2) cycle

                                        	x3 = gb_xyz_list(1,i3)
                                        	y3 = gb_xyz_list(2,i3)
                                        	z3 = gb_xyz_list(3,i3)

	                                	vr1 = [(x-x2), (y-y2), (z-z2)]
	                                	vr2 = [(x-x3), (y-y3), (z-z3)]
	                                	vr3 = [(x3-x2), (y3-y2), (z3-z2)]
						r1 = dsqrt(vr1(1)**2.0_dl + vr1(2)**2.0_dl + vr1(3)**2.0_dl)
						r2 = dsqrt(vr2(1)**2.0_dl + vr2(2)**2.0_dl + vr2(3)**2.0_dl)
						r3 = dsqrt(vr3(1)**2.0_dl + vr3(2)**2.0_dl + vr3(3)**2.0_dl)
                                        	if(r1 .ge. rmax .or. r1 .le. rmin .or. r2 .ge. rmax .or. r2 .le. rmin &
							.or. r3 .ge. rmax .or. r3 .le. rmin) cycle
						ir1bin = ceiling((r1-rmin)*odeltar)
						ir2bin = ceiling((r2-rmin)*odeltar)
						ir3bin = ceiling((r3-rmin)*odeltar)
						if(ir1bin .ne. ir2bin .or. ir1bin.ne.ir3bin) cycle
						irbin = ir1bin
                                        	counts(irbin,itbin) = counts(irbin,itbin) + gb_mass_list(i1)*gb_mass_list(i2)*gb_mass_list(i3)
					enddo
					enddo
					enddo
					enddo


                                        !counts(irbin,itbin) = counts(irbin,itbin) + gb_mass_list(i1)*gb_mass_list(i2)
                                enddo
                                !if(i1>10) stop
                          elseif (decomp .eq. 1) then
				print *, "ERROR (3pCF)! decomp = 1 not supporte!"
				stop
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
                                        if(.false.) then

                                         v4 = [(x-x2), (y-y2), (z-z2)]
                                         sep  = dsqrt(v4(1)**2.0_dl + v4(2)**2.0_dl + v4(3)**2.0_dl)
                                         if(sep .ge. 1.415_dl*rmax .or. sep .le. rmin) cycle

                                         v3 = [(x+x2)*0.5_dl, (y+y2)*0.5_dl, (z+z2)*0.5_dl]
                                         dist = dsqrt(v3(1)**2.0_dl + v3(2)**2.0_dl + v3(3)**2.0_dl)
                                         theta = abs(v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)) / dist / sep
                                         sig = sep*theta; pi = dsqrt(sep*sep - sig*sig)

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
                datanorm = sum(gb_mass_list(1:gb_numdata))**3.0  ! * (gb_numdata-1)/dble(gb_numdata)

                counts = counts/datanorm

                print*, 'result: '
                print*, counts                          
        end subroutine Threepcf
end module LSS_3pcf
