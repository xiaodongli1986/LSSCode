
module cic
implicit none
contains
        
        integer function check_blocks(block1,block2)
                integer :: block1, block2
                if(block1.eq.block2) then
                      !print *, '  (check_blocks) consistent blocks:', block1, block2
                      check_blocks=1
                else
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      check_blocks=0
                endif
        end function check_blocks

!      for i in range(ndat):
!        X, V = pos(i,:), vel(i,:)
!        ix = (X * dx_inv).astype('int32')
!        w = 1. - (X*dx_inv - ix)
!        ix0 = (ix+nc)%nc; ix1 = (ix+1+nc)%nc 
!        if i%50000 == 0: print('\t\t',i,'-th particle...') 
!        w1, w2, w3, w4, w5, w6, w7, w8 = w(0)*w(1)*w(2), w(0)*(1-w(1))*w(2), w(0)*w(1)*(1-w(2)), w(0)*(1-w(1))*(1-w(2)), \
!            (1-w(0))*w(1)*w(2), (1-w(0))*(1-w(1))*w(2), (1-w(0))*w(1)*(1-w(2)), (1-w(0))*(1-w(1))*(1-w(2))
!        for grid, wei in (
!                (rhogrid, 1),
!                (vxgrid, V(0)),
!                (vygrid, V(1)),
!                (vzgrid, V(2)),
!                ):
!          grid(ix0(0),  ix0(1), ix0(2)) += w1*wei
!          grid(ix0(0),  ix1(1), ix0(2)) += w2*wei
!          grid(ix0(0),  ix0(1), ix1(2)) += w3*wei
!          grid(ix0(0),  ix1(1), ix1(2)) += w4*wei
!          grid(ix1(0),  ix0(1), ix0(2)) += w5*wei
!          grid(ix1(0),  ix1(1), ix0(2)) += w6*wei
!          grid(ix1(0),  ix0(1), ix1(2)) += w7*wei
!          grid(ix1(0),  ix1(1), ix1(2)) += w8*wei


        subroutine cic_rgrid(rgrid, xyzmin,xyzmax, nc)
                integer, intent(in) :: nc
                real, intent(in) :: xyzmin,xyzmax
                real, intent(out) :: rgrid(nc,nc,nc)
                integer :: i1,i2,i3 
                real :: dx,dx_inv
                dx = (xyzmax-xyzmin)/real(nc); dx_inv = 1./dx
                write(*,'(A,i6,2f15.7)') '   (cic_xyz_vs) nc, dx, dx_inv = ',nc, dx, dx_inv
                do i1=1,nc
                do i2=1,nc 
                do i3=1,nc
                        rgrid(i3,i2,i1) =  ((xyzmin+dx*(i1-0.5))**2. + (xyzmin+dx*(i2-0.5))**2. + (xyzmin+dx*(i3-0.5))**2.)**0.5
                enddo
                enddo
                enddo
        end subroutine cic_rgrid

                
                

        subroutine cic_xyz_vs(xyzs,vs,rhogrid,vxgrid,vygrid,vzgrid,xyzmin,xyzmax,ndat,nc)
                integer, intent(in) :: ndat, nc
                real, intent(in) :: xyzs(3,ndat), vs(3,ndat),xyzmin,xyzmax
                real, intent(out):: rhogrid(nc,nc,nc), vxgrid(nc,nc,nc), vygrid(nc,nc,nc), vzgrid(nc,nc,nc)
                integer :: i,ix(3),ix0(3),ix1(3)
                real :: w1,w2,w3,w4,w5,w6,w7,w8, x,y,z,vx,vy,vz, dx,dx_inv, xyz(3), w(3),wb(3), wei
                dx = (xyzmax-xyzmin)/real(nc); dx_inv = 1./dx
                write(*,'(A,i6,2f15.7)') '   (cic_xyz_vs) nc, dx, dx_inv = ',nc, dx, dx_inv
                do i = 1, ndat
                        xyz = xyzs(:,i); x=xyzs(1,i);y=xyzs(2,i);z=xyzs(3,i)
                        if(x<xyzmin.or.x>xyzmax.or.y<xyzmin.or.y>xyzmax.or.z<xyzmin.or.z>xyzmax) cycle
                        vx=vs(1,i);   vy=vs(2,i);   vz=vs(3,i); 
                        ix=int(xyz*dx_inv); 
                        !print *, 'xyz= ', xyz
                        w=1. - (xyz*dx_inv -ix); wb=1.-w
                        ix0 = mod((ix+nc),nc)+1; ix1 = mod((ix+1+nc),nc)+1
                        !print *, 'ix0, ix1 = ', ix0, ix1
                        if(mod(i,1000000).eq.0) print *, '  (cic_xyz_vs) processing', i,'-th particle...'
                        w1 = w(1)*w(2)*w(3); w2 = w(1)*wb(2)*w(3); w3 = w(1)*w(2)*wb(3); w4 = w(1)*wb(2)*wb(3);
                        w5 = wb(1)*w(2)*w(3); w6 = wb(1)*wb(2)*w(3); w7 = wb(1)*w(2)*wb(3); w8 = wb(1)*wb(2)*wb(3);
                        wei = 1.0
                        !print *, 'weights = ', w1,w2,w3,w4,w5,w6,w7,w8
                        !if(i.eq.3) stop
                        rhogrid(ix0(3),ix0(2), ix0(1)) = rhogrid(ix0(3),ix0(2), ix0(1)) + w1*wei
                        rhogrid(ix0(3),ix1(2), ix0(1)) = rhogrid(ix0(3),ix1(2), ix0(1)) + w2*wei
                        rhogrid(ix1(3),ix0(2), ix0(1)) = rhogrid(ix1(3),ix0(2), ix0(1)) + w3*wei
                        rhogrid(ix1(3),ix1(2), ix0(1)) = rhogrid(ix1(3),ix1(2), ix0(1)) + w4*wei
                        rhogrid(ix0(3),ix0(2), ix1(1)) = rhogrid(ix0(3),ix0(2), ix1(1)) + w5*wei
                        rhogrid(ix0(3),ix1(2), ix1(1)) = rhogrid(ix0(3),ix1(2), ix1(1)) + w6*wei
                        rhogrid(ix1(3),ix0(2), ix1(1)) = rhogrid(ix1(3),ix0(2), ix1(1)) + w7*wei
                        rhogrid(ix1(3),ix1(2), ix1(1)) = rhogrid(ix1(3),ix1(2), ix1(1)) + w8*wei
                        wei = vx
                        vxgrid(ix0(3),ix0(2), ix0(1)) = vxgrid(ix0(3),ix0(2), ix0(1)) + w1*wei
                        vxgrid(ix0(3),ix1(2), ix0(1)) = vxgrid(ix0(3),ix1(2), ix0(1)) + w2*wei
                        vxgrid(ix1(3),ix0(2), ix0(1)) = vxgrid(ix1(3),ix0(2), ix0(1)) + w3*wei
                        vxgrid(ix1(3),ix1(2), ix0(1)) = vxgrid(ix1(3),ix1(2), ix0(1)) + w4*wei
                        vxgrid(ix0(3),ix0(2), ix1(1)) = vxgrid(ix0(3),ix0(2), ix1(1)) + w5*wei
                        vxgrid(ix0(3),ix1(2), ix1(1)) = vxgrid(ix0(3),ix1(2), ix1(1)) + w6*wei
                        vxgrid(ix1(3),ix0(2), ix1(1)) = vxgrid(ix1(3),ix0(2), ix1(1)) + w7*wei
                        vxgrid(ix1(3),ix1(2), ix1(1)) = vxgrid(ix1(3),ix1(2), ix1(1)) + w8*wei
                        wei = vy
                        vygrid(ix0(3),ix0(2), ix0(1)) = vygrid(ix0(3),ix0(2), ix0(1)) + w1*wei
                        vygrid(ix0(3),ix1(2), ix0(1)) = vygrid(ix0(3),ix1(2), ix0(1)) + w2*wei
                        vygrid(ix1(3),ix0(2), ix0(1)) = vygrid(ix1(3),ix0(2), ix0(1)) + w3*wei
                        vygrid(ix1(3),ix1(2), ix0(1)) = vygrid(ix1(3),ix1(2), ix0(1)) + w4*wei
                        vygrid(ix0(3),ix0(2), ix1(1)) = vygrid(ix0(3),ix0(2), ix1(1)) + w5*wei
                        vygrid(ix0(3),ix1(2), ix1(1)) = vygrid(ix0(3),ix1(2), ix1(1)) + w6*wei
                        vygrid(ix1(3),ix0(2), ix1(1)) = vygrid(ix1(3),ix0(2), ix1(1)) + w7*wei
                        vygrid(ix1(3),ix1(2), ix1(1)) = vygrid(ix1(3),ix1(2), ix1(1)) + w8*wei
                        wei = vz
                        vzgrid(ix0(3),ix0(2), ix0(1)) = vzgrid(ix0(3),ix0(2), ix0(1)) + w1*wei
                        vzgrid(ix0(3),ix1(2), ix0(1)) = vzgrid(ix0(3),ix1(2), ix0(1)) + w2*wei
                        vzgrid(ix1(3),ix0(2), ix0(1)) = vzgrid(ix1(3),ix0(2), ix0(1)) + w3*wei
                        vzgrid(ix1(3),ix1(2), ix0(1)) = vzgrid(ix1(3),ix1(2), ix0(1)) + w4*wei
                        vzgrid(ix0(3),ix0(2), ix1(1)) = vzgrid(ix0(3),ix0(2), ix1(1)) + w5*wei
                        vzgrid(ix0(3),ix1(2), ix1(1)) = vzgrid(ix0(3),ix1(2), ix1(1)) + w6*wei
                        vzgrid(ix1(3),ix0(2), ix1(1)) = vzgrid(ix1(3),ix0(2), ix1(1)) + w7*wei
                        vzgrid(ix1(3),ix1(2), ix1(1)) = vzgrid(ix1(3),ix1(2), ix1(1)) + w8*wei

                enddo
        end subroutine cic_xyz_vs

end module cic

program main_cic

use LSS_cosmo_funs
use cic

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputname, printstr, inputfiles(10000), readin_fmt_str
        integer :: i,j,k, i1,i2,j1,j2,k1,k2, ifile, nfile, nc, &
		tmpint(64), block1, block2, ntotal, nnow, numarg, chbk, nsplit, n, readin_fmt
        real :: xyzmin, xyzmax, tmpfloats(64), xyz_rescale, xyz(3), vxvyvz(3)  !, xshift=0., yshift=0., zshift=0.
        real, allocatable:: rhogrid(:,:,:), vxgrid(:,:,:), vygrid(:,:,:), vzgrid(:,:,:) , xyzs(:,:), vs(:,:), rgrid(:,:,:)
        integer, parameter :: fmt_gadget = 1, fmt_ascii = 2
        logical :: do_rgrid = .false.

        type head
                integer :: npart(6)
                double precision :: mass(6)
                double precision :: time
                double precision :: redshift
                integer :: flag_sfr, flag_feedback, npartTotal(6), flag_cooling, num_files
                double precision :: BoxSize, Omega0, OmegaLambda, HubbleParam;
                integer :: fill(24)
        end type head
        type(head) :: headinfo

	!headinfo%BoxSize=1.0

	
	printstr = "Usage: LSS_cic -input inputfile   -nc n-cells   -xyzmin xyzmin   -xyzmax xyzmax   "//&
		"-output outputname -xyz_rescale xyz_rescale   -nsplit nsplit    -readin_fmt  ascii/gadget   -do_grid T/F. "//&
		"#### Example: LSS_cic -input snp04000e.\? -nc 512 -output snp04000e_cic"
        numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') printstr
		stop
	endif

        xyz_rescale = 1; xyzmin=0; xyzmax=0; nsplit=1

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-nc') then
			read(tmpstr2,*) nc
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmin') then
			read(tmpstr2,*) xyzmin
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmax') then
			read(tmpstr2,*) xyzmax
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
			read(tmpstr2,'(A)') outputname
		elseif(trim(adjustl(tmpstr1)).eq.'-xyz_rescale') then
			read(tmpstr2,*) xyz_rescale
		elseif(trim(adjustl(tmpstr1)).eq.'-nsplit') then
			read(tmpstr2,*) nsplit
!		elseif(trim(adjustl(tmpstr1)).eq.'-xshift') then
!			read(tmpstr2,*) xshift
!		elseif(trim(adjustl(tmpstr1)).eq.'-yshift') then
!			read(tmpstr2,*) yshift
!		elseif(trim(adjustl(tmpstr1)).eq.'-zshift') then
!			read(tmpstr2,*) zshift
		elseif(trim(adjustl(tmpstr1)).eq.'-readin_fmt') then
			read(tmpstr2,'(A)') readin_fmt_str
                        if(trim(adjustl(readin_fmt_str)).eq.'gadget') then
                                readin_fmt = fmt_gadget
                        elseif(trim(adjustl(readin_fmt_str)).eq.'ascii') then
                                readin_fmt = fmt_ascii
                        else
                                print *, ' (LSS_cic) ERROR! wrong readin_fmt_str: ', trim(adjustl(readin_fmt_str))
                                stop
                        endif
                elseif(trim(adjustl(tmpstr1)).eq.'-do_rgrid') then
                        read(tmpstr2,*) do_rgrid
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	print *, 'LSS_cic, CIC density and velocity fields'
        write(*,'(A,A)')  '  inputfile  = ', trim(adjustl(inputfile))
        write(*,'(A,A)')  '  outputname = ', trim(adjustl(outputname))
        print *, ' nc (celss) = ', nc
        print *, ' xyz_rescale= ', xyz_rescale
        print *, ' xyz range  = ', xyzmin, xyzmax
        print *, ' do_rgrid   = ', do_rgrid
        print *, ' readin_fmt = ', trim(adjustl(readin_fmt_str))
        write(*,'(A,i3,A,i8,A)'), '  nsplit^3   = ', nsplit, '^3 = ', nsplit**3, ' sets of outupt files.'


        call system("ls "//trim(adjustl(inputfile))//' > cic_tmp_filelists.tmp ');

        open(file="cic_tmp_filelists.tmp",action="read",unit=100)
        nfile = 0
        do while(.true.)
                nfile = nfile+1; read(100,'(A)',end=100) inputfiles(nfile)
!                print *, nfile, trim(adjustl(inputfiles(nfile)))
                cycle
100             exit
        enddo
        close(100); nfile = nfile-1
        print *, ' found ', nfile, 'files to read-in...'
        do i = 1, nfile
                write(*,'(i5,5x,A)') i, trim(adjustl(inputfiles(i)))
        enddo

        write(*,*) 'allocating rho/px/py/pz fields...'
        allocate(rhogrid(nc,nc,nc),vxgrid(nc,nc,nc),vygrid(nc,nc,nc),vzgrid(nc,nc,nc))
        ntotal=0

        rhogrid = 0.; vxgrid=0.; vygrid=0.; vzgrid=0.
        do ifile = 1, nfile
                print *
                write(*,'(A)') ' ##########################'
                write(*,'(A,A,A)'), ' open ', trim(adjustl(inputfiles(ifile))), '...'
                if(readin_fmt .eq. fmt_gadget) then
                        !open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='binary',access='sequential') !gfortran fmt
			! warning!  binary sequential -> unformatted stream
                        open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='unformatted',access='stream') 
                        read(10001) block1, headinfo, block2
			! goftran fmt !headinfo%boxsize = (headinfo%boxsize) * xyz_rescale
			headinfo%boxsize = (headinfo%boxsize) * xyz_rescale
                        chbk = check_blocks(block1,block2)
                        if(ifile.eq.1) then
                                write(*,'(A,3f14.5)')'   om /ol / h         = ', &
					headinfo%Omega0, headinfo%OmegaLambda, headinfo%HubbleParam
                                write(*,'(A,3f14.5)')'   redshift / boxsize = ', &
					headinfo%redshift, headinfo%boxsize
                                if(xyzmax.eq.0.) then
                                        xyzmax = headinfo%boxsize
                                        print *, '  xyz range (reset)  = ', xyzmin, xyzmax
                                endif
                        endif
                        nnow = 0;
                        do i = 1, 6
                                if(headinfo%npart(i).ne.0) then
                                        write(*,'("   In total",i12,A,i3,A,e14.7,A,i12,A)') &
						headinfo%npart(i), ' particles with type ', i, ', mass',&
                                                headinfo%mass(i),' (',headinfo%nparttotal(i),') in total'
                                        nnow = nnow+headinfo%npart(i)
                                endif
                        enddo
                        ntotal = ntotal+nnow
        
                        allocate(xyzs(3,nnow),vs(3,nnow))
        
                        read(10001) block1
                        read(10001) xyzs; xyzs=xyzs*xyz_rescale
                        read(10001) block2
                        chbk = check_blocks(block1,block2)
                        write(*,'(A,3f10.3,5x,3f10.3)') '     begin and end of xyzs: ', xyzs(:,1), xyzs(:,nnow)

                        read(10001) block1
                        read(10001) vs
                        read(10001) block2
                        chbk = check_blocks(block1,block2)
                        write(*,'(A,3f10.3,5x,3f10.3)') '     begin and end of   vs: ', vs(:,1), vs(:,nnow)
                        close(10001)
                elseif(readin_fmt .eq. fmt_ascii) then
                        call count_line_number(inputfiles(ifile), nnow)
                        open(unit=10001,file=trim(adjustl(inputfiles(ifile))),action='read')
                        allocate(xyzs(3,nnow),vs(3,nnow))
                        do i = 1, nnow
                                read(10001,*) xyzs(1:3,i), vs(1:3,i)
                        enddo
                        close(10001)
                        ntotal = ntotal + nnow
                endif


                ! run cic
                call cic_xyz_vs(xyzs,vs,rhogrid,vxgrid,vygrid,vzgrid,xyzmin,xyzmax,nnow,nc)

                deallocate(xyzs,vs)
                        
        enddo
        print * 
        print *, '#########################'
        write(*,*) 'Finishing cic. In total we read-in', ntotal, 'particles in ', nfile, 'files.'
        if(readin_fmt .eq. fmt_gadget) then
                if(ntotal.eq.sum(headinfo%nparttotal)) then
                        print *, 'Consistent. ntotal, sum(headinfo%nparttotal) = ', ntotal, sum(headinfo%nparttotal)
                else
                        print *, ' (WARNING!!!) Inconsistent: ntotal, sum(headinfo%nparttotal) = ', ntotal, sum(headinfo%nparttotal)
                        print *, ' (WARNING!!!) Inconsistent: ntotal, sum(headinfo%nparttotal) = ', ntotal, sum(headinfo%nparttotal)
                        print *, ' (WARNING!!!) Inconsistent: ntotal, sum(headinfo%nparttotal) = ', ntotal, sum(headinfo%nparttotal)
                endif
        endif
        do i = 1, 10
                print *, i, rhogrid(i,i,i)
        enddo

        if(do_rgrid) then
                allocate(rgrid(nc,nc,nc))
                call cic_rgrid(rgrid, xyzmin,xyzmax, nc)
        endif


        print *, '#########################'
        write(*,*) 'output...'
        if(nsplit.eq.1) then
          open(file=trim(adjustl(outputname))//'.rhogrid',unit=2001,form='unformatted',action='write')
          open(file=trim(adjustl(outputname))//'.pxgrid',unit=2002,form='unformatted',action='write')
          open(file=trim(adjustl(outputname))//'.pygrid',unit=2003,form='unformatted',action='write')
          open(file=trim(adjustl(outputname))//'.pzgrid',unit=2004,form='unformatted',action='write')
          write(2001) rhogrid
          write(2002) vxgrid
          write(2003) vygrid
          write(2004) vzgrid
          close(2001); close(2002); close(2003); close(2004)
          if(do_rgrid)  then
                open(file=trim(adjustl(outputname))//'.rgrid',unit=2005,form='unformatted',action='write')
                write(2005) rgrid; close(2005)
          endif
        else
          n = nc / nsplit
          print *, ' nc, nsplit, n = ', nc, nsplit, n
          do i=1,nsplit
          do j=1,nsplit
          do k=1,nsplit
             ifile = k+(j-1)*nsplit+(i-1)*nsplit*nsplit
             i1 = (i-1)*n+1; i2=i1+n-1
             j1 = (j-1)*n+1; j2=j1+n-1
             k1 = (k-1)*n+1; k2=k1+n-1
             write(tmpstr1,*) ifile
             write(*,'(3x,A,A,A,2i4,3x,2i4,3x,2i4)') 'write to ', trim(adjustl(outputname)), '. i/j/k range = ', i1,i2, j1,j2, k1,k2
             open(file=trim(adjustl(outputname))//'.ifile'//trim(adjustl(tmpstr1))//'.rhogrid',&
	        unit=2001,form='unformatted',action='write')
             open(file=trim(adjustl(outputname))//'.ifile'//trim(adjustl(tmpstr1))//'.pxgrid',&
	         unit=2002,form='unformatted',action='write')
             open(file=trim(adjustl(outputname))//'.ifile'//trim(adjustl(tmpstr1))//'.pygrid',&
	         unit=2003,form='unformatted',action='write')
             open(file=trim(adjustl(outputname))//'.ifile'//trim(adjustl(tmpstr1))//'.pzgrid',&
	         unit=2004,form='unformatted',action='write')
             write(2001) rhogrid(k1:k2,j1:j2,i1:i2)
             write(2002) vxgrid(k1:k2,j1:j2,i1:i2)
             write(2003) vygrid(k1:k2,j1:j2,i1:i2)
             write(2004) vzgrid(k1:k2,j1:j2,i1:i2)
             close(2001); close(2002); close(2003); close(2004)
             if(do_rgrid)  then
                open(file=trim(adjustl(outputname))//'.ifile'//trim(adjustl(tmpstr1))//'.rgrid',&
		    unit=2005,form='unformatted',action='write')
                write(2005) rgrid(k1:k2,j1:j2,i1:i2)
                close(2005)
             endif
           enddo
           enddo
           enddo

            
        endif


end program
