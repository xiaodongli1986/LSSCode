
module illustris_cic
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

        subroutine cic_features(xyzs,features,nfeature,xyzmin,xyzmax,ndat,nc,featuregrid)
            ! output featuregrid; first nfeatures columns are the different features; last column is the number density 
                integer, intent(in) :: ndat, nc, nfeature
                real, intent(in) :: xyzs(3,ndat), features(nfeature,ndat),xyzmin,xyzmax
                real, intent(out):: featuregrid(nc,nc,nc,nfeature+1)  ! the last column has wei=1 (i.e. number density)
                integer :: i,ix(3),ix0(3),ix1(3), ifeature
                real :: w1,w2,w3,w4,w5,w6,w7,w8, x,y,z, dx,dx_inv, xyz(3), w(3),wb(3), wei
                dx = (xyzmax-xyzmin)/real(nc); dx_inv = 1./dx
                write(*,'(A,i6,2f15.7)') '   (cic_features) nc, dx, dx_inv = ',nc, dx, dx_inv
                do i = 1, ndat
                        xyz = xyzs(:,i); x=xyzs(1,i);y=xyzs(2,i);z=xyzs(3,i)
                        if(x<xyzmin.or.x>xyzmax.or.y<xyzmin.or.y>xyzmax.or.z<xyzmin.or.z>xyzmax) cycle
                        ix=int(xyz*dx_inv); 
                        !print *, 'xyz= ', xyz
                        w=1. - (xyz*dx_inv -ix); wb=1.-w
                        ix0 = mod((ix+nc),nc)+1; ix1 = mod((ix+1+nc),nc)+1
                        !print *, 'ix0, ix1 = ', ix0, ix1
                        if(mod(i,1000000).eq.0) print *, '  (cic_features) processing', i,'-th particle...'
                        w1 = w(1)*w(2)*w(3); w2 = w(1)*wb(2)*w(3); w3 = w(1)*w(2)*wb(3); w4 = w(1)*wb(2)*wb(3);
                        w5 = wb(1)*w(2)*w(3); w6 = wb(1)*wb(2)*w(3); w7 = wb(1)*w(2)*wb(3); w8 = wb(1)*wb(2)*wb(3);
                        !print *, 'weights = ', w1,w2,w3,w4,w5,w6,w7,w8
                        !if(i.eq.3) stop
                        do ifeature = 1, nfeature+1
                             if(ifeature .le. nfeature) then
                                wei = features(ifeature, i)
                             else
                                wei = 1.0
                             endif
                             !print *, 'i, ifeature, wei = ', i, ifeature, wei
                             featuregrid(ix0(3),ix0(2),ix0(1),ifeature) = featuregrid(ix0(3),ix0(2),ix0(1),ifeature) + w1*wei
                             featuregrid(ix0(3),ix1(2),ix0(1),ifeature) = featuregrid(ix0(3),ix1(2),ix0(1),ifeature) + w2*wei
                             featuregrid(ix1(3),ix0(2),ix0(1),ifeature) = featuregrid(ix1(3),ix0(2),ix0(1),ifeature) + w3*wei
                             featuregrid(ix1(3),ix1(2),ix0(1),ifeature) = featuregrid(ix1(3),ix1(2),ix0(1),ifeature) + w4*wei
                             featuregrid(ix0(3),ix0(2),ix1(1),ifeature) = featuregrid(ix0(3),ix0(2),ix1(1),ifeature) + w5*wei
                             featuregrid(ix0(3),ix1(2),ix1(1),ifeature) = featuregrid(ix0(3),ix1(2),ix1(1),ifeature) + w6*wei
                             featuregrid(ix1(3),ix0(2),ix1(1),ifeature) = featuregrid(ix1(3),ix0(2),ix1(1),ifeature) + w7*wei
                             featuregrid(ix1(3),ix1(2),ix1(1),ifeature) = featuregrid(ix1(3),ix1(2),ix1(1),ifeature) + w8*wei
                        enddo
                enddo
        end subroutine cic_features

        subroutine cic_features_arbitary_position(xyzs,features,nfeature,xmin,ymin,zmin,boxsize,ndat,nc,featuregrid)
            ! output featuregrid; first nfeatures columns are the different features; last column is the number density 
                integer, intent(in) :: ndat, nc, nfeature
                real, intent(in) :: xyzs(3,ndat), features(nfeature,ndat),xmin,ymin,zmin,boxsize
                real, intent(out):: featuregrid(nc,nc,nc,nfeature+1)  ! the last column has wei=1 (i.e. number density)
                integer :: i,ix(3),ix0(3),ix1(3), ifeature
                real :: w1,w2,w3,w4,w5,w6,w7,w8, x,y,z, dx,dx_inv, w(3),wb(3), wei, dxyz(3), xmax,ymax,zmax
                xmax=xmin+boxsize; ymax=ymin+boxsize; zmax=zmin+boxsize
                dx = boxsize/real(nc); dx_inv = 1./dx
                write(*,'(A,i6,2f15.7)') '   (cic_features) nc, dx, dx_inv = ',nc, dx, dx_inv
                do i = 1, ndat
                        x=xyzs(1,i);y=xyzs(2,i);z=xyzs(3,i)
                        if(x<xmin.or.x>xmax.or.y<ymin.or.y>ymax.or.z<zmin.or.z>zmax) cycle
                        dxyz(1)=x-xmin; dxyz(2)=y-ymin; dxyz(3)=z-zmin 
                        ix=int(dxyz*dx_inv); 
                        !print *, 'xyz= ', xyz
                        w=1. - (dxyz*dx_inv -ix); wb=1.-w

                        !ix0 = mod((ix+nc),nc)+1; ix1 = mod((ix+1+nc),nc)+1 !! periodic boundary condition

                        ix0 = min(ix+1,nc); ix1 = min((ix+2),nc) !! no periodic boundary condition

                        !ix0(0) = ix; ix1 = ix
                        !print *, 'ix0, ix1 = ', ix0, ix1
                        if(mod(i,1000000).eq.0) print *, '  (cic_features) processing', i,'-th particle...'
                        w1 = w(1)*w(2)*w(3); w2 = w(1)*wb(2)*w(3); w3 = w(1)*w(2)*wb(3); w4 = w(1)*wb(2)*wb(3);
                        w5 = wb(1)*w(2)*w(3); w6 = wb(1)*wb(2)*w(3); w7 = wb(1)*w(2)*wb(3); w8 = wb(1)*wb(2)*wb(3);
                        !print *, 'weights = ', w1,w2,w3,w4,w5,w6,w7,w8
                        !if(i.eq.3) stop
                        do ifeature = 1, nfeature+1
                             if(ifeature .le. nfeature) then
                                wei = features(ifeature, i)
                             else
                                wei = 1.0
                             endif
                             !print *, 'i, ifeature, wei = ', i, ifeature, wei
                             featuregrid(ix0(3),ix0(2),ix0(1),ifeature) = featuregrid(ix0(3),ix0(2),ix0(1),ifeature) + w1*wei
                             featuregrid(ix0(3),ix1(2),ix0(1),ifeature) = featuregrid(ix0(3),ix1(2),ix0(1),ifeature) + w2*wei
                             featuregrid(ix1(3),ix0(2),ix0(1),ifeature) = featuregrid(ix1(3),ix0(2),ix0(1),ifeature) + w3*wei
                             featuregrid(ix1(3),ix1(2),ix0(1),ifeature) = featuregrid(ix1(3),ix1(2),ix0(1),ifeature) + w4*wei
                             featuregrid(ix0(3),ix0(2),ix1(1),ifeature) = featuregrid(ix0(3),ix0(2),ix1(1),ifeature) + w5*wei
                             featuregrid(ix0(3),ix1(2),ix1(1),ifeature) = featuregrid(ix0(3),ix1(2),ix1(1),ifeature) + w6*wei
                             featuregrid(ix1(3),ix0(2),ix1(1),ifeature) = featuregrid(ix1(3),ix0(2),ix1(1),ifeature) + w7*wei
                             featuregrid(ix1(3),ix1(2),ix1(1),ifeature) = featuregrid(ix1(3),ix1(2),ix1(1),ifeature) + w8*wei
                        enddo
                enddo
        end subroutine cic_features_arbitary_position
end module illustris_cic

program main_illustris_cic

use LSS_cosmo_funs
use illustris_cic

implicit none

    character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputname, outputfile, printstr, inputfiles(10000)
    integer :: i,j,k, i1,i2,j1,j2,k1,k2, ifile, nfile, nc, tmpint(64), &
       block1, block2, ntotal, nnow, numarg, chbk, nsplit, n
    real :: xmin,ymin,zmin, boxsize, tmpfloats(64), xyz_rescale, xyz(3)! , xshift=0., yshift=0., zshift=0.
    real, allocatable:: rgrid(:,:,:), featuregrid(:,:,:,:), all_data(:,:)
    logical :: do_rgrid = .false.

    character(len=char_len) :: filelistfile
    real(dl) :: tmpx, tmpy
    integer :: size=1, tmpi, ifeature
    integer, allocatable :: seed(:)

    integer :: nfeature, npar

    printstr = "Usage: LSS_illustris_cic -input inputfile   -nc n-cells   -xmin xmin   -ymin ymin   -zmin zmin  "//&
        " -boxsize boxsize   -npar npar   -nfeature nfeature  "//&
        "-output outputname -xyz_rescale xyz_rescale   -nsplit nsplit    -readin_fmt  ascii/gadget   "//&
        "#### This Code Does not use periodic BC!   ## Example: LSS_illustris_cic ... "
    numarg = iargc()
    if(numarg.le.1) then
        write(*,'(A)') printstr
        stop
    endif

    xyz_rescale = 1; xmin=0; ymin=0; zmin=0; nsplit=1

    do i = 1, numarg
        if(mod(i,2).eq.0) cycle
        call getarg(i,tmpstr1)
        call getarg(i+1,tmpstr2)
        if(trim(adjustl(tmpstr1)).eq.'-input') then
            read(tmpstr2,'(A)') inputfile
        elseif(trim(adjustl(tmpstr1)).eq.'-nc') then
            read(tmpstr2,*) nc
        elseif(trim(adjustl(tmpstr1)).eq.'-xmin') then
            read(tmpstr2,*) xmin
        elseif(trim(adjustl(tmpstr1)).eq.'-ymin') then
            read(tmpstr2,*) ymin
        elseif(trim(adjustl(tmpstr1)).eq.'-zmin') then
            read(tmpstr2,*) zmin
        elseif(trim(adjustl(tmpstr1)).eq.'-boxsize') then
            read(tmpstr2,*) boxsize
        elseif(trim(adjustl(tmpstr1)).eq.'-output') then
            read(tmpstr2,'(A)') outputname
        elseif(trim(adjustl(tmpstr1)).eq.'-xyz_rescale') then
            read(tmpstr2,*) xyz_rescale
        elseif(trim(adjustl(tmpstr1)).eq.'-nsplit') then
            read(tmpstr2,*) nsplit
        elseif(trim(adjustl(tmpstr1)).eq.'-npar') then
            read(tmpstr2,*) npar
        elseif(trim(adjustl(tmpstr1)).eq.'-nfeature') then
            read(tmpstr2,*) nfeature
!        elseif(trim(adjustl(tmpstr1)).eq.'-yshift') then
!            read(tmpstr2,*) yshift
!        elseif(trim(adjustl(tmpstr1)).eq.'-zshift') then
!            read(tmpstr2,*) zshift
        else
            print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
            write(*,'(A)') trim(adjustl(printstr))
            stop
        endif
    enddo

    print *, '(LSS_illustris_cic) CIC feature fields for illustris simulation'
    write(*,'(A,A)')  '    inputfile  = ', trim(adjustl(inputfile))
    write(*,'(A,A)')  '    outputname = ', trim(adjustl(outputname))
    print *, '   nc (celss) = ', nc
    print *, '   xyz_rescale= ', xyz_rescale
    print *, '   x,y,z min= ', xmin, ymin, zmin
    print *, '   boxsize= ', boxsize
    print *, '   do_rgrid   = ', do_rgrid
    print *, '   npar, nfeature = ', npar, nfeature
    write(*,'(A,i3,A,i8,A)'), '    nsplit^3   = ', nsplit, '^3 = ', nsplit**3, ' sets of outupt files.'
    allocate(all_data( nfeature, npar))

    call random_seed(size=size); allocate(seed(size)); call random_seed(put=seed)
    call random_number(tmpx); call random_number(tmpy); tmpi = int(tmpx*100000000 + tmpy * 100000)
    write(filelistfile, *) tmpi; filelistfile= 'illustris_subboxsplit_filelist_'//trim(adjustl(filelistfile))//'.tmp'
    call system("ls "//trim(adjustl(inputfile))//' > '//trim(adjustl(filelistfile)));
    !call system("ls "//trim(adjustl(inputfile))//' > cic_tmp_filelists.tmp ');

    open(file=trim(adjustl(filelistfile)),action="read",unit=100)
    nfile = 0
    do while(.true.)
        nfile = nfile+1; read(100,'(A)',end=100) inputfiles(nfile)
!       print *, nfile, trim(adjustl(inputfiles(nfile)))
        cycle
100     exit
    enddo
    close(100); nfile = nfile-1
    print *, '  (LSS_illustris_cic) found ', nfile, 'files to read-in...'
    do i = 1, nfile
        write(*,'(i5,5x,A)') i, trim(adjustl(inputfiles(i)))
    enddo

    write(*,*) ' (LSS_illustris_cic) allocating grid fields (ncol= ',nfeature-3+1,')...'
    allocate(featuregrid(nc,nc,nc,nfeature-3+1))
    ntotal=0

    featuregrid=0.
    do ifile = 1, nfile
        print *
        write(*,'(A,A,A)'), ' (LSS_illustris_cic) open ', trim(adjustl(inputfiles(ifile))), '...'
        if(.true.) then
            open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='unformatted',access='stream') 
            read(10001) all_data
            ntotal = ntotal+nnow

            print *, '(LSS_illustris_cic) finish read in all_data'
        
            !call cic_features(all_data(1:3,:),all_data(4:nfeature,:),nfeature-3,xyzmin,xyzmax,npar,nc,featuregrid)

            call cic_features_arbitary_position(all_data(1:3,:),all_data(4:nfeature,:),nfeature-3,&
                    xmin,ymin,zmin,boxsize,npar,nc,featuregrid)
        endif
                        
    enddo
    print * 
    write(*,*) '(LSS_illustris_cic) Finishing cic. In total we read-in', npar, 'particles in ', nfile, 'files.'
    do i = 1, 10
        print *, i, featuregrid(i,i,i,:)
    enddo

!    if(do_rgrid) then
!        allocate(rgrid(nc,nc,nc))
!        call cic_rgrid(rgrid, xyzmin,xyzmax, nc)
!    endif

    write(*,*) '(LSS_illustris_cic) output...'
    if(nsplit.eq.1) then
         do ifeature = 1, nfeature -3+1
             write(tmpstr1,*) nc; outputfile = trim(adjustl(outputname))//'_'//trim(adjustl(tmpstr1))//'grid' 
             write(tmpstr1,*) ifeature; outputfile = trim(adjustl(outputfile))//'_feature'//trim(adjustl(tmpstr1))
             write(*,'(A,A)') ' (LSS_illustris_cic) saving result to ', trim(adjustl(outputfile))
             open(file=trim(adjustl(outputfile)),unit=2001,form='unformatted',action='write')
             write(2001) featuregrid(:,:,:,ifeature); close(2001)
         enddo


         if(do_rgrid)  then
               write(tmpstr1,*) nc; outputfile = trim(adjustl(outputname))//'_'//trim(adjustl(tmpstr1))//'grid' 
               outputfile = trim(adjustl(outputfile))//'_rgrid'
               open(file=trim(adjustl(outputfile)),unit=2001,form='unformatted',action='write')
               write(2001) rgrid; close(2001)
         endif
            
    endif

    print *, '(LSS_illustris_cic) rm tmp file: ', trim(adjustl(filelistfile))
    call system("rm "//trim(adjustl(filelistfile)))

end program
