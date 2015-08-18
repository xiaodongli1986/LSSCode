
        


program main

use LSS_tools

implicit none

        character(len=char_len) :: catname, basename, printstr, tmpstr1, tmpstr2, filename
        integer :: i,j
        type savedhalotype 
                integer*8 :: mbp
                integer*8 :: nowhid, majorglobalid
                double precision :: x, y, z, hx, hy, hz
                real :: mass, fofhmass
                real :: vx, vy, vz, hvx, hvy, hvz
                integer*2 :: statusflag, linkflag
        end type savedhalotype
        type(savedhalotype) :: nh


!typedef struct SavedHaloType{
!        INDXTYPE mbp; /* mbp ID */
!        IDTYPE nowhid,majorglobalid; /* nowhid is halo id at the redshift and 
!majorglobalid is global array id to the major subhalo 
!                                                                        (or mbp !) */
!        POSTYPE x,y,z,hx,hy,hz; 
!                /* in unit of h^-1 Mpc in comoving coordinate 
!                                                         (x,y,z) is position of !the mbp and (hx,hy,hz) is the CoM position of the hosting FoF halo */
!        float mass,fofhmass; 
!                /* in unit of M_sun/h 
!                                                  mass is the mass of MBP or halo (if one mbp is in the halo) just before entering larger halo region.
!                                                  fofhmass is the hosting halo mass. */
!        float vx,vy,vz,hvx,hvy,hvz; 
!                /* in unit of km/second 
!                                                                 (vx,vy,vz) is the peculiar velocity of the mbp &
!                                                                 (hvx,hvy,hvz) is the peculiar velocity of the hosting FoF halo. */
!        short int statusflag, linkflag;
!} SavedHaloType;

        ! values
        integer, parameter :: nval=14
        character(len=char_len), parameter :: &
                name1 = 'x', &
                name2 = 'y', &
                name3 = 'z', &
                name4 = 'hx', &
                name5 = 'hy', &
                name6 = 'hz', &
                name7 = 'vx', &
                name8 = 'vy', &
                name9 = 'vz', &
                name10= 'hvx', &
                name11= 'hvy', &
                name12= 'hvz', &
                name13= 'mass', &
                name14= 'fofhmass' 
     
        character(len=char_len) :: valnames(nval) = &
                (/name1, name2, name3, name4, name5, name6,name7, name8, name9, &
                        name10, name11, name12, name13, name14/)
        real(dl) :: x,y,z,hx,hy,hz,vx,vy,vz,hvx,hvy,hvz,m,fm, vals(nval)

        ! histogram
        integer, parameter :: ncell = 40, nmasscell = ncell*100
        real(dl), parameter :: xyzmin=-4000.0, xyzmax=4000.0, vmin=-800, vmax=800, massmin = 1.0e0, massmax = 1.0e20
        integer :: xyzhist(ncell,ncell,ncell)=0, vxyzhist(ncell,ncell,ncell)=0, &
                hxyzhist(ncell,ncell,ncell)=0, hvxyzhist(ncell,ncell,ncell)=0, &
                mhist(nmasscell)=0, fmhist(nmasscell)=0, ix,iy,iz,im, &
                nvout, nhvout

        ! means, sq-means
        real*16 :: means(nval), sqmeans(nval)
        integer*8 :: num_halo, nm0, nfm0
 
        means = 0; sqmeans = 0; num_halo = 0

        ! print information
        print *, 'This is a sample of reading HR4 LC data... Please use it on baekdu!'
        printstr = 'Usage: EXE catnme'
        
        if(iargc().ne.1) then
                print *, printstr
                stop
        endif

        call getarg(1, catname)

        ! Name of directory
!        basename ='/home/eostm/mergerTime/fullCatalogue/LC93/lightcone/'//trim(adjustl(catname))//'.dat.LC.'
        basename ='/home/eostm/mergerTime/fullCatalogue/'//trim(adjustl(catname))//'/lightcone/data/'//trim(adjustl(catname))//'.dat.LC.'
        num_halo = 0; nm0=0; nfm0=0
        nvout = 0; nhvout = 0;
        do i = 1, 16
                write(tmpstr1, *) i
                if (i .le. 9) then
                        filename = trim(adjustl(basename))//'0'//trim(adjustl(tmpstr1))//'_16'
                else
                        filename = trim(adjustl(basename))//trim(adjustl(tmpstr1))//'_16'
                endif
                print *, '##################################'
                print *, 'opening ', trim(adjustl(filename))
                open(unit=100,file=filename,form='binary',action='read')!,action='read',access='stream',status='old')!,recl=8)!,iostat=ok)
                do while (.true.)
                        read(100,end=100) nh
                        if(.false.) then
                                print *, '################'
                                print *, num_halo, 'th halo:'
                                print *, nh
                        endif
                        num_halo = num_halo + 1
                        ! vals
                        x=nh%x;y=nh%y;z=nh%z;hx=nh%hx;hy=nh%hy;hz=nh%hz;
                        vx=nh%vx;vy=nh%vy;vz=nh%vz;hvx=nh%hvx;hvy=nh%hvy;hvz=nh%hvz;
                        m=nh%mass;fm=nh%fofhmass
                        vals = (/ x,y,z,hx,hy,hz,vx,vy,vz,hvx,hvy,hvz,m,fm/)
                        ! means, sqmeans
                        do j=1, nval
                                means(j)=means(j)+vals(j)
                                sqmeans(j)=sqmeans(j)+vals(j)**2.0
                        enddo
                        ! hists
                        call  xi3d(x,y,z,ix,iy,iz,xyzmin,xyzmax,ncell)
                        xyzhist(ix,iy,iz)=xyzhist(ix,iy,iz)+1 
                        call  xi3d(hx,hy,hz,ix,iy,iz,xyzmin,xyzmax,ncell)
                        hxyzhist(ix,iy,iz)=hxyzhist(ix,iy,iz)+1 
                        if( vmin<vx.and.vx<vmax .and. vmin<vy.and.vy<vmax .and. vmin<vz.and.vz<vmax) then
                                call  xi3d(vx,vy,vz,ix,iy,iz,vmin,vmax,ncell)
                                vxyzhist(ix,iy,iz)=vxyzhist(ix,iy,iz)+1 
                        else
                                nvout = nvout + 1
                        endif
                        if( vmin<hvx.and.hvx<vmax .and. vmin<hvy.and.hvy<vmax .and. vmin<hvz.and.hvz<vmax) then
                                call  xi3d(hvx,hvy,hvz,ix,iy,iz,vmin,vmax,ncell)
                                hvxyzhist(ix,iy,iz)=hvxyzhist(ix,iy,iz)+1 
                        else
                                nhvout = nhvout + 1
                        endif
                        if(m.ne.0.0_dl) then
                                im=logmi(m, massmin, massmax, nmasscell)
                                mhist(im)=mhist(im)+1
                        else
                                nm0 = nm0+1
                        endif
                        if(fm.ne.0.0_dl) then
                                im=logmi(fm, massmin, massmax, nmasscell)
                                fmhist(im)=fmhist(im)+1
                        else
                                nfm0 = nfm0+1
                        endif
!                       if(num_halo > 1000000) exit
                        cycle
100                     exit
                enddo                
                close(100)
                print *, 'current #: ', num_halo, '...'
        enddo

        do j = 1, nval
                means(j) = means(j) / dble(num_halo)
                sqmeans(j) = sqmeans(j) / dble(num_halo)
        enddo
        print *, 'Finishing processing ', num_halo, 'halos...'


        ! means
        write(tmpstr1, *) ncell
        catname = trim(adjustl(catname))//'.ncell'//trim(adjustl(tmpstr1))
        filename = trim(adjustl(catname))//'.meansqmeans'
        open(file=filename,unit=100)
        write(100,'(A,i10,A)') '# Name   mean   mean-square   RMS   var #',&
                num_halo, ' halos processed'
        do j = 1, nval
                write(100,'(A,3x,4(e15.7,1x))') trim(adjustl(valnames(j))), means(j), &
                        sqmeans(j), sqrt(sqmeans(j)), sqmeans(j)-means(j)**2.0
        enddo
        close(100)

        ! hists
        filename = trim(adjustl(catname))//'.xyzhist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,2e15.7)') '# i j k num #',&
                num_halo, ' halos processed; min, max = ', xyzmin, xyzmax
        open(file=filename,unit=100)
        do ix = 1,ncell
        do iy = 1,ncell
        do iz = 1,ncell
                write(100,'(3i5,i12)'), ix,iy,iz, xyzhist(ix,iy,iz)
        enddo
        enddo
        enddo
        close(100)

        
        filename = trim(adjustl(catname))//'.vxyzhist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,2e15.7,i10,A)') '# i j k num #',&
                num_halo, ' halos processed; min, max = ', vmin, vmax, &
                nvout, ' lying out of range. '
        open(file=filename,unit=100)
        do ix = 1,ncell
        do iy = 1,ncell
        do iz = 1,ncell
                write(100,'(3i5,i12)'), ix,iy,iz, vxyzhist(ix,iy,iz)
        enddo
        enddo
        enddo
        close(100)

        filename = trim(adjustl(catname))//'.hxyzhist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,2e15.7)') '# i j k num #',&
                num_halo, ' halos processed; min, max = ', xyzmin, xyzmax
        open(file=filename,unit=100)
        do ix = 1,ncell
        do iy = 1,ncell
        do iz = 1,ncell
                write(100,'(3i5,i12)'), ix,iy,iz, hxyzhist(ix,iy,iz)
        enddo
        enddo
        enddo
        close(100)

        
        filename = trim(adjustl(catname))//'.hvxyzhist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,2e15.7,i10,A)') '# i j k num #',&
                num_halo, ' halos processed; min, max = ', vmin, vmax, &
                nhvout, ' lying out of range. '
        open(file=filename,unit=100)
        do ix = 1,ncell
        do iy = 1,ncell
        do iz = 1,ncell
                write(100,'(3i5,i12)'), ix,iy,iz, hvxyzhist(ix,iy,iz)
        enddo
        enddo
        enddo
        close(100)

        filename = trim(adjustl(catname))//'.masshist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,i15,A)') '# m_min, m_max, num #',&
                num_halo, ' halos processed; ', nm0, 'has mass=0'
        do im = 1,nmasscell
                write(100,'(2e20.7,3x,i12)'), 10.0**(log10(massmin) + (log10(massmax)-log10(massmin))*(im-1)/dble(nmasscell)), &
                        10.0**(log10(massmin) + (log10(massmax)-log10(massmin))*(im)/dble(nmasscell)), mhist(im)
        enddo
        close(100)

        filename = trim(adjustl(catname))//'.fofhmasshist'
        open(file=filename,unit=100)
        write(100,'(A,i15,A,i15,A)') '# m_min, m_max, num #',&
                num_halo, ' halos processed; ', nfm0, 'has fofhmass=0'
        do im = 1,nmasscell
                write(100,'(2e20.7,3x,i12)'), 10.0**(log10(massmin) + (log10(massmax)-log10(massmin))*(im-1)/dble(nmasscell)), &
                        10.0**(log10(massmin) + (log10(massmax)-log10(massmin))*(im)/dble(nmasscell)), fmhist(im)
        enddo
        close(100)

end program main
