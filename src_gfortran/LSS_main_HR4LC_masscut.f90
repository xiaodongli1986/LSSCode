
! Applying minimal-mass and maximal-r cut to HR4 LC

program main

use LSS_tools

implicit none

        character(len=char_len) :: catname, basename, printstr, tmpstr1, tmpstr2, inputfilename, outputfilename
        integer :: i
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
        
        real(dl) :: masscut, rcut, maxrsq, rsq
        integer :: num_halo, num_acpt


        print *, 'This apply the masscut and distance cut to HR4 LC'
        printstr = 'Usage: EXE catname masscut rcut'
        
        if(iargc().ne.3) then
                print *, printstr
                stop
        endif

        call getarg(1, catname)
        call getarg(2, tmpstr1); read(tmpstr1, *) masscut
        call getarg(3, tmpstr2); read(tmpstr2, *) rcut

        ! Name of directory
        !print *, trim(adjustl(catname))
	! gfortran fmt
        basename ='/home/eostm/mergerTime/fullCatalogue/'//trim(adjustl(catname))//&
		'/lightcone/data/'//trim(adjustl(catname))//'.dat.LC.'
        outputfilename = &
                '/home/xiaodongli/SparseFilaments/data/input/HR4/LC/' &
                //trim(adjustl(catname)) &
                //'_massge'//trim(adjustl(tmpstr1))// &
                '_rle'//trim(adjustl(tmpstr2))
        print *, '##################################'
        print *, 'Applying minimal-mass and maximal-r cut to HR4 LC data; output = ', trim(adjustl(outputfilename))

        num_halo = 0; num_acpt = 0; maxrsq = rcut**2.0;
        open(unit=200,file=outputfilename)
        do i = 1, 16
                write(tmpstr1, *) i
                if (i .le. 9) then
                        inputfilename = trim(adjustl(basename))//'0'//trim(adjustl(tmpstr1))//'_16'
                else
                        inputfilename = trim(adjustl(basename))//trim(adjustl(tmpstr1))//'_16'
                endif
                print *, 'opening ', trim(adjustl(inputfilename))
                ! gfortran fmt
		!open(unit=100,file=inputfilename,form='binary',action='read')
		open(unit=100,file=inputfilename,form='unformatted',access='stream',action='read')
                do while(.true.)
                        read(100,end=1000) nh
                        num_halo = num_halo + 1
                        if(nh%mass < masscut) cycle
                        rsq = nh%x **2.0 + nh%y**2.0 + nh%z**2.0
                        if(rsq > maxrsq) cycle
                        write(200, '(7(e15.7))') nh%x, nh%y, nh%z, nh%vx, nh%vy, nh%vz, nh%mass
                        num_acpt = num_acpt + 1
!                        if(num_acpt > 10000) exit
                        cycle
1000                    exit
                enddo
                close(100)
        enddo
        close(200)

        write(*, '(A,i12,A,i12,A)') 'Among ', num_halo, ' halos, ', num_acpt, ' accepted satisfying our condition.' 
        

end program main
