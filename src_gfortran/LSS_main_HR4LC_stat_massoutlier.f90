
program main

use LSS_tools

implicit none

        character(len=char_len) :: basename, printstr, tmpstr1, tmpstr2, filename, outputfile1, outputfile2
        integer :: i
        type savedhalotype 
                integer*8 :: mbp
                integer*8 :: nowhid, majorglobalid
                double precision :: x, y, z, hx, hy, hz
                real :: mass, fofhmass
                real :: vx, vy, vz, hvx, hvy, hvz
                integer*2 :: statusflag, linkflag
        end type savedhalotype
        type(savedhalotype) :: nowhalo

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

        integer :: i1,i2
        real(dl) :: x

        print *, 'This searches for the outlier-mass HR4 LC data... Study the mass/r distribution... Please use it on baekdu!'
        printstr = 'Usage: EXE catnme'
        
        if(iargc().ne.1) then
                print *, printstr
                stop
        endif

        call getarg(1, basename)

        ! Name of directory
        outputfile1 = trim(adjustl(basename))//'.massoutlier'
        outputfile2 = trim(adjustl(basename))//'.masswithin'
        basename ='/home/eostm/mergerTime/fullCatalogue/'//trim(adjustl(basename))//&
		'/lightcone/data/'//trim(adjustl(basename))//'.dat.LC.'
        open(unit=201,file=outputfile1)
        open(unit=202,file=outputfile2)
        do i = 1, 16
                write(tmpstr1, *) i
                if (i .le. 9) then
                        filename = trim(adjustl(basename))//'0'//trim(adjustl(tmpstr1))//'_16'
                else
                        filename = trim(adjustl(basename))//trim(adjustl(tmpstr1))//'_16'
                endif
                print *, '##################################'
                print *, 'opening ', trim(adjustl(filename))
                !open(unit=100,file=filename,form='binary',action='read')!,action='read',access='stream',status='old')!,recl=8)!,iostat=ok) ! gfortran fmt
                open(unit=100,file=filename,form='unformatted',access='stream',action='read')!,action='read',access='stream',status='old')!,recl=8)!,iostat=ok)
                i1=0; i2=0
                do while(.true.)
                        read(100,end=100) nowhalo
                        if (nowhalo%mass .le. 0.27E+12) then
                                call random_number(x)
                                if(x.le.0.02) then
                                        write(201,'(7e15.7)') nowhalo%x, nowhalo%y, nowhalo%z, &
                                                nowhalo%vx, nowhalo%vy, nowhalo%vz, nowhalo%mass
                                        i1 = i1+1
                                endif
                        else
                                call random_number(x)
                                if(x.le.0.002) then
                                        write(202,'(7e15.7)') nowhalo%x, nowhalo%y, nowhalo%z, &
                                                nowhalo%vx, nowhalo%vy, nowhalo%vz, nowhalo%mass
                                        i2 = i2+1
                                endif
                        endif
                        cycle
100                     exit 
                enddo               
                close(100)
                close(201)
                close(202)
        enddo



        

end program main
