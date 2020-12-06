
program main

use LSS_tools

implicit none

        character(len=char_len) :: basename, printstr, tmpstr1, tmpstr2, filename
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


        print *, 'This is a sample of reading HR4 LC data... Please use it on baekdu!'
        printstr = 'Usage: EXE catnme'
        
        if(iargc().ne.1) then
                print *, printstr
                stop
        endif

        call getarg(1, basename)

        ! Name of directory
        basename ='/home/eostm/mergerTime/fullCatalogue/LC93/lightcone/data/'//trim(adjustl(basename))//'.dat.LC.'
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
                read(100) nowhalo
                print *, nowhalo
                print *
                close(100)
        enddo



        

end program main
