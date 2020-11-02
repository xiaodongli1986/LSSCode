
program main

use LSS_tools

implicit none

        character(len=char_len) :: basename, printstr, tmpstr1, tmpstr2, filename, outputfile1, outputfile2
        integer :: i, j
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

        real(dl) :: masscut1, masscut2, masscutmin, x0,y0,z0, x,y,z, rmaxsq
        integer * 8 :: nhalo, nwrite1, nwrite2, nmasscut1

        print *, 'This reads J08 Slice data, do mass and distance cut'
        printstr = 'Usage: EXE basename'
        
        if(iargc().ne.1) then
                print *, printstr
                stop
        endif

        call getarg(1, basename)

        ! Name of Files
        masscut1 = 4.05d12; outputfile1 = trim(adjustl(basename))//'.mge4.05e12.00'
        masscut2 = 8.0d12; outputfile2 = trim(adjustl(basename))//'.mge8e12.00'

        ! Spherical Region
        rmaxsq = 1808.4_dl **2.0; 

        print *, 'Only save halos within a sphere with maximal radius ', real(sqrt(rmaxsq))

        filename='/home/eostm/mergerTime/fullCatalogue/J08/slice/'//trim(adjustl(basename))
        print *, 'Reading from ', trim(adjustl(filename))
        print *, 'Doing masscut radius cut; mass cutted halos within spherical region written to files:'
        print *, masscut1, trim(adjustl(outputfile1))
        print *, masscut2, trim(adjustl(outputfile2))

        nhalo = 0; nwrite1 = 0; nwrite2 = 0; nmasscut1 = 0;
        masscutmin = min(masscut1, masscut2)
        !open(unit=100,file=filename,form='binary',action='read') ! gfortran fmt
        open(unit=100,file=filename,form='unformatted',access='stream',action='read')
        open(unit=201,file=outputfile1)
        open(unit=202,file=outputfile2)
        print *
        print *, 'opening ', trim(adjustl(filename))
        do while(.true.)
                read(100,end=100) nowhalo
                nhalo = nhalo + 1
                if(nowhalo%mass < masscutmin) cycle
                if(nowhalo%mass > masscut1) nmasscut1 = nmasscut1 + 1
                
                x0=nowhalo%x; y0=nowhalo%y; z0=nowhalo%z;

                do j = 1, 8
                        if(j.eq.1) then; x=x0;y=y0;z=z0; endif
                        if(j.eq.2) then; x=x0-3150;y=y0;z=z0; endif
                        if(j.eq.3) then; x=x0;y=y0-3150;z=z0; endif
                        if(j.eq.4) then; x=x0;y=y0;z=z0-3150; endif
                        if(j.eq.5) then; x=x0;y=y0-3150;z=z0-3150; endif
                        if(j.eq.6) then; x=x0-3150;y=y0;z=z0-3150; endif
                        if(j.eq.7) then; x=x0-3150;y=y0-3150;z=z0; endif
                        if(j.eq.8) then; x=x0-3150;y=y0-3150;z=z0-3150; endif

                        if(x*x + y*y + z*z > rmaxsq) cycle

                        if(nowhalo%mass > masscut1) then
                                nwrite1 = nwrite1 +1
                                write(201, '(7e15.7)') x,y,z, nowhalo%vx, nowhalo%vy, nowhalo%vz, nowhalo%mass
                        endif
                        if(nowhalo%mass > masscut2) then
                                nwrite2 = nwrite2 +1
                                write(202, '(7e15.7)') x,y,z, nowhalo%vx, nowhalo%vy, nowhalo%vz, nowhalo%mass
                        endif
                enddo
                cycle
100             exit        
        enddo

        close(100)
        close(201)
        close(202)
        print *, 'Finishing processing ', nhalo, 'halos'
        print *, nmasscut1, 'halos pass masscut1'
        print *, nwrite1, 'halos written to ', trim(adjustl(outputfile1))
        print *, nwrite2, 'halos written to ', trim(adjustl(outputfile2))

        

end program main
