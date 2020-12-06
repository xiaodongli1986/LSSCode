

program main

implicit none


        integer*8 :: ID, CentralID(15), fofhaloid
        double precision :: x,y,z
        real :: vx,vy,vz
        real :: subhalomass, galmass
        integer :: idummy


character(*), parameter :: dir = '/home/eostm/galaxyCatalog/fullCatalog/2LPT2048p1024s/'
character(len=100) :: cosmostrs(5), cosmostrs2(5)
character(len=500) :: nowfile, nowoutfile
character(len=4) :: snapIDstrs(5,4)
character(len=3) :: zchars(4) = (/ '2.0', '1.0', '0.5', '0.0' /) 

integer:: icosmo, iss, igal

cosmostrs(1)  = 'Om0.26w-1';
cosmostrs2(1) = 'om0.2600_w-1.0000';
cosmostrs(2)  = 'Om0.26w-0.5';
cosmostrs2(2) = 'om0.2600_w-0.5000';
cosmostrs(3)  = 'Om0.26w-1.5';
cosmostrs2(3) = 'om0.2600_w-1.5000';
cosmostrs(4)  = 'Om0.21w-1';
cosmostrs2(4) = 'om0.2100_w-1.0000';
cosmostrs(5)  = 'Om0.31w-1';
cosmostrs2(5) = 'om0.3100_w-1.0000';

!Om0.21w-1/data:
!J08.new.noNorm.0031  J08.new.noNorm.0059  J08.new.noNorm.0088
!J08.new.noNorm.0139

!Om0.26w-0.5/data:
!J08.new.noNorm.0026  J08.new.noNorm.0048  J08.new.noNorm.0069
!J08.new.noNorm.0115

!Om0.26w-1.5/data:
!J08.new.noNorm.0030  J08.new.noNorm.0058  J08.new.noNorm.0088
!J08.new.noNorm.0142

!Om0.26w-1/data:
!J08.new.noNorm.0028  J08.new.noNorm.0054  J08.new.noNorm.0080
!J08.new.noNorm.0131

!Om0.31w-1/data:
!J08.new.noNorm.0027  J08.new.noNorm.0051  J08.new.noNorm.0077
!J08.new.noNorm.0126


snapIDstrs(1,:) = (/'0028', '0054', '0080', '0131' /)
snapIDstrs(2,:) = (/'0026', '0048', '0069', '0115' /)
snapIDstrs(3,:) = (/'0030', '0058', '0088', '0142' /)
snapIDstrs(4,:) = (/'0031', '0059', '0088', '0139' /)
snapIDstrs(5,:) = (/'0027', '0051', '0077', '0126' /)

!xmin=1.0e30;ymin=xmin; zmin=xmin;
!xmax=-xmin; ymax=xmax; zmax=xmax;

do icosmo = 1, 5
 do iss = 1, 4
   igal = 0
   nowfile = trim(adjustl(dir))//trim(adjustl(cosmostrs(icosmo)))//'/data/J08.new.noNorm.'//trim(adjustl(snapIDstrs(icosmo,iss)))
   nowoutfile = trim(adjustl(cosmostrs2(icosmo)))//'_red'//zchars(iss)
   write(*,'(A)') 'Opening  ', trim(adjustl(nowfile))
   write(*,'(A)') 'Output to', trim(adjustl(nowoutfile))
   call system('ls -alh '//trim(adjustl(nowfile)))
   open(100987,file=nowfile,form='binary',action='read')
   open(100988,file=nowoutfile)
   do while(.true.)
     if(icosmo.eq.4) then
           read(100987,end=10) ID, CentralID(1:15), fofhaloid, x,y,z, vx,vy,vz,  subhalomass, galmass, idummy
     else
           read(100987,end=10) ID, CentralID(1:10), fofhaloid, x,y,z, vx,vy,vz,  subhalomass, galmass, idummy
     endif
     write(100988,'(7e17.9)') x,y,z, vx,vy,vz, galmass
     igal = igal + 1
     cycle
10   exit
   enddo
   close(100987); close(100988)
   print *, 'Finishing reading / writing ', igal, 'galaxies.'
 enddo
enddo



end program
