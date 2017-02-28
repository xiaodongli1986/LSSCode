
! 1. Plans

! 1. The file system: can find all files;
! 2. Loadin & ximu: can load in 2pcf files, and compute ximu
! 3. Cosmo convert: convert 2pcf from 1 to another
! 4. Covmat: compute covmat from all covmat mocks
! 5. Chisq: load in files, compute chisqs!



module LSS_ximu_tools

!use LSS_cosmo_funs
implicit none

  integer, parameter :: rt = kind(1.0d0)
  integer, parameter :: charlen = 300
  type :: par 
    real(rt) :: omegam, w
  end type

  character(len=charlen), parameter :: filedir = '/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/'
 
!----------------------------
! 1. 

contains


  real(rt) function DAofz(nowpar)
    type(par) :: nowpar
  end function DAofz
  real(rt) function Hofz(nowpar)
    type(par) :: nowpar
  end function Hofz



  character(25) function omwstr(omegam, w)
    real(rt) :: omegam, w
    character(10) :: str1, str2
    write(str1, '(f8.4)') omegam
    write(str2, '(f8.4)') w
    omwstr = 'om'//trim(adjustl(str1))//'_w'//trim(adjustl(str2))
    !omstr = str3
  end function omwstr

!1. Name of 2-point correlation function (2pcf) files 

!-----------------------------------
! data (baseline)    
  character(len=charlen) function data2pcffile_base(i_redshiftbin)
    integer :: i_redshiftbin
    character(len=charlen) :: filename
    if (i_redshiftbin .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    endif
    data2pcffile_base = filename
  end function data2pcffile_base

!-----------------------------------
! data (in general, i.e. in different cosmologies)
  character(len=charlen) function data2pcffile_general(i_redshiftbin, nowpar)
    ! args
    integer :: i_redshiftbin
    type(par) :: nowpar
    ! local
    character(20) :: nowomwstr
    character(len=charlen) :: filename

    nowomwstr = omwstr(nowpar%omegam, nowpar%w)
    if (i_redshiftbin .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    endif
    data2pcffile_general = filename
  end function data2pcffile_general

!-----------------------------------
! mock for systematic correction
  character(len=charlen) function syscor2pcffile(i_redshiftbin, imock)
    integer :: i_redshiftbin, imock
    character(len=charlen) :: filename
    character(len=15) :: str1
     
    write(str1, '(i3)') imock-1 
    str1 = 'J08.RSD.00'//trim(adjustl(str1))

    if (i_redshiftbin .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    endif
    syscor2pcffile = filename
  end function syscor2pcffile

!-----------------------------------
! mock for systematic correction
  character(len=charlen) function covmat2pcffile(i_redshiftbin, imock)
    integer :: i_redshiftbin, imock
    character(len=charlen) :: filename
    character(30) :: str1
     
    write(str1, '(i5)') imock-1 
    if (imock .le. 10) then
       str1 = 'PatchyV6C.RSD.000'//trim(adjustl(str1))
    elseif (imock .le. 100) then
       str1 = 'PatchyV6C.RSD.00'//trim(adjustl(str1))
    elseif (imock .le. 1000) then
       str1 = 'PatchyV6C.RSD.0'//trim(adjustl(str1))
    elseif (imock .le. 10000) then
       str1 = 'PatchyV6C.RSD.'//trim(adjustl(str1))
    endif

    if (i_redshiftbin .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (i_redshiftbin .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    endif
    covmat2pcffile = filename
  end function covmat2pcffile


  subroutine ximu_loadsmufile(filename, smutab, smax, nummubin)
    ! argument
    character(len=charlen), intent(in) :: filename
    integer, intent(in) :: smax, nummubin
    real(rt), intent(out) :: smutab(smax,nummubin,4) ! Four values at each s/mu: DD, DR, RR, xi

    smutab = 0.0_rt
  end subroutine

    

end module LSS_ximu_tools




program main_LSS_ximu

use LSS_ximu_tools

implicit none

    character(len=charlen) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, nowfile
    integer :: i, i_redshiftbin, imock
    logical :: logvar
    type(par) :: nowpar

    

    printstr = "Now it is empty!"

    print *, '  (LSS_ximu) checking existence of data base files...'
    do i_redshiftbin = 1, 6
        nowfile=data2pcffile_base(i_redshiftbin)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, i_redshiftbin, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
    enddo

    nowpar%omegam=0.26; nowpar%w=0.0
    do i_redshiftbin = 1, 6
        nowfile=data2pcffile_general(i_redshiftbin,nowpar)
        inquire(file=nowfile, exist=logvar)
        if (.not.logvar) then 
            print *, i_redshiftbin, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for systematic correction)...'
    do i_redshiftbin = 1, 6
     do imock = 1, 4
        nowfile=syscor2pcffile(i_redshiftbin,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, i_redshiftbin, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
     enddo
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for covmat estimation)...'
    do i_redshiftbin = 1, 6
     do imock = 1, 2000
        nowfile=covmat2pcffile(i_redshiftbin,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, i_redshiftbin, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
     enddo
    enddo


    if(iargc().le.1) then
        print *, printstr
        stop
    endif

    outputfile = ""
    do i = 1, iargc()
        if(mod(i,2).eq.0) cycle
        call getarg(i,tmpstr1)
        call getarg(i+1,tmpstr2)
        if(trim(adjustl(tmpstr1)).eq."-inputfile") then
            read(tmpstr2,"(A)") inputfile
        elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
            read(tmpstr2,"(A)") outputfile
        else
            print *, "Unkown argument: ", trim(adjustl(tmpstr1))
            write(*,"(A)") trim(adjustl(printstr))
            stop
        endif
    enddo

    if(trim(adjustl(outputfile)).eq."") then
        
    endif

    print *, 'This is an empty program LSS_ximu!!'

end program
