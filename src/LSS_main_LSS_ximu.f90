
! 1. Plans

! 1. The file system: can find all files;
! 2. Loadin & ximu: can load in 2pcf files, and compute ximu
! 3. Cosmo convert: convert 2pcf from 1 to another
! 4. Covmat: compute covmat from all covmat mocks
! 5. Chisq: load in files, compute chisqs!


!###############################################################
!## Module: Types, constants, tools
module types_constants
implicit none

  integer, parameter  :: rt =  kind(1.0d0)
  integer, parameter  :: charlen=1000
  real(rt), parameter :: CONST_C = 299792.458 !unit: km/s
  
  
  ! Setting: The binning schemes (number of bins in the mu space )
  !  integer, parameter :: N1=1,N2=1,mubins(N1)=(/10/); real(rt),parameter::mucuts(N2)=(/0.99_rt/)
  integer, parameter :: N1=36
  integer, parameter :: mubins(N1) = (/ 5,6,7,8,9,10,11,12,13,14,&
  		15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30, &
  		31,32,33,34,35,36,37,38,39,40 /)

! Setting: The binning schemes (maximal cuts on mu )
  integer, parameter :: N2=15
  real(rt), parameter :: mucuts(N2) = (/ 0.99_rt, 0.98_rt, 0.97_rt, 0.96_rt, 0.95_rt, 0.94_rt, &
     0.93_rt, 0.92_rt, 0.91_rt, 0.90_rt, 0.89_rt, 0.88_rt, 0.87_rt, 0.86_rt, 0.85_rt /)
! In total we have N1*N2 schemes with different (#-bin, mucut)

! Setting: Integration limits of s
  real(rt), parameter :: ints1 = 6.0_rt, ints2 = 40.0_rt

! Setting: directory for covmat files, chisq files
  character(len=charlen), parameter :: covmatdir = '/home/xiaodongli/LSS/2PCF_AP/covmats/', &
    chisqdir = '/home/xiaodongli/LSS/2PCF_AP/chisqs/'
  
! Setting: directory of the 2pCF (2-point correlation funtion) files 
  character(len=charlen), parameter :: filedir = '/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/'

! Settings: info of the 2pCF files. Including: maximal s (rmax); number of bins in s and mu (nbins, mubins)
!  integer, parameter :: smax_database=150, nbins_database=750, mubins_database=600, & ! 2pCF of observational data, in "baseline" cosmologies
!  integer, parameter :: smax_database=150, nbins_database=1200, mubins_database=960, & ! 2pCF of observational data, in "baseline" cosmologies
  integer, parameter :: smax_database=150, nbins_database=600, mubins_database=480, & ! 2pCF of observational data, in "baseline" cosmologies
			smax_data=51, nbins_data=51, mubins_data=120, & ! 2pCF of observational data, in general cosmologies
                        smax_sysmock=150, nbins_sysmock=150, mubins_sysmock=120, & ! 2pCF from Horizon Run 4 mocks, used for estimation of covmat 
                        smax_covmock=51,  nbins_covmock=51, mubins_covmock=120, &  ! 2pCF from Patchy mocks, used for systematic correction
                        ncovmocks = 2000, nsysmocks = 4 ! number of mocks used for covmat-estimation and systematic-correction
!                        ncovmocks = 50, nsysmocks = 4 ! number of mocks used for covmat-estimation and systematic-correction

! Settings: polynomical fitting degeree for dintxi_sys
  integer, parameter :: polyfitdeg = -1
                        
! Settings: effectove redshifts of the redshift bins
  integer, parameter :: nz = 6                        
  real(rt), parameter :: zeffs(nz) = (/ 0.2154098242_rt, 0.3156545036_rt, 0.3869159269_rt, &
                                        0.4785988874_rt, 0.5408209467_rt, 0.6186561000_rt  /)
!  integer, parameter :: nz = 2
!  real(rt), parameter :: zeffs(nz) = (/ 0.2154098242_rt, 0.3156545036_rt /)

contains
  subroutine nizhen(aa,b,n)
    ! Arguments
    real(rt), intent(in) :: aa(n,n)
    integer, intent(in) :: n
    real(rt), intent(out) :: b(n,n)
    ! Local
    integer :: i,j,k
    real(rt) :: a(n,n)
    a=aa
    b=0.0_rt
    do i=1,n
      b(i,i)=1
    enddo
    do i=1,n
      b(i,:)=b(i,:)/a(i,i)
      a(i,i:n)=a(i,i:n)/a(i,i) 
      do j=i+1,n     
        do k=1,n   
          b(j,k)=b(j,k)-b(i,k)*a(j,i)
        enddo
        a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
      enddo
    enddo
    do i=n,1,-1
      do j=i-1,1,-1
        do k=1,n
          b(j,k)=b(j,k)-b(i,k)*a(j,i)
        enddo
      enddo
    enddo
  end subroutine nizhen
  !------------------------------------------
  ! Polynomial Regression:
  !  polynomial fitting to data points
  ! Y(i) = A(1) + A(2)*X(1) + A(3)*X(2)^2 + ... 
  !        + A(n+1)*X(n)^n
  !------------------------------------------  	
	subroutine poly_fit(X,Y,A,ndat,n)
		! Dummy
		real(rt), intent(in) :: X(ndat),Y(ndat)
		integer, intent(in) :: ndat,n
		real(rt), intent(out) :: A(n+1)
		! Local
		real(rt) :: CapX(ndat,n+1), CapMatA(n+1,n+1), CapMatB(n+1,n+1), CapMatC(n+1,ndat)
		integer :: i,j,k
		! CapX(i,j) = X(i) ^ (j-1)
		do i = 1, ndat
		do j = 1, n+1
			CapX(i,j) = X(i)**(j-1)
		enddo
		enddo
		! CapMatA = ( CapX^T CapX)
		do i = 1,n+1
		do j = 1,n+1
			CapMatA(i,j) = 0.0_rt
			do k = 1,ndat
				CapMatA(i,j) = CapMatA(i,j) + CapX(k,i)*CapX(k,j)
			enddo
		enddo
		enddo
		! CapMatB = ( CapX^T CapX )^(-1)
		call nizhen(CapMatA,CapMatB,n+1)
		! CapMatC = ( CapX^T CapX )^(-1) CapX^T
		do i = 1, n+1
		do j = 1, ndat
			CapMatC(i,j) = 0.0_rt
			do k = 1, n+1
				CapMatC(i,j) = CapMatC(i,j) + CapMatB(i,k)*CapX(j,k)
			enddo
		enddo
		enddo
		! A = ( CapX^T CapX )^(-1) CapX^T Y
		do i = 1, n+1
			A(i) = 0.0_rt
			do j = 1, ndat
				A(i) = A(i) + CapMatC(i,j)*Y(j)
			enddo
		enddo
  	end subroutine poly_fit
  !------------------------------------------
  ! Value of a n-th polynomial 
  !  at some value of x
  !------------------------------------------   	
  	real(rt) function poly(x,A,n)
  		! Dummy
  		real(rt), intent(in) :: x, A(n+1)
  		integer, intent(in) :: n
  		! local
  		integer :: i
  		poly = A(1)
  		do i = 1, n
  			poly = poly + A(i+1)*x**(i)
  		enddo
  	end function poly
  !------------------------------------------
  ! polynomial regression of curve Y
  !------------------------------------------
   function polyfitY(X,Y,n,polyfitdeg)
     integer, intent(in) :: n, polyfitdeg
     real(rt), intent(in) :: X(n), Y(n)
     real(rt) :: polyfitY(n), coeff(polyfitdeg+1)
     integer :: i
     call poly_fit(X,Y,coeff,n,polyfitdeg)
     !print *, coeff
     do i = 1, n
       polyfitY(i) = poly(X(i),coeff,polyfitdeg)
     enddo
   end function polyfitY

end module types_constants
!###############################################################

!###############################################################
!## Module: Cosmological functions (DA and H)
module cosmology

use types_constants

implicit none

  type :: par 
    real(rt) :: omegam, w
  end type


contains

!---------------------------------------------------------------
! Simpson integration for some cosmological function 
  real(rt) function CosmoSimpson(nowpar,fun,xleft,xright,N)
    real(rt), EXTERNAL :: fun
    real(rt), intent(in) :: xleft, xright
    integer, intent(in) :: N
    real(rt) :: x1,x2,BC,f1,f2
    integer :: i
    type(par) :: nowpar
    BC=(xright-xleft)/dble(N)
    x1=xleft;x2=x1+BC;
    f1=fun(nowpar,x1);f2=fun(nowpar,x2);CosmoSimpson=(f1+fun(nowpar,(x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0;
    do i = 2,N
      x1=x2;f1=f2;x2=x2+BC;
      f2=fun(nowpar,x2);CosmoSimpson=CosmoSimpson+(f1+fun(nowpar,(x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0
    enddo
  end function CosmoSimpson
!----------------------------------------------------------------
! inv_ez = 1/E(z); E(z) = H(z)/H0
  real(rt) function inv_ez(nowpar,z)
    type(par) :: nowpar
    real(rt) :: z
    inv_ez =  1. / dsqrt(dble(nowpar%omegam*(1.+z)**3.  &
      + (1.-nowpar%omegam)*(1.0+z)**(3.*(1.+nowpar%w))))
  end function inv_ez 
!----------------------------------------------------------------
! rz = \int 1/H(z) dz, in unit of Mpc/h
  real(rt) function rz(nowpar,z)
    type(par) :: nowpar
    real(rt), intent(in) :: z
    ! use at least 128 bins in redshift space; 128 bins for every redshift interval of 0.25
    rz = CosmoSimpson(nowpar,inv_ez,0.0_rt,z,N=max(128,128*ceiling(z/0.25_rt)))*CONST_C/100.d0
  end function rz
!----------------------------------------------------------------
! DA(z), in unit of Mpc/h
  real(rt) function DAofz(nowpar,z)
    type(par) :: nowpar
    real(rt), intent(in) :: z
    DAofz = rz(nowpar,z) / (1.0+z)
  end function DAofz
!----------------------------------------------------------------
! H(z), in unit of km/s/(Mpc/h)
  real(rt) function Hofz(nowpar,z)
    type(par) :: nowpar
    real(rt), intent(in) :: z
    Hofz=100./inv_ez(nowpar,z)
  end function Hofz
end module cosmology


module LSS_ximu_tools
use cosmology
implicit none
  


! Structure storing covmats
  type :: M2
    real(rt), allocatable :: A(:,:)
    integer :: nA=-1
  end type
  type(M2) :: covmats(N1,N2,nz-1)
  real(rt) :: dintxi_syscor(maxval(mubins),N1,N2,nz-1)
     
!----------------------------
! 1. 

contains

  character(25) function omwstr(omegam, w)
    real(rt) :: omegam, w
    character(10) :: str1, str2
    write(str1, '(f8.4)') omegam
    write(str2, '(f8.4)') w
    omwstr = 'om'//trim(adjustl(str1))//'_w'//trim(adjustl(str2))
  end function omwstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Settings: Name of 2pCF files
!-----------------------------------
! data (baseline)    
  character(len=charlen) function data2pcffile_base(iz, ombase, wbase)
    integer :: iz
    real(rt) :: ombase, wbase
    character(len=charlen) :: filename, filestr, tmpstr1,tmpstr2,tmpstr3
    write(tmpstr1,*) nbins_database; write(tmpstr2,*) mubins_database; write(tmpstr3,*) smax_database
    filestr = '.rmax'//trim(adjustl(tmpstr3))//'.'//trim(adjustl(tmpstr1))//'rbins.'&
      //trim(adjustl(tmpstr2))//'mubins.'
    if (iz .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    endif
    data2pcffile_base = filename
  end function data2pcffile_base
!-----------------------------------
! data (in general, i.e. in different cosmologies)
  character(len=charlen) function data2pcffile_general(iz, nowpar)
    ! args
    integer :: iz
    type(par) :: nowpar
    ! local
    character(20) :: nowomwstr
    character(len=charlen) :: filename
    nowomwstr = omwstr(nowpar%omegam, nowpar%w)
    if (iz .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    endif
    data2pcffile_general = filename
  end function data2pcffile_general
!-----------------------------------
! mock (for systematic correction)
  character(len=charlen) function syscor2pcffile(iz, imock)
    integer :: iz, imock
    character(len=charlen) :: filename
    character(len=15) :: str1
    write(str1, '(i3)') imock-1 
    str1 = 'J08.RSD.00'//trim(adjustl(str1))
    if (iz .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    endif
    syscor2pcffile = filename
  end function syscor2pcffile
!-----------------------------------
! mock for covariance matrix estimation
  character(len=charlen) function cov2pcffile(iz, imock)
    integer :: iz, imock
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
    if (iz .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    endif
    cov2pcffile = filename
  end function cov2pcffile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Load in in the xi(s,mu) table
!   All files are in format of : smin, smax, minimal value of 1-mu, maximal value of 1-mu, DD, DR, RR, three estimators of xi
!   First row in comment
  subroutine ximu_loadsmufile(filename, smutab, smax, nbins, mubins)
    ! argument
    character(len=charlen), intent(in) :: filename
    integer, intent(in) :: smax, nbins, mubins
    real(rt), intent(out) :: smutab(nbins,mubins,3) ! Four values at each s/mu: DD, DR, RR, xi
    ! variables
    character(len=charlen) :: nowstr
    integer :: i, j
    real :: tmpx

    smutab = 0.0_rt
    
    open(unit=44817,file=filename)
    print *, ' (ximu_loadsmufile) Load in 2pCF file: ', trim(adjustl(filename))
    read(44817,*) nowstr
    !print *, nowstr
    do i = 1, nbins
    do j = 1, mubins
    !    print *, i,j
        read(44817,*) tmpx, tmpx, tmpx, tmpx, smutab(i,j,1:3) !, tmpx, tmpx, smutab(i,j,4)
    enddo
    enddo
    close(44817)
  end subroutine ximu_loadsmufile

! Compute covmats
  subroutine calc_covmats()
    integer :: i,j,i1,i2,n, iz,imock, maxmubin
    real(rt), allocatable :: intxis(:,:,:,:,:), tmpX(:), dintxi(:,:)
    real(rt) :: smutab_covmock(nbins_covmock,mubins_covmock), deltas, tmpx1, tmpx2, tmpx3
    character(charlen) :: nowfile, tmpstr
    logical :: logvar
    ! 1. Initialize the structure "covmats"
    do i = 1, N1
    do j = 1, N2
    do iz= 2, nz
      if (allocated(covmats(i,j,iz-1)%A)) deallocate(covmats(i,j,iz-1)%A)
      n=mubins(i); covmats(i,j,iz-1)%nA = n-1
      allocate(covmats(i,j,iz-1)%A(n-1,n-1))
      covmats(i,j,iz-1)%A = 0.0_rt
    enddo
    enddo
    enddo
    ! 2. Compute intxi, for nz redshift bins, all 2000 mocks
    ! 2.1. Initialize the structure intxis
    maxmubin = maxval(mubins)
    print *, 'maxmubin = ', maxmubin
    allocate(intxis(maxmubin,ncovmocks,nz,N1,N2))
    allocate(dintxi(maxmubin,ncovmocks),tmpX(maxmubin))
    intxis = 0.0_rt; dintxi = 0.0_rt
    ! 2.2 Compute intxis
    do iz = 1, nz
      do imock = 1, ncovmocks
        ! 2.2.1 Check existence of 2pCF file
        nowfile=cov2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_covmats) file not found: ', trim(adjustl(nowfile)); stop
        endif
        ! 2.2.2 Load in 2pCF file
        call ximu_loadsmufile(nowfile, smutab_covmock, smax_covmock, nbins_covmock, mubins_covmock)
        ! 2.2.3 Compute intxi in the N1*N2 binning schemes
        deltas = smax_covmock / dble(nbins_covmock)
        do i = 1, N1
        do j = 1, N2
          call smuintxi(smutab_covmock, deltas, nbins_covmock, mubins_covmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)))
            call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
        enddo
        enddo
      enddo
    enddo
    ! 2.3 Compute covariance matrices
    ! 2.3.1 Loop of schemes & redshifts
    do i = 1, N1
    do j = 1, N2
    do iz= 2, nz
      print *, 'mubin, mucut, iz = ', mubins(i), real(mucuts(j)), iz
      ! 2.3.2 dintxi = intxi@high_z - intxi@low_z; compute it for all covmocks
      n = mubins(i)
      do imock= 1, ncovmocks
       dintxi(1:n-1,imock) = intxis(1:n-1,imock,iz,i,j)-intxis(1:n-1,imock,1,i,j)
      enddo
      ! 2.3.3 compute covmat using dintxi
      covmats(i,j,iz-1)%A = 0.0_rt
      do i1= 1,n-1
      do i2 = i1,n-1
        ! mean(X*Y), mean(X), mean(Y)
        tmpx1 = 0.0_rt; tmpx2 = 0.0_rt; tmpx3 = 0.0_rt
        do imock = 1, ncovmocks
          tmpx1 = tmpx1 + dintxi(i1,imock)*dintxi(i2,imock)
          tmpx2 = tmpx2 + dintxi(i1,imock)
          tmpx3 = tmpx3 + dintxi(i2,imock)
        enddo
        tmpx1 = tmpx1 / dble(ncovmocks); tmpx2 = tmpx2 / dble(ncovmocks); tmpx3 = tmpx3 / dble(ncovmocks)
        ! covariance = mean(X*Y) - mean(X)*mean(Y)
        covmats(i,j,iz-1)%A(i1,i2) = (tmpx1 - tmpx2*tmpx3) * dble(ncovmocks) / dble(ncovmocks-1)
      enddo
      enddo
      ! symmetric matrix
      do i1=1,n-1
      do i2=1,i1-1
        covmats(i,j,iz-1)%A(i1,i2) = covmats(i,j,iz-1)%A(i2,i1)
      enddo 
      enddo
      if(i.eq.1 .and. j.eq.1) print *, covmats(i,j,iz-1)%A(1:n-1,1:n-1)
    enddo
    enddo
    enddo
  end subroutine calc_covmats
! Name of files storing covariance matrix
  character(charlen) function covmatfilename(mubin,mucut,iz,suffixstr)
    integer, intent(in) :: mubin, iz
    real(rt), intent(in) :: mucut
    character(*), intent(in), optional :: suffixstr  
    character(charlen) :: tmpstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    write(tmpstr1,*) mubin
    write(tmpstr2,'(f5.2)') mucut
    write(tmpstr3,*) iz
    write(tmpstr4,*) ncovmocks
    write(tmpstr5,'(f5.1)') ints1
    write(tmpstr6,'(f5.1)') ints2
    covmatfilename = trim(adjustl(covmatdir))//''//trim(adjustl(tmpstr1)) &
        //'mubins.mumax' // trim(adjustl(tmpstr2))//'.iz'//trim(adjustl(tmpstr3)) &
        //'.CovMock_'//trim(adjustl(tmpstr4)) &
        //'.s'//trim(adjustl(tmpstr5))//'to'//trim(adjustl(tmpstr6))//'.covmat'
  end function covmatfilename
! Output covmats to files
  subroutine output_covmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,iz,k,l,nA
    character(len=charlen) :: nowfile, tmpstr, tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
    do iz= 2, nz
      nowfile = covmatfilename(mubins(i),mucuts(j),iz,tmpstr)
      print *, 'Write covmat to : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile)
      nA = covmats(i,j,iz-1)%nA
      do k = 1,nA
       write(tmpstr1, '(e15.7)') covmats(i,j,iz-1)%A(k,1)
       do l = 2,nA
         write(tmpstr2,'(e15.7)') covmats(i,j,iz-1)%A(k,l)
         tmpstr1 = trim(adjustl(tmpstr1))//' '//tmpstr2
       enddo
       write(7733,'(A)') trim(adjustl(tmpstr1))
      enddo
      close(7733)
    enddo
    enddo
    enddo
  end subroutine output_covmats
! Load covmats from files
  subroutine load_covmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,iz,k,l,nA
    character(len=charlen) :: nowfile, tmpstr, tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
    do iz= 2, nz
      nowfile = covmatfilename(mubins(i),mucuts(j),iz,tmpstr)
!      print *, 'Write covmat to : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile,action='read')
      nA = mubins(i)-1
      if (allocated(covmats(i,j,iz-1)%A) .and. covmats(i,j,iz-1)%nA .ne. nA) deallocate(covmats(i,j,iz-1)%A)
      if (.not. allocated(covmats(i,j,iz-1)%A)) allocate(covmats(i,j,iz-1)%A(nA,nA))
      covmats(i,j,iz-1)%nA = nA
      do k = 1, nA
       read(7733,*) covmats(i,j,iz-1)%A(:,k)
      enddo
      close(7733)
    enddo
    enddo
    enddo
  end subroutine load_covmats
! Invert all covmats
  subroutine invert_covmats()
    integer :: i,j,iz, n
    real(rt), allocatable :: B(:,:)
    do i = 1, n1
      allocate(B(mubins(i)-1,mubins(i)-1))
      do j = 1, n2
        do iz= 2, nz
          call nizhen(covmats(i,j,iz-1)%A, B, mubins(i)-1)
          covmats(i,j,iz-1)%A = B
        enddo
      enddo
      deallocate(B)
    enddo
  end subroutine invert_covmats

  subroutine calc_syscor()
    integer :: i,j,i1,i2,n, iz,imock, maxmubin
    real(rt) :: smutab_sysmock(nbins_sysmock,mubins_sysmock), &
      intxis(maxval(mubins),nsysmocks,nz,N1,N2),&
      intxi0(maxval(mubins)), intxi1(maxval(mubins)),&
      tmpX(maxval(mubins)),X(maxval(mubins)),mumids(maxval(mubins),N1,N2),&
      deltas, tmpx1, tmpx2, tmpx3
    character(charlen) :: nowfile, tmpstr
    logical :: logvar
    intxis = 0.0_rt
    ! 1. Compute intxi, for nz z-bins, all syscor mocks, N1*N2 schemes
    do iz = 1, nz
      do imock = 1, nsysmocks
        ! 1.1 Check existence of 2pCF file
        nowfile=syscor2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_covmats) file not found: ', trim(adjustl(nowfile)); stop
        endif
        ! 1.2 Load in 2pCF file
        call ximu_loadsmufile(nowfile, smutab_sysmock, smax_sysmock, nbins_sysmock, mubins_sysmock)
        ! 1.3 Compute intxi in the N1*N2 binning schemes
        deltas = smax_sysmock / dble(nbins_sysmock)
        do i = 1, N1
        do j = 1, N2
          call smuintxi(smutab_sysmock, deltas, nbins_sysmock, mubins_sysmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)),&
            mumids=mumids(1:mubins(i),i,j))
          call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
!          print *, real(mumids(1:mubins(i)))
        enddo
        enddo
!        stop
      enddo
    enddo

    ! 2. Redshift evolution of intxi, for nz-1 z-bins, all syscor mocks, N1*N2 schemes
    do i = 1, N1
    do j = 1, N2
      ! 2.1 reference intxi at first bin
      intxi0 = 0.0_rt
      do imock = 1, nsysmocks
        intxi0(1:mubins(i)-1) = intxi0(1:mubins(i)-1) + intxis(1:mubins(i)-1,imock,1,i,j)
      enddo
      intxi0 = intxi0 / dble(nsysmocks)
      ! 2.2 intxi at a higher redshift bin
      do iz = 2, nz
        intxi1 = 0.0_rt
        do imock = 1, nsysmocks
          intxi1(1:mubins(i)-1) = intxi1(1:mubins(i)-1) + intxis(1:mubins(i)-1,imock,iz,i,j)
        enddo
        intxi1 = intxi1 / dble(nsysmocks)
        ! 2.3 Redshift evolution
        dintxi_syscor(1:mubins(i)-1,i,j,iz-1) = &
          intxi1(1:mubins(i)-1) - intxi0(1:mubins(i)-1)
      enddo
    enddo
    enddo
    
    ! polynomial fit 
    if (polyfitdeg .ge. 1) then
      do i = 1, maxval(mubins)
        X(i) = i-1
      enddo
      do i = 1, N1
        do j = 1, N2
          do iz = 2, nz
            dintxi_syscor(1:mubins(i)-1,i,j,iz-1) = &
              polyfitY(mumids(1:mubins(i)-1,i,j),dintxi_syscor(1:mubins(i)-1,i,j,iz-1),mubins(i)-1,polyfitdeg) ! using exact values of mumids for polyfit; almost no effect on the resulted chisq values
!              polyfitY(X(1:mubins(i)-1),dintxi_syscor(1:mubins(i)-1,i,j,iz-1),mubins(i)-1,polyfitdeg)
          enddo
        enddo
      enddo
   endif
  end subroutine calc_syscor
  
! Coordiante transformation of (s,mu) : from one cosmology to another
  subroutine smu__CosmoConvert(s,mu,DA1,DA2,H1,H2,s_prime,mu_prime)
    real(rt), intent(in) :: s,mu,DA1,DA2,H1,H2
    real(rt), intent(out) :: s_prime,mu_prime
    real(rt) :: alpha1, alpha2, s1, s2  ! s1: angular direction; s2: LOS direction 
    s2 = s*mu;
    s1 = sqrt(s*s - s2*s2)
    alpha1 = DA2 / DA1
    alpha2 = H1 / H2
    s_prime  =  sqrt((alpha1*s1)**2 + (alpha2*s2)**2)
    mu_prime =  alpha2*s2 / s_prime
  end subroutine smu__CosmoConvert


! Mapping xi(s,mu) from baseline cosmology to another 
! Dense grid in baseline cosmology, sparse grid in another
  subroutine DSMapping(smutabstd, nums1, nummu1, smutab2, &
	  nums2, nummu2, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  &
	  smin_mapping, smax_mapping ) ! range of s considered in the coordinate transformation:  
	                               ! basically, smin_mapping < s1 < s2 < smax_mapping.
    ! argument
    real(rt), intent(in) :: smutabstd(nums1,nummu1,3), DAstd, DAnew, Hstd, Hnew, &
      deltas1, deltas2, smin_mapping, smax_mapping
    integer, intent(in) :: nums1, nummu1, nums2, nummu2
    real(rt), intent(out):: smutab2(nums2,nummu2,3)
    ! variables
    real(rt) :: deltamu1,deltamu2, mubound1,mubound2,sbound1,sbound2, &
      scenter,anglecenter,mucenter,scenter2,anglecenter2,mucenter2, &
      scenter2_bound,anglecenter2_bound,anglecenter3,scenter3, countDD,countDR,countRR,&
      s,angle,ds,dangle,d1,rat,rat_bound,rats,ratangle,rat1,rat2,rat3,rat4
    real(rt) :: smutabstd_centers(nums1,nummu1,2)
    integer :: maxs1, mins1, is1,iangle1,is2,iangle2, is_nearby,imu_nearby, &
      is2_bound, is2_b, iangle2_b, iangle2_bound, &
      is1_bound, iangle1_bound, iangle1_nearby, is1_nearby, &
      smutabstd_ismus(nums1,nummu1,2)
    logical :: sboundflag, muboundflag
    
    deltamu1 = 1.0d0/real(nummu1); deltamu2 = 1.0d0/real(nummu2)
    smutab2=0.0d0; smutabstd_ismus=10000
    ! range of s for smutab1, smutab2
    call smu__CosmoConvert(smax_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smax_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    maxs1 = floor( max(sbound1,sbound2) / deltas1 + 0.5 )
    if(maxs1>nums1) print *, ' WARNING (mapping_smudata_to_another_cosmology_DenseToSparse)!'&
      ' Outflow of s: nums1, maxs1 = ', nums1, maxs1
    call smu__CosmoConvert(smin_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smin_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    mins1 = floor( min(sbound1,sbound2) / deltas1 + 0.5 )
    mins1=max(mins1,0); maxs1=min(maxs1,nums1-1)
    
!    print *, 'deltas1, deltamu1 ', deltas1, deltamu1
!    print *, 'deltas2, deltamu2 ', deltas2, deltamu2
!    print *, 'sbound1, sbound2  ', sbound1, sbound2
    !print *, mins1, maxs1
    
!    maxs2 = min(floor(smax_mapping / deltas2 + 0.5), nums2)
  !      deltas1=0.2, deltamu1=1.0/600.0, deltas2=1.0, deltamu2=1.0/120.0, smin_mapping=1,smax_mapping=51,

!    ''' method can be simple_bin or divided_pixel'''
!    nummu1 = int(1.0/deltamu1 + 0.5)
!    sbound1, mubound1 = smu__CosmoConvert(smax_mapping,0,DAnew,DAstd,Hnew,Hstd)
!    sbound2, mubound2 = smu__CosmoConvert(smax_mapping,1,DAnew,DAstd,Hnew,Hstd)
!    sbound = max(sbound1, sbound2); nums1 = int(sbound / deltas1 + 0.5)
!    
!    sbound1, mubound1 = smu__CosmoConvert(smin_mapping,0,DAnew,DAstd,Hnew,Hstd)
!    sbound2, mubound2 = smu__CosmoConvert(smin_mapping,1,DAnew,DAstd,Hnew,Hstd)
!    sbound = min(sbound1, sbound2); mins1 = int(sbound / deltas1)
!    
!    nums2, nummu2 = int(smax_mapping / deltas2 + 0.5), int(1.0/deltamu2 + 0.5)
!    numrow3=len(smutabstd[0][0])
!    smutab2 = [[[0 for row3 in range(numrow3+1)] for row2 in range(nummu2)] for row1 in range(nums2)]

!    open(unit=9231489,file='smutab_database.txt')! haha

    do is1 = mins1, maxs1-1
      do iangle1 = 0, nummu1-1
        scenter=(is1+0.5_rt)*deltas1; anglecenter=(iangle1+0.5)*deltamu1
        mucenter=1.0_rt-anglecenter
        call smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,scenter2,mucenter2)
        anglecenter2 = 1.0_rt - mucenter2
        is2 = floor(scenter2 / deltas2); iangle2 = floor(anglecenter2/deltamu2)
        smutabstd_centers(is1+1,iangle1+1,1) = scenter2
        smutabstd_centers(is1+1,iangle1+1,2) = anglecenter2
        smutabstd_ismus(is1+1,iangle1+1,1) = is2
        smutabstd_ismus(is1+1,iangle1+1,2) = iangle2
!        write(9231489,'(i4,i4,3(e14.7))') is1, iangle1, smutabstd(is1+1,iangle1+1,1:3)! haha
      enddo
    enddo
!    close(9231489) ! haha
    
    
    !print *, 'Initialization of tabs done!'

!    
!    #print mins1, nums1, nummu1, nums2, nummu2
!    if method == 'divided_pixel':
!        smutabstd_centers = [[0 for row2 in range(nummu1)] for row1 in range(nums1)]
!    for is1 in range(mins1, nums1):
!        for iangle1 in range(nummu1):
!            scenter, anglecenter = (is1+0.5)*deltas1, (iangle1+0.5)*deltamu1
!            mucenter = 1.0 - anglecenter
!            scenter2, mucenter2 = smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,)
!            anglecenter2 = 1.0 - mucenter2
!            is2 = int(scenter2  / deltas2  )
!            iangle2 = int(anglecenter2 / deltamu2)
!            if method == 'simple_bin':
!                if is2 < nums2 and iangle2 < nummu2:
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
!                    if save_counts_row != None:
!                        smutab2[is2][iangle2][save_counts_row] += 1
!            elif method == 'divided_pixel':
!                smutabstd_centers[is1][iangle1] = [scenter2, anglecenter2, is2, iangle2]
!                
    do is1 = mins1, maxs1-1
      do iangle1 = 0, nummu1-1
        scenter2=smutabstd_centers(is1+1,iangle1+1,1)
        anglecenter2=smutabstd_centers(is1+1,iangle1+1,2)
	is2=smutabstd_ismus(is1+1,iangle1+1,1)
	iangle2=smutabstd_ismus(is1+1,iangle1+1,2)
        if (is2.ge.nums2 .or. iangle2.ge.nummu2) cycle
        
 !    if method == 'divided_pixel':
!        for is1 in range(mins1, nums1):
!            for iangle1 in range(nummu1):
!                scenter2, anglecenter2, is2, iangle2 = smutabstd_centers[is1][iangle1]
!                if not (is2 < nums2 and iangle2 < nummu2):
!                    continue       
        
        sboundflag=.false.; muboundflag=.false.
        do is1_nearby = max(is1-1,mins1), min(is1+1,maxs1-1)
          is2_b = smutabstd_ismus(is1_nearby+1,iangle1+1,1)
          iangle2_b = smutabstd_ismus(is1_nearby+1,iangle1+1,2)
          if (is2_b .ne. is2) then 
            sboundflag=.true.; is1_bound=is1_nearby; is2_bound=is2_b;
            scenter2_bound=smutabstd_centers(is1_nearby+1,iangle1+1,1)
          endif
        enddo
        do iangle1_nearby = max(iangle1-1,0), min(iangle1+1,nummu1-1)        
          is2_b = smutabstd_ismus(is1+1,iangle1_nearby+1,1)
          iangle2_b = smutabstd_ismus(is1+1,iangle1_nearby+1,2)
          if (iangle2_b .ne. iangle2) then 
            muboundflag=.true.; iangle1_bound=iangle1_nearby; iangle2_bound=iangle2_b;
            anglecenter2_bound=smutabstd_centers(is1+1,iangle1_nearby+1,2)
          endif
        enddo

            !seriously possible bug found in the python code! Decide to stop and checking the python code for a while.
          
            
!                ### firstly, check boundary:
!                sboundflag, muboundflag = False, False
!                for is1_nearby in [max(is1-1,mins1), min(is1+1,nums1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1_nearby][iangle1][2], smutabstd_centers[is1_nearby][iangle1][3]
!                    if is2_b != is2:
!                        sboundflag=True; is1_bound = is1_nearby; is2_bound = is2_b; 
!                        scenter2_bound = smutabstd_centers[is1_nearby][iangle1][0]
!                for iangle1_nearby in [max(iangle1-1,0), min(iangle1+1,nummu1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1][iangle1_nearby][2], smutabstd_centers[is1][iangle1_nearby][3]
!                    if iangle2_b != iangle2:
!                        muboundflag=True; iangle1_bound = iangle1_nearby; iangle2_bound = iangle2_b
!                        anglecenter2_bound = smutabstd_centers[is1][iangle1_nearby][1]
         
!                    
          countDD = smutabstd(is1+1,iangle1+1,1); 
          countDR = smutabstd(is1+1,iangle1+1,2); 
          countRR = smutabstd(is1+1,iangle1+1,3)
          if ( .not.sboundflag .and. .not. muboundflag) then 
            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR
          endif

!                
!                ### Then, treat them case by case...
!                ## s, mu are all not near the boundary of tab 2
!                if ((not sboundflag)and(not muboundflag)):
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += 1


          if ( sboundflag .and. .not. muboundflag) then 
            s = (is2 + is2_bound + 1) * 0.5 * deltas2
            scenter3 = (scenter2+scenter2_bound) / 2.0
            ds = scenter3 -scenter2
            d1 = s-scenter2+ds
            rat = min(d1 / (2.0_rt*ds),1.0_rt)
            rat_bound = 1.0-rat

            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat
            
            if (is2_bound < nums2) then
             smutab2(is2_bound+1,iangle2+1,1) = smutab2(is2_bound+1,iangle2+1,1)+countDD*rat_bound
             smutab2(is2_bound+1,iangle2+1,2) = smutab2(is2_bound+1,iangle2+1,2)+countDR*rat_bound
             smutab2(is2_bound+1,iangle2+1,3) = smutab2(is2_bound+1,iangle2+1,3)+countRR*rat_bound
            endif
          endif

!                            
!                ## s is near the boundary of tab2
!                if sboundflag and (not muboundflag):
!                    s = (is2 + is2_bound+1) * 0.5 * deltas2
!                    if False:
!                        rat = (s-scenter2) / (scenter2_bound-scenter2)
!                    else:
!                        scenter3=(scenter2+scenter2_bound) * 0.5
!                        ds = scenter3-scenter2
!                        d1 = s-scenter2+ds
!                        rat = d1 / (2*ds)
!                        rat = min(rat, 1)
!                        
!                    rat_bound = 1-rat
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
!                    if is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat_bound
!                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat_bound
!              
          if ( muboundflag .and. .not. sboundflag) then 
            angle = (iangle2 + iangle2_bound + 1) * 0.5 * deltamu2
            anglecenter3 = (anglecenter2+anglecenter2_bound) / 2.0
            dangle = anglecenter3 -anglecenter2
            d1 = angle-anglecenter2+dangle
            rat = min(d1 / (2.0_rt*dangle),1.0_rt)
            rat_bound = 1.0-rat

            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat
            
            if (iangle2_bound < nummu2) then
             smutab2(is2+1,iangle2_bound+1,1) = smutab2(is2+1,iangle2_bound+1,1)+countDD*rat_bound
             smutab2(is2+1,iangle2_bound+1,2) = smutab2(is2+1,iangle2_bound+1,2)+countDR*rat_bound
             smutab2(is2+1,iangle2_bound+1,3) = smutab2(is2+1,iangle2_bound+1,3)+countRR*rat_bound
            endif
          endif              

!                ## mu is near the boundary of tab2
!                if muboundflag and (not sboundflag):
!                    angle = (iangle2 + iangle2_bound+1) * 0.5 * deltamu2
!!                    if False:
 !                       rat = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
 !!                   else:
!!                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
!                        dangle = anglecenter3-anglecenter2
!                        d1 = angle-anglecenter2+dangle
!                        rat = d1 / (2*dangle)
!!                        rat = min(rat, 1)
                        
!                    rat_bound = 1-rat
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
!                    if iangle2_bound < nummu2:
!                      for row3 in compute_rows:
!                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat_bound
 !                     if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat_bound
 !               

          if ( muboundflag .and. sboundflag) then 
            s = (is2 + is2_bound + 1) * 0.5 * deltas2
            angle = (iangle2 + iangle2_bound + 1) * 0.5 * deltamu2
            scenter3 = (scenter2+scenter2_bound) / 2.0
            ds = scenter3 -scenter2
            d1 = s-scenter2+ds
            rats = min(d1 / (2.0_rt*ds),1.0_rt)

            anglecenter3 = (anglecenter2+anglecenter2_bound) / 2.0
            dangle = anglecenter3 -anglecenter2
            d1 = angle-anglecenter2+dangle
            ratangle = min(d1 / (2.0_rt*dangle),1.0_rt)
  !              ## both s, mu are near the boundy...
!                if muboundflag and sboundflag:
!                    s = (is2 + is2_bound+1) * 0.5 * deltas2
!                    angle = (iangle2 + iangle2_bound+1) * 0.5 * deltamu2
!                    if False:
!                        rats = (s-scenter2) / (scenter2_bound-scenter2)
!                        ratangle = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
!                    else:
!                        scenter3=(scenter2+scenter2_bound) * 0.5
!                        ds = scenter3-scenter2
!                        d1 = s-scenter2+ds
!                        rats = d1 / (2*ds)
!                        rats = min(rats, 1)
!!                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
 !                       dangle = anglecenter3-anglecenter2
 !                       d1 = angle-anglecenter2+dangle
 !                       ratangle = d1 / (2*dangle)
 !                       ratangle = min(ratangle, 1)            
            rat1 = rats*ratangle
            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat1
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat1
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat1

 !                   # original pixel
 !                   rat1 = rats*ratangle
 !                   for row3 in compute_rows:
 !                       smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat1
 !                   if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat1

            rat2 = (1-rats) * ratangle
            if (is2_bound < nums2) then
             smutab2(is2_bound+1,iangle2+1,1) = smutab2(is2_bound+1,iangle2+1,1)+countDD*rat2
             smutab2(is2_bound+1,iangle2+1,2) = smutab2(is2_bound+1,iangle2+1,2)+countDR*rat2
             smutab2(is2_bound+1,iangle2+1,3) = smutab2(is2_bound+1,iangle2+1,3)+countRR*rat2
            endif
            
            rat3=rats*(1-ratangle)
            if (iangle2_bound < nummu2) then
             smutab2(is2+1,iangle2_bound+1,1) = smutab2(is2+1,iangle2_bound+1,1)+countDD*rat3
             smutab2(is2+1,iangle2_bound+1,2) = smutab2(is2+1,iangle2_bound+1,2)+countDR*rat3
             smutab2(is2+1,iangle2_bound+1,3) = smutab2(is2+1,iangle2_bound+1,3)+countRR*rat3
            endif
            
            rat4=(1-rats)*(1-ratangle)
            if (iangle2_bound < nummu2 .and. is2_bound<nums2) then
             smutab2(is2_bound+1,iangle2_bound+1,1) = smutab2(is2_bound+1,iangle2_bound+1,1)+countDD*rat4
             smutab2(is2_bound+1,iangle2_bound+1,2) = smutab2(is2_bound+1,iangle2_bound+1,2)+countDR*rat4
             smutab2(is2_bound+1,iangle2_bound+1,3) = smutab2(is2_bound+1,iangle2_bound+1,3)+countRR*rat4
            endif
          endif              
          cycle ! haha	
!                    # diff s
 !                   rat2 = (1-rats)*ratangle
!                    if is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat2
!                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat2
!                    # diff angle
!                    rat3 = rats*(1-ratangle)
!                    if iangle2_bound < nummu2:
!                      for row3 in compute_rows:
!                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat3
!                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat3
!                    # diff s and diff angle
!                    rat4 = (1-rats)*(1-ratangle)
!                    if iangle2_bound < nummu2 and is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat4
!                      if save_counts_row != None: smutab2[is2_bound][iangle2_bound][save_counts_row] += rat4

      enddo
    enddo

!                            
!    if div_counts_rows != [] and save_counts_row != None:
!        for is2 in range(nums2):
!            for iangle2 in range(nummu2):
!                for row3 in div_counts_rows:
!                    smutab2[is2][iangle2][row3] /= smutab2[is2][iangle2][save_counts_row]
!    return smutab2
  end subroutine DSMapping
  
  subroutine smuintxi(smutab, deltas, nbins, mubins, anglemin, anglemax, smin, smax, nummuedge, intxi, mumids)
    ! argument
    integer, intent(in)  :: nbins, mubins, nummuedge
    real(rt), intent(in) :: smutab(nbins,mubins,3), deltas,  &
      anglemin,anglemax, smin,smax ! range of s,mu to be integrated
    real(rt), intent(out) :: intxi(nummuedge-1)
    real(rt), intent(out), optional :: mumids(nummuedge-1) ! middle value of mu in each mu-bin 
    ! variable
    real(rt) :: deltamu, angleindices(nummuedge), dangle, dangleindex, &
      DD,DR,RR,xi, dk,k1,k2
    integer :: i,j,k,sindex1, sindex2
    logical :: testprint=.false.

    intxi = 0.0_rt
    deltamu = 1.0 / dble(mubins)
    
    ! first of all, decide the fractional bins in angular space
    dangle = (anglemax-anglemin) / dble(nummuedge-1)
    dangleindex = dangle/deltamu
    angleindices(1) = anglemin / deltamu
    do j = 2, nummuedge
      angleindices(j) = angleindices(j-1) + dangleindex
    enddo
    if(testprint) print *, 'real(angleindices) = ', real(angleindices)
    do i = 1, nummuedge
      if(testprint) print *, i, angleindices(i) * deltamu
    enddo
    
    ! This only holds for deltas eq 1!
    sindex1 = int(smin / deltas+1.5)
    sindex2 = int(smax / deltas+0.5)
    if(abs(deltas-1.0).gt.1.0e-5) then
      print *, ' (smuintxi) Warning! deltas not equal to 1. May have significant error in the range of integral!'
      print *, ' deltas = ' , deltas
    endif
    if(testprint) print *, 'sindex1, sindex2 = ', sindex1, sindex2
        
    if(sindex2 > nbins) then
      print *, ' (smuintxi) WARNING!! index of s outflow: sindex, nbins = ', sindex2, nbins, '; forcing sindex2 equal to nbins'
      stop
      sindex2 = nbins
    endif

    if(present(mumids)) then
      do j = 1,nummuedge-1
	k1 = floor(angleindices(j)+1+0.00001) 
	k2 = floor(angleindices(j+1)+0.00001)
	mumids(j) = dble((k1-1) + k2)/2.0_rt*deltamu
      enddo
    endif

    do i = sindex1, sindex2    
      if(testprint) print *, 'Doing integration for s = ', i
    if (i.eq.sindex1.and.testprint) then
      do j =1,nummuedge-1
        k1 = angleindices(j)
        k2 = angleindices(j+1)
        print *,  floor(k1+1+0.00001), floor(k2+0.00001)
      enddo
    endif
    do j = 1, nummuedge-1
      DD=0.0_rt; DR=0.0_rt; RR=0.0_rt;
      k1 = angleindices(j)
      k2 = angleindices(j+1)
      if(i.eq.sindex1.and.testprint) print *, j
      ! summation of counts: the whole bins
      if (.false.) then
       if(i.eq.sindex1 ) print *, 'k1,k2 = ', k1,k2
       if(i.eq.sindex1 ) print *, 'ceiling(k1),floor(k2) = ', ceiling(k1),floor(k2)
       do k = ceiling(k1),floor(k2)
         if(i.eq.sindex1.and.testprint) print *, 'Integrating: ', k
         DD = DD + smutab(i,k,1)
         DR = DR + smutab(i,k,2)
         RR = RR + smutab(i,k,3)
       enddo
       ! summatin of counts: the fractional bins (tails)
       dk = (ceiling(k1)-k1)
       if(i.eq.sindex1.and.testprint ) print *, 'ceiling(k1)-k1 = ', dk
       DD = DD + smutab(i,ceiling(k1)-1,1)*dk
       DR = DR + smutab(i,ceiling(k1)-1,2)*dk
       RR = RR + smutab(i,ceiling(k1)-1,3)*dk
       dk = (k2-floor(k2))
       if(i.eq.sindex1.and. testprint) print *, 'k2-floor(k2) = ', dk
       DD = DD + smutab(i,floor(k2)+1,1)*dk
       DR = DR + smutab(i,floor(k2)+1,2)*dk
       RR = RR + smutab(i,floor(k2)+1,3)*dk
      else
      ! suppose k1 = 3.5; k2=4.5; then floor(k1)+1 = 4; floor(k-1)+1 = 4
      ! suppose k = 3; then floor(k)+1=4; floor(k-1)+1=3
       do k =  floor(k1+1+0.00001),floor(k2+0.00001)
         DD = DD + smutab(i,k,1)
         DR = DR + smutab(i,k,2)
         RR = RR + smutab(i,k,3)
       enddo
      endif
	      
      xi = (DD-2.0_rt*DR) / RR + 1.0_rt
      intxi(j) = intxi(j) + xi 
    enddo
    enddo
  end subroutine smuintxi

! normalize an array: amplitude shifted to 1; skip the LAST element
  subroutine normfun(A,nA,Anormed)
    integer, intent(in) :: nA
    real(rt), intent(in) :: A(nA)
    real(rt), intent(out) :: Anormed(nA-1)
    real(rt) :: avg
    integer :: imid, i
    avg = sum(A) / dble(nA)
    imid = nA / 2 + 1
!    do i = 1, imid-1
    do i = 1, nA-1
      Anormed(i) = A(i) / avg
    enddo
!    do i = imid, nA
!      Anormed(i-1) = A(i) / avg
!    enddo
  end subroutine normfun
  
  subroutine smu_ximu_calcchisqs(&
    omlist, numom, wlist, numw, & ! List of omegam, w
    outputdir, baseoutputfile, & ! Basic name of the outputfile
    omstd, wstd & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base
    )
    !arguments
    integer, intent(in) :: numom, numw
    real(rt), intent(in) :: omlist(numom), wlist(numw),  &
       omstd, wstd
    character(*), intent(in) :: outputdir, baseoutputfile
    
    ! variables
    real(rt) :: omwlist(2,numw,numom), chisqs_nosyscor(N1,N2), chisqs_syscorn(N1,N2), &
      smutabstds(nbins_database,mubins_database,3,nz), &
      smutab_data(nbins_data,mubins_data,3)
    real(rt) :: smin_mapping=1.0_rt, smax_mapping=50_rt, DAstd,DAnew,Hstd,Hnew, deltas1,deltas2,t0,t1,t2,dt,&
      intxis(maxval(mubins),N1,N2,nz), intxi(maxval(mubins)), dintxi(maxval(mubins)), &
      chisq_nosyscor(n1,n2,nz-1), chisq_syscor(n1,n2,nz-1), &
      chisq_nosyscor_allschemes(nz-1), chisq_syscor_allschemes(nz-1)
    type(par) :: parstd, parnew
    character(len=1000) :: tmpstr,tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,tmpstr9,&
      nowchisqstr,filenames(N1,N2), nowfile, filename_allschemes
    integer :: nowfileunit, fileunits(N1,N2), fileunit_allschemes, i,j,k,iz,i1,i2, n,iom,iw,iomw
    integer, parameter :: basefileunit = 58483
    
    ! 1. Initialization
    parstd%omegam = omstd; parstd%w = wstd
    ! 1.1 Check covmat files ready or not
    do i=1,N1
    do j=1,N2
    do iz=2,nz
      if (.not.allocated(covmats(i,j,iz-1)%A).or.covmats(i,j,iz-1)%nA.ne.mubins(i)-1) then
        print *, ' (smu_ximu_calcchisqs) ERROR! Covmats not ready: iz,i,j,mubin,mucut = ', iz,i,j,mubins(i),mucuts(j)
        stop
      endif
    enddo
    enddo
    enddo
    ! 1.2 Load in 2pCF computed in standard cosmology
    do iz = 1, nz
      call ximu_loadsmufile(data2pcffile_base(iz,omstd,wstd), smutabstds(:,:,:,iz), smax_database, nbins_database, mubins_database)
    enddo

    ! 2. Names of files storing chisq values
    write(tmpstr1,*) ncovmocks
    nowchisqstr = trim(adjustl(baseoutputfile)) // '.CovMock_'//trim(adjustl(tmpstr1)) 
    nowfileunit = basefileunit
    do i=1,N1
    do j=1,N2
      write(tmpstr1, *) mubins(i) ! how many bins in mu space
      write(tmpstr2, '(f5.2)') mucuts(j) ! maximal value of mu
      write(tmpstr3, '(f5.1)') ints1
      write(tmpstr4, '(f5.1)') ints2
      write(tmpstr5, *) ncovmocks
      tmpstr6=''; if(polyfitdeg.ge.0) write(tmpstr6,'(A,i1)') '.polyfitdeg',polyfitdeg
      tmpstr7=''; tmpstr8=''; tmpstr9=''; write(tmpstr7,*) nbins_database; write(tmpstr8,*) mubins_database; 
      tmpstr7='.baseline_nbins'//trim(adjustl(tmpstr7))//'_mubins'//trim(adjustl(tmpstr8))
      nowfile = trim(adjustl(outputdir))//'/'//trim(adjustl(nowchisqstr))//'__'//trim(adjustl(tmpstr1)) &
        //'mubins.mumax' // trim(adjustl(tmpstr2))  &
        //'.CovMock_'//trim(adjustl(tmpstr5)) &
        //'.s'//trim(adjustl(tmpstr3))//'to'//trim(adjustl(tmpstr4))//trim(adjustl(tmpstr6)) &
        //trim(adjustl(tmpstr7)) &
        //'.txt'
        
      filenames(i,j) = nowfile; fileunits(i,j) = nowfileunit
      print *, ' (smu_ximu_calcchisqs) Opening file for output: ', &
        trim(adjustl(nowfile))
      open(unit=fileunits(i,j),file=filenames(i,j))
      write(fileunits(i,j),'(A)') '### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor'
      nowfileunit = nowfileunit+1
    enddo
    enddo
    tmpstr3=''; write(tmpstr3,*) N1*N2
    fileunit_allschemes=fileunits(N1,N2)+1; 
    filename_allschemes = trim(adjustl(filenames(1,1)))//'.'//trim(adjustl(tmpstr3))//'schemes'
    open(unit=fileunit_allschemes,file=filename_allschemes)
    print *, ' (smu_ximu_calcchisqs) Opening file for output: ', trim(adjustl(filename_allschemes))
    write(fileunit_allschemes,'(A)') '### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor'    
    
    ! 3. Compute chisq values
    call cpu_time(t0); t1=t0
    dt = 10.0; iomw=0;    
    do iom = 1, numom
    do iw  = 1, numw
      chisq_nosyscor_allschemes=0.0_rt; chisq_syscor_allschemes = 0.0_rt;
      ! 3.1 Compute intxi
      iomw=iomw+1; call cpu_time(t2)
      if (t2-t1.gt.dt) then
        write(*,'(f10.1,A,i5,A,f4.1,A)') t2-t0, ' seconds passed.   #-(om,w) = ', &
           iomw, ' (',100*float(iomw)/float(numom*numw),'%)'
        t1=t2
      endif
      do iz  = 1, nz
        parnew%omegam = omlist(iom); parnew%w=wlist(iw)
        DAstd = DAofz(parstd,zeffs(iz)); Hstd = Hofz(parstd,zeffs(iz))
        DAnew = DAofz(parnew,zeffs(iz)); Hnew = Hofz(parnew,zeffs(iz))
        ! 3.1 smutabdata obtained from DSMapping
        deltas1 = smax_database / float(nbins_database) 
        deltas2 = smax_data / float(nbins_data)
!       print *, ' (smu_ximu_calcchisqs) Do DSMapping... iz = ', iz
        call DSMapping(smutabstds(:,:,:,iz), nbins_database, mubins_database, smutab_data, &
	  nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
!       print *, ' (smu_ximu_calcchisqs) DSMapping done: iz = ', iz
      
        do i1 = 1, N1
        do i2 = 1, N2
!         print *, ' (smu_ximu_calcchisqs) Compute \int xi...'
	  call smuintxi(smutab_data, deltas2, nbins_data, mubins_data, &
	    anglemin=1.0_rt-mucuts(i2), anglemax=1.0_rt, &
	    smin=ints1, smax=ints2, &
	    nummuedge=mubins(i1)+1, &
	    intxi=intxi(1:mubins(i1)))
!         print *, ' (smu_ximu_calcchisqs) Normalize \int xi...'
	  call normfun(intxi(1:mubins(i1)),mubins(i1),intxis(1:mubins(i)-1,i1,i2,iz)) 
        enddo
        enddo
!       print *, ' (smu_ximu_calcchisqs) Compute intxi done: iz = ', iz
      enddo
      ! 3.2 Compute chisqs
      do i1 = 1,N1
      do i2 = 1,N2
        do iz = 2, nz
          n = mubins(i1)
          dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz) - intxis(1:n-1,i1,i2,1)
          chisq_nosyscor(i1,i2,iz-1) = chisq_cov_xbar(dintxi, covmats(i1,i2,iz-1)%A, n-1)
          dintxi(1:n-1) = dintxi(1:n-1) - dintxi_syscor(1:n-1,i1,i2,iz-1)
          chisq_syscor(i1,i2,iz-1) = chisq_cov_xbar(dintxi, covmats(i1,i2,iz-1)%A, n-1)
        enddo
        chisq_nosyscor_allschemes(:) = chisq_nosyscor_allschemes(:) + chisq_nosyscor(i1,i2,:)/ dble(N1*N2)
        chisq_syscor_allschemes(:)   = chisq_syscor_allschemes(:) + chisq_syscor(i1,i2,:) / dble(N1*N2)
        ! write chisq values to files...
        tmpstr=''; tmpstr1=''
        write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
          sum(chisq_nosyscor(i1,i2,1:nz-1)), sum(chisq_syscor(i1,i2,1:nz-1))
        do k = 1,nz-1
          write(tmpstr1,'(e15.7)') chisq_nosyscor(i1,i2,k)
          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
        enddo
        do k = 1,nz-1
          write(tmpstr1,'(e15.7)') chisq_syscor(i1,i2,k)
          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
        enddo
        write(fileunits(i1,i2),'(A)') trim(adjustl(tmpstr))
      enddo
      enddo
      ! output allscheme chisq value
      tmpstr=''; tmpstr1=''
      write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
        sum(chisq_nosyscor_allschemes(1:nz-1)), sum(chisq_syscor_allschemes(1:nz-1))
!     print *, tmpstr
      do k = 1,nz-1
        write(tmpstr1,'(e15.7)') chisq_nosyscor_allschemes(k)
        tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
      enddo
      do k = 1,nz-1
        write(tmpstr1,'(e15.7)') chisq_syscor_allschemes(k)
        tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
      enddo
      write(fileunit_allschemes,'(A)') trim(adjustl(tmpstr))
    enddo
    enddo
    
    do i1 = 1,n1
    do i2 = 1,n2
      close(fileunits(i1,i2))
    enddo
    enddo
    close(fileunit_allschemes)
!    chisqs_nosyscor(N1,N2)
!    chisqs_syscorn(A1,N2)
  end subroutine smu_ximu_calcchisqs
    
  real(rt) function chisq_cov_xbar(xbar, invcov, n)
    integer, intent(in) :: n
    real(rt), intent(in) :: xbar(n), invcov(n,n)
    integer :: i,j
    chisq_cov_xbar = 0.0_rt
    do i = 1, n
    do j = 1, n
      chisq_cov_xbar = chisq_cov_xbar+xbar(i)*invcov(i,j)*xbar(j)
    enddo
    enddo
  end function chisq_cov_xbar
    
end module LSS_ximu_tools


module LSS_ximu_tests
use LSS_ximu_tools
implicit none
contains

subroutine check_load_files()
    character(len=charlen) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, nowfile
    integer :: i, j, iz, imock
    logical :: logvar
    type(par) :: parstd, parnew, nowpar
    real(rt) :: smutab_database(nbins_database, mubins_database, 3), &
                smutab_sysmock(nbins_sysmock, mubins_sysmock, 3), &
                smutab_covmock(nbins_covmock, mubins_covmock, 3)

    real(rt) :: smutab_data(nbins_data, mubins_data, 3)
    real(rt) :: smin_mapping=1.0_rt, smax_mapping=50_rt, DAstd, DAnew, Hstd, Hnew, deltas1, deltas2, intxi(1000), avg, &
      tmpx1,tmpx2,tmpx3,tmpx4,tmpx5,tmpx6,xi



    printstr = "Now it is empty!"
    intxi = 0.0_rt
    parstd%omegam = 0.26_rt; parstd%w = -1.0_rt
    parnew%omegam = 0.26_rt; parnew%w = -1.5_rt

    do iz = 1, nz
        print *, 'Redshift = ', zeffs(iz), '...'
        print *, 
    enddo

    print *, '  (LSS_ximu) checking existence of data base files...'
    do iz = 1, nz
        nowfile=data2pcffile_base(iz, 0.26_rt, -1.0_rt)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif

        call ximu_loadsmufile(nowfile, smutab_database, smax_database, nbins_database, mubins_database)
        print *, trim(adjustl(nowfile))
        !call DSMapping(smutabstd, nums1, nummu1, smutab2, &
	!  nums2, nummu2, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
	DAstd=DAofz(parstd,zeffs(iz))
	DAnew=DAofz(parnew,zeffs(iz))
	Hstd=Hofz(parstd,zeffs(iz))
	Hnew=Hofz(parnew,zeffs(iz))
	deltas1 = smax_database / float(nbins_database) 
	deltas2 = smax_data / float(nbins_data)
	call DSMapping(smutab_database, nbins_database, mubins_database, smutab_data, &
	  nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
	do i = 30, 50
	  j = 100
	  print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
	enddo
!	stop

	open(unit=37802,file='smutab_data.txt')
	do i = 1, nbins_data
	do j = 1, mubins_data
		write(37802,'(i6,i6,e14.7,e14.7,e14.7)') i-1, j-1, smutab_data(i,j,1:3)
	enddo
	enddo
	close(37802)
	
!	smuintxi(smutab, deltas, nbins, mubins, anglemin, anglemax, smin, smax, nummuedge, intxi)
	call smuintxi(smutab_data, deltas2, nbins_data, mubins_data, 0.05_rt, 1.0_rt, 6.0_rt, 40.0_rt, 26, intxi(1:25))
	print *, 'intxi = ', real(intxi(1:25))
	avg = sum(intxi(1:25)) / 25.0
	intxi(1:25) = intxi(1:25) / avg
	print *, 'intxi = ', real(intxi(1:25))
	stop
    enddo

    nowpar%omegam=0.26; nowpar%w=0.0
    do iz = 1, nz
        nowfile=data2pcffile_general(iz,nowpar)
        inquire(file=nowfile, exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for systematic correction)...'
    do iz = 1, nz
     do imock = 1, nsysmocks
        nowfile=syscor2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
        call ximu_loadsmufile(nowfile, smutab_sysmock, smax_sysmock, nbins_sysmock, mubins_sysmock)
     enddo
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for covmat estimation)...'
    do iz = 1, nz
     do imock = 1, ncovmocks
        nowfile=cov2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
        call ximu_loadsmufile(nowfile, smutab_covmock, smax_covmock, nbins_covmock, mubins_covmock)
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
end subroutine check_load_files
end module LSS_ximu_tests

program main
use LSS_ximu_tests
  integer, parameter :: numom=25, numw=51
  integer :: i,j,k,i1,i2,iz
  real(rt) :: omlist(numom), wlist(numw), ommin, ommax, wmin, wmax
  character(charlen) :: outputdir

!  call check_load_files()
!  stop

  ommin= 0.06_rt; ommax= 0.66_rt;
  wmin =-2.5_rt; wmax  = 0.0_rt
  do i = 1, numom
    omlist(i) = ommin + (ommax-ommin)/dble(numom-1.0_rt)*dble(i-1)
  enddo
  do j = 1, numw
    wlist(j) = wmin + (wmax-wmin)/dble(numw-1.0_rt)*dble(j-1)
  enddo

  print *,  'List of omegam to be computed:'
  do i =1, numom
   write(*,'(f8.4)',advance='no'), real(omlist(i)) 
  enddo
  print *
  print *,  'List of w to be computed:'
  do i =1, numw
   write(*,'(f8.4)',advance='no'), real(wlist(i)) 
  enddo
  print *

!  stop
!  print *, 'Compute covmats...'
!  call calc_covmats()
!  print *, 'Output covmats...'
!  call output_covmats('.170411')

  print *, 'Load in covmats...'
  call load_covmats('.170411')
  print *, 'Invert covmats...'
  call invert_covmats()
  print *, 'Compute systematic correction...'
  call calc_syscor()
!  do iz = 2, nz
!    print *, real(dintxi_syscor(1:mubins(1)-1,1,5,iz-1))
!  enddo
!  stop
!  call check_load_files()
  if(.true.) &
   call smu_ximu_calcchisqs(&
    omlist=omlist, numom=numom, &
    wlist=wlist, numw=numw, & ! List of omegam, w
    outputdir='/home/xiaodongli/LSS/2PCF_AP/chisqs', &
    baseoutputfile='TestScan_170410.nopolyfit', & ! Basic name of the outputfile
    omstd=0.26_rt, wstd=-1.0_rt & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base)
    )
end program
