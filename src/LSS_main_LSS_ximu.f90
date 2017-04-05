
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

  integer, parameter :: rmax_database=150, nbins_database=750, mubins_database=600, &
                        rmax_sysmock=150, nbins_sysmock=150, mubins_sysmock=120, &
                        rmax_covmock=51,  nbins_covmock=51, mubins_covmock=120
 
!----------------------------
! 1. 

contains


  real(rt) function DAofz(nowpar)
    type(par) :: nowpar
    DAofz=0.0
  end function DAofz
  real(rt) function Hofz(nowpar)
    type(par) :: nowpar
    Hofz=0.0
  end function Hofz



  character(25) function omwstr(omegam, w)
    real(rt) :: omegam, w
    character(10) :: str1, str2
    write(str1, '(f8.4)') omegam
    write(str2, '(f8.4)') w
    omwstr = 'om'//trim(adjustl(str1))//'_w'//trim(adjustl(str2))
    !omstr = str3
  end function omwstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1. Name of 2-point correlation function (2pcf) files 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------
! data (baseline)    
  character(len=charlen) function data2pcffile_base(i_redshiftbin, ombase, wbase)
    integer :: i_redshiftbin
    real(rt) :: ombase, wbase
    character(len=charlen) :: filename, filestr
    filestr = '.rmax150.750rbins.600mubins.'
    if (i_redshiftbin .eq. 1) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (i_redshiftbin .eq. 2) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (i_redshiftbin .eq. 3) then
      filename = trim(adjustl(filedir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (i_redshiftbin .eq. 4) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (i_redshiftbin .eq. 5) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (i_redshiftbin .eq. 6) then
      filename = trim(adjustl(filedir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'//&
	trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function for loading in the xi(s,mu) table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ximu_loadsmufile(filename, smutab, smax, nbins, mubins)
    ! argument
    character(len=charlen), intent(in) :: filename
    integer, intent(in) :: smax, nbins, mubins
    real(rt), intent(out) :: smutab(nbins,mubins,4) ! Four values at each s/mu: DD, DR, RR, xi
    ! variables
    character(len=charlen) :: nowstr
    integer :: i, j
    real :: tmpx

    smutab = 0.0_rt
    
    open(unit=44817,file=filename)
    print *, filename
    read(44817,*) nowstr
    !print *, nowstr
    do i = 1, nbins
    do j = 1, mubins
    !    print *, i,j
        read(44817,*) tmpx, tmpx, tmpx, tmpx, smutab(i,j,1:3) !, tmpx, tmpx, smutab(i,j,4)
    enddo
    enddo
    close(44817)
  end subroutine

!def smu__CosmoConvert(s,mu,DA1,DA2,H1,H2):
!    ''' s1: angular direction; s2: LOS direction '''
!    s2 = s*mu;
!    s1 = np.sqrt(s*s - s2*s2)
!    alpha1 = DA2 / DA1
!    alpha2 = H1 / H2
!    s_prime  =  np.sqrt((alpha1*s1)**2 + (alpha2*s2)**2)
!    mu_prime =  alpha2*s2 / s_prime
!    return s_prime, mu_prime

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function for mapping xi(s,mu) from dense to sparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  subroutine mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd, nums1, nummu1, smutab2, &
	  nums2, nummu2, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping )
  !      compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[], method='simple_bin'):
    ! argument
    real(rt), intent(in) :: smutabstd(nums1,nummu1,3), DAstd, DAnew, Hstd, Hnew, &
      deltas1, deltas2, smin_mapping, smax_mapping
    integer, intent(in) :: nums1, nummu1, nums2, nummu2
    real(rt), intent(out):: smutab2(nums2,nummu2,3)
    ! variables
    real(rt) :: deltamu1,deltamu2, mubound1,mubound2,sbound1,sbound2, &
      scenter,anglecenter,mucenter,scenter2,anglecenter2,mucenter2
    real(rt) :: smutabstd_centers(nums2,nummu2,2)
    integer :: maxs1, mins1, is1,iangle1,is2,iangle2, smutabstd_ismus(nums2,nummu2,2)
    
    deltamu1 = 1.0d0/real(nummu1); deltamu2 = 1.0d0/real(nummu2)
    smutab2=0.0d0
    ! range of s for smutab1, smutab2
    call smu__CosmoConvert(smax_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smax_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    maxs1 = floor( max(sbound1,sbound2) / deltas1 + 0.5 )
    call smu__CosmoConvert(smin_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smin_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    mins1 = floor( min(sbound1,sbound2) / deltas1 + 0.5 )
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

    do is1 = mins1, maxs1-1
      do iangle1 = 0, nummu1-1
        scenter=(is1+0.5_rt)*deltas1; anglecenter=(iangle1+0.5)*deltamu1
        mucenter=1.0_rt-anglecenter
        call smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,scenter2,mucenter2)
        anglecenter2 = 1.0_rt - mucenter2
        is2 = floor(scenter2 / deltas2); iangle2 = floor(anglecenter2/deltamu2)
        if ((is2 < nums2) .and. (iangle2 < nummu2)) then
          smutabstd_centers(is1+1,iangle1+1,1) = scenter2
          smutabstd_centers(is1+1,iangle1+1,2) = anglecenter2
          smutabstd_ismus(is1+1,iangle1+1,1) = is2
          smutabstd_ismus(is1+1,iangle1+1,2) = iangle2
        endif
      enddo
    enddo

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
    do is1 = mins1, maxs1
      do iangle1 = 0, nummu1-1
        scenter2=smutabstd_centers(is1+1,iangle1+1,1)
        anglecenter2=smutabstd_centers(is1+1,iangle1+1,2)
	is2=smutabstd_ismu(is1+1,iangle1+1,1)
	iangle2=smutabstd_ismu(is1+1,iangle1+1,2)
        if (is2.ge.nums2 .or. iangle2.ge.nummu2) cycle
        
 !    if method == 'divided_pixel':
!        for is1 in range(mins1, nums1):
!            for iangle1 in range(nummu1):
!                scenter2, anglecenter2, is2, iangle2 = smutabstd_centers[is1][iangle1]
!                if not (is2 < nums2 and iangle2 < nummu2):
!                    continue       
        
        sboundflag=.false.; muboundflag=.false.
        do is1_nearby = max(is1-1,mins1), min(is1+1,nums1-1)
          is2_b = smutabstd_ismu(is1_nearby+1,iangle1+1,1)
          iangle2_b = smutabstd_ismu(is1_nearby+1,iangle1+1,2)
          if (is2_b .ne. is2) then 
            sboundflag=.true.; is1_bound=is1_nearby; is2_bound=is2_b;
            scenter2_bound=smutabstd_centers(is1_nearby,iangle1,1)
            !seriously possible bug found in the python code! Decide to stop and checking the python code for a while.
          
            
!                ### firstly, check boundary:
!                sboundflag, muboundflag = False, False
!                for is1_nearby in [max(is1-1,mins1), min(is1+1,nums1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1_nearby][iangle1][2], smutabstd_centers[is1_nearby][iangle1][3]
!                    if is2_b != is2:
!                        sboundflag=True; is1_bound = is1_nearby; is2_bound = is2_b; 
!                        scenter2_bound = smutabstd_centers[is1_nearby][iangle1][1]
!                for iangle1_nearby in [max(iangle1-1,0), min(iangle1+1,nummu1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1][iangle1_nearby][2], smutabstd_centers[is1][iangle1_nearby][3]
!                    if iangle2_b != iangle2:
!                        muboundflag=True; iangle1_bound = iangle1_nearby; iangle2_bound = iangle2_b
!                        anglecenter2_bound = smutabstd_centers[is1][iangle1_nearby][2]
         
      enddo
    enddo

!                    

!                
!                ### Then, treat them case by case...
!                ## s, mu are all not near the boundary of tab 2
!                if ((not sboundflag)and(not muboundflag)):
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += 1
!                            
!                ## s is near the boundary of tab2
!                if sboundflag and (not muboundflag):
!                    s = (is2 + is2_bound) * 0.5 * deltas2
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
!                ## mu is near the boundary of tab2
!                if muboundflag and (not sboundflag):
!                    angle = (iangle2 + iangle2_bound) * 0.5 * deltamu2
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
  !              ## both s, mu are near the boundy...
!                if muboundflag and sboundflag:
!                    s = (is2 + is2_bound) * 0.5 * deltas2
!                    angle = (iangle2 + iangle2_bound) * 0.5 * deltamu2
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
 !                   # original pixel
 !                   rat1 = rats*ratangle
 !                   for row3 in compute_rows:
 !                       smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat1
 !                   if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat1
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
!                            
!    if div_counts_rows != [] and save_counts_row != None:
!        for is2 in range(nums2):
!            for iangle2 in range(nummu2):
!                for row3 in div_counts_rows:
!                    smutab2[is2][iangle2][row3] /= smutab2[is2][iangle2][save_counts_row]
!    return smutab2
  end subroutine mapping_smudata_to_another_cosmology_DenseToSparse
end module LSS_ximu_tools




program main_LSS_ximu

use LSS_ximu_tools

implicit none

    character(len=charlen) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, nowfile
    integer :: i, j, i_redshiftbin, imock
    logical :: logvar
    type(par) :: nowpar
    real(rt) :: smutab_database(nbins_database, mubins_database, 4), &
                smutab_sysmock(nbins_sysmock, mubins_sysmock, 4), &
                smutab_covmock(nbins_covmock, mubins_covmock, 4)

    

    printstr = "Now it is empty!"

    print *, '  (LSS_ximu) checking existence of data base files...'
    do i_redshiftbin = 1, 6
        nowfile=data2pcffile_base(i_redshiftbin, 0.26_rt, -1.0_rt)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, i_redshiftbin, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
        call ximu_loadsmufile(nowfile, smutab_database, rmax_database, nbins_database, mubins_database)
    enddo
    do i_redshiftbin = 1, 100
     do i = 1, nbins_database
     do j = 1, mubins_database
     smutab_database(i,j,:) = smutab_database(i,j,:)+0.1
     enddo
     enddo
     print *, i_redshiftbin
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
        call ximu_loadsmufile(nowfile, smutab_sysmock, rmax_sysmock, nbins_sysmock, mubins_sysmock)
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
        call ximu_loadsmufile(nowfile, smutab_covmock, rmax_covmock, nbins_covmock, mubins_covmock)
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
