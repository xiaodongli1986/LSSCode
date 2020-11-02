module cic_nuisance

use LSS_constants_types
implicit none
contains


!void cic_mass_assignment(Snapshot const * const snapshot, float* const grid, const int nc)
!{
!  const int np= snapshot->np_local;
!  ParticleMinimum* const p= snapshot->p;
!  const float dx_inv= nc/snapshot->boxsize;
!
!  for(int i=0; i<np; i++) {
!    int ix[3], ix0[3], ix1[3];
!    float w[3];

!    for(int j=0; j<3; ++j) {
!      ix[j]= (int) floor(p[i].x[j]*dx_inv);
!      w[j]= 1.0f - (p[i].x[j]*dx_inv - ix[j]);  // CIC weight for left point
!      ix0[j]= (ix[j] + nc) % nc;               // left grid (periodic)
!      ix1[j]= (ix[j] + 1 + nc) % nc;           // right grid (periodic)
!    }
!    
!    grid[(ix0[0]*nc + ix0[1])*nc + ix0[2]] += w[0]*w[1]*w[2];
!    grid[(ix0[0]*nc + ix1[1])*nc + ix0[2]] += w[0]*(1-w[1])*w[2];
!    grid[(ix0[0]*nc + ix0[1])*nc + ix1[2]] += w[0]*w[1]*(1-w[2]);
!    grid[(ix0[0]*nc + ix1[1])*nc + ix1[2]] += w[0]*(1-w[1])*(1-w[2]);
!
!    grid[(ix1[0]*nc + ix0[1])*nc + ix0[2]] += (1-w[0])*w[1]*w[2];
!    grid[(ix1[0]*nc + ix1[1])*nc + ix0[2]] += (1-w[0])*(1-w[1])*w[2];
!    grid[(ix1[0]*nc + ix0[1])*nc + ix1[2]] += (1-w[0])*w[1]*(1-w[2]);
!    grid[(ix1[0]*nc + ix1[1])*nc + ix1[2]] += (1-w[0])*(1-w[1])*(1-w[2]);
!  }
!}

end module cic_nuisance

program main_cic

use cic_nuisance

implicit none

	
	print *, " This program is empty!!"

end program

