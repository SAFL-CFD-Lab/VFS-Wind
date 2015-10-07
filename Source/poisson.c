/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscmg.h"
#include <stdlib.h>

extern int pseudo_periodic, rans, implicit;
extern PetscInt block_number, freesurface, immersed, ti, tistart, les, tiout;
extern char path[256];

#define lidx2(i,j,k,user)	(int)(gid[k][j][i])

//PetscErrorCode VolumeFlux(UserCtx *user, Vec lUcor, PetscReal *ibm_Flux, PetscReal *ibm_Area);

double time_coeff()
{
  if(levelset || rans) {
	if(ti==tistart) return 1.;
	else return 1.5;
  }
  else return 1.;
	
};
		
PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnom, void *dummy)
{
  UserCtx *user = (UserCtx*)dummy;
  Vec x;
  PetscInt tttt, mi, j, k;
  PetscReal norm;


  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  VecDuplicate(user->P, &x);
  KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &x);
  VecMax(x, &tttt, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "KSP Max %d %le\n", tttt, norm);
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (lidx(mi,j,k,user) ==tttt) {
	  PetscPrintf(PETSC_COMM_SELF, "KspMMax %d %d %d %d %le\n", lidx(mi,j,k,user),mi,j, k, norm);
	}
      }
    }
  }

  VecMin(x, &tttt, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "KSP Max %d %le\n", tttt, norm);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (lidx(mi,j,k,user) ==tttt) {
	  PetscPrintf(PETSC_COMM_SELF, "KspMMin %d %d %d %d %le\n", lidx(mi,j,k,user),mi,j, k, norm);
	}
      }
    }
  }

  VecDestroy(x);

  return 0;
}

PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)
{
	DALocalInfo	info = user->info;
	PetscInt	gxs, gxe, gys, gye, gzs, gze;

	gxs = info.gxs; gxe = gxs + info.gxm;
	gys = info.gys; gye = gys + info.gym;
	gzs = info.gzs; gze = gzs + info.gzm;


	if (!(user->aotopetsc)) {
		user->aotopetsc = PETSC_TRUE;
		DAGetGlobalIndices(user->da, PETSC_NULL, &user->idx_from);
	}
	
	return (user->idx_from[(k-gzs) * (info.gxm*info.gym) + (j-gys)*(info.gxm) + (i-gxs)]);
  
  //  return (i + j * mx + k * mx * my);
}

void Convert_Phi2_Phi(UserCtx *user)
{
	DALocalInfo	info = user->info;
	DA		da = user->da;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, ***phi, ***lphi, *phi2;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	VecSet(user->Phi,0);
	DAVecGetArray(da, user->Phi, &phi);
	DAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(user->Phi2, &phi2);
	
	int pos=0;
	for(k=lzs; k<lze; k++)
	for(j=lys; j<lye; j++)
	for(i=lxs; i<lxe; i++) {
		/*if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
			phi[k][j][i]=0;
		}
		else */
		if( (int) nvert[k][j][i]>poisson_threshold ) {
			/*if(movefsi || rotatefsi) {
				phi[k][j][i] = phi2[pos++];
			}
			else*/ phi[k][j][i]=0;
		}
		else {
			phi[k][j][i] = phi2[pos++];
		}
	}
	
	DAVecRestoreArray(da, user->Phi, &phi);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(user->Phi2, &phi2);
	
	DAGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
	DAGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
	
	DAVecGetArray(da, user->lPhi, &lphi);
	DAVecGetArray(da, user->Phi, &phi);
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		if(flag) {
			phi[k][j][i] = lphi[c][b][a];
		}
		
		//if(k==1 && i!=0)	printf("%d,%d,%d, %f %f %f\n", i,j,k, lphi[-2][j][i], lphi[k][j][i], lphi[k+1][j][i]);
	}
	DAVecRestoreArray(da, user->lPhi, &lphi);
	DAVecRestoreArray(da, user->Phi, &phi);
	//exit(0);
}

PetscInt setup_lidx2(UserCtx *user)
{
	DALocalInfo	info = user->info;
	DA		da = user->da;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, ***gid, ***lid;

	Vec	Lid;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	if(ti==tistart) {
		VecDuplicate(user->lNvert, &user->Gid);
	}
	VecDuplicate(user->lNvert, &Lid);
	
	VecSet(user->Gid, -1);
	VecSet(Lid, -1);
	
	DAVecGetArray(da, user->Gid, &gid);
	DAVecGetArray(da, Lid, &lid);
	DAVecGetArray(da, user->lNvert, &nvert);

	int r, myrank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
		
	std::vector<int> ndof_node(size), ndof_node_tmp(size);	// # of pressure dof for processors
	
	int ndof_node_accu;
	
	for(r=0; r<size; r++) {
		ndof_node_tmp[r] = 0;
	}
	
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		
		if(poisson==-1) {
			lid[k][j][i] = (PetscReal)ndof_node_tmp[myrank]++;
		}
		else {
			if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) { }
			else if( (int)nvert[k][j][i]>poisson_threshold ) {
				/*if(movefsi || rotatefsi) {
					lid[k][j][i] = (PetscReal)ndof_node_tmp[myrank]++;
				}
				else {}*/
			}
			else {
				lid[k][j][i] = (PetscReal)ndof_node_tmp[myrank]++;
			}
		}
	}
	
	MPI_Allreduce( &ndof_node_tmp[0], &ndof_node[0], size, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		
	ndof_node_accu = 0;
	for(r=0; r<myrank; r++) ndof_node_accu += ndof_node[r];
		
	int n;
	user->p_global_begin = ndof_node_accu;
		
	VecGetSize(user->Phi,&n);
	if(myrank==size-1) {
		printf("\n\n********* %d %d ***********\n\n", ndof_node_accu + ndof_node[myrank], n);
		user->reduced_p_size = ndof_node_accu + ndof_node[myrank];
	}
		
	MPI_Bcast(&user->reduced_p_size, 1, MPI_INT, size-1, PETSC_COMM_WORLD);
		
	PetscBarrier(PETSC_NULL);
		
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if((int)(lid[k][j][i])>=0) {
			gid[k][j][i] = lid[k][j][i] + ndof_node_accu;	// gid is double, be careful
		}
	}
		
	user->local_Phi2_size = ndof_node[myrank];
	
	if(ti!=tistart) VecDestroy (user->Phi2);
	
	VecCreateMPI(PETSC_COMM_WORLD, ndof_node[myrank], PETSC_DETERMINE, &user->Phi2);
		
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->Gid, &gid);
	DAVecRestoreArray(da, Lid, &lid);
	
	VecDestroy(Lid);

        DALocalToLocalBegin(da, user->Gid, INSERT_VALUES, user->Gid);
        DALocalToLocalEnd(da, user->Gid, INSERT_VALUES, user->Gid);
	
	if(periodic) {
		DAVecGetArray(da, user->Gid, &gid);
		for(k=zs; k<ze; k++)
		for(j=ys; j<ye; j++)
		for(i=xs; i<xe; i++) {
			int flag=0, a=i, b=j, c=k;
			
			if(i_periodic && i==0) a=mx-2, flag=1;
			else if(i_periodic && i==mx-1) a=1, flag=1;
			
			if(j_periodic && j==0) b=my-2, flag=1;
			else if(j_periodic && j==my-1) b=1, flag=1;
			
			if(k_periodic && k==0) c=mz-2, flag=1;
			else if(k_periodic && k==mz-1) c=1, flag=1;
			
			if(ii_periodic && i==0) a=-2, flag=1;
			else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
			
			if(jj_periodic && j==0) b=-2, flag=1;
			else if(jj_periodic && j==my-1) b=my+1, flag=1;
			
			if(kk_periodic && k==0) c=-2, flag=1;
			else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
			
			if(flag) gid[k][j][i] = gid[c][b][a];
			
			//if(i==1 && k==mz-2) printf("%d, %d, %d, %d <= %d %d %d \n", (int)gid[k][j][i], k,j,i, c,b,a);
		}
		DAVecRestoreArray(da, user->Gid, &gid);
		
		DALocalToLocalBegin(da, user->Gid, INSERT_VALUES, user->Gid);
		DALocalToLocalEnd(da, user->Gid, INSERT_VALUES, user->Gid);
	}	
	
	return 0;
}

PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user)
{
  if (*(user->isc)) {
    *ih = i;
  }
  else {
    *ih = 2 * i;
  }

  if (*(user->jsc)) {
    *jh = j;
  }
  else {
    *jh = 2 * j;
  }

  if (*(user->ksc)) {
    *kh = k;
  }
  else {
    *kh = 2 * k;
  }

  return 0;
}

PetscErrorCode MyRestriction(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);


  
  DA	da = user->da;

  DA	da_f = *user->da_f;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  
  PetscReal ***f, ***x, ***nvert;
  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  DAVecGetArray(da,   F, &f);

  Vec lX;
  //  DAGetLocalVector(da_f, &lX);
  DACreateLocalVector(da_f, &lX);
  DAGlobalToLocalBegin(da_f, X, INSERT_VALUES, lX);
  DAGlobalToLocalEnd(da_f, X, INSERT_VALUES, lX);  
  DAVecGetArray(da_f, lX, &x);

  DAVecGetArray(da, user->lNvert, &nvert);

  PetscReal ***nvert_f;
  DAVecGetArray(da_f, user->user_f->lNvert, &nvert_f);

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (k==0) {
	  f[k][j][i] = 0.;
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;
	}
	else if (j==0) {
	  f[k][j][i] = 0.;
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;
	}
	else if (i==0) {
	  f[k][j][i] = 0.;
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;
	}
	else {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user);
	  f[k][j][i] = 0.125 *
	    (x[kh   ][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) +
	     x[kh   ][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) +
	     x[kh   ][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) +
	     x[kh   ][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) +
	     x[kh-ka][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) +
	     x[kh-ka][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia]));

/* 	  PetscReal scale; */
/* 	  scale = (PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) + */
/* 		   PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia])); */
/* 	  if (scale > 0) f[k][j][i] *= (8./scale); */

	  if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
	}
      }
    }
  }


  DAVecRestoreArray(da_f, user->user_f->lNvert, &nvert_f);

  DAVecRestoreArray(da_f, lX, &x);
  VecDestroy(lX);
  //  DARestoreLocalVector(da_f, &lX);
  DAVecRestoreArray(da,   F,  &f);
  DAVecRestoreArray(da, user->lNvert, &nvert);
/*   PetscReal norm; */
/*   VecNorm(F, NORM_2, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Restriction Norm %le\n", norm); */
  return 0;
}

#define GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user) \
  if (*(user->isc)) { \
    ic = i; \
    ia = 0; \
  } \
  else { \
    ic = (i+1) / 2; \
    ia = (i - 2 * (ic)) == 0 ? 1 : -1; \
    if (i==1 || i==mx-2) ia = 0; \
  }\
  if (*(user->jsc)) { \
    jc = j; \
    ja = 0; \
  } \
  else { \
    jc = (j+1) / 2; \
    ja = (j - 2 * (jc)) == 0 ? 1 : -1; \
    if (j==1 || j==my-2) ja = 0; \
  } \
  if (*(user->ksc)) { \
    kc = k; \
    ka = 0; \
  } \
  else { \
    kc = (k+1) / 2; \
    ka = (k - 2 * (kc)) == 0 ? 1 : -1; \
    if (k==1 || k==mz-2) ka = 0; \
  } \
  if (ka==-1 && nvert_c[kc-1][jc][ic] > 0.1) ka = 0; \
  else if (ka==1 && nvert_c[kc+1][jc][ic] > 0.1) ka = 0; \
  if (ja==-1 && nvert_c[kc][jc-1][ic] > 0.1) ja = 0; \
  else if (ja==1 && nvert_c[kc][jc+1][ic] > 0.1) ja = 0; \
  if (ia==-1 && nvert_c[kc][jc][ic-1] > 0.1) ia = 0; \
  else if (ia==1 && nvert_c[kc][jc][ic+1] > 0.1) ia = 0;
 

/* PetscErrorCode MyCoarseBC(Mat A, Vec X) */
/* { */
/*   UserCtx *user; */

/*   MatShellGetContext(A, (void*)&user); */

/*   UserCtx *user_c; */
/*   user_c = user->user_c; */

/*   DALocalInfo	info = user_c->info; */
/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*   PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
/*   PetscInt	lxs, lxe, lys, lye, lzs, lze; */

/* } */
PetscErrorCode MyInterpolation(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);

  
  
  DA	da = user->da;

  DA	da_c = *user->da_c;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert, ***nvert_c;
  PetscInt i, j, k, ic, jc, kc, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


  DAVecGetArray(da,   F, &f);


  Vec lX;
  DACreateLocalVector(da_c, &lX);
  //  DAGetLocalVector(da_c, &lX);
  DAGlobalToLocalBegin(da_c, X, INSERT_VALUES, lX);
  DAGlobalToLocalEnd(da_c, X, INSERT_VALUES, lX);  
  DAVecGetArray(da_c, lX, &x);

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);
/* 	  kc = (k + 1) / 2; */
/* 	  jc = (j + 1) / 2; */
/* 	  ic = (i + 1) / 2; */
/* 	  ka = (k - 2 * kc)==0 ? 1 : -1; */
/* 	  ja = (j - 2 * jc)==0 ? 1 : -1; */
/* 	  ia = (i - 2 * ic)==0 ? 1 : -1; */
/* 	  if (ka==-1 &&(k==1 || nvert_c[kc-1][jc][ic]>0.1)) ka = 0; */
/* 	  else if (ka==1 && (k==mz-2 || nvert_c[kc+1][jc][ic]>0.1)) ka=0; */

/* 	  if (ja==-1 &&(j==1 || nvert_c[kc][jc-1][ic]>0.1)) ja = 0; */
/* 	  else if (ja==1 &&(j==my-2 || nvert_c[kc][jc+1][ic]>0.1)) ja=0; */

/* 	  if (ia==-1 &&(i==1 || nvert_c[kc][jc][ic-1]>0.1)) ia = 0; */
/* 	  else if (ia==1 && (i==mx-2 || nvert_c[kc][jc][ic+1]>0.1)) ia=0; */

/* 	f[k][j][i] = x[kc][jc][ic]; */
	  f[k][j][i] = (x[kc   ][jc   ][ic   ] * 9 +
			x[kc   ][jc+ja][ic   ] * 3 +
			x[kc   ][jc   ][ic+ia] * 3 +
			x[kc   ][jc+ja][ic+ia]) * 3./64. +
	    (x[kc+ka][jc   ][ic   ] * 9 +
	     x[kc+ka][jc+ja][ic   ] * 3 +
	     x[kc+ka][jc   ][ic+ia] * 3 +
	     x[kc+ka][jc+ja][ic+ia]) /64.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0) {
	  f[k][j][i] = 0.;//-f[k][j][i+1];
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;//-f[k][j][i-1];
	}
	else if (j==0) {
	  f[k][j][i] = 0.;//-f[k][j+1][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;//-f[k][j-1][i];
	}
	else if (k==0) {
	  f[k][j][i] = 0.;//-f[k+1][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;//-f[k-1][j][i];
	}
	if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
/* 	  f[k][j][i] = 0.125 * */
/* 	    (x[kc  ][jh  ][ih  ] + */
/* 	     x[kc  ][jh  ][ih-1] + */
/* 	     x[kh  ][jh-1][ih  ] + */
/* 	     x[kh-1][jh  ][ih  ] + */
/* 	     x[kh  ][jh-1][ih-1] + */
/* 	     x[kh-1][jh-1][ih  ] + */
/* 	     x[kh-1][jh  ][ih-1] + */
/* 	     x[kh-1][jh-1][ih-1]); */
/* 	} */
      }
    }
  }

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);

  DAVecRestoreArray(da_c, lX, &x);
  //  DARestoreLocalVector(da_c, &lX);
  VecDestroy(lX);
  DAVecRestoreArray(da,   F,  &f);

/*   PetscReal sum; */
/*   PetscInt N; */
/*   VecSum(F, &sum); */
/*   VecGetSize(F, &N); */
/*   sum = sum / (-1.*N); */
/*   VecShift(F, sum); */

/*   PetscReal max, min; */
/*   PetscInt  maxi, mini; */
/*   VecMax(X, &maxi, &max); */
/*   VecMin(X, &mini, &min); */
/*   PetscPrintf(PETSC_COMM_WORLD, "MM %d %le %d %le\n", maxi, max, mini, min); */
/*   PetscReal norm; */
/*   VecNorm(F, NORM_2, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Interpolation Norm %le\n", norm); */

  return 0;
  //  MatNullSpaceRemove(da.nullsp, F, PETSC_NULL);
}

PetscErrorCode PoissonNullSpaceFunction(Vec X, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  DA da = user->da;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	/* ***x,*/ ***nvert;
  PetscInt	i, j, k;
	
	PetscReal	***aj, ***gid;
	DAVecGetArray(user->da, user->lAj, &aj);
	DAVecGetArray(user->da, user->Gid, &gid);
	

/*   /\* First remove a constant from the Vec field X *\/ */
/*   MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp); */
/*   MatNullSpaceRemove(nullsp, X, PETSC_NULL); */
/*   MatNullSpaceDestroy(nullsp); */

  /* Then apply boundary conditions */
  
  DAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal lsum, sum;
  PetscReal  lnum, num;

  if (user->multinullspace) PetscPrintf(PETSC_COMM_WORLD, "MultiNullSpace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  if (!user->multinullspace) {
    lsum = 0;
    lnum = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    //lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    //PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    VecSum(X, &sum);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum = sum / (-1.0 * num);
/*     PetscPrintf(PETSC_COMM_WORLD, "NullSpace %e\n", sum); */
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if((int)nvert[k][j][i]==0 && i!=0 && i!=mx-1 && j!=0 && j!=my-1 && k!=0 && k!=mz-1) {
		//x[k][j][i] +=sum;
		double val=sum;
		VecSetValue(X, (int)gid[k][j][i], val, ADD_VALUES);
	  }
	}
      }
    }
  }
  VecAssemblyBegin(X);
  VecAssemblyEnd(X);
  
  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(user->da, user->lAj, &aj);
  DAVecRestoreArray(user->da, user->Gid, &gid);
  return 0;
}

PetscErrorCode PoissonNullSpaceFunction_original(MatNullSpace nsp, Vec X, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  DA da = user->da;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	***x, ***nvert;
  PetscInt	i, j, k;

/*   /\* First remove a constant from the Vec field X *\/ */
/*   MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp); */
/*   MatNullSpaceRemove(nullsp, X, PETSC_NULL); */
/*   MatNullSpaceDestroy(nullsp); */

  /* Then apply boundary conditions */
  DAVecGetArray(da, X, &x);
  DAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal lsum, sum;
  PetscReal  lnum, num;

  if (user->multinullspace) PetscPrintf(PETSC_COMM_WORLD, "MultiNullSpace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  if (!user->multinullspace) {
    lsum = 0;
    lnum = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum = sum / (-1.0 * num);
/*     PetscPrintf(PETSC_COMM_WORLD, "NullSpace %e\n", sum); */
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    x[k][j][i] +=sum;
	  }
	}
      }
    }
  }
  else {
    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }

    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }
			   
  }
  if (zs == 0) {
    k = 0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k = mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] > 0.1)
	  x[k][j][i] = 0.;
      }
    }
  }
  DAVecRestoreArray(da, X, &x);
  DAVecRestoreArray(da, user->lNvert, &nvert);

  return 0;
}

PetscErrorCode mymatmultadd(Mat mat, Vec v1, Vec v2, Vec v3)
{

  Vec vt;
  VecDuplicate(v3, &vt);
  MatMult(mat, v1, vt);
  VecWAXPY(v3, 1., v2, vt);
  VecDestroy(vt);
  return(0);
}


#define CP 0
#define EP 1
#define WP 2
#define NP 3
#define SP 4
#define TP 5
#define BP 6
#define NE 7
#define SE 8
#define NW 9
#define SW 10
#define TN 11
#define BN 12
#define TS 13
#define BS 14
#define TE 15
#define BE 16
#define TW 17
#define BW 18

PetscErrorCode PoissonLHSNew(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;

	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;

	PetscReal	***aj, ***iaj, ***jaj, ***kaj;

	PetscInt lxs, lxe, lys, lye, lzs, lze;
	PetscInt gxs, gxe, gys, gye, gzs, gze;

	Vec		G11, G12, G13, G21, G22, G23, G31, G32, G33;
	PetscReal	***g11, ***g12, ***g13, ***g21, ***g22, ***g23;
	PetscReal	***g31, ***g32, ***g33;
	PetscReal	***rho;

	PetscReal	***nvert, ***nvert_o, ***gid, ***level;
	PetscScalar	vol[27];
	PetscInt	idx[27], row;

	PetscInt	i, j, k, N;
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	gxs = info.gxs; gxe = gxs + info.gxm;
	gys = info.gys; gye = gys + info.gym;
	gzs = info.gzs; gze = gzs + info.gzm;
	
	if (!user->assignedA) {
		user->assignedA = PETSC_TRUE;
		
		setup_lidx2(user);
		N = mx * my * mz;
		PetscInt M;
		
		MatCreate(PETSC_COMM_WORLD, &(user->A));
		VecGetLocalSize(user->Phi2, &M);
		MatSetSizes(user->A,M,M,PETSC_DETERMINE,PETSC_DETERMINE);
		MatSetType(user->A,MATMPIAIJ);
		MatMPIAIJSetPreallocation(user->A, 19, PETSC_NULL, 19, PETSC_NULL);
		MatSetFromOptions(user->A);
	}

	MatZeroEntries(user->A);

	if (levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lLevelset, &level);
	}
	
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);

	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lIEta, &ieta);
	DAVecGetArray(fda, user->lIZet, &izet);

	DAVecGetArray(fda, user->lJCsi, &jcsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lJZet, &jzet);

	DAVecGetArray(fda, user->lKCsi, &kcsi);
	DAVecGetArray(fda, user->lKEta, &keta);
	DAVecGetArray(fda, user->lKZet, &kzet);

	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);

	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lNvert_o, &nvert_o);

	VecDuplicate(user->lAj, &G11);
	VecDuplicate(user->lAj, &G12);
	VecDuplicate(user->lAj, &G13);
	VecDuplicate(user->lAj, &G21);
	VecDuplicate(user->lAj, &G22);
	VecDuplicate(user->lAj, &G23);
	VecDuplicate(user->lAj, &G31);
	VecDuplicate(user->lAj, &G32);
	VecDuplicate(user->lAj, &G33);
/*	
	VecSet(G11,1.e10);
	VecSet(G12,1.e10);
	VecSet(G13,1.e10);
	VecSet(G21,1.e10);
	VecSet(G22,1.e10);
	VecSet(G23,1.e10);
	VecSet(G31,1.e10);
	VecSet(G32,1.e10);
	VecSet(G33,1.e10);
*/
	DAVecGetArray(da, G11, &g11);
	DAVecGetArray(da, G12, &g12);
	DAVecGetArray(da, G13, &g13);
	DAVecGetArray(da, G21, &g21);
	DAVecGetArray(da, G22, &g22);
	DAVecGetArray(da, G23, &g23);
	DAVecGetArray(da, G31, &g31);
	DAVecGetArray(da, G32, &g32);
	DAVecGetArray(da, G33, &g33);
	
	/*for (k=gzs; k<gze; k++)
	for (j=gys; j<gye; j++)
	for (i=gxs; i<gxe; i++) */
	for (k=lzs-1; k<lze+1; k++)
	for (j=lys-1; j<lye+1; j++)
	for (i=lxs-1; i<lxe+1; i++) 
	{
		int a=i, b=j, c=k;
		
		
			int i_flag=0, j_flag=0, k_flag=0;
			
			if(i_periodic && i==0) a=mx-2, i_flag=1;
			else if(i_periodic && i==mx-1) a=1, i_flag=1;
			
			if(j_periodic && j==0) b=my-2, j_flag=1;
			else if(j_periodic && j==my-1) b=1, j_flag=1;
			
			if(k_periodic && k==0) c=mz-2, k_flag=1;
			else if(k_periodic && k==mz-1) c=1, k_flag=1;
			
			if(ii_periodic && i==0) a=-2, i_flag=1;
			else if(ii_periodic && i==mx-1) a=mx+1, i_flag=1;
			
			if(jj_periodic && j==0) b=-2, j_flag=1;
			else if(jj_periodic && j==my-1) b=my+1, j_flag=1;
			
			if(kk_periodic && k==0) c=-2, k_flag=1;
			else if(kk_periodic && k==mz-1) c=mz+1, k_flag=1;
			
		g11[k][j][i] = (icsi[c][b][a].x * icsi[c][b][a].x + icsi[c][b][a].y * icsi[c][b][a].y + icsi[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g12[k][j][i] = (ieta[c][b][a].x * icsi[c][b][a].x + ieta[c][b][a].y * icsi[c][b][a].y + ieta[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g13[k][j][i] = (izet[c][b][a].x * icsi[c][b][a].x + izet[c][b][a].y * icsi[c][b][a].y + izet[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g21[k][j][i] = (jcsi[c][b][a].x * jeta[c][b][a].x + jcsi[c][b][a].y * jeta[c][b][a].y + jcsi[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g22[k][j][i] = (jeta[c][b][a].x * jeta[c][b][a].x + jeta[c][b][a].y * jeta[c][b][a].y + jeta[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g23[k][j][i] = (jzet[c][b][a].x * jeta[c][b][a].x + jzet[c][b][a].y * jeta[c][b][a].y + jzet[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g31[k][j][i] = (kcsi[c][b][a].x * kzet[c][b][a].x + kcsi[c][b][a].y * kzet[c][b][a].y + kcsi[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
		g32[k][j][i] = (keta[c][b][a].x * kzet[c][b][a].x + keta[c][b][a].y * kzet[c][b][a].y + keta[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
		g33[k][j][i] = (kzet[c][b][a].x * kzet[c][b][a].x + kzet[c][b][a].y * kzet[c][b][a].y + kzet[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];	
	}
	
	PetscInt m;
	DAVecGetArray(da, user->Gid, &gid);	// for macro gid2()

	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		double one=1.0;
		/*
		if(poisson==-1) row = lidx(i, j, k, user);
		else*/ 
		row = lidx2(i, j, k, user);
		
		if ( i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {
			if(poisson==-1) {
				if(periodic) row = lidx(i, j, k, user);	// tricky!!!!!!!!!!!!!!!!!!!!!!!! seokkoo
				MatSetValues(user->A, 1, &row, 1, &row, &one, INSERT_VALUES);
			}
		}
		else if ( nvert[k][j][i]>0.1 && row>=0 ) {	// for fsi
			MatSetValues(user->A, 1, &row, 1, &row, &one, INSERT_VALUES);
		}
		else {
			if (nvert[k][j][i] > poisson_threshold) { // i, j, k is not a fluid point
				continue;
			}
			else { // i, j, k is a fluid point
				for (m=0; m<19; m++) vol[m] = 0.;
			    
				/* Contribution from i+1 - i */
				if (nvert[k][j][i+1] < poisson_threshold && (i != mx-2 || i_periodic || ii_periodic) ) { // i+1, j, k is a fluid point
					double r = 1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k][j][i+1] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						r = mean ( rho[k][j][i+1], rho[k][j][i] );
					}
					
					/* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
					vol[CP] -= g11[k][j][i] / r; //i, j, k
					vol[EP] += g11[k][j][i] / r; // i+1, j, k
					
					/* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1}) * 0.25 * g12[k][j][i] */
					if ( (j == my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i] + nvert[k][j+1][i+1] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < poisson_threshold && j!=1 ) {
							vol[CP] += g12[k][j][i] * 0.5 / r; //i, j, k
							vol[EP] += g12[k][j][i] * 0.5 / r; // i+1, j, k
							vol[SP] -= g12[k][j][i] * 0.5 / r; //i, j-1, k
							vol[SE] -= g12[k][j][i] * 0.5 / r; //i+1, j-1, k
						}
					}
					else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < poisson_threshold) {
							vol[NP] += g12[k][j][i] * 0.5 / r;  //i, j+1, k
							vol[NE] += g12[k][j][i] * 0.5 / r; //i+1, j+1, k
							vol[CP] -= g12[k][j][i] * 0.5 / r; //i, j, k
							vol[EP] -= g12[k][j][i] * 0.5 / r; //i+1, j, k
						}
					}
					else {
						vol[NP] += g12[k][j][i] * 0.25 / r; // i, j+1, k
						vol[NE] += g12[k][j][i] * 0.25 / r; // i+1, j+1, k
						vol[SP] -= g12[k][j][i] * 0.25 / r; // i, j-1, k
						vol[SE] -= g12[k][j][i] * 0.25 / r; // i+1, j-1, k
						
					}
					/* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1}) * 0.25 / r * g13[k][j][i] */
					if ( (k == mz-2 && !k_periodic  && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < poisson_threshold && k!=1) {
							vol[CP] += g13[k][j][i] * 0.5 / r; // i, j, k
							vol[EP] += g13[k][j][i] * 0.5 / r; // i+1, j, k
							vol[BP] -= g13[k][j][i] * 0.5 / r; // i, j, k-1
							vol[BE] -= g13[k][j][i] * 0.5 / r; // i+1, j, k-1
						
						}
					}
					else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j][i+1] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < poisson_threshold) {
							vol[TP] += g13[k][j][i] * 0.5 / r;  // i, j, k+1
							vol[TE] += g13[k][j][i] * 0.5 / r; // i+1, j, k+1
							vol[CP] -= g13[k][j][i] * 0.5 / r;  // i, j, k
							vol[EP] -= g13[k][j][i] * 0.5 / r;  // i+1, j, k
						}
					}
					else {
						vol[TP] += g13[k][j][i] * 0.25 / r; //i, j, k+1
						vol[TE] += g13[k][j][i] * 0.25 / r; //i+1, j, k+1
						vol[BP] -= g13[k][j][i] * 0.25 / r; //i, j, k-1
						vol[BE] -= g13[k][j][i] * 0.25 / r; //i+1, j, k-1
						
					}
				}  // end of i+1 - i

				/* Contribution from i - i-1 */
				if (nvert[k][j][i-1] < poisson_threshold && (i != 1 || i_periodic || ii_periodic) ) { // i-1, j, k is a fluid point
					double r=1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k][j][i-1] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						r = mean ( rho[k][j][i-1], rho[k][j][i] );
					}
					
					/* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
					vol[CP] -= g11[k][j][i-1] / r;  //i, j, k
					vol[WP] += g11[k][j][i-1] / r;  //i-1, j, k
					
					/* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1}) * 0.25 / r * g12[k][j][i-1] */
					if ( (j == my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i] + nvert[k][j+1][i-1] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < poisson_threshold && j!=1 ) {
							vol[CP] -= g12[k][j][i-1] * 0.5 / r; // i, j, k
							vol[WP] -= g12[k][j][i-1] * 0.5 / r; // i-1, j, k
							vol[SP] += g12[k][j][i-1] * 0.5 / r; //i, j-1, k
							vol[SW] += g12[k][j][i-1] * 0.5 / r; // i-1, j-1, k
						}
					}
					else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k][j-1][i-1] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < poisson_threshold) {
							vol[NP] -= g12[k][j][i-1] * 0.5 / r; // i, j+1, k
							vol[NW] -= g12[k][j][i-1] * 0.5 / r; // i-1, j+1, k
							vol[CP] += g12[k][j][i-1] * 0.5 / r; // i, j, k
							vol[WP] += g12[k][j][i-1] * 0.5 / r; // i-1, j, k
						}
					}
					else {
						vol[NP] -= g12[k][j][i-1] * 0.25 / r; // i, j+1, k
						vol[NW] -= g12[k][j][i-1] * 0.25 / r; //i-1, j+1, k
						vol[SP] += g12[k][j][i-1] * 0.25 / r; // i, j-1, k
						vol[SW] += g12[k][j][i-1] * 0.25 / r; // i-1, j-1, k
					}

					/* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1}) * 0.25 / r * g13[k][j][i] */
					if ( (k == mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j][i-1] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < poisson_threshold && k!=1 ) {
							vol[CP] -= g13[k][j][i-1] * 0.5 / r; // i, j, k
							vol[WP] -= g13[k][j][i-1] * 0.5 / r; // i-1, j, k
							vol[BP] += g13[k][j][i-1] * 0.5 / r; // i, j, k-1
							vol[BW] += g13[k][j][i-1] * 0.5 / r; // i-1, j, k-1
						}
					}
					else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j][i-1] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < poisson_threshold) {
							vol[TP] -= g13[k][j][i-1] * 0.5 / r; // i, j, k+1
							vol[TW] -= g13[k][j][i-1] * 0.5 / r; //i-1, j, k+1
							vol[CP] += g13[k][j][i-1] * 0.5 / r; // i, j, k
							vol[WP] += g13[k][j][i-1] * 0.5 / r; //i-1, j, k
						}
					}
					else {
						vol[TP] -= g13[k][j][i-1] * 0.25 / r;  // i, j, k+1
						vol[TW] -= g13[k][j][i-1] * 0.25 / r; // i-1, j, k+1
						vol[BP] += g13[k][j][i-1] * 0.25 / r;  // i, j, k-1
						vol[BW] += g13[k][j][i-1] * 0.25 / r; // i-1, j, k-1
					}
				} // end of i - i-1

				/* Contribution from j+1 - j */
				if (nvert[k][j+1][i] < poisson_threshold && (j != my-2 || j_periodic || jj_periodic || (levelset && user->bctype[3]==4 && j==my-2) ) ) {
					double r=1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k][j+1][i] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						r = mean ( rho[k][j+1][i], rho[k][j][i] );
					}
										
					/* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) * 0.25 / r */
					if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] += g21[k][j][i] * 0.5 / r; // i, j, k
							vol[NP] += g21[k][j][i] * 0.5 / r; // i, j+1, k
							vol[WP] -= g21[k][j][i] * 0.5 / r; // i-1, j, k
							vol[NW] -= g21[k][j][i] * 0.5 / r; // i-1, j+1, k
						}
					}
					else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < poisson_threshold) {
							vol[EP] += g21[k][j][i] * 0.5 / r; // i+1, j, k
							vol[NE] += g21[k][j][i] * 0.5 / r; // i+1, j+1, k
							vol[CP] -= g21[k][j][i] * 0.5 / r; // i, j, k
							vol[NP] -= g21[k][j][i] * 0.5 / r; // i, j+1, k
						}
					}
					else {
						vol[EP] += g21[k][j][i] * 0.25 / r; //i+1, j, k
						vol[NE] += g21[k][j][i] * 0.25 / r; //i+1, j+1, k
						vol[WP] -= g21[k][j][i] * 0.25 / r; //i-1, j, k
						vol[NW] -= g21[k][j][i] * 0.25 / r; //i-1, j+1, k
					}
					
					// zero Dirichlet pressure
					if(levelset && user->bctype[3]==4 && j==my-2) {
						vol[CP] -= g22[k][j][i] / 0.5 / r;
					}
					else {
						/* dpde{j} = (p{j+1} - p{j}) * g22[k][j][i] */
						vol[CP] -= g22[k][j][i] / r;
						vol[NP] += g22[k][j][i] / r;
					}
					
					/* dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25 / r*/
					if ( (k == mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j+1][i] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < poisson_threshold && k!=1 ) {
							vol[CP] += g23[k][j][i] * 0.5 / r; //i,j,k
							vol[NP] += g23[k][j][i] * 0.5 / r; //i, j+1, k
							vol[BP] -= g23[k][j][i] * 0.5 / r;//i, j, k-1
							vol[BN] -= g23[k][j][i] * 0.5 / r;//i, j+1, k-1
						}
					}
					else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j+1][i] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < poisson_threshold) {
							vol[TP] += g23[k][j][i] * 0.5 / r; //i, j, k+1
							vol[TN] += g23[k][j][i] * 0.5 / r;//i, j+1, k+1
							vol[CP] -= g23[k][j][i] * 0.5 / r;//i, j, k
							vol[NP] -= g23[k][j][i] * 0.5 / r;//i, j+1, k
						}
					}
					else {
						vol[TP] += g23[k][j][i] * 0.25 / r; // i, j, k+1
						vol[TN] += g23[k][j][i] * 0.25 / r; // i, j+1, k+1
						vol[BP] -= g23[k][j][i] * 0.25 / r; // i, j, k-1
						vol[BN] -= g23[k][j][i] * 0.25 / r; // i, j+1, k-1
					}
				} // End of j+1 - j

				/* Contribution j - j-1 */
				if (nvert[k][j-1][i] < poisson_threshold && (j!=1 || j_periodic || jj_periodic) ) {
					double r=1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k][j-1][i] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						r = mean ( rho[k][j-1][i], rho[k][j][i] );
					}
					
					/* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) * 0.25 / r * g21[k][j-1][i] */
					if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k][j-1][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] -= g21[k][j-1][i] * 0.5 / r;// i, j, k
							vol[SP] -= g21[k][j-1][i] * 0.5 / r;// i, j-1, k
							vol[WP] += g21[k][j-1][i] * 0.5 / r;// i-1, j, k
							vol[SW] += g21[k][j-1][i] * 0.5 / r;// i-1, j-1, k
						}
					}
					else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k][j-1][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < poisson_threshold) {
							vol[EP] -= g21[k][j-1][i] * 0.5 / r;//i+1, j, k
							vol[SE] -= g21[k][j-1][i] * 0.5 / r;//i+1, j-1, k
							vol[CP] += g21[k][j-1][i] * 0.5 / r;//i, j, k
							vol[SP] += g21[k][j-1][i] * 0.5 / r;//i, j-1, k
						}
					}
					else {
						vol[EP] -= g21[k][j-1][i] * 0.25 / r;// i+1, j, k
						vol[SE] -= g21[k][j-1][i] * 0.25 / r;// i+1, j-1, k
						vol[WP] += g21[k][j-1][i] * 0.25 / r;// i-1, j, k
						vol[SW] += g21[k][j-1][i] * 0.25 / r;// i-1, j-1, k
					}
			      
					/* -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i] */
					vol[CP] -= g22[k][j-1][i] / r;
					vol[SP] += g22[k][j-1][i] / r;

					/* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) * 0.25 / r * g23[k][j-1][i] */
					if ( (k == mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j-1][i] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < poisson_threshold && k!=1 ) {
							vol[CP] -= g23[k][j-1][i] * 0.5 / r;//i, j, k
							vol[SP] -= g23[k][j-1][i] * 0.5 / r;//i, j-1, k
							vol[BP] += g23[k][j-1][i] * 0.5 / r;//i, j, k-1
							vol[BS] += g23[k][j-1][i] * 0.5 / r;//i, j-1, k-1
						}
					}
					else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j-1][i] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < poisson_threshold) {
							vol[TP] -= g23[k][j-1][i] * 0.5 / r;//i, j, k+1
							vol[TS] -= g23[k][j-1][i] * 0.5 / r;//i, j-1, k+1
							vol[CP] += g23[k][j-1][i] * 0.5 / r;//i, j, k
							vol[SP] += g23[k][j-1][i] * 0.5 / r;//i, j-1, k
						}
					}
					else {
						vol[TP] -= g23[k][j-1][i] * 0.25 / r;//i, j, k+1
						vol[TS] -= g23[k][j-1][i] * 0.25 / r;//i, j-1, k+1
						vol[BP] += g23[k][j-1][i] * 0.25 / r;//i, j, k-1
						vol[BS] += g23[k][j-1][i] * 0.25 / r;//i, j-1, k-1
					}
				} // End of j - j-1

				/* contribution from k+1 - k */
				if (nvert[k+1][j][i] < poisson_threshold && (k != mz-2 || k_periodic || kk_periodic || (levelset && user->bctype[5]==4 && k==mz-2) ) ) {
					double r=1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k+1][j][i] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						r = mean ( rho[k+1][j][i], rho[k][j][i] );
					}
					
					/* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) * 0.25 / r * g31[k][j][i] */
					if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k+1][j][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] += g31[k][j][i] * 0.5 / r;//i, j, k
							vol[TP] += g31[k][j][i] * 0.5 / r;//i, j, k+1
							vol[WP] -= g31[k][j][i] * 0.5 / r;//i-1, j, k
							vol[TW] -= g31[k][j][i] * 0.5 / r;//i-1, j, k+1
						}
					}
					else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k+1][j][i-1] > poisson_threshold) {
							if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < poisson_threshold) {
								vol[EP] += g31[k][j][i] * 0.5 / r;//i+1, j, k
								vol[TE] += g31[k][j][i] * 0.5 / r;//i+1, j, k+1
								vol[CP] -= g31[k][j][i] * 0.5 / r;//i, j, k
								vol[TP] -= g31[k][j][i] * 0.5 / r;//i, j, k+1
							}
					}
					else {
						vol[EP] += g31[k][j][i] * 0.25 / r;//i+1, j, k
						vol[TE] += g31[k][j][i] * 0.25 / r;//i+1, j, k+1
						vol[WP] -= g31[k][j][i] * 0.25 / r;//i-1, j, k
						vol[TW] -= g31[k][j][i] * 0.25 / r;//i-1, j, k+1
					}

					/* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 0.25 / r * g32[k][j][i] */
					if ( (j == my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i] + nvert[k+1][j+1][i] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < poisson_threshold && j!=1 ) {
							vol[CP] += g32[k][j][i] * 0.5 / r;//i, j,k
							vol[TP] += g32[k][j][i] * 0.5 / r;//i, j, k+1
							vol[SP] -= g32[k][j][i] * 0.5 / r;//i, j-1, k
							vol[TS] -= g32[k][j][i] * 0.5 / r;//i, j-1, k+1
						}
					}
					else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k+1][j-1][i] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < poisson_threshold) {
							vol[NP] += g32[k][j][i] * 0.5 / r;//i, j+1, k
							vol[TN] += g32[k][j][i] * 0.5 / r;//i, j+1, k+1
							vol[CP] -= g32[k][j][i] * 0.5 / r;//i, j, k
							vol[TP] -= g32[k][j][i] * 0.5 / r;//i, j, k+1
						}
					}
					else {
						vol[NP] += g32[k][j][i] * 0.25 / r;//i, j+1, k
						vol[TN] += g32[k][j][i] * 0.25 / r;//i, j+1, k+1
						vol[SP] -= g32[k][j][i] * 0.25 / r;//i, j-1, k
						vol[TS] -= g32[k][j][i] * 0.25 / r;//i, j-1, k+1
					}
					
					// zero Dirichlet pressure
					if( levelset && user->bctype[5] == 4 && k==mz-2) {
						vol[CP] -= g33[k][j][i] / 0.5 / r;
					}
					else
					{
						/* dpdz{k} = p{k+1} - p{k} */
						vol[CP] -= g33[k][j][i] / r; //i, j, k
						vol[TP] += g33[k][j][i] / r; //i, j, k+1
					}
				} // End of k+1 - k

				/* Contribution from k - k-1 */
				if (nvert[k-1][j][i] < poisson_threshold && (k != 1 || k_periodic || kk_periodic || (levelset && user->bctype[4]==5000 && k==1 && level[k][j][i]<0)) ) {
					double r=1.0;
					if (levelset) {
						double lm = 0.5 * ( level[k-1][j][i] + level[k][j][i] );
						//r = rho_air + (rho_water - rho_air) * H ( lm );
						if(k==1) r = rho[k][j][i];
						else r = mean ( rho[k-1][j][i], rho[k][j][i] );
					}
					
					/* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) * 0.25 / r * g31[k-1][j][i] */
					if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k-1][j][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] -= g31[k-1][j][i] * 0.5 / r;//i, j, k
							vol[BP] -= g31[k-1][j][i] * 0.5 / r;//i, j, k-1
							vol[WP] += g31[k-1][j][i] * 0.5 / r;//i-1, j, k
							vol[BW] += g31[k-1][j][i] * 0.5 / r;//i-1, j, k-1
						}
					}
					else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k-1][j][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < poisson_threshold) {
							vol[EP] -= g31[k-1][j][i] * 0.5 / r;//i+1, j, k
							vol[BE] -= g31[k-1][j][i] * 0.5 / r;//i+1, j, k-1
							vol[CP] += g31[k-1][j][i] * 0.5 / r;//i, j, k
							vol[BP] += g31[k-1][j][i] * 0.5 / r;//i, j, k-1
						}
					}
					else {
						vol[EP] -= g31[k-1][j][i] * 0.25 / r;//i+1, j, k
						vol[BE] -= g31[k-1][j][i] * 0.25 / r;//i+1, j, k-1
						vol[WP] += g31[k-1][j][i] * 0.25 / r;//i-1, j, k
						vol[BW] += g31[k-1][j][i] * 0.25 / r;//i-1, j, k-1
					}
			      
					/* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) *  0.25 / r * g32[k-1][j][i] */
					// ( p{i, j+1,k-1/2} - p{i, j-1,k-1/2} ) / (2eta)
					if ( (j == my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i] + nvert[k-1][j+1][i] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < poisson_threshold && j!=1 ) {
							vol[CP] -= g32[k-1][j][i] * 0.5 / r;//i, j,k
							vol[BP] -= g32[k-1][j][i] * 0.5 / r;//i, j, k-1
							vol[SP] += g32[k-1][j][i] * 0.5 / r;//i, j-1, k 
							vol[BS] += g32[k-1][j][i] * 0.5 / r;//i, j-1, k-1
						}
					}
					else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k-1][j-1][i] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < poisson_threshold) {
							vol[NP] -= g32[k-1][j][i] * 0.5 / r;//i, j+1, k
							vol[BN] -= g32[k-1][j][i] * 0.5 / r;//i, j+1, k-1
							vol[CP] += g32[k-1][j][i] * 0.5 / r;//i, j, k
							vol[BP] += g32[k-1][j][i] * 0.5 / r;//i, j, k-1
						}
					}
					else {
						vol[NP] -= g32[k-1][j][i] * 0.25 / r;//i, j+1, k
						vol[BN] -= g32[k-1][j][i] * 0.25 / r;//i, j+1, k-1
						vol[SP] += g32[k-1][j][i] * 0.25 / r;//i, j-1, k
						vol[BS] += g32[k-1][j][i] * 0.25 / r;//i, j-1, k-1
					}
					// zero Dirichlet pressure
                                        if( levelset && user->bctype[4] == 5 && k==1) {
						vol[CP] -= g33[k-1][j][i] / 0.5 / r;
                                        }
                                        else
					{
					  /* -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i] */
					  vol[CP] -= g33[k-1][j][i] / r; // i, j, k
					  vol[BP] += g33[k-1][j][i] / r; //i, j, k-1
					}
				} // End of k - k-1
				
				if(poisson==-1) for (m=0; m<19; m++)  vol[m] *= aj[k][j][i];
				/*
				if(poisson==-1) {
					idx[CP] = lidx(i  , j  , k  , user);
					idx[EP] = lidx(i+1, j  , k  , user);
					idx[WP] = lidx(i-1, j  , k  , user);
					idx[NP] = lidx(i  , j+1, k  , user);
					idx[SP] = lidx(i  , j-1, k  , user);
					idx[TP] = lidx(i  , j  , k+1, user);
					idx[BP] = lidx(i  , j  , k-1, user);
					idx[NE] = lidx(i+1, j+1, k  , user);
					idx[SE] = lidx(i+1, j-1, k  , user);
					idx[NW] = lidx(i-1, j+1, k  , user);
					idx[SW] = lidx(i-1, j-1, k  , user);
					idx[TN] = lidx(i  , j+1, k+1, user);
					idx[BN] = lidx(i  , j+1, k-1, user);
					idx[TS] = lidx(i  , j-1, k+1, user);
					idx[BS] = lidx(i  , j-1, k-1, user);
					idx[TE] = lidx(i+1, j  , k+1, user);
					idx[BE] = lidx(i+1, j  , k-1, user);
					idx[TW] = lidx(i-1, j  , k+1, user);
					idx[BW] = lidx(i-1, j  , k-1, user);
				}
				else*/ {
					idx[CP] = lidx2(i  , j  , k  , user);
					idx[EP] = lidx2(i+1, j  , k  , user);
					idx[WP] = lidx2(i-1, j  , k  , user);
					idx[NP] = lidx2(i  , j+1, k  , user);
					idx[SP] = lidx2(i  , j-1, k  , user);
					idx[TP] = lidx2(i  , j  , k+1, user);
					idx[BP] = lidx2(i  , j  , k-1, user);
					idx[NE] = lidx2(i+1, j+1, k  , user);
					idx[SE] = lidx2(i+1, j-1, k  , user);
					idx[NW] = lidx2(i-1, j+1, k  , user);
					idx[SW] = lidx2(i-1, j-1, k  , user);
					idx[TN] = lidx2(i  , j+1, k+1, user);
					idx[BN] = lidx2(i  , j+1, k-1, user);
					idx[TS] = lidx2(i  , j-1, k+1, user);
					idx[BS] = lidx2(i  , j-1, k-1, user);
					idx[TE] = lidx2(i+1, j  , k+1, user);
					idx[BE] = lidx2(i+1, j  , k-1, user);
					idx[TW] = lidx2(i-1, j  , k+1, user);
					idx[BW] = lidx2(i-1, j  , k-1, user);
				}
				
				//MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);
				//for(m=0; m<19; m++) if( fabs(vol[m])>1.e-70 || m==CP ) MatSetValues(user->A, 1, &row, 1, &idx[m], &vol[m], INSERT_VALUES);
				for(m=0; m<19; m++) {
					if( (fabs(vol[m])>1.e-10 && idx[m]>=0) || m==CP ) {
						MatSetValues(user->A, 1, &row, 1, &idx[m], &vol[m], INSERT_VALUES);
					}
				}

				/*
				if( (i==1 && j==2 && k==1) || (i==1 && j==2 && k==2)  || (i==1 && j==2 && k==3)  || (i==1 && j==2 && k==59) || (i==1 && j==2 && k==60) ) {
					printf("%d %d %d, my_idx = %d\n", i, j, k, idx[0]);
					for(m=0; m<19; m++) printf("%d\t", m);
					printf("\n");
					for(m=0; m<19; m++) printf("%d\t", idx[m]);
					printf("\n");
					for(m=0; m<19; m++) printf("%.2f\t", vol[m]);
					printf("\n\n");
				}*/
			} // End of fluid point
		} // End of interial points
	}
    //exit(0);
	DAVecRestoreArray(da, user->Gid, &gid);

	//kangsk  
	PetscPrintf(PETSC_COMM_WORLD, "MatAssemblyBegin...\n");
	MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);
	PetscPrintf(PETSC_COMM_WORLD, "MatAssemblyEnd...\n");

	DAVecRestoreArray(da, G11, &g11);
	DAVecRestoreArray(da, G12, &g12);
	DAVecRestoreArray(da, G13, &g13);
	DAVecRestoreArray(da, G21, &g21);
	DAVecRestoreArray(da, G22, &g22);
	DAVecRestoreArray(da, G23, &g23);
	DAVecRestoreArray(da, G31, &g31);
	DAVecRestoreArray(da, G32, &g32);
	DAVecRestoreArray(da, G33, &g33);
  
	VecDestroy(G11);
	VecDestroy(G12);
	VecDestroy(G13);
	VecDestroy(G21);
	VecDestroy(G22);
	VecDestroy(G23);
	VecDestroy(G31);
	VecDestroy(G32);
	VecDestroy(G33);

	if (levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lLevelset, &level);
	}
	
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);

	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);

	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);

	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);

	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);

	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lNvert_o, &nvert_o);
	
	return 0;
}


PetscErrorCode PoissonRHS(UserCtx *user, Vec B)
{
  DALocalInfo info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscReal	***nvert, ***aj, ***rb, ***gid, dt = user->dt;
  Cmpnts	***ucont;

 // DAVecGetArray(user->da, B, &rb);
  DAVecGetArray(user->fda, user->lUcont, &ucont);
  DAVecGetArray(user->da, user->lNvert, &nvert);
  DAVecGetArray(user->da, user->lAj, &aj);
  
	DAVecGetArray(user->da, user->Gid, &gid);
	

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	      double val;
	if (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
		val = 0.;
	}
	else if (nvert[k][j][i] >= poisson_threshold) {
		val = 0;
	}
	else {
		double coeff=time_coeff();
		#ifdef DIRICHLET
		//if(freesurface && j==user->free_surface_j[i][k]) coeff *= user->vof[i][k];
		#endif
		val =- (ucont[k][j][i].x - ucont[k][j][i-1].x + ucont[k][j][i].y - ucont[k][j-1][i].y + ucont[k][j][i].z - ucont[k-1][j][i].z ) / dt * user->st * coeff;// * aj[k][j][i]; 
	}
	VecSetValue(B, lidx(i,j,k,user), val, INSERT_VALUES);
      }
    }
  }
  
  VecAssemblyBegin(B);
  VecAssemblyEnd(B);

  DAVecRestoreArray(user->da, user->Gid, &gid);
  
 // DAVecRestoreArray(user->da, B, &rb);

#ifndef DIRICHLET  
	int N;
	double sum0;
	VecGetSize(B,&N);
	VecSum(B,&sum0);
	sum0  = sum0/(-1.0*N);
	VecShift(B,sum0);
        user->multinullspace = PETSC_FALSE;
	PoissonNullSpaceFunction(B, user);
#endif 
  
  DAVecGetArray(user->da, B, &rb);
  
   PetscReal lsum=0, sum=0; 
   for (k=zs; k<ze; k++) { 
     for (j=ys; j<ye; j++) { 
       for (i=xs; i<xe; i++) { 
 	if (i==0 || i==mx-1 || j==0 || j==my-1 || 
 	    k==0 || k==mz-1) { 
 	  lsum += 0.; 
 	} 
 	else { 
 	  lsum += rb[k][j][i];// / aj[k][j][i]; 
 	} 
       } 
     } 
   }

   PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); 
   PetscPrintf(PETSC_COMM_WORLD, "Summation RHS %e\n", sum); 
	
  DAVecRestoreArray(user->fda, user->lUcont, &ucont);
  DAVecRestoreArray(user->da, user->lNvert, &nvert);
  DAVecRestoreArray(user->da, user->lAj, &aj);
  DAVecRestoreArray(user->da, B, &rb);

  return 0;
}


PetscErrorCode PoissonRHS2(UserCtx *user, Vec B)
{
	DALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
	PetscReal	***nvert, ***aj, ***gid, dt = user->dt;
	Cmpnts	***ucont, ***lucont;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  	
	DAVecGetArray(user->fda, user->lUcont, &ucont);
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->da, user->lAj, &aj);
  	DAVecGetArray(user->da, user->Gid, &gid);
	
	int lcount=0;
		
	VecSet(B, 0);
		
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double val;
		if (nvert[k][j][i] >= poisson_threshold) {
			if( (int) (gid[k][j][i]) >=0 ) val = 0;	// for fsi
			else continue;
		}
		else {
			double coeff=time_coeff();
			#ifdef DIRICHLET
			//if(freesurface && j==user->free_surface_j[i][k]) coeff *= user->vof[i][k];
			#endif

			val=0;
			
			val -= ucont[k][j][i].x;

			if(i==1 && i_periodic) val += ucont[k][j][mx-2].x;
			else if(i==1 && ii_periodic)  val += ucont[k][j][-2].x;
			else val += ucont[k][j][i-1].x;

			val -= ucont[k][j][i].y;

			if(j==1 && j_periodic) val += ucont[k][my-2][i].y;
			else if(j==1 && jj_periodic) val += ucont[k][-2][i].y;
			else val += ucont[k][j-1][i].y;

			val -= ucont[k][j][i].z;

			if(k==1 && k_periodic) val += ucont[mz-2][j][i].z;
			else if(k==1 && kk_periodic) val += ucont[-2][j][i].z;
			else val += ucont[k-1][j][i].z;
			
			
		/*
			if(  (i==1 && j==2) ) {
				printf("%d %d %d, my_idx = %d, val=%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", i, j, k, 
					(int)gid[k][j][i], ucont[k][j][i-1].x, ucont[k][j][i].x, ucont[k][j-1][i].y, ucont[k][j][i].y, ucont[k-1][j][i].z, ucont[k][j][i].z);
				printf("\n");
			}*/
		
			val *=  -1.0 / dt * user->st * coeff;
			if(poisson==-1) val *= aj[k][j][i];
	
			lcount++;
		}
		VecSetValue(B, (int)gid[k][j][i], val, INSERT_VALUES);
		
		
	}
  //exit(0);
	VecAssemblyBegin(B);
	VecAssemblyEnd(B);

	
  
#ifndef DIRICHLET  
	
	
	double sum, sum1;
	VecSum(B, &sum);
	/*
	int N;
	VecGetSize(B, &N);
	sum  = sum/(-1.0*N);
	VecShift(B,sum);
	user->multinullspace = PETSC_FALSE;
	*/
	MPI_Allreduce( &lcount, &user->rhs_count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
	
	double val = -sum/(double) (user->rhs_count);	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]<0.1) VecSetValue(B, (int)gid[k][j][i], val, ADD_VALUES);
	}
	
	VecAssemblyBegin(B);
	VecAssemblyEnd(B);
	VecSum(B, &sum1);
#endif 
  
	PetscPrintf(PETSC_COMM_WORLD, "Summation RHS %e %e\n", sum, sum1); 
	
	DAVecRestoreArray(user->da, user->Gid, &gid);
	DAVecRestoreArray(user->fda, user->lUcont, &ucont);
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, user->lAj, &aj);
	
  

  return 0;
}


PetscErrorCode FullyBlocked(UserCtx *user)
{
  DA da = user->da;
  Vec nNvert;
  DALocalInfo info = user->info;
/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;

  PetscInt *KSKE = user->KSKE;
  PetscReal ***nvert;
  PetscTruth *Blocked;

  DACreateNaturalVector(da, &nNvert);
  DAGlobalToNaturalBegin(da, user->Nvert, INSERT_VALUES, nNvert);
  DAGlobalToNaturalEnd(da, user->Nvert, INSERT_VALUES, nNvert);

  VecScatter ctx;
  Vec Zvert;
  VecScatterCreateToZero(nNvert, &ctx, &Zvert);

  VecScatterBegin(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);

  VecScatterDestroy(ctx);
  VecDestroy(nNvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {

    VecGetArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    PetscMalloc(mx*my*sizeof(PetscTruth), &Blocked);
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	Blocked[j*mx+i] = PETSC_FALSE;
	for (k=0; k<mz; k++) {
	  if (nvert[k][j][i] > 0.1) {
	    if (!Blocked[j*mx+i]) {
	      KSKE[2*(j*mx+i)] = k;
	      Blocked[j*mx+i] = PETSC_TRUE;
	    }
	    else {
	      KSKE[2*(j*mx+i)] = PetscMin(KSKE[2*(j*mx+i)], k);
	    }
	  }
	}
      }
    }


    user->multinullspace = PETSC_TRUE;
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	if (!Blocked[j*mx+i]) {
	  user->multinullspace = PETSC_FALSE;
	  break;
	}
      }
    }
    PetscFree(Blocked);
    VecRestoreArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);

    }
  }
  else {
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }

/*   DACreateNaturalVector(da, &nNvert); */
/*   VecDestroy(nNvert); */

  VecDestroy(Zvert);
  return 0;
}

PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c)
{
  DALocalInfo	info = user_c->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal ***nvert, ***nvert_h;
  DAVecGetArray(user_h->da, user_h->lNvert, &nvert_h);
  DAVecGetArray(user_c->da, user_c->Nvert, &nvert);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  if (*(user_c->isc)) ia = 0;
  else ia = 1;
  if (*(user_c->jsc)) ja = 0;
  else ja = 1;
  if (*(user_c->ksc)) ka = 0;
  else ka = 1;
  VecSet(user_c->Nvert, 0.);
  if (user_c->thislevel > 0) {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
				for (i=lxs; i<lxe; i++) {
					GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
					if (nvert_h[kh   ][jh   ][ih   ] *
							nvert_h[kh   ][jh   ][ih-ia] *
							nvert_h[kh   ][jh-ja][ih   ] *
							nvert_h[kh-ka][jh   ][ih   ] *
							nvert_h[kh   ][jh-ja][ih-ia] *
							nvert_h[kh-ka][jh   ][ih-ia] *
							nvert_h[kh-ka][jh-ja][ih   ] *
							nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
						nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
					}
				}
      }
    }
  }
  else {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
				for (i=lxs; i<lxe; i++) {
					GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
					if (nvert_h[kh   ][jh   ][ih   ] *
							nvert_h[kh   ][jh   ][ih-ia] *
							nvert_h[kh   ][jh-ja][ih   ] *
							nvert_h[kh-ka][jh   ][ih   ] *
							nvert_h[kh   ][jh-ja][ih-ia] *
							nvert_h[kh-ka][jh   ][ih-ia] *
							nvert_h[kh-ka][jh-ja][ih   ] *
							nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
						nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
					}
				}
      }
    }
  }
  DAVecRestoreArray(user_h->da, user_h->lNvert, &nvert_h);
  DAVecRestoreArray(user_c->da, user_c->Nvert, &nvert);
  DAGlobalToLocalBegin(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DAGlobalToLocalEnd(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DAVecGetArray(user_c->da, user_c->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i] < 0.1) {
					if (nvert[k][j][i+1] + nvert[k][j][i-1] > 1.1 &&
							nvert[k][j+1][i] + nvert[k][j-1][i] > 1.1 &&
							nvert[k+1][j][i] + nvert[k-1][j][i] > 1.1) {
						nvert[k][j][i] = 1.;
					}
				}
      }
    }
  }
  DAVecRestoreArray(user_c->da, user_c->lNvert, &nvert);
  DALocalToGlobal(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert);
  return 0;
}

PetscErrorCode MyKSPMonitor1(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
  //PetscPrintf(PETSC_COMM_WORLD,"     (%D) KSP Residual norm %14.12e \n",n,rnorm);
	KSPMonitorTrueResidualNorm(ksp, n, rnorm, dummy);
	return 0;
}

// Seokkoo Kang, August 2008
// Poisson solver based on algebraic multigrid
KSP ksp_amg[100];
int was_ksp_set=0;
int was_lhs_set=0;

PetscErrorCode PoissonSolver_MG(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo)
{
	PetscInt l;
	PetscInt bi;

	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	
	MGCtx *mgctx = usermg->mgctx;
	UserCtx	*user;

	l = usermg->mglevels-1;
	user = mgctx[l].user;
  
	extern PetscInt movefsi, rotatefsi;
	if(movefsi || rotatefsi) {
		for (bi=0; bi<block_number; bi++) PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp);
	}
	else if(!was_lhs_set) {
		was_lhs_set=1;
		
		for (bi=0; bi<block_number; bi++) {
			#ifdef DIRICHLET
			VecDuplicate(user[bi].P, &user[bi].rhsD);
			PoissonLHSNew_seokkoo(&(user[bi]), ibm, user[bi].ibm_intp);
			//PoissonLHSNew_seokkoo2(&(user[bi]), ibm, user[bi].ibm_intp);
			#else
			PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp);
			#endif
		}
		
	}
	
	for (bi=0; bi<block_number; bi++) {
		user[bi].multinullspace = PETSC_FALSE;
		//VecDuplicate(user[bi].P, &user[bi].B);
		VecDuplicate(user[bi].Phi2, &user[bi].B2);
		
		//if(poisson_threshold<1)
		//double ibm_Area, ibm_Flux;
		//VolumeFlux(&user[bi], user[bi].Ucont, &ibm_Flux, &ibm_Area);
		//if(immersed==2) VolumeFlux_zero(&user[bi], user[bi].Ucont, &ibm_Flux, &ibm_Area);
		//    #ifdef PRIMITIVE_PRESSURE
		//PoissonRHS(&(user[bi]), user[bi].B);
		
		PoissonRHS2(&(user[bi]), user[bi].B2);
	}
	

  l = usermg->mglevels-1;
  user = mgctx[l].user;
  
  for (bi=0; bi<block_number; bi++) {
	if(was_ksp_set != 10) {
		KSPCreate(PETSC_COMM_WORLD, &ksp_amg[bi]);
		//KSPSetNullSpace(ksp_amg[bi], user[bi].nullsp);
		PC pc;
		KSPGetPC(ksp_amg[bi],&pc);
		PCSetType(pc, PCHYPRE);
		
		PCHYPRESetType(pc, "boomeramg");
		//PCHYPRESetType(pc, "parasails");	// works
		//PCHYPRESetType(pc, "pilut");
		
		PetscPrintf(PETSC_COMM_WORLD, "KSPSetType...\n");
		MatStructure flag;
		if(movefsi || rotatefsi) flag=DIFFERENT_NONZERO_PATTERN;
		else flag=SAME_NONZERO_PATTERN;
		
		KSPSetOperators(ksp_amg[bi], user[bi].A, user[bi].A, flag);
		KSPSetUp(ksp_amg[bi]);
		
		//
		
		KSPSetType(ksp_amg[bi],KSPGMRES);
		//KSPSetType(ksp_amg[bi],KSPCG);
		//KSPSetNormType(ksp_amg[bi], KSP_NORM_UNPRECONDITIONED);
		
		//KSPSetType(ksp_amg[bi],KSPRICHARDSON);	// To use boomeramg without PETSC Kryliv solver, # V-cycles are determined by KSPMAXIT
		KSPGMRESSetRestart(ksp_amg[bi],51);
		
		KSPMonitorSet(ksp_amg[bi],MyKSPMonitor1,PETSC_NULL,0);
		
		int size;
		MPI_Comm_size(PETSC_COMM_WORLD, &size);
		
		//PetscOptionsInsertString("-pc_hypre_boomeramg_relax_type_all SOR/Jacobi");
		PetscOptionsInsertString("-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_relax_weight_all 0");	// 07.10.2009
		//PetscOptionsInsertString("-pc_hypre_boomeramg_max_iter 3");	// # cycles
		//PetscOptionsInsertString("-pc_hypre_boomeramg_cycle_type W "); 
		PetscOptionsInsertString("-pc_hypre_boomeramg_print_statistics 1");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_nodal_relaxation 1");
		
		//PetscOptionsInsertString("-pc_hypre_boomeramg_agg_num_paths 10");
		PetscOptionsInsertString("-pc_hypre_boomeramg_max_levels 15");
		
		PetscOptionsInsertString("-pc_hypre_euclid_levels 5");	// Number of levels of fill ILU(k)
		//PetscOptionsInsertString("-pc_hypre_euclid_bj");	// Use block Jacobi, crash why?
		PetscOptionsInsertString("-pc_hypre_euclid_print_statistics");

		PetscOptionsInsertString("-pc_hypre_pilut_maxiter 5");
		PetscOptionsInsertString("-pc_hypre_pilut_tol 1.e-5");
		
		PetscOptionsInsertString("-pc_hypre_parasails_nlevels 1");
		PetscOptionsInsertString("-pc_hypre_parasails_thresh 0.1");
		PetscOptionsInsertString("-pc_hypre_parasails_filter -0.9");	// automatic
		PetscOptionsInsertString("-pc_hypre_parasails_sym SPD");
		PetscOptionsInsertString("-pc_hypre_parasails_reuse 1");
		PetscOptionsInsertString("-pc_hypre_parasails_logging 1");
		PetscOptionsInsertString("-pc_hypre_parasails_loadbal 1");
		/*
		if(size<=84) {
			PetscOptionsInsertString("-pc_hypre_boomeramg_strong_threshold 0.25");
			PetscOptionsInsertString("-pc_hypre_boomeramg_coarsen_type Falgout");	// HMIS, PMIS, modifiedRuge-Stueben, Ruge-Stueben, Falgout, CLJP
		}
		else */{
			//PetscOptionsInsertString("-pc_hypre_boomeramg_strong_threshold 0.1");
			PetscOptionsInsertString("-pc_hypre_boomeramg_strong_threshold 0.08");	//0.12 is good
			
			PetscOptionsInsertString("-pc_hypre_boomeramg_coarsen_type PMIS");	// HMIS, PMIS, modifiedRuge-Stueben, Ruge-Stueben, Falgout, CLJP
			//PetscOptionsInsertString("-pc_hypre_boomeramg_truncfactor 0.25");
		}
		// HMIS is better than Falgout when CPU>100 'cause uses less levels
		
		
		
		if(user->reduced_p_size<22300000) {
			PetscOptionsInsertString("-pc_hypre_boomeramg_P_max 5");
			PetscOptionsInsertString("-pc_hypre_boomeramg_interp_type ext+i");
		}
		else {
			PetscOptionsInsertString("-pc_hypre_boomeramg_interp_type FF1");
			PetscOptionsInsertString("-pc_hypre_boomeramg_agg_nl 1");
		}// FF1 ext+i
		
		
		//PetscOptionsInsertString("-pc_hypre_boomeramg_cycle_type V");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_truncfactor");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_grid_sweeps_coarse 5");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_grid_sweeps_up 3");
		//PetscOptionsInsertString("-pc_hypre_boomeramg_grid_sweeps_down 3");
		
		PCSetOperators(pc, user[bi].A, user[bi].A, flag);
		PCSetFromOptions(pc);
		PCSetUp(pc);
		was_ksp_set = 10;
		
		KSPSetInitialGuessNonzero(ksp_amg[bi], PETSC_TRUE);	// put guess into Phi2
		//KSPSetInitialGuessKnoll(ksp_amg[bi], PETSC_TRUE);
		
		//Create_Hypre_Solver();	// 07.09.2009 Seokkoo
		//Create_Hypre_P_Matrix_Vector(user);
	}
	PetscBarrier(PETSC_NULL);
	
	extern PetscInt tistart;
	extern double poisson_tol;
	
	if(ti==tistart) KSPSetTolerances(ksp_amg[bi], 1.e-11, 1.e-12, PETSC_DEFAULT, 50);
	else if(ti==tistart+1) KSPSetTolerances(ksp_amg[bi], 1.e-20, poisson_tol, PETSC_DEFAULT, 50);
	//KSPSolve(ksp_amg[bi], user[bi].B, user[bi].Phi);
		
	KSPSolve(ksp_amg[bi], user[bi].B2, user[bi].Phi2);
	
	#ifndef DIRICHLET
	int N;
	double sum0;
	VecGetSize(user[bi].Phi2,&N);
	VecSum(user[bi].Phi2,&sum0);
	sum0  = sum0/(-1.0*N);
	VecShift(user[bi].Phi2,sum0);
	/*
	user->multinullspace = PETSC_FALSE;
	PoissonNullSpaceFunction(user[bi].Phi2, user);
	*/
	#endif
	
	Convert_Phi2_Phi(&user[bi]);
	
	//
	
	#ifndef DIRICHLET
	/*
	int N;
	double sum0;
	VecGetSize(user[bi].Phi,&N);
	VecSum(user[bi].Phi,&sum0);
	sum0  = sum0/(-1.0*N);
	VecShift(user[bi].Phi,sum0);
        user->multinullspace = PETSC_FALSE;
	PoissonNullSpaceFunction(user[bi].Phi, user);
	*/
	#endif
	#ifdef DIRICHLET
		PetscPrintf(PETSC_COMM_WORLD, "Dirichlet BC !!! \n"); 
	#endif
  }
	
  
  for (bi=0; bi<block_number; bi++) {
    DAGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
    DAGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
  }

  for (bi=0; bi<block_number; bi++) {
    //VecDestroy(mgctx[usermg->mglevels-1].user[bi].B);
    VecDestroy(mgctx[usermg->mglevels-1].user[bi].B2);
  }
  
	PetscGetTime(&te);
	cput=te-ts;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {
		FILE *f;
		char filen[80];
		sprintf(filen, "%s/Converge_dU", path);
		f = fopen(filen, "a");
		//PetscFPrintf(PETSC_COMM_WORLD, f, "%d(poisson)  %le %.2e(s) %d(iter)\n", ti, poisson::Poisson_Solver[bi]->TrueResidual(), cput, poisson::Poisson_Solver[bi]->NumIters());
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d(poisson)  %.2e(s)", ti, cput);
		fclose(f);
	}
	
  return 0;
}

PetscErrorCode PoissonSolver_MG_original(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo)
{
  PetscInt l;
  PetscInt bi;

	PetscReal ts,te,cput;
	
	
  MGCtx *mgctx = usermg->mgctx;

  KSP	*mgksp, subksp, csksp;
  PC	*mgpc, subpc;
  UserCtx	*user;

  PetscInt	m_c, m_f, M_c, M_f;

  PetscMalloc(block_number*sizeof(KSP), &mgksp);
  PetscMalloc(block_number*sizeof(PC), &mgpc);

  //PetscPrintf(PETSC_COMM_WORLD, "TEST         !\n");  
  //  MyNFaceInit(usermg);


  if (immersed) {
    for (l=usermg->mglevels-1; l>0; l--) {
      for (bi=0; bi<block_number; bi++) {
	mgctx[l].user[bi].multinullspace = PETSC_FALSE;
	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]);
      }
    }
    /* At the corsest level, check whether the grid is separated into sevearal 
       blockes by the immersed body */
    l = 0;
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      PetscMalloc(user[bi].info.mx*user[bi].info.my*2*sizeof(PetscInt), &user[bi].KSKE);
      FullyBlocked(&user[bi]);

    }
  }

  if( ti==tistart || movefsi || rotatefsi || levelset )
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      if (l==usermg->mglevels-1)
	PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp);
      else
	PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp);
    }
  }
  
  l = usermg->mglevels-1;
  user = mgctx[l].user;
  for (bi=0; bi<block_number; bi++) {
    VecDuplicate(user[bi].P, &user[bi].B);
    PetscReal ibm_Flux, ibm_Area;
    VolumeFlux(&user[bi], user[bi].Ucont, &ibm_Flux, &ibm_Area);

    PoissonRHS2(&(user[bi]), user[bi].B);
/*     PetscReal sum; */
/*     Vec Temp, Aj; */
/*     VecDuplicate(user[bi].B, &Temp); */
/*     VecDuplicate(user[bi].B, &Aj); */
/*     DALocalToGlobal(user[bi].da, user[bi].lAj, INSERT_VALUES, Aj); */
    
/*     VecPointwiseDivide(Temp, user[bi].B, Aj); */

/*     VecSum(Temp, &sum); */
/*     VecDestroy(Temp); */
/*     VecDestroy(Aj); */
		    
/*     PetscPrintf(PETSC_COMM_WORLD, "Summation RHS %e\n", sum); */
  }

  

/*   l = usermg->mglevels-1; */
/*   user = mgctx[l].user; */
/*   for (bi=0; bi<block_number; bi++) { */
/*     KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]); */
/*     KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN); */
/*     KSPSetFromOptions(mgksp[bi]); */
/*     KSPSetUp(mgksp[bi]); */
/*     KSPSetTolerances(mgksp[bi], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10); */

/*     KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); */

/*     KSPDestroy(mgksp[bi]); */
/*   } */
  for (bi=0; bi<block_number; bi++) {
    /* Create ksp for Multigrid */
    KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]);
    KSPAppendOptionsPrefix(mgksp[bi], "ps_");
    
/*     KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); */
    
    /* Use multigrid as preconditioner */
    KSPGetPC(mgksp[bi], &mgpc[bi]);	
    PCSetType(mgpc[bi], PCMG);	
    /* Setup MG levels from usergm->mglevels */
    PCMGSetLevels(mgpc[bi], usermg->mglevels, PETSC_NULL);	
    /* V cycle */
    //PCMGSetCycles(mgpc[bi], 1);
    PCMGSetCycleType(mgpc[bi],PC_MG_CYCLE_V);	// kangsk, v 2.3.3	
    
    PCMGSetType(mgpc[bi], PC_MG_MULTIPLICATIVE);	
/*     PCMGSetType(mgpc[bi], PC_MG_FULL); */

    /* Create Restriction and Interpolate schemes
       This is needed for all levels other than the coarsest one */
    for (l=usermg->mglevels-1; l>0; l--) {
      user = mgctx[l].user;
      m_c = (usermg->mgctx[l-1].user[bi].info.xm *
	     usermg->mgctx[l-1].user[bi].info.ym *
	     usermg->mgctx[l-1].user[bi].info.zm);

      m_f = (usermg->mgctx[l].user[bi].info.xm *
	     usermg->mgctx[l].user[bi].info.ym *
	     usermg->mgctx[l].user[bi].info.zm);

      M_c = (usermg->mgctx[l-1].user[bi].info.mx *
	     usermg->mgctx[l-1].user[bi].info.my *
	     usermg->mgctx[l-1].user[bi].info.mz);

      M_f = (usermg->mgctx[l].user[bi].info.mx *
	     usermg->mgctx[l].user[bi].info.my *
	     usermg->mgctx[l].user[bi].info.mz);

      MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)&mgctx[l-1].user[bi], &user[bi].MR);
      MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)&user[bi], &user[bi].MP);

      PCMGSetRestriction(mgpc[bi], l, user[bi].MR);
      PCMGSetInterpolation(mgpc[bi], l, user[bi].MP);

      /* Use subroutine MyRestriction and MyInterpolation for
	 Mat * Vec operation */
      MatShellSetOperation(user[bi].MR, MATOP_MULT, (void(*)(void))MyRestriction);
      MatShellSetOperation(user[bi].MP, MATOP_MULT, (void(*)(void))MyInterpolation);

      MatShellSetOperation(user[bi].MR, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);
      MatShellSetOperation(user[bi].MP, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);

    }
  }

PetscPrintf(PETSC_COMM_WORLD, "\nAAA\n");
  
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      
      if (l) { /* Not the coarset grid level */
	PCMGGetSmoother(mgpc[bi], l, &subksp);
	/* Set the left hand side for every KSP at each grid level */
	KSPSetOperators(subksp, user[bi].A, user[bi].A,
			DIFFERENT_NONZERO_PATTERN);
      
	KSPGMRESSetRestart(subksp,20);//seokkoo 
	KSPSetType(subksp,KSPFGMRES);//seokkoo: 
	      
      	//KSPSetFromOptions(subksp);// removed by seokkoo
	KSPSetUp(subksp);

	KSPGetPC(subksp, &subpc);
/* 	PCSetType(subpc, PCASM); */
	KSP *subsubksp;
	PC subsubpc;
	PetscInt abi, nlocal;
	PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
	for (abi = 0; abi<nlocal; abi++) {
	  KSPGetPC(subsubksp[abi], &subsubpc);
	  //PCFactorSetShiftNonzero(subsubpc, 1.e-10);//removed seokkoo
		PCFactorSetShiftType(subsubpc, MAT_SHIFT_NONZERO); //seokkoo                                                                                                            
                PCFactorSetShiftAmount(subsubpc, 1.e-10); //seokkoo 
	}

	//PCFactorSetShiftNonzero(subpc, 1.e-10);//removed seokkoo
	PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO); //seokkoo                                                                                                               
	PCFactorSetShiftAmount(subpc, 1.e-10); //seokkoo 
      }
      else {  /* Coarsest grid */

	/* The default solver for the coarset grid is
	   KSPPreonly and PCLU.
	   One can choose other solvers, such as a KSP solver by changing the
	   following lines. */
	PCMGGetCoarseSolve(mgpc[bi], &subksp);
/* 	KSPSetType(subksp, KSPBCGS); */
//	KSPSetType(subksp, KSPPREONLY);//seokkoo

	KSPGetPC(subksp, &subpc);
	KSPSetOperators(subksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);

	PCSetType(subpc, PCBJACOBI);
	//KSPGMRESSetRestart(subksp,200);//seokkoo 2010.1.5
	KSPSetUp(subksp);

	KSP *subsubksp;
	PC subsubpc;
	PetscInt abi, nlocal;
	PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
	for (abi = 0; abi<nlocal; abi++) {
	  KSPGetPC(subsubksp[abi], &subsubpc);
	  //PCFactorSetShiftNonzero(subsubpc, 1.e-10);//removed seokkoo
	}

	if(usermg->mglevels>1) {
		KSPSetTolerances(subksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, 20);
		KSPGMRESSetRestart(subksp,20);//seokkoo 
		KSPSetType(subksp,KSPFGMRES);//seokkoo: 
		//KSPSetType(subksp,KSPCG);//seokkoo: 
	}
	
	

/* 	PC *subpc; */
/* 	KSPGetPC(subksp,  */
/* 	KSPGetPC(subksp, &subpc); */
/* 	PCFactorSetShiftNonzero(subpc, 1.e-9); */
/* 	PCSetUp(subpc); */
      }

      /* The Poisson equation has Neumann boundary conditions, thus
	 need to use NullSpace to solve the equation.
	 We use the subroutine PoissonNullSpaceFunction to
	 (1) reduce a constant value from the solution
	 (2) Set the solution values at boundaries and blanking nodes to
	 zero */
      MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp);
      MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction_original, &user[bi]);
      
    
      KSPSetNullSpace(subksp, user[bi].nullsp);

      PCMGSetResidual(mgpc[bi], l, PCMGDefaultResidual, user[bi].A);

      KSPSetUp(subksp);

      if (l<usermg->mglevels-1) {
	MatGetVecs(user[bi].A, &user[bi].R, PETSC_NULL);
	//	VecDuplicate(user[bi].P, &user[bi].R);
	PCMGSetRhs(mgpc[bi], l, user[bi].R);
      }

    }
  }
PetscPrintf(PETSC_COMM_WORLD, "\nBBB Poisson.c\n");
  l = usermg->mglevels-1;
  user = mgctx[l].user;
  
  for (bi=0; bi<block_number; bi++) {
    KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
    KSPSetNullSpace(mgksp[bi], user[bi].nullsp);

   
	KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); //seokkoo
	KSPSetType(mgksp[bi],KSPFGMRES);//seokkoo
    //KSPSetType(mgksp[bi],KSPPREONLY);//seokkoo
    //KSPSetType(mgksp[bi],KSPCG);//seokkoo
    KSPGMRESSetRestart(mgksp[bi],20);//seokkoo
    KSPSetTolerances(mgksp[bi], 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 200);//seokkoo
    

    KSPSetFromOptions(mgksp[bi]);
    KSPSetUp(mgksp[bi]);

    /* PetscReal norm; */
/*     VecNorm(user[bi].B, NORM_2, &norm); */
/*     PetscPrintf(PETSC_COMM_WORLD, "KSP RHS Norm %e\n", norm); */
    KSPMonitorSet(mgksp[bi], MyKSPMonitor1, PETSC_NULL, 0);
    
    PetscGetTime(&ts);
    KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi);

/*     Vec vv; */
/*     KSPBuildResidual(mgksp[bi], PETSC_NULL, PETSC_NULL, &vv); */

    /*
    if (user->main->prfield) {
      
      if (ti/user->main->prio * user->main->prio == ti) {
	Vec Temp;
	VecDuplicate(user[bi].P, &Temp);
	MatMult(user[bi].A, user[bi].Phi, Temp);
	VecAXPY(Temp, -1., user[bi].B);


	PetscViewer viewer;
	char filen[90];
	PetscOptionsClearValue("-vecload_block_size");

	sprintf(filen, "prfield%5.5d_%1.1d.dat", ti, user->this);

	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(Temp, viewer);
	//     VecLoadIntoVector(viewer, Temp);
	PetscViewerDestroy(viewer);
    
	VecDestroy(Temp);
      }
    } */


/*     MatMult(user[bi].A, user[bi].Phi, user[bi].X); */
/*     VecAXPY(user[bi].X, -1., user[bi].B); */
/*     VecCopy(vv, user[bi].Phi); */

/*     VecNorm(user[bi].X, NORM_2, &norm); */
/*     PetscPrintf(PETSC_COMM_WORLD, "Poisson Norm %e\n", norm); */

/*     PetscTruth flg = PETSC_FALSE; */
/*     PetscOptionsGetTruth(PETSC_NULL, "-after_mg", &flg, PETSC_NULL); */
/*     if (flg) { */
/*       PCSetType(mgpc[bi], PCBJACOBI); */
/*       KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); */
	 
/*       KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); */
/*     } */
  }
  
	PetscGetTime(&te);
	cput=te-ts;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {
		FILE *f;
		char filen[80];
		sprintf(filen, "Converge_dU");
		f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d(poisson)  %.2e(s)", ti, cput);
		fclose(f);
	}
		

  for (bi=0; bi<block_number; bi++) {
    DAGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
    DAGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
  }

  //  MyNFaceFinalize(usermg);

  /* Release the allocated spaces */
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      MatNullSpaceDestroy(user[bi].nullsp);

      if( movefsi || rotatefsi || levelset ) MatDestroy(user[bi].A);
      user[bi].assignedA = PETSC_FALSE;
      if (l) { /* Grid level > 0 (coarsest) */
	MatDestroy(user[bi].MR);
	MatDestroy(user[bi].MP);
      }
      else {
	PetscFree(user[bi].KSKE);
      }
    }
  }
  for (bi=0; bi<block_number; bi++) {
    //    PCDestroy(mgpc[bi]);
    
    for (l=0; l<usermg->mglevels-1; l++) {
      VecDestroy(mgctx[l].user[bi].R);
    }
    
    KSPDestroy(mgksp[bi]);
    VecDestroy(mgctx[usermg->mglevels-1].user[bi].B);
  }
  PetscFree(mgksp);
  PetscFree(mgpc);

  return 0;
}

PetscErrorCode Projection(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;
  PetscReal	***phi, ***rho, ***p;

  Cmpnts	***ucont, ***lucont;
  PetscReal	***nvert, ***level;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

	if (levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lLevelset, &level);
	}
	
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);

	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lIEta, &ieta);
	DAVecGetArray(fda, user->lIZet, &izet);

	DAVecGetArray(fda, user->lJCsi, &jcsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lJZet, &jzet);

	DAVecGetArray(fda, user->lKCsi, &kcsi);
	DAVecGetArray(fda, user->lKEta, &keta);
	DAVecGetArray(fda, user->lKZet, &kzet);

	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);

	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lPhi, &phi);
	DAVecGetArray(da, user->lP, &p);
  
	DAVecGetArray(fda, user->Ucont, &ucont);

	PetscReal dpdc, dpde, dpdz;
  
	Cmpnts ***cent;
	DAVecGetArray(fda, user->lCent, &cent);

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	      
		double coeff = time_coeff();
		double r1=1.0, r2=1.0, r3=1.0;
		
		if(levelset) {
			if(i==mx-2) r1=rho[k][j][i];
			else r1 = mean ( rho[k][j][i], rho[k][j][i+1] );
			
			if(j==my-2) r2=rho[k][j][i];
			else r2 = mean ( rho[k][j][i], rho[k][j+1][i] );
			
			if(k==mz-2) r3=rho[k][j][i];
			else r3 = mean ( rho[k][j][i], rho[k+1][j][i] );
		}
		
		if (i<mx-2 || ( (i_periodic || ii_periodic) && i==mx-2) ) {
			dpdc = phi[k][j][i+1] - phi[k][j][i];
			if( i==mx-2 && i_periodic) dpdc = phi[k][j][1] - phi[k][j][i];
			else if( i==mx-2 && ii_periodic) dpdc = phi[k][j][mx+1] - phi[k][j][i];

			dpde = 0.;
			dpdz = 0.;
			 
			// seokkoo - take into account Dirichlet pressure here !!
			if ( (j==my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i]+nvert[k][j+1][i+1] > poisson_threshold) {
				if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < poisson_threshold && (j!=1) ) {
					dpde = (phi[k][j  ][i] + phi[k][j  ][i+1] - phi[k][j-1][i] - phi[k][j-1][i+1]) * 0.5;
				}
			}
			else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > poisson_threshold) {
				if (nvert[k][j+1][i] + nvert[k][j+1][i+1] <poisson_threshold) {
					dpde = (phi[k][j+1][i] + phi[k][j+1][i+1] - phi[k][j  ][i] - phi[k][j  ][i+1]) * 0.5;
				}
			}
			else {
				dpde = (phi[k][j+1][i] + phi[k][j+1][i+1] - phi[k][j-1][i] - phi[k][j-1][i+1]) * 0.25;
			}

			if ( (k == mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j][i+1] >poisson_threshold) {
				if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < poisson_threshold && (k!=1) ) {
					dpdz = (phi[k  ][j][i] + phi[k  ][j][i+1] - phi[k-1][j][i] - phi[k-1][j][i+1]) * 0.5;
				}
			}
			else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j][i+1] >poisson_threshold) {
				if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < poisson_threshold) {
					dpdz = (phi[k+1][j][i] + phi[k+1][j][i+1] - phi[k  ][j][i] - phi[k  ][j][i+1]) * 0.5;
				}
			}
			else {
				dpdz = (phi[k+1][j][i] + phi[k+1][j][i+1] - phi[k-1][j][i] - phi[k-1][j][i+1]) * 0.25;
			}

			if ((nvert[k][j][i] + nvert[k][j][i+1])<poisson_threshold) {	// how to take account ib interface velocity?
				ucont[k][j][i].x -= 
					(dpdc * (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y + icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
					dpde * (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y + ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] + 
					dpdz * (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y + izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i]) * user->dt * user->st / coeff / r1;
			}
		}
		
		if (j<my-2 || ( (j_periodic || jj_periodic) && j==my-2) ||
		    (levelset && user->bctype[3]==4 && j==my-2)
		    ) {
			dpdc = 0.;
			dpdz = 0.;
			if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k][j+1][i+1] >poisson_threshold) {	
				if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < poisson_threshold && (i!=1) ) {
					dpdc = (phi[k][j][i  ] + phi[k][j+1][i  ] - phi[k][j][i-1] - phi[k][j+1][i-1]) * 0.5;
				}
			}
			else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > poisson_threshold) {
				if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < poisson_threshold) {
					dpdc = (phi[k][j][i+1] + phi[k][j+1][i+1] - phi[k][j][i  ] - phi[k][j+1][i  ]) * 0.5;
				}
			}
			else {
				dpdc = (phi[k][j][i+1] + phi[k][j+1][i+1] - phi[k][j][i-1] - phi[k][j+1][i-1]) * 0.25;
			}

			if ( levelset && user->bctype[3]==4 && j==my-2 ) dpde = phi[k][j][i] - phi[k][j-1][i];
			else if( j==my-2 && j_periodic) dpde = phi[k][1][i] - phi[k][j][i];
			else if( j==my-2 && jj_periodic) dpde = phi[k][my+1][i] - phi[k][j][i];
			else {
			  dpde = phi[k][j+1][i] - phi[k][j][i];
			}

			if ( (k == mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i] + nvert[k+1][j+1][i] > poisson_threshold) {
				if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < poisson_threshold && (k!=1) ) {
					dpdz = (phi[k  ][j][i] + phi[k  ][j+1][i] - phi[k-1][j][i] - phi[k-1][j+1][i]) * 0.5;
				}
			}
			else if ( (k == 1 && !k_periodic && !kk_periodic) || nvert[k-1][j][i] + nvert[k-1][j+1][i] > poisson_threshold) {
				if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < poisson_threshold) {
					dpdz = (phi[k+1][j][i] + phi[k+1][j+1][i] - phi[k  ][j][i] - phi[k  ][j+1][i]) * 0.5;
				}
			}
			else {
				dpdz = (phi[k+1][j][i] + phi[k+1][j+1][i] - phi[k-1][j][i] - phi[k-1][j+1][i]) * 0.25;
			}
			
			
			if ((nvert[k][j][i] + nvert[k][j+1][i])<poisson_threshold || (int)(nvert[k][j][i]+nvert[k][j+1][i])==5) {//seokkoo
				ucont[k][j][i].y -=
					(dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
					dpde * (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
					dpdz * (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i]) * user->dt * user->st / coeff / r2;
			}
		}
		
		if (	k < mz-2 || 
			( (kk_periodic || k_periodic) && k==mz-2 )  ||
			(levelset && user->bctype[5]==4 && k==mz-2) ) {
			dpdc = 0.;
			dpde = 0.;
			if ( (i == mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1] + nvert[k+1][j][i+1] > poisson_threshold) {
				if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < poisson_threshold && (i!=1) ) {
					dpdc = (phi[k][j][i  ] + phi[k+1][j][i  ] - phi[k][j][i-1] - phi[k+1][j][i-1]) * 0.5;
				}
			}
			else if ( (i == 1 && !i_periodic && !ii_periodic) || nvert[k][j][i-1] + nvert[k+1][j][i-1] > poisson_threshold) {
				if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < poisson_threshold) {
					dpdc = (phi[k][j][i+1] + phi[k+1][j][i+1] - phi[k][j][i  ] - phi[k+1][j][i  ]) * 0.5;
				}
			}
			else {
				dpdc = (phi[k][j][i+1] + phi[k+1][j][i+1] - phi[k][j][i-1] - phi[k+1][j][i-1]) * 0.25;
			}
		  
			if ( (j == my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i] + nvert[k+1][j+1][i] > poisson_threshold) {
				if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < poisson_threshold && (j!=1) ) {
					dpde = (phi[k][j  ][i] + phi[k+1][j  ][i] - phi[k][j-1][i] - phi[k+1][j-1][i]) * 0.5;
				}
			}
			else if ( (j == 1 && !j_periodic && !jj_periodic) || nvert[k][j-1][i] + nvert[k+1][j-1][i] > poisson_threshold) {
				if (nvert[k][j+1][i] + nvert[k+1][j+1][i] <poisson_threshold) {
					dpde = (phi[k][j+1][i] + phi[k+1][j+1][i] - phi[k][j  ][i] - phi[k+1][j  ][i]) * 0.5;
				}
			}
			else {
				dpde = (phi[k][j+1][i] + phi[k+1][j+1][i] - phi[k][j-1][i] - phi[k+1][j-1][i]) * 0.25;
			}
			
			
			if( k==mz-2 && k_periodic) dpdz = phi[1][j][i] - phi[k][j][i];
			else if( k==mz-2 && kk_periodic) dpdz = phi[mz+1][j][i] - phi[k][j][i];
			else if (k==mz-2 && levelset && user->bctype[5]==4) {
				double pBC_old = 1.5 * p[k][j][i] - 0.5 * p[k-1][j][i];
				double pBC_new = p[k+1][j][i];
				double phiBC_new = pBC_new - pBC_old;
				
				//dpdz = -	(0 - phi[k][j][i])/0.5;
				dpdz = phi[k][j][i] - phi[k-1][j][i];
			}
			else dpdz = phi[k+1][j][i] - phi[k][j][i];
			
			if ((nvert[k][j][i] + nvert[k+1][j][i])<poisson_threshold) {
				/*
				if(k_periodic || kk_periodic) {
					double dz = 1./ kaj[k][j][i] / sqrt(kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);

					if ( fabs(mean_pressure_gradient)>1.e-5) {
						dpdz += mean_pressure_gradient * dz;
					}
					else if(inletprofile!=17) {
						extern double inlet_flux;
						dpdz += dz * (user->mean_k_flux-inlet_flux) / user->dt / user->mean_k_area * coeff;//sqrt(kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);
					}
				}
				*/
				
				ucont[k][j][i].z -=
						(dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
						dpde * (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
						dpdz * (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) * user->dt *user->st / coeff / r3;
				
				if(levelset && k==1 && user->bctype[4]==5000 && level[k][j][i]<0) {
					dpdc=0, dpde=0;
					r3 = rho[k][j][i];
					ucont[0][j][i].z -=
						(dpdz * (kzet[0][j][i].x * kzet[0][j][i].x + kzet[0][j][i].y * kzet[0][j][i].y + kzet[0][j][i].z * kzet[0][j][i].z) * kaj[0][j][i]) * user->dt *user->st / coeff / r3;
				}
			}
		}
	}
    
	if (levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lLevelset, &level);
	}
  
	DAVecRestoreArray(fda, user->lCent, &cent);
  
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);

	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);

	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);

	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);

	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);
	DAVecRestoreArray(da, user->lPhi, &phi);
	DAVecRestoreArray(da, user->lP, &p);
	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(da, user->lNvert, &nvert);

	DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	
	if(periodic) {
	  DAVecGetArray(user->fda, user->Ucont, &ucont);
	  DAVecGetArray(user->fda, user->lUcont, &lucont);
	
	  for (k=zs; k<ze; k++)
	    for (j=ys; j<ye; j++)
	      for (i=xs; i<xe; i++) {
		int i_flag=0, j_flag=0, k_flag=0;
		int a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, i_flag=1;
		else if(i_periodic && i==mx-1) a=1, i_flag=1;
		
		if(j_periodic && j==0) b=my-2, j_flag=1;
		else if(j_periodic && j==my-1) b=1, j_flag=1;
		
		if(k_periodic && k==0) c=mz-2, k_flag=1;
		else if(k_periodic && k==mz-1) c=1, k_flag=1;
		
		if(ii_periodic && i==0) a=-2, i_flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, i_flag=1;
		
		if(jj_periodic && j==0) b=-2, j_flag=1;
		else if(jj_periodic && j==my-1) b=my+1, j_flag=1;
		
		if(kk_periodic && k==0) c=-2, k_flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, k_flag=1;
		
		if(i_flag) lucont[k][j][i].x = lucont[c][b][a].x;
		if(j_flag) lucont[k][j][i].y = lucont[c][b][a].y;
		if(k_flag) lucont[k][j][i].z = lucont[c][b][a].z;
		
		ucont[k][j][i] = lucont[k][j][i];
	      }
	
	  DAVecRestoreArray(user->fda, user->Ucont, &ucont);
	  DAVecRestoreArray(user->fda, user->lUcont, &lucont);
	
	  DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); //101229
	  DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	}
	/*
	
	DAVecGetArray(fda, user->lUcont, &ucont);
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i_periodic && i==0) ucont[k][j][0].x=ucont[k][j][mx-2].x;
		if(j_periodic && j==0) ucont[k][0][i].y=ucont[k][my-2][i].y;
		if(k_periodic && k==0) ucont[0][j][i].z=ucont[mz-2][j][i].z;
		
		if(ii_periodic && i==0) 	ucont[k][j][0].x=ucont[k][j][-2].x;
		if(jj_periodic && j==0) 	ucont[k][0][i].y=ucont[k][-2][i].y;
		if(kk_periodic && k==0) 	ucont[0][j][i].z=ucont[-2][j][i].z;
	}
	DAVecRestoreArray(fda, user->lUcont, &ucont);
	
	DALocalToGlobal(fda, user->lUcont, INSERT_VALUES, user->Ucont);
	DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  */
	
	//PetscBarrier(PETSC_NULL);

	Contra2Cart(user);
	
  //extern GhostNodeVelocity2(UserCtx *);
  //GhostNodeVelocity(user);
  //GhostNodeVelocity2(user);
  return(0);
}
/*
void SetDirichletPressure(UserCtx *user)
{
	if(!levelset) return;
	
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;

	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt	lxs, lys, lzs, lxe, lye, lze;
	PetscInt	i, j, k;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	PetscReal ***p, ***lp, ***level;
	Cmpnts ***cent;
	
	DAVecGetArray(fda, user->lCent, &cent); 
	DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, user->lP, &lp);
	DAVecGetArray(da, user->P, &p);
	
	if(ze==mz && levelset && user->bctype[5]==4) {
		k=mz-1;
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
			double L = level[k][j][i];
			
			if(L<0) {
				if(inlet_y_flag) {
					p[k][j][i] = rho_air*gravity_y*(inlet_y - y);
				}
				else if(inlet_z_flag) {
					p[k][j][i] = rho_air*gravity_z*(inlet_z - z);
				}
			}
			else {
				if(inlet_y_flag) {
					p[k][j][i] = rho_water*gravity_y*(inlet_y - y);
				}
				else if(inlet_z_flag) {
					p[k][j][i] = rho_water*gravity_z*(inlet_z - z);
				}
			}
			lp[k][j][i] = p[k][j][i];
		}
	}
	
	DAVecRestoreArray(fda, user->lCent, &cent); 
	DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, user->lP, &lp);
	DAVecRestoreArray(da, user->P, &p);
}
*/

PetscErrorCode UpdatePressure(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;

	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt	lxs, lys, lzs, lxe, lye, lze;
	PetscInt	i, j, k;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	Vec coords;	// seokkoo
	Cmpnts ***ucont, ***coor, ***ucat;
	Cmpnts ***csi, ***eta, ***zet;
	
	PetscReal ***p, ***phi, ***lp, ***lphi, ***aj, ***nvert, ***lnu_t;

	DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
    
	
	DAVecGetArray(da, user->P, &p);
	DAVecGetArray(da, user->Phi, &phi);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, user->lUcont, &ucont); 
	DAVecGetArray(fda, user->lUcat, &ucat); 
	DAVecGetArray(da, user->lAj, &aj); 
   
	DAGetGhostedCoordinates(da, &coords);
	DAVecGetArray(fda, coords, &coor);

	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
  
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		#ifdef	PRIMITIVE_PRESSURE
		p[k][j][i] = phi[k][j][i];
		
		#ifdef DIRICHLET
		double upper_y = 0.25 * ( coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y + coor[k  ][j  ][i-1].y +coor[k-1][j  ][i-1].y );
		double lower_y  = 0.25 * ( coor[k  ][j-1][i  ].y + coor[k-1][j-1][i  ].y + coor[k  ][j-1][i-1].y + coor[k-1][j-1][i-1].y );
		/*
			p (cent[k][j+1][i].y - upper_y ) + p_j+1 ( upper_y - cent[k][j][i].y ) = pD
			p ( lower_y - cent[k][j-1][i].y ) + p_j-1 ( cent[k][j][i].y - lower_y ) = pD
		*/
		      /*
		if( freesurface && (int)nvert[k][j][i]==5 && (int)nvert[k][j-1][i]==0 ) {
			double AB = cent[k][j][i].y - cent[k][j-1][i].y;
			double AS = free_surface_y[i][k] - cent[k][j-1][i].y;
			assert(AS>=0);
			
			double D = cent[k][j][i].y - cent[k][j-1][i].y;
			p[k][j][i] = (0 - phi[k][j-1][i]) * ( cent[k][j][i].y - lower_y ) / ( lower_y - cent[k][j-1][i].y );
		}
		*/
		#endif
		
		#else
		if( nvert[k][j][i]<0.1) p[k][j][i] += phi[k][j][i];
		
		if( nvert[k][j][i]>1.1 ) p[k][j][i] = 0;
		//p[k][j][i] -= (ucont[k][j][i].x - ucont[k][j][i-1].x + ucont[k][j][i].y - ucont[k][j-1][i].y + ucont[k][j][i].z - ucont[k-1][j][i].z ) * aj[k][j][i] / user->ren; 
		/*
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double ajc = aj[k][j][i];
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		p[k][j][i] -= (du_dx + dv_dy + dw_dz) / user->ren; 
		*/
		#endif
	}
    
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	
	DAVecRestoreArray(da, user->Phi, &phi);
	DAVecRestoreArray(da, user->P, &p);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, user->lUcont, &ucont); 
	DAVecRestoreArray(fda, user->lUcat, &ucat); 
	DAVecRestoreArray(da, user->lAj, &aj); 
	DAVecRestoreArray(fda, coords, &coor);
  
	DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
	
	DAVecGetArray(da, user->lP, &lp);
	DAVecGetArray(da, user->lPhi, &lphi);
	DAVecGetArray(da, user->P, &p);
	DAVecGetArray(da, user->Phi, &phi);
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		if(flag) {
			p[k][j][i] = lp[c][b][a];
			phi[k][j][i] = lphi[c][b][a];
		}
		
	}
	DAVecRestoreArray(da, user->lP, &lp);
	DAVecRestoreArray(da, user->lPhi, &lphi);
	DAVecRestoreArray(da, user->P, &p);
	DAVecRestoreArray(da, user->Phi, &phi);
	/*
	DALocalToGlobal(da, user->lP, INSERT_VALUES, user->P);
	DALocalToGlobal(da, user->lPhi, INSERT_VALUES, user->Phi);
	*/
	DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
		
	DAGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
	DAGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
	
	
  return 0;
}

PetscErrorCode VolumeFlux(UserCtx *user, Vec lUcor, PetscReal *ibm_Flux, PetscReal *ibm_Area)
{
  DA	da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal ***nvert;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DAVecGetArray(fda, lUcor, &ucor);
  DAVecGetArray(fda, user->lCsi, &csi);
  DAVecGetArray(fda, user->lEta, &eta);
  DAVecGetArray(fda, user->lZet, &zet);
  DAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area;
  
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < poisson_threshold) {
	  if (nvert[k][j][i+1] > poisson_threshold && i <= mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x + csi[k][j][i].y * csi[k][j][i].y + csi[k][j][i].z * csi[k][j][i].z);
	  }
	  if (nvert[k][j+1][i] > poisson_threshold && j <= my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y +  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > poisson_threshold && k <= mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > poisson_threshold ) {
	  if (nvert[k][j][i+1] < poisson_threshold && i <= mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x + csi[k][j][i].y * csi[k][j][i].y + csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < poisson_threshold && j <= my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y + eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < poisson_threshold && k <= mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }
  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD);
  PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux %le %le\n", *ibm_Flux, *ibm_Area);

  PetscReal correction;

  if (*ibm_Area > 1.e-15) {
    correction = *ibm_Flux / *ibm_Area;
  }
  else {
    correction = 0;
  }

	if(immersed!=2)
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < poisson_threshold) {
	  if (nvert[k][j][i+1] > poisson_threshold && i <= mx-2) {
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x + csi[k][j][i].y * csi[k][j][i].y + csi[k][j][i].z * csi[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].x=0;
	  }
	  if (nvert[k][j+1][i] > poisson_threshold && j <= my-2) {
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y + eta[k][j][i].z * eta[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].y=0;
	  }
	  if (nvert[k+1][j][i] > poisson_threshold && k <= mz-2) {
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].z=0;
	  }
	}

	if (nvert[k][j][i] > poisson_threshold ) {
	  if (nvert[k][j][i+1] < poisson_threshold && i <= mx-2) {
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x + csi[k][j][i].y * csi[k][j][i].y + csi[k][j][i].z * csi[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].x=0;
	  }
	  if (nvert[k][j+1][i] < poisson_threshold && j <= my-2) {
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y + eta[k][j][i].z * eta[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].y=0;
	  }
	  if (nvert[k+1][j][i] < poisson_threshold && k <= mz-2) {
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z) * correction;
		if(immersed==3 || !immersed) ucor[k][j][i].z=0;
	  }
	}

      }
    }
  }
  


  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < poisson_threshold) {
	  if (nvert[k][j][i+1] > poisson_threshold && i <= mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > poisson_threshold && j <= my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > poisson_threshold && k <= mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > poisson_threshold ) {
	  if (nvert[k][j][i+1] < poisson_threshold && i <= mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < poisson_threshold && j <= my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < poisson_threshold && k <= mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }
  
  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD);
  PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux22 %le %le\n", *ibm_Flux, *ibm_Area);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(fda, lUcor, &ucor);

  DAGlobalToLocalBegin(fda, lUcor, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, lUcor, INSERT_VALUES, user->lUcont);
  return 0;
}

void Add_IBMFlux_to_Outlet(UserCtx *user, PetscReal ibm_Flux)
{
	if(user->bctype[5]!=4) return;
		
	DA	da = user->da, fda = user->fda;

	DALocalInfo	info = user->info;

	PetscInt xs = info.xs, xe = info.xs + info.xm;
	PetscInt ys = info.ys, ye = info.ys + info.ym;
	PetscInt zs = info.zs, ze = info.zs + info.zm;
	PetscInt mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
	PetscInt lxs, lys, lzs, lxe, lye, lze;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
		
	PetscReal ***nvert;
	Cmpnts ***ucont, ***csi, ***eta, ***zet;
	
	DAVecGetArray(fda, user->Ucont, &ucont);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lNvert, &nvert);
  
	double AreaSum=0, lArea=0;
	if(ze==mz) {
		k = ze-1;
		for (j=lys; j<lye; j++) 
		for (i=lxs; i<lxe; i++) {
			if (nvert[k-1][j][i] < 0.1) {
				lArea += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
			}
		}
	}
	PetscGlobalSum(&lArea, &AreaSum, PETSC_COMM_WORLD);
	
	if(ze==mz) {
		k = ze-1;
		for (j=lys; j<lye; j++) 
		for (i=lxs; i<lxe; i++) {
			if (nvert[k-1][j][i] < 0.1) {
				double A = sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
				ucont[k-1][j][i].z += ibm_Flux * A / AreaSum;
			}
		}
	}

	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(da, user->lNvert, &nvert);
}


PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg)
{
  DA	da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal epsilon=1.e-8;
  PetscReal ***nvert, ibmval=1.1;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DAVecGetArray(fda, user->Ucont, &ucor);
  DAVecGetArray(fda, user->lCsi, &csi);
  DAVecGetArray(fda, user->lEta, &eta);
  DAVecGetArray(fda, user->lZet, &zet);
  DAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area, libm_Flux_abs=0., ibm_Flux_abs;
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux += ucor[k][j][i].x;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							  csi[k][j][i].y * csi[k][j][i].y +
							  csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
	    } else 
	      ucor[k][j][i].x=0.;
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux += ucor[k][j][i].y;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/sqrt(eta[k][j][i].x * eta[k][j][i].x +
							  eta[k][j][i].y * eta[k][j][i].y +
							  eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    } else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux += ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux -= ucor[k][j][i].x;
	    if (flg==3)
	    libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							csi[k][j][i].y * csi[k][j][i].y +
							csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	    }else 
	      ucor[k][j][i].x=0.;
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux -= ucor[k][j][i].y;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/ sqrt(eta[k][j][i].x * eta[k][j][i].x +
							   eta[k][j][i].y * eta[k][j][i].y +
							   eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    }else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux -= ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

      }
    }
  }
  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD);
  PetscGlobalSum(&libm_Flux_abs, &ibm_Flux_abs, PETSC_COMM_WORLD);
  PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD);

  PetscReal correction;

  if (*ibm_Area > 1.e-15) {
    if (flg>1) 
      correction = (*ibm_Flux + user->FluxIntpSum)/ ibm_Flux_abs;
    else if (flg)
      correction = (*ibm_Flux + user->FluxIntpSum) / *ibm_Area;
    else
      correction = *ibm_Flux / *ibm_Area;
  }
  else {
    correction = 0;
  }

  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux %le %le %le\n", *ibm_Flux, *ibm_Area, correction);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] <ibmval && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon){
	    if (flg==3) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x + 
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < 1.1) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x +
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

      }
    }
  }
  


  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }
  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD);
  PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux22 %le %le\n", *ibm_Flux, *ibm_Area);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(fda, user->Ucont, &ucor);

  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  return 0;
}

PetscErrorCode MaxPosition(UserCtx *user, PetscInt pos)
{
  
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt i, j, k;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (lidx(i,j,k,user) == pos) {
	  PetscPrintf(PETSC_COMM_SELF, "Position %i %i %i\n", i, j, k);
	}
      }
    }
  }

  return 0;
}
#define Epsilon_Eq 1.e-6
//#define PartFloat_Eq(a, b) (fabs(a-b)>Epislon_Eq) ? 0:1;
#define Float_Eq(a, b) (a==b) ? PETSC_TRUE : (((a-b)<Epsilon_Eq) && (a-b) > -Epsilon_Eq)

PetscErrorCode MyNFaceFine(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscReal ***nvert;
  Cmpnts ***nface;

  PetscInt      i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(fda, user->lNFace, &nface);

  VecSet(user->lNFace, 0.);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (Float_Eq(nvert[k][j][i], 0)) {
	  if (Float_Eq(nvert[k][j][i+1], 1)) {
	    nface[k][j][i].x = 1.;
	  }
	  if (Float_Eq(nvert[k][j+1][i], 1)) {
	    nface[k][j][i].y = 1.;
	  }
	  if (Float_Eq(nvert[k+1][j][i], 1)) {
	    nface[k][j][i].z = 1.;
	  }
	}
	else {
	  nface[k][j][i].x = 1.;
	  nface[k][j][i].y = 1.;
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }
  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lNFace, &nface);

  DALocalToLocalBegin(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  DALocalToLocalEnd(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  return 0;
}

PetscErrorCode MyNFaceRestrict(UserCtx *user)
{
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  
  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  Cmpnts ***nface_h, ***nface;
  PetscReal ***nvert;

  UserCtx *user_h = user->user_f;
  DA	fda_h = user_h->fda, fda = user->fda, da = user->da;
  DAVecGetArray(fda_h, user_h->lNFace, &nface_h);
  DAVecGetArray(fda,   user->lNFace, &nface);

  VecSet(user->lNFace, 0.);

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);
	if (Float_Eq(nface_h[kh   ][jh   ][ih].x, 1) &&
	    Float_Eq(nface_h[kh   ][jh-ja][ih].x, 1) &&
	    Float_Eq(nface_h[kh-ka][jh   ][ih].x, 1) &&
	    Float_Eq(nface_h[kh-ka][jh-ja][ih].x, 1)) {
	  nface[k][j][i].x = 1.;
	}

	if (Float_Eq(nface_h[kh   ][jh][ih   ].y, 1) &&
	    Float_Eq(nface_h[kh   ][jh][ih-ia].y, 1) &&
	    Float_Eq(nface_h[kh-ka][jh][ih   ].y, 1) &&
	    Float_Eq(nface_h[kh-ka][jh][ih-ia].y, 1)) {
	  nface[k][j][i].y = 1.;
	}

	if (Float_Eq(nface_h[kh][jh   ][ih   ].z, 1) &&
	    Float_Eq(nface_h[kh][jh-ja][ih   ].z, 1) &&
	    Float_Eq(nface_h[kh][jh   ][ih-ia].z, 1) &&
	    Float_Eq(nface_h[kh][jh-ja][ih-ia].z, 1)) {
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }

  DAVecGetArray(da, user->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (Float_Eq(nvert[k][j][i], 0)) {
	  if (Float_Eq(nvert[k][j][i+1], 1)) {
	    nface[k][j][i].x = 1.;
	  }
	  if (Float_Eq(nvert[k][j+1][i], 1)) {
	    nface[k][j][i].y = 1.;
	  }
	  if (Float_Eq(nvert[k+1][j][i], 1)) {
	    nface[k][j][i].z = 1.;
	  }
	}
	else {
	  nface[k][j][i].x = 1.;
	  nface[k][j][i].y = 1.;
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }


  
  DAVecRestoreArray(da, user->lNvert, &nvert);

  DAVecRestoreArray(fda, user->lNFace, &nface);
  DAVecRestoreArray(fda_h, user_h->lNFace, &nface_h);

  DALocalToLocalBegin(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  DALocalToLocalEnd(fda, user->lNFace, INSERT_VALUES, user->lNFace);

  //  VecSet(user->lNFace, 0.);
  return 0;
}

PetscErrorCode MyNFaceInit(UserMG *usermg)
{
  PetscInt l;

  PetscInt bi;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user;
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      VecDuplicate(user[bi].lCsi, &user[bi].lNFace);
      if (l == usermg->mglevels-1) {
	MyNFaceFine(&user[bi]);
      }
      else {
	MyNFaceRestrict(&user[bi]);
      }
    }
  }
  return 0;
}

PetscErrorCode MyNFaceFinalize(UserMG *usermg)
{
  PetscInt l;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user;
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    VecDestroy(user->lNFace);
  }
  return 0;
}


PetscErrorCode Do_averaging(UserCtx *user)
{
	DALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
  
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Cmpnts ***ucat, ***u_cross_sum, ***csi, ***eta, ***zet;
	Cmpnts ***u2sum, ***vort_sum, ***vort2_sum, ***uuusum, ***du2sum;
	PetscReal ***p2sum, ***p, ***k_sum, ***lnu_t, ***nut_sum,  ***aj, ***nvert, ***udpsum, ***taussum;
	Cmpnts2 ***ko;
	
	VecAXPY(user->Ucat_sum, 1., user->Ucat);
	VecAXPY(user->P_sum, 1., user->P);
	
	double max_norm;
	VecMax(user->Ucat_sum, &i, &max_norm);
	PetscPrintf(PETSC_COMM_WORLD, "\n*** Max Ucat_avg = %e \n", max_norm / ti);
	
	DAVecGetArray(user->fda, user->lUcat, &ucat);
	DAVecGetArray(user->da, user->lAj, &aj);
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->fda,user->lCsi, &csi);
	DAVecGetArray(user->fda,user->lEta, &eta);
	DAVecGetArray(user->fda,user->lZet, &zet);
	
	DAVecGetArray(user->fda, user->Ucat_square_sum, &u2sum);
	DAVecGetArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
	DAVecGetArray(user->da, user->lP, &p);
	
	if(rans) {
		DAVecGetArray(user->fda2, user->K_Omega, &ko);
		DAVecGetArray(user->da, user->K_sum, &k_sum);
	}
	
	if(les) {
		DAVecGetArray(user->da, user->lNu_t, &lnu_t);
		DAVecGetArray(user->da, user->Nut_sum, &nut_sum);
	}

	if(averaging>=2) {
		DAVecGetArray(user->da, user->P_square_sum, &p2sum);
		//DAVecGetArray(user->fda, user->P_cross_sum, &p_cross_sum);
	}
	if(averaging>=3) {
		if(les) {
			DAVecGetArray(user->da, user->tauS_sum, &taussum);
		}
		
		DAVecGetArray(user->da, user->Udp_sum, &udpsum);
		DAVecGetArray(user->fda, user->dU2_sum, &du2sum);
		DAVecGetArray(user->fda, user->UUU_sum, &uuusum);

		DAVecGetArray(user->fda, user->Vort_sum, &vort_sum);
		DAVecGetArray(user->fda, user->Vort_square_sum, &vort2_sum);
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  double U = ucat[k][j][i].x, V = ucat[k][j][i].y, W = ucat[k][j][i].z;
		u2sum[k][j][i].x += ucat[k][j][i].x * ucat[k][j][i].x;
		u2sum[k][j][i].y += ucat[k][j][i].y * ucat[k][j][i].y;
		u2sum[k][j][i].z += ucat[k][j][i].z * ucat[k][j][i].z;

		if(rans) {
			k_sum[k][j][i] += ko[k][j][i].x;
		}

		if(les) {
			nut_sum[k][j][i] += lnu_t[k][j][i];
		}

		if(averaging>=2) {
			p2sum[k][j][i] += p[k][j][i] * p[k][j][i];
			/*
			p_cross_sum[k][j][i].x += p[k][j][i] * ucat[k][j][i].x;
			p_cross_sum[k][j][i].y += p[k][j][i] * ucat[k][j][i].y;
			p_cross_sum[k][j][i].z += p[k][j][i] * ucat[k][j][i].z;
			*/
		}
		
		if(averaging>=3) {
			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			double dpdc, dpde, dpdz;
			double dp_dx, dp_dy, dp_dz;

			double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
			double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
			double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
			double ajc = aj[k][j][i];
				
			Compute_du_center ( i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_dscalar_center ( i, j, k, mx, my, mz, p, nvert, &dpdc, &dpde, &dpdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);

			double vort_x = dw_dy - dv_dz, vort_y = du_dz - dw_dx, vort_z = dv_dx - du_dy;
			
			vort_sum[k][j][i].x += vort_x;
			vort_sum[k][j][i].y += vort_y;
			vort_sum[k][j][i].z += vort_z;
			
			vort2_sum[k][j][i].x += vort_x*vort_x;
			vort2_sum[k][j][i].y += vort_y*vort_y;
			vort2_sum[k][j][i].z += vort_z*vort_z;

			udpsum[k][j][i] += U * dp_dx + V * dp_dy + W * dp_dz;

			du2sum[k][j][i].x += pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.);
			du2sum[k][j][i].y += pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.);
			du2sum[k][j][i].z += pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.);
			
			uuusum[k][j][i].x += (U*U + V*V + W*W) * U;
			uuusum[k][j][i].y += (U*U + V*V + W*W) * V;
			uuusum[k][j][i].z += (U*U + V*V + W*W) * W;
			
			if(les) {
				double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
				double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
				double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
				double SS = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz;
				taussum[k][j][i] += 2. * lnu_t[k][j][i] * SS;
			}
		}
		
		u_cross_sum[k][j][i].x += ucat[k][j][i].x * ucat[k][j][i].y;	// uv
		u_cross_sum[k][j][i].y += ucat[k][j][i].y * ucat[k][j][i].z;	// vw
		u_cross_sum[k][j][i].z += ucat[k][j][i].z * ucat[k][j][i].x;	// wu	
	}
	
	DAVecRestoreArray(user->fda, user->lUcat, &ucat);
	DAVecRestoreArray(user->da, user->lAj, &aj);
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->fda,user->lCsi, &csi);
	DAVecRestoreArray(user->fda,user->lEta, &eta);
	DAVecRestoreArray(user->fda,user->lZet, &zet);
	
	DAVecRestoreArray(user->fda, user->Ucat_square_sum, &u2sum);
	DAVecRestoreArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
	DAVecRestoreArray(user->da, user->lP, &p);
	
	if(rans) {
          DAVecRestoreArray(user->fda2, user->K_Omega, &ko);
          DAVecRestoreArray(user->da, user->K_sum, &k_sum);
        }
	if(les) {
          DAVecRestoreArray(user->da, user->lNu_t, &lnu_t);
          DAVecRestoreArray(user->da, user->Nut_sum, &nut_sum);
        }

	if(averaging>=2) {
		DAVecRestoreArray(user->da, user->P_square_sum, &p2sum);
		//DAVecRestoreArray(user->fda, user->P_cross_sum, &p_cross_sum);
	}
	if(averaging>=3) {
		if(les) {
			DAVecRestoreArray(user->da, user->tauS_sum, &taussum);
		}
		DAVecRestoreArray(user->da, user->Udp_sum, &udpsum);
		DAVecRestoreArray(user->fda, user->dU2_sum, &du2sum);
		DAVecRestoreArray(user->fda, user->UUU_sum, &uuusum);
		DAVecRestoreArray(user->fda, user->Vort_sum, &vort_sum);
		DAVecRestoreArray(user->fda, user->Vort_square_sum, &vort2_sum);
	}
	
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
			
	if (ti!=0 && ti == (ti/tiout) * tiout && periodic && inletprofile==13 ) {
		Cmpnts ***u_sum, ***uu_sum, ***cent;
		double N=ti;
		double buffer[20][1000];	// [var][points]
		
		for(i=0; i<20; i++)
		for(j=0; j<1000; j++) buffer[i][j]=0;
		
		DAVecGetArray(user->fda, user->lCent, &cent);
		DAVecGetArray(user->fda, user->Ucat_sum, &u_sum);
		DAVecGetArray(user->fda, user->Ucat_square_sum, &uu_sum);
		DAVecGetArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
		
		std::vector<int> count (my), total_count(my);
		std::vector<double> Y_sum(my), U_sum(my), V_sum(my), W_sum(my), UU_sum(my), VV_sum(my), WW_sum(my), UV_sum(my), VW_sum(my), WU_sum(my);
		std::vector<double> Y_sum_tmp(my), U_sum_tmp(my), V_sum_tmp(my), W_sum_tmp(my), UU_sum_tmp(my), VV_sum_tmp(my), WW_sum_tmp(my), UV_sum_tmp(my), VW_sum_tmp(my), WU_sum_tmp(my);
	
		for (j=0; j<my; j++) Y_sum_tmp[j] = 0, U_sum_tmp[j]=0, V_sum_tmp[j]=0, W_sum_tmp[j]=0, UU_sum_tmp[j]=0, VV_sum_tmp[j]=0, WW_sum_tmp[j]=0, UV_sum_tmp[j]=0, VW_sum_tmp[j]=0, WU_sum_tmp[j]=0;
	
		std::fill( count.begin(), count.end(), 0 );
		
			
		for (j=ys; j<ye; j++) {
			
			for (k=lzs; k<lze; k++)
			for (i=lxs; i<lxe; i++) {
				{
					count[j] ++;
					
					double U = u_sum[k][j][i].x / N, V = u_sum[k][j][i].y / N, W = u_sum[k][j][i].z / N;
					double uu = uu_sum[k][j][i].x/N - U*U;
					double vv = uu_sum[k][j][i].y/N - V*V;
					double ww = uu_sum[k][j][i].z/N - W*W;
					double uv = u_cross_sum[k][j][i].x/N - U*V;
					double vw = u_cross_sum[k][j][i].y/N - V*W;
					double wu = u_cross_sum[k][j][i].z/N - W*U;
					
					Y_sum_tmp[j] += cent[k][j][i].y;
					
					U_sum_tmp[j] += U;
					V_sum_tmp[j] += V;
					W_sum_tmp[j] += W;
						
					UU_sum_tmp[j] += uu;
					VV_sum_tmp[j] += vv;
					WW_sum_tmp[j] += ww;
						
					UV_sum_tmp[j] += uv;
					VW_sum_tmp[j] += vw;
					WU_sum_tmp[j] += wu;
				}
			}
		}
		
		DAVecRestoreArray(user->fda, user->lCent, &cent);
		DAVecRestoreArray(user->fda, user->Ucat_sum, &u_sum);
		DAVecRestoreArray(user->fda, user->Ucat_square_sum, &uu_sum);
		DAVecRestoreArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
		
		MPI_Reduce( &count[0], &total_count[0], my, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &Y_sum_tmp[0], &Y_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &U_sum_tmp[0], &U_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &V_sum_tmp[0], &V_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &W_sum_tmp[0], &W_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &UU_sum_tmp[0], &UU_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &VV_sum_tmp[0], &VV_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &WW_sum_tmp[0], &WW_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &UV_sum_tmp[0], &UV_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &VW_sum_tmp[0], &VW_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		MPI_Reduce( &WU_sum_tmp[0], &WU_sum[0], my, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
		
		if(!rank) {
			char filen[80];
			
			sprintf(filen, "%s/Channel_Profile_%06d.dat", path, ti);
			FILE *fp = fopen(filen, "w");
			if(fp) {
				fprintf(fp, "VARIABLES = \"y\" \"y+\" \"w\" \"w+\" \"uu+\" \"vv+\" \"ww+\" \"urms+\" \"vrms+\" \"wrms+\" \"vw+\" \"log\"\n");
				fprintf(fp, "ZONE T=\"Channel\"\n\n");
				for (j=1; j<my-1; j++) {
					
					double Re_tau=user->ustar_now[2]*user->ren;
					double U_star=user->ustar_now[2];
					double Y = Y_sum[j] / total_count[j];
					double Yp = Re_tau*Y;
					double ustar = U_star;
					double ustar2 = ustar * ustar;
					double W = W_sum[j] / total_count[j];
					double uu = UU_sum[j] / total_count[j];
					double vv = VV_sum[j] / total_count[j];
					double ww = WW_sum[j] / total_count[j];
					double vw = VW_sum[j] / total_count[j];
					
					
					fprintf(fp, "%e ", Y);	// y
					fprintf(fp, "%e ", Yp);	// y+
					fprintf(fp, "%e ", W);	// w
					fprintf(fp, "%e ", W/ustar);	// w+
					fprintf(fp, "%e ", uu/ustar2);	// uu+
					fprintf(fp, "%e ", vv/ustar2);	// vv+
					fprintf(fp, "%e ", ww/ustar2);	// ww+
					fprintf(fp, "%e ", sqrt(uu/ustar2));	// urms
					fprintf(fp, "%e ", sqrt(vv/ustar2));	// vrms
					fprintf(fp, "%e ", sqrt(ww/ustar2));	// wrms
					fprintf(fp, "%e ", vw/ustar2);	// vw+
					fprintf(fp, "%e ", 1/.41*log(Yp)+5.2);
					fprintf(fp, "\n");
				}
				fclose(fp);
			}
			
		}
	}
	
	PetscBarrier(PETSC_NULL);
	
	return 0;
}

PetscErrorCode KE_Output(UserCtx *user)
{
	DALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
	PetscInt	lxs, lys, lzs, lxe, lye, lze;
  
	Cmpnts	***ucat;
	PetscReal ***aj;
	PetscScalar ***level;
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	DAVecGetArray(user->fda, user->lUcat, &ucat);
	DAVecGetArray(user->da, user->lAj, &aj);
	if(levelset) {
		DAVecGetArray(user->da, user->lLevelset, &level);
	}
	
	double local_sum=0, sum=0;
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(!air_flow_levelset_periodic){
			local_sum += 0.5 * ucat[k][j][i].x * ucat[k][j][i].x / aj[k][j][i];
			local_sum += 0.5 * ucat[k][j][i].y * ucat[k][j][i].y / aj[k][j][i];
			local_sum += 0.5 * ucat[k][j][i].z * ucat[k][j][i].z / aj[k][j][i];
		}
		if(air_flow_levelset_periodic){
			if(level[k][j][i]<-dthick){
				local_sum += 0.5 * ucat[k][j][i].x * ucat[k][j][i].x / aj[k][j][i];
				local_sum += 0.5 * ucat[k][j][i].y * ucat[k][j][i].y / aj[k][j][i];
				local_sum += 0.5 * ucat[k][j][i].z * ucat[k][j][i].z / aj[k][j][i];	
			}
		}
	}
	PetscGlobalSum(&local_sum, &sum, PETSC_COMM_WORLD);
	
	DAVecRestoreArray(user->fda, user->lUcat, &ucat);
	DAVecRestoreArray(user->da, user->lAj, &aj);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {	
		char filen[80];
		sprintf(filen, "%s/Kinetic_Energy.dat", path);
		FILE *f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d\t%.7e\n", ti, sum);
		fclose(f);
	}
	if(levelset) {
		DAVecRestoreArray(user->da, user->lLevelset, &level);
	}	
	return 0;
}


PetscInt setup_lidx3(UserCtx *user)	// with component, 1,2,3, for momentum
{
	DALocalInfo	info = user->info;
	DA		da = user->da, fda = user->fda;
	PetscInt	gxs, gxe, gys, gye, gzs, gze;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, ***lidm;
	Cmpnts ***gidm;

	Vec	Lid;
	VecDuplicate(user->lUcont, &user->Gidm);
	VecDuplicate(user->lNvert, &Lid);
	
	DAVecGetArray(fda, user->Gidm, &gidm);
	DAVecGetArray(da, Lid, &lidm);
	DAVecGetArray(da, user->lNvert, &nvert);

	gxs = info.gxs; gxe = gxs + info.gxm;
	gys = info.gys; gye = gys + info.gym;
	gzs = info.gzs; gze = gzs + info.gzm;
		
	//int nblank_node[4096], nblank_node_tmp[4096];	// # of blank nodes for processors
	int ndof_node[4096], ndof_node_tmp[4096];	// # of pressure dof for processors
	int ndof_node_accu;
	
	int r, myrank, size;
	
		MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
		MPI_Comm_size(PETSC_COMM_WORLD,&size);
		
		for(r=0; r<size; r++) {
			ndof_node_tmp[r] = 0;
			//nblank_node_tmp[r] = 0;
		}
	
		for(k=zs; k<ze; k++)
		for(j=ys; j<ye; j++)
		for(i=xs; i<xe; i++) {
			if(! (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) &&
					(nvert[k][j][i]+nvert[k][j][i+1]>1.1 || nvert[k][j][i]+nvert[k][j+1][i]>1.1 || nvert[k][j][i]+nvert[k+1][j][i]>1.1) ) {
				//nblank_node_tmp[myrank]+=3;
				lidm[k][j][i]=0;
			}
			else {
				lidm[k][j][i] = (PetscReal)ndof_node_tmp[myrank];
				ndof_node_tmp[myrank] += 3;	// vector size
			}
		}
	
		//MPI_Allreduce( &nblank_node_tmp, &nblank_node, size, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		MPI_Allreduce( &ndof_node_tmp, &ndof_node, size, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		
		ndof_node_accu = 0;
		for(r=0; r<myrank; r++) ndof_node_accu += ndof_node[r];
		
		int n;
		VecGetSize(user->Ucont,&n);
		if(myrank==size-1) {
			printf("\n\n********* momentum: %d %d ***********\n\n", ndof_node_accu + ndof_node[myrank], n);
		}
		
		//MPI_Bcast(&user->reduced_p_size, 1, MPI_INT, size-1, PETSC_COMM_WORLD);
		//PetscPrintf(PETSC_COMM_WORLD, "%d***********************\n", reduced_p_size);
		
		
		PetscBarrier(PETSC_NULL);
		
		for(k=zs; k<ze; k++)
		for(j=ys; j<ye; j++)
		for(i=xs; i<xe; i++) {
			gidm[k][j][i].x = lidm[k][j][i] + ndof_node_accu;	// gidm is double, be careful
			gidm[k][j][i].y = gidm[k][j][i].x + 1;
			gidm[k][j][i].z = gidm[k][j][i].x + 2;
		}
		
		
	VecCreateMPI(PETSC_COMM_WORLD, ndof_node[myrank], PETSC_DETERMINE, &user->Ucont2);
	
		
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, user->Gidm, &gidm);
	DAVecRestoreArray(da, Lid, &lidm);
	
	
	VecDestroy(Lid);
	
	DALocalToLocalBegin(da, user->Gidm, INSERT_VALUES, user->Gidm);
	DALocalToLocalEnd(da, user->Gidm, INSERT_VALUES, user->Gidm);
	return 0;
}


void Convert_Ucont2_Ucont(UserCtx *user)
{
	DALocalInfo	info = user->info;
	DA		da = user->da, fda = user->fda;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, *ucont2;
	Cmpnts ***ucont;
	
	DAVecGetArray(fda, user->Ucont, &ucont);
	DAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(user->Ucont2, &ucont2);
	
	int pos=0;
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if(! (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) &&
					(nvert[k][j][i]+nvert[k][j][i+1]>1.1 || nvert[k][j][i]+nvert[k][j+1][i]>1.1 || nvert[k][j][i]+nvert[k+1][j][i]>1.1) ) {
			// do nothing
		}
		else {
			ucont[k][j][i].x = ucont2[pos++];
			ucont[k][j][i].y = ucont2[pos++];
			ucont[k][j][i].z = ucont2[pos++];
		}
	}
	
	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(user->Ucont2, &ucont2);
	
	DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
}

void Convert_Ucont_Ucont2(UserCtx *user)
{
	DALocalInfo	info = user->info;
	DA		da = user->da, fda = user->fda;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, *ucont2;
	Cmpnts ***ucont;
	
	DAVecGetArray(fda, user->Ucont, &ucont);
	DAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(user->Ucont2, &ucont2);
	
	int pos=0;
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if(! (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) &&
					(nvert[k][j][i]+nvert[k][j][i+1]>1.1 || nvert[k][j][i]+nvert[k][j+1][i]>1.1 || nvert[k][j][i]+nvert[k+1][j][i]>1.1) ) {
			// do nothing
		}
		else {
			ucont2[pos++] = ucont[k][j][i].x;
			ucont2[pos++] = ucont[k][j][i].y;
			ucont2[pos++] = ucont[k][j][i].z;
		}
	}
	
	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(user->Ucont2, &ucont2);
	
	DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
}

void Convert_RHS_RHS2(UserCtx *user, Vec RHS, Vec RHS2)
{
	DALocalInfo	info = user->info;
	DA		da = user->da, fda = user->fda;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, *rhs2;
	Cmpnts ***rhs;
	
	DAVecGetArray(fda, RHS, &rhs);
	DAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(RHS2, &rhs2);
	
	int pos=0;
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if(! (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) &&
					(nvert[k][j][i]+nvert[k][j][i+1]>1.1 || nvert[k][j][i]+nvert[k][j+1][i]>1.1 || nvert[k][j][i]+nvert[k+1][j][i]>1.1) ) {
			// do nothing
		}
		else {
			rhs2[pos++] = rhs[k][j][i].x;
			rhs2[pos++] = rhs[k][j][i].y;
			rhs2[pos++] = rhs[k][j][i].z;
		}
	}
	
	DAVecRestoreArray(fda, RHS, &rhs);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(RHS2, &rhs2);
}
