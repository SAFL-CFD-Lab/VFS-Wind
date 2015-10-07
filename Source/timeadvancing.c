/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

extern PetscReal FluxInSum, FluxOutSum;
extern PetscInt immersed;

PetscErrorCode MyFieldRestriction(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  DA	da_f = *user->da_f;

  UserCtx *user_f = user->user_f;
  DA	fda_f = user_f->fda;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont, ***ucont_f;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DAVecGetArray(fda, user->Ucont, &ucont);
  DAVecGetArray(fda_f, user_f->lUcont, &ucont_f);

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);
	
	ucont[k][j][i].x = (ucont_f[kh   ][jh   ][ih  ].x +
			    ucont_f[kh-ka][jh   ][ih  ].x +
			    ucont_f[kh   ][jh-ja][ih  ].x +
			    ucont_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	ucont[k][j][i].y = (ucont_f[kh   ][jh  ][ih   ].y +
			    ucont_f[kh-ka][jh  ][ih   ].y +
			    ucont_f[kh   ][jh  ][ih-ia].y +
			    ucont_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);
      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	ucont[k][j][i].z = (ucont_f[kh  ][jh   ][ih   ].z +
			    ucont_f[kh  ][jh   ][ih-ia].z +
			    ucont_f[kh  ][jh-ja][ih   ].z +
			    ucont_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);
      }
    }
  }

  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(fda_f, user_f->lUcont, &ucont_f);
  
  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  VecCopy(user->Ucont, user->Ucont_MG);

  PetscReal ***p, ***p_f;

  DAVecGetArray(da, user->P, &p);
  DAVecGetArray(da_f, user_f->lP, &p_f);

  Cmpnts ***ucat, ***ucat_f;
  DAVecGetArray(user_f->fda, user_f->lUcat, &ucat_f);
  DAVecGetArray(user->fda, user->Ucat, &ucat); 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	p[k][j][i] = (p_f[kh   ][jh   ][ih   ] +
		      p_f[kh   ][jh   ][ih-ia] +
		      p_f[kh   ][jh-ja][ih   ] +
		      p_f[kh   ][jh-ja][ih-ia] +
		      p_f[kh-ka][jh   ][ih   ] +
		      p_f[kh-ka][jh   ][ih-ia] +
		      p_f[kh-ka][jh-ja][ih   ] +
		      p_f[kh-ka][jh-ja][ih-ia]) * 0.125;

/* 	ucat[k][j][i].x = (ucat_f[kh   ][jh   ][ih   ].x + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].x + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].x + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].x + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].x + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].x + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].x + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].x) * 0.125; */

/* 	ucat[k][j][i].y = (ucat_f[kh   ][jh   ][ih   ].y + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].y + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].y + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].y + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].y + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].y + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].y + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].y) * 0.125; */

/* 	ucat[k][j][i].z = (ucat_f[kh   ][jh   ][ih   ].z + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].z + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].z + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].z + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].z + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].z + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].z + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].z) * 0.125; */
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user_f->fda, user_f->lUcat, &ucat_f);
  
  DAVecRestoreArray(da, user->P, &p);
  DAVecRestoreArray(da_f, user_f->lP, &p_f);
  
  DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
		      
  return 0;
}

PetscErrorCode GhostNodeVelocity(UserCtx *user)
{

  DA da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec		Coor;
  Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet, ***ucont_o;
  PetscScalar	FluxIn, FluxOut, r, uin, ratio, xc, yc;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
/*   DAVecGetArray(fda, user->Ucont, &ucont); */
  DAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DAVecGetArray(fda, user->Ucat,  &ucat);


  DAVecGetArray(fda, user->lCsi,  &csi);
  DAVecGetArray(fda, user->lEta,  &eta);
  DAVecGetArray(fda, user->lZet,  &zet);

/*   Contra2Cart(user, user->lUcont, user->Ucat); */
/*   if (user->thislevel == user->mglevels-1) { */
/*     if (user->bctype[5]==4) { */
/*       if (ze == mz) { */
/* 	k = ze-2; */
/* 	FluxOut = 0; */
/* 	for (j=lys; j<lye; j++) { */
/* 	  for (i=lxs; i<lxe; i++) { */
/* 	    FluxOut += (ucat[k][j][i].x * (zet[k][j][i].x + zet[k-1][j][i].x) + */
/* 			ucat[k][j][i].y * (zet[k][j][i].y + zet[k-1][j][i].y) + */
/* 			ucat[k][j][i].z * (zet[k][j][i].z + zet[k-1][j][i].z)) * 0.5; */
/* 	  } */
/* 	} */
/*       } */
/*       else { */
/* 	FluxOut = 0.; */
/*       } */
  
/*       PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD); */
/*       ratio = FluxInSum / FluxOutSum; */
/*       if (FluxOutSum < 1.e-6) ratio = 1.; */
/*       //    PetscPrintf(PETSC_COMM_SELF, "Ratio %d %le %le %le %d %d\n", ti, ratio, FluxInSum, FluxOutSum, zs, ze); */

/*       if (ze==mz && fabs(FluxOutSum) > 1.e-6) { */
/* 	k = ze-1; */
/* 	for (j=lys; j<lye; j++) { */
/* 	  for (i=lxs; i<lxe; i++) {   */
/* 	    ubcs[k][j][i].x = ucat[k-1][j][i].x; */
/* 	    ubcs[k][j][i].y = ucat[k-1][j][i].y; */
/* 	    ubcs[k][j][i].z = ucat[k-1][j][i].z * ratio; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*     else if (user->bctype[5]==0) { */
/*       if (ze==mz) { */
/* 	k = ze-1; */
/* 	for (j=lys; j<lye; j++) { */
/* 	  for (i=lxs; i<lxe; i++) {   */
/* 	    ubcs[k][j][i].x = ucat[k-1][j][i].x; */
/* 	    ubcs[k][j][i].y = ucat[k-1][j][i].y; */
/* 	    ubcs[k][j][i].z = ucat[k-1][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
  if (user->bctype[5]==2) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }

/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }


  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }

/* ==================================================================================             */
/*   INTERFACE BC */
/* ==================================================================================             */
  if (user->bctype[0]==0) {
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i+1].x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==0) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i-1].y;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }

  if (user->bctype[2]==0) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = ucat[k][j+1][i].y;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==0) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = ucat[k][j-1][i].y;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }


  if (user->bctype[4]==0) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==0) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	}
      }
    }
  }

/*   if (user->bctype[4]==5) { */
/*     if (zs==0) { */
/*       k = 0; */
/*       for (j=lys; j<lye; j++) { */
/* 	for (i=lxs; i<lxe; i++) {   */
/* 	  ubcs[k][j][i].x = 0.; */
/* 	  ubcs[k][j][i].y = 0.; */
/* 	  ubcs[k][j][i].z = 1.; */
/* 	} */
/*       } */
/*     }     */
/*   } */
  // boundary conditions on ghost nodes
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;

/* 	// w = 0 for cavity */
/* 	ucat[k][j][i].x = ucat[k][j+1][i].x; */
/* 	ucat[k][j][i].y = ucat[k][j+1][i].y; */
/* 	ucat[k][j][i].z = -ucat[k][j+1][i].z; */
      }
    }
  }

  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
/* 	// w = 0 for cavity */
/* 	ucat[k][j][i].x = ucat[k][j-1][i].x; */
/* 	ucat[k][j][i].y = ucat[k][j-1][i].y; */
/* 	ucat[k][j][i].z = -ucat[k][j-1][i].z; */
      }
    }
  }

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }


  // 0 velocity on the corner point
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }

/*   DAVecRestoreArray(fda, user->Ucont, &ucont); */
  DAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DAVecRestoreArray(fda, user->Ucat,  &ucat);
  DAVecRestoreArray(fda, Coor, &coor);

  DAVecRestoreArray(fda, user->lCsi,  &csi);
  DAVecRestoreArray(fda, user->lEta,  &eta);
  DAVecRestoreArray(fda, user->lZet,  &zet);

  //  DAVecRestoreArray(fda, user->Ucont_o, &ucont_o);

  DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  return(0);
  
}

PetscErrorCode MyRKRHSRestriction(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  UserCtx *user_f = user->user_f;
  DA	fda_f = user_f->fda;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***f, ***rhs_f;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  VecSet(user->Forcing, 0.);
  DAVecGetArray(fda, user->Forcing, &f);
  Vec lRhs;
  VecDuplicate(user_f->lUcont, &lRhs);
  DAGlobalToLocalBegin(fda_f, user_f->Rhs, INSERT_VALUES, lRhs);
  DAGlobalToLocalEnd(fda_f, user_f->Rhs, INSERT_VALUES, lRhs);
  DAVecGetArray(fda_f, lRhs, &rhs_f);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (i==0 || i==mx-2) {
	  f[k][j][i].x = 0.;
	}
	else {
	  f[k][j][i].x = ((rhs_f[kh   ][jh   ][ih   ].x +
			   rhs_f[kh-ka][jh   ][ih   ].x +
			   rhs_f[kh   ][jh-ja][ih   ].x +
			   rhs_f[kh-ka][jh-ja][ih   ].x) * (1.+ka) * (1.+ja) +
			  (rhs_f[kh   ][jh   ][ih-ia].x +
			   rhs_f[kh-ka][jh   ][ih-ia].x +
			   rhs_f[kh   ][jh-ja][ih-ia].x +
			   rhs_f[kh-ka][jh-ja][ih-ia].x) * (1.+ka) * (1.+ja) / 2.+
			  (rhs_f[kh   ][jh   ][ih+ia].x +
			   rhs_f[kh-ka][jh   ][ih+ia].x +
			   rhs_f[kh   ][jh-ja][ih+ia].x +
			   rhs_f[kh-ka][jh-ja][ih+ia].x) * (1.+ka) * (1.+ja) / 2.)
	    * 0.125;
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (j==0 || j==my-2) {
	  f[k][j][i].y = 0.;
	}
	else {
	  f[k][j][i].y = ((rhs_f[kh   ][jh   ][ih   ].y +
			   rhs_f[kh-ka][jh   ][ih   ].y +
			   rhs_f[kh   ][jh   ][ih-ia].y +
			   rhs_f[kh-ka][jh   ][ih-ia].y) * (1.+ka) * (1.+ia) +
			  (rhs_f[kh   ][jh-ja][ih   ].y +
			   rhs_f[kh-ka][jh-ja][ih   ].y +
			   rhs_f[kh   ][jh-ja][ih-ia].y +
			   rhs_f[kh-ka][jh-ja][ih-ia].y) * (1.+ka) * (1.+ia) / 2 +
			  (rhs_f[kh   ][jh+ja][ih   ].y +
			   rhs_f[kh-ka][jh+ja][ih   ].y +
			   rhs_f[kh   ][jh+ja][ih-ia].y +
			   rhs_f[kh-ka][jh+ja][ih-ia].y) * (1.+ka) * (1.+ia) / 2)
	    * 0.125;
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (k==0 || k==mz-2) {
	  f[k][j][i].z = 0.;
	}
	else {
	  f[k][j][i].z = ((rhs_f[kh   ][jh   ][ih   ].z +
			   rhs_f[kh   ][jh   ][ih-ia].z +
			   rhs_f[kh   ][jh-ja][ih   ].z +
			   rhs_f[kh   ][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) +
			  (rhs_f[kh-ka][jh   ][ih   ].z +
			   rhs_f[kh-ka][jh   ][ih-ia].z +
			   rhs_f[kh-ka][jh-ja][ih   ].z +
			   rhs_f[kh-ka][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) / 2 +
			  (rhs_f[kh+ka][jh   ][ih   ].z +
			   rhs_f[kh+ka][jh   ][ih-ia].z +
			   rhs_f[kh+ka][jh-ja][ih   ].z +
			   rhs_f[kh+ka][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) / 2)
	    * 0.125;
	}
      }
    }
  }

  DAVecRestoreArray(fda, user->Forcing, &f);
/*   VecScale(user->Forcing, -1); */
  DAVecRestoreArray(fda_f, lRhs, &rhs_f);
  VecDestroy(lRhs);
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


PetscErrorCode MyRKRHSInterpolation(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  DA	da_c = *user->da_c;

  UserCtx *user_c = user->user_c;

  DA	fda_c = user_c->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
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

  Vec dU;
  Cmpnts	***du, ***f;

  PetscInt	i, j, k, ic, jc, kc, ia, ja, ka;
  PetscReal	***nvert_c;

  PetscInt ii, jj, kk;
  if (*(user->isc)) ii = 0;
  else ii = 1;

  if (*(user->jsc)) jj = 0;
  else jj = 1;

  if (*(user->ksc)) kk = 0;
  else kk = 1;

  VecDuplicate(user->Ucont, &dU);
  VecSet(dU, 0.);

  Vec ldU;
  DAGetLocalVector(fda, &ldU);
  VecSet(ldU, 0.);
  DAVecGetArray(fda, ldU, &du);


  VecWAXPY(user_c->Forcing, -1., user_c->Ucont_MG, user_c->Ucont);

  Vec lForcing;
  DACreateLocalVector(fda_c, &lForcing);
  DAGlobalToLocalBegin(fda_c, user_c->Forcing, INSERT_VALUES, lForcing);
  DAGlobalToLocalEnd(fda_c, user_c->Forcing, INSERT_VALUES, lForcing);

  DAVecGetArray(fda_c, lForcing, &f);

  DAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);

/*   totally wrong! */
/*     1) du is flux not velocity! therefore surface or aj should be */
/*        come into calcutations. */

  /* i component */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	if (!(i%2)) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].x = (f[kc   ][jc   ][ic   ].x * 9 +
			   f[kc   ][jc+ja][ic   ].x * 3 +
			   f[kc+ka][jc   ][ic   ].x * 3 +
			   f[kc+ka][jc+ja][ic   ].x) / 16. / ((1.+kk)*(1.+jj));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {      
      for (i=lxs; i<lxe; i++) {

	if ((i%2)) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].x = (du[k][j][i-1].x + du[k][j][i+1].x) * 0.5;
	}

      }
    }
  }

  /* j component */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      if (!(j%2)) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].y = (f[kc   ][jc   ][ic   ].y * 9 +
			   f[kc   ][jc   ][ic+ia].y * 3 +
			   f[kc+ka][jc   ][ic   ].y * 3 +
			   f[kc+ka][jc   ][ic+ia].y) / 16. / ((1.+kk)*(1.+ii));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      if (j%2) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].y = (du[k][j-1][i].y + du[k][j+1][i].y) * 0.5;

	}
      }
    }
  }

  /* k component */
  for (k=lzs; k<lze; k++) {
    if (!(k%2)) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].z = (f[kc   ][jc   ][ic   ].z * 9 +
			   f[kc   ][jc   ][ic+ia].z * 3 +
			   f[kc   ][jc+ja][ic   ].z * 3 +
			   f[kc   ][jc+ja][ic+ia].z) / 16. / ((1.+jj)*(1.+ii));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);

  for (k=lzs; k<lze; k++) {
    if ((k%2)) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].z = (du[k-1][j][i].z + du[k+1][j][i].z) * 0.5;
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DAVecRestoreArray(fda_c, lForcing, &f);

  DALocalToGlobal(fda, ldU, INSERT_VALUES, dU);


  VecAXPY(user->Ucont, 1., dU);

  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DARestoreLocalVector(fda, &ldU);
  VecDestroy(dU);
  VecDestroy(lForcing);

  DAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);
  
  return 0;
}

#define MG_MAXIT 12
/*
PetscErrorCode TSMG(UserCtx *user)
{
  PetscInt l;

  PetscInt mglevels = user->mglevels;

  PetscInt mgit, ps, *PseudoSteps;

  PetscMalloc(mglevels*sizeof(PetscInt), &PseudoSteps);
  for (l=0; l<mglevels; l++) {
    PseudoSteps[l] = (l||mglevels==1) ? 1:10;
  }

  UserCtx *user_mg;
  user_mg = user;
  for (l=mglevels-1; l>=0; l--) {
    if (user->thislevel == user->mglevels-1) {
      InflowFlux(user);
      OutflowFlux(user);
      
      FormBCS(user);
      if (immersed) {
	ibm_interpolation_advanced(user, user->ibm);
      }
    }

    VecDuplicate(user_mg->Ucont, &user_mg->Ucont_o);
    VecCopy(user_mg->Ucont, user_mg->Ucont_o);
    if (l) { // Not on the coarest grid level
      MyFieldRestriction(user_mg->user_c);
      user_mg = user_mg->user_c;
    }
  }

  user_mg = user;
  PetscTruth MG_Converged = PETSC_FALSE;
  PetscReal norm, norm_old, temp;

  norm_old=0.;
  for (mgit = 0; mgit < MG_MAXIT; mgit++) {
    for (l=mglevels-1; l>=0; l--) {
      for (ps=0; ps<PseudoSteps[l]; ps++) {
	RungeKutta_MG(user_mg, ps);
      }

      if (l==mglevels-1) {
	VecWAXPY(user_mg->Rhs, -1., user_mg->Ucont_o, user_mg->Ucont);
	VecNorm(user_mg->Rhs, NORM_2, &norm);
	temp = norm;
	norm = norm-norm_old;
	norm_old = temp;
	PetscPrintf(PETSC_COMM_WORLD, "TS MG Norm %e\n", norm);
	if (norm < 1.e-20 && norm > -1.e-20) {
	  PetscPrintf(PETSC_COMM_WORLD, "TSConverged %i\n", mgit);
	  MG_Converged = PETSC_TRUE;
	  break;
	}
      }

      if (l>0) {
	//	user_mg->istage = 0;
	ComputeRHS(user_mg, 1);
	MyFieldRestriction(user_mg->user_c);
	MyRKRHSRestriction(user_mg->user_c);
	user_mg = user_mg->user_c;
      }
      if (MG_Converged) break;
    }

    if (MG_Converged) break;
    if (mglevels>1) user_mg = user_mg->user_f;
    for (l=1; l<mglevels; l++) {
      MyRKRHSInterpolation(user_mg);
      if (l<mglevels-1) user_mg = user_mg->user_f;
    }
  }
  PetscFree(PseudoSteps);

  user_mg = user;
  for (l=mglevels-1; l>=0; l--) {
    VecDestroy(user_mg->Ucont_o);
    user_mg = user_mg->user_c;
  }

  return 0;
}*/
