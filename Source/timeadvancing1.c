/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
#include "petscsnes.h"

PetscErrorCode FormFunctionSNES(SNES snes, Vec X, Vec Rhs, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  Cmpnts ***rhs;

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;


  Vec		ICsi = user->lICsi, IEta = user->lIEta, IZet = user->lIZet;
  Vec		JCsi = user->lJCsi, JEta = user->lJEta, JZet = user->lJZet;
  Vec		KCsi = user->lKCsi, KEta = user->lKEta, KZet = user->lKZet;
  Vec		IAj = user->lIAj, JAj = user->lJAj, KAj = user->lKAj;

  Vec dUcont;

  Vec	Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  PetscReal	***p, ***iaj, ***jaj, ***kaj, ***aj;

  PetscInt i, j, k; 
  PetscReal ***nvert, ***nvert_o;

  PetscReal dpdc, dpde, dpdz;

  VecDuplicate(Rhs, &dUcont);
  
  DAVecGetArray(user->fda, Rhs, &rhs);

  Vec lUcont;
  DAGetLocalVector(user->fda, &lUcont);  


  DAGlobalToLocalBegin(user->fda, X, INSERT_VALUES, lUcont);
  DAGlobalToLocalEnd(user->fda, X, INSERT_VALUES, lUcont);

  PetscBarrier(PETSC_NULL);
  Contra2Cart(user);

/*   OutflowVelocity(user, X); */

/*   DAGlobalToLocalBegin(user->fda, X, INSERT_VALUES, lUcont); */

/*   DAGlobalToLocalEnd(user->fda, X, INSERT_VALUES, lUcont); */
  GhostNodeVelocity(user);

  Vec Conv, Visc, Rc, Rct;
  VecDuplicate(user->Ucont, &Conv);
  VecDuplicate(user->Ucont, &Visc);
  VecDuplicate(user->Ucont, &Rc);
  VecDuplicate(user->lUcont, &Rct);

  Convection(user, lUcont, user->lUcat, Conv);
  Viscous(user, lUcont, user->lUcat, Visc);

  VecWAXPY(Rc, -1., Conv, Visc);

  VecDestroy(Conv);
  VecDestroy(Visc);

  Cmpnts ***rc, ***rct;
  DAVecGetArray(user->fda, Rc, &rc);

  DAVecGetArray(user->fda, Rct, &rct);
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  DAVecGetArray(fda, Csi, &csi);
  DAVecGetArray(fda, Eta, &eta);
  DAVecGetArray(fda, Zet, &zet);

  DAVecGetArray(da, user->lAj, &aj);

  /* Calculate the contravariant rhs from the cartesian rhs */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	rct[k][j][i].x = aj[k][j][i] * 
	  (0.5 * (csi[k][j][i].x + csi[k][j][i-1].x) * rc[k][j][i].x +
	   0.5 * (csi[k][j][i].y + csi[k][j][i-1].y) * rc[k][j][i].y +
	   0.5 * (csi[k][j][i].z + csi[k][j][i-1].z) * rc[k][j][i].z);
	rct[k][j][i].y = aj[k][j][i] *
	  (0.5 * (eta[k][j][i].x + eta[k][j-1][i].x) * rc[k][j][i].x +
	   0.5 * (eta[k][j][i].y + eta[k][j-1][i].y) * rc[k][j][i].y +
	   0.5 * (eta[k][j][i].z + eta[k][j-1][i].z) * rc[k][j][i].z);
	rct[k][j][i].z = aj[k][j][i] * 
	  (0.5 * (zet[k][j][i].x + zet[k-1][j][i].x) * rc[k][j][i].x +
	   0.5 * (zet[k][j][i].y + zet[k-1][j][i].y) * rc[k][j][i].y +
	   0.5 * (zet[k][j][i].z + zet[k-1][j][i].z) * rc[k][j][i].z);
      }
    }
  }

  DAVecRestoreArray(da, user->lAj, &aj);
  DAVecRestoreArray(fda, Rct, &rct);
  DAVecRestoreArray(fda, Rc, &rc);

  VecDestroy(Rc);

  DALocalToLocalBegin(fda, Rct, INSERT_VALUES, Rct);
  DALocalToLocalEnd(fda, Rct, INSERT_VALUES, Rct);

  DAVecGetArray(fda, Rct, &rct);

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(da, user->lNvert_o, &nvert_o);

  DAVecGetArray(fda, ICsi, &icsi);
  DAVecGetArray(fda, IEta, &ieta);
  DAVecGetArray(fda, IZet, &izet);
           
  DAVecGetArray(fda, JCsi, &jcsi);
  DAVecGetArray(fda, JEta, &jeta);
  DAVecGetArray(fda, JZet, &jzet);
           
  DAVecGetArray(fda, KCsi, &kcsi);
  DAVecGetArray(fda, KEta, &keta);
  DAVecGetArray(fda, KZet, &kzet);

  DAVecGetArray(da, IAj, &iaj);
  DAVecGetArray(da, JAj, &jaj);
  DAVecGetArray(da, KAj, &kaj);

  DAVecGetArray(da, user->lP, &p);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	dpdc = p[k][j][i+1] - p[k][j][i];

	if ((int)(nvert[k][j+1][i]+0.5)==1 ||
	    (int)(nvert[k][j+1][i+1]+0.5)==1 || j==my-2) {
	  dpde = (p[k][j  ][i  ] - p[k][j-1][i  ] +
		  p[k][j  ][i+1] - p[k][j-1][i+1]) * 0.5;
	}
	else if ((int)(nvert[k][j-1][i]+0.5)==1 ||
		 (int)(nvert[k][j-1][i+1]+0.5)==1 || j==1) {
	  dpde = (p[k][j+1][i  ] - p[k][j  ][i  ] +
		  p[k][j+1][i+1] - p[k][j  ][i+1]) * 0.5;
	}
	else {
	  dpde = (p[k][j+1][i] - p[k][j-1][i] +
		  p[k][j+1][i+1] - p[k][j-1][i+1]) * 0.25;
	}

	if ((int)(nvert[k+1][j][i]+0.5)==1 ||
	    (int)(nvert[k+1][j][i+1]+0.5)==1 || k==mz-2) {
	  dpdz = (p[k][j][i  ] - p[k-1][j][i  ] +
		  p[k][j][i+1] - p[k-1][j][i+1]) * 0.5;
	}
	else if ((int)(nvert[k-1][j][i]+0.5)==1 ||
		 (int)(nvert[k-1][j][i+1]+0.5)==1 || k==1) {
	  dpdz = (p[k+1][j][i  ] - p[k][j][i  ] +
		  p[k+1][j][i+1] - p[k][j][i+1]) * 0.5;
	}
	else {
	  dpdz = (p[k+1][j][i] - p[k-1][j][i] +
		  p[k+1][j][i+1] - p[k-1][j][i+1]) * 0.25;
	}
	    
/* 	if (i>1) { */
/* 	  rhs[k][j][i].x = 0.125 * (rct[k][j][i+1].x * 3. + */
/* 				    rct[k][j][i  ].x * 6. - */
/* 				    rct[k][j][i-1].x); */
/* 	} */
/* 	else { */
/* 	rhs[k][j][i].x = (rct[k][j][i].x * gs[k][j][i+1].x + */
/* 			  rct[k][j][i+1].x * gs[k][j][i].x) / */
/* 	  (gs[k][j][i+1].x + gs[k][j][i].x); */
	rhs[k][j][i].x = 0.5 * (rct[k][j][i].x + rct[k][j][i+1].x);
	
	
/*  	if (ucont[k][j][i].x > 0) { */
/* 	  if (i>1 && (int)(nvert[k][j][i-1]+0.5)==0) { */
/* /\* 	    c1 = 0.5 * gs[k][j][i].x; *\/ */
/* /\* 	    c2 = 0.5 * gs[k][j][i].x + gs[k][j][i-1].x; *\/ */
/* /\* 	    c3 = gs[k][j][i].x; *\/ */
/* /\* 	    c4 = gs[k][j][i].x + gs[k][j][i-1].x; *\/ */

/* /\* 	    g1 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	    c2 = gs[k][j][i].x * 0.5; *\/ */
/* /\* 	    c3 = gs[k][j][i-1].x; *\/ */

/* /\* 	    g2 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	    rhs[k][j][i].x = (rct[k][j][i-1].x + *\/ */
/* /\* 		    (rct[k][j][i].x - rct[k][j][i-1].x) * g1 + *\/ */
/* /\* 		    (rct[k][j][i-1].x - rct[k][j][i-2].x) * g2); *\/ */
/* 	    rhs[k][j][i].x = coef * (-    rct[k][j][i-1].x - */
/* 				     2. * rct[k][j][i  ].x + */
/* 				     3. * rct[k][j][i+1].x) + */
/* 	      rct[k][j][i].x; */
/* 	  } */
/* 	  else { */
/* 	    rhs[k][j][i].x = 0.5 * (rct[k][j][i+1].x + rct[k][j][i].x); */
/* 	  } */
/* 	} */
/* 	else { */
/* 	  if (i < mx-3&&(int)(nvert[k][j][i+2]+0.5)==0) { */
/* /\* 	    c1 = -0.5 * gs[k][j][i].x; *\/ */
/* /\* 	    c2 = -(0.5 * gs[k][j][i].x + gs[k][j][i+1].x); *\/ */
/* /\* 	    c3 = -gs[k][j][i+1].x; *\/ */
/* /\* 	    c4 = -(gs[k][j][i].x + gs[k][j][i+1].x); *\/ */

/* /\* 	    g3 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	    c2 = -(0.5 * gs[k][j][i].x); *\/ */
/* /\* 	    c3 = -gs[k][j][i+1].x; *\/ */

/* /\* 	    g4 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	    rhs[k][j][i].x = (rct[k][j][i].x + *\/ */
/* /\* 		    (rct[k][j][i-1].x - rct[k][j][i].x) * g3 + *\/ */
/* /\* 		    (rct[k][j][i].x   - rct[k][j][i+1].x) * g4); *\/ */
/* 	    rhs[k][j][i].x = coef * (-    rct[k][j][i+2].x - */
/* 				     2. * rct[k][j][i+1].x + */
/* 				     3. * rct[k][j][i  ].x) + */
/* 	      rct[k][j][i+1].x; */
/* 	  } */
/* 	  else { */
/* 	    rhs[k][j][i].x = 0.5 * (rct[k][j][i+1].x + rct[k][j][i].x); */
/* 	  } */
/* 	} */
	  
/* 	} */
/* 	if ((int)(nvert_o[k][j][i  ]+0.5)==0 && */
/* 	    (int)(nvert_o[k][j][i+1]+0.5)==0) { */
	  rhs[k][j][i].x -=
	    (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		     icsi[k][j][i].y * icsi[k][j][i].y +
		     icsi[k][j][i].z * icsi[k][j][i].z) +
	     dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		     ieta[k][j][i].y * icsi[k][j][i].y +
		     ieta[k][j][i].z * icsi[k][j][i].z) +
	     dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		     izet[k][j][i].y * icsi[k][j][i].y +
		     izet[k][j][i].z * icsi[k][j][i].z)) * iaj[k][j][i];
/* 	} */
/* 	PetscPrintf(PETSC_COMM_WORLD, "x%le %le %d %d %d\n", icsi[k][j][i].x, csi[k][j][i].x, i, j, k); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "y%le %le\n", icsi[k][j][i].y, csi[k][j][i].y); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "z%le %le\n", icsi[k][j][i].z, csi[k][j][i].z); */

	if ((int)(nvert[k][j][i+1]+0.5)==1 ||
	    (int)(nvert[k][j+1][i+1]+0.5)==1 || i==mx-2) {
	  dpdc = (p[k][j  ][i] - p[k][j  ][i-1] +
		  p[k][j+1][i] - p[k][j+1][i-1]) * 0.5;
	}
	else if ((int)(nvert[k][j][i-1]+0.5)==1 ||
		 (int)(nvert[k][j+1][i-1]+0.5)==1 || i==1) {
	  dpdc = (p[k][j  ][i+1] - p[k][j  ][i] +
		  p[k][j+1][i+1] - p[k][j+1][i]) * 0.5;
	}
	else {
	  dpdc = (p[k][j  ][i+1] - p[k][j  ][i-1] +
		  p[k][j+1][i+1] - p[k][j+1][i-1]) * 0.25;
	}

	dpde = p[k][j+1][i] - p[k][j][i];
	
	if ((int)(nvert[k+1][j][i]+0.5)==1  ||
	    (int)(nvert[k+1][j+1][i]+0.5)==1 || k==mz-2) {
	  dpdz = (p[k][j  ][i] - p[k-1][j  ][i] +
		  p[k][j+1][i] - p[k-1][j+1][i]) * 0.5;
	}
	else if ((int)(nvert[k-1][j][i]+0.5)==1  ||
		 (int)(nvert[k-1][j+1][i]+0.5)==1 || k==1) {
	  dpdz = (p[k+1][j  ][i] - p[k][j  ][i] +
		  p[k+1][j+1][i] - p[k][j+1][i]) * 0.5;
	}
	else {
	  dpdz = (p[k+1][j  ][i] - p[k-1][j  ][i] +
		  p[k+1][j+1][i] - p[k-1][j+1][i]) * 0.25;
	}

/* 	if (j>1) { */
/* 	  rhs[k][j][i].y = 0.125 * (rct[k][j+1][i].y * 3. + */
/* 				    rct[k][j][i].y * 6. - */
/* 				    rct[k][j-1][i].y); */
/* 	} */
/* 	else { */
/* 	  rhs[k][j][i].y = (rct[k][j][i].y * gs[k][j+1][i].y + */
/* 			    rct[k][j+1][i].y * gs[k][j][i].y) / */
/* 	    (gs[k][j][i].y + gs[k][j+1][i].y); */
	  rhs[k][j][i].y = 0.5 * (rct[k][j][i].y + rct[k][j+1][i].y);

/* 	  if (ucont[k][j][i].y > 0) { */
/* 	    if (j>1&&(int)(nvert[k][j-1][i]+0.5)==0) { */
/* /\* 	      c1 = 0.5 * gs[k][j][i].y; *\/ */
/* /\* 	      c2 = 0.5 * gs[k][j][i].y + gs[k][j-1][i].y; *\/ */
/* /\* 	      c3 = gs[k][j][i].y; *\/ */
/* /\* 	      c4 = gs[k][j][i].y + gs[k][j-1][i].y; *\/ */

/* /\* 	      g1 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      c2 = gs[k][j][i].y * 0.5; *\/ */
/* /\* 	      c3 = gs[k][j-1][i].y; *\/ */

/* /\* 	      g2 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      rhs[k][j][i].y = (rct[k][j-1][i].y + *\/ */
/* /\* 		      (rct[k][j][i].y - rct[k][j-1][i].y) * g1 + *\/ */
/* /\* 		      (rct[k][j-1][i].y - rct[k][j-2][i].y) * g2); *\/ */
/* 	    rhs[k][j][i].y = coef * (-    rct[k][j-1][i].y - */
/* 				     2. * rct[k][j  ][i].y + */
/* 				     3. * rct[k][j+1][i].y) + */
/* 	      rct[k][j][i].y; */

/* 	    } */
/* 	    else { */
/* 	      rhs[k][j][i].y = 0.5 * (rct[k][j+1][i].y + rct[k][j][i].y); */
/* 	    } */
/* 	  } */
/* 	  else { */
/* 	    if (j < my-3&&(int)(nvert[k][j+2][i]+0.5)==0) { */
/* /\* 	      c1 = -0.5 * gs[k][j][i].y; *\/ */
/* /\* 	      c2 = -(0.5 * gs[k][j][i].y + gs[k][j+1][i].y); *\/ */
/* /\* 	      c3 = -gs[k][j+1][i].y; *\/ */
/* /\* 	      c4 = -(gs[k][j][i].y + gs[k][j+1][i].y); *\/ */

/* /\* 	      g3 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      c2 = -(0.5 * gs[k][j][i].y); *\/ */
/* /\* 	      c3 = -gs[k][j+1][i].y; *\/ */

/* /\* 	      g4 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      rhs[k][j][i].y = (rct[k][j][i].y + *\/ */
/* /\* 		      (rct[k][j-1][i].y - rct[k][j][i].y) * g3 + *\/ */
/* /\* 		      (rct[k][j][i].y   - rct[k][j+1][i].y) * g4); *\/ */
/* 	      rhs[k][j][i].y = coef * (-    rct[k][j+2][i].y - */
/* 				       2. * rct[k][j+1][i].y + */
/* 				       3. * rct[k][j  ][i].y) + */
/* 		rct[k][j+1][i].y; */

/* 	    } */
/* 	    else { */
/* 	      rhs[k][j][i].y = 0.5 * (rct[k][j+1][i].y + rct[k][j][i].y); */
/* 	    } */
/* 	  } */


/* 	} */
/* 	if ((int)(nvert_o[k][j  ][i]+0.5)==0 && */
/* 	    (int)(nvert_o[k][j+1][i]+0.5)==0) { */
	  rhs[k][j][i].y -=
	    (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		     jcsi[k][j][i].y * jeta[k][j][i].y +
		     jcsi[k][j][i].z * jeta[k][j][i].z) +
	     dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		     jeta[k][j][i].y * jeta[k][j][i].y +
		     jeta[k][j][i].z * jeta[k][j][i].z) +
	     dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		     jzet[k][j][i].y * jeta[k][j][i].y +
		     jzet[k][j][i].z * jeta[k][j][i].z)) * jaj[k][j][i];
/*  	} */
	if ((int)(nvert[k][j][i+1]+0.5)==1 ||
	    (int)(nvert[k+1][j][i+1]+0.5)==1 || i==mx-2) {
	  dpdc = (p[k  ][j][i] - p[k  ][j][i-1] +
		  p[k+1][j][i] - p[k+1][j][i-1]) * 0.5;
	}
	else if ((int)(nvert[k][j][i-1]+0.5)==1 ||
		 (int)(nvert[k+1][j][i-1]+0.5)==1 || i==1) {
	  dpdc = (p[k  ][j][i+1] - p[k  ][j][i] +
		  p[k+1][j][i+1] - p[k+1][j][i]) * 0.5;
	}
	else {
	  dpdc = (p[k  ][j][i+1] - p[k  ][j][i-1] +
		  p[k+1][j][i+1] - p[k+1][j][i-1]) * 0.25;
	}

	if ((int)(nvert[k][j+1][i]+0.5) ==1 ||
	    (int)(nvert[k+1][j+1][i]+0.5)==1 || j==my-2) {
	  dpde = (p[k  ][j][i] - p[k  ][j-1][i] +
		  p[k+1][j][i] - p[k+1][j-1][i]) * 0.5;
	}
	else if ((int)(nvert[k][j-1][i]+0.5) ==1 ||
		 (int)(nvert[k+1][j-1][i]+0.5)==1 || j==1) {
	  dpde = (p[k  ][j+1][i] - p[k  ][j][i] +
		  p[k+1][j+1][i] - p[k+1][j][i]) * 0.5;
	}
	else {
	  dpde = (p[k  ][j+1][i] - p[k  ][j-1][i] +
		  p[k+1][j+1][i] - p[k+1][j-1][i]) * 0.25;
	}
	dpdz = (p[k+1][j][i] - p[k][j][i]);


/* 	if (k>1) { */
/* 	  rhs[k][j][i].z = 0.125 * (rct[k+1][j][i].z * 3. + */
/* 				    rct[k][j][i].z * 6. - */
/* 				    rct[k-1][j][i].x); */
/* 	} */
/* 	else { */
/* 	  rhs[k][j][i].z = (rct[k][j][i].z * gs[k+1][j][i].z + */
/* 			    rct[k+1][j][i].z * gs[k][j][i].z) / */
/* 	    (gs[k+1][j][i].z + gs[k][j][i].z); */

	  rhs[k][j][i].z = 0.5 * (rct[k][j][i].z + rct[k+1][j][i].z);

/* 	  if (ucont[k][j][i].z > 0) { */
/* 	    if (k>1&&(int)(nvert[k-1][j][i]+0.5)==0) { */
/* /\* 	      c1 = 0.5 * gs[k][j][i].z; *\/ */
/* /\* 	      c2 = 0.5 * gs[k][j][i].z + gs[k-1][j][i].z; *\/ */
/* /\* 	      c3 = gs[k][j][i].z; *\/ */
/* /\* 	      c4 = gs[k][j][i].z + gs[k-1][j][i].z; *\/ */

/* /\* 	      g1 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      c2 = gs[k][j][i].z * 0.5; *\/ */
/* /\* 	      c3 = gs[k-1][j][i].z; *\/ */

/* /\* 	      g2 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      rhs[k][j][i].z = (rct[k-1][j][i].z + *\/ */
/* /\* 		      (rct[k][j][i].z - rct[k-1][j][i].z) * g1 + *\/ */
/* /\* 		      (rct[k-1][j][i].z - rct[k-2][j][i].z) * g2); *\/ */
/* 	      rhs[k][j][i].z = coef * (-    rct[k-1][j][i].z - */
/* 				       2. * rct[k  ][j][i].z + */
/* 				       3. * rct[k+1][j][i].z) + */
/* 		rct[k][j][i].z; */
/* 	    } */
/* 	    else { */
/* 	      rhs[k][j][i].z = 0.5 * (rct[k+1][j][i].z + rct[k][j][i].z); */
/* 	    } */
/* 	  } */
/* 	  else { */
/* 	    if (k < mz-3&&(int)(nvert[k+2][j][i]+0.5)==0) { */
/* /\* 	      c1 = -0.5 * gs[k][j][i].z; *\/ */
/* /\* 	      c2 = -(0.5 * gs[k][j][i].z + gs[k+1][j][i].z); *\/ */
/* /\* 	      c3 = -gs[k+1][j][i].z; *\/ */
/* /\* 	      c4 = -(gs[k][j][i].z + gs[k+1][j][i].z); *\/ */

/* /\* 	      g3 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      c2 = -(0.5 * gs[k][j][i].z); *\/ */
/* /\* 	      c3 = -gs[k+1][j][i].z; *\/ */

/* /\* 	      g4 = c1 * c2 / (c3 * c4); *\/ */

/* /\* 	      rhs[k][j][i].z = (rct[k][j][i].z + *\/ */
/* /\* 		      (rct[k-1][j][i].z - rct[k][j][i].z) * g3 + *\/ */
/* /\* 		      (rct[k][j][i].z   - rct[k+1][j][i].z) * g4); *\/ */
/* 	      rhs[k][j][i].z = coef * (-    rct[k+2][j][i].z - */
/* 				       2. * rct[k+1][j][i].z + */
/* 				       3. * rct[k  ][j][i].z) + */
/* 		rct[k+1][j][i].z; */
/* 	    } */
/* 	    else { */
/* 	      rhs[k][j][i].z = 0.5 * (rct[k+1][j][i].z + rct[k][j][i].z); */
/* 	    } */
/* 	  } */

/* 	  if (j==6 &&(i==1 || i==2)) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD, "rhs%e %e %d\n", conv[k][j][i].z, rhs[k][j][i].z, i); */
/* 	  } */
/* 	} */
/* 	if ((int)(nvert_o[k  ][j][i]+0.5)==0 && */
/* 	    (int)(nvert_o[k+1][j][i]+0.5)==0) { */
	  rhs[k][j][i].z -=
	    (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		     kcsi[k][j][i].y * kzet[k][j][i].y +
		     kcsi[k][j][i].z * kzet[k][j][i].z) +
	     dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		     keta[k][j][i].y * kzet[k][j][i].y +
		     keta[k][j][i].z * kzet[k][j][i].z) +
	     dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		     kzet[k][j][i].y * kzet[k][j][i].y +
		     kzet[k][j][i].z * kzet[k][j][i].z)) * kaj[k][j][i];
/* 	} */
/* 	if (i==10 && j==10) { */
/* 	  PetscPrintf(PETSC_COMM_WORLD, "P%le\n", (dpdc*kcsi[k][j][i].z+dpde*keta[k][j][i].z+dpdz*kzet[k][j][i].z) * kaj[k][j][i]); */
/* 	} */


	//rhs[k][j][i].x = 0.; rhs[k][j][i].y = 0.; rhs[k][j][i].z= 0.;
	if ((int)(nvert[k][j][i]+0.5)!=0) {
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	if ((int)(nvert[k][j][i+1]+0.5)!=0) {
	  rhs[k][j][i].x=0;
	}
	if ((int)(nvert[k][j+1][i]+0.5)!=0) {
	  rhs[k][j][i].y=0;
	}
	if ((int)(nvert[k+1][j][i]+0.5)!=0) {
	  rhs[k][j][i].z=0;
	}
	
      }
    }
  }

  /* i direction boundary conditions*/
  if (xs ==0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (xe == mx) {
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	i = mx-2;
	rhs[k][j][i].x = 0;
	i = mx-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  /* j direction boundary conditions */
  if (ys == 0) {
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	j=0;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (ye == my) {
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	j=my-2;
	rhs[k][j][i].y = 0;
	j=my-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }
  /* k direction boundary conditions */
  if (zs == 0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (ze == mz) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	k=mz-2;
	rhs[k][j][i].z = 0.;//rhs[k-1][j][i].z;
	k=mz-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  DAVecRestoreArray(fda, Csi, &csi);
  DAVecRestoreArray(fda, Eta, &eta);
  DAVecRestoreArray(fda, Zet, &zet);

  DAVecRestoreArray(fda, ICsi, &icsi);
  DAVecRestoreArray(fda, IEta, &ieta);
  DAVecRestoreArray(fda, IZet, &izet);
           
  DAVecRestoreArray(fda, JCsi, &jcsi);
  DAVecRestoreArray(fda, JEta, &jeta);
  DAVecRestoreArray(fda, JZet, &jzet);
           
  DAVecRestoreArray(fda, KCsi, &kcsi);
  DAVecRestoreArray(fda, KEta, &keta);
  DAVecRestoreArray(fda, KZet, &kzet);

  DAVecRestoreArray(fda, Rct, &rct);
  DAVecRestoreArray(fda, Rhs, &rhs);

  DAVecRestoreArray(da, IAj, &iaj);
  DAVecRestoreArray(da, JAj, &jaj);
  DAVecRestoreArray(da, KAj, &kaj);

  DAVecRestoreArray(da, user->lP, &p);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(da, user->lNvert_o, &nvert_o);

  VecWAXPY(dUcont, -1., X, user->Ucont_o);
  
  VecAXPY(Rhs, 1./user->dt, dUcont);

/*   DAVecGetArray(fda, Rhs, &rhs); */
/*   if (ze == mz) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	k=mz-2; */
/* 	rhs[k][j][i].z = rhs[k-1][j][i].z; */
/* 	k=mz-1; */
/* 	rhs[k][j][i].x = 0; */
/* 	rhs[k][j][i].y = 0; */
/* 	rhs[k][j][i].z = 0; */
/*       } */
/*     } */
/*   } */
/*   DAVecRestoreArray(fda, Rhs, &rhs); */
  VecDestroy(dUcont);
  VecDestroy(Rct);
  DARestoreLocalVector(fda, &lUcont);
  return 0;
}


PetscErrorCode Prediction(UserCtx *user)
{
  SNES snes;
  //  KSP ksp;

  VecDuplicate(user->Ucont, &user->Ucont_o);
  VecDuplicate(user->Ucont, &user->Rhs);
  VecCopy(user->Ucont, user->Ucont_o);
  SNESCreate(PETSC_COMM_WORLD, &snes);

  SNESSetApplicationContext(snes, (void*)user);

  InflowFlux(user);
  OutflowVelocity(user, user->Ucont);
  GhostNodeVelocity(user);

  VecCopy(user->Ucont, user->Ucont_o);

  SNESSetFunction(snes, user->Rhs, FormFunctionSNES, (void*)user);

  SNESSetFromOptions(snes);

  PetscInt snesit;
  for (snesit=0; snesit < 1; snesit++) {
    SNESSolve(snes, PETSC_NULL, user->Ucont);
    DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
    Contra2Cart(user);
    OutflowVelocity(user, user->Ucont);
    GhostNodeVelocity(user);
  }
  FormFunctionSNES(snes, user->Ucont, user->Rhs, (void*)user);
  PetscInt itn;
  SNESGetIterationNumber(snes, &itn);
  PetscPrintf(PETSC_COMM_WORLD, "Iteration Number %i\n", itn);

  VecDestroy(user->Ucont_o);
  VecDestroy(user->Rhs);
  //  OutflowVelocity(user);


  DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  InflowFlux(user);
  PetscPrintf(PETSC_COMM_WORLD, "Inflow Flux %e\n", user->FluxInSum);
  OutflowFlux(user);

  OutflowVelocity(user, user->Ucont);

/*   GhostNodeVelocity(user); */
  Contra2Cart(user);


  SNESDestroy(snes);
  return 0;
}
