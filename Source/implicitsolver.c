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

extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed;
extern PetscInt ti,tistart;
extern PetscInt imp_MAX_IT;
extern PetscReal imp_atol, imp_rtol, imp_stol;
extern PetscInt mg_idx, mg_preItr, mg_poItr, mg_MAX_IT;
extern char path[256];
extern IBMNodes	*ibm_ptr;
extern  FSInfo        *fsi_ptr;


extern double time_coeff();
extern void Convert_Phi2_Phi(UserCtx *user);
extern void Convert_Ucont2_Ucont(UserCtx *user);
extern void Convert_RHS_RHS2(UserCtx *user, Vec RHS, Vec RHS2);
extern void Convert_Ucont_Ucont2(UserCtx *user);
extern void IB_BC(UserCtx *user);
extern PetscInt les;


//int momentum_option=0;	//seokkoo

PetscErrorCode ImplicitSolverLHSnew03(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***ucont, ***dtow, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.01;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateMPIAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DAVecGetArray(fda, user->lICsi, &csi);
    DAVecGetArray(fda, user->lIEta, &eta);
    DAVecGetArray(fda, user->lIZet, &zet);
    DAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecGetArray(fda, user->lJCsi, &csi);
    DAVecGetArray(fda, user->lJEta, &eta);
    DAVecGetArray(fda, user->lJZet, &zet);
    DAVecGetArray(da, user->lJAj, &aj);
  } else {
    DAVecGetArray(fda, user->lKCsi, &csi);
    DAVecGetArray(fda, user->lKEta, &eta);
    DAVecGetArray(fda, user->lKZet, &zet);
    DAVecGetArray(da, user->lKAj, &aj);
  }
  //DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(fda, user->psuedot, &dtow);
  DAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = lidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 || nvert[k][j][i] + nvert[k][j][i+1] < 0.1) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].x)>1e-10)
		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;			
	      }	  
	    }//csi0

	    if (fabs(csi1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].y)>1e-10)
		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;			
	      }	  
	    }//csi1

	    if (fabs(csi2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].z)>1e-10)
		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;			
	      }	  
	    }//csi2

	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 || nvert[k][j][i] + nvert[k][j+1][i] < 0.1) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].x)>1e-10)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].y)>1e-10)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)	
		  /eta[k+1][j][i].y;			
	      }	  
	    }//eta1

	    if (fabs(eta2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].z)>1e-10)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;			
	      }	  
	    }//eta2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 || nvert[k][j][i] + nvert[k+1][j][i] < 0.1) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2/Re;

		if (fabs(zet[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2/Re;
		if (fabs(zet[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2/Re;
		if (fabs(zet[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2/Re;
		if (fabs(zet[k][j+1][i].x)>1e-10)
		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2/Re;
		if (fabs(zet[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2/Re;
		if (fabs(zet[k+1][j][i].x)>1e-10)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].y)>1e-10)
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].y)>1e-10)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;			
	      }	  
	    }//zet1

	    if (fabs(zet2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].z)>1e-10)
		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].z)>1e-10)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	      }	  
	    }//zet2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DAVecRestoreArray(fda, user->lICsi, &csi);
    DAVecRestoreArray(fda, user->lIEta, &eta);
    DAVecRestoreArray(fda, user->lIZet, &zet);
    DAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecRestoreArray(fda, user->lJCsi, &csi);
    DAVecRestoreArray(fda, user->lJEta, &eta);
    DAVecRestoreArray(fda, user->lJZet, &zet);
    DAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DAVecRestoreArray(fda, user->lKCsi, &csi);
    DAVecRestoreArray(fda, user->lKEta, &eta);
    DAVecRestoreArray(fda, user->lKZet, &zet);
    DAVecRestoreArray(da, user->lKAj, &aj);
  }
  DAVecRestoreArray(fda, user->lUcont, &ucont);  
  DAVecRestoreArray(fda, user->psuedot, &dtow);
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode ImplicitSolverLHSnew04(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.01;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateMPIAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DAVecGetArray(fda, user->lICsi, &csi);
    DAVecGetArray(fda, user->lIEta, &eta);
    DAVecGetArray(fda, user->lIZet, &zet);
    DAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecGetArray(fda, user->lJCsi, &csi);
    DAVecGetArray(fda, user->lJEta, &eta);
    DAVecGetArray(fda, user->lJZet, &zet);
    DAVecGetArray(da, user->lJAj, &aj);
  } else {
    DAVecGetArray(fda, user->lKCsi, &csi);
    DAVecGetArray(fda, user->lKEta, &eta);
    DAVecGetArray(fda, user->lKZet, &zet);
    DAVecGetArray(da, user->lKAj, &aj);
  }
  //DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(fda, user->psuedot, &dtow);
  DAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = lidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].x)>1e-10)
		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;			
	      }	  
	    }//csi0

	    if (fabs(csi1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1|| (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].y)>1e-10)
		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;			
	      }	  
	    }//csi1

	    if (fabs(csi2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ( i==1 || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].z)>1e-10)
		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;			
	      }	  
	    }//csi2

	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].y+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].x)>1e-10)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].y)>1e-10)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)	
		  /eta[k+1][j][i].y;			
	      }	  
	    }//eta1

	    if (fabs(eta2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].z)>1e-10)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.; 
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;			
	      }	  
	    }//eta2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].z+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2/Re;

		if (fabs(zet[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2/Re;
		if (fabs(zet[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2/Re;
		if (fabs(zet[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2/Re;
		if (fabs(zet[k][j+1][i].x)>1e-10)
		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2/Re;
		if (fabs(zet[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2/Re;
		if (fabs(zet[k+1][j][i].x)>1e-10)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].y)>1e-10)
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].y)>1e-10)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;			
	      }	  
	    }//zet1

	    if (fabs(zet2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].z)>1e-10)
		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].z)>1e-10)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	      }	  
	    }//zet2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DAVecRestoreArray(fda, user->lICsi, &csi);
    DAVecRestoreArray(fda, user->lIEta, &eta);
    DAVecRestoreArray(fda, user->lIZet, &zet);
    DAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecRestoreArray(fda, user->lJCsi, &csi);
    DAVecRestoreArray(fda, user->lJEta, &eta);
    DAVecRestoreArray(fda, user->lJZet, &zet);
    DAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DAVecRestoreArray(fda, user->lKCsi, &csi);
    DAVecRestoreArray(fda, user->lKEta, &eta);
    DAVecRestoreArray(fda, user->lKZet, &zet);
    DAVecRestoreArray(da, user->lKAj, &aj);
  }
  DAVecRestoreArray(fda, user->lUcont, &ucont);  
  DAVecRestoreArray(fda, user->psuedot, &dtow);
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode ImplicitSolverLHSnew05(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir)
{ // works for Stretched grid also
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.01;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateMPIAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DAVecGetArray(fda, user->lICsi, &csi);
    DAVecGetArray(fda, user->lIEta, &eta);
    DAVecGetArray(fda, user->lIZet, &zet);
    DAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecGetArray(fda, user->lJCsi, &csi);
    DAVecGetArray(fda, user->lJEta, &eta);
    DAVecGetArray(fda, user->lJZet, &zet);
    DAVecGetArray(da, user->lJAj, &aj);
  } else {
    DAVecGetArray(fda, user->lKCsi, &csi);
    DAVecGetArray(fda, user->lKEta, &eta);
    DAVecGetArray(fda, user->lKZet, &zet);
    DAVecGetArray(da, user->lKAj, &aj);
  }
  //DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = lidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].x)>1e-10)
		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;			
	      }	  
	    }//csi0

	    if (fabs(csi1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1|| (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].y)>1e-10)
		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;			
	      }	  
	    }//csi1

	    if (fabs(csi2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ( i==1 || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].z)>1e-10)
		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;			
	      }	  
	    }//csi2

	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].x)>1e-10)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].y)>1e-10)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)	
		  /eta[k+1][j][i].y;			
	      }	  
	    }//eta1

	    if (fabs(eta2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].z)>1e-10)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.; 
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;			
	      }	  
	    }//eta2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (user->bctype[5]==8 && k==mz-2) {
	    val[0]=COEF_TIME_ACCURACY/dt + 0.8; 
	    val[1]= -.8;
	    idx[0]=row;
	    idx[1]=lidx(i  , j ,k-1, user);

	    MatSetValues(user->A, 1, &row, 2, idx, val, INSERT_VALUES);	  	    

	  } else if (k<mz-2 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2/Re;

		if (fabs(zet[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2/Re;
		if (fabs(zet[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2/Re;
		if (fabs(zet[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2/Re;
		if (fabs(zet[k][j+1][i].x)>1e-10)
		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2/Re;
		if (fabs(zet[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2/Re;
		if (fabs(zet[k+1][j][i].x)>1e-10)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].y)>1e-10)
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].y)>1e-10)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;			
	      }	  
	    }//zet1

	    if (fabs(zet2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].z)>1e-10)
		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].z)>1e-10)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	      }	  
	    }//zet2


	    idx[0]=lidx(i  , j , k,  user);           
	    idx[1]=lidx(i-1, j , k,  user);
	    idx[2]=lidx(i+1, j , k,  user);
	    idx[3]=lidx(i  ,j-1, k , user);
	    idx[4]=lidx(i  ,j+1, k , user);
	    idx[5]=lidx(i  , j ,k-1, user);
	    idx[6]=lidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DAVecRestoreArray(fda, user->lICsi, &csi);
    DAVecRestoreArray(fda, user->lIEta, &eta);
    DAVecRestoreArray(fda, user->lIZet, &zet);
    DAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DAVecRestoreArray(fda, user->lJCsi, &csi);
    DAVecRestoreArray(fda, user->lJEta, &eta);
    DAVecRestoreArray(fda, user->lJZet, &zet);
    DAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DAVecRestoreArray(fda, user->lKCsi, &csi);
    DAVecRestoreArray(fda, user->lKEta, &eta);
    DAVecRestoreArray(fda, user->lKZet, &zet);
    DAVecRestoreArray(da, user->lKAj, &aj);
  }
  DAVecRestoreArray(fda, user->lUcont, &ucont);  
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}


PetscErrorCode GetPsuedoTime(UserCtx *user)
{ 
  DA da = user->da, fda = user->fda; 
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

  PetscInt      i, j, k;
  Cmpnts        ***pdt, ***ucont, ru, absUcont;
  Cmpnts	***csi, ***eta, ***zet;
  PetscReal     ***iaj,***jaj,***kaj;
  PetscReal     g11,g22,g33;
  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal     cfl=user->cfl, vnn=user->vnn*user->ren;

  DAVecGetArray(fda, user->psuedot, &pdt);
  DAVecGetArray(fda, user->lICsi, &csi);
  DAVecGetArray(fda, user->lJEta, &eta);
  DAVecGetArray(fda, user->lKZet, &zet);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	csi0 = csi[k][j][i].x*iaj[k][j][i];
	csi1 = csi[k][j][i].y*iaj[k][j][i];
	csi2 = csi[k][j][i].z*iaj[k][j][i];
	
	eta0 = eta[k][j][i].x*jaj[k][j][i];
	eta1 = eta[k][j][i].y*jaj[k][j][i];
	eta2 = eta[k][j][i].z*jaj[k][j][i];
	
	zet0 = zet[k][j][i].x*kaj[k][j][i];
	zet1 = zet[k][j][i].y*kaj[k][j][i];
	zet2 = zet[k][j][i].z*kaj[k][j][i];
	
	g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;	   
	g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;	   	    
	g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
	
	absUcont.x=fabs(ucont[k][j][i].x*iaj[k][j][i]);
	absUcont.y=fabs(ucont[k][j][i].y*jaj[k][j][i]);
	absUcont.z=fabs(ucont[k][j][i].z*kaj[k][j][i]);

	ru.x=absUcont.x+sqrt(absUcont.x*absUcont.x+g11);
	ru.y=absUcont.y+sqrt(absUcont.y*absUcont.y+g22);
	ru.z=absUcont.z+sqrt(absUcont.z*absUcont.z+g33);
	
/* 	pdt[k][j][i].x= PetscMin(cfl/ru.x/iaj[k][j][i],vnn/g11/iaj[k][j][i]); */
/* 	pdt[k][j][i].y= PetscMin(cfl/ru.y/jaj[k][j][i],vnn/g22/jaj[k][j][i]); */
/* 	pdt[k][j][i].z= PetscMin(cfl/ru.z/kaj[k][j][i],vnn/g33/kaj[k][j][i]); */

	pdt[k][j][i].x= PetscMin(cfl/ru.x,vnn/g11);
	pdt[k][j][i].x= PetscMin(pdt[k][j][i].x,user->dt);
	pdt[k][j][i].y= PetscMin(cfl/ru.y,vnn/g22);
	pdt[k][j][i].y= PetscMin(pdt[k][j][i].y,user->dt);
	pdt[k][j][i].z= PetscMin(cfl/ru.z,vnn/g33);
	pdt[k][j][i].z= PetscMin(pdt[k][j][i].z,user->dt);
	//PetscPrintf(PETSC_COMM_WORLD, "Psuedo Time %le %le %le %le %d %d %d\n", pdt[k][j][i].x, absUcont.x, csi0 ,ru.x,i,j,k);
	
      }
    }
  }
  DAVecRestoreArray(fda, user->psuedot, &pdt);
  DAVecRestoreArray(fda, user->lICsi, &csi);
  DAVecRestoreArray(fda, user->lJEta, &eta);
  DAVecRestoreArray(fda, user->lKZet, &zet);
  DAVecRestoreArray(fda, user->lUcont, &ucont);  
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);

  return(0);
}

PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm,
				      FSInfo *fsi)
{
  DA            da, fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];  alfa[0]=1.;//0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec           B;
  Vec           Ucont_i, RB;
  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ***nvert;
  PetscReal     ts,te,cput,cfl_i;
  //Vec    dUcont, pUcont;
  PetscInt	bi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;

  if (block_number>1) {
    Block_Interface_U(user);
  }

  for (bi=0; bi<block_number; bi++) {

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;
    ksp = user[bi].ksp;
    cfl_i= user[bi].cfl;
    
    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]),&fsi[0]);

    //ibm_interpolation_advanced2(&user[bi], ibm);
    if (immersed && ti>0 ) {
      /*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
	ibm_interpolation_advanced(&user[bi]);
      }
    }

/*     if (!CGSolver) { */
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "CGS!\n" ); */
/*     } */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  

  //if (ti>0) {
    //Get the direction 
    //0==i dir, 1==j dir, 2==k dir
/*     for (dir=0;dir<3;dir++){ */


    PetscInt pseudot;
    PetscReal normdU=10.,normdU1,reldU=1.,normdT, normF;
    // pseudo time iteration
    //for(pseudot=0; pseudot<5;pseudot++) {
    PetscGetTime(&ts);
    pseudot=0;
    while (( (normdU>imp_atol && reldU>imp_rtol) || pseudot<1) && pseudot<imp_MAX_IT) {
      pseudot++;
      GetPsuedoTime(&(user[bi]));
      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	VecSet(user[bi].Rhs,0);
	Formfunction_2(&user[bi], user[bi].Rhs, 1.0);
	//FormFunction1(user[bi].Ucont, user[bi].Rhs, &(user[bi]));
	// Calculate du/dt
	// du = 1.5 u - 2 u_o + 0.5 u_rm1
	if (COEF_TIME_ACCURACY>1.) {
	  VecCopy(user[bi].Ucont, user[bi].dUcont);
	  VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	  VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	  VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);
	} else {
	  VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);
	}

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 
	//VecWAXPY(B, alfa[istage]* user[bi].dt * user[bi].st  , user[bi].Rhs, user[bi].Ucont_o);
	// Calc B=RHS
	VecCopy(user[bi].Rhs, B);
	//VecScale(B,alfa[istage]*user[bi].dt* user[bi].st);
	VecSet(user[bi].dUcont, 0.);
	PetscPrintf(PETSC_COMM_WORLD, "RHS done!\n" );

	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){

	  // Set the LHS	 
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);

	  // set rhs of solver
	  DAVecGetArray(fda, B, &b);
	  DAVecGetArray(da, RB, &rb);
	  DAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;			  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || i>mx-3 || i<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].x=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || j>my-3 || j<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].y=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || k>mz-3 || k<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].z=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;	
		  }
		}
	      }
	    }
	  }

	  DAVecRestoreArray(fda, B, &b);
	  DAVecRestoreArray(da, RB, &rb);
	  DAVecRestoreArray(da, user[bi].lNvert, &nvert);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
	  KSPSetType(ksp, KSPBCGS);
/* 	  KSPSetType(ksp, KSPGMRES); */

/* 	  VecNorm(RB, NORM_INFINITY, &normF); */
/* 	  /\* normF=PetscMin(normF,.9); *\/ */
/* 	  normF=PetscMax(normF,1.e-16); */
/* 	  rtol=normF*normF;//PetscMin(normF*normF,1.e-1); */
/* 	  KSPSetTolerances(ksp,normF*normF,rtol,PETSC_DEFAULT,400); */
	  
	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  PetscBarrier(PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);
/* 	  VecNorm(Ucont_i, NORM_INFINITY, &normdU); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, "norm dU %le dir %d!\n",normdU,dir); */

	  // Restore
	  DAVecGetArray(da, Ucont_i, &ucont_i);
	  //DAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DAVecGetArray(fda, user[bi].dUcont, &ucont);
	  
	  //PetscPrintf(PETSC_COMM_SELF, "du.i=dui!\n");
	  //PetscBarrier(PETSC_NULL);
	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    
/* 		  if ((j==97||j==98)  && (k==33 || k==32 || k==31 ||k==30 || k==29 || k==10) ){ */
/* 		    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
/* 		  } */
		  
		}
	      }
	    }
	  }
	  //PetscPrintf(PETSC_COMM_SELF, "End Restore!\n");
	  //DAVecRestoreArray(fda, B, &b);
	  
	  DAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  //DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  //DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  
	  //VecAssemblyBegin(user->Ucont);
	  //VecAssemblyEnd(user->Ucont);
	  
	  //MatNullSpaceDestroy(nullsp);
	  PetscTruth temp;
	  MatValid(user[bi].A, &temp);
	  if (temp) {
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(user[bi].A);
	    //PetscPrintf(PETSC_COMM_WORLD, "Destroy A!\n");
	  }
	  
	}//dir
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	PetscBarrier(PETSC_NULL);
	DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU1=normdU;
	if (pseudot>1) reldU=normdU/normdU1;
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	VecNorm(B, NORM_INFINITY, &normF);
	PetscGetTime(&te);
	cput=te-ts;
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le %le %le\n", pseudot, normdU, reldU,normF, normdT, cput);

	if (!rank) {
	  FILE *f;
	  char filen[80];
	  
	  sprintf(filen, "Converge_dU");
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU, reldU, normF,normdT, cput);
	  fclose(f);
	}

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum);
	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum);
	
	
	//DAVecRestoreArray(fda, Ucont_o, &ucont);
	//VecCopy(Ucont_o, user->Ucont);
	/*   DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	/*   DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	
	// Apply BC
	       
	FormBCS(&(user[bi]),&fsi[0]);
      } // istage
            
      if (immersed && ti>0) {
	/*for (ibi=0;ibi<NumberOfBodies;ibi++)*/ {
	  ibm_interpolation_advanced(&user[bi]);
	}

      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      

    } //psuedo time
  //  } // ti>0

/*     }//dir */

  // Destroy   
    KSPDestroy(user[bi].ksp);
    VecDestroy(user[bi].Rhs);

    VecDestroy(B);
    VecDestroy(Ucont_i);
    VecDestroy(RB);
    VecDestroy(user[bi].dUcont);
    VecDestroy(user[bi].pUcont);

  } //bi

  if (block_number>1) {
    Block_Interface_U(user);
  }

  // Destroy
/*   for (bi=0; bi<block_number; bi++) { */
/*     PetscTruth temp; */
/*     MatValid(user[bi].A, &temp); */
/*     if (temp) { */
/*       user[bi].assignedA = PETSC_FALSE; */
/*       MatDestroy(user[bi].A); */
/*     }     */
/*   } */

  return(0);
}

PetscErrorCode ImplicitMomentumSolver1(UserCtx *user, IBMNodes *ibm,
				      FSInfo *fsi)
{
  DA            da, fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];  alfa[0]=1.;//0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec           B;
  Vec           Ucont_i, RB;
  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ***nvert;
  PetscReal     ts,te,cput,cfl_i;
  //Vec    dUcont, pUcont;
  PetscInt	bi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;

  if (block_number>1) {
    Block_Interface_U(user);
  }

  for (bi=0; bi<block_number; bi++) {

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;
    ksp = user[bi].ksp;
    cfl_i= user[bi].cfl;
    
    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]),&fsi[0]);

    if (immersed && ti>0 ) {
      {
	ibm_interpolation_advanced(&user[bi]);
      }
    }

/*     if (!CGSolver) { */
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "CGS!\n" ); */
/*     } */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  

  //if (ti>0) {
    //Get the direction 
    //0==i dir, 1==j dir, 2==k dir
/*     for (dir=0;dir<3;dir++){ */


    PetscInt pseudot,ls_itr;
    PetscReal normdU=10.,normdU_old,normdU1,reldU=1.,normdT, normF;
    PetscScalar lambda=.5;

    // pseudo time iteration
    //for(pseudot=0; pseudot<5;pseudot++) {
    PetscGetTime(&ts);
    pseudot=0;
    ls_itr=0;
    while (( (normdU>imp_atol && reldU>imp_rtol) || pseudot<1) && pseudot<imp_MAX_IT) {
      pseudot++;
      GetPsuedoTime(&(user[bi]));
      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	VecSet(user[bi].Rhs,0);
	Formfunction_2(&user[bi], user[bi].Rhs, 1.0);
	//FormFunction1(user[bi].Ucont, user[bi].Rhs, &(user[bi]));
	// Calculate du/dt
	// du = 1.5 u - 2 u_o + 0.5 u_rm1
	if (COEF_TIME_ACCURACY>1.) {
	  VecCopy(user[bi].Ucont, user[bi].dUcont);
	  VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	  VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	  VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);
	} else {
	  VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);
	}

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 
	//VecWAXPY(B, alfa[istage]* user[bi].dt * user[bi].st  , user[bi].Rhs, user[bi].Ucont_o);
	// Calc B=RHS
	VecCopy(user[bi].Rhs, B);
	//VecScale(B,alfa[istage]*user[bi].dt* user[bi].st);
	VecSet(user[bi].dUcont, 0.);
	PetscPrintf(PETSC_COMM_WORLD, "RHS done!\n" );

	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){

	  // Set the LHS	 
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);

	  // set rhs of solver
	  DAVecGetArray(fda, B, &b);
	  DAVecGetArray(da, RB, &rb);
	  DAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;			  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || i>mx-3 || i<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].x=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || j>my-3 || j<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].y=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || k>mz-3 || k<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].z=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;	
		  }
		}
	      }
	    }
	  }

	  DAVecRestoreArray(fda, B, &b);
	  DAVecRestoreArray(da, RB, &rb);
	  DAVecRestoreArray(da, user[bi].lNvert, &nvert);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
	  KSPSetType(ksp, KSPBCGS);
/* 	  KSPSetType(ksp, KSPGMRES); */

/* 	  VecNorm(RB, NORM_INFINITY, &normF); */
/* 	  /\* normF=PetscMin(normF,.9); *\/ */
/* 	  normF=PetscMax(normF,1.e-16); */
/* 	  rtol=normF*normF;//PetscMin(normF*normF,1.e-1); */
/* 	  KSPSetTolerances(ksp,normF*normF,rtol,PETSC_DEFAULT,400); */
	  
	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  PetscBarrier(PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);
/* 	  VecNorm(Ucont_i, NORM_INFINITY, &normdU); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, "norm dU %le dir %d!\n",normdU,dir); */

	  // Restore
	  DAVecGetArray(da, Ucont_i, &ucont_i);
	  //DAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DAVecGetArray(fda, user[bi].dUcont, &ucont);
	  
	  //PetscPrintf(PETSC_COMM_SELF, "du.i=dui!\n");
	  //PetscBarrier(PETSC_NULL);
	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    
/* 		  if ((j==97||j==98)  && (k==33 || k==32 || k==31 ||k==30 || k==29 || k==10) ){ */
/* 		    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
/* 		  } */
		  
		}
	      }
	    }
	  }
	  //PetscPrintf(PETSC_COMM_SELF, "End Restore!\n");
	  //DAVecRestoreArray(fda, B, &b);
	  
	  DAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  //DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  //DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  
	  //VecAssemblyBegin(user->Ucont);
	  //VecAssemblyEnd(user->Ucont);
	  
	  //MatNullSpaceDestroy(nullsp);
	  PetscTruth temp;
	  MatValid(user[bi].A, &temp);
	  if (temp) {
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(user[bi].A);
	    //PetscPrintf(PETSC_COMM_WORLD, "Destroy A!\n");
	  }
	  
	}//dir
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	PetscBarrier(PETSC_NULL);
	DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	normdU_old=normdU;
	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU1=normdU;
	if (pseudot>1) reldU=normdU/normdU1;
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	VecNorm(B, NORM_INFINITY, &normF);
	PetscGetTime(&te);
	cput=te-ts;
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le %le %le\n", pseudot, normdU, reldU,normF, normdT, cput);

	if (!rank) {
	  FILE *f;
	  char filen[80];
	  //sprintf(filen, "Converge_dU%1.1d",dir);
	  sprintf(filen, "Converge_dU");
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU, reldU, normF,normdT, cput);
	  fclose(f);
	}

	lambda=0.5;
	if (normdU>normdU_old && ls_itr<10 && pseudot>1) {
	  ls_itr++ ;
	  pseudot--;
	  normdU=normdU_old;

	  user[bi].cfl *=lambda;

	  VecCopy(user[bi].pUcont, user[bi].Ucont);
	  
	  DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  	  
	  PetscPrintf(PETSC_COMM_WORLD, "LineSearch  %d |F| %le |Fold| %le %le %le\n",ls_itr,normdU,normdU_old,normF,lambda);
	} else if (pseudot>15 && normdU<imp_atol*5. && ls_itr>2) {
	   user[bi].cfl /=lambda;
	}

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum);
	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum);
	
	
	//DAVecRestoreArray(fda, Ucont_o, &ucont);
	//VecCopy(Ucont_o, user->Ucont);
	/*   DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	/*   DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	
	// Apply BC
	       
	FormBCS(&(user[bi]),&fsi[0]);
      } // istage
            
      if (immersed && ti>0) {
	/*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
	  ibm_interpolation_advanced(&user[bi]);
	}

	//ibm_interpolation_advanced(&user[bi], ibm, fsi);
	//ibm_interpolation_advanced2(&user[bi], ibm);
/* 	ibm_interpolation(ibminfo, user, ibm); */
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      

    } //psuedo time
  //  } // ti>0

/*     }//dir */


    user[bi].cfl=cfl_i;

    // Destroy   
    KSPDestroy(user[bi].ksp);
    VecDestroy(user[bi].Rhs);

    VecDestroy(B);
    VecDestroy(Ucont_i);
    VecDestroy(RB);
    VecDestroy(user[bi].dUcont);
    VecDestroy(user[bi].pUcont);

  } //bi

  if (block_number>1) {
    Block_Interface_U(user);
  }

  return(0);
}

PetscErrorCode CalcRHS(UserCtx *user, int dudt)
{
  PetscInt      i, j, k;
  DALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Cmpnts        ***rhs;
  PetscReal     ***nvert;
  Vec           dUcont;

  DA            da = user->da,fda = user->fda;
  info = user->info;
  
  xs = info.xs; xe = info.xs + info.xm;
  ys = info.ys; ye = info.ys + info.ym;
  zs = info.zs; ze = info.zs + info.zm;
  mx = info.mx; my = info.my; mz = info.mz;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscInt     mz_end;
  if (user->bctype[5]==8) 
    mz_end=mz-2;
  else
    mz_end=mz-3;
   
  
	double threshold=0.1;   //03.29
	
	
  
	VecDuplicate(user->Ucont, &dUcont);

  // Calculate du/dt
  // du = 1.5 u - 2 u_o + 0.5 u_rm1
  
  double coeff = time_coeff();
  
  if (coeff>1.4 && coeff<1.6) {
	VecCopy(user->Ucont, dUcont);
	VecScale(dUcont, 1.5);
	VecAXPY(dUcont, -2.,user->Ucont_o );
	VecAXPY(dUcont, 0.5, user->Ucont_rm1);
  }
  /*
  else if (COEF_TIME_ACCURACY>1.8 && COEF_TIME_ACCURACY<1.9) {
		VecCopy(user->Ucont, dUcont);
		VecScale(dUcont, 11./6.);
		VecAXPY(dUcont, -3.,user->Ucont_o );
		VecAXPY(dUcont, 1.5, user->Ucont_rm1);
		VecAXPY(dUcont, -1./3., user->Ucont_rm2);
  }*/

  else {
    VecWAXPY(dUcont, -1., user->Ucont_o, user->Ucont);
  }

      // Calc the RHS
	VecSet(user->Rhs,0);
	Formfunction_2(user, user->Rhs, 1.0);
  //FormFunction1(user->Ucont, user->Rhs, user);
  //extern PetscErrorCode Formfunction_2(UserCtx *user, Vec Rhs);
  //Formfunction_2(user, user->Rhs);
    
  // Add -du/dt to right hand side
  VecAXPY(user->Rhs, -1./user->dt, dUcont);

  // set rhs BC
  DAVecGetArray(fda, user->Rhs, &rhs);//
  DAVecGetArray(da, user->lNvert, &nvert);//
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
		
		if ( (i==mx-2 && !i_periodic && !ii_periodic) || (nvert[k][j][i]+nvert[k][j][i+1])>threshold ) {
			rhs[k][j][i].x=0.;
		} 
		
		if ( (j==my-2 && !j_periodic && !jj_periodic) || (nvert[k][j][i]+nvert[k][j+1][i])>threshold ) {
			rhs[k][j][i].y=0.;
		} 
		
		if ( (k==mz-2 && !k_periodic && !kk_periodic) || (nvert[k][j][i]+nvert[k+1][j][i])>threshold ) {
			rhs[k][j][i].z=0.;
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
	if(!i_periodic && !ii_periodic) {
		i = mx-2;
		rhs[k][j][i].x = 0;
	}
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
	if(!j_periodic && !jj_periodic) {
		j=my-2;
		rhs[k][j][i].y = 0;
	}
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
	if(!k_periodic && !kk_periodic) {
		k=mz-2;
		rhs[k][j][i].z = 0;
	}
	k=mz-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  DAVecRestoreArray(fda, user->Rhs, &rhs);//
  DAVecRestoreArray(da, user->lNvert, &nvert);//
  

  VecDestroy(dUcont);

  return(0);
}

PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscPrintf(PETSC_COMM_WORLD,"     (%D) SNES Residual norm %14.12e \n",n,rnorm);
	return 0;
}

PetscErrorCode FormFunction_seokkoo(SNES snes, Vec Ucont, Vec Rhs, void *ptr) //UserCtx *user)
{
	UserCtx *user = (UserCtx*)ptr;
	VecCopy(Ucont, user->Ucont);
	
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	
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
	
	
	DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	
	Contra2Cart_2(user);
	IB_BC(user);
	
	/*
			0      1      2                     mx-2  mx-1
		|	|   1	|   2	|   3	|  ...	| mx-2| mx-1 |
	*/
			
	CalcRHS(user, 1);
	VecCopy(user->Rhs, Rhs);
	
	//VecScale(Rhs,user->dt);	// 11.9.2008
	//VecScale(Rhs,user->ren);	// 12.9.2008
	return(0);
}

PetscErrorCode FormFunction_seokkoo2(SNES snes, Vec Ucont2, Vec Rhs2, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;
	
	VecCopy(Ucont2, user->Ucont2);
	Convert_Ucont2_Ucont(user);
	
	CalcRHS(user, 1);
	Convert_RHS_RHS2(user, user->Rhs, Rhs2);
	
	return(0);
}


int snes_created=0;

PetscErrorCode Implicit_MatrixFree(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{

	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	
	PetscInt	bi;
	

	int rank;
	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	
	for (bi=0; bi<block_number; bi++) {
		
		KSP ksp;
		PC pc;
		
		VecDuplicate(user[bi].Ucont, &user[bi].rhs);
		if(!snes_created) {
			SNESCreate(PETSC_COMM_WORLD,&user[bi].snes);
			SNESMonitorSet (user[bi].snes, MySNESMonitor, PETSC_NULL, PETSC_NULL);
		}
		SNESSetFunction(user[bi].snes,user[bi].rhs,FormFunction_SNES,(void *)&user[bi]);
		
		#if defined(PETSC_HAVE_ADIC___)
		DAGetMatrix(user->da,MATAIJ,&user[bi].J);
		ISColoring             iscoloring;
		DAGetColoring(user->da,IS_COLORING_GHOSTED,&iscoloring);
		MatSetColoring(J[bi],iscoloring);
		ISColoringDestroy(iscoloring);
		SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,SNESDAComputeJacobianWithAdic,(void *)&user[bi]);
		#else
		MatCreateSNESMF(user[bi].snes, &user[bi].J);
		SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,MatMFFDComputeJacobian,(void *)&user[bi]);
		#endif
		
		extern double imp_free_tol;
		SNESSetType(user[bi].snes, SNESTR);			//SNESTR,SNESLS	: SNESLS is better for stiff PDEs such as the one including IB but slower

		SNESSetMaxLinearSolveFailures(user[bi].snes,10000);
		SNESSetMaxNonlinearStepFailures(user[bi].snes,10000);		
		SNESKSPSetUseEW(user[bi].snes, PETSC_TRUE);
		SNESKSPSetParametersEW(user[bi].snes,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		//SNESSetTolerances(user[bi].snes,imp_free_tol,PETSC_DEFAULT,PETSC_DEFAULT,50,50000);
		SNESSetTolerances(user[bi].snes, /*5.e-4*/PETSC_DEFAULT,imp_free_tol,PETSC_DEFAULT,50,50000);
		
		SNESGetKSP(user[bi].snes, &ksp);
		KSPGetPC(ksp,&pc);
		
		if(!snes_created) {
			KSPSetType(ksp, KSPGMRES);
			//KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);	//2009.09.22 poor performance
			//KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);	//2009.09.22
			
			//KSPFischerGuess itg;
			//KSPFischerGuessCreate(ksp,1,100,&itg);
			//KSPSetFischerGuess(ksp, itg);	//2009.09.22
			
			//KSPGMRESSetPreAllocateVectors(ksp);	--> crazy thing consumes memory 
		}
		#if defined(PETSC_HAVE_ADIC___)
		PCSetType(pc,PCBJACOBI);
		#else
		PCSetType(pc,PCNONE);
		#endif
		int maxits=1000;	
		//double rtol=PETSC_DEFAULT, atol=imp_free_tol, dtol=PETSC_DEFAULT;
		double rtol=imp_free_tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;

		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
		
		snes_created=1;
	}
	
	bi =0;
	//int it=0;
		
	Vec U;
	VecDuplicate(user[bi].Ucont, &U);
	VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
	extern int outflow_scale;
		
		
	InflowFlux(&(user[bi]));
	outflow_scale=0;
	FormBCS(&(user[bi]),&fsi[0]);
		
	VecCopy(user[bi].Ucont, U);
		
	SNESSolve(user[bi].snes, PETSC_NULL, U);
		
	PetscPrintf(PETSC_COMM_WORLD, "\nMomentum eqs computed ...\n");
	
	double norm;
	SNESGetFunctionNorm(user[bi].snes, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "\nSNES residual norm=%.5e\n\n", norm);
		
	VecCopy(U, user[bi].Ucont);
		
	DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	
	Vec Coor;
	Cmpnts ***ucont, ***lucont, ***coor, ***csi, ***eta, ***ucat, ***lucat, ***zet;
	PetscReal ***level;
	int i, j, k;
	
	DAGetGhostedCoordinates(user[0].da, &Coor);
	if(levelset) DAVecGetArray(user[0].da, user[0].lLevelset, &level);
	DAVecGetArray(user[0].fda, Coor, &coor);
	DAVecGetArray(user[0].fda, user[0].Ucont, &ucont);
	DAVecGetArray(user[0].fda, user[0].Ucat, &ucat);
	DAVecGetArray(user[0].fda, user[0].lUcont, &lucont);
	DAVecGetArray(user[0].fda, user[0].lUcat, &lucat);
	DAVecGetArray(user[0].fda, user[0].lCsi, &csi);
	DAVecGetArray(user[0].fda, user[0].lEta, &eta);
	DAVecGetArray(user[0].fda, user[0].lZet, &zet);
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
	  
		if ( (user[0].bctype[0]==1 || user[0].bctype[0]==-1 || user[0].bctype[0]==-2 || user[0].bctype[0]==10) && i==0) ucont[k][j][i].x = 0;
		if ( (user[0].bctype[1]==1 || user[0].bctype[1]==-1 || user[0].bctype[1]==-2 || user[0].bctype[1]==10) && i==mx-2) ucont[k][j][i].x = 0;
		if ( (user[0].bctype[2]==1 || user[0].bctype[2]==-1 || user[0].bctype[2]==-2 || user[0].bctype[2]==10) && j==0) ucont[k][j][i].y = 0;
		if ( (user[0].bctype[3]==1 || user[0].bctype[3]==-1 || user[0].bctype[3]==-2 || user[0].bctype[3]==10 || user[0].bctype[3]==2) && j==my-2) ucont[k][j][i].y = 0;
		if ( (user[0].bctype[4]==1 || user[0].bctype[4]==-1 || user[0].bctype[4]==-2 || user[0].bctype[4]==10) && k==0) ucont[k][j][i].z = 0;
		if ( (user[0].bctype[5]==1 || user[0].bctype[5]==-1 || user[0].bctype[5]==-2 || user[0].bctype[5]==10) && k==mz-2) ucont[k][j][i].z = 0;

//add (toni)
	//for air_flow_levelset
		if(air_flow_levelset && k==0 && level[k][j][i]>-dthick)ucont[k][j][i].z = 0;
		if(air_flow_levelset && k==mz-2 && level[k][j][i]>-dthick)ucont[k][j][i].z = 0;
//end (toni)
		
		if ( (user[0].bctype[0]==-1 || user[0].bctype[0]==-2) && i==1) ucont[k][j][i].x = 0;
		if ( (user[0].bctype[1]==-1 || user[0].bctype[1]==-2) && i==mx-3) ucont[k][j][i].x = 0;
		if ( (user[0].bctype[2]==-1 || user[0].bctype[2]==-2) && j==1) ucont[k][j][i].y = 0;
		if ( (user[0].bctype[3]==-1 || user[0].bctype[3]==-2) && j==my-3) ucont[k][j][i].y = 0;
		if ( (user[0].bctype[4]==-1 || user[0].bctype[4]==-2) && k==1) ucont[k][j][i].z = 0;
		if ( (user[0].bctype[5]==-1 || user[0].bctype[5]==-2) && k==mz-3) ucont[k][j][i].z = 0;
		
		if ( user[0].bctype[4]==5 && k==0 ) {
		  if(levelset) {
		    //if(level[k][j][i]>0) ucont[k][j][i].y = 0;
		    ucat[k][j][i].x = - lucat[k+1][j][i].x;
		    ucat[k][j][i].y = - lucat[k+1][j][i].y;
                    ucat[k][j][i].z = - lucat[k+1][j][i].z;
		    lucat[k][j][i] = ucat[k][j][i];
		  }
		}

		if ( user[0].bctype[3]==4 && j==my-2 ) {
		  ucat[k][j+1][i] = ucat[k][j][i];
		  lucat[k][j+1][i] = ucat[k][j+1][i];
		  ucont[k][j][i].y = 0.5*(lucat[k][j][i].x+lucat[k][j+1][i].x)*eta[k][j][i].x + 0.5*(lucat[k][j][i].y+lucat[k][j+1][i].y)*eta[k][j][i].y + 0.5*(lucat[k][j][i].z+lucat[k][j+1][i].z)*eta[k][j][i].z;
                }
		
		if ( user[0].bctype[5]==4 && k==mz-2) {
		  ucat[k+1][j][i] = ucat[k][j][i]; 
		  lucat[k+1][j][i] = ucat[k+1][j][i];
		  ucont[k][j][i].z = 0.5*(lucat[k][j][i].x+lucat[k+1][j][i].x)*zet[k][j][i].x + 0.5*(lucat[k][j][i].y+lucat[k+1][j][i].y)*zet[k][j][i].y + 0.5*(lucat[k][j][i].z+lucat[k+1][j][i].z)*zet[k][j][i].z;
		  if(/*!levelset &&*/ ucont[k][j][i].z<0) {
		    ucont[k][j][i].z = 0;
		    ucat[k][j][i].x = ucat[k][j][i].y = ucat[k][j][i].z = 0;
		    ucat[k+1][j][i] = ucat[k][j][i];
		  }
			if (air_flow_levelset && level[k][j][i]>-dthick) {
				//ucat[k+1][j][i] = - ucat[k][j][i]; 
				//lucat[k+1][j][i] =  ucat[k+1][j][i];		
				ucont[k][j][i].z = 0.;
			}			
		}
		
		if (user->bctype[0]==11 && i==0 && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) {
			double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
			if( zc > 0 ) {
				ucont[k][j][i].x = 0 * csi[k][j][i].x + 0 * csi[k][j][i].y + lucat[k][j][i].z * csi[k][j][i].z;
			}
		}
	}
	if(levelset) DAVecRestoreArray(user[0].da, user[0].lLevelset, &level);
	DAVecRestoreArray(user[0].fda, Coor, &coor);
	DAVecRestoreArray(user[0].fda, user[0].Ucont, &ucont);
	DAVecRestoreArray(user[0].fda, user[0].Ucat, &ucat);
	DAVecRestoreArray(user[0].fda, user[0].lUcont, &lucont);
	DAVecRestoreArray(user[0].fda, user[0].lUcat, &lucat);
	DAVecRestoreArray(user[0].fda, user[0].lCsi, &csi);
	DAVecRestoreArray(user[0].fda, user[0].lEta, &eta);
	DAVecRestoreArray(user[0].fda, user[0].lZet, &zet);
	
	outflow_scale=1;
	FormBCS(&(user[bi]),&fsi[0]);
		
	// haha
	
	if (immersed && ti>0 )
	  ibm_interpolation_advanced(&user[0]);
	
	//

	if (block_number>1) Block_Interface_U(user);

	PetscGetTime(&te);
	cput=te-ts;
	if (!rank) {
		FILE *f;
		char filen[80];
  		sprintf(filen, "%s/Converge_dU", path);
		f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d(momentum) %.2e(s) %le\n", ti, cput, norm);
		fclose(f);
	}
	VecDestroy(U);
	VecDestroy(user[bi].Rhs);
	
	
	
	for (bi=0; bi<block_number; bi++) {
		DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
		DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
		VecWAXPY(user[bi].DUold, -1., user[bi].Ucont_o, user[bi].Ucont);
		/*
		
		SNESDestroy(user[bi].snes);
		*/
		VecDestroy(user[bi].rhs);
		MatDestroy(user[bi].J);
	}
	
	double max_norm;
	VecMax(user->Ucat, &i, &max_norm);
	PetscPrintf(PETSC_COMM_WORLD, "\n*** Max Ucat = %e \n", max_norm);
	
	
	Contra2Cart(user);
	return(0);
}

PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, 
		     FSInfo *fsi)
{

  PetscReal alfa[4],ts,te,cput;
  DALocalInfo	info = user->info;

  PetscInt	bi, rank;
  PetscInt      pseudot;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


    alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  PetscReal normdU=10.,normdU0,normdU1=1.,reldU=1.,normF, normF0, relF=1.,normF1=1.; 
  PetscGetTime(&ts);
  
  pseudot=0;

  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));

    FormBCS(&(user[bi]),&fsi[0]);
    if (immersed) 
    /*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
      ibm_interpolation_advanced(&user[bi]);
    }
    
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    //VecCopy(user[bi].Ucont, user[bi].dUcont);
    VecCopy(user[bi].Ucont, user[bi].pUcont);


    // Calc RHS
    CalcRHS(&user[bi], 1);
    
    VecNorm(user[bi].Rhs, NORM_INFINITY, &normF0);
    normF1=normF0;
  }

  PetscScalar lambda=.2;
    
  //  if (ti>0) {

    // pseudo time iteration
    //for(pseudot=0; pseudot<1;pseudot++) {
    while ( ( (normF>imp_atol && relF>imp_rtol && normdU>imp_stol) || pseudot<1) && pseudot<imp_MAX_IT) {
      pseudot++;
      for (bi=0; bi<block_number; bi++) {
	for (int istage=0; istage<4; istage++) {
	  
	  // Advanced in time using RK scheme
	  VecWAXPY(user[bi].Ucont, lambda*alfa[istage] * user[bi].dt * user[bi].st , user[bi].Rhs, user[bi].pUcont);
	  
	  DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	  InflowFlux(&(user[bi]));
	  OutflowFlux(&(user[bi]));

	  PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum);

	  FormBCS(&(user[bi]),&fsi[0]);

	  // Calc the RHS
	  CalcRHS(&user[bi], 1);

	}//istage

	if (immersed) 
	  /*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
	    ibm_interpolation_advanced(&user[bi]);
	  }	

	VecWAXPY(user[bi].dUcont, -1.  , user[bi].Ucont, user[bi].pUcont);  

	VecNorm(user[bi].Rhs, NORM_INFINITY, &normF);
	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU0=normdU;
	if (pseudot>1) reldU=normdU/normdU0;
	if (pseudot==1) normF0=normF;
	if (pseudot>1) relF=normF/normF0;

	if (normF>normF1 && lambda>0.01) {
	  lambda=0.5*lambda;
	  pseudot--;
	  VecCopy(user[bi].pUcont, user[bi].Ucont);
	  
	  DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);	  
	}else {
	  VecCopy(user[bi].Ucont, user[bi].pUcont);
	  normF1=normF;
	  normdU1=normdU;
	}

	PetscGetTime(&te);
	cput=te-ts;
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le %le %le\n", pseudot, normdU, reldU,normF,relF, cput);
	if (!rank) {
	  FILE *f;
	  char filen[80];
	  //sprintf(filen, "Converge_dU%1.1d",dir);
	  sprintf(filen, "Converge_dU");
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU, reldU, normF,relF, cput);
	  fclose(f);
	}
	
      }
      if (block_number>1) {
	Block_Interface_U(user);
      }    
      
    }
    // }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(user[bi].Rhs);
    VecWAXPY(user[bi].DUold, -1., user[bi].Ucont_o, user[bi].Ucont);
    VecDestroy(user[bi].dUcont);
    VecDestroy(user[bi].pUcont);
  }

  return(0);
}


PetscErrorCode ImplicitSmoother(UserCtx *user, IBMNodes *ibm, FSInfo *fsi,PetscInt ItrNo)
{
  DA da = user->da, fda = user->fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];alfa[0] = 1.;// 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

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

  Vec           B;
  Vec           Ucont_i, RB;
  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ***nvert;
  //Vec    dUcont, pUcont;
  PetscInt	bi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp = user[0].ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;


  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));

    FormBCS(&(user[bi]),&fsi[0]);

    //ibm_interpolation_advanced2(&user[bi], ibm);
    /*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
      ibm_interpolation_advanced(&user[bi]);
    }

    //    ibm_interpolation_advanced(&user[bi], ibm, fsi);

/*     if (!CGSolver) { */
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "CGS!\n" ); */
/*     } */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  }

  //if (ti>0) {

  PetscInt pseudot;
  PetscReal normdU=10.,normdU1,reldU=1.,normdT;
  // pseudo time iteration
  //for(pseudot=0; pseudot<5;pseudot++) {
  pseudot=0;
  for (bi=0; bi<block_number; bi++) {    
    while ((normdU>1e-7 && reldU>1e-3) && pseudot<ItrNo) {
      pseudot++;
      GetPsuedoTime(&(user[bi]));    
      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	VecSet(user[bi].Rhs,0);
	Formfunction_2(&user[bi], user[bi].Rhs, 1.0);
	//FormFunction1(user[bi].Ucont, user[bi].Rhs, &(user[bi]));

	// Calculate du/dt
	VecCopy(user[bi].Ucont, user[bi].dUcont);
	VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);

	//VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 
	//VecWAXPY(B, alfa[istage]* user[bi].dt * user[bi].st  , user[bi].Rhs, user[bi].Ucont_o);
	// Calc B=RHS
	VecCopy(user[bi].Rhs, B);
	//VecScale(B,alfa[istage]*user[bi].dt* user[bi].st);
	VecSet(user[bi].dUcont, 0.);
	//PetscPrintf(PETSC_COMM_WORLD, "RHS done!\n" );
	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){
	  DAVecGetArray(fda, B, &b);
	  DAVecGetArray(da, RB, &rb);
	  DAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;			  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;	
		  }
		}
	      }
	    }
	  }

	  DAVecRestoreArray(fda, B, &b);
	  DAVecRestoreArray(da, RB, &rb);
	  DAVecRestoreArray(da, user[bi].lNvert, &nvert);

	  // Set the LHS
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
	  KSPSetType(ksp, KSPBCGS);

	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  PetscBarrier(PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);
	  
	  // Restore
	  DAVecGetArray(da, Ucont_i, &ucont_i);
	  //DAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DAVecGetArray(fda, user[bi].dUcont, &ucont);
	  
	  //PetscPrintf(PETSC_COMM_SELF, "du.i=dui!\n");
	  //PetscBarrier(PETSC_NULL);
	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    
/* 		  if ((j==97||j==98)  && (k==33 || k==32 || k==31 ||k==30 || k==29 || k==10) ){ */
/* 		    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
/* 		  } */
		  
		}
	      }
	    }
	  }
	  //PetscPrintf(PETSC_COMM_SELF, "End Restore!\n");
	  //DAVecRestoreArray(fda, B, &b);
	  
	  DAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  //DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  //DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  
	  //VecAssemblyBegin(user->Ucont);
	  //VecAssemblyEnd(user->Ucont);
	  
	  //MatNullSpaceDestroy(nullsp);
	  PetscTruth temp;
	  MatValid(user[bi].A, &temp);
	  if (temp) {
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(user[bi].A);
	    //PetscPrintf(PETSC_COMM_WORLD, "Destroy A!\n");
	  }
	  
	}
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	PetscBarrier(PETSC_NULL);
	DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU1=normdU;
	if (pseudot>1) reldU=normdU/normdU1;
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le\n", pseudot, normdU, reldU, normdT);

/* 	if (!rank) { */
/* 	  FILE *f; */
/* 	  char filen[80]; */
/* 	  sprintf(filen, "Converge_dU"); */
/* 	  f = fopen(filen, "a"); */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le\n",ti, pseudot, normdU, reldU, normdT); */
/* 	  fclose(f); */
/* 	} */

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
/* 	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum); */
	
	
	//DAVecRestoreArray(fda, Ucont_o, &ucont);
	//VecCopy(Ucont_o, user->Ucont);
	/*   DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	/*   DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	
	// Apply BC
	       
	FormBCS(&(user[bi]),&fsi[0]);
      } // istage
            
      if (immersed && ti>0) {
	/*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
	  ibm_interpolation_advanced(&user[bi]);
	}

	//ibm_interpolation_advanced(&user[bi], ibm, fsi);
	//ibm_interpolation_advanced2(&user[bi], ibm);
/* 	ibm_interpolation(ibminfo, user, ibm); */
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      
    } //bi
  } //psuedo time
  //  } // ti>0
  
  // Destroy
  for (bi=0; bi<block_number; bi++) {
    PetscTruth temp;
    MatValid(user[bi].A, &temp);
    if (temp) {
      user[bi].assignedA = PETSC_FALSE;
      MatDestroy(user[bi].A);
    }
    
    KSPDestroy(user[bi].ksp);
   
    VecDestroy(user[bi].Rhs);

/*     if (!CGSolver) */
/*       VecDestroy(user[bi].Ucont_o); */

    VecDestroy(B);
    VecDestroy(Ucont_i);
    VecDestroy(RB);
    VecDestroy(user[bi].dUcont);
    VecDestroy(user[bi].pUcont);
  }

  return(0);
}



PetscReal LinearInterpolation(PetscReal host[2], PetscReal p, PetscReal val[2])
{
  if (fabs(host[0]-host[1])>1e-9) {

  return(val[0]*(p-host[1])/(host[0]-host[1])
       + val[1]*(p-host[0])/(host[1]-host[0]));
  } else {
    PetscPrintf(PETSC_COMM_SELF, "Linear intrp Failed!!!!!!!!!!!!!\n"); 
  }
  return(0.);
}

  
PetscReal BilinearInterpolation(Cpt2D host[4], Cpt2D p, PetscReal val[4])
{
  /* Needs special ordering of host
             vkji 
     host[0]=v?00
     host[1]=v?10
     host[2]=v?01
     host[3]=v?11

  */
  PetscReal  v,x23,x01,v01,v23;

  // bilinear interpolation
  // in y
  if (fabs(host[0].y-host[1].y)>1e-9) {
  v01=val[0]*(p.y-host[1].y)/(host[0].y-host[1].y)
    + val[1]*(p.y-host[0].y)/(host[1].y-host[0].y);

  x01=host[0].x*(p.y-host[1].y)/(host[0].y-host[1].y)
    + host[1].x*(p.y-host[0].y)/(host[1].y-host[0].y);

  v23=val[2]*(p.y-host[3].y)/(host[2].y-host[3].y)
    + val[3]*(p.y-host[2].y)/(host[3].y-host[2].y);

  x23=host[2].x*(p.y-host[3].y)/(host[2].y-host[3].y)
    + host[3].x*(p.y-host[2].y)/(host[3].y-host[2].y);
  
  // in x
  if (fabs(x01-x23)>1e-9) {
  v = v01*(p.x-x23)/(x01-x23) +
      v23*(p.x-x01)/(x23-x01);
  return(v);
  } else {
   PetscPrintf(PETSC_COMM_SELF, "Bilinear intrp Failed x!!!!!!!!!!!!!\n"); 
  }

  } else {
    PetscPrintf(PETSC_COMM_SELF, "Bilinear intrp Failed y!!!!!!!!!!!!!\n"); 
  }
  return(0.);
}

PetscReal TrilinearInterpolation(Cmpnts host[8], Cmpnts p, PetscReal val[8])
{
  /* Needs special ordering of host
             vkji 
     host[0]=v000
     host[1]=v100
     host[2]=v010
     host[3]=v110
     host[4]=v001
     host[5]=v101
     host[6]=v011
     host[7]=v111
  */

  Cpt2D      bih[4], bip ;
  PetscReal  v, bval[4] ;

  bip.x=p.x;
  bip.y=p.y;
  // in z
  if (fabs(host[0].z-host[1].z)>1e-9) {
    bval[0]=val[0]*(p.z-host[1].z)/(host[0].z-host[1].z)
          + val[1]*(p.z-host[0].z)/(host[1].z-host[0].z);    
    bih[0].x=host[0].x*(p.z-host[1].z)/(host[0].z-host[1].z)
           + host[1].x*(p.z-host[0].z)/(host[1].z-host[0].z);
    bih[0].y=host[0].y*(p.z-host[1].z)/(host[0].z-host[1].z)
           + host[1].y*(p.z-host[0].z)/(host[1].z-host[0].z);
    
    bval[1]=val[2]*(p.z-host[3].z)/(host[2].z-host[3].z)
          + val[3]*(p.z-host[2].z)/(host[3].z-host[2].z);    
    bih[1].x=host[2].x*(p.z-host[3].z)/(host[2].z-host[3].z)
           + host[3].x*(p.z-host[2].z)/(host[3].z-host[2].z);
    bih[1].y=host[2].y*(p.z-host[3].z)/(host[2].z-host[3].z)
           + host[3].y*(p.z-host[2].z)/(host[3].z-host[2].z);

    bval[2]=val[4]*(p.z-host[5].z)/(host[4].z-host[5].z)
          + val[5]*(p.z-host[4].z)/(host[5].z-host[4].z);    
    bih[2].x=host[4].x*(p.z-host[5].z)/(host[4].z-host[5].z)
           + host[5].x*(p.z-host[4].z)/(host[5].z-host[4].z);
    bih[2].y=host[4].y*(p.z-host[5].z)/(host[4].z-host[5].z)
           + host[5].y*(p.z-host[4].z)/(host[5].z-host[4].z);

    bval[3]=val[6]*(p.z-host[7].z)/(host[6].z-host[7].z)
          + val[7]*(p.z-host[6].z)/(host[7].z-host[6].z);    
    bih[3].x=host[6].x*(p.z-host[7].z)/(host[6].z-host[7].z)
          + host[7].x*(p.z-host[6].z)/(host[7].z-host[6].z);
    bih[3].y=host[6].y*(p.z-host[7].z)/(host[6].z-host[7].z)
          + host[7].y*(p.z-host[6].z)/(host[7].z-host[6].z);
    
    // in y,x
    v = BilinearInterpolation(bih,bip,bval);
    return(v);
  } else {
    PetscPrintf(PETSC_COMM_SELF, "Trilinear intrp Failed!!!!!!!!!!!!!\n"); 
  }
 return(0.);
}
 
PetscErrorCode MgRestrictionP(UserCtx *user,PetscInt flevel)
{ 
  /* 6/28/06 Iman
     Restrics Pressure from the fine grid to coarse grid
     flevel   finer grid level
     clevel   coarser grid level
  */

  DA	da = user->da, fda = user->fda;

  UserCtx *user_f = user->user_f;
  DA	da_f =  user_f->da, fda_f = user_f->fda;  	
  DALocalInfo	info;
  DAGetLocalInfo(da, &info);

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

  PetscReal	***p_c, ***p_f, pRhost[8];
  Cmpnts        ***coor_c, ***coor_f, Rhost[8], R;

  PetscInt      i,j,k;

  // Get the vectors
  DAVecGetArray(fda, user->lCent,&coor_f);
  DAVecGetArray(fda_f, user_f->lCent,&coor_c);
  DAVecGetArray(da, user->lP, &p_f);
  DAVecGetArray(da_f, user_f->lP, &p_c);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	R=coor_c[k][j][i];

      /*  trilinear interpolation
	Vxyz = 	V000 (1 - x) (1 - y) (1 - z) +
	V100 x (1 - y) (1 - z) +
	V010 (1 - x) y (1 - z) +
	V001 (1 - x) (1 - y) z +
	V101 x (1 - y) z +
	V011 (1 - x) y z +
	V110 x y (1 - z) +
	V111 x y z
      */ 
       	
	Rhost[0]=coor_f[2*k-1][2*j-1][2*i-1];
	Rhost[1]=coor_f[2*k  ][2*j-1][2*i-1];
	Rhost[2]=coor_f[2*k-1][2*j  ][2*i-1];
	Rhost[3]=coor_f[2*k  ][2*j  ][2*i-1];
	Rhost[4]=coor_f[2*k-1][2*j-1][2*i  ];
	Rhost[5]=coor_f[2*k  ][2*j-1][2*i  ];
	Rhost[6]=coor_f[2*k-1][2*j  ][2*i  ];
	Rhost[7]=coor_f[2*k  ][2*j  ][2*i  ];

	pRhost[0]=p_f[2*k-1][2*j-1][2*i-1];
	pRhost[1]=p_f[2*k  ][2*j-1][2*i-1];
	pRhost[2]=p_f[2*k-1][2*j  ][2*i-1];
	pRhost[3]=p_f[2*k  ][2*j  ][2*i-1];
	pRhost[4]=p_f[2*k-1][2*j-1][2*i  ];
	pRhost[5]=p_f[2*k  ][2*j-1][2*i  ];
	pRhost[6]=p_f[2*k-1][2*j  ][2*i  ];
	pRhost[7]=p_f[2*k  ][2*j  ][2*i  ];

	p_c[k][j][i]= TrilinearInterpolation(Rhost,R,pRhost);

/* 	// Calc Coeff using lagrange formula */
/* 	for (n=0;n<8;n++) { */
/* 	  Rcoeff[n]=1.; */
/* 	  for (l=0;l<8;l++) { */
/* 	    if (l!=n) { */
/* 	      if (fabs(Rhost[n].x-Rhost[l].x)>1e-8)  */
/* 		Rcoeff[n] *=(R.x-Rhost[l].x)/(Rhost[n].x-Rhost[l].x); */
/* 	      if (fabs(Rhost[n].y-Rhost[l].y)>1e-8)  */
/* 		Rcoeff[n] *=(R.y-Rhost[l].y)/(Rhost[n].y-Rhost[l].y); */
/* 	      if (fabs(Rhost[n].z-Rhost[l].z)>1e-8)  */
/* 		Rcoeff[n] *=(R.z-Rhost[l].z)/(Rhost[n].z-Rhost[l].z); */
/* 	    }//if */
/* 	  }//l */
/* 	}//n */	     
/* 	p_c[k][j][i]=0.; */
/* 	for (n=0;n<8;n++) { */
/* 	  Rcoeff[n]=1/8.; */
/* 	  p_c[k][j][i]+=Rcoeff[n]*pRhost[n]; */
/* 	} */

      }
    }
  }

  // Restore vectors
  DAVecRestoreArray(fda, user->lCent,&coor_f);
  DAVecRestoreArray(fda_f, user_f->lCent,&coor_c);
  DAVecRestoreArray(da, user->lP, &p_f);
  DAVecRestoreArray(da_f, user_f->lP, &p_c);

  return(0);
}

PetscErrorCode MgGridInterpolation(PetscInt i, PetscInt j, PetscInt k,
				   PetscInt *ih, PetscInt *jh, PetscInt *kh,
				   PetscInt dir, UserCtx *user)
{
  //PetscPrintf(PETSC_COMM_WORLD, "grid intp/n");
  if (*(user->isc)) {
    *ih = i;
  }
  else if (dir==0 || i%2==0){
    *ih = i/2;
  }
  else {
    *ih = i/2 + 1;
  }

  if (*(user->jsc)) {
    *jh = j;
  }
  else if (dir==1 || j%2==0){
    *jh = j/2;
  }
  else {
    *jh = j/2 + 1;
  }

  if (*(user->ksc)) {
    *kh = k;
  }
  else if (dir==2 || k%2==0){
    *kh = k/2;
  }
  else {
    *kh = k/2 + 1;
  }

  return 0;
}

PetscErrorCode MgInterpolationdU_old(UserCtx *user)
{
/*  Input: Coarsegrid dU (user) */
/*  Output: Finegrid dU   */
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
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka, dir;
  Cmpnts ***dU, ***dU_f;

  if (*(user->isc)) ia = 1;
  else ia = 0;

  if (*(user->jsc)) ja = 1;
  else ja = 0;

  if (*(user->ksc)) ka = 1;
  else ka = 0;

  DAVecGetArray(fda, user->dUcont, &dU);
  DAVecGetArray(fda_f, user_f->dUcont, &dU_f);

  // x-dir
  dir=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);

	if (i%2==0) {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih  ].x )
		            / (2.-ka) / (2.-ja);
	} else {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih     ].x +
			   dU[kh   ][jh   ][ih+1-ia].x )
		            / (2.-ka) / (2.-ja) / 2.;
	}
      }
    }
  }

  // y-dir
  dir=1;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (j%2==0) {
	dU_f[k][j][i].y = (dU[kh   ][jh   ][ih  ].y )
		            / (2.-ka) / (2.-ia);
	} else {
	dU_f[k][j][i].y = (dU[kh   ][jh     ][ih  ].y +
			   dU[kh   ][jh+1-ja][ih  ].y )
		            / (2.-ka) / (2.-ia) / 2.;
	}
      }
    }
  }

  // z-dir
  dir=2;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (k%2==0) {
	dU_f[k][j][i].z = (dU[kh   ][jh   ][ih  ].z )
		            / (2.-ia) / (2.-ja);
	} else {
	dU_f[k][j][i].z = (dU[kh     ][jh   ][ih  ].z +
			   dU[kh+1-ka][jh   ][ih  ].z )
		            / (2.-ia) / (2.-ja) / 2.;
	}
      }
    }
  }

  DAVecRestoreArray(fda, user->dUcont, &dU);
  DAVecRestoreArray(fda_f, user_f->dUcont, &dU_f);

  return(0);
}
/*
PetscErrorCode MgInterpolationdU(UserCtx *user)
{
  DA	fda = user->fda;

  UserCtx *user_f = user->user_f;
  DA	da_f =  user_f->da, fda_f = user_f->fda;  	
  DALocalInfo	info;
  DAGetLocalInfo(da_f, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze, llze,llye,llxe;
  lxs = xs; lxe = xe; llxe=xe;
  lys = ys; lye = ye; llye=ye;
  lzs = zs; lze = ze; llze=ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if (xe==mx) llxe = xe-2;
  if (ye==my) llye = ye-2;
  if (ze==mz) llze = ze-2;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka, dir;
  Cmpnts ***dU, ***dU_f;
  Cmpnts ***dA, ***dA_f;

  if (*(user->isc)) ia = 1;
  else ia = 0;

  if (*(user->jsc)) ja = 1;
  else ja = 0;

  if (*(user->ksc)) ka = 1;
  else ka = 0;

  DAVecGetArray(fda, user->dUcont, &dU);
  DAVecGetArray(fda_f, user_f->dUcont, &dU_f);
  DAVecGetArray(fda, user->lArea, &dA);
  DAVecGetArray(fda_f, user_f->lArea, &dA_f);

  // x-dir
  dir=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<llxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);

	if (i%2==0) {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih  ].x /
			   dA[kh   ][jh   ][ih  ].x)
		         * dA_f[k][j][i].x;
	} else {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih     ].x /
			   dA[kh   ][jh   ][ih     ].x +
			   dU[kh   ][jh   ][ih+1-ia].x /
			   dA[kh   ][jh   ][ih+1-ia].x)
		            * dA_f[k][j][i].x / 2.;
	}
       
	if (fabs(dA_f[k][j][i].x)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_x is zero!!!! %d %d %d \n", i, j, k);

      }
    }
  }

  // y-dir
  dir=1;
  for (k=lzs; k<lze; k++) {
    for (j=ys; j<llye; j++) {
      for (i=lxs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (j%2==0) {
	dU_f[k][j][i].y = (dU[kh   ][jh   ][ih  ].y /
			   dA[kh   ][jh   ][ih  ].y)
		         * dA_f[k][j][i].y;
	} else {
	dU_f[k][j][i].y = (dU[kh   ][jh     ][ih  ].y /
			   dA[kh   ][jh     ][ih  ].y +
			   dU[kh   ][jh+1-ja][ih  ].y /
			   dA[kh   ][jh+1-ja][ih  ].y)
		         * dA_f[k][j][i].y   / 2.;
	}

	if (fabs(dA_f[k][j][i].y)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_y is zero!!!! %d %d %d \n", i, j, k);
      }
    }
  }

  // z-dir
  dir=2;
  for (k=zs; k<llze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (k%2==0) {
	dU_f[k][j][i].z = (dU[kh   ][jh   ][ih  ].z /
			   dA[kh   ][jh   ][ih  ].z)
		          * dA_f[k][j][i].z;
	} else {
	dU_f[k][j][i].z = (dU[kh     ][jh   ][ih  ].z /
			   dA[kh     ][jh   ][ih  ].z +
			   dU[kh+1-ka][jh   ][ih  ].z /
			   dA[kh+1-ka][jh   ][ih  ].z)
		           * dA_f[k][j][i].z / 2.;
	}

	if (fabs(dA_f[k][j][i].z)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_z is zero!!!! %d %d %d \n", i, j, k);

      }
    }
  }

  DAVecRestoreArray(fda, user->dUcont, &dU);
  DAVecRestoreArray(fda_f, user_f->dUcont, &dU_f);
  DAVecRestoreArray(fda, user->lArea, &dA);
  DAVecRestoreArray(fda_f, user_f->lArea, &dA_f);

  return(0);
}
*/
/*
PetscErrorCode MgFieldRestriction(UserCtx *user)
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
  Cmpnts ***ucont_o, ***ucont_o_f;
  Cmpnts ***ucont_rm1, ***ucont_rm1_f;

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
  DAVecGetArray(fda, user->Ucont_o, &ucont_o);
  DAVecGetArray(fda_f, user_f->lUcont_o, &ucont_o_f);
  DAVecGetArray(fda, user->Ucont_rm1, &ucont_rm1);
  DAVecGetArray(fda_f, user_f->lUcont_rm1, &ucont_rm1_f);

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

	ucont_o[k][j][i].x = (ucont_o_f[kh   ][jh   ][ih  ].x +
			      ucont_o_f[kh-ka][jh   ][ih  ].x +
			      ucont_o_f[kh   ][jh-ja][ih  ].x +
			      ucont_o_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);

	ucont_rm1[k][j][i].x = (ucont_rm1_f[kh   ][jh   ][ih  ].x +
				ucont_rm1_f[kh-ka][jh   ][ih  ].x +
				ucont_rm1_f[kh   ][jh-ja][ih  ].x +
				ucont_rm1_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);
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

	ucont_o[k][j][i].y = (ucont_o_f[kh   ][jh  ][ih   ].y +
			      ucont_o_f[kh-ka][jh  ][ih   ].y +
			      ucont_o_f[kh   ][jh  ][ih-ia].y +
			      ucont_o_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);

	ucont_rm1[k][j][i].y = (ucont_rm1_f[kh   ][jh  ][ih   ].y +
				ucont_rm1_f[kh-ka][jh  ][ih   ].y +
				ucont_rm1_f[kh   ][jh  ][ih-ia].y +
				ucont_rm1_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);
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

	ucont_o[k][j][i].z = (ucont_o_f[kh  ][jh   ][ih   ].z +
			      ucont_o_f[kh  ][jh   ][ih-ia].z +
			      ucont_o_f[kh  ][jh-ja][ih   ].z +
			      ucont_o_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);

	ucont_rm1[k][j][i].z = (ucont_rm1_f[kh  ][jh   ][ih   ].z +
				ucont_rm1_f[kh  ][jh   ][ih-ia].z +
				ucont_rm1_f[kh  ][jh-ja][ih   ].z +
				ucont_rm1_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);
      }
    }
  }

  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(fda_f, user_f->lUcont, &ucont_f);
  DAVecRestoreArray(fda, user->Ucont_o, &ucont_o);
  DAVecRestoreArray(fda_f, user_f->lUcont_o, &ucont_o_f);
  DAVecRestoreArray(fda, user->Ucont_rm1, &ucont_rm1);
  DAVecRestoreArray(fda_f, user_f->lUcont_rm1, &ucont_rm1_f);
  
  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DAGlobalToLocalBegin(fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
  DAGlobalToLocalEnd(fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);

  DAGlobalToLocalBegin(fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1);
  DAGlobalToLocalEnd(fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1);

  VecCopy(user->Ucont, user->Ucont_MG);

  PetscReal ***p, ***p_f, ***v_f, v_tot;

  DAVecGetArray(da, user->P, &p);
  DAVecGetArray(da_f, user_f->lP, &p_f);
  DAVecGetArray(da_f, user_f->lVolume, &v_f);


  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	v_tot =(      v_f[kh   ][jh   ][ih   ] +		
		      v_f[kh   ][jh   ][ih-ia] +		
		      v_f[kh   ][jh-ja][ih   ] +		
		      v_f[kh   ][jh-ja][ih-ia] +		
		      v_f[kh-ka][jh   ][ih   ] +		
		      v_f[kh-ka][jh   ][ih-ia] +		
		      v_f[kh-ka][jh-ja][ih   ] +		
		      v_f[kh-ka][jh-ja][ih-ia]    );

	p[k][j][i] = (p_f[kh   ][jh   ][ih   ] *
		      v_f[kh   ][jh   ][ih   ] +
		      p_f[kh   ][jh   ][ih-ia] *
		      v_f[kh   ][jh   ][ih-ia] +
		      p_f[kh   ][jh-ja][ih   ] *
		      v_f[kh   ][jh-ja][ih   ] +
		      p_f[kh   ][jh-ja][ih-ia] *
		      v_f[kh   ][jh-ja][ih-ia] +
		      p_f[kh-ka][jh   ][ih   ] *
		      v_f[kh-ka][jh   ][ih   ] +
		      p_f[kh-ka][jh   ][ih-ia] *
		      v_f[kh-ka][jh   ][ih-ia] +
		      p_f[kh-ka][jh-ja][ih   ] *
		      v_f[kh-ka][jh-ja][ih   ] +
		      p_f[kh-ka][jh-ja][ih-ia] *
		      v_f[kh-ka][jh-ja][ih-ia])/v_tot;

      }
    }
  }

  
  DAVecRestoreArray(da, user->P, &p);
  DAVecRestoreArray(da_f, user_f->lP, &p_f);
  DAVecRestoreArray(da_f, user_f->lVolume, &v_f);
  
  DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
		      
  return 0;
}
*/

/* ==================================================================================             */
/*      Multi-Grid Cycle  */
/*
PetscErrorCode MGCYC(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi,PetscInt level, PetscInt cycl_idx, PetscInt preItr, PetscInt postItr,PetscInt CGItr)
{

  PetscInt i;
  Vec      pUcont, dUcont;

  UserCtx *user, *user_c;

  user   = usermg->mgctx[level].user;
  user_c = usermg->mgctx[level-1].user;

  DA    fda=user->fda;

//     Presmoothing fG 
  ImplicitSmoother(user,ibm, fsi,preItr);
  PetscPrintf(PETSC_COMM_WORLD, "PRE SMOOTHING!\n");


//    CGC: 
//    1) Restrict U,p to CG
    
    //MgRestrictionP(user, level, isc, jsc, ksc);
    MgFieldRestriction(user_c);
    //MyNvertRestriction(user, user_c);
    PetscPrintf(PETSC_COMM_WORLD, "RESTRICTION!!!\n");

//      save Ucont  
    VecDuplicate(user_c->lUcont, &pUcont);   
    VecDuplicate(user_c->lUcont, &dUcont);   
    VecCopy(user_c->lUcont, pUcont);

//      2) Solve Ucont on CG 
    if (level-1==0) { // Coarsest level: Solve to machine zero!
      //ImplicitMomentumSolver(user_c, ibm, fsi);
      ImplicitSmoother(user_c,ibm, fsi, CGItr);
      PetscPrintf(PETSC_COMM_WORLD, "CG Solver!!!\n");

    } else { // not the coarsest level: Solve by MGM
      for (i=0; i<cycl_idx; i++) {
	MGCYC(usermg, ibm, fsi,level-1,cycl_idx, preItr, postItr, CGItr);
      }
    }
  

//      calc dUcont  
    VecWAXPY(dUcont, -1., pUcont, user_c->lUcont);
    VecDuplicate(user_c->lUcont, &(user_c->dUcont));
    VecCopy(dUcont, user_c->dUcont);
    VecDuplicate(user->Ucont, &(user->dUcont));


      //3) Interpolate dU from CG to finer Grid
    MgInterpolationdU(user_c);  
    // Ucont=Ucont+dU on fine Grid
    VecDuplicate(user->Ucont, &(user->pUcont));
    VecCopy(user->Ucont, user->pUcont);
    VecWAXPY(user->Ucont, +1., user->dUcont, user->pUcont);  
    
    DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

    PetscPrintf(PETSC_COMM_WORLD, "INTERPOLATION!!!\n");
    

//      Destroy 
    VecDestroy(dUcont);
    VecDestroy(pUcont);
    VecDestroy(user->pUcont);
    VecDestroy(user_c->dUcont);
    VecDestroy(user->dUcont);

//    Postsmoothing fG 
  ImplicitSmoother(user,ibm,fsi, postItr);
  PetscPrintf(PETSC_COMM_WORLD, "POST SMOOTHING!\n");

  
  return(0);
}
*/
/* ==================================================================================             */
/* ==================================================================================             */
/*
PetscErrorCode MGMomentumSolver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi)
{
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt finest_level=usermg->mglevels - 1;
  UserCtx *user, *user_l;
  user = usermg->mgctx[finest_level].user;
  
  Vec pUcont, dUcont;
  
  PetscInt bi, level;
  PetscInt pseudot;
  PetscReal normdU=10.,normdU1=1.,reldU=1.,normdT, cput,ts,te;
  pseudot=0;
  
  for (bi=0; bi<block_number; bi++) { 

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]),&fsi[0]);

    if (immersed)
	{
      ibm_interpolation_advanced(&user[bi]);
    }

    for (level=finest_level-1; level>-1; level--) {
      user_l=usermg->mgctx[level].user;
      MyFieldRestriction(user_l);
    }

    VecDuplicate(user[bi].Ucont, &pUcont);   
    VecDuplicate(user[bi].Ucont, &dUcont);

    PetscGetTime(&ts);
    
    while (( (normdU>1e-7 && reldU>1e-3) || pseudot<3) && pseudot<mg_MAX_IT) {
      pseudot++;
      VecCopy(user[bi].Ucont, pUcont);

      MGCYC(usermg, ibm, fsi, finest_level, mg_idx, mg_preItr,mg_poItr,imp_MAX_IT);
      
      VecWAXPY(dUcont, -1., pUcont, user[bi].Ucont);


      VecNorm(dUcont, NORM_INFINITY, &normdU);
      if (pseudot==1) normdU1=normdU;
      if (pseudot>1) reldU=normdU/normdU1;
      VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
      PetscGetTime(&te);
      cput=te-ts;
      PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU MG %d  %le %le %le %le\n", pseudot, normdU, reldU, normdT, cput);

      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      
      InflowFlux(&(user[bi]));
      OutflowFlux(&(user[bi]));            

      FormBCS(&(user[bi]),&fsi[0]);

      if (immersed) 
	{
	  ibm_interpolation_advanced(&user[bi]);
	}

      if (!rank) {
	FILE *f;
	char filen[80];
	sprintf(filen, "Converge_dU_MG");
	f = fopen(filen, "a");
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le\n",ti,pseudot, normdU, reldU, normdT, cput);
	fclose(f);
      }
    }

    VecDestroy(pUcont);
    VecDestroy(dUcont);
  }
  return(0);
}

*/

void initial_guess_for_snes(UserCtx *user, Vec U)
{
	DA	da = user->da, fda = user->fda;
	
	DALocalInfo	info = user->info;
	int	xs = info.xs, xe = info.xs + info.xm;
	int	ys = info.ys, ye = info.ys + info.ym;
	int	zs = info.zs, ze = info.zs + info.zm;
	int	mx = info.mx, my = info.my, mz = info.mz;

	int  i,j,k;
	PetscReal	***nvert;
	Cmpnts ***u, ***rhs, ***ucont;

	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	VecCopy(user->Ucont, U);
	
	PetscPrintf (PETSC_COMM_WORLD, "Generaing initial guess ...\n");
	VecSet(user->Rhs,0);
	Formfunction_2(user, user->Rhs, 1.0);
	
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, U, &u);
	DAVecGetArray(fda, user->Rhs, &rhs);
	DAVecGetArray(fda, user->Ucont, &ucont);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if((int)nvert[k][j][i]==0 && (int)nvert[k][j][i+1]==0) u[k][j][i].x = ucont[k][j][i].x + user->dt * rhs[k][j][i].x;
		else u[k][j][i].x = ucont[k][j][i].x;
		
		if((int)nvert[k][j][i]==0 && (int)nvert[k][j+1][i]==0) u[k][j][i].y = ucont[k][j][i].y + user->dt * rhs[k][j][i].y;
		else u[k][j][i].y = ucont[k][j][i].y;
		
		if((int)nvert[k][j][i]==0 && (int)nvert[k+1][j][i]==0) u[k][j][i].z = ucont[k][j][i].z + user->dt * rhs[k][j][i].z;
		else u[k][j][i].z = ucont[k][j][i].z;
	}
   	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, U, &u);
	DAVecRestoreArray(fda, user->Rhs, &rhs);
	DAVecRestoreArray(fda, user->Ucont, &ucont);
}
