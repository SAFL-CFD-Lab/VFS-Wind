/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
#define NEWMETRIC

PetscErrorCode FormMetrics(UserCtx *user)
{
  DA		cda;
  Cmpnts	***csi, ***eta, ***zet;
  PetscScalar	***aj;
  Vec		coords;
  Cmpnts	***coor;

  DA		da = user->da, fda = user->fda;
  Vec		Csi = user->Csi, Eta = user->Eta, Zet = user->Zet;
  Vec		Aj = user->Aj;
  Vec		ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet;
  Vec		JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet;
  Vec		KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet;
  Vec		IAj = user->IAj, JAj = user->JAj, KAj = user->KAj;

  
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  Cmpnts	***gs;
  PetscReal	***iaj, ***jaj, ***kaj;

  Vec		Cent = user->Cent; //local working array for storing cell center geometry

  Vec		Centx, Centy, Centz, lCoor;
  Cmpnts	***cent, ***centx, ***centy, ***centz;

  PetscInt	xs, ys, zs, xe, ye, ze;
  DALocalInfo	info;

  PetscInt	mx, my, mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt	i, j, k;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscErrorCode	ierr;

  PetscReal	xcp, ycp, zcp, xcm, ycm, zcm;
  DAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DAGetCoordinateDA(da, &cda);
  DAVecGetArray(fda, Csi, &csi);
  DAVecGetArray(fda, Eta, &eta);
  DAVecGetArray(fda, Zet, &zet);
  ierr = DAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DAGetGhostedCoordinates(da, &coords);
  DAVecGetArray(fda, coords, &coor);

  DAGetLocalVector(fda, &Centx);

  VecDuplicate(Centx, &Centy);
  VecDuplicate(Centx, &Centz);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


	DAVecGetArray(fda, user->Cent, &cent);
  	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++) 
	for (i=lxs; i<lxe; i++) {
		cent[k][j][i].x = 0.125 *
			(coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
			coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
			coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
		cent[k][j][i].y = 0.125 *
			(coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
			coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
			coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
		cent[k][j][i].z = 0.125 *
			(coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
			coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
			coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);
	}
	DAVecRestoreArray(fda, user->Cent, &cent);
	
	DAGlobalToLocalBegin(fda, user->Cent, INSERT_VALUES, user->lCent);
	DAGlobalToLocalEnd(fda, user->Cent, INSERT_VALUES, user->lCent);
	
	DAVecGetArray(fda, user->lCent, &cent);
	
	#if !defined (NEWMETRIC)
  /* Calculating transformation metrics in i direction */
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=xs; i<lxe; i++) {
		/* csi = de X dz */
		double dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
			      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
		double dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
			      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
		double dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
			      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);

		double dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
			      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
		double dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
			      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
		double dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
			      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);
		
		csi[k][j][i].x = dyde * dzdz - dzde * dydz;
		csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
		csi[k][j][i].z = dxde * dydz - dyde * dxdz;
		
	}
    
	// Need more work -- lg65
	/* calculating j direction metrics */
	for (k=lzs; k<lze; k++)
	for (j=ys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		/* eta = dz X de */
		double dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
			      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
		double dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
			      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
		double dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
			      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);
											 
		double dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
			      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
		double dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
			      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
		double dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
			      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);
		
		eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
		eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
		eta[k][j][i].z = dxdz * dydc - dydz * dxdc;
	}

	/* calculating k direction metrics */
	for (k=zs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
			      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
		double dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
			      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
		double dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
			      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);
									 
		double dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
			      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
		double dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
			      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
		double dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
			      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);
		
		zet[k][j][i].x = dydc * dzde - dzdc * dyde;
		zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
		zet[k][j][i].z = dxdc * dyde - dydc * dxde;
		
	}
	
	DAVecRestoreArray(fda, user->lCent, &cent);
	#endif
	
  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	if (i>0 && j>0 && k>0) {
	  double dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
			 coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  double dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
			 coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  double dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
			 coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

	  double dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x - 
			 coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
	  double dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y - 
			 coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
	  double dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z - 
			 coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

	  double dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  double dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  double dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);
	  
		aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
				dydc * (dxde * dzdz - dzde * dxdz) +
				dzdc * (dxde * dydz - dyde * dxdz);
		aj[k][j][i] = 1./aj[k][j][i];
	  
		#ifdef NEWMETRIC
		csi[k][j][i].x = dyde * dzdz - dzde * dydz;
		csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
		csi[k][j][i].z = dxde * dydz - dyde * dxdz;
		
		eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
		eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
		eta[k][j][i].z = dxdz * dydc - dydz * dxdc;
	  
		zet[k][j][i].x = dydc * dzde - dzdc * dyde;
		zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
		zet[k][j][i].z = dxdc * dyde - dydc * dxde;
		#endif
	}
      }
    }
  }

  // mirror grid outside the boundary
	if (xs==0) {
		i = xs;
		for (k=zs; k<ze; k++) 
		for (j=ys; j<ye; j++) {
			#ifdef NEWMETRIC
			csi[k][j][i] = csi[k][j][i+1];
			#endif
			eta[k][j][i] = eta[k][j][i+1];
			zet[k][j][i] = zet[k][j][i+1];
			aj[k][j][i] = aj[k][j][i+1];
		}
	}

	if (xe==mx) {
		i = xe-1;
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++) {
			#ifdef NEWMETRIC
			csi[k][j][i] = csi[k][j][i-1];
			#endif
			eta[k][j][i] = eta[k][j][i-1];
			zet[k][j][i] = zet[k][j][i-1];
			aj[k][j][i] = aj[k][j][i-1];
		}
	}
  

	if (ys==0) {
		j = ys;
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			eta[k][j][i] = eta[k][j+1][i];
			#endif
			csi[k][j][i] = csi[k][j+1][i];
			zet[k][j][i] = zet[k][j+1][i];
			aj[k][j][i] = aj[k][j+1][i];
		}
	}
  

	if (ye==my) {
		j = ye-1;
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			eta[k][j][i] = eta[k][j-1][i];
			#endif
			csi[k][j][i] = csi[k][j-1][i];
			zet[k][j][i] = zet[k][j-1][i];
			aj[k][j][i] = aj[k][j-1][i];
		}
	}
  

	if (zs==0) {
		k = zs;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			zet[k][j][i] = zet[k+1][j][i];
			#endif
			eta[k][j][i] = eta[k+1][j][i];
			csi[k][j][i] = csi[k+1][j][i];
			aj[k][j][i] = aj[k+1][j][i];
		}
	}
	

	if (ze==mz) {
		k = ze-1;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			zet[k][j][i] = zet[k-1][j][i];
			#endif
			eta[k][j][i] = eta[k-1][j][i];
			csi[k][j][i] = csi[k-1][j][i];
			aj[k][j][i] = aj[k-1][j][i];
		}
	}
	
	DAVecRestoreArray(fda, Csi, &csi);
	DAVecRestoreArray(fda, Eta, &eta);
	DAVecRestoreArray(fda, Zet, &zet);
	DAVecRestoreArray(da, Aj,  &aj);
	
	DAGlobalToLocalBegin(fda, user->Csi, INSERT_VALUES, user->lCsi);
	DAGlobalToLocalEnd(fda, user->Csi, INSERT_VALUES, user->lCsi);

	DAGlobalToLocalBegin(fda, user->Eta, INSERT_VALUES, user->lEta);
	DAGlobalToLocalEnd(fda, user->Eta, INSERT_VALUES, user->lEta);

	DAGlobalToLocalBegin(fda, user->Zet, INSERT_VALUES, user->lZet);
	DAGlobalToLocalEnd(fda, user->Zet, INSERT_VALUES, user->lZet);
  
	DAGlobalToLocalBegin(da, user->Aj, INSERT_VALUES, user->lAj);
	DAGlobalToLocalEnd(da, user->Aj, INSERT_VALUES, user->lAj);
	
		DAVecGetArray(fda, user->lCsi, &csi);
		DAVecGetArray(fda, user->lEta, &eta);
		DAVecGetArray(fda, user->lZet, &zet);
		DAVecGetArray(da, user->lAj,  &aj);
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			int flag=0, a=i, b=j, c=k;
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
			
			flag = i_flag + j_flag + k_flag;

			if(flag) {
				csi[k][j][i] = csi[c][b][a];
				eta[k][j][i] = eta[c][b][a];
				zet[k][j][i] = zet[c][b][a];
				aj[k][j][i] = aj[c][b][a];
			}
		}
		DAVecRestoreArray(fda, user->lCsi, &csi);
		DAVecRestoreArray(fda, user->lEta, &eta);
		DAVecRestoreArray(fda, user->lZet, &zet);
		DAVecRestoreArray(da, user->lAj,  &aj);
	
	Cmpnts	***lcsi, ***leta, ***lzet;
	PetscScalar	***laj;
	
	
	DAVecGetArray(fda, user->lCsi, &lcsi);
	DAVecGetArray(fda, user->lEta, &leta);
	DAVecGetArray(fda, user->lZet, &lzet);
	DAVecGetArray(da, user->lAj,  &laj);
	
  DAVecGetArray(fda, ICsi, &icsi);
  DAVecGetArray(fda, IEta, &ieta);
  DAVecGetArray(fda, IZet, &izet);
  DAVecGetArray(da, IAj,  &iaj);

	
	
	DAVecGetArray(fda, user->GridSpace, &gs);	
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	xcp = 0.25 *
	  (coor[k][j][i].x + coor[k][j-1][i].x +
	   coor[k-1][j-1][i].x + coor[k-1][j][i].x);
	ycp = 0.25 *
	  (coor[k][j][i].y + coor[k][j-1][i].y +
	   coor[k-1][j-1][i].y + coor[k-1][j][i].y);
	zcp = 0.25 *
	  (coor[k][j][i].z + coor[k][j-1][i].z +
	   coor[k-1][j-1][i].z + coor[k-1][j][i].z);

	xcm = 0.25 *
	  (coor[k][j][i-1].x + coor[k][j-1][i-1].x +
	   coor[k-1][j-1][i-1].x + coor[k-1][j][i-1].x);
	ycm = 0.25 *
	  (coor[k][j][i-1].y + coor[k][j-1][i-1].y +
	   coor[k-1][j-1][i-1].y + coor[k-1][j][i-1].y);
	zcm = 0.25 *
	  (coor[k][j][i-1].z + coor[k][j-1][i-1].z +
	   coor[k-1][j-1][i-1].z + coor[k-1][j][i-1].z);

	gs[k][j][i].x = sqrt((xcp-xcm) * (xcp-xcm) +
			     (ycp-ycm) * (ycp-ycm) +
			     (zcp-zcm) * (zcp-zcm));

	xcp = 0.25 *
	  (coor[k][j][i].x + coor[k][j][i-1].x +
	   coor[k-1][j][i].x + coor[k-1][j][i-1].x);
	ycp = 0.25 *
	  (coor[k][j][i].y + coor[k][j][i-1].y +
	   coor[k-1][j][i].y + coor[k-1][j][i-1].y);
	zcp = 0.25 *
	  (coor[k][j][i].z + coor[k][j][i-1].z +
	   coor[k-1][j][i].z + coor[k-1][j][i-1].z);

	xcm = 0.25 *
	  (coor[k][j-1][i].x + coor[k][j-1][i-1].x +
	   coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x);
	ycm = 0.25 *
	  (coor[k][j-1][i].y + coor[k][j-1][i-1].y +
	   coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y);
	zcm = 0.25 *
	  (coor[k][j-1][i].z + coor[k][j-1][i-1].z +
	   coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z);

	gs[k][j][i].y = sqrt((xcp-xcm) * (xcp-xcm) +
			     (ycp-ycm) * (ycp-ycm) +
			     (zcp-zcm) * (zcp-zcm));

	xcp = 0.25 *
	  (coor[k][j][i].x + coor[k][j][i-1].x +
	   coor[k][j-1][i].x + coor[k][j-1][i-1].x);
	ycp = 0.25 *
	  (coor[k][j][i].y + coor[k][j][i-1].y +
	   coor[k][j-1][i].y + coor[k][j-1][i-1].y);
	zcp = 0.25 *
	  (coor[k][j][i].z + coor[k][j][i-1].z +
	   coor[k][j-1][i].z + coor[k][j-1][i-1].z);

	xcm = 0.25 *
	  (coor[k-1][j][i].x + coor[k-1][j][i-1].x +
	   coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x);
	ycm = 0.25 *
	  (coor[k-1][j][i].y + coor[k-1][j][i-1].y +
	   coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y);
	zcm = 0.25 *
	  (coor[k-1][j][i].z + coor[k-1][j][i-1].z +
	   coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z);

	gs[k][j][i].z = sqrt((xcp-xcm) * (xcp-xcm) +
			     (ycp-ycm) * (ycp-ycm) +
			     (zcp-zcm) * (zcp-zcm));

      }
    }
  }

  DAVecRestoreArray(fda, user->GridSpace, &gs);

  DAVecGetArray(fda, Centx, &centx);
  for(k=gzs+1; k<gze; k++) {
    for (j=gys+1; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	centx[k][j][i].x = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x +
			    coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;
	centx[k][j][i].y = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y +
			    coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;
	centx[k][j][i].z = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z +
			    coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25;

      }
    }
  }


  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {  
	PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
	
	if (i==0) {
	  dxdc = centx[k][j][i+1].x - centx[k][j][i].x;
	  dydc = centx[k][j][i+1].y - centx[k][j][i].y;
	  dzdc = centx[k][j][i+1].z - centx[k][j][i].z;
	}
	else if (i==mx-2) {
	  dxdc = centx[k][j][i].x - centx[k][j][i-1].x;
	  dydc = centx[k][j][i].y - centx[k][j][i-1].y;
	  dzdc = centx[k][j][i].z - centx[k][j][i-1].z;
	}
	else {
	  dxdc = (centx[k][j][i+1].x - centx[k][j][i-1].x) * 0.5;
	  dydc = (centx[k][j][i+1].y - centx[k][j][i-1].y) * 0.5;
	  dzdc = (centx[k][j][i+1].z - centx[k][j][i-1].z) * 0.5;
	}

	
	if (j==1) {
	  dxde = centx[k][j+1][i].x - centx[k][j][i].x;
	  dyde = centx[k][j+1][i].y - centx[k][j][i].y;
	  dzde = centx[k][j+1][i].z - centx[k][j][i].z;
	}
	else if (j==my-2) {
	  dxde = centx[k][j][i].x - centx[k][j-1][i].x;
	  dyde = centx[k][j][i].y - centx[k][j-1][i].y;
	  dzde = centx[k][j][i].z - centx[k][j-1][i].z;
	}
	else {
	  dxde = (centx[k][j+1][i].x - centx[k][j-1][i].x) * 0.5;
	  dyde = (centx[k][j+1][i].y - centx[k][j-1][i].y) * 0.5;
	  dzde = (centx[k][j+1][i].z - centx[k][j-1][i].z) * 0.5;
	}
	
	if (k==1) {
	  dxdz = (centx[k+1][j][i].x - centx[k][j][i].x);
	  dydz = (centx[k+1][j][i].y - centx[k][j][i].y);
	  dzdz = (centx[k+1][j][i].z - centx[k][j][i].z);
	}
	else if (k==mz-2) {
	  dxdz = (centx[k][j][i].x - centx[k-1][j][i].x);
	  dydz = (centx[k][j][i].y - centx[k-1][j][i].y);
	  dzdz = (centx[k][j][i].z - centx[k-1][j][i].z);
	}
	else {
	  dxdz = (centx[k+1][j][i].x - centx[k-1][j][i].x) * 0.5;
	  dydz = (centx[k+1][j][i].y - centx[k-1][j][i].y) * 0.5;
	  dzdz = (centx[k+1][j][i].z - centx[k-1][j][i].z) * 0.5;
	}
	
	icsi[k][j][i].x = dyde * dzdz - dzde * dydz;
	icsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	icsi[k][j][i].z = dxde * dydz - dyde * dxdz;

	ieta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	ieta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	ieta[k][j][i].z = dxdz * dydc - dydz * dxdc;

	izet[k][j][i].x = dydc * dzde - dzdc * dyde;
	izet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	izet[k][j][i].z = dxdc * dyde - dydc * dxde;

	iaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	  dydc * (dxde * dzdz - dzde * dxdz) +
	  dzdc * (dxde * dydz - dyde * dxdz);
	iaj[k][j][i] = 1./iaj[k][j][i];

	#ifdef NEWMETRIC
	AxByC ( 0.5, lcsi[k][j][i], 0.5, lcsi[k][j][i+1], &icsi[k][j][i]);
	AxByC ( 0.5, leta[k][j][i], 0.5, leta[k][j][i+1], &ieta[k][j][i]);
	AxByC ( 0.5, lzet[k][j][i], 0.5, lzet[k][j][i+1], &izet[k][j][i]);
	iaj[k][j][i] = 2. / ( 1./laj[k][j][i] + 1./laj[k][j][i+1] );
	#endif

      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "test\n");

  DAVecRestoreArray(fda, ICsi, &icsi);
  DAVecRestoreArray(fda, IEta, &ieta);
  DAVecRestoreArray(fda, IZet, &izet);
  DAVecRestoreArray(da, IAj,  &iaj);

  // j direction
  DAVecGetArray(fda, JCsi, &jcsi);
  DAVecGetArray(fda, JEta, &jeta);
  DAVecGetArray(fda, JZet, &jzet);
  DAVecGetArray(da, JAj,  &jaj);

  DAVecGetArray(fda, Centy, &centy);
  for(k=gzs+1; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs+1; i<gxe; i++) {
	centy[k][j][i].x = (coor[k  ][j][i  ].x + coor[k-1][j][i  ].x +
			    coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;
	centy[k][j][i].y = (coor[k  ][j][i  ].y + coor[k-1][j][i  ].y +
			    coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;
	centy[k][j][i].z = (coor[k  ][j][i  ].z + coor[k-1][j][i  ].z +
			    coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25;
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
	if (i==1) {
	  dxdc = centy[k][j][i+1].x - centy[k][j][i].x;
	  dydc = centy[k][j][i+1].y - centy[k][j][i].y;
	  dzdc = centy[k][j][i+1].z - centy[k][j][i].z;
	}
	else if (i==mx-2) {
	  dxdc = centy[k][j][i].x - centy[k][j][i-1].x;
	  dydc = centy[k][j][i].y - centy[k][j][i-1].y;
	  dzdc = centy[k][j][i].z - centy[k][j][i-1].z;
	}
	else {
	  dxdc = (centy[k][j][i+1].x - centy[k][j][i-1].x) * 0.5;
	  dydc = (centy[k][j][i+1].y - centy[k][j][i-1].y) * 0.5;
	  dzdc = (centy[k][j][i+1].z - centy[k][j][i-1].z) * 0.5;
	}

	if (j==0) {
	  dxde = centy[k][j+1][i].x - centy[k][j][i].x;
	  dyde = centy[k][j+1][i].y - centy[k][j][i].y;
	  dzde = centy[k][j+1][i].z - centy[k][j][i].z;
	}
	else if (j==my-2) {
	  dxde = centy[k][j][i].x - centy[k][j-1][i].x;
	  dyde = centy[k][j][i].y - centy[k][j-1][i].y;
	  dzde = centy[k][j][i].z - centy[k][j-1][i].z;
	}
	else {
	  dxde = (centy[k][j+1][i].x - centy[k][j-1][i].x) * 0.5;
	  dyde = (centy[k][j+1][i].y - centy[k][j-1][i].y) * 0.5;
	  dzde = (centy[k][j+1][i].z - centy[k][j-1][i].z) * 0.5;
	}

	if (k==1) {
	  dxdz = (centy[k+1][j][i].x - centy[k][j][i].x);
	  dydz = (centy[k+1][j][i].y - centy[k][j][i].y);
	  dzdz = (centy[k+1][j][i].z - centy[k][j][i].z);
	}
	else if (k==mz-2) {
	  dxdz = (centy[k][j][i].x - centy[k-1][j][i].x);
	  dydz = (centy[k][j][i].y - centy[k-1][j][i].y);
	  dzdz = (centy[k][j][i].z - centy[k-1][j][i].z);
	}
	else {
	  dxdz = (centy[k+1][j][i].x - centy[k-1][j][i].x) * 0.5;
	  dydz = (centy[k+1][j][i].y - centy[k-1][j][i].y) * 0.5;
	  dzdz = (centy[k+1][j][i].z - centy[k-1][j][i].z) * 0.5;
	}

	jcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
	jcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	jcsi[k][j][i].z = dxde * dydz - dyde * dxdz;
	
	jeta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	jeta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	jeta[k][j][i].z = dxdz * dydc - dydz * dxdc;

	jzet[k][j][i].x = dydc * dzde - dzdc * dyde;
	jzet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	jzet[k][j][i].z = dxdc * dyde - dydc * dxde;

	
	jaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	  dydc * (dxde * dzdz - dzde * dxdz) +
	  dzdc * (dxde * dydz - dyde * dxdz);
	jaj[k][j][i] = 1./jaj[k][j][i];
	
	#ifdef NEWMETRIC
	AxByC ( 0.5, lcsi[k][j][i], 0.5, lcsi[k][j+1][i], &jcsi[k][j][i]);
	AxByC ( 0.5, leta[k][j][i], 0.5, leta[k][j+1][i], &jeta[k][j][i]);
	AxByC ( 0.5, lzet[k][j][i], 0.5, lzet[k][j+1][i], &jzet[k][j][i]);
	jaj[k][j][i] = 2. / ( 1./laj[k][j][i] + 1./laj[k][j+1][i] );
	#endif
      }
    }
  }

  DAVecRestoreArray(fda, JCsi, &jcsi);
  DAVecRestoreArray(fda, JEta, &jeta);
  DAVecRestoreArray(fda, JZet, &jzet);
  DAVecRestoreArray(da, JAj,  &jaj);

  // k direction
  DAVecGetArray(fda, KCsi, &kcsi);
  DAVecGetArray(fda, KEta, &keta);
  DAVecGetArray(fda, KZet, &kzet);
  DAVecGetArray(da, KAj,  &kaj);

  DAVecGetArray(fda, Centz, &centz);
  for(k=gzs; k<gze; k++) {
    for (j=gys+1; j<gye; j++) {
      for (i=gxs+1; i<gxe; i++) {
	centz[k][j][i].x = (coor[k  ][j][i  ].x + coor[k][j-1][i  ].x +
			    coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
	centz[k][j][i].y = (coor[k  ][j][i  ].y + coor[k][j-1][i  ].y +
			    coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
	centz[k][j][i].z = (coor[k  ][j][i  ].z + coor[k][j-1][i  ].z +
			    coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;
      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {  
	PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
	
	if (i==1) {
	  dxdc = centz[k][j][i+1].x - centz[k][j][i].x;
	  dydc = centz[k][j][i+1].y - centz[k][j][i].y;
	  dzdc = centz[k][j][i+1].z - centz[k][j][i].z;
	}
	else if (i==mx-2) {
	  dxdc = centz[k][j][i].x - centz[k][j][i-1].x;
	  dydc = centz[k][j][i].y - centz[k][j][i-1].y;
	  dzdc = centz[k][j][i].z - centz[k][j][i-1].z;
	}
	else {
	  dxdc = (centz[k][j][i+1].x - centz[k][j][i-1].x) * 0.5;
	  dydc = (centz[k][j][i+1].y - centz[k][j][i-1].y) * 0.5;
	  dzdc = (centz[k][j][i+1].z - centz[k][j][i-1].z) * 0.5;
	}

	
	if (j==1) {
	  dxde = centz[k][j+1][i].x - centz[k][j][i].x;
	  dyde = centz[k][j+1][i].y - centz[k][j][i].y;
	  dzde = centz[k][j+1][i].z - centz[k][j][i].z;
	}
	else if (j==my-2) {
	  dxde = centz[k][j][i].x - centz[k][j-1][i].x;
	  dyde = centz[k][j][i].y - centz[k][j-1][i].y;
	  dzde = centz[k][j][i].z - centz[k][j-1][i].z;
	}
	else {
	  dxde = (centz[k][j+1][i].x - centz[k][j-1][i].x) * 0.5;
	  dyde = (centz[k][j+1][i].y - centz[k][j-1][i].y) * 0.5;
	  dzde = (centz[k][j+1][i].z - centz[k][j-1][i].z) * 0.5;
	}

	
	if (k==0) {
	  dxdz = (centz[k+1][j][i].x - centz[k][j][i].x);
	  dydz = (centz[k+1][j][i].y - centz[k][j][i].y);
	  dzdz = (centz[k+1][j][i].z - centz[k][j][i].z);
	}
	else if (k==mz-2) {
	  dxdz = (centz[k][j][i].x - centz[k-1][j][i].x);
	  dydz = (centz[k][j][i].y - centz[k-1][j][i].y);
	  dzdz = (centz[k][j][i].z - centz[k-1][j][i].z);
	}
	else {
	  dxdz = (centz[k+1][j][i].x - centz[k-1][j][i].x) * 0.5;
	  dydz = (centz[k+1][j][i].y - centz[k-1][j][i].y) * 0.5;
	  dzdz = (centz[k+1][j][i].z - centz[k-1][j][i].z) * 0.5;
	}

	kcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
	kcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	kcsi[k][j][i].z = dxde * dydz - dyde * dxdz;

	keta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	keta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	keta[k][j][i].z = dxdz * dydc - dydz * dxdc;

	kzet[k][j][i].x = dydc * dzde - dzdc * dyde;
	kzet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	kzet[k][j][i].z = dxdc * dyde - dydc * dxde;


	kaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	  dydc * (dxde * dzdz - dzde * dxdz) +
	  dzdc * (dxde * dydz - dyde * dxdz);
	kaj[k][j][i] = 1./kaj[k][j][i];
	
	#ifdef NEWMETRIC
	AxByC ( 0.5, lcsi[k][j][i], 0.5, lcsi[k+1][j][i], &kcsi[k][j][i]);
	AxByC ( 0.5, leta[k][j][i], 0.5, leta[k+1][j][i], &keta[k][j][i]);
	AxByC ( 0.5, lzet[k][j][i], 0.5, lzet[k+1][j][i], &kzet[k][j][i]);
	kaj[k][j][i] = 2. / ( 1/laj[k][j][i] + 1/laj[k+1][j][i] );
	#endif

      }
    }
  }
  
	DAVecRestoreArray(fda, user->lCsi, &lcsi);
	DAVecRestoreArray(fda, user->lEta, &leta);
	DAVecRestoreArray(fda, user->lZet, &lzet);
	DAVecRestoreArray(da, user->lAj,  &laj);

/* ==================================================================================             */
/*  Calculating the Area of each side of the grid */
/*  the numbering is the same as flux (Ucont) */
  
  /*
  PetscReal  x13,x24,y13,y24,z13,z24;
  PetscReal  Ar1,Ar2,Ar3;
  Cmpnts     ***area;
  
  DAVecGetArray(fda, user->Area, &area);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	x13= coor[k][j][i].x-coor[k-1][j-1][i].x;
	y13= coor[k][j][i].y-coor[k-1][j-1][i].y;
	z13= coor[k][j][i].z-coor[k-1][j-1][i].z;

	x24= coor[k][j-1][i].x-coor[k-1][j][i].x;
	y24= coor[k][j-1][i].y-coor[k-1][j][i].y;
	z24= coor[k][j-1][i].z-coor[k-1][j][i].z;

	Ar1 =  y13*z24 - z13*y24 ;
	Ar2 =-(x13*z24 - z13*x24);
	Ar3 =  x13*y24 - y13*x24 ;
	
	area[k][j][i].x = 0.5*sqrt( Ar1*Ar1+Ar2*Ar2+Ar3*Ar3 );
	      
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	x13= coor[k][j][i].x-coor[k-1][j][i-1].x;
	y13= coor[k][j][i].y-coor[k-1][j][i-1].y;
	z13= coor[k][j][i].z-coor[k-1][j][i-1].z;

	x24= coor[k][j][i-1].x-coor[k-1][j][i].x;
	y24= coor[k][j][i-1].y-coor[k-1][j][i].y;
	z24= coor[k][j][i-1].z-coor[k-1][j][i].z;

	Ar1 =  y13*z24 - z13*y24 ;
	Ar2 =-(x13*z24 - z13*x24);
	Ar3 =  x13*y24 - y13*x24 ;

	area[k][j][i].y = 0.5*sqrt( Ar1*Ar1+Ar2*Ar2+Ar3*Ar3 );
	
      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	x13= coor[k][j][i].x-coor[k][j-1][i-1].x;
	y13= coor[k][j][i].y-coor[k][j-1][i-1].y;
	z13= coor[k][j][i].z-coor[k][j-1][i-1].z;

	x24= coor[k][j-1][i].x-coor[k][j][i-1].x;
	y24= coor[k][j-1][i].y-coor[k][j][i-1].y;
	z24= coor[k][j-1][i].z-coor[k][j][i-1].z;

	Ar1 =  y13*z24 - z13*y24 ;
	Ar2 =-(x13*z24 - z13*x24);
	Ar3 =  x13*y24 - y13*x24 ;

	area[k][j][i].z = 0.5*sqrt( Ar1*Ar1+Ar2*Ar2+Ar3*Ar3 );
	
      }
    }
  }

  DAVecRestoreArray(fda, user->Area, &area);

  DAGlobalToLocalBegin(fda, user->Area, INSERT_VALUES, user->lArea);
  DAGlobalToLocalEnd(fda, user->Area, INSERT_VALUES, user->lArea);

  VecDestroy(user->Area);
*/

/* ==================================================================================             */
/*  Calculating the Volume of each grid cell*/
/*  the numbering is the same as cell center (nvert, ucat) */
/*
  Cmpnts      v1,v2,v3;
  PetscReal   ***vol;

  DAVecGetArray(da, user->Volume, &vol);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	v1.x = centx[k][j][i].x - centx[k][j][i-1].x;
	v1.y = centx[k][j][i].y - centx[k][j][i-1].y;
	v1.z = centx[k][j][i].z - centx[k][j][i-1].z;

	v2.x = centy[k][j][i].x - centy[k][j-1][i].x;
	v2.y = centy[k][j][i].y - centy[k][j-1][i].y;
	v2.z = centy[k][j][i].z - centy[k][j-1][i].z;

	v3.x = centz[k][j][i].x - centz[k-1][j][i].x;
	v3.y = centz[k][j][i].y - centz[k-1][j][i].y;
	v3.z = centz[k][j][i].z - centz[k-1][j][i].z;

	// Volume = v1.(v2xv3)
	vol[k][j][i] = v1.x*(v2.y*v3.z-v3.y*v2.z)
	              -v1.y*(v2.x*v3.z-v3.x*v2.z) +
           	       v1.z*(v2.x*v3.y-v3.x*v2.y);

	vol[k][j][i] = fabs(vol[k][j][i]);

      }
    }
  }

  DAVecRestoreArray(da, user->Volume, &vol);

  DAGlobalToLocalBegin(da, user->Volume, INSERT_VALUES, user->lVolume);
  DAGlobalToLocalEnd(da, user->Volume, INSERT_VALUES, user->lVolume);

  VecDestroy(user->Volume);
  
  */
  
/* ==================================================================================             */

  DAVecRestoreArray(fda, Centz, &centz);
  DAVecRestoreArray(fda, Centy, &centy);
  DAVecRestoreArray(fda, Centx, &centx);
/*
  DAVecRestoreArray(cda, Csi, &csi);
  DAVecRestoreArray(cda, Eta, &eta);
  DAVecRestoreArray(cda, Zet, &zet);
  DAVecRestoreArray(da, Aj,  &aj);
*/
  DAVecRestoreArray(fda, KCsi, &kcsi);
  DAVecRestoreArray(fda, KEta, &keta);
  DAVecRestoreArray(fda, KZet, &kzet);
  DAVecRestoreArray(da, KAj,  &kaj);
  

  DAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);

  VecAssemblyBegin(user->Cent);
  VecAssemblyEnd(user->Cent);

  VecAssemblyBegin(user->ICsi);
  VecAssemblyEnd(user->ICsi);
  VecAssemblyBegin(user->IEta);
  VecAssemblyEnd(user->IEta);
  VecAssemblyBegin(user->IZet);
  VecAssemblyEnd(user->IZet);
  VecAssemblyBegin(user->IAj);
  VecAssemblyEnd(user->IAj);

  VecAssemblyBegin(user->JCsi);
  VecAssemblyEnd(user->JCsi);
  VecAssemblyBegin(user->JEta);
  VecAssemblyEnd(user->JEta);
  VecAssemblyBegin(user->JZet);
  VecAssemblyEnd(user->JZet);
  VecAssemblyBegin(user->JAj);
  VecAssemblyEnd(user->JAj);

  VecAssemblyBegin(user->KCsi);
  VecAssemblyEnd(user->KCsi);
  VecAssemblyBegin(user->KEta);
  VecAssemblyEnd(user->KEta);
  VecAssemblyBegin(user->KZet);
  VecAssemblyEnd(user->KZet);
  VecAssemblyBegin(user->KAj);
  VecAssemblyEnd(user->KAj);

  DARestoreLocalVector(fda, &Centx);

  VecDestroy(Centy);
  VecDestroy(Centz);

  DAGlobalToLocalBegin(fda, user->Csi, INSERT_VALUES, user->lCsi);
  DAGlobalToLocalEnd(fda, user->Csi, INSERT_VALUES, user->lCsi);

  DAGlobalToLocalBegin(fda, user->Eta, INSERT_VALUES, user->lEta);
  DAGlobalToLocalEnd(fda, user->Eta, INSERT_VALUES, user->lEta);

  DAGlobalToLocalBegin(fda, user->Zet, INSERT_VALUES, user->lZet);
  DAGlobalToLocalEnd(fda, user->Zet, INSERT_VALUES, user->lZet);

  DAGlobalToLocalBegin(fda, user->ICsi, INSERT_VALUES, user->lICsi);
  DAGlobalToLocalEnd(fda, user->ICsi, INSERT_VALUES, user->lICsi);

  DAGlobalToLocalBegin(fda, user->IEta, INSERT_VALUES, user->lIEta);
  DAGlobalToLocalEnd(fda, user->IEta, INSERT_VALUES, user->lIEta);

  DAGlobalToLocalBegin(fda, user->IZet, INSERT_VALUES, user->lIZet);
  DAGlobalToLocalEnd(fda, user->IZet, INSERT_VALUES, user->lIZet);

  DAGlobalToLocalBegin(fda, user->JCsi, INSERT_VALUES, user->lJCsi);
  DAGlobalToLocalEnd(fda, user->JCsi, INSERT_VALUES, user->lJCsi);

  DAGlobalToLocalBegin(fda, user->JEta, INSERT_VALUES, user->lJEta);
  DAGlobalToLocalEnd(fda, user->JEta, INSERT_VALUES, user->lJEta);

  DAGlobalToLocalBegin(fda, user->JZet, INSERT_VALUES, user->lJZet);
  DAGlobalToLocalEnd(fda, user->JZet, INSERT_VALUES, user->lJZet);

  DAGlobalToLocalBegin(fda, user->KCsi, INSERT_VALUES, user->lKCsi);
  DAGlobalToLocalEnd(fda, user->KCsi, INSERT_VALUES, user->lKCsi);

  DAGlobalToLocalBegin(fda, user->KEta, INSERT_VALUES, user->lKEta);
  DAGlobalToLocalEnd(fda, user->KEta, INSERT_VALUES, user->lKEta);

  DAGlobalToLocalBegin(fda, user->KZet, INSERT_VALUES, user->lKZet);
  DAGlobalToLocalEnd(fda, user->KZet, INSERT_VALUES, user->lKZet);

  DAGlobalToLocalBegin(da, user->Aj, INSERT_VALUES, user->lAj);
  DAGlobalToLocalEnd(da, user->Aj, INSERT_VALUES, user->lAj);

  DAGlobalToLocalBegin(da, user->IAj, INSERT_VALUES, user->lIAj);
  DAGlobalToLocalEnd(da, user->IAj, INSERT_VALUES, user->lIAj);

  DAGlobalToLocalBegin(da, user->JAj, INSERT_VALUES, user->lJAj);
  DAGlobalToLocalEnd(da, user->JAj, INSERT_VALUES, user->lJAj);

  DAGlobalToLocalBegin(da, user->KAj, INSERT_VALUES, user->lKAj);
  DAGlobalToLocalEnd(da, user->KAj, INSERT_VALUES, user->lKAj);

  DAGlobalToLocalBegin(fda, user->GridSpace, INSERT_VALUES, user->lGridSpace);
  DAGlobalToLocalEnd(fda, user->GridSpace, INSERT_VALUES, user->lGridSpace);

  DAGlobalToLocalBegin(fda, user->Cent, INSERT_VALUES, user->lCent);
  DAGlobalToLocalEnd(fda, user->Cent, INSERT_VALUES, user->lCent);

  VecDestroy(user->Csi);
  VecDestroy(user->Eta);
  VecDestroy(user->Zet);

  VecDestroy(user->ICsi);
  VecDestroy(user->IEta);
  VecDestroy(user->IZet);

  VecDestroy(user->JCsi);
  VecDestroy(user->JEta);
  VecDestroy(user->JZet);

  VecDestroy(user->KCsi);
  VecDestroy(user->KEta);
  VecDestroy(user->KZet);

  VecDestroy(user->Aj);
  VecDestroy(user->IAj);
  VecDestroy(user->JAj);
  VecDestroy(user->KAj);

  PetscBarrier(PETSC_NULL);
	//
	if(periodic) {
		
		DAVecGetArray(fda, user->lCsi, &csi);
		DAVecGetArray(fda, user->lEta, &eta);
		DAVecGetArray(fda, user->lZet, &zet);
		DAVecGetArray(da, user->lAj,  &aj);
	  
		DAVecGetArray(fda, user->lICsi, &icsi);
		DAVecGetArray(fda, user->lIEta, &ieta);
		DAVecGetArray(fda, user->lIZet, &izet);
		DAVecGetArray(da, user->lIAj,  &iaj);
		
		DAVecGetArray(fda, user->lJCsi, &jcsi);
		DAVecGetArray(fda, user->lJEta, &jeta);
		DAVecGetArray(fda, user->lJZet, &jzet);
		DAVecGetArray(da, user->lJAj,  &jaj);
		
		DAVecGetArray(fda, user->lKCsi, &kcsi);
		DAVecGetArray(fda, user->lKEta, &keta);
		DAVecGetArray(fda, user->lKZet, &kzet);
		DAVecGetArray(da, user->lKAj,  &kaj);
		
		DAVecGetArray(fda, user->lCent, &cent);
		/*
		if(xs==0 || xe==mx) {
			int from, to;
			for (k=lzs; k<lze; k++)
			for (j=lys; j<lye; j++) {
				if(xs==0) {
					i = 1, from = i, to = 0;
					
					if(i_periodic) from = mx-2;
					else if(ii_periodic) from = -2;
				}
				
				if(xe==mx) {
					i = mx-2, from = i, to = mx-1;
					
					if(i_periodic) from = 1;
					else if(ii_periodic) from = mx+1;
				}
				csi[k][j][to] = csi[k][j][from];
				eta[k][j][to] = eta[k][j][from];
				zet[k][j][to] = zet[k][j][from];
				aj[k][j][to] = aj[k][j][from];
				
				icsi[k][j][to] = icsi[k][j][from];
				ieta[k][j][to] = ieta[k][j][from];
				izet[k][j][to] = izet[k][j][from];
				iaj[k][j][to] = iaj[k][j][from];
				
				jcsi[k][j][to] = jcsi[k][j][from];
				jeta[k][j][to] = jeta[k][j][from];
				jzet[k][j][to] = jzet[k][j][from];
				jaj[k][j][to] = jaj[k][j][from];
				
				kcsi[k][j][to] = kcsi[k][j][from];
				keta[k][j][to] = keta[k][j][from];
				kzet[k][j][to] = kzet[k][j][from];
				kaj[k][j][to] = kaj[k][j][from];
			}
		}
		
		if(ys==0 || ye==my) {
			int from, to;
					
			for (k=lzs; k<lze; k++)
			for (i=lxs; i<lxe; i++) {
				if(ys==0) {
					j = 1, from = j, to = 0;
					
					if(j_periodic) from = my-2;
					else if(jj_periodic) from = -2;
				}
				
				if(ye==my) {
					j = my-2, from = j, to = my-1;
					
					if(j_periodic) from = 1;
					else if(jj_periodic) from = my+1;
				}
				
				csi[k][to][i] = csi[k][from][i];
				eta[k][to][i] = eta[k][from][i];
				zet[k][to][i] = zet[k][from][i];
				aj[k][to][i] = aj[k][from][i];
				
				icsi[k][to][i] = icsi[k][from][i];
				ieta[k][to][i] = ieta[k][from][i];
				izet[k][to][i] = izet[k][from][i];
				iaj[k][to][i] = iaj[k][from][i];
				
				jcsi[k][to][i] = jcsi[k][from][i];
				jeta[k][to][i] = jeta[k][from][i];
				jzet[k][to][i] = jzet[k][from][i];
				jaj[k][to][i] = jaj[k][from][i];
				
				kcsi[k][to][i] = kcsi[k][from][i];
				keta[k][to][i] = keta[k][from][i];
				kzet[k][to][i] = kzet[k][from][i];
				kaj[k][to][i] = kaj[k][from][i];				
			}
		}
		
		if(zs==0 || ze==mz) {
			int from, to;
			
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if(zs==0) {
					k = 1, from = k, to = 0;
					
					if(k_periodic) from = mz-2;
					else if(kk_periodic) from = -2;
				}
				
				if(ze==mz) {
					k = mz-2, from = k, to = mz-1;
					
					if(k_periodic) from = 1;
					else if(kk_periodic) from = mz+1;
				}
				
				csi[to][j][i] = csi[from][j][i];
				eta[to][j][i] = eta[from][j][i];
				zet[to][j][i] = zet[from][j][i];
				aj[to][j][i] = aj[from][j][i];
				
				icsi[to][j][i] = icsi[from][j][i];
				ieta[to][j][i] = ieta[from][j][i];
				izet[to][j][i] = izet[from][j][i];
				iaj[to][j][i] = iaj[from][j][i];
				
				jcsi[to][j][i] = jcsi[from][j][i];
				jeta[to][j][i] = jeta[from][j][i];
				jzet[to][j][i] = jzet[from][j][i];
				jaj[to][j][i] = jaj[from][j][i];
				
				kcsi[to][j][i] = kcsi[from][j][i];
				keta[to][j][i] = keta[from][j][i];
				kzet[to][j][i] = kzet[from][j][i];
				kaj[to][j][i] = kaj[from][j][i];			
			}
		}// needs to check using periodic open duct
		*/
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
				
				csi[k][j][i] = csi[c][b][a];
				eta[k][j][i] = eta[c][b][a];
				zet[k][j][i] = zet[c][b][a];
				aj[k][j][i] = aj[c][b][a];
				
				icsi[k][j][i] = icsi[c][b][a];
				ieta[k][j][i] = ieta[c][b][a];
				izet[k][j][i] = izet[c][b][a];
				iaj[k][j][i] = iaj[c][b][a];
				
				jcsi[k][j][i] = jcsi[c][b][a];
				jeta[k][j][i] = jeta[c][b][a];
				jzet[k][j][i] = jzet[c][b][a];
				jaj[k][j][i] = jaj[c][b][a];
				
				kcsi[k][j][i] = kcsi[c][b][a];
				keta[k][j][i] = keta[c][b][a];
				kzet[k][j][i] = kzet[c][b][a];
				kaj[k][j][i] = kaj[c][b][a];
				
				//cent[k][j][i] = cent[c][b][a];
			}
		}
		
		DAVecRestoreArray(fda, user->lCsi, &csi);
		DAVecRestoreArray(fda, user->lEta, &eta);
		DAVecRestoreArray(fda, user->lZet, &zet);
		DAVecRestoreArray(da, user->lAj,  &aj);
		
		DAVecRestoreArray(fda, user->lICsi, &icsi);
		DAVecRestoreArray(fda, user->lIEta, &ieta);
		DAVecRestoreArray(fda, user->lIZet, &izet);
		DAVecRestoreArray(da, user->lIAj,  &iaj);
		
		DAVecRestoreArray(fda, user->lJCsi, &jcsi);
		DAVecRestoreArray(fda, user->lJEta, &jeta);
		DAVecRestoreArray(fda, user->lJZet, &jzet);
		DAVecRestoreArray(da, user->lJAj,  &jaj);
		
		DAVecRestoreArray(fda, user->lKCsi, &kcsi);
		DAVecRestoreArray(fda, user->lKEta, &keta);
		DAVecRestoreArray(fda, user->lKZet, &kzet);
		DAVecRestoreArray(da, user->lKAj,  &kaj);
		
		DAVecRestoreArray(fda, user->lCent, &cent);
		
		// calling localtolocal is essential
		DALocalToLocalBegin(fda, user->lCsi, INSERT_VALUES, user->lCsi);
		DALocalToLocalEnd(fda, user->lCsi, INSERT_VALUES, user->lCsi);
		
		DALocalToLocalBegin(fda, user->lEta, INSERT_VALUES, user->lEta);
		DALocalToLocalEnd(fda, user->lEta, INSERT_VALUES, user->lEta);
		
		DALocalToLocalBegin(fda, user->lZet, INSERT_VALUES, user->lZet);
		DALocalToLocalEnd(fda, user->lZet, INSERT_VALUES, user->lZet);
		
		DALocalToLocalBegin(da, user->lAj, INSERT_VALUES, user->lAj);
		DALocalToLocalEnd(da, user->lAj, INSERT_VALUES, user->lAj);
		
		DALocalToLocalBegin(fda, user->lICsi, INSERT_VALUES, user->lICsi);
		DALocalToLocalEnd(fda, user->lICsi, INSERT_VALUES, user->lICsi);
		
		DALocalToLocalBegin(fda, user->lIEta, INSERT_VALUES, user->lIEta);
		DALocalToLocalEnd(fda, user->lIEta, INSERT_VALUES, user->lIEta);
		
		DALocalToLocalBegin(fda, user->lIZet, INSERT_VALUES, user->lIZet);
		DALocalToLocalEnd(fda, user->lIZet, INSERT_VALUES, user->lIZet);
		
		DALocalToLocalBegin(da, user->lIAj, INSERT_VALUES, user->lIAj);
		DALocalToLocalEnd(da, user->lIAj, INSERT_VALUES, user->lIAj);
		
		DALocalToLocalBegin(fda, user->lJCsi, INSERT_VALUES, user->lJCsi);
		DALocalToLocalEnd(fda, user->lJCsi, INSERT_VALUES, user->lJCsi);
		
		DALocalToLocalBegin(fda, user->lJEta, INSERT_VALUES, user->lJEta);
		DALocalToLocalEnd(fda, user->lJEta, INSERT_VALUES, user->lJEta);
		
		DALocalToLocalBegin(fda, user->lJZet, INSERT_VALUES, user->lJZet);
		DALocalToLocalEnd(fda, user->lJZet, INSERT_VALUES, user->lJZet);
		
		DALocalToLocalBegin(da, user->lJAj, INSERT_VALUES, user->lJAj);
		DALocalToLocalEnd(da, user->lJAj, INSERT_VALUES, user->lJAj);
		
		DALocalToLocalBegin(fda, user->lKCsi, INSERT_VALUES, user->lKCsi);
		DALocalToLocalEnd(fda, user->lKCsi, INSERT_VALUES, user->lKCsi);
		
		DALocalToLocalBegin(fda, user->lKEta, INSERT_VALUES, user->lKEta);
		DALocalToLocalEnd(fda, user->lKEta, INSERT_VALUES, user->lKEta);
		
		DALocalToLocalBegin(fda, user->lKZet, INSERT_VALUES, user->lKZet);
		DALocalToLocalEnd(fda, user->lKZet, INSERT_VALUES, user->lKZet);
		
		DALocalToLocalBegin(da, user->lKAj, INSERT_VALUES, user->lKAj);
		DALocalToLocalEnd(da, user->lKAj, INSERT_VALUES, user->lKAj);
	}	
	return 0;
}

