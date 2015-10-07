/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

extern PetscInt ti;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, inviscid;
extern PetscInt TwoD;
extern int wallfunction;
extern double mean_pressure_gradient;

extern int momentum_option;
extern PetscInt les, rans, mixed;


extern double integrate_gridfilter(double val[3][3][3]);
extern void Calculate_Covariant_metrics(double g[3][3], double G[3][3]);

extern void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3]);
extern void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3]);
/* Reconstruct Cartesian velocity components at cell centers from
   volume fluxes defined at cell surface centers
   Input: 
         *user
	 Ucont --------- Local Vec of surface volume fluxes
	 Ucat  --------- Global Vec of Cartesian velocity
*/

void Contra2Cart(UserCtx *user)
{
	Contra2Cart_2(user);
	//DALocalToLocalBegin(user->da, user->lUstar, INSERT_VALUES, user->lUstar);
	//DALocalToLocalEnd(user->da, user->lUstar, INSERT_VALUES, user->lUstar);
}

void Contra2Cart_single(Cmpnts &csi, Cmpnts &eta, Cmpnts &zet, Cmpnts &ucont, Cmpnts *ucat)
{
	double det = csi.x * (eta.y * zet.z - eta.z * zet.y) -
		csi.y * (eta.x * zet.z - eta.z * zet.x) +
		csi.z * (eta.x * zet.y - eta.y * zet.x);

	double det0 = ucont.x * (eta.y * zet.z - eta.z * zet.y) -
		ucont.y * (csi.y * zet.z - csi.z * zet.y) +
		ucont.z * (csi.y * eta.z - csi.z * eta.y);

	double det1 = -ucont.x * (eta.x * zet.z - eta.z * zet.x) +
		ucont.y * (csi.x * zet.z - csi.z * zet.x) -
		ucont.z * (csi.x * eta.z - csi.z * eta.x);

	double det2 = ucont.x * (eta.x * zet.y - eta.y * zet.x) -
		ucont.y * (csi.x * zet.y - csi.y * zet.x) +
		ucont.z * (csi.x * eta.y - csi.y * eta.x);

	(*ucat).x = det0 / det;
	(*ucat).y = det1 / det;
	(*ucat).z = det2 / det;		
};

void Contra2Cart_2(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	PetscReal	mat[3][3], det, det0, det1, det2;
	
	Vec		Aj = user->lAj;
	Vec Coor;
	
	PetscReal	***aj, ***rho, ***mu, ***llevel;
	Cmpnts	***ucat, ***lucat_o;
	Cmpnts ***lucat, ***lucont;
	
	PetscReal	***nvert;
	PetscReal	q[3]; //local working array
	PetscInt	i, j, k;
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

	if (lxs==0) lxs++;
	if (lxe==mx) lxe--;
	if (lys==0) lys++;
	if (lye==my) lye--;
	if (lzs==0) lzs++;
	if (lze==mz) lze--;

	Cmpnts	***icsi, ***jeta, ***kzet;
	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***coor;
	
	DAGetGhostedCoordinates(da, &Coor);
	DAVecGetArray(fda, Coor, &coor);
	
	if(levelset) {
		DAVecGetArray(da, user->lMu, &mu);
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lLevelset, &llevel);
	}

	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lKZet, &kzet);
	
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	
	DAVecGetArray(da,  Aj,  &aj);
	DAVecGetArray(da, user->lNvert, &nvert);

	PetscInt rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	DAVecGetArray(fda, user->lUcont, &lucont);
	DAVecGetArray(fda, user->Ucat,  &ucat);
	DAVecGetArray(fda, user->lUcat_old,  &lucat_o);
	
	// important
	if(periodic)
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
			lucont[k][j][i] = lucont[c][b][a];
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ( nvert[k][j][i] < 0.1 ) {
			mat[0][0] = (csi[k][j][i].x);
			mat[0][1] = (csi[k][j][i].y);
			mat[0][2] = (csi[k][j][i].z);
					      
			mat[1][0] = (eta[k][j][i].x);
			mat[1][1] = (eta[k][j][i].y);
			mat[1][2] = (eta[k][j][i].z);
					      
			mat[2][0] = (zet[k][j][i].x);
			mat[2][1] = (zet[k][j][i].y);
			mat[2][2] = (zet[k][j][i].z);
			
			
			int iLL =i-2,  iL =i-1,  iR= i,  iRR =i+1;
			int jLL =j-2,  jL =j-1,  jR= j,  jRR =j+1;
			int kLL=k-2, kL=k-1, kR=k, kRR=k+1;
			/*
			if(i==1) {
				if(i_periodic) iLL = mx-3, iL = mx-2;
				else iLL = iL, iRR = iR;
			}
			else if(i==mx-2) {
				if(i_periodic) iRR = 1;
				else iLL = iL, iRR = iR;
			}
			else if( nvert[k][j][iLL]+nvert[k][j][iRR] > 0.1 ) iLL = iL, iRR = iR;
			
			if(j==1) {
				if(j_periodic) jLL = my-3, jL = my-2;
				else jLL = jL, jRR = jR;
			}
			else if(j==my-2) {
				if(j_periodic) jRR = 1;
				else jLL = jL, jRR = jR;
			}
			else if( nvert[k][jLL][i]+nvert[k][jRR][i] > 0.1 ) jLL = jL, jRR = jR;
			
			if(k==1) {
				if(k_periodic) kLL = mz-3, kL = mz-2;
				else if(kk_periodic) kLL = -3, kL = -2;
				else kLL = kL, kRR = kR;
			}
			else if(k==mz-2) {
				if(k_periodic) kRR = 1;
				else if(kk_periodic) kRR = mz+1;
				else kLL = kL, kRR = kR;
			}
			else if( nvert[kLL][j][i]+nvert[kRR][j][i] > 0.1 ) kLL = kL, kRR = kR;
			*/
			q[0] = 0.5 * ( lucont[k][j][iL].x + lucont[k][j][iR].x);
			q[1] = 0.5 * ( lucont[k][jL][i].y + lucont[k][jR][i].y);
			q[2] = 0.5 * ( lucont[kL][j][i].z + lucont[kR][j][i].z);
			
			/*
			q[0] = 0.0625 * ( -ucont[k][j][iLL].x + 9.*ucont[k][j][iL].x + 9.*ucont[k][j][iR].x - ucont[k][j][iRR].x );
			q[1] = 0.0625 * ( -ucont[k][jLL][i].y + 9.*ucont[k][jL][i].y + 9.*ucont[k][jR][i].y - ucont[k][jRR][i].y );
			q[2] = 0.0625 * ( -ucont[kLL][j][i].z + 9.*ucont[kL][j][i].z + 9.*ucont[kR][j][i].z - ucont[kRR][j][i].z );
			*/
					
			det = 	mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
					mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
					mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

			det0 =	 q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
					q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
					q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

			det1 = 	-q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
					q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
					q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

			det2 = 	q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
					q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
					q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

			ucat[k][j][i].x = det0 / det;
			ucat[k][j][i].y = det1 / det;
			ucat[k][j][i].z = det2 / det;
			/*
			if (levelset && k==1 && user->bctype[4]==5 && 
			    (llevel[k][j][i]<0 || llevel[k][j][i]*llevel[k][j-1][i]<0 || llevel[k][j][i]*llevel[k][j+1][i]<0) ) {
			  Set(&ucat[k][j][i], 0);
			}
			*/
		}
	}
	
	DAVecRestoreArray(fda, user->Ucat,  &ucat);
	
	DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	
	if(periodic) {
		DAVecGetArray(fda, user->Ucat,  &ucat);
		DAVecGetArray(fda, user->lUcat,  &lucat);
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
				lucat[k][j][i] = lucat[c][b][a];
				ucat[k][j][i] = lucat[c][b][a];
			}
		}
		DAVecRestoreArray(fda, user->Ucat,  &ucat);
		DAVecRestoreArray(fda, user->lUcat,  &lucat);
	}
	
	//second call is needed
	DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
        DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

	//DALocalToGlobal(fda, user->lUcat, INSERT_VALUES, user->Ucat);
	
	DAVecGetArray(fda, user->lUcat,  &lucat);
	DAVecGetArray(fda, user->Ucat,  &ucat);
    
	PetscReal ***ustar;
	DAVecGetArray(da, user->lUstar, &ustar);
	
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
	  int solid_flag=0;
		if ( !movefsi && !rotatefsi && (int)(nvert[k][j][i]+0.1)==3 ) {
			Set(&ucat[k][j][i], 0);
			continue;
		}
		
		// wall function for boundary
		if( j!=0 && j!=my-1 && k!=0 && k!=mz-1 && 
		    ( ( user->bctype[0]==-1 && i==1) || (user->bctype[1]==-1 &&  i==mx-2) ) ) {
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double sb, sc; 
			double ni[3], nj[3], nk[3];
			Cmpnts Uc, Ua, Ub;
			
			Ua.x = Ua.y = Ua.z = 0;
			sb = 0.5/aj[k][j][i]/area;
			
			if(i==1) {
				sc = 2*sb + 0.5/aj[k][j][i+1]/area;
				Uc = lucat[k][j][i+1];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j][i-1]/area;
				Uc = lucat[k][j][i-1];
			}
			
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
			
			double ren = user->ren;
			
			wall_function (1./ren, sc, sb, Ua, Uc, &Ub, &ustar[k][j][i], ni[0], ni[1], ni[2]);
			
			ucat[k][j][i] = Ub;
			
			//if(user->bctype[2]<=1 && j==1) Set(&ucat[k][j][i],0), ustar[k][j][i]=0;
                        //if(user->bctype[3]<=1 && j==my-2) Set(&ucat[k][j][i],0), ustar[k][j][i]=0;
			//if(user->bctype[3]==10 && j==my-2) Set(&ucat[k][j][i],0), ustar[k][j][i]=0;
			solid_flag=1;
		}
		
		// rough wall function for boundary
		if( j!=0 && j!=my-1 && k!=0 && k!=mz-1 && 
		    ( ( user->bctype[0]==-2 && i==1) || (user->bctype[1]==-2 &&  i==mx-2) ) ) {
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double sb, sc; 
			double ni[3], nj[3], nk[3];
			Cmpnts Uc, Ua, Ub;
			
			Ua.x = Ua.y = Ua.z = 0;
			sb = 0.5/aj[k][j][i]/area;
			
			if(i==1) {
				sc = 2*sb + 0.5/aj[k][j][i+1]/area;
				Uc = lucat[k][j][i+1];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j][i-1]/area;
				Uc = lucat[k][j][i-1];
			}
			
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
			
			double ren = user->ren;
			
			wall_function_roughness_loglaw (1./ren, roughness_size, sc, sb, Ua, Uc, &Ub, &ustar[k][j][i], ni[0], ni[1], ni[2]);
			
			ucat[k][j][i] = Ub;
			solid_flag=1;
		}
		
		// wall function for boundary
		if( i!=0 && i!=mx-1 && k!=0 && k!=mz-1 && ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 &&  j==my-2) ) ) {
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double sb, sc; 
			double ni[3], nj[3], nk[3];
			Cmpnts Uc, Ua, Ub;
			
			Ua.x = Ua.y = Ua.z = 0;
			sb = 0.5/aj[k][j][i]/area;
			
			if(j==1) {
				sc = 2*sb + 0.5/aj[k][j+1][i]/area;
				Uc = lucat[k][j+1][i];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j-1][i]/area;
				Uc = lucat[k][j-1][i];
			}
			
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
			
			double ren = user->ren;
			
			wall_function (1./ren, sc, sb, Ua, Uc, &Ub, &ustar[k][j][i], nj[0], nj[1], nj[2]);
						
			ucat[k][j][i] = Ub;
			
			
			//if(user->bctype[0]<=1 && i==1) Set(&ucat[k][j][i],0), ustar[k][j][i]=0;
			//if(user->bctype[1]<=1 && i==mx-2) Set(&ucat[k][j][i],0), ustar[k][j][i]=0;
			solid_flag=1;
		}
		
		// rough wall function for boundary
		if( i!=0 && i!=mx-1 && k!=0 && k!=mz-1 && ( (user->bctype[2]==-2 && j==1) || (user->bctype[3]==-2 &&  j==my-2) ) ) {
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double sb, sc; 
			double ni[3], nj[3], nk[3];
			Cmpnts Uc, Ua, Ub;
			
			Ua.x = Ua.y = Ua.z = 0;
			sb = 0.5/aj[k][j][i]/area;
			
			if(j==1) {
				sc = 2*sb + 0.5/aj[k][j+1][i]/area;
				Uc = lucat[k][j+1][i];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j-1][i]/area;
				Uc = lucat[k][j-1][i];
			}
			
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
			
			double ren = user->ren;
			
			wall_function_roughness_loglaw (1./ren, roughness_size, sc, sb, Ua, Uc, &Ub, &ustar[k][j][i], nj[0], nj[1], nj[2]);
			
			ucat[k][j][i] = Ub;
			solid_flag=1;
		}
			
		//int solid_flag=0;
		
		if(user->bctype[3]==13 && j==my-1) { //slip top wall
		  ucat[k][j][i].x = ucat[k][j-1][i].x;
		  ucat[k][j][i].z = ucat[k][j-1][i].z;
		  ucat[k][j][i].y = -ucat[k][j-1][i].y;
		}
		if(user->bctype[3]==14 && j==my-1) { //slip top wall
                  ucat[k][j][i].x = ucat[k][j-1][i].x;
                  ucat[k][j][i].y = ucat[k][j-1][i].y;
                  ucat[k][j][i].z = -ucat[k][j-1][i].z;
                }
		/*slip BC*/
		if (user->bctype[0]==10 && i==0 && (j!=0 /*&& j!=my-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			double g[3][3], G[3][3];
			g[0][0]=csi[k][j][i+1].x, g[0][1]=csi[k][j][i+1].y, g[0][2]=csi[k][j][i+1].z;
			g[1][0]=eta[k][j][i+1].x, g[1][1]=eta[k][j][i+1].y, g[1][2]=eta[k][j][i+1].z;
			g[2][0]=zet[k][j][i+1].x, g[2][1]=zet[k][j][i+1].y, g[2][2]=zet[k][j][i+1].z;
			
			Calculate_Covariant_metrics(g, G);
			double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
			double nx = - xcsi, ny = - ycsi, nz = - zcsi;
			double sum=sqrt(nx*nx+ny*ny+nz*nz);
			nx /= sum, ny /= sum, nz /= sum;
			
			double un = lucat_o[k][j][i+1].x*nx + lucat_o[k][j][i+1].y*ny + lucat_o[k][j][i+1].z*nz;
			ucat[k][j][i].x = lucat_o[k][j][i+1].x - 2 * un * nx;
			ucat[k][j][i].y = lucat_o[k][j][i+1].y - 2 * un * ny;
			ucat[k][j][i].z = lucat_o[k][j][i+1].z - 2 * un * nz;
			
			if( nvert[k][j][i+1]>0.1 ) Set(&ucat[k][j][i],0);
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		if (user->bctype[1]==10 && i==mx-1 && (j!=0 /*&& j!=my-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			double g[3][3], G[3][3];
			g[0][0]=csi[k][j][i-1].x, g[0][1]=csi[k][j][i-1].y, g[0][2]=csi[k][j][i-1].z;
			g[1][0]=eta[k][j][i-1].x, g[1][1]=eta[k][j][i-1].y, g[1][2]=eta[k][j][i-1].z;
			g[2][0]=zet[k][j][i-1].x, g[2][1]=zet[k][j][i-1].y, g[2][2]=zet[k][j][i-1].z;
			
			Calculate_Covariant_metrics(g, G);
			double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
			double nx = xcsi, ny = ycsi, nz = zcsi;
			double sum=sqrt(nx*nx+ny*ny+nz*nz);
			nx /= sum, ny /= sum, nz /= sum;
			
			double un = lucat_o[k][j][i-1].x*nx + lucat_o[k][j][i-1].y*ny + lucat_o[k][j][i-1].z*nz;
			ucat[k][j][i].x = lucat_o[k][j][i-1].x - 2 * un * nx;
			ucat[k][j][i].y = lucat_o[k][j][i-1].y - 2 * un * ny;
			ucat[k][j][i].z = lucat_o[k][j][i-1].z - 2 * un * nz;
		
			if( nvert[k][j][i-1]>0.1 ) Set(&ucat[k][j][i],0);
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		if (user->bctype[2]==10 && j==0 && (i!=0/* && i!=mx-1*/ && k!=0/* && k!=mz-1*/) ) {
			double g[3][3], G[3][3];
			g[0][0]=csi[k][j+1][i].x, g[0][1]=csi[k][j+1][i].y, g[0][2]=csi[k][j+1][i].z;
			g[1][0]=eta[k][j+1][i].x, g[1][1]=eta[k][j+1][i].y, g[1][2]=eta[k][j+1][i].z;
			g[2][0]=zet[k][j+1][i].x, g[2][1]=zet[k][j+1][i].y, g[2][2]=zet[k][j+1][i].z;
			
			Calculate_Covariant_metrics(g, G);
			double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
			double nx = - xeta, ny = - yeta, nz = - zeta;
			double sum=sqrt(nx*nx+ny*ny+nz*nz);
			nx /= sum, ny /= sum, nz /= sum;
			
			double un = lucat_o[k][j+1][i].x*nx + lucat_o[k][j+1][i].y*ny + lucat_o[k][j+1][i].z*nz;
			ucat[k][j][i].x = lucat_o[k][j+1][i].x - 2 * un * nx;
			ucat[k][j][i].y = lucat_o[k][j+1][i].y - 2 * un * ny;
			ucat[k][j][i].z = lucat_o[k][j+1][i].z - 2 * un * nz;
			
			if( nvert[k][j+1][i]>0.1 ) Set(&ucat[k][j][i],0);
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		if (std::abs(user->bctype[3])==10 && j==my-1 && (i!=0 /*&& i!=mx-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			/*
			double ni[3], nj[3], nk[3];
			double nx, ny, nz;
			Calculate_normal(csi[k][j-1][i], eta[k][j-1][i], zet[k][j-1][i], ni, nj, nk);
			nx = nj[0];
			ny = nj[1];
			nz = nj[2];
		*/
			double g[3][3], G[3][3];
			g[0][0]=csi[k][j-1][i].x, g[0][1]=csi[k][j-1][i].y, g[0][2]=csi[k][j-1][i].z;
			g[1][0]=eta[k][j-1][i].x, g[1][1]=eta[k][j-1][i].y, g[1][2]=eta[k][j-1][i].z;
			g[2][0]=zet[k][j-1][i].x, g[2][1]=zet[k][j-1][i].y, g[2][2]=zet[k][j-1][i].z;
			
			Calculate_Covariant_metrics(g, G);
			double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
			double nx = xeta, ny = yeta, nz = zeta;
			double sum=sqrt(nx*nx+ny*ny+nz*nz);
			nx /= sum, ny /= sum, nz /= sum;
			
			double un = lucat_o[k][j-1][i].x*nx + lucat_o[k][j-1][i].y*ny + lucat_o[k][j-1][i].z*nz;
			ucat[k][j][i].x = lucat_o[k][j-1][i].x - 2 * un * nx;
			ucat[k][j][i].y = lucat_o[k][j-1][i].y - 2 * un * ny;
			ucat[k][j][i].z = lucat_o[k][j-1][i].z - 2 * un * nz;
		
			if( nvert[k][j-1][i]>0.1 ) Set(&ucat[k][j][i],0);
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		/* noslip BC */
		if(i==0 && user->bctype[0]==1 && (j!=0 /*&& j!=my-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k][j][i+1], &ucat[k][j][i]);
			//else AxC(-1, lucat[k][j][i+1], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		if(i==mx-1 && user->bctype[1]==1 && (j!=0 /*&& j!=my-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k][j][i-1], &ucat[k][j][i]);
			//else AxC(-1, lucat[k][j][i-1], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		
		if(j==0 && user->bctype[2]==1 && (i!=0 /*&& i!=mx-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k][j+1][i], &ucat[k][j][i]);
			//else AxC(-1, lucat[k][j+1][i], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		if(j==my-1 && user->bctype[3]==1 && (i!=0 /*&& i!=mx-1*/ && k!=0 /*&& k!=mz-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k][j-1][i], &ucat[k][j][i]);
			//else AxC(-1, lucat[k][j-1][i], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		
		if(k==0 && user->bctype[4]==1 && (i!=0 /*&& i!=mx-1*/ && j!=0 /*&& j!=my-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k+1][j][i], &ucat[k][j][i]);
			//else AxC(-1, lucat[k+1][j][i], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		if(k==mz-1 && user->bctype[5]==1  && (i!=0 /*&& i!=mx-1*/ && j!=0 /*&& j!=my-1*/) ) {
			/*if(solid_flag)*/ AxC(-1, lucat[k-1][j][i], &ucat[k][j][i]);
			//else AxC(-1, lucat[k-1][j][i], &ucat[k][j][i]);
			//if(solid_flag) Set(&ucat[k][j][i], 0);
			//else solid_flag=1;
			solid_flag=1;
		}
		
		//cavity problem 
		if (j==my-1 && user->bctype[3]==2) {
			if(solid_flag) {
				ucat[k][j][i].x=2.0-lucat_o[k][j-1][i].x;
				ucat[k][j][i].y=-lucat_o[k][j-1][i].y;
				ucat[k][j][i].z=-lucat_o[k][j-1][i].z;
				//Set(&ucat[k][j][i], 0);
			}
			else {
				ucat[k][j][i].x=2.0-lucat[k][j-1][i].x;
				ucat[k][j][i].y=-lucat[k][j-1][i].y;
				ucat[k][j][i].z=-lucat[k][j-1][i].z;
			}
		}
		
		// couette flow j=0
		if (j==0 && user->bctype[2]==12) {
			ucat[k][j][i].x=-lucat[k][j+1][i].x;
			ucat[k][j][i].y=-lucat[k][j+1][i].y;
			ucat[k][j][i].z=2.0-lucat[k][j+1][i].z;
			
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		// couette flow j=my-1
		if (j==my-1 && user->bctype[3]==12) {
			ucat[k][j][i].x=-lucat[k][j-1][i].x;
			ucat[k][j][i].y=-lucat[k][j-1][i].y;
			ucat[k][j][i].z=2.0-lucat[k][j-1][i].z;
			
			if(solid_flag) Set(&ucat[k][j][i], 0);
		}
		
		// body-fitted cylinder : inflow & outflow
		if (user->bctype[0]==11 && i==0 && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) {
			double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
			if( zc <= 0 ) {
				ucat[k][j][i].x = - lucat[k][j][i+1].x;
				ucat[k][j][i].y = - lucat[k][j][i+1].y;
				ucat[k][j][i].z = 2.0 - lucat[k][j][i+1].z;
				
				if(solid_flag) Set(&ucat[k][j][i], 0);
			}
		}
		
		

		/* outflow */
		if (user->bctype[3]==4 && j==my-1 && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1) ) {
		  ucat[k][j][i] = lucat[k][j-1][i];
		}
		/*
		if (user->bctype[4]==5 && k==0 && (i!=0 && i!=mx-1 && j!=0 && j!=my-1)){
		  Cmpnts u, ucon=lucont[k][j][i];
		  ucon.x = ucon.y = 0;
		  Contra2Cart_single(csi[k+1][j][i], eta[k+1][j][i], zet[k+1][j][i], ucon, &u);

		  ucat[k][j][i].x = - lucat[k+1][j][i].x + 2*u.x;
		  ucat[k][j][i].y = - lucat[k+1][j][i].y + 2*u.y;
		  ucat[k][j][i].z = - lucat[k+1][j][i].z + 2*u.z;
		}
		*/
		if (user->bctype[5]==4 && k==mz-1 && (i!=0 && i!=mx-1 && j!=0 && j!=my-1) ) {
		  /*
		  double ni[3], nj[3], nk[3];
		  double nx, ny, nz;
		  Calculate_normal(csi[k-1][j][i], eta[k-1][j][i], zet[k-1][j][i], ni, nj, nk);
		  nx = nk[0];
		  ny = nk[1];
		  nz = nk[2];
		  
		  double un = lucat[k-1][j][i].x * nx + lucat[k-1][j][i].y * ny + lucat[k-1][j][i].z * nz;
		  
		  ucat[k][j][i].x = - lucat[k-1][j][i].x + 2 * un * nx;
		  ucat[k][j][i].y = - lucat[k-1][j][i].y + 2 * un * ny;
		  ucat[k][j][i].z = - lucat[k-1][j][i].z + 2 * un * nz;
		  */
		  
		  if ( nvert[k-1][j][i]>0.1 ) Set(&ucat[k][j][i],0);
		  if (solid_flag) Set(&ucat[k][j][i], 0);
		}

		
		// 101220
		if (i!=0 && i!=mx-1 && j!=0 && j!=my-2 && k!=0 && k!=mz-2) {
		  if ( user->bctype[0]<=1 && user->bctype[2]<=1 && (i==1 && j==1) ) Set(&ucat[k][j][i],0);
		  if ( user->bctype[1]<=1 && user->bctype[2]<=1 && (i==mx-2 && j==1) ) Set(&ucat[k][j][i],0);
		  if ( user->bctype[0]<=1 && user->bctype[3]<=1 &&(i==1 && j==my-2) ) Set(&ucat[k][j][i],0);
		  if ( user->bctype[1]<=1 && user->bctype[3]<=1 &&(i==mx-2 && j==my-2) ) Set(&ucat[k][j][i],0);
		}
	}
	
	DAVecRestoreArray(da, user->lUstar, &ustar);
	DAVecRestoreArray(fda, user->lUcont, &lucont);
	DAVecRestoreArray(fda, user->lUcat,  &lucat);
	DAVecRestoreArray(fda, user->Ucat,  &ucat);
	DAVecRestoreArray(fda, user->lUcat_old,  &lucat_o);
	
	DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	
	if(levelset) {
		DAVecRestoreArray(da, user->lMu, &mu);
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lLevelset, &llevel);
        }

	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	
	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);
  
	DAVecRestoreArray(da,  Aj,  &aj);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	
	DAVecRestoreArray(fda, Coor, &coor);
	
	if(periodic) {
		DAVecGetArray(fda, user->Ucat,  &ucat);
		DAVecGetArray(fda, user->lUcat,  &lucat);
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
				lucat[k][j][i] = lucat[c][b][a];
				ucat[k][j][i] = lucat[c][b][a];
			}
		}
		DAVecRestoreArray(fda, user->Ucat,  &ucat);
		DAVecRestoreArray(fda, user->lUcat,  &lucat);

		DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat); //101229
		DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	}
};

PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv)
{
  //  Vec		Ucont = user->lUcont, Ucat = user->lUcat;
  
  Cmpnts	***ucont, ***ucat;
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***conv;
  //  Cmpnts	***zet;
  PetscReal	ucon, up, um;
  PetscReal	coef = 0.125;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	***nvert;

  DAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  DAVecGetArray(fda, Ucont, &ucont);
  DAVecGetArray(fda, Ucat,  &ucat);
  DAVecGetArray(fda, Conv,  &conv);

  VecDuplicate(Ucont, &Fp1);
  VecDuplicate(Ucont, &Fp2);
  VecDuplicate(Ucont, &Fp3);

  DAVecGetArray(fda, Fp1, &fp1);
  DAVecGetArray(fda, Fp2, &fp2);
  DAVecGetArray(fda, Fp3, &fp3);

  DAVecGetArray(da, user->lNvert, &nvert);

  /*  for (j=ys; j<ye; j++) {
    for (i=xs; i<xe; i++) {
      k=0;
      ucont[k][j][i].z = Flux_in * zet[k][j][i].z;
    }
  }
  if (ti==0) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	k=ze-2;
	ucont[k][j][i].z = Flux_in * zet[k][j][i].z;
      }
    }
  }
  else {
    PetscReal Flux_out=0;
    PetscReal fin = 0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	k=ze-3;
	Flux_out += ucont[k][j][i].z;
	k=0;
	fin += ucont[k][j][i].z;
      }
    }
    k=ze-2;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucont[k][j][i].z = ucont[k-1][j][i].z * fin / Flux_out;
      }
    }
    }*/
  /* We have two different sets of node: 1. grid node, the physical points
     where grid lines intercross; 2. storage node, where we store variables.
     All node without explicitly specified as "grid node" refers to
     storage node.

     The integer node is defined at cell center while half node refers to
     the actual grid node. (The reason to choose this arrangement is we need
     ghost node, which is half node away from boundaries, to specify boundary
     conditions. By using this storage arrangement, the actual storage need
     is (IM+1) * (JM + 1) * (KM+1) where IM, JM, & KM refer to the number of
     grid nodes along i, j, k directions.)

     DA, the data structure used to define the storage of 3D arrays, is defined
     as mx * my * mz. mx = IM+1, my = JM+1, mz = KM+1.

     Staggered grid arrangement is used in this solver.
     Pressure is stored at interger node (hence the cell center) and volume
     fluxes defined on the center of each surface of a given control volume
     is stored on the cloest upper integer node. */

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test a "); */
  /* Calculating the convective terms on cell centers.
     First calcualte the contribution from i direction
     The flux is evaluated by QUICK scheme */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){


	ucon = ucont[k][j][i].x * 0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (i>0 && i<mx-2 &&
	    (nvert[k][j][i+1]) < 0.1 &&
	    (nvert[k][j][i-1]) < 0.1) { // interial nodes
		fp1[k][j][i].x = 
			um * (coef * (-    ucat[k][j][i+2].x -  2. * ucat[k][j][i+1].x +  3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
			up * (coef * (-    ucat[k][j][i-1].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) + ucat[k][j][i  ].x);
		fp1[k][j][i].y = 
			um * (coef * (-    ucat[k][j][i+2].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) + ucat[k][j][i+1].y) +
			up * (coef * (-    ucat[k][j][i-1].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) + ucat[k][j][i  ].y);
		fp1[k][j][i].z = 
			um * (coef * (-    ucat[k][j][i+2].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
			up * (coef * (-    ucat[k][j][i-1].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) + ucat[k][j][i  ].z);
	}
	else if (i==0 || (nvert[k][j][i-1]) > 0.1) {
		fp1[k][j][i].x = 
			um * (coef * (-    ucat[k][j][i+2].x - 2. * ucat[k][j][i+1].x + 3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
			up * (coef * (-    ucat[k][j][i  ].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) +  ucat[k][j][i  ].x);
		fp1[k][j][i].y = 
			um * (coef * (-    ucat[k][j][i+2].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) +  ucat[k][j][i+1].y) +
			up * (coef * (-    ucat[k][j][i  ].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) +  ucat[k][j][i  ].y);
		fp1[k][j][i].z = 
			um * (coef * (-    ucat[k][j][i+2].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
			up * (coef * (-    ucat[k][j][i  ].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) +  ucat[k][j][i  ].z);
	}
	else if (i==mx-2 || (nvert[k][j][i+1]) > 0.1) {
		fp1[k][j][i].x = 
			um * (coef * (-    ucat[k][j][i+1].x -  2. * ucat[k][j][i+1].x +  3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
			up * (coef * (-    ucat[k][j][i-1].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) + ucat[k][j][i  ].x);
		fp1[k][j][i].y = 
			um * (coef * (-    ucat[k][j][i+1].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) + ucat[k][j][i+1].y) +
			up * (coef * (-    ucat[k][j][i-1].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) + ucat[k][j][i  ].y);
		fp1[k][j][i].z = 
			um * (coef * (-    ucat[k][j][i+1].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
			up * (coef * (-    ucat[k][j][i-1].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) + ucat[k][j][i  ].z);
	}
      }
    }
  }

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test aa "); */

  /* j direction */
  for (k=lzs; k<lze; k++) {
    for(j=lys-1; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].y * 0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (j>0 && j<my-2 &&
	    (nvert[k][j+1][i]) < 0.1 &&
	    (nvert[k][j-1][i]) < 0.1) {
		fp2[k][j][i].x = 
			um * (coef * (-    ucat[k][j+2][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) + ucat[k][j+1][i].x) +		     
			up * (coef * (-    ucat[k][j-1][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) + ucat[k][j][i  ].x);		     
		fp2[k][j][i].y = 				     
			um * (coef * (-    ucat[k][j+2][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
			up * (coef * (-    ucat[k][j-1][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) + ucat[k][j][i  ].y);		     
		fp2[k][j][i].z = 				     
			um * (coef * (-    ucat[k][j+2][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) + ucat[k][j+1][i].z) +		     
			up * (coef * (-    ucat[k][j-1][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) +  ucat[k][j][i  ].z);
	}
	else if (j==0 || (nvert[k][j-1][i]) > 0.1) {
		fp2[k][j][i].x = 
			um * (coef * (-    ucat[k][j+2][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) +  ucat[k][j+1][i].x) +		     
			up * (coef * (-    ucat[k][j  ][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) +  ucat[k][j][i  ].x);		     
		fp2[k][j][i].y = 				     
			um * (coef * (-    ucat[k][j+2][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
			up * (coef * (-    ucat[k][j  ][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) +  ucat[k][j][i  ].y);		     
		fp2[k][j][i].z = 				     
			um * (coef * (-    ucat[k][j+2][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +		     
			up * (coef * (-    ucat[k][j  ][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) + ucat[k][j][i  ].z);
	}
	else if (j==my-2 || (nvert[k][j+1][i]) > 0.1) {
		fp2[k][j][i].x = 
			um * (coef * (-    ucat[k][j+1][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) + ucat[k][j+1][i].x) +		     
			up * (coef * (-    ucat[k][j-1][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) + ucat[k][j][i  ].x);		     
		fp2[k][j][i].y = 				     
			um * (coef * (-    ucat[k][j+1][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
			up * (coef * (-    ucat[k][j-1][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) + ucat[k][j][i  ].y);		     
		fp2[k][j][i].z = 				     
			um * (coef * (-    ucat[k][j+1][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +		     
			up * (coef * (-    ucat[k][j-1][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) +  ucat[k][j][i  ].z);
	}
      }
    }
  }

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test aaa "); */

  /* k direction */
  for (k=lzs-1; k<lze; k++) {
    for(j=lys; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].z * 0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (k>0 && k<mz-2 &&
	    (nvert[k+1][j][i]) < 0.1 &&
	    (nvert[k-1][j][i]) < 0.1) {
		fp3[k][j][i].x = 
			um * (coef * (-    ucat[k+2][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
			up * (coef * (-    ucat[k-1][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) + ucat[k][j][i  ].x);  		     
		fp3[k][j][i].y = 	       			     
			um * (coef * (-    ucat[k+2][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) + ucat[k+1][j][i].y)   +		     
			up * (coef * (-    ucat[k-1][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) + ucat[k][j][i  ].y);  		     
		fp3[k][j][i].z = 	       			     
			um * (coef * (-    ucat[k+2][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
			up * (coef * (-    ucat[k-1][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);

	}
	else if (k<mz-2 && (k==0 || (nvert[k-1][j][i]) > 0.1)) {
		fp3[k][j][i].x = 
			um * (coef * (-    ucat[k+2][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
			up * (coef * (-    ucat[k  ][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);  		     
		fp3[k][j][i].y = 	       			     
			um * (coef * (-    ucat[k+2][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) + ucat[k+1][j][i].y)   +		     
			up * (coef * (-    ucat[k  ][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) +  ucat[k][j][i  ].y);  		     
		fp3[k][j][i].z = 	       			     
			um * (coef * (-    ucat[k+2][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
			up * (coef * (-    ucat[k  ][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
	}
	else if (k>0 && (k==mz-2 || (nvert[k+1][j][i]) > 0.1)) {
		fp3[k][j][i].x = 
			um * (coef * (-    ucat[k+1][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
			up * (coef * (-    ucat[k-1][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);  		     
		fp3[k][j][i].y = 	       			     
			um * (coef * (-    ucat[k+1][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) +  ucat[k+1][j][i].y)   +		     
			up * (coef * (-    ucat[k-1][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) + ucat[k][j][i  ].y);  		     
		fp3[k][j][i].z = 	       			     
			um * (coef * (-    ucat[k+1][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
			up * (coef * (-    ucat[k-1][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
	}
      }
    }
  }


  /*for (j=lys; j<lye; j++) {
    for (i=lxs; i<lxe; i++) {
      k=lze-1;
      fp3[k][j][i].x = fp3[k-1][j][i].x;
      fp3[k][j][i].y = fp3[k-1][j][i].y;
      fp3[k][j][i].z = fp3[k-1][j][i].z;

      
    }
    }*/

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test aaaa "); */

  /* Calculate the convective terms under cartesian coordinates */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	conv[k][j][i].x = 
	  fp1[k][j][i].x - fp1[k][j][i-1].x +
	  fp2[k][j][i].x - fp2[k][j-1][i].x +
	  fp3[k][j][i].x - fp3[k-1][j][i].x;
	conv[k][j][i].y = 
	  fp1[k][j][i].y - fp1[k][j][i-1].y +
	  fp2[k][j][i].y - fp2[k][j-1][i].y +
	  fp3[k][j][i].y - fp3[k-1][j][i].y;
	conv[k][j][i].z = 
	  fp1[k][j][i].z - fp1[k][j][i-1].z +
	  fp2[k][j][i].z - fp2[k][j-1][i].z +
	  fp3[k][j][i].z - fp3[k-1][j][i].z;

      }
    }
  }

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test aaaaa "); */

  DAVecRestoreArray(fda, Ucont, &ucont);
  DAVecRestoreArray(fda, Ucat,  &ucat);
  DAVecRestoreArray(fda, Conv,  &conv);

  DAVecRestoreArray(fda, Fp1, &fp1);
  DAVecRestoreArray(fda, Fp2, &fp2);
  DAVecRestoreArray(fda, Fp3, &fp3);
  DAVecRestoreArray(da, user->lNvert, &nvert);

  VecDestroy(Fp1);
  VecDestroy(Fp2);
  VecDestroy(Fp3);


  return (0);
}

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc)
{
	Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;

	Cmpnts	***ucont, ***ucat;//, ***ucat_old;

	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;

	PetscReal	***nvert;

	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	Vec		Fp1, Fp2, Fp3;
	Cmpnts	***fp1, ***fp2, ***fp3;
	Cmpnts	***visc;
	PetscReal	***aj, ***iaj, ***jaj, ***kaj;//, ***volume;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	PetscReal	ucon, up, um, ajc;

	PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
	PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
	PetscReal	g11, g21, g31;
	PetscReal	r11, r21, r31, r12, r22, r32, r13, r23, r33;

	PetscScalar	solid;

	solid = 0.5;
  
	DAGetLocalInfo(da, &info);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	DAVecGetArray(fda, Ucont, &ucont);
	DAVecGetArray(fda, Ucat,  &ucat);
	//DAVecGetArray(fda, user->lUcat_old,  &lucat_old);
	
	DAVecGetArray(fda, Visc,  &visc);

	DAVecGetArray(fda, Csi, &csi);
	DAVecGetArray(fda, Eta, &eta);
	DAVecGetArray(fda, Zet, &zet);

	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lIEta, &ieta);
	DAVecGetArray(fda, user->lIZet, &izet);

	DAVecGetArray(fda, user->lJCsi, &jcsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lJZet, &jzet);

	DAVecGetArray(fda, user->lKCsi, &kcsi);
	DAVecGetArray(fda, user->lKEta, &keta);
	DAVecGetArray(fda, user->lKZet, &kzet);

	DAVecGetArray(da, user->lNvert, &nvert);

	VecDuplicate(Ucont, &Fp1);
	VecDuplicate(Ucont, &Fp2);
	VecDuplicate(Ucont, &Fp3);

	DAVecGetArray(fda, Fp1, &fp1);
	DAVecGetArray(fda, Fp2, &fp2);
	DAVecGetArray(fda, Fp3, &fp3);

	DAVecGetArray(da, user->lAj, &aj);
	//DAVecGetArray(da, user->lVolume, &volume);
    
	PetscReal ***distance, ***lnu_t;	//***lS, 
	Cmpnts2 ***K_Omega;
	//Cmpnts ***lnu_t_IJK;
	
	//seokkoo
	if(les) {
		//DAVecGetArray(da, user->lCs, &Cs);
		//DAVecGetArray(fda, user->lNu_t_IJK, &lnu_t_IJK);
		DAVecGetArray(da, user->lNu_t, &lnu_t);
	}
	else if (rans) {
		DAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
		DAVecGetArray(user->da, user->Distance, &distance);
		//DAVecGetArray(user->da, user->lSrans, &lS);
		DAVecGetArray(da, user->lNu_t, &lnu_t);
	}
    
  /* The visc flux on each surface center is stored at previous integer node */

  // i direction
	DAVecGetArray(da, user->lIAj, &iaj);
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs-1; i<lxe; i++) {
		
		dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
		dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
		dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;

		if ((nvert[k][j+1][i])> solid || (nvert[k][j+1][i+1])> solid) {
			dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
			dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
			dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
		}
		else if  ((nvert[k][j-1][i])> solid || (nvert[k][j-1][i+1])> solid) {
			dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
			dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
			dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
		}
		else {
			dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
			dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
			dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
		}	  

		if ((nvert[k+1][j][i])> solid || (nvert[k+1][j][i+1])> solid) {
			dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
			dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
			dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
		}
		else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j][i+1])> solid) {
			dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
			dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
			dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
		}
		else {
			dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
			dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
			dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
		}
		
		csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;

		g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

		r11 = dudc * csi0 + dude * eta0 + dudz * zet0;	//du_dx
		r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;	//dv_dx
		r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;	//dw_dx

		r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
		r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
		r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

		r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
		r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
		r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

		ajc = iaj[k][j][i];
		
		double dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
				
		Compute_du_i (i, j, k, mx, my, mz, ucat, nvert, &dudc_o, &dvdc_o, &dwdc_o, &dude_o, &dvde_o, &dwde_o, &dudz_o, &dvdz_o, &dwdz_o);
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o,&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );	
		
		double nu = 1./user->ren, nu_t=0;
		
		if(les) {
			//nu_t = 0.5 * ( lnu_t_IJK[k][j][i].x + lnu_t_IJK[k][j][i+1].x );
			nu_t = 0.5 * ( lnu_t[k][j][i] + lnu_t[k][j][i+1] );
		}
		else if(rans && ti>0) {
			nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
		}
		
		fp1[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu+nu_t);
		fp1[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu+nu_t);
		fp1[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu+nu_t);
	}
  
	DAVecRestoreArray(da, user->lIAj, &iaj);  

  
  // j direction
	DAVecGetArray(da, user->lJAj, &jaj);
	for (k=lzs; k<lze; k++)
	for (j=lys-1; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ((nvert[k][j][i+1])> solid || (nvert[k][j+1][i+1])> solid) {
		  dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
		  dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
		  dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
		}
		else if ((nvert[k][j][i-1])> solid || (nvert[k][j+1][i-1])> solid) {
		  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
		  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
		  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
		}
		else {
		  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
		  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
		  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
		}

		dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
		dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
		dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;

		if ((nvert[k+1][j][i])> solid || (nvert[k+1][j+1][i])> solid) {
		  dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
		  dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
		  dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
		}
		else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j+1][i])> solid) {
		  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
		  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
		  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
		}
		else {
		  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
		  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
		  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
		}
		
		csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;

		g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
		g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

		r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
		r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
		r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

		r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
		r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
		r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

		r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
		r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
		r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

		ajc = jaj[k][j][i];
		
		double dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
				
		Compute_du_j (i, j, k, mx, my, mz, ucat, nvert, &dudc_o, &dvdc_o, &dwdc_o, &dude_o, &dvde_o, &dwde_o, &dudz_o, &dvdz_o, &dwdz_o);
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o,&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );	
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
		double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
		
		double nu = 1./user->ren, nu_t = 0;
		
		if(les) {
			//nu_t = 0.5 * ( lnu_t_IJK[k][j][i].y + lnu_t_IJK[k][j+1][i].y );
			nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
		}
		else if(rans && ti>0) {
			nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
		}
		
		fp2[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu+nu_t);
		fp2[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu+nu_t);
		fp2[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu+nu_t);
	}


	DAVecRestoreArray(da, user->lJAj, &jaj);
  
  
  // k direction
	DAVecGetArray(da, user->lKAj, &kaj);
	
	for (k=lzs-1; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ((nvert[k][j][i+1])> solid || (nvert[k+1][j][i+1])> solid) {
			dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
			dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
			dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
		}
		else if ((nvert[k][j][i-1])> solid || (nvert[k+1][j][i-1])> solid) {
			dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
			dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
			dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
		}
		else {
			dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
			dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
			dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
		}

		if ((nvert[k][j+1][i])> solid || (nvert[k+1][j+1][i])> solid) {
			dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
			dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
			dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
		}
		else if ((nvert[k][j-1][i])> solid || (nvert[k+1][j-1][i])> solid) {
			dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
			dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
			dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
		}
		else {
			dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
			dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
			dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
		}

		dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
		dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
		dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;

		csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		
		g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
		g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
		g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

		r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
		r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
		r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

		r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
		r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
		r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

		r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
		r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
		r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

		ajc = kaj[k][j][i];
	
		double dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		Compute_du_k (i, j, k, mx, my, mz, ucat, nvert, &dudc_o, &dvdc_o, &dwdc_o, &dude_o, &dvde_o, &dwde_o, &dudz_o, &dvdz_o, &dwdz_o);
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc_o, dvdc_o, dwdc_o, dude_o, dvde_o, dwde_o, dudz_o, dvdz_o, dwdz_o,&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );	
		
		double nu = 1./user->ren, nu_t =0;
		
		if(les) {
			//nu_t = 0.5 * ( lnu_t_IJK[k][j][i].z + lnu_t_IJK[k+1][j][i].z );
			nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
		}
		else if(rans && ti>0) {
			nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
		}
		
		fp3[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu+nu_t);
		fp3[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu+nu_t);
		fp3[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu+nu_t);
	}
	DAVecRestoreArray(da, user->lKAj, &kaj);

	DAVecRestoreArray(fda, Fp1, &fp1);
	DAVecRestoreArray(fda, Fp2, &fp2);
	DAVecRestoreArray(fda, Fp3, &fp3);
	
	DALocalToLocalBegin(fda, Fp1, INSERT_VALUES, Fp1);
	DALocalToLocalEnd(fda, Fp1, INSERT_VALUES, Fp1);
	DALocalToLocalBegin(fda, Fp2, INSERT_VALUES, Fp2);
	DALocalToLocalEnd(fda, Fp2, INSERT_VALUES, Fp2);
	DALocalToLocalBegin(fda, Fp3, INSERT_VALUES, Fp3);
	DALocalToLocalEnd(fda, Fp3, INSERT_VALUES, Fp3);
	
	DAVecGetArray(fda, Fp1, &fp1);
	DAVecGetArray(fda, Fp2, &fp2);
	DAVecGetArray(fda, Fp3, &fp3);
	
        DAVecGetArray(da, user->lIAj, &iaj);
        DAVecGetArray(da, user->lJAj, &jaj);
        DAVecGetArray(da, user->lKAj, &kaj);

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		visc[k][j][i].x = (fp1[k][j][i].x - fp1[k][j][i-1].x + fp2[k][j][i].x - fp2[k][j-1][i].x + fp3[k][j][i].x - fp3[k-1][j][i].x) ;	// for x-momentum
		visc[k][j][i].y = (fp1[k][j][i].y - fp1[k][j][i-1].y + fp2[k][j][i].y - fp2[k][j-1][i].y + fp3[k][j][i].y - fp3[k-1][j][i].y) ;	// for y-momentum
		visc[k][j][i].z = (fp1[k][j][i].z - fp1[k][j][i-1].z + fp2[k][j][i].z - fp2[k][j-1][i].z + fp3[k][j][i].z - fp3[k-1][j][i].z) ;		// for z-momentum
	}

	if(les) {
		//DAVecRestoreArray(da, user->lCs, &Cs);
		//DAVecRestoreArray(fda, user->lNu_t_IJK, &lnu_t_IJK);
		DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	}
	else if (rans) {
		DAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
		DAVecRestoreArray(user->da, user->Distance, &distance);
		//DAVecRestoreArray(user->da, user->lSrans, &lS);
		DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	}
	
        DAVecRestoreArray(da, user->lIAj, &iaj);
        DAVecRestoreArray(da, user->lJAj, &jaj);
        DAVecRestoreArray(da, user->lKAj, &kaj);

 
  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }

  if (ze==mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }
  if (ye==my) {
      j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	visc[k][j][i].x = 0; visc[k][j][i].y = 0; visc[k][j][i].z = 0;
      }
    }
  }

			   
	DAVecRestoreArray(fda, Ucont, &ucont);
	DAVecRestoreArray(fda, Ucat,  &ucat);
	//DAVecRestoreArray(fda, user->lUcat_old,  &lucat_old);
	DAVecRestoreArray(fda, Visc,  &visc);

	DAVecRestoreArray(fda, Csi, &csi);
	DAVecRestoreArray(fda, Eta, &eta);
	DAVecRestoreArray(fda, Zet, &zet);
	  
	DAVecRestoreArray(fda, Fp1, &fp1);
	DAVecRestoreArray(fda, Fp2, &fp2);
	DAVecRestoreArray(fda, Fp3, &fp3);

	DAVecRestoreArray(da, user->lAj, &aj);
	//DAVecRestoreArray(da, user->lVolume, &volume);
  
	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);

	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);

	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);

	DAVecRestoreArray(da, user->lNvert, &nvert);

	VecDestroy(Fp1);
	VecDestroy(Fp2);
	VecDestroy(Fp3);

  
	return(0);
};

//int time_marching=BDF;

PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
//  DA da = user->da, fda = user->fda;
  PetscInt istage;
  PetscReal alfa[4];

  PetscInt	bi;

  alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;
  // Create local working vector for Ucont and Ucat components

/*   DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*   DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
	
	COEF_TIME_ACCURACY=1.0;	//seokkoo
	//time_marching=0;

extern int outflow_scale;

  for (bi=0; bi<block_number; bi++) {
//    InflowFlux(&(user[bi]));
    //OutflowFlux(&(user[bi]));
    //FormBCS(&(user[bi]),&fsi[0]);

    //   PetscPrintf(PETSC_COMM_WORLD, "RK-1\n");

	InflowFlux(&(user[bi]));
	//	PetscPrintf(PETSC_COMM_WORLD, "RK-2\n");

	outflow_scale=0;
	FormBCS(&(user[bi]),&fsi[0]);
	//PetscPrintf(PETSC_COMM_WORLD, "RK-3\n");

    if (immersed) 
    /*for (ibi=0;ibi<NumberOfBodies;ibi++) */{
      ibm_interpolation_advanced(&user[bi]);
    }
    
    ///PetscPrintf(PETSC_COMM_WORLD, "RK-4\n");

    //    ibm_interpolation_advanced(&user[bi], ibm, fsi);
    //ibm_interpolation_advanced2(&user[bi], ibm);
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].dUcont);
    VecCopy(user[bi].Ucont, user[bi].pUcont);
  }

  // PetscPrintf(PETSC_COMM_WORLD, "RK-5\n");

  extern PetscErrorCode CalcRHS(UserCtx *user, int dudt);
  
    PetscInt pseudot;
    // pseudo time iteration
    for(pseudot=0; pseudot<1;pseudot++) {
      for (bi=0; bi<block_number; bi++) {
	for (istage=0; istage<4; istage++) {

	  //  PetscPrintf(PETSC_COMM_WORLD, "RK-6\n");

		Contra2Cart(user);
		//PetscPrintf(PETSC_COMM_WORLD, "RK-7\n");

		VecSet(user[bi].Rhs,0);	// seokkoo
		Formfunction_2(&user[bi], user[bi].Rhs, 1.0);
		VecAXPY(user[bi].Rhs, -1, user->dP);

		//PetscPrintf(PETSC_COMM_WORLD, "RK-8\n");

	  VecWAXPY(user[bi].Ucont, alfa[istage] * user[bi].dt * user[bi].st, user[bi].Rhs, user[bi].pUcont);

	  DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

		Contra2Cart(user);

		//	PetscPrintf(PETSC_COMM_WORLD, "RK-8-1\n");

	}
	VecCopy(user[bi].Ucont, user[bi].pUcont);
	//PetscPrintf(PETSC_COMM_WORLD, "RK-8-2\n");

      }
      //PetscPrintf(PETSC_COMM_WORLD, "RK-8-3\n");

      if (block_number>1) {
	Block_Interface_U(user);
      }
    }
    
    //PetscPrintf(PETSC_COMM_WORLD, "RK-9\n");

	DAGlobalToLocalBegin(user->fda, user[0].Ucont, INSERT_VALUES, user[0].lUcont);
	DAGlobalToLocalEnd(user->fda, user[0].Ucont, INSERT_VALUES, user[0].lUcont);
				
	//extern void update_ibm_ucont(UserCtx *user);
	//update_ibm_ucont(&user[0]);
		
	DAGlobalToLocalBegin(user[0].fda, user[0].Ucat, INSERT_VALUES, user[0].lUcat);
	DAGlobalToLocalEnd(user[0].fda, user[0].Ucat, INSERT_VALUES, user[0].lUcat);
    
	//PetscPrintf(PETSC_COMM_WORLD, "RK-10\n");

	outflow_scale=1;
	FormBCS(&(user[0]),&fsi[0]);
	
    
    //  }
  for (bi=0; bi<block_number; bi++) {
    VecDestroy(user[bi].Rhs);

    VecWAXPY(user[bi].DUold, -1., user[bi].Ucont_o, user[bi].Ucont);
/*     VecDestroy(user[bi].Ucont_o); */
    VecDestroy(user[bi].dUcont);
    VecDestroy(user[bi].pUcont);
  }
  //PetscPrintf(PETSC_COMM_WORLD, "RK-100\n");

  return 0;
}

PetscErrorCode Spectral(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscReal cfl = 0.5;
  PetscReal vnn = 0.5;

  PetscOptionsGetReal(PETSC_NULL, "-cfl", &cfl, PETSC_NULL);
  PetscInt i, j, k;
  PetscReal abu, abv, abw, ren = user->ren;
  Cmpnts ***ucont, ***csi, ***eta, ***zet;
  PetscReal ***dt, ***aj;
  PetscReal g1, g2, g3, eigen1, eigen2, eigen3, temp, dtii, dtij, dtik;
  PetscReal dtvi, dtvj, dtvk;
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(fda, user->lCsi, &csi);
  DAVecGetArray(fda, user->lEta, &eta);
  DAVecGetArray(fda, user->lZet, &zet);
  DAVecGetArray(da, user->lAj, &aj);
  DAVecGetArray(da, user->Dt, &dt);
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	g1 = (csi[k][j][i].x * csi[k][j][i].x +
	      csi[k][j][i].y * csi[k][j][i].y +
	      csi[k][j][i].z * csi[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g2 = (eta[k][j][i].x * eta[k][j][i].x +
	      eta[k][j][i].y * eta[k][j][i].y +
	      eta[k][j][i].z * eta[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g3 = (zet[k][j][i].x * zet[k][j][i].x +
	      zet[k][j][i].y * zet[k][j][i].y +
	      zet[k][j][i].z * zet[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	abu = fabs(ucont[k][j][i].x + ucont[k][j][i-1].x) * 0.5 * aj[k][j][i];
	abv = fabs(ucont[k][j][i].y + ucont[k][j][i-1].y) * 0.5 * aj[k][j][i];
	abw = fabs(ucont[k][j][i].z + ucont[k][j][i-1].z) * 0.5 * aj[k][j][i];

	if (g1<1.e-10) g1 = 1.;
	if (g2<1.e-10) g2 = 1.;
	if (g3<1.e-10) g3 = 1.;
	eigen1 = abu + sqrt(abu*abu + g1);
	eigen2 = abv + sqrt(abv*abv + g2);
	eigen3 = abw + sqrt(abw*abw + g3);

	dtii = cfl / eigen1;
	dtij = cfl / eigen2;
	dtik = cfl / eigen3;

	temp = vnn * ren;
	dtvi = temp / g1;
	dtvj = temp / g2;
	dtvk = temp / g3;

	temp = PetscMin(dtii, dtij);
	temp = PetscMin(temp, dtik);
	temp = PetscMin(temp, dtvi);
	temp = PetscMin(temp, dtvj);
	dt[k][j][i] = PetscMin(temp, dtvk);
      }
    }
  }
  DAVecRestoreArray(fda, user->lUcont, &ucont);
  DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(da, user->lAj, &aj);
  DAVecRestoreArray(da, user->Dt, &dt);
  return 0;
}
