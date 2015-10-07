/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

//#define SECOND_ORDER

extern PetscInt ti;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, inviscid;
extern PetscInt TwoD;
extern char path[256];

extern PetscInt les, rans, mixed;
extern IBMNodes	*ibm_ptr;
extern FSInfo *fsi_ptr;
//extern PetscErrorCode ibm_interpolation_advanced(UserCtx *user, /*IBMNodes *ibm, FSInfo *fsi,*//*PetscInt ibi, */PetscInt Add_dUndt);

// (A - B)*a = C
void Subtract_Scale_AddTo ( Cmpnts &A, Cmpnts &B, double a, Cmpnts *C) 
{
	(*C).x += (A.x - B.x) * a;
	(*C).y += (A.y - B.y) * a;
	(*C).z += (A.z - B.z) * a;
}

void Subtract_Scale_Set ( Cmpnts &A, Cmpnts &B, double a, Cmpnts *C) 
{
	(*C).x = (A.x - B.x) * a;
	(*C).y = (A.y - B.y) * a;
	(*C).z = (A.z - B.z) * a;
}

// C=aX+bY
void AxByC ( double a, Cmpnts &X, double b, Cmpnts &Y, Cmpnts *C) 
{
	(*C).x = a*X.x + b*Y.x;
	(*C).y = a*X.y + b*Y.y;
	(*C).z = a*X.z + b*Y.z;
}

// C=aX
void AxC ( double a, Cmpnts &X, Cmpnts *C) 
{
	(*C).x = a*X.x;
	(*C).y = a*X.y;
	(*C).z = a*X.z;
}

void Scale ( Cmpnts *A, double a )
{
	(*A).x *= a;
	(*A).y *= a;
	(*A).z *= a;
};

void Set ( Cmpnts *A, double a )
{
	(*A).x = a;
	(*A).y = a;
	(*A).z = a;
};

void Average(Cmpnts &a1, Cmpnts &a2, Cmpnts &a3, Cmpnts &a4, Cmpnts &a5, Cmpnts &a6, Cmpnts &a7, Cmpnts &a8, Cmpnts *A)
{
	(*A).x = 0.125 * ( a1.x + a2.x + a3.x + a4.x + a5.x + a6.x + a7.x + a8.x );
	(*A).y = 0.125 * ( a1.y + a2.y + a3.y + a4.y + a5.y + a6.y + a7.y + a8.y );
	(*A).z = 0.125 * ( a1.z + a2.z + a3.z + a4.z + a5.z + a6.z + a7.z + a8.z );
};

void Average(Cmpnts &a1, Cmpnts &a2, Cmpnts &a3, Cmpnts &a4, Cmpnts *A)
{
	(*A).x = 0.25 * ( a1.x + a2.x + a3.x + a4.x );
	(*A).y = 0.25 * ( a1.y + a2.y + a3.y + a4.y );
	(*A).z = 0.25 * ( a1.z + a2.z + a3.z + a4.z );
};

void Average(Cmpnts &a1, Cmpnts &a2, Cmpnts *A)
{
	(*A).x = 0.5 * ( a1.x + a2.x );
	(*A).y = 0.5 * ( a1.y + a2.y );
	(*A).z = 0.5 * ( a1.z + a2.z );
};

void Average(double &a1, double &a2, double &a3, double &a4, double &a5, double &a6, double &a7, double &a8, double *A)
{
	(*A) = 0.125 * ( a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 );
};

void Average(double &a1, double &a2, double &a3, double &a4, double *A)
{
	(*A) = 0.25 * ( a1 + a2 + a3 + a4 );
};

void Weighted_Average(double &a1, double &a2, double &a3, double &a4, double &w1, double &w2, double &w3, double &w4, double *A)
{
	(*A) = ( a1*w1 + a2*w2 + a3*w3 + a4*w4 ) / ( w1 + w2 + w3 + w4 );
};

void Average(double &a1, double &a2, double *A)
{
	(*A) = 0.5 * ( a1 + a2 );
};

double Calc_Minimum_dt (UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscReal	***aj;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	Cmpnts ***ucont, ***cent;
	Cmpnts	***csi, ***eta, ***zet;
	PetscReal ***nvert;

	DAGetLocalInfo(da, &info);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	DAVecGetArray(fda, user->lUcont,  &ucont);
	DAVecGetArray(fda, user->lCent,  &cent);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lNvert, &nvert);
	
	double ldt=1.e7, dt=0;
	double ldx=1.e7;
	double ldi_min=1.e7, ldj_min=1.e7, ldk_min=1.e7;
	double ldi_max=0, ldj_max=0, ldk_max=0;
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {

		ldt = PetscMin ( fabs(1./ucont[k][j][i].x/aj[k][j][i]), ldt );
		ldt = PetscMin ( fabs(1./ucont[k][j][i].y/aj[k][j][i]), ldt );
		ldt = PetscMin ( fabs(1./ucont[k][j][i].z/aj[k][j][i]), ldt );
		
		if(aj[k][j][i]<0) printf("Negative jacobian %d,%d,%d\n", i,j,k);

		double ldi = 1./aj[k][j][i]/sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
		double ldj = 1./aj[k][j][i]/sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
		double ldk = 1./aj[k][j][i]/sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
		
		ldi_min = PetscMin ( ldi_min, ldi );
		ldj_min = PetscMin ( ldj_min, ldj );
		ldk_min = PetscMin ( ldk_min, ldk );
		
		ldi_max = PetscMax ( ldi_max, ldi );
		ldj_max = PetscMax ( ldj_max, ldj );
		ldk_max = PetscMax ( ldk_max, ldk );
		
		ldx = PetscMin ( ldi_min, PetscMin ( ldj_min, ldk_min ) );
	}
	
	PetscGlobalMin(&ldt, &dt, PETSC_COMM_WORLD);
	PetscGlobalMin(&ldx, &dx_min, PETSC_COMM_WORLD);
	PetscGlobalMin(&ldi_min, &di_min, PETSC_COMM_WORLD);
	PetscGlobalMin(&ldj_min, &dj_min, PETSC_COMM_WORLD);
	PetscGlobalMin(&ldk_min, &dk_min, PETSC_COMM_WORLD);
	
	PetscGlobalMax(&ldi_max, &di_max, PETSC_COMM_WORLD);
	PetscGlobalMax(&ldj_max, &dj_max, PETSC_COMM_WORLD);
	PetscGlobalMax(&ldk_max, &dk_max, PETSC_COMM_WORLD);
	
	DAVecRestoreArray(fda, user->lUcont,  &ucont);
	DAVecRestoreArray(fda, user->lCent,  &cent);
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	
	PetscPrintf(PETSC_COMM_WORLD, "CFL 1.0 time step=%.6f, dx_min=%.5f\n", dt, dx_min);
	
	return dt;
}

void Pressure_Gradient(UserCtx *user, Vec dP)
{
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet, ***dp;
	PetscScalar ***p, ***level, ***rho;

	PetscReal	***nvert;

	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscReal	***iaj, ***jaj, ***kaj;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

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
	
	
	DAGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);
	
	if(periodic) {
		DAVecGetArray(da, user->lP, &p);
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
				p[k][j][i] = p[c][b][a];
			}
		}
		DAVecRestoreArray(da, user->lP, &p);
		
		DALocalToLocalBegin(da, user->lP, INSERT_VALUES, user->lP);
		DALocalToLocalEnd(da, user->lP, INSERT_VALUES, user->lP);
		DALocalToGlobal(da, user->lP, INSERT_VALUES, user->P);
	}
	
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
	DAVecGetArray(da, user->lP, &p);

	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);
	
	if(levelset) {
		DAVecGetArray(user->da, user->lLevelset, &level);
		DAVecGetArray(da, user->lDensity, &rho);
	}
	
	VecSet(dP, 0.);
	
	DAVecGetArray(fda, dP,  &dp);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double g11_i = (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y + icsi[k][j][i].z * icsi[k][j][i].z);
		double g12_i = (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y + ieta[k][j][i].z * icsi[k][j][i].z);
		double g13_i = (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y + izet[k][j][i].z * icsi[k][j][i].z);
		double g21_j = (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z);
		double g22_j = (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z);
		double g23_j = (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z);
		double g31_k = (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z);
		double g32_k = (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z);
		double g33_k = (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);
		double r1=1.0, r2=1.0, r3=1.0;
		
		if(levelset) {
			/*
			double lm;
			
			lm = 0.5 * ( level[k][j][i] + level[k][j][i+1] );
			r1 = rho_air + (rho_water - rho_air) * H ( lm );
			
			lm = 0.5 * ( level[k][j][i] + level[k][j+1][i] );
			r2 = rho_air + (rho_water - rho_air) * H ( lm );
			
			lm = 0.5 * ( level[k][j][i] + level[k+1][j][i] );
			r3 = rho_air + (rho_water - rho_air) * H ( lm );
			*/
			r1 = mean ( rho[k][j][i], rho[k][j][i+1] );
			r2 = mean ( rho[k][j][i], rho[k][j+1][i] );
			r3 = mean ( rho[k][j][i], rho[k+1][j][i] );
		}
		
		double dpdc, dpde, dpdz;
		if( (i==0 || i==mx-2) && i_periodic) dpdc = p[k][j][1] - p[k][j][i];
		else if( i==mx-2 && ii_periodic) dpdc = p[k][j][mx+1] - p[k][j][i];
		else dpdc = p[k][j][i+1] - p[k][j][i];
		
		if ((int)(nvert[k][j+1][i]+0.5)==1 || (int)(nvert[k][j+1][i+1]+0.5)==1 || (j==my-2 && !j_periodic && !jj_periodic)) dpde = (p[k][j  ][i  ] - p[k][j-1][i  ] + p[k][j  ][i+1] - p[k][j-1][i+1]) * 0.5;
		else if ((int)(nvert[k][j-1][i]+0.5)==1 || (int)(nvert[k][j-1][i+1]+0.5)==1 || (j==1 && !j_periodic && !jj_periodic)) dpde = (p[k][j+1][i  ] - p[k][j  ][i  ] + p[k][j+1][i+1] - p[k][j  ][i+1]) * 0.5;
		else dpde = (p[k][j+1][i] - p[k][j-1][i] + p[k][j+1][i+1] - p[k][j-1][i+1]) * 0.25;

		if ((int)(nvert[k+1][j][i]+0.5)==1 || (int)(nvert[k+1][j][i+1]+0.5)==1 || (k==mz-2 && !k_periodic && !kk_periodic)) dpdz = (p[k][j][i  ] - p[k-1][j][i  ] + p[k][j][i+1] - p[k-1][j][i+1]) * 0.5;
		else if ((int)(nvert[k-1][j][i]+0.5)==1 || (int)(nvert[k-1][j][i+1]+0.5)==1 || (k==1 && !k_periodic && !kk_periodic)) dpdz = (p[k+1][j][i  ] - p[k][j][i  ] + p[k+1][j][i+1] - p[k][j][i+1]) * 0.5;
		else dpdz = (p[k+1][j][i] - p[k-1][j][i] + p[k+1][j][i+1] - p[k-1][j][i+1]) * 0.25;
			
		dp[k][j][i].x = (dpdc * g11_i + dpde *  g12_i + dpdz * g13_i ) * iaj[k][j][i] / r1;
		
		if ((int)(nvert[k][j][i+1]+0.5)==1 || (int)(nvert[k][j+1][i+1]+0.5)==1 || (i==mx-2&& !i_periodic && !ii_periodic)) dpdc = (p[k][j  ][i] - p[k][j  ][i-1] + p[k][j+1][i] - p[k][j+1][i-1]) * 0.5;
		else if ((int)(nvert[k][j][i-1]+0.5)==1 || (int)(nvert[k][j+1][i-1]+0.5)==1 || (i==1&& !i_periodic && !ii_periodic)) dpdc = (p[k][j  ][i+1] - p[k][j  ][i] + p[k][j+1][i+1] - p[k][j+1][i]) * 0.5;
		else dpdc = (p[k][j  ][i+1] - p[k][j  ][i-1] + p[k][j+1][i+1] - p[k][j+1][i-1]) * 0.25;
		
		if( (j==0 || j==my-2) && j_periodic) dpde = p[k][1][i] - p[k][j][i];
		else if( j==my-2 && jj_periodic) dpde = p[k][my+1][i] - p[k][j][i];
		else dpde = p[k][j+1][i] - p[k][j][i];
		
		if ((int)(nvert[k+1][j][i]+0.5)==1  || (int)(nvert[k+1][j+1][i]+0.5)==1 || (k==mz-2 && !k_periodic && !kk_periodic)) dpdz = (p[k][j  ][i] - p[k-1][j  ][i] + p[k][j+1][i] - p[k-1][j+1][i]) * 0.5;
		else if ((int)(nvert[k-1][j][i]+0.5)==1  || (int)(nvert[k-1][j+1][i]+0.5)==1 || (k==1 && !k_periodic && !kk_periodic)) dpdz = (p[k+1][j  ][i] - p[k][j  ][i] + p[k+1][j+1][i] - p[k][j+1][i]) * 0.5;
		else dpdz = (p[k+1][j  ][i] - p[k-1][j  ][i] + p[k+1][j+1][i] - p[k-1][j+1][i]) * 0.25;
		
		dp[k][j][i].y = (dpdc * g21_j + dpde * g22_j + dpdz * g23_j ) * jaj[k][j][i] / r2;
		
		if ((int)(nvert[k][j][i+1]+0.5)==1 || (int)(nvert[k+1][j][i+1]+0.5)==1 || (i==mx-2 && !i_periodic && !ii_periodic)) dpdc = (p[k  ][j][i] - p[k  ][j][i-1] + p[k+1][j][i] - p[k+1][j][i-1]) * 0.5;
		else if ((int)(nvert[k][j][i-1]+0.5)==1 || (int)(nvert[k+1][j][i-1]+0.5)==1 || (i==1 && !i_periodic && !ii_periodic)) dpdc = (p[k  ][j][i+1] - p[k  ][j][i] + p[k+1][j][i+1] - p[k+1][j][i]) * 0.5;
		else dpdc = (p[k  ][j][i+1] - p[k  ][j][i-1] + p[k+1][j][i+1] - p[k+1][j][i-1]) * 0.25;

		if ((int)(nvert[k][j+1][i]+0.5) ==1 || (int)(nvert[k+1][j+1][i]+0.5)==1 || (j==my-2 && !j_periodic && !jj_periodic)) dpde = (p[k  ][j][i] - p[k  ][j-1][i] + p[k+1][j][i] - p[k+1][j-1][i]) * 0.5;
		else if ((int)(nvert[k][j-1][i]+0.5) ==1 || (int)(nvert[k+1][j-1][i]+0.5)==1 || (j==1 && !j_periodic && !jj_periodic)) dpde = (p[k  ][j+1][i] - p[k  ][j][i] + p[k+1][j+1][i] - p[k+1][j][i]) * 0.5;
		else dpde = (p[k  ][j+1][i] - p[k  ][j-1][i] + p[k+1][j+1][i] - p[k+1][j-1][i]) * 0.25;
		
		if( (k==0 || k==mz-2) && k_periodic) dpdz = p[1][j][i] - p[k][j][i];
		else if( k==mz-2 && kk_periodic) dpdz = p[mz+1][j][i] - p[k][j][i];
		else dpdz = (p[k+1][j][i] - p[k][j][i]);
		
		if(k_periodic || kk_periodic) {
			double dz = 1./ kaj[k][j][i] / sqrt(kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);
			if ( dpdz_set ) {
				dpdz += mean_pressure_gradient * dz;
			}
			else if(air_flow_levelset_periodic && level[k][j][i]<-dthick) {
				extern double inlet_flux;
				dpdz += dz * (user->mean_k_flux-inlet_flux) / user->dt / user->mean_k_area;
			}
			else if(inletprofile!=17) {
				extern double inlet_flux;
				dpdz += dz * (user->mean_k_flux-inlet_flux) / user->dt / user->mean_k_area;
			}
		}
		
		dp[k][j][i].z = (dpdc * g31_k + dpde * g32_k + dpdz * g33_k ) * kaj[k][j][i] / r3;
				
		if( i==0 || nvert[k][j][i]+nvert[k][j][i+1]>0.1 || (!i_periodic && !ii_periodic && i==mx-2) ) {
			dp[k][j][i].x = 0;
		}
		if( j==0 || nvert[k][j][i]+nvert[k][j+1][i]>0.1 || (!j_periodic && !jj_periodic && j==my-2) ) {
			dp[k][j][i].y = 0;
		}
		if( k==0 || nvert[k][j][i]+nvert[k+1][j][i]>0.1 || (!k_periodic && !kk_periodic && k==mz-2) ) {
			dp[k][j][i].z = 0;
		}
	}

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
	DAVecRestoreArray(da, user->lP, &p);

	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);
		
	if(levelset) {
		DAVecRestoreArray(user->da, user->lLevelset, &level);
		DAVecRestoreArray(da, user->lDensity, &rho);
	}
	
	DAVecRestoreArray(fda, dP,  &dp);
};

//add Toni
	// for wave_momentum _source
double epscos (double x, double eps)
{

	if ( x < -eps ) return 0;
	else if ( fabs(x) <= eps ) return 0.5 / eps * ( 1.0 + cos ( M_PI * x / eps ) );
	else return 0;
};	
//end Toni



PetscErrorCode Formfunction_2(UserCtx *user, Vec Rhs, double scale)
{
	Vec	Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
	Vec Ucont=user->lUcont, Ucat=user->lUcat;
	
	Cmpnts	***ucont, ***ucont_o, ***ucat, ***cent;

	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscScalar ***p;

	PetscReal	***nvert;

	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	Vec	Div1, Div2, Div3, Fp;
	Vec	Adv1, Adv2, Adv3;
	Vec	Visc1, Visc2, Visc3;
	
	Cmpnts	***div1, ***div2, ***div3, ***fp;
	Cmpnts	***adv1, ***adv2, ***adv3;
	Cmpnts	***visc1, ***visc2, ***visc3;
	Cmpnts	***rhs, ***stension;
	PetscReal	***aj, ***iaj, ***jaj, ***kaj, ***level, ***rho, ***mu;//, ***vol;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

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

	DAVecGetArray(fda, user->lCent, &cent);
	DAVecGetArray(fda, Ucont, &ucont);
	DAVecGetArray(fda, user->lUcont_o, &ucont_o);
	DAVecGetArray(fda, Ucat,  &ucat);
	DAVecGetArray(fda, Rhs,  &rhs);
//	DAVecGetArray(fda, user->RHS_o,  &rhs_o);
	//DAVecGetArray(fda, user->RHS_rm1,  &rhs_rm1);

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
	DAVecGetArray(da, user->lP, &p);
	
	if(levelset) {
	  DAVecGetArray(fda, user->lST, &stension);
		DAVecGetArray(user->da, user->lLevelset, &level);
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	
// add (Toni)
//	for wave_momentum_source
	PetscReal ***ustar;
	DAVecGetArray(da, user->lUstar, &ustar);
	Cmpnts ***WAVE_fp;
	double t=user->dt*(double)ti;
	if(wave_momentum_source) DAVecGetArray(fda, user->lWAVE_fp, &WAVE_fp);	
// end (Toni)

	Fp = user->Fp;
	Div1 = user->Div1;
	Div2 = user->Div2;
	Div3 = user->Div3;
	Visc1 = user->Visc1;
	Visc2 = user->Visc2;
	Visc3 = user->Visc3;
	
	if(user->bctype[0]==11) 
	{
		user->lA_cyl=0;
		user->lA_cyl_x=0;
		user->lA_cyl_z=0;
		user->lFvx_cyl=0;
		user->lFvz_cyl=0;
		user->lFpx_cyl=0;
		user->lFpz_cyl=0;
	}
	
	
	/*
	VecDuplicate(Ucont, &Fp);
	VecDuplicate(Ucont, &Div1);
	VecDuplicate(Ucont, &Div2);
	VecDuplicate(Ucont, &Div3);
	VecDuplicate(Ucont, &Visc1);
	VecDuplicate(Ucont, &Visc2);
	VecDuplicate(Ucont, &Visc3);
	*/
	
	VecSet(Div1, 0);
	VecSet(Div2, 0);
	VecSet(Div3, 0);
	VecSet(Visc1, 0);
	VecSet(Visc2, 0);
	VecSet(Visc3, 0);
	
	DAVecGetArray(fda, Div1, &div1);
	DAVecGetArray(fda, Div2, &div2);
	DAVecGetArray(fda, Div3, &div3);
	
	DAVecGetArray(fda, Visc1, &visc1);
	DAVecGetArray(fda, Visc2, &visc2);
	DAVecGetArray(fda, Visc3, &visc3);
	
	if(skew) {
		VecDuplicate(Ucont, &Adv1);
		VecDuplicate(Ucont, &Adv2);
		VecDuplicate(Ucont, &Adv3);
		
		DAVecGetArray(fda, Adv1, &adv1);
		DAVecGetArray(fda, Adv2, &adv2);
		DAVecGetArray(fda, Adv3, &adv3);
	}

	
	DAVecGetArray(da, user->lAj, &aj);
	//DAVecGetArray(da, user->lVolume, &vol);
    
	PetscReal ***lnu_t;
	Cmpnts2 ***K_Omega;
	
	//seokkoo
	if(les) {
		DAVecGetArray(da, user->lNu_t, &lnu_t);
	}
	else if (rans) {
		DAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
		//DAVecGetArray(user->da, user->Distance, &distance);
		DAVecGetArray(da, user->lNu_t, &lnu_t);
	}
      
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);

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
			ucont[k][j][i] = ucont[c][b][a];
		}
	}
	
	// i direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(j==0 || k==0) continue;
		
		PetscReal ajc = iaj[k][j][i];
		csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;

		Compute_du_i (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		
		g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

		r11 = dudc * csi0 + dude * eta0 + dudz * zet0;	//du_dx * J
		r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;	//dv_dx * J
		r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;	//dw_dx * J

		r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
		r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
		r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

		r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
		r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
		r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		int iL=i-1, iR=i+2;
		double denom=3.;
		double lm=1000.;
		
		if(levelset) lm = 0.5 * ( level[k][j][i] + level[k][j][i+1] );
		
		
		if(i==0 || i==mx-2) {
			if(i_periodic) iL = mx-3, iR = 2;
			else if(ii_periodic && i==mx-2) iR=mx+2;
			else if(ii_periodic && i==0) iL=-3;
			else iL = i, iR=i+1, denom=1.;
		}
		else if( nvert[k][j][iL]+nvert[k][j][iR] > 0.1 ) iL = i, iR=i+1;
		
		//#ifdef SECOND_ORDER
		if(second_order) {
		  iL = i, iR=i+1;
		  denom=1.;
		}
		//#endif
		
		double ucon = ucont[k][j][i].x;
		
		if(i_periodic && (i==0 || i==mx-2) )  ucon = ucont[k][j][mx-2].x;
		else if(ii_periodic && i==0)  ucon = ucont[k][j][-2].x;
		
		double up = - 0.5 * ( ucon + fabs(ucon) );
		double um = - 0.5 * ( ucon - fabs(ucon) );
		
		if ( immersed && i!=mx-2 && nvert[k][j][i]>0.1) {
			div1[k][j][i].x = 
				um * (0.125 * (-    ucat[k][j][i+2].x - 2. * ucat[k][j][i+1].x + 3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
				up * (0.125 * (-    ucat[k][j][i  ].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) +  ucat[k][j][i  ].x);
			div1[k][j][i].y = 
				um * (0.125 * (-    ucat[k][j][i+2].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) +  ucat[k][j][i+1].y) +
				up * (0.125 * (-    ucat[k][j][i  ].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) +  ucat[k][j][i  ].y);
			div1[k][j][i].z = 
				um * (0.125 * (-    ucat[k][j][i+2].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
				up * (0.125 * (-    ucat[k][j][i  ].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) +  ucat[k][j][i  ].z);
		}
		else if ( immersed &&  i!=0 && nvert[k][j][i+1]>0.1 ) {
			div1[k][j][i].x = 
				um * (0.125 * (-    ucat[k][j][i+1].x -  2. * ucat[k][j][i+1].x +  3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
				up * (0.125 * (-    ucat[k][j][i-1].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) + ucat[k][j][i  ].x);
			div1[k][j][i].y = 
				um * (0.125 * (-    ucat[k][j][i+1].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) + ucat[k][j][i+1].y) +
				up * (0.125 * (-    ucat[k][j][i-1].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) + ucat[k][j][i  ].y);
			div1[k][j][i].z = 
				um * (0.125 * (-    ucat[k][j][i+1].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
				up * (0.125 * (-    ucat[k][j][i-1].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) + ucat[k][j][i  ].z);
		}
		else {
		  if( inviscid || (rans && !central) || (les && (   (levelset_weno==1) ||  (levelset_weno==2 && (fabs(lm) <= dthick))   ||  (levelset_weno==3 && (lm > 0.))  ||   (levelset_weno==4 && (lm > -dthick))   ))  || (levelset_weno==5) ) {
			  if(1==0) {
				div1[k][j][i].x = 
					um * (0.125 * (-    ucat[k][j][iR].x -  2. * ucat[k][j][i+1].x +  3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
					up * (0.125 * (-    ucat[k][j][iL].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) + ucat[k][j][i  ].x);
				div1[k][j][i].y = 
					um * (0.125 * (-    ucat[k][j][iR].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) + ucat[k][j][i+1].y) +
					up * (0.125 * (-    ucat[k][j][iL].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) + ucat[k][j][i  ].y);
				div1[k][j][i].z = 
					um * (0.125 * (-    ucat[k][j][iR].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
					up * (0.125 * (-    ucat[k][j][iL].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) + ucat[k][j][i  ].z);
			  }
			  else {
				div1[k][j][i].x = - ucont[k][j][i].x * weno3(ucat[k][j][iL].x, ucat[k][j][i].x, ucat[k][j][i+1].x, ucat[k][j][iR].x, ucont[k][j][i].x);
				div1[k][j][i].y = - ucont[k][j][i].x * weno3(ucat[k][j][iL].y, ucat[k][j][i].y, ucat[k][j][i+1].y, ucat[k][j][iR].y, ucont[k][j][i].x);
				div1[k][j][i].z = - ucont[k][j][i].x * weno3(ucat[k][j][iL].z, ucat[k][j][i].z, ucat[k][j][i+1].z, ucat[k][j][iR].z, ucont[k][j][i].x);
			  }
		  }
			else {
			  //#ifdef SECOND_ORDER
			  if(second_order) {
				div1[k][j][i].x = - ucont[k][j][i].x * 0.5 * ( ucat[k][j][i].x + ucat[k][j][i+1].x);
				div1[k][j][i].y = - ucont[k][j][i].x * 0.5 * ( ucat[k][j][i].y + ucat[k][j][i+1].y);
				div1[k][j][i].z = - ucont[k][j][i].x * 0.5 * ( ucat[k][j][i].z + ucat[k][j][i+1].z);
			  }
			  else {
				//#else
				div1[k][j][i].x = -ucont[k][j][i].x * 0.0625 * ( -ucat[k][j][iL].x + 9.*ucat[k][j][i].x + 9.*ucat[k][j][i+1].x - ucat[k][j][iR].x );
				div1[k][j][i].y = -ucont[k][j][i].x * 0.0625 * ( -ucat[k][j][iL].y + 9.*ucat[k][j][i].y + 9.*ucat[k][j][i+1].y - ucat[k][j][iR].y );
				div1[k][j][i].z = -ucont[k][j][i].x * 0.0625 * ( -ucat[k][j][iL].z  + 9.*ucat[k][j][i].z + 9.*ucat[k][j][i+1].z - ucat[k][j][iR].z );
				//#endif
			  }
			}
		}
		
		if(skew) {
			div1[k][j][i].x *= 0.5;
			div1[k][j][i].y *= 0.5;
			div1[k][j][i].z *= 0.5;
			
			Cmpnts ducat_dc4;
			
			Subtract_Scale_Set(ucat[k][j][iR], ucat[k][j][iL], 1./denom, &ducat_dc4);
			
			adv1[k][j][i].x = - 0.5 * ucont[k][j][i].x * ( 9./8. * dudc - 1./8. * ducat_dc4.x );
			adv1[k][j][i].y = - 0.5 * ucont[k][j][i].x * ( 9./8. * dvdc - 1./8. * ducat_dc4.y );
			adv1[k][j][i].z = - 0.5 * ucont[k][j][i].x * ( 9./8. * dwdc - 1./8. * ducat_dc4.z );
		}
		
		if(nvert[k][j][i]+nvert[k][j][i+1]>0.1) {
			if(immersed==3 || !immersed) {
				div1[k][j][i].x = 0;
				div1[k][j][i].y = 0;
				div1[k][j][i].z = 0;
				if(skew) {
				  adv1[k][j][i].x = 0;
				  adv1[k][j][i].y = 0;
				  adv1[k][j][i].z = 0;
				}
			}
			/*
			if(skew) {
				adv1[k][j][i].x = 0;
				adv1[k][j][i].y = 0;
				adv1[k][j][i].z = 0;
				}*/
		}
		
		if(user->bctype[0]==11 && user->bctype[1]==1  && i==mx-2) {
			double Sxx = 0.5*( du_dx + du_dx ), Sxz = 0.5*(du_dz + dw_dx);
			double Szx = Sxz, Szz = 0.5*(dw_dz + dw_dz);
						
			double A = sqrt ( icsi[k][j][i].x*icsi[k][j][i].x + icsi[k][j][i].y*icsi[k][j][i].y + icsi[k][j][i].z*icsi[k][j][i].z );
			
			double ni[3], nj[3], nk[3];
			double Ai, Aj, Ak;
			
			Calculate_normal_and_area(icsi[k][j][i], ieta[k][j][i], izet[k][j][i], ni, nj, nk, &Ai, &Aj, &Ak);
			double nx = ni[0], ny = ni[1], nz = ni[2];
			nx *= -1;
			ny *= -1;
			nz *= -1;
			
			Ai=A;
			//printf("%f \n", user->A_cyl);
			user->lA_cyl += Ai;
			user->lA_cyl_x += fabs(icsi[k][j][i].x);//fabs( Ai * nx );
			user->lA_cyl_z += fabs(icsi[k][j][i].z);//fabs( Ai * nz );
			
			double P = 2*p[k][j][i-1] - p[k][j][i];
			user->lFpx_cyl += - P * nx * Ai;
			user->lFpz_cyl += - P * nz * Ai;
			user->lFvx_cyl += 2.0 * ( Sxx * nx + 0 * ny + Sxz * nz ) * Ai / user->ren;
			user->lFvz_cyl += 2.0 * ( Szx * nx + 0 * ny + Szz * nz ) * Ai / user->ren;
		}
		
		double nu = 1./user->ren, nu_t=0;
		
		if(levelset) {
			//double lm = 0.5 * ( level[k][j][i] + level[k][j][i+1] );
			//nu = mu_air + (mu_water - mu_air) * H ( lm );
			
			//nu = mean ( mu[k][j][i], mu[k][j][i+1] );
			
			if ( (i==0 && !i_periodic && !ii_periodic) || nvert[k][j][i]>0.1 ) nu = mu[k][j][i+1];
			else if ( (i==mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1]>0.1 ) nu = mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j][i+1] );
		}
			
		if( les || (rans && ti>0) ) {

						
			if ( (i==0 && !i_periodic && !ii_periodic) || nvert[k][j][i]>0.1 ) nu_t = lnu_t[k][j][i+1];
			else if( (i==mx-2 && !i_periodic && !ii_periodic) || nvert[k][j][i+1]>0.1 ) nu_t = lnu_t[k][j][i];
			else nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

			//if (i==0 || i==mx-2) nu_t=0;
			//if (j==0 || j==my-2) nu_t=nu_tau_wall;
			//if (k==0 || k==mz-2) nu_t=0;
			
			
			if(les && levelset) {
			  if(  ((fabs(lm) >= dthick) && (levelset_weno==2))  ||  (levelset_weno==0)  ||  ((levelset_weno==3)&&(lm < 0.))  ||  ((levelset_weno==4)&&(lm < -dthick))) {
					div1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu_t);
					div1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu_t);
					div1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu_t);
			  }
			}
			else {
				visc1[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu_t);
				visc1[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu_t);
				visc1[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu_t);
			}
		}
		else {
			visc1[k][j][i].x = 0;
			visc1[k][j][i].y = 0;
			visc1[k][j][i].z = 0;
		}
		
		if(laplacian) {
			r11=r21=r31=0.;
			r12=r22=r32=0.;
			r13=r23=r33=0.;
		}
		
		visc1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu);
		visc1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu);
		visc1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu);
		
		if(clark) {
			double dx, dy, dz;
			Calculate_dxdydz (ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
			double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
			
			double t11 = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 );
			double t12 = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 );
			double t13 = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 );
			double t21 = t12;
			double t22 = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 );
			double t23 = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 );
			double t31 = t13;
			double t32 = t23;
			double t33 = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 );
			
			visc1[k][j][i].x -= ( t11 * csi0 + t12 * csi1 + t13 * csi2 ) / 12.;
			visc1[k][j][i].y -= ( t21 * csi0 + t22 * csi1 + t23 * csi2 ) / 12.;
			visc1[k][j][i].z -= ( t31 * csi0 + t32 * csi1 + t33 * csi2 ) / 12.;
		}
		/*
		if( nvert[k][j][i]+nvert[k][j][i+1]>1.1 ) {
			Set (&div1[k][j][i], 0.);
			Set (&visc1[k][j][i], 0.);
			}*/
	}
  
	// j direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(i==0 || k==0) continue;

		PetscReal ajc = jaj[k][j][i];
		csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;

		Compute_du_j (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		
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

		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		int jL=j-1, jR=j+2;
		double denom=3.;
		
		double lm=1000.;
		
		if(levelset) lm = 0.5 * ( level[k][j][i] + level[k][j+1][i] );
		
		if(j==0 || j==my-2) {
			if(j_periodic) jL = my-3, jR = 2;
			else if(jj_periodic && j==my-2) jR=my+2;
			else if(jj_periodic && j==0) jL=-3;
			else jL = j, jR=j+1, denom=1.;
		}
		else if( nvert[k][jL][i]+nvert[k][jR][i] > 0.1 ) jL = j, jR=j+1;
		
		//#ifdef SECOND_ORDER
		if(second_order) {
		  jL = j, jR=j+1;
		  denom=1.;
		}
		//#endif
		
		double ucon = ucont[k][j][i].y;
		if(j_periodic && (j==0 || j==my-2) )  ucon = ucont[k][my-2][i].y;
		else if(jj_periodic && j==0)  ucon = ucont[k][-2][i].y;
		
		double up = - 0.5 * ( ucon + fabs(ucon) );
		double um = - 0.5 * ( ucon - fabs(ucon) );
		
		if ( immersed &&  j!=my-2 && nvert[k][j][i]>0.1 ) {
			div2[k][j][i].x = 
				um * (0.125 * (-    ucat[k][j+2][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) +  ucat[k][j+1][i].x) +		     
				up * (0.125 * (-    ucat[k][j  ][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) +  ucat[k][j][i  ].x);		     
			div2[k][j][i].y = 				     
				um * (0.125 * (-    ucat[k][j+2][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
				up * (0.125 * (-    ucat[k][j  ][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) +  ucat[k][j][i  ].y);		     
			div2[k][j][i].z = 				     
				um * (0.125 * (-    ucat[k][j+2][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +		     
				up * (0.125 * (-    ucat[k][j  ][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) + ucat[k][j][i  ].z);
		}
		else if ( immersed &&  j!=0 && nvert[k][j+1][i]>0.1 ) {
			div2[k][j][i].x = 
				um * (0.125 * (-    ucat[k][j+1][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) + ucat[k][j+1][i].x) +		     
				up * (0.125 * (-    ucat[k][j-1][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) + ucat[k][j][i  ].x);		     
			div2[k][j][i].y = 				     
				um * (0.125 * (-    ucat[k][j+1][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
				up * (0.125 * (-    ucat[k][j-1][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) + ucat[k][j][i  ].y);		     
			div2[k][j][i].z = 				     
				um * (0.125 * (-    ucat[k][j+1][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +		     
				up * (0.125 * (-    ucat[k][j-1][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) +  ucat[k][j][i  ].z);
		}
		else {
			if( inviscid || (rans && !central) || (les && (   (levelset_weno==1) ||  (levelset_weno==2 && (fabs(lm) <= dthick))   ||  (levelset_weno==3 && (lm > 0.))  ||   (levelset_weno==4 && (lm > -dthick))   ))  || (levelset_weno==5) ) {
			  if(1==0) {
				div2[k][j][i].x = 
					um * (0.125 * (-    ucat[k][jR][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) + ucat[k][j+1][i].x) +		     
					up * (0.125 * (-    ucat[k][jL][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) + ucat[k][j][i  ].x);		     
				div2[k][j][i].y = 				     
					um * (0.125 * (-    ucat[k][jR][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +		     
					up * (0.125 * (-    ucat[k][jL][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) + ucat[k][j][i  ].y);		     
				div2[k][j][i].z = 				     
					um * (0.125 * (-    ucat[k][jR][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) + ucat[k][j+1][i].z) +		     
					up * (0.125 * (-    ucat[k][jL][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) +  ucat[k][j][i  ].z);
			  }
			  else {
				
				div2[k][j][i].x = - ucont[k][j][i].y * weno3(ucat[k][jL][i].x, ucat[k][j][i].x, ucat[k][j+1][i].x, ucat[k][jR][i].x, ucont[k][j][i].y);
				div2[k][j][i].y = - ucont[k][j][i].y * weno3(ucat[k][jL][i].y, ucat[k][j][i].y, ucat[k][j+1][i].y, ucat[k][jR][i].y, ucont[k][j][i].y);
				div2[k][j][i].z = - ucont[k][j][i].y * weno3(ucat[k][jL][i].z, ucat[k][j][i].z, ucat[k][j+1][i].z, ucat[k][jR][i].z, ucont[k][j][i].y);
			  }
		  }
			else {
			  //#ifdef SECOND_ORDER
			  if(second_order) {
				div2[k][j][i].x = - ucont[k][j][i].y * 0.5 * (ucat[k][j][i].x + ucat[k][j+1][i].x);
				div2[k][j][i].y = - ucont[k][j][i].y * 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
				div2[k][j][i].z = - ucont[k][j][i].y * 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);
			  }
			  else {
				//#else
				div2[k][j][i].x = -ucont[k][j][i].y * 0.0625 * ( -ucat[k][jL][i].x + 9.*ucat[k][j][i].x + 9.*ucat[k][j+1][i].x - ucat[k][jR][i].x );
				div2[k][j][i].y = -ucont[k][j][i].y * 0.0625 * ( -ucat[k][jL][i].y + 9.*ucat[k][j][i].y + 9.*ucat[k][j+1][i].y - ucat[k][jR][i].y );
				div2[k][j][i].z = -ucont[k][j][i].y * 0.0625 * ( -ucat[k][jL][i].z  + 9.*ucat[k][j][i].z + 9.*ucat[k][j+1][i].z - ucat[k][jR][i].z );
			  }
				//#endif
			}
		}
		
		if(skew) {
			div2[k][j][i].x *= 0.5;
			div2[k][j][i].y *= 0.5;
			div2[k][j][i].z *= 0.5;
			
			Cmpnts ducat_de4;
			Subtract_Scale_Set(ucat[k][jR][i], ucat[k][jL][i], 1./denom, &ducat_de4);
			
			adv2[k][j][i].x = - 0.5 * ucont[k][j][i].y * ( 9./8. * dude - 1./8. * ducat_de4.x );
			adv2[k][j][i].y = - 0.5 * ucont[k][j][i].y * ( 9./8. * dvde - 1./8. * ducat_de4.y );
			adv2[k][j][i].z = - 0.5 * ucont[k][j][i].y * ( 9./8. * dwde - 1./8. * ducat_de4.z );
		}
		
		
		if( nvert[k][j][i]+nvert[k][j+1][i]>0.1) {
			if(immersed==3 || !immersed) {
				div2[k][j][i].x = 0;
				div2[k][j][i].y = 0;
				div2[k][j][i].z = 0;
				if(skew) {
				  adv2[k][j][i].x = 0;
				  adv2[k][j][i].y = 0;
				  adv2[k][j][i].z = 0;
				}
			}
			/*
			if(skew) {
				adv2[k][j][i].x = 0;
				adv2[k][j][i].y = 0;
				adv2[k][j][i].z = 0;
			}
			*/
		}
		
		double nu = 1./user->ren, nu_t = 0;
		
		if(levelset) {
			//double lm = 0.5 * ( level[k][j][i] + level[k][j+1][i] );
			//nu = mu_air + (mu_water - mu_air) * H ( lm );
			
			//nu = mean ( mu[k][j][i], mu[k][j+1][i] );
						
			if ( (j==0 && !j_periodic && !jj_periodic) || nvert[k][j][i]>0.1 ) nu = mu[k][j+1][i];
			else if ( (j==my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i]>0.1 ) nu = mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j+1][i] );
		}
		
		if( les || (rans && ti>0) ) {
			if ( (j==0 && !j_periodic && !jj_periodic) || nvert[k][j][i]>0.1 ) nu_t = lnu_t[k][j+1][i];
			else if( (j==my-2 && !j_periodic && !jj_periodic) || nvert[k][j+1][i]>0.1 ) nu_t = lnu_t[k][j][i];
			else nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
			
						//if (i==0 || i==mx-2) nu_t=0;
			//if (j==0 && k==25 && i==10) {
				//nu_t=nu_tau_wall;
			
			//}
			//if (k==0 || k==mz-2) nu_t=0;
			
			// wall function for boundary
			if(i!=0 && i!=mx-1 && k!=0 && k!=mz-1 && freesurface_wallmodel && level[k][j][i]>0.0 && level[k][j+1][i]<0.0) {
				double area = sqrt( eta[k][j+1][i].x*eta[k][j+1][i].x + eta[k][j+1][i].y*eta[k][j+1][i].y + eta[k][j+1][i].z*eta[k][j+1][i].z );
				double sb, sc; 
				double ni[3], nj[3], nk[3];
				Cmpnts Uc, Ua, Ub;
				Ua.x = Ua.y = Ua.z = 0;
				sb = 0.5/aj[k][j+1][i]/area;
				sc = 1.5/aj[k][j+1][i]/area;
				Ub = ucat[k][j+1][i];
				Ua.x = .5*ucat[k][j][i].x+.5*ucat[k][j+1][i].x;
				Ua.y = .5*ucat[k][j][i].y+.5*ucat[k][j+1][i].y;
				Ua.z = .5*ucat[k][j][i].z+.5*ucat[k][j+1][i].z;
				Uc= ucat[k][j+2][i];
				double nnu=mu[k][j+2][i]/rho[k][j+2][i];
				Calculate_normal(csi[k][j+1][i], eta[k][j+1][i], zet[k][j+1][i], ni, nj, nk);
				if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
				double ren = user->ren;
				double ustar1, ustar2, nu_t_1, nu_t_2;
				wall_function_freesurface(nnu, sc, Ua, Uc, &ustar2, nj[0], nj[1], nj[2]);
				wall_function_freesurface(nnu, sb, Ua, Ub, &ustar1, nj[0], nj[1], nj[2]);
				nu_t_1=ustar1*ustar1/((Ub.z-Ua.z)/sb)-nnu;
				nu_t_2=ustar2*ustar2/((Uc.z-Ua.z)/sc)-nnu;
				if (i==mx/2 && k==mz/2) printf("nu_t=%f,nu_t_2=%f,ustar=%f,ustar2=%f\n",nu_t_1,nu_t_2,ustar1,ustar2);
				nu_t=nu_t_2;
				//ustar[k][j+1][i]=ustar1;
				if(nu_t<0.0)nu_t=0.;
				//if(k>31 || i<1 || i>31 || k<1)nu_t=0.;
			}
			else if(i!=0 && i!=mx-1 && k!=0 && k!=mz-1 && j==0 && viscosity_wallmodel) {
				double area = sqrt( eta[k][j+1][i].x*eta[k][j+1][i].x + eta[k][j+1][i].y*eta[k][j+1][i].y + eta[k][j+1][i].z*eta[k][j+1][i].z );
				double sb, sc; 
				double ni[3], nj[3], nk[3];
				Cmpnts Uc, Ua, Ub;
				Ua.x = Ua.y = Ua.z = 0;
				sb = 0.5/aj[k][j+1][i]/area;
				Ub = ucat[k][j+1][i];
				Calculate_normal(csi[k][j+1][i], eta[k][j+1][i], zet[k][j+1][i], ni, nj, nk);
				if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
				double ren = user->ren;
				wall_function_freesurface(1./ren, sb, Ua, Ub, &ustar[k][j+1][i], nj[0], nj[1], nj[2]);
				nu_t=ustar[k][j+1][i]*ustar[k][j+1][i]/(Ub.z/sb)-1./ren;
				if(nu_t<0.0)nu_t=0.;
				if (i==mx/2 && k==mz/2) printf("nu_t=%f,ustar=%f \n",nu_t,ustar[k][j+1][i]);
			}					
			
			
			if(les && levelset) {
			  if(  ((fabs(lm) >= dthick) && (levelset_weno==2))  ||  (levelset_weno==0)  ||  ((levelset_weno==3)&&(lm < 0.))  ||  ((levelset_weno==4)&&(lm < -dthick))) {
					div2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu_t);
					div2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu_t);
					div2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu_t);
			  }
			}
			else {
				visc2[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu_t);
				visc2[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu_t);
				visc2[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu_t);
			}
		}
		else {
			visc2[k][j][i].x = 0;
			visc2[k][j][i].y = 0;
			visc2[k][j][i].z = 0;
		}
		
		if(laplacian) {
			r11=r21=r31=0.;
			r12=r22=r32=0.;
			r13=r23=r33=0.;
		}
		
		visc2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu);
		visc2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu);
		visc2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu);
		
		if(clark) {
			double dx, dy, dz;
			Calculate_dxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
			double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
			
			double t11 = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 );
			double t12 = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 );
			double t13 = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 );
			double t21 = t12;
			double t22 = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 );
			double t23 = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 );
			double t31 = t13;
			double t32 = t23;
			double t33 = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 );
			
			visc2[k][j][i].x -= ( t11 * eta0 + t12 * eta1 + t13 * eta2 ) / 12.;
			visc2[k][j][i].y -= ( t21 * eta0 + t22 * eta1 + t23 * eta2 ) / 12.;
			visc2[k][j][i].z -= ( t31 * eta0 + t32 * eta1 + t33 * eta2 ) / 12.;
		}
		/*
		if( nvert[k][j][i]+nvert[k][j+1][i]>1.1 ) {
			Set (&div2[k][j][i], 0.);
			Set (&visc2[k][j][i], 0.);
			}*/
	}
  
	// k direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(i==0 || j==0) continue;
		
		PetscReal ajc = kaj[k][j][i];
		csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		
		Compute_du_k (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		
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

		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		int kL=k-1, kR=k+2;
		double denom=3.;
		double lm=1000.;
		
		if(levelset) lm = 0.5 * ( level[k][j][i] + level[k+1][j][i] );
		
		if(k==0 || k==mz-2) {
			if(k_periodic) kL = mz-3, kR = 2;
			else if(kk_periodic && k==mz-2) kR=mz+2;
			else if(kk_periodic && k==0) kL=-3;
			else kL = k, kR=k+1, denom=1.;
		}
		else if( nvert[kL][j][i]+nvert[kR][j][i] > 0.1 ) kL = k, kR=k+1, denom=1.;
		
		//#ifdef SECOND_ORDER
		if(second_order) {
		  kL = k, kR=k+1;
		  denom=1.;
		}
		//#endif
		
		double ucon = ucont[k][j][i].z;
		
		if(k_periodic && (k==0 || k==mz-2) )  ucon = ucont[mz-2][j][i].z;
		else if(kk_periodic && k==0)  ucon = ucont[-2][j][i].z;
		
		if(k==mz-2 && user->bctype[5]==4 && (int)nvert[k][j][i]==0) {
			ucon = ucont[k-1][j][i].z;
		}

		double up = - 0.5 * ( ucon + fabs(ucon) );
		double um = - 0.5 * ( ucon - fabs(ucon) );

		if ( immersed && k!=mz-2 && nvert[k][j][i]>0.1 ) {
			div3[k][j][i].x = 
				um * (0.125 * (-    ucat[k+2][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
				up * (0.125 * (-    ucat[k  ][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);  		     
			div3[k][j][i].y = 	       			     
				um * (0.125 * (-    ucat[k+2][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) + ucat[k+1][j][i].y)   +		     
				up * (0.125 * (-    ucat[k  ][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) +  ucat[k][j][i  ].y);  		     
			div3[k][j][i].z = 	       			     
				um * (0.125 * (-    ucat[k+2][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
				up * (0.125 * (-    ucat[k  ][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
		}
		else if ( immersed &&  k!=0 && nvert[k+1][j][i]>0.1 ) {
			div3[k][j][i].x = 
				um * (0.125 * (-    ucat[k+1][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
				up * (0.125 * (-    ucat[k-1][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);  		     
			div3[k][j][i].y = 	       			     
				um * (0.125 * (-    ucat[k+1][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) +  ucat[k+1][j][i].y)   +		     
				up * (0.125 * (-    ucat[k-1][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) + ucat[k][j][i  ].y);  		     
			div3[k][j][i].z = 	       			     
				um * (0.125 * (-    ucat[k+1][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
				up * (0.125 * (-    ucat[k-1][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
		}
		else {
			if( inviscid || (rans && !central) || (les && (   (levelset_weno==1) ||  (levelset_weno==2 && (fabs(lm) <= dthick))   ||  (levelset_weno==3 && (lm > 0.))  ||   (levelset_weno==4 && (lm > -dthick))   ))  || (levelset_weno==5) ) {
		    if(1==0) {
				div3[k][j][i].x = 
					um * (0.125 * (-    ucat[kR][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +		     
					up * (0.125 * (-    ucat[kL][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) + ucat[k][j][i  ].x);  		     
				div3[k][j][i].y = 	       			     
					um * (0.125 * (-    ucat[kR][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) + ucat[k+1][j][i].y)   +		     
					up * (0.125 * (-    ucat[kL][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) + ucat[k][j][i  ].y);  		     
				div3[k][j][i].z = 	       			     
				  um * (0.125 * (-    ucat[kR][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +		     
					up * (0.125 * (-    ucat[kL][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
		    }
				
		    else {
				div3[k][j][i].x = - ucont[k][j][i].z * weno3(ucat[kL][j][i].x, ucat[k][j][i].x, ucat[k+1][j][i].x, ucat[kR][j][i].x, ucont[k][j][i].z);
				div3[k][j][i].y = - ucont[k][j][i].z * weno3(ucat[kL][j][i].y, ucat[k][j][i].y, ucat[k+1][j][i].y, ucat[kR][j][i].y, ucont[k][j][i].z);
				div3[k][j][i].z = - ucont[k][j][i].z * weno3(ucat[kL][j][i].z, ucat[k][j][i].z, ucat[k+1][j][i].z, ucat[kR][j][i].z, ucont[k][j][i].z);
		    }
		  }
			else {
			  //#ifdef SECOND_ORDER
			  if(second_order) {
				div3[k][j][i].x = - ucont[k][j][i].z * 0.5 * (ucat[k][j][i].x + ucat[k+1][j][i].x);
				div3[k][j][i].y = - ucont[k][j][i].z * 0.5 * (ucat[k][j][i].y + ucat[k+1][j][i].y);
				div3[k][j][i].z = - ucont[k][j][i].z * 0.5 * (ucat[k][j][i].z + ucat[k+1][j][i].z);
			  }				
			  //#else
			  else {
				div3[k][j][i].x = -ucont[k][j][i].z * 0.0625 * ( -ucat[kL][j][i].x + 9.*ucat[k][j][i].x + 9.*ucat[k+1][j][i].x - ucat[kR][j][i].x );
				div3[k][j][i].y = -ucont[k][j][i].z * 0.0625 * ( -ucat[kL][j][i].y + 9.*ucat[k][j][i].y + 9.*ucat[k+1][j][i].y - ucat[kR][j][i].y );
				div3[k][j][i].z = -ucont[k][j][i].z * 0.0625 * ( -ucat[kL][j][i].z  + 9.*ucat[k][j][i].z + 9.*ucat[k+1][j][i].z - ucat[kR][j][i].z );
				//#endif
			  }
			}
		}		
		
		if(skew) {
			div3[k][j][i].x *= 0.5;
			div3[k][j][i].y *= 0.5;
			div3[k][j][i].z *= 0.5;
			
			Cmpnts ducat_dz4;
			Subtract_Scale_Set(ucat[kR][j][i], ucat[kL][j][i], 1./denom, &ducat_dz4);
			
			adv3[k][j][i].x = - 0.5 * ucont[k][j][i].z * ( 9./8. * dudz - 1./8. * ducat_dz4.x );
			adv3[k][j][i].y = - 0.5 * ucont[k][j][i].z * ( 9./8. * dvdz - 1./8. * ducat_dz4.y );
			adv3[k][j][i].z = - 0.5 * ucont[k][j][i].z * ( 9./8. * dwdz - 1./8. * ducat_dz4.z );
		}
		
		
		if( nvert[k][j][i]+nvert[k+1][j][i]>0.1) {
			if(immersed==3 || !immersed) {
				div3[k][j][i].x = 0;
				div3[k][j][i].y = 0;
				div3[k][j][i].z = 0;
				if(skew) {
				  adv3[k][j][i].x = 0;
				  adv3[k][j][i].y = 0;
				  adv3[k][j][i].z = 0;
				}
			}
			/*
			if(skew) {
				adv3[k][j][i].x = 0;
				adv3[k][j][i].y = 0;
				adv3[k][j][i].z = 0;
			}
			*/
		}
		/*
		div3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu+nu_t);
		div3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu+nu_t);
		div3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu+nu_t);
		*/
		
		double nu = 1./user->ren, nu_t =0;
		
		if(levelset) {
			//double lm = 0.5 * ( level[k][j][i] + level[k+1][j][i] );
			//nu = mu_air + (mu_water - mu_air) * H ( lm );
			
			//nu = mean ( mu[k][j][i], mu[k+1][j][i] );
			
			if ( (k==0 && !k_periodic && !kk_periodic) || nvert[k][j][i]>0.1 ) nu = mu[k+1][j][i];
			else if ( (k==mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i]>0.1 ) nu = mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k+1][j][i] );
		}
		
		if( les || (rans && ti>0) ) {
			if ( (k==0 && !k_periodic && !kk_periodic) || nvert[k][j][i]>0.1 ) nu_t = lnu_t[k+1][j][i];
			else if( (k==mz-2 && !k_periodic && !kk_periodic) || nvert[k+1][j][i]>0.1 ) nu_t = lnu_t[k][j][i];
			else nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
			
						//if (i==0 || i==mx-2) nu_t=0;
			//if (j==0 || j==my-2) nu_t=nu_tau_wall;
			//if (k==0 || k==mz-2) nu_t=0;
			
			if(les && levelset) {
			  if(  ((fabs(lm) >= dthick) && (levelset_weno==2))  ||  (levelset_weno==0)  ||  ((levelset_weno==3)&&(lm < 0.))  ||  ((levelset_weno==4)&&(lm < -dthick))) {
					div3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu_t);
					div3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu_t);
					div3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu_t);
			  }
			}
			else {
				visc3[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu_t);
				visc3[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu_t);
				visc3[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu_t);
			}
		}
		else {
			visc3[k][j][i].x = 0;
			visc3[k][j][i].y = 0;
			visc3[k][j][i].z = 0;
		}
		
		if(laplacian) {
			r11=r21=r31=0.;
			r12=r22=r32=0.;
			r13=r23=r33=0.;
		}
		
		visc3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu);
		visc3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu);
		visc3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu);
		
		if(clark) {
			double dx, dy, dz;
			Calculate_dxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
			double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
			
			double t11 = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 );
			double t12 = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 );
			double t13 = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 );
			double t21 = t12;
			double t22 = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 );
			double t23 = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 );
			double t31 = t13;
			double t32 = t23;
			double t33 = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 );
			
			visc3[k][j][i].x -= ( t11 * zet0 + t12 * zet1 + t13 * zet2 ) / 12.;
			visc3[k][j][i].y -= ( t21 * zet0 + t22 * zet1 + t23 * zet2 ) / 12.;
			visc3[k][j][i].z -= ( t31 * zet0 + t32 * zet1 + t33 * zet2 ) / 12.;
		}
		/*
		if( nvert[k][j][i]+nvert[k+1][j][i]>1.1 ) {
			Set (&div3[k][j][i], 0.);
			Set (&visc3[k][j][i], 0.);
			}*/
	}
	
	
		DAVecRestoreArray(fda, Div1, &div1);
		DAVecRestoreArray(fda, Div2, &div2);
		DAVecRestoreArray(fda, Div3, &div3);
		
		DALocalToLocalBegin(fda, Div1, INSERT_VALUES, Div1);
		DALocalToLocalEnd(fda, Div1, INSERT_VALUES, Div1);
		DALocalToLocalBegin(fda, Div2, INSERT_VALUES, Div2);
		DALocalToLocalEnd(fda, Div2, INSERT_VALUES, Div2);
		DALocalToLocalBegin(fda, Div3, INSERT_VALUES, Div3);
		DALocalToLocalEnd(fda, Div3, INSERT_VALUES, Div3);
		
		DAVecGetArray(fda, Div1, &div1);
		DAVecGetArray(fda, Div2, &div2);
		DAVecGetArray(fda, Div3, &div3);
	
	
	if(skew) {
		DAVecRestoreArray(fda, Adv1, &adv1);
		DAVecRestoreArray(fda, Adv2, &adv2);
		DAVecRestoreArray(fda, Adv3, &adv3);
		
		DALocalToLocalBegin(fda, Adv1, INSERT_VALUES, Adv1);
		DALocalToLocalEnd(fda, Adv1, INSERT_VALUES, Adv1);
		DALocalToLocalBegin(fda, Adv2, INSERT_VALUES, Adv2);
		DALocalToLocalEnd(fda, Adv2, INSERT_VALUES, Adv2);
		DALocalToLocalBegin(fda, Adv3, INSERT_VALUES, Adv3);
		DALocalToLocalEnd(fda, Adv3, INSERT_VALUES, Adv3);
		
		DAVecGetArray(fda, Adv1, &adv1);
		DAVecGetArray(fda, Adv2, &adv2);
		DAVecGetArray(fda, Adv3, &adv3);
	}
	
	DAVecRestoreArray(fda, Visc1, &visc1);
	DAVecRestoreArray(fda, Visc2, &visc2);
	DAVecRestoreArray(fda, Visc3, &visc3);
	
	DALocalToLocalBegin(fda, Visc1, INSERT_VALUES, Visc1);
	DALocalToLocalEnd(fda, Visc1, INSERT_VALUES, Visc1);
	DALocalToLocalBegin(fda, Visc2, INSERT_VALUES, Visc2);
	DALocalToLocalEnd(fda, Visc2, INSERT_VALUES, Visc2);
	DALocalToLocalBegin(fda, Visc3, INSERT_VALUES, Visc3);
	DALocalToLocalEnd(fda, Visc3, INSERT_VALUES, Visc3);
	
	DAVecGetArray(fda, Visc1, &visc1);
	DAVecGetArray(fda, Visc2, &visc2);
	DAVecGetArray(fda, Visc3, &visc3);

	
	
	DAVecGetArray(fda, Fp, &fp);
	
	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int a=i, b=j, c=k;

		int flag=0;
		
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
			div1[k][j][i] = div1[c][b][a];
			div2[k][j][i] = div2[c][b][a];
			div3[k][j][i] = div3[c][b][a];
			visc1[k][j][i] = visc1[c][b][a];
			visc2[k][j][i] = visc2[c][b][a];
			visc3[k][j][i] = visc3[c][b][a];
			
			if(skew) {
				adv1[k][j][i] = adv1[c][b][a];
				adv2[k][j][i] = adv2[c][b][a];
				adv3[k][j][i] = adv3[c][b][a];
			}
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		Cmpnts div, div_4th;
		
		div.x = (div1[k][j][i].x - div1[k][j][i-1].x + div2[k][j][i].x - div2[k][j-1][i].x + div3[k][j][i].x - div3[k-1][j][i].x);
		div.y = (div1[k][j][i].y - div1[k][j][i-1].y + div2[k][j][i].y - div2[k][j-1][i].y + div3[k][j][i].y - div3[k-1][j][i].y);
		div.z = (div1[k][j][i].z - div1[k][j][i-1].z + div2[k][j][i].z - div2[k][j-1][i].z + div3[k][j][i].z - div3[k-1][j][i].z);
		
		div_4th.x = div_4th.y = div_4th.z = 0;
		 
		#ifdef SECOND_ORDER
		
		fp[k][j][i] = div;
		
		#else
		
		int iR=i+1, iL=i-2;
		int jR=j+1, jL=j-2;
		int kR=k+1, kL=k-2;
		double denom_i=3.;
		double denom_j=3.;
		double denom_k=3.;
		
		if (i==1) {
			if( i_periodic ) iL=mx-3;
			else if( ii_periodic ) iL=-3;
			else iR=i, iL=i-1, denom_i=1.;
		}
		else if(i==2 || i==mx-3) {
			if( i_periodic ) {}
			else if( ii_periodic ) {}
			else iR=i, iL=i-1, denom_i=1.;
		}
		else if(i==mx-2) {
			if( i_periodic) iR=1;
			else if( ii_periodic) iR=mx+1;
			else iR=i, iL=i-1, denom_i=1.;
		}
		
		if (j==1) {
			if( j_periodic ) jL=my-3;
			else if( jj_periodic ) jL=-3;
			else jR=j, jL=j-1, denom_j=1.;
		}
		else if(j==2 || j==my-3) {
			if( j_periodic ) {}
			else if( jj_periodic ) {}
			else jR=j, jL=j-1, denom_j=1.;
		}
		else if(j==my-2) {
			if( j_periodic) jR=1;
			else if( jj_periodic) jR=my+1;
			else jR=j, jL=j-1, denom_j=1.;
		}
		
		if (k==1) {
			if( k_periodic ) kL=mz-3;
			else if( kk_periodic ) kL=-3;
			else kR=k, kL=k-1, denom_k=1.;
		}
		else if(k==2 || k==mz-3) {
			if( k_periodic ) {}
			else if( kk_periodic ) {}
			else kR=k, kL=k-1, denom_k=1.;
		}
		else if(k==mz-2) {
			if( k_periodic) kR=1;
			else if( kk_periodic) kR=mz+1;
			else kR=k, kL=k-1, denom_k=1.;
		}
		
		if (nvert[k][j][i-1]+nvert[k][j][i]+nvert[k][j][i+1]>0.1) iR=i, iL=i-1, denom_i=1.;
		if (nvert[k][j-1][i]+nvert[k][j][i]+nvert[k][j+1][i]>0.1) jR=j, jL=j-1, denom_j=1.;
		if (nvert[k-1][j][i]+nvert[k][j][i]+nvert[k+1][j][i]>0.1) kR=k, kL=k-1, denom_k=1.;
			
		if(inviscid || (rans && !central)) fp[k][j][i] = div;
		else {
		  if(second_order) {
		    fp[k][j][i] = div;
		    iR=i, iL=i-1;
		    jR=j, jL=j-1;
		    kR=k, kL=k-1;
		  }
		  else {
			Subtract_Scale_AddTo(div1[k][j][iR], div1[k][j][iL], 1./denom_i, &div_4th); 
			Subtract_Scale_AddTo(div2[k][jR][i], div2[k][jL][i], 1./denom_j, &div_4th); 
			Subtract_Scale_AddTo(div3[kR][j][i], div3[kL][j][i], 1./denom_k, &div_4th); 
			AxByC ( 9./8., div, -1./8., div_4th, &fp[k][j][i]);
		  }
			if(skew) {
				fp[k][j][i].x += 9./8. * 0.5 * (adv1[k][j][i].x+adv1[k][j][i-1].x) - 1./8. * 0.5 * ( adv1[k][j][iR].x + adv1[k][j][iL].x );
				fp[k][j][i].x += 9./8. * 0.5 * (adv2[k][j][i].x+adv2[k][j-1][i].x) - 1./8. * 0.5 * ( adv2[k][jR][i].x + adv2[k][jL][i].x );
				fp[k][j][i].x += 9./8. * 0.5 * (adv3[k][j][i].x+adv3[k-1][j][i].x) - 1./8. * 0.5 * ( adv3[kR][j][i].x + adv3[kL][j][i].x );
				
				fp[k][j][i].y += 9./8. * 0.5 * (adv1[k][j][i].y+adv1[k][j][i-1].y) - 1./8. * 0.5 * ( adv1[k][j][iR].y + adv1[k][j][iL].y );
				fp[k][j][i].y += 9./8. * 0.5 * (adv2[k][j][i].y+adv2[k][j-1][i].y) - 1./8. * 0.5 * ( adv2[k][jR][i].y + adv2[k][jL][i].y );
				fp[k][j][i].y += 9./8. * 0.5 * (adv3[k][j][i].y+adv3[k-1][j][i].y) - 1./8. * 0.5 * ( adv3[kR][j][i].y + adv3[kL][j][i].y );
				
				fp[k][j][i].z += 9./8. * 0.5 * (adv1[k][j][i].z+adv1[k][j][i-1].z) - 1./8. * 0.5 * ( adv1[k][j][iR].z + adv1[k][j][iL].z );
				fp[k][j][i].z += 9./8. * 0.5 * (adv2[k][j][i].z+adv2[k][j-1][i].z) - 1./8. * 0.5 * ( adv2[k][jR][i].z + adv2[k][jL][i].z );
				fp[k][j][i].z += 9./8. * 0.5 * (adv3[k][j][i].z+adv3[k-1][j][i].z) - 1./8. * 0.5 * ( adv3[kR][j][i].z + adv3[kL][j][i].z );
			}
		}
		#endif

		double r=1.0;
		
		if(levelset) {
			double lm = level[k][j][i];
			double dx=pow(1./aj[k][j][i],1./3.);
			//r = rho_air + (rho_water - rho_air) * H ( lm, dx );//100305
			r = rho[k][j][i];
		}
		
		if ( inviscid ) {}
		//else if ( user->bctype[4]==5 && levelset && k<=inlet_buffer_k ) {} // neglect gas viscosity
		else {
			fp[k][j][i].x += (visc1[k][j][i].x - visc1[k][j][i-1].x + visc2[k][j][i].x - visc2[k][j-1][i].x + visc3[k][j][i].x - visc3[k-1][j][i].x) / r;
			fp[k][j][i].y += (visc1[k][j][i].y - visc1[k][j][i-1].y + visc2[k][j][i].y - visc2[k][j-1][i].y + visc3[k][j][i].y - visc3[k-1][j][i].y) / r;
			fp[k][j][i].z += (visc1[k][j][i].z - visc1[k][j][i-1].z + visc2[k][j][i].z - visc2[k][j-1][i].z + visc3[k][j][i].z - visc3[k-1][j][i].z) / r;
		}
		
		/*
		fp[k][j][i].x *= aj[k][j][i];
		fp[k][j][i].y *= aj[k][j][i];
		fp[k][j][i].z *= aj[k][j][i];
		*/
		//
		
	}
	
	DAVecRestoreArray(fda, Fp, &fp);
	
	DALocalToLocalBegin(fda, Fp, INSERT_VALUES, Fp);
	DALocalToLocalEnd(fda, Fp, INSERT_VALUES, Fp);
	
	DAVecGetArray(fda, Fp, &fp);
	
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
		
		if(flag) fp[k][j][i] = fp[c][b][a];

	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
				
		/*#ifdef SECOND_ORDER
		rhs[k][j][i].x += 	( icsi[k][j][i].x * ( fp[k][j][i].x/aj[k][j][i+1] + fp[k][j][i+1].x/aj[k][j][i] ) +
						icsi[k][j][i].y * ( fp[k][j][i].y/aj[k][j][i+1] + fp[k][j][i+1].y/aj[k][j][i] ) +
						icsi[k][j][i].z * ( fp[k][j][i].z/aj[k][j][i+1] + fp[k][j][i+1].z/aj[k][j][i] ) ) / (1./aj[k][j][i] + 1./aj[k][j][i+1]) * iaj[k][j][i];
		
		rhs[k][j][i].y +=	( jeta[k][j][i].x * ( fp[k][j][i].x/aj[k][j+1][i] + fp[k][j+1][i].x/aj[k][j][i] ) +
						jeta[k][j][i].y * ( fp[k][j][i].y/aj[k][j+1][i] + fp[k][j+1][i].y/aj[k][j][i] ) +
						jeta[k][j][i].z * ( fp[k][j][i].z/aj[k][j+1][i] + fp[k][j+1][i].z/aj[k][j][i] ) ) / (1./aj[k][j][i] + 1./aj[k][j+1][i]) * jaj[k][j][i];
		
		rhs[k][j][i].z +=	( kzet[k][j][i].x * ( fp[k][j][i].x/aj[k+1][j][i] + fp[k+1][j][i].x/aj[k][j][i] ) +
						kzet[k][j][i].y * ( fp[k][j][i].y/aj[k+1][j][i] + fp[k+1][j][i].y/aj[k][j][i] ) +
						kzet[k][j][i].z * ( fp[k][j][i].z/aj[k+1][j][i] + fp[k+1][j][i].z/aj[k][j][i] ) ) / (1./aj[k][j][i] + 1./aj[k+1][j][i]) * kaj[k][j][i];
		#else*/
		
		rhs[k][j][i].x += scale * ( 0.5 * ( csi[k][j][i].x * fp[k][j][i].x + csi[k][j][i].y * fp[k][j][i].y + csi[k][j][i].z * fp[k][j][i].z) + 0.5 * ( csi[k][j][i+1].x * fp[k][j][i+1].x + csi[k][j][i+1].y * fp[k][j][i+1].y + csi[k][j][i+1].z * fp[k][j][i+1].z) ) * iaj[k][j][i];
		rhs[k][j][i].y += scale * ( 0.5 * ( eta[k][j][i].x * fp[k][j][i].x + eta[k][j][i].y * fp[k][j][i].y + eta[k][j][i].z * fp[k][j][i].z) + 0.5 * ( eta[k][j+1][i].x * fp[k][j+1][i].x + eta[k][j+1][i].y * fp[k][j+1][i].y + eta[k][j+1][i].z * fp[k][j+1][i].z) ) * jaj[k][j][i];
		rhs[k][j][i].z += scale * ( 0.5 * ( zet[k][j][i].x * fp[k][j][i].x + zet[k][j][i].y * fp[k][j][i].y + zet[k][j][i].z * fp[k][j][i].z) + 0.5 * ( zet[k+1][j][i].x * fp[k+1][j][i].x + zet[k+1][j][i].y * fp[k+1][j][i].y + zet[k+1][j][i].z * fp[k+1][j][i].z) ) * kaj[k][j][i];
		
		/*
		rhs[k][j][i].x +=	0.5 * ( icsi[k][j][i].x * ( fp[k][j][i].x + fp[k][j][i+1].x ) + icsi[k][j][i].y * ( fp[k][j][i].y + fp[k][j][i+1].y ) + icsi[k][j][i].z * ( fp[k][j][i].z + fp[k][j][i+1].z ) ) * iaj[k][j][i];
		rhs[k][j][i].y +=	0.5 * ( jeta[k][j][i].x * ( fp[k][j][i].x + fp[k][j+1][i].x ) + jeta[k][j][i].y * ( fp[k][j][i].y + fp[k][j+1][i].y ) +jeta[k][j][i].z * ( fp[k][j][i].z + fp[k][j+1][i].z ) ) * jaj[k][j][i];
		rhs[k][j][i].z +=	0.5 * ( kzet[k][j][i].x * ( fp[k][j][i].x + fp[k+1][j][i].x ) + kzet[k][j][i].y * ( fp[k][j][i].y + fp[k+1][j][i].y ) + kzet[k][j][i].z * ( fp[k][j][i].z + fp[k+1][j][i].z ) ) * kaj[k][j][i];
		*/
		/*#endif*/
		
		if(levelset) {
			/*double dH1 = dH ( 0.5 * ( level[k][j][i] + level[k][j][i+1] ) );
			double dH2 = dH ( 0.5 * ( level[k][j][i] + level[k][j+1][i] ) );
			double dH3 = dH ( 0.5 * ( level[k][j][i] + level[k+1][j][i] ) );
*/
			double r1 = mean ( rho[k][j][i], rho[k][j][i+1] );
			double r2 = mean ( rho[k][j][i], rho[k][j+1][i] );
			double r3 = mean ( rho[k][j][i], rho[k+1][j][i] );
			
			//double gx = 0;
			//double gy = 0;
			//double gz = -1./(Fr*Fr);
			
			
			  
			rhs[k][j][i].x += scale * ( gravity_x * icsi[k][j][i].x + gravity_y * icsi[k][j][i].y + gravity_z * icsi[k][j][i].z );
			rhs[k][j][i].y += scale * ( gravity_x * jeta[k][j][i].x + gravity_y * jeta[k][j][i].y + gravity_z * jeta[k][j][i].z );
			rhs[k][j][i].z += scale * ( gravity_x * kzet[k][j][i].x + gravity_y * kzet[k][j][i].y + gravity_z * kzet[k][j][i].z );
			  
// add (Toni)
	// For wave_momentum_source
			if (wave_momentum_source){
				double sx, sy, sz;
				double al=0.5, ar=0.5;				
				sx = (al*WAVE_fp[k][j][i].x+ar*WAVE_fp[k][j][i+1].x);
				sy = (al*WAVE_fp[k][j][i].y+ar*WAVE_fp[k][j][i+1].y);
				sz = (al*WAVE_fp[k][j][i].z+ar*WAVE_fp[k][j][i+1].z);		
				rhs[k][j][i].x += scale * ( sx * icsi[k][j][i].x + sy * icsi[k][j][i].y + sz * icsi[k][j][i].z );	
				sx = (al*WAVE_fp[k][j][i].x+ar*WAVE_fp[k][j+1][i].x);
				sy = (al*WAVE_fp[k][j][i].y+ar*WAVE_fp[k][j+1][i].y);
				sz = (al*WAVE_fp[k][j][i].z+ar*WAVE_fp[k][j+1][i].z);
				rhs[k][j][i].y += scale * ( sx * jeta[k][j][i].x + sy * jeta[k][j][i].y + sz * jeta[k][j][i].z );
				sx = (al*WAVE_fp[k][j][i].x+ar*WAVE_fp[k+1][j][i].x);
				sy = (al*WAVE_fp[k][j][i].y+ar*WAVE_fp[k+1][j][i].y);
				sz = (al*WAVE_fp[k][j][i].z+ar*WAVE_fp[k+1][j][i].z);
				rhs[k][j][i].z += scale * ( sx * kzet[k][j][i].x + sy * kzet[k][j][i].y + sz * kzet[k][j][i].z );
			}					  
	
// end (Toni) 
			if(surface_tension ) {
			  
				double sx, sy, sz;
				double al=0.5, ar=0.5;
				
				if ( nvert[k][j][i]>0.1 ) al=0., ar=1.;
				else if ( nvert[k][j][i+1]>0.1 ) al=1., ar=0.;
				else al=0.5, ar=0.5;
						
				sx = (al*stension[k][j][i].x+ar*stension[k][j][i+1].x);
				sy = (al*stension[k][j][i].y+ar*stension[k][j][i+1].y);
				sz = (al*stension[k][j][i].z+ar*stension[k][j][i+1].z);
					
				rhs[k][j][i].x += scale * ( sx * icsi[k][j][i].x + sy * icsi[k][j][i].y + sz * icsi[k][j][i].z );
					
				if ( nvert[k][j][i]>0.1 ) al=0., ar=1.;
				else if ( nvert[k][j+1][i]>0.1 ) al=1., ar=0.;
				else al=0.5, ar=0.5;
					
				sx = (al*stension[k][j][i].x+ar*stension[k][j+1][i].x);
				sy = (al*stension[k][j][i].y+ar*stension[k][j+1][i].y);
				sz = (al*stension[k][j][i].z+ar*stension[k][j+1][i].z);
					
				rhs[k][j][i].y += scale * ( sx * jeta[k][j][i].x + sy * jeta[k][j][i].y + sz * jeta[k][j][i].z );
					
				if ( nvert[k][j][i]>0.1 ) al=0., ar=1.;
				else if ( nvert[k+1][j][i]>0.1 ) al=1., ar=0.;
				else al=0.5, ar=0.5;
					
				sx = (al*stension[k][j][i].x+ar*stension[k+1][j][i].x);
				sy = (al*stension[k][j][i].y+ar*stension[k+1][j][i].y);
				sz = (al*stension[k][j][i].z+ar*stension[k+1][j][i].z);
					
				rhs[k][j][i].z += scale * ( sx * kzet[k][j][i].x + sy * kzet[k][j][i].y + sz * kzet[k][j][i].z );
				
				/*
				rhs[k][j][i].x += scale * ( 0.5 * (stension[k][j][i].x+stension[k][j][i+1].x) * icsi[k][j][i].x + 0.5 * (stension[k][j][i].y+stension[k][j][i+1].y) * icsi[k][j][i].y + 0.5 * (stension[k][j][i].z+stension[k][j][i+1].z) * icsi[k][j][i].z );
				rhs[k][j][i].y += scale * ( 0.5 * (stension[k][j][i].x+stension[k][j+1][i].x) * jeta[k][j][i].x + 0.5 * (stension[k][j][i].y+stension[k][j+1][i].y) * jeta[k][j][i].y + 0.5 * (stension[k][j][i].z+stension[k][j+1][i].z) * jeta[k][j][i].z );
				rhs[k][j][i].z += scale * ( 0.5 * (stension[k][j][i].x+stension[k+1][j][i].x) * kzet[k][j][i].x + 0.5 * (stension[k][j][i].y+stension[k+1][j][i].y) * kzet[k][j][i].y + 0.5 * (stension[k][j][i].z+stension[k+1][j][i].z) * kzet[k][j][i].z );
				*/
					/*
				rhs[k][j][i].x += scale * ( 0.5 * ( csi[k][j][i].x * stension[k][j][i].x + csi[k][j][i].y * stension[k][j][i].y + csi[k][j][i].z * stension[k][j][i].z) + 0.5 * ( csi[k][j][i+1].x * stension[k][j][i+1].x + csi[k][j][i+1].y * stension[k][j][i+1].y + csi[k][j][i+1].z * stension[k][j][i+1].z) );// * dH1;
				rhs[k][j][i].y += scale * ( 0.5 * ( eta[k][j][i].x * stension[k][j][i].x + eta[k][j][i].y * stension[k][j][i].y + eta[k][j][i].z * stension[k][j][i].z) + 0.5 * ( eta[k][j+1][i].x * stension[k][j+1][i].x + eta[k][j+1][i].y * stension[k][j+1][i].y + eta[k][j+1][i].z * stension[k][j+1][i].z) );// * dH2;
				rhs[k][j][i].z += scale * ( 0.5 * ( zet[k][j][i].x * stension[k][j][i].x + zet[k][j][i].y * stension[k][j][i].y + zet[k][j][i].z * stension[k][j][i].z) + 0.5 * ( zet[k+1][j][i].x * stension[k+1][j][i].x + zet[k+1][j][i].y * stension[k+1][j][i].y + zet[k+1][j][i].z * stension[k+1][j][i].z) );// * dH3;
					*/
				}
			
		}
		
		
		if(nvert[k][j][i]+nvert[k][j][i+1]>0.1 || (!i_periodic && !ii_periodic && i==mx-2) ) {
			rhs[k][j][i].x = 0;
		}
		if(nvert[k][j][i]+nvert[k][j+1][i]>0.1 || (!j_periodic && !jj_periodic && j==my-2) ) {
			rhs[k][j][i].y = 0;
		}
		if(nvert[k][j][i]+nvert[k+1][j][i]>0.1 || (!k_periodic && !kk_periodic && k==mz-2) ) {
			rhs[k][j][i].z = 0;
		}
		/*
		if ( (user->bctype[0]==-1 && i==1) ||( user->bctype[1]==-1 && i==mx-2) ) {
		  rhs[k][j][i].y = 0;
		  rhs[k][j][i].z = 0;
		}
		if ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 && j==my-2) ) {
		  rhs[k][j][i].x = 0;
                  rhs[k][j][i].z = 0;
		  }*/
	}

	if(les) {
		DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	}
	else if (rans) {
		DAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
		//DAVecRestoreArray(user->da, user->Distance, &distance);
		DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	}
	
        DAVecRestoreArray(da, user->lIAj, &iaj);
        DAVecRestoreArray(da, user->lJAj, &jaj);
        DAVecRestoreArray(da, user->lKAj, &kaj);

	if (xs ==0) {
		i = 0;
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++) {
			rhs[k][j][i].x = 0;
			rhs[k][j][i].y = 0;
			rhs[k][j][i].z = 0;
		}
	}

	if (xe == mx) {
		for (k=zs; k<ze; k++)
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


	if (ys == 0) {
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			j=0;
			rhs[k][j][i].x = 0;
			rhs[k][j][i].y = 0;
			rhs[k][j][i].z = 0;
		}
	}
  
	if (ye == my) {
		for (k=zs; k<ze; k++) 
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
	
	
	if (zs == 0) {
		k=0;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			rhs[k][j][i].x = 0;
			rhs[k][j][i].y = 0;
			rhs[k][j][i].z = 0;
		}
	}
  
	if (ze == mz) {
		for (j=ys; j<ye; j++)
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
	
	DAVecRestoreArray(fda, user->lCent, &cent);
	DAVecRestoreArray(fda, Ucont, &ucont);
	DAVecRestoreArray(fda, user->lUcont_o, &ucont_o);
	DAVecRestoreArray(fda, Ucat,  &ucat);
	DAVecRestoreArray(fda, Rhs,  &rhs);
	//DAVecRestoreArray(fda, user->RHS_o,  &rhs_o);
	//DAVecRestoreArray(fda, user->RHS_rm1,  &rhs_rm1);

	DAVecRestoreArray(fda, Csi, &csi);
	DAVecRestoreArray(fda, Eta, &eta);
	DAVecRestoreArray(fda, Zet, &zet);
	  
	DAVecRestoreArray(fda, Fp, &fp);
	DAVecRestoreArray(fda, Div1, &div1);
	DAVecRestoreArray(fda, Div2, &div2);
	DAVecRestoreArray(fda, Div3, &div3);
	
	DAVecRestoreArray(fda, Visc1, &visc1);
	DAVecRestoreArray(fda, Visc2, &visc2);
	DAVecRestoreArray(fda, Visc3, &visc3);
	
	DAVecRestoreArray(da, user->lAj, &aj);
	//DAVecRestoreArray(da, user->lVolume, &vol);
  
	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);

	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);

	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);
	
// add (Toni)
	//for wave_momentum_source
	if(wave_momentum_source) DAVecRestoreArray(fda, user->lWAVE_fp, &WAVE_fp);	 
	DAVecRestoreArray(da, user->lUstar, &ustar);	
// end (Toni)  
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lP, &p);
	
	if(levelset) {
	  DAVecRestoreArray(fda, user->lST, &stension);
		DAVecRestoreArray(user->da, user->lLevelset, &level);
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}

	/*
	VecDestroy(Fp);
	VecDestroy(Div1);
	VecDestroy(Div2);
	VecDestroy(Div3);
	VecDestroy(Visc1);
	VecDestroy(Visc2);
	VecDestroy(Visc3);
	*/

	if(skew) {
		DAVecRestoreArray(fda, Adv1, &adv1);
		DAVecRestoreArray(fda, Adv2, &adv2);
		DAVecRestoreArray(fda, Adv3, &adv3);
		
		VecDestroy(Adv1);
		VecDestroy(Adv2);
		VecDestroy(Adv3);
	}
	
	
	return(0);
};


void IB_BC(UserCtx *user)
{
	PetscInt      i, j, k;
	DALocalInfo	info ;
	PetscInt	xs, xe, ys, ye, zs, ze;
	PetscInt  	mx,my,mz;	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts        ***ucont;
	PetscReal     ***nvert;
	PetscScalar   ***level;
	
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
  
	Cmpnts ***ucat, ***icsi, ***jeta, ***kzet;
	
	if(!immersed && ti==tistart) {
          DAVecGetArray(da, user->lNvert, &nvert);
          for (k=lzs; k<lze; k++) {
            for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
                if(user->bctype[0]==-1 && i==1) nvert[k][j][i]=1;
                if(user->bctype[1]==-1 && i==mx-2) nvert[k][j][i]=1;
                if(user->bctype[2]==-1 && j==1) nvert[k][j][i]=1;
                if(user->bctype[3]==-1 && j==my-2) nvert[k][j][i]=1;
                if(user->bctype[4]==-1 && k==1) nvert[k][j][i]=1;
                if(user->bctype[5]==-1 && k==mz-2) nvert[k][j][i]=1;

                if(user->bctype[0]==-2 && i==1) nvert[k][j][i]=1;
                if(user->bctype[1]==-2 && i==mx-2) nvert[k][j][i]=1;
                if(user->bctype[2]==-2 && j==1) nvert[k][j][i]=1;
		if(user->bctype[3]==-2 && j==my-2) nvert[k][j][i]=1;
                if(user->bctype[4]==-2 && k==1) nvert[k][j][i]=1;
                if(user->bctype[5]==-2 && k==mz-2) nvert[k][j][i]=1;
              }
	    }
	  }
          DAVecRestoreArray(da, user->lNvert, &nvert);
	  DALocalToLocalBegin(da, user->lNvert, INSERT_VALUES, user->lNvert);
	  DALocalToLocalEnd(da, user->lNvert, INSERT_VALUES, user->lNvert);
	  DALocalToGlobal(da, user->lNvert, INSERT_VALUES, user->Nvert);
	}
	
	DAVecGetArray(fda, user->lUcat, &ucat);
	DAVecGetArray(fda, user->lUcont, &ucont);//
	DAVecGetArray(fda, user->lICsi, &icsi);//
	DAVecGetArray(fda, user->lJEta, &jeta);//
	DAVecGetArray(fda, user->lKZet, &kzet);//
	DAVecGetArray(da, user->lNvert, &nvert);//

	if(levelset) DAVecGetArray(da, user->lLevelset,  &level);	
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
			ucat[k][j][i] = ucat[c][b][a];
		}
	}

	double  ucx, ucy, ucz;
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  
		if(immersed) {
			double f = 1.0;

			if(immersed==3) f = 0;
			
			if ((int)nvert[k][j][i]==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
				ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z) * f;
					
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z) * f;
					
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z) * f;
			}
				
			if ((int)(nvert[k][j][i+1])==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
				ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z) * f;
			}
						
			if ((int)(nvert[k][j+1][i])==1) {
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z) * f;
			}
						
			if ((int)(nvert[k+1][j][i])==1) {
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z) * f;
			}
			
			if(!movefsi && !rotatefsi && immersed==3) {
				if ((nvert[k][j][i+1]+nvert[k][j][i])>1.1) ucont[k][j][i].x = 0;
				if ((nvert[k][j+1][i]+nvert[k][j][i])>1.1) ucont[k][j][i].y = 0;
				if ((nvert[k+1][j][i]+nvert[k][j][i])>1.1) ucont[k][j][i].z = 0;
			}
		}
	  
		// wall func
		
		if ( (user->bctype[0]==-1 && i==1) ||( user->bctype[1]==-1 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[0]==-2 && i==1) ||( user->bctype[1]==-2 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
		if ( (user->bctype[2]==-2 && j==1) || (user->bctype[3]==-2 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
	}
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if ( (user->bctype[0]==10) && i==1) ucont[k][j][0].x = 0;
		if ( (user->bctype[1]==10) && i==mx-2) ucont[k][j][mx-2].x = 0;
		if ( (user->bctype[2]==10) && j==1) ucont[k][0][i].y = 0;
		if ( (user->bctype[3]==10) && j==my-2) ucont[k][my-2][i].y = 0;
		if ( (user->bctype[4]==10) && k==1) ucont[0][j][i].z = 0;
		if ( (user->bctype[5]==10) && k==mz-2) ucont[mz-2][j][i].z = 0;
//add (Toni)
			//for air_flow_levelset
		if(air_flow_levelset && k==1 && level[k][j][i]>-dthick)ucont[0][j][i].z = 0;
		if(air_flow_levelset && k==mz-2 && level[k][j][i]>-dthick)ucont[mz-2][j][i].z = 0;
//end (Toni)		
		if ( i_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][mx-2].x;
		if ( i_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][1].x;
		if ( j_periodic && j==0 ) ucont[k][0][i].y = ucont[k][my-2][i].y;
		if ( j_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][1][i].y;
		if ( k_periodic && k==0 ) ucont[0][j][i].z = ucont[mz-2][j][i].z;
		if ( k_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[1][j][i].z;
		
		if ( ii_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][-2].x;
		if ( ii_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][mx+1].x;
		
		if ( jj_periodic && j==0 ) ucont[k][0][i].y = ucont[k][-2][i].y;
		if ( jj_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][my+1][i].y;
		
		if ( kk_periodic && k==0 ) ucont[0][j][i].z = ucont[-2][j][i].z;
		if ( kk_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[mz+1][j][i].z;
	}
	DAVecRestoreArray(fda, user->lUcat, &ucat);//
	DAVecRestoreArray(fda, user->lICsi, &icsi);//
	DAVecRestoreArray(fda, user->lJEta, &jeta);//
	DAVecRestoreArray(fda, user->lKZet, &kzet);//
	DAVecRestoreArray(fda, user->lUcont, &ucont);//
	DAVecRestoreArray(da, user->lNvert, &nvert);//
	
	if(levelset) DAVecRestoreArray(da, user->lLevelset,  &level);	
	DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	
	return;
}

PetscErrorCode FormFunction_SNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;
	VecCopy(Ucont, user->Ucont);
	
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscInt	i, j, k;
	
	Cmpnts ***ucont;
	PetscReal ***nvert;
	
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

	if (lxs==0) lxs++;
	if (lxe==mx) lxe--;
	if (lys==0) lys++;
	if (lye==my) lye--;
	if (lzs==0) lzs++;
	if (lze==mz) lze--;
	
	DAVecGetArray(user->fda, user->Ucont, &ucont);
	DAVecGetArray(user->da, user->lNvert, &nvert);
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		// noslip BC 
		if(i==0 && user->bctype[0]==1) ucont[k][j][i].x = 0;
		if(i==mx-1 && user->bctype[1]==1) ucont[k][j][i-1].x = 0;
		if(j==0 && user->bctype[2]==1) ucont[k][j][i].y = 0;
		if(j==my-1 && user->bctype[3]==1) ucont[k][j-1][i].y = 0;
		if(k==0 && user->bctype[4]==1) ucont[k][j][i].z = 0;
		if(k==mz-1 && user->bctype[5]==1) ucont[k-1][j][i].z = 0;
		
		//cavity problem 
		if (j==my-1 && user->bctype[3]==2) ucont[k][j-1][i].y = 0;
		
		// couette flow j=0
		if (j==0 && user->bctype[2]==12) ucont[k][j][i].y = 0;
		
		// couette flow j=my-1
		if (j==my-1 && user->bctype[3]==12) ucont[k][j-1][i].y = 0;
		
		//slip BC
		if (user->bctype[0]==10 && i==0 && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) ucont[k][j][i].x = 0;
		if (user->bctype[1]==10 && i==mx-1 && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) ucont[k][j][i-1].x = 0;
		if (user->bctype[2]==10 && j==0 && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1) ) ucont[k][j][i].y = 0;
		if (std::abs(user->bctype[3])==10 && j==my-1 && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1) ) ucont[k][j-1][i].y = 0;
	}	
	DAVecRestoreArray(user->fda, user->Ucont, &ucont);
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	
	DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

	Contra2Cart_2(user);
	IB_BC(user);
	
	VecSet(Rhs,0);
	
	const double dt = user->dt;
	double coeff = time_coeff();
	if( coeff>0.9 && coeff<1.1 ) {
		VecAXPY(Rhs, -1./dt, user->Ucont);
		VecAXPY(Rhs, 1./dt, user->Ucont_o);
	}
	else/* if( coeff > 1.4 && coeff < 1.6 )*/ {
		VecAXPY(Rhs, -1.5/dt, user->Ucont);
		VecAXPY(Rhs, 2./dt, user->Ucont_o);
		VecAXPY(Rhs, -0.5/dt, user->Ucont_rm1);
	}/*
	else {
		VecAXPY(Rhs, -11./6./dt, user->Ucont);
		VecAXPY(Rhs, 3./dt, user->Ucont_o);
		VecAXPY(Rhs, -1.5/dt, user->Ucont_rm1);
		VecAXPY(Rhs, 1./3./dt, user->Ucont_rm2);
	}*/
	
	//Formfunction_2(user, Rhs, 1.0);	// careful ! adding values to Rhs
	if( coeff>0.9 && coeff<1.1 ) {
	  Formfunction_2(user, Rhs, 0.5);	// careful ! adding values to Rhs
	  VecAXPY(Rhs, 0.5, user->RHS_o);
	}
	else Formfunction_2(user, Rhs, 1.0);

	VecAXPY(Rhs, -1, user->dP);

	// begin add (xiaolei)
	if (nacelle_model || rotor_model || IB_delta) {
		VecAXPY(Rhs, 1, user->F_eul); 
	}
	// end add (xiaolei)	
	
	//if( time_coeff()>1.1 && time_coeff()<2.0 ) VecScale(Rhs, 1./1.5);
	return 0;
}

