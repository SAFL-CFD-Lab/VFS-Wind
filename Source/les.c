/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern double integrate_hat(int np, double *val, double *w);
extern double integrate_testfilter_i(double val[3][3][3], double vol[3][3][3]);
extern double integrate_testfilter_j(double val[3][3][3], double vol[3][3][3]);
extern double integrate_testfilter_k(double val[3][3][3], double vol[3][3][3]);
extern double integrate_testfilter_simpson(double val[3][3][3], double w[3][3][3]);

double integrate_gridfilter(double val[3][3][3], double vol[3][3][3]);
extern void Calculate_Covariant_tensor(double g[3][3], double G[3][3]);
extern void Calculate_Covariant_metrics(double g[3][3], double G[3][3]);

extern int mixed;
extern int les, ti, tistart;
extern PetscTruth rstart_flg;

const double les_eps=1.e-4;// osl 1.e-7
const double wall_cs=0.001;

void get_weight ( int i, int j, int k, int mx, int my, int mz, PetscReal ***aj, PetscReal ***nvert, PetscReal nv, double weight[3][3][3])
{
		for(int p=-1; p<=1; p++)
		for(int q=-1; q<=1; q++)
		for(int r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;

			if( nvert[K][J][I]>nv ) weight[R][Q][P]=0;
			else weight[R][Q][P] = 1./aj[K][J][I];
			/*
			if( i_periodic ) {
				if( I==0 ) I=mx-2;
				else if( I==mx-1 ) I=1;
			}
			else if( ii_periodic ) {
				if( I==0 ) I=-2;
				else if( I==mx-1 ) I=mx+1;
			}
			else if( I==0 || I==mx-1) weight[R][Q][P]=0;
			
			if( j_periodic ) {
				if( J==0 ) J=my-2;
				else if( J==my-1 ) J=1;
			}
			else if( jj_periodic ) {
				if( J==0 ) J=-2;
				else if( J==my-1 ) J=my+1;
			}
			else if(J==0 || J==my-1) weight[R][Q][P]=0;
			
			if( k_periodic ) {
				if( K==0 ) K=mz-2;
				else if( K==mz-1 ) K=1;
			}
			else if( kk_periodic ) {
				if( K==0 ) K=-2;
				else if( K==mz-1 ) K=mz+1;
			}
			else if( K==0 || K==mz-1 ) weight[R][Q][P]=0;
			*/
		}
}

void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat)
{
	if(ti<2 && tistart==0 && !rstart_flg){
		VecSet(user->lCs, 0.0);
		return;
	}
	
	if(les==1) {
		VecSet(user->lCs, 0.01);
		return;
	}
	
	Vec Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;

	Cmpnts	***ucont, ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
	
	PetscReal	***nvert, ***Cs;//, ***lnu_t;

	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
  
	PetscReal	***aj;
	Cmpnts ***Ax, ***Ay, ***Az, ***ucat_f, ***cent;//, ;
	double ***LM, ***MM;//, ***NM;
	PetscReal ***Sabs;//, ***Cs, ***lCs_o;//, ***nu_t;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	PetscReal	ajc;

	PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
	
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

	Vec lSx;		// Sx :	 Sxx, Sxy, Sxz
	Vec lSy;		// Sy :	 Syx, Syy, Syz
	Vec lSz;		// Sz :	 Szx, Szy, Szz
	Vec lS;
	Vec lUcat_f;
	
	DAVecGetArray(fda, Ucont, &ucont);
	DAVecGetArray(fda, Ucat,  &ucat);
	DAVecGetArray(fda, Csi, &csi);
	DAVecGetArray(fda, Eta, &eta);
	DAVecGetArray(fda, Zet, &zet);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lAj, &aj); 
	
	DAVecGetArray(da, user->lCs, &Cs);
	//DAVecGetArray(da, user->lNu_t, &lnu_t);
	
	DAVecGetArray(fda, user->lCent, &cent);

	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal	***iaj, ***jaj, ***kaj;

	DAVecGetArray(da, user->lIAj, &iaj);  
	DAVecGetArray(da, user->lJAj, &jaj);  
	DAVecGetArray(da, user->lKAj, &kaj);  
	
	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lIEta, &ieta);
	DAVecGetArray(fda, user->lIZet, &izet);

	DAVecGetArray(fda, user->lJCsi, &jcsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lJZet, &jzet);

	DAVecGetArray(fda, user->lKCsi, &kcsi);
	DAVecGetArray(fda, user->lKEta, &keta);
	DAVecGetArray(fda, user->lKZet, &kzet);
	//
  
	VecDuplicate(user->lP, &user->lLM);
	VecDuplicate(user->lP, &user->lMM);
	//VecDuplicate(user->lP, &user->lNM);
	
	VecSet(user->lLM, 0);
	VecSet(user->lMM, 0);
	//VecSet(user->lNM, 0);
	
	VecDuplicate(user->lUcont, &lSx);
	VecDuplicate(user->lUcont, &lSy);
	VecDuplicate(user->lUcont, &lSz);
	VecDuplicate(user->lNvert, &lS);
	VecDuplicate(user->lUcont, &lUcat_f);

	VecSet(lSx, 0);  
	VecSet(lSy, 0);  
	VecSet(lSz, 0);	
	VecSet(lS, 0);
	
	DAVecGetArray(da, user->lLM, &LM);//
	DAVecGetArray(da, user->lMM, &MM);//
	//DAVecGetArray(da, user->lNM, &NM);//
  	
	DAVecGetArray(fda, lSx, &Ax);//
	DAVecGetArray(fda, lSy, &Ay);//
	DAVecGetArray(fda, lSz, &Az);//
	DAVecGetArray(da, lS, &Sabs);//
	DAVecGetArray(fda, lUcat_f, &ucat_f);//
  	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if( nvert[k][j][i]>1.1) continue;
		
		ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

				
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	
		Sabs[k][j][i] = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );
		
		Ax[k][j][i].x = du_dx;	Ax[k][j][i].y = du_dy;	Ax[k][j][i].z = du_dz;
		Ay[k][j][i].x = dv_dx;	Ay[k][j][i].y = dv_dy;	Ay[k][j][i].z = dv_dz;
		Az[k][j][i].x = dw_dx;	Az[k][j][i].y = dw_dy;	Az[k][j][i].z = dw_dz;
		
		double weight[3][3][3];
		double u[3][3][3];
		double v[3][3][3];
		double w[3][3][3];
		get_weight (i, j, k, mx, my, mz, aj, nvert, 0.1, weight);
		
		for(int p=-1; p<=1; p++)
		for(int q=-1; q<=1; q++)
		for(int r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;
		}
		
		ucat_f[k][j][i].x = integrate_testfilter_simpson(u, weight);
		ucat_f[k][j][i].y = integrate_testfilter_simpson(v, weight);
		ucat_f[k][j][i].z = integrate_testfilter_simpson(w, weight);
	}
	
	DAVecRestoreArray(fda, lSx, &Ax);//
	DAVecRestoreArray(fda, lSy, &Ay);//
	DAVecRestoreArray(fda, lSz, &Az);//
	DAVecRestoreArray(da, lS, &Sabs);//
	DAVecRestoreArray(fda, lUcat_f, &ucat_f);//
		
	DALocalToLocalBegin(fda, lSx, INSERT_VALUES, lSx);
	DALocalToLocalEnd(fda, lSx, INSERT_VALUES, lSx);
 
	DALocalToLocalBegin(fda, lSy, INSERT_VALUES, lSy);
	DALocalToLocalEnd(fda, lSy, INSERT_VALUES, lSy);
 
	DALocalToLocalBegin(fda, lSz, INSERT_VALUES, lSz);
	DALocalToLocalEnd(fda, lSz, INSERT_VALUES, lSz);
 
	DALocalToLocalBegin(da, lS, INSERT_VALUES, lS);
	DALocalToLocalEnd(da, lS, INSERT_VALUES, lS);
	
	DALocalToLocalBegin(fda, lUcat_f, INSERT_VALUES, lUcat_f);
	DALocalToLocalEnd(fda, lUcat_f, INSERT_VALUES, lUcat_f);
	
	DAVecGetArray(fda, lSx, &Ax);//
	DAVecGetArray(fda, lSy, &Ay);//
	DAVecGetArray(fda, lSz, &Az);//
	DAVecGetArray(da, lS, &Sabs);//
	DAVecGetArray(fda, lUcat_f, &ucat_f);//
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int c=k, b=j, a=i, flag=0;
		
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
			Sabs[k][j][i] = Sabs[c][b][a];
			Ax[k][j][i] = Ax[c][b][a];
			Ay[k][j][i] = Ay[c][b][a];
			Az[k][j][i] = Az[c][b][a];
			ucat_f[k][j][i] = ucat_f[c][b][a];
		}
	}
	
 	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			LM[k][j][i]=MM[k][j][i]=0;//NM[k][j][i]=0;
			continue;
		}
		ajc = aj[k][j][i];
	
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	
		int a, b;
		double Lij[3][3], Sij_hat[3][3], SSij_hat[3][3], Mij[3][3], Nij[3][3], Nij_cat[3][3];
		double Lij_cat[3][3], Mij_cat[3][3];
		
		double filter, test_filter;
		double S[3][3][3], S_hat;
		double u[3][3][3], v[3][3][3], w[3][3][3];
		
		double U[3][3][3], V[3][3][3], W[3][3][3];
		double Uu[3][3][3], Uv[3][3][3], Uw[3][3][3];
		double Vu[3][3][3], Vv[3][3][3], Vw[3][3][3];
		double Wu[3][3][3], Wv[3][3][3], Ww[3][3][3];
		
		double uu[3][3][3], uv[3][3][3], uw[3][3][3];
		double vv[3][3][3], vw[3][3][3], ww[3][3][3];
		double S11[3][3][3], S12[3][3][3], S13[3][3][3], S21[3][3][3], S22[3][3][3], S23[3][3][3], S31[3][3][3], S32[3][3][3], S33[3][3][3];
		double SS11[3][3][3], SS12[3][3][3], SS13[3][3][3], SS21[3][3][3], SS22[3][3][3], SS23[3][3][3], SS31[3][3][3], SS32[3][3][3], SS33[3][3][3];
		
		double dU_dx[3][3][3], dU_dy[3][3][3], dU_dz[3][3][3];
		double dV_dx[3][3][3], dV_dy[3][3][3], dV_dz[3][3][3];
		double dW_dx[3][3][3], dW_dy[3][3][3], dW_dz[3][3][3];
		double T11[3][3][3], T12[3][3][3], T13[3][3][3], T21[3][3][3], T22[3][3][3], T23[3][3][3], T31[3][3][3], T32[3][3][3], T33[3][3][3];	// Clark model term
		
		double weight[3][3][3];
		
		double dx,dy,dz;
		Calculate_dxdydz(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
		double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
		

		get_weight (i, j, k, mx, my, mz, aj, nvert, 0.1, weight);
		
		int p,q,r;
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;

			u[R][Q][P] = v[R][Q][P] = w[R][Q][P] = 0;
			uu[R][Q][P] = uv[R][Q][P] = uw[R][Q][P] = 0;
			vv[R][Q][P] = vw[R][Q][P] = ww[R][Q][P] = 0;
			
			S11[R][Q][P] = S12[R][Q][P] = S13[R][Q][P] = 0;
			S21[R][Q][P] = S22[R][Q][P] = S23[R][Q][P] = 0;
			S31[R][Q][P] = S32[R][Q][P] = S33[R][Q][P] = 0;
						
			S[R][Q][P] = 0;
				
			SS11[R][Q][P] = SS12[R][Q][P] = SS13[R][Q][P] = 0;
			SS21[R][Q][P] = SS22[R][Q][P] = SS23[R][Q][P] = 0;
			SS31[R][Q][P] = SS32[R][Q][P] = SS33[R][Q][P] = 0;
						
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;
			
			// metric tensors are also test-filtered : Big difference
			
			/*if(I==0)*/ U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
			//else U[R][Q][P] = 0.5 * ( ucont[K][J][I].x + ucont[K][J][I-1].x );
			
			/*if(J==0)*/ V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
			//else V[R][Q][P] = 0.5 * ( ucont[K][J][I].y + ucont[K][J-1][I].y );
			
			/*if(K==0)*/ W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;
			//else W[R][Q][P] = 0.5 * ( ucont[K][J][I].z + ucont[K-1][J][I].z );
			
			Uu[R][Q][P] = U[R][Q][P]*u[R][Q][P];
			Uv[R][Q][P] = U[R][Q][P]*v[R][Q][P];
			Uw[R][Q][P] = U[R][Q][P]*w[R][Q][P];
			
			Vu[R][Q][P] = V[R][Q][P]*u[R][Q][P];
			Vv[R][Q][P] = V[R][Q][P]*v[R][Q][P];
			Vw[R][Q][P] = V[R][Q][P]*w[R][Q][P];
			
			Wu[R][Q][P] = W[R][Q][P]*u[R][Q][P];
			Wv[R][Q][P] = W[R][Q][P]*v[R][Q][P];
			Ww[R][Q][P] = W[R][Q][P]*w[R][Q][P];
			
			uu[R][Q][P] = u[R][Q][P]*u[R][Q][P];	
			uv[R][Q][P] = u[R][Q][P]*v[R][Q][P];	
			uw[R][Q][P] = u[R][Q][P]*w[R][Q][P];
			vv[R][Q][P] = v[R][Q][P]*v[R][Q][P];	
			vw[R][Q][P] = v[R][Q][P]*w[R][Q][P];	
			ww[R][Q][P] = w[R][Q][P]*w[R][Q][P];
			
			const double du_dx = Ax[K][J][I].x, du_dy = Ax[K][J][I].y, du_dz = Ax[K][J][I].z;
			const double dv_dx = Ay[K][J][I].x, dv_dy = Ay[K][J][I].y, dv_dz = Ay[K][J][I].z;
			const double dw_dx = Az[K][J][I].x, dw_dy = Az[K][J][I].y, dw_dz = Az[K][J][I].z;
		
			const double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
			const double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
			const double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
			
			dU_dx[R][Q][P] = du_dx, dU_dy[R][Q][P] = du_dy, dU_dz[R][Q][P] = du_dz;
			dV_dx[R][Q][P] = dv_dx, dV_dy[R][Q][P] = dv_dy, dV_dz[R][Q][P] = dv_dz;
			dW_dx[R][Q][P] = dw_dx, dW_dy[R][Q][P] = dw_dy, dW_dz[R][Q][P] = dw_dz;
			
			T11[R][Q][P] = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 ) / 12.;
			T12[R][Q][P] = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 ) / 12.;
			T13[R][Q][P] = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 ) / 12.;
			T21[R][Q][P] = T12[R][Q][P];
			T22[R][Q][P] = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 ) / 12.;
			T23[R][Q][P] = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 ) / 12.;
			T31[R][Q][P] = T13[R][Q][P];
			T32[R][Q][P] = T23[R][Q][P];
			T33[R][Q][P] = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 ) / 12.;
		
			S11[R][Q][P] = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
			S21[R][Q][P] = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
			S31[R][Q][P] = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;
			
			S[R][Q][P] = Sabs[K][J][I];
			
			SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
			SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
			SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];
		}
		
		double sum_weight=0;
		double coef[3][3][3]={
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250, 
			0.125, 0.250, 0.125, 
				
			0.250, 0.500, 0.250,
			0.500, 1.000, 0.500,
			0.250, 0.500, 0.250,
				
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250,
			0.125, 0.250, 0.125
		};
		
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			sum_weight += weight[r+1][q+1][p+1] * coef[r+1][q+1][p+1];
		}
		
		filter = pow( 1./aj[k][j][i], 1./3. );
		
		if(testfilter_ik) test_filter = pow(5.0, 1./3.) * filter;
		else {
			//test_filter = 2.0 * filter;
			test_filter = pow( sum_weight, 1./3. );
		}
		/*		
		double _U=integrate_testfilter_simpson(U, weight);
		double _V=integrate_testfilter_simpson(V, weight);
		double _W=integrate_testfilter_simpson(W, weight);
		
		double _u=integrate_testfilter_simpson(u, weight);
		double _v=integrate_testfilter_simpson(v, weight);
		double _w=integrate_testfilter_simpson(w, weight);
		*/
		double _u=ucat_f[k][j][i].x;
		double _v=ucat_f[k][j][i].y;
		double _w=ucat_f[k][j][i].z;
		double _U = _u * csi[k][j][i].x + _v * csi[k][j][i].y + _w * csi[k][j][i].z;
		double _V = _u * eta[k][j][i].x + _v * eta[k][j][i].y + _w * eta[k][j][i].z;
		double _W = _u * zet[k][j][i].x + _v * zet[k][j][i].y + _w * zet[k][j][i].z;
		
		double _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz;
		double _du_dx, _du_dy, _du_dz, _dv_dx, _dv_dy, _dv_dz, _dw_dx, _dw_dy, _dw_dz;
		Compute_du_center ( i, j, k, mx, my, mz, ucat_f, nvert, &_dudc, &_dvdc, &_dwdc, &_dude, &_dvde, &_dwde, &_dudz, &_dvdz, &_dwdz);        // derivative of ucat_f
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz, &_du_dx, &_dv_dx, &_dw_dx, &_du_dy, &_dv_dy, &_dw_dy, &_du_dz, &_dv_dz, &_dw_dz );
		double _Sxx = 0.5*( _du_dx + _du_dx ), _Sxy = 0.5*(_du_dy + _dv_dx), _Sxz = 0.5*(_du_dz + _dw_dx);
		double _Syx = _Sxy, _Syy = 0.5*(_dv_dy + _dv_dy), _Syz = 0.5*(_dv_dz + _dw_dy);
		double _Szx = _Sxz, _Szy=_Syz, _Szz = 0.5*(_dw_dz + _dw_dz);
		Sij_hat[0][0] = _Sxx, Sij_hat[0][1] = _Sxy, Sij_hat[0][2] = _Sxz;
		Sij_hat[1][0] = _Syx, Sij_hat[1][1] = _Syy, Sij_hat[1][2] = _Syz;
		Sij_hat[2][0] = _Szx, Sij_hat[2][1] = _Szy, Sij_hat[2][2] = _Szz;
		S_hat = sqrt( 2.0 *( _Sxx*_Sxx + _Sxy*_Sxy + _Sxz*_Sxz + _Syx*_Syx + _Syy*_Syy + _Syz*_Syz + _Szx*_Szx + _Szy*_Szy + _Szz*_Szz ) );

		double _T11=integrate_testfilter_simpson(T11, weight);
		double _T12=integrate_testfilter_simpson(T12, weight);
		double _T13=integrate_testfilter_simpson(T13, weight);
		double _T21=integrate_testfilter_simpson(T21, weight);
		double _T22=integrate_testfilter_simpson(T22, weight);
		double _T23=integrate_testfilter_simpson(T23, weight);
		double _T31=integrate_testfilter_simpson(T31, weight);
		double _T32=integrate_testfilter_simpson(T32, weight);
		double _T33=integrate_testfilter_simpson(T33, weight);
		/*
		double _dU_dx=integrate_testfilter_simpson(dU_dx, weight);
		double _dV_dx=integrate_testfilter_simpson(dV_dx, weight);
		double _dW_dx=integrate_testfilter_simpson(dW_dx, weight);
		double _dU_dy=integrate_testfilter_simpson(dU_dy, weight);
		double _dV_dy=integrate_testfilter_simpson(dV_dy, weight);
		double _dW_dy=integrate_testfilter_simpson(dW_dy, weight);
		double _dU_dz=integrate_testfilter_simpson(dU_dz, weight);
		double _dV_dz=integrate_testfilter_simpson(dV_dz, weight);
		double _dW_dz=integrate_testfilter_simpson(dW_dz, weight);
		*/
		double _R11 = ( _du_dx * _du_dx * dx2 + _du_dy * _du_dy * dy2 + _du_dz * _du_dz * dz2 ) / 12.;
		double _R12 = ( _du_dx * _dv_dx * dx2 + _du_dy * _dv_dy * dy2 + _du_dz * _dv_dz * dz2 ) / 12.;
		double _R13 = ( _du_dx * _dw_dx * dx2 + _du_dy * _dw_dy * dy2 + _du_dz * _dw_dz * dz2 ) / 12.;
		double _R21 = _R12;
		double _R22 = ( _dv_dx * _dv_dx * dx2 + _dv_dy * _dv_dy * dy2 + _dv_dz * _dv_dz * dz2 ) / 12.;
		double _R23 = ( _dv_dx * _dw_dx * dx2 + _dv_dy * _dw_dy * dy2 + _dv_dz * _dw_dz * dz2 ) / 12.;
		double _R31 = _R13;
		double _R32 = _R23;
		double _R33 = ( _dw_dx * _dw_dx * dx2 + _dw_dy * _dw_dy * dy2 + _dw_dz * _dw_dz * dz2 ) / 12.;
		
		Lij[0][0] = integrate_testfilter_simpson(Uu, weight) - _U*_u;
		Lij[0][1] = integrate_testfilter_simpson(Uv, weight) - _U*_v;
		Lij[0][2] = integrate_testfilter_simpson(Uw, weight) - _U*_w;
		Lij[1][0] = integrate_testfilter_simpson(Vu, weight) - _V*_u;
		Lij[1][1] = integrate_testfilter_simpson(Vv, weight) - _V*_v;
		Lij[1][2] = integrate_testfilter_simpson(Vw, weight) - _V*_w;
		Lij[2][0] = integrate_testfilter_simpson(Wu, weight) - _W*_u;
		Lij[2][1] = integrate_testfilter_simpson(Wv, weight) - _W*_v;
		Lij[2][2] = integrate_testfilter_simpson(Ww, weight) - _W*_w;
		
		
		Nij_cat[0][0] = 4.0 * _R11  - _T11;
		Nij_cat[0][1] = 4.0 * _R12  - _T12;
		Nij_cat[0][2] = 4.0 * _R13  - _T13;
		Nij_cat[1][0] = 4.0 * _R21  - _T21;
		Nij_cat[1][1] = 4.0 * _R22  - _T22;
		Nij_cat[1][2] = 4.0 * _R23  - _T23;
		Nij_cat[2][0] = 4.0 * _R31  - _T31;
		Nij_cat[2][1] = 4.0 * _R32  - _T32;
		Nij_cat[2][2] = 4.0 * _R33  - _T33;
		
		Nij[0][0] = Nij_cat[0][0] * csi0 + Nij_cat[0][1] * csi1 + Nij_cat[0][2] * csi2;
		Nij[0][1] = Nij_cat[0][0] * eta0 + Nij_cat[0][1] * eta1 + Nij_cat[0][2] * eta2;
		Nij[0][2] = Nij_cat[0][0] * zet0 + Nij_cat[0][1] * zet1 + Nij_cat[0][2] * zet2;
		Nij[1][0] = Nij_cat[1][0] * csi0 + Nij_cat[1][1] * csi1 + Nij_cat[1][2] * csi2;
		Nij[1][1] = Nij_cat[1][0] * eta0 + Nij_cat[1][1] * eta1 + Nij_cat[1][2] * eta2;
		Nij[1][2] = Nij_cat[1][0] * zet0 + Nij_cat[1][1] * zet1 + Nij_cat[1][2] * zet2;
		Nij[2][0] = Nij_cat[2][0] * csi0 + Nij_cat[2][1] * csi1 + Nij_cat[2][2] * csi2;
		Nij[2][1] = Nij_cat[2][0] * eta0 + Nij_cat[2][1] * eta1 + Nij_cat[2][2] * eta2;
		Nij[2][2] = Nij_cat[2][0] * zet0 + Nij_cat[2][1] * zet1 + Nij_cat[2][2] * zet2;
		
		/*
		double _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz;
		double _du_dx, _du_dy, _du_dz, _dv_dx, _dv_dy, _dv_dz, _dw_dx, _dw_dy, _dw_dz;
		Compute_du_center ( i, j, k, mx, my, mz, ucat_f, nvert, &_dudc, &_dvdc, &_dwdc, &_dude, &_dvde, &_dwde, &_dudz, &_dvdz, &_dwdz);	// derivative of ucat_f
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz, &_du_dx, &_dv_dx, &_dw_dx, &_du_dy, &_dv_dy, &_dw_dy, &_du_dz, &_dv_dz, &_dw_dz );
		double _Sxx = 0.5*( _du_dx + _du_dx ), _Sxy = 0.5*(_du_dy + _dv_dx), _Sxz = 0.5*(_du_dz + _dw_dx);
		double _Syx = _Sxy, _Syy = 0.5*(_dv_dy + _dv_dy), _Syz = 0.5*(_dv_dz + _dw_dy);
		double _Szx = _Sxz, _Szy=_Syz, _Szz = 0.5*(_dw_dz + _dw_dz);
		Sij_hat[0][0] = _Sxx, Sij_hat[0][1] = _Sxy, Sij_hat[0][2] = _Sxz;
		Sij_hat[1][0] = _Syx, Sij_hat[1][1] = _Syy, Sij_hat[1][2] = _Syz;
		Sij_hat[2][0] = _Szx, Sij_hat[2][1] = _Szy, Sij_hat[2][2] = _Szz;
		S_hat = sqrt( 2.0 *( _Sxx*_Sxx + _Sxy*_Sxy + _Sxz*_Sxz + _Syx*_Syx + _Syy*_Syy + _Syz*_Syz + _Szx*_Szx + _Szy*_Szy + _Szz*_Szz ) );
		*/
		//
		
		/*
		Sij_hat[0][0] = integrate_testfilter_simpson(S11, weight);	
		Sij_hat[0][1] = integrate_testfilter_simpson(S12, weight);	
		Sij_hat[0][2] = integrate_testfilter_simpson(S13, weight);
		Sij_hat[1][0] = integrate_testfilter_simpson(S21, weight);	
		Sij_hat[1][1] = integrate_testfilter_simpson(S22, weight);	
		Sij_hat[1][2] = integrate_testfilter_simpson(S23, weight);
		Sij_hat[2][0] = integrate_testfilter_simpson(S31, weight);	
		Sij_hat[2][1] = integrate_testfilter_simpson(S32, weight);	
		Sij_hat[2][2] = integrate_testfilter_simpson(S33, weight);
		
		S_hat=0;
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			S_hat += pow( Sij_hat[a][b], 2. );
		}
		S_hat = sqrt ( 2 * S_hat );
		*/
		
		SSij_hat[0][0] = integrate_testfilter_simpson(SS11, weight);	
		SSij_hat[0][1] = integrate_testfilter_simpson(SS12, weight);	
		SSij_hat[0][2] = integrate_testfilter_simpson(SS13, weight);
		SSij_hat[1][0] = integrate_testfilter_simpson(SS21, weight);
		SSij_hat[1][1] = integrate_testfilter_simpson(SS22, weight);	
		SSij_hat[1][2] = integrate_testfilter_simpson(SS23, weight);
		SSij_hat[2][0] = integrate_testfilter_simpson(SS31, weight);	
		SSij_hat[2][1] = integrate_testfilter_simpson(SS32, weight);	
		SSij_hat[2][2] = integrate_testfilter_simpson(SS33, weight);
		
		
		
		//S_hat = integrate_testfilter_simpson(S, weight);
		
		
		double gg[3][3], ggc[3][3], G[3][3];
		double xcsi, xeta, xzet, ycsi, yeta, yzet, zcsi, zeta, zzet;

		gg[0][0]=csi0, gg[0][1]=csi1, gg[0][2]=csi2;
		gg[1][0]=eta0, gg[1][1]=eta1, gg[1][2]=eta2;
		gg[2][0]=zet0, gg[2][1]=zet1, gg[2][2]=zet2;
		Calculate_Covariant_metrics(gg, ggc);
		xcsi=ggc[0][0], xeta=ggc[0][1], xzet=ggc[0][2];
		ycsi=ggc[1][0], yeta=ggc[1][1], yzet=ggc[1][2];
		zcsi=ggc[2][0], zeta=ggc[2][1], zzet=ggc[2][2];
		G[0][0] = xcsi * xcsi + ycsi * ycsi + zcsi * zcsi;
		G[1][1] = xeta * xeta + yeta * yeta + zeta * zeta;
		G[2][2] = xzet * xzet + yzet * yzet + zzet * zzet;
		G[0][1] = G[1][0] = xeta * xcsi + yeta * ycsi + zeta * zcsi;
		G[0][2] = G[2][0] = xzet * xcsi + yzet * ycsi + zzet * zcsi;
		G[1][2] = G[2][1] = xeta * xzet + yeta * yzet + zeta * zzet;
		
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
		  	Mij_cat[a][b] = -pow( test_filter, 2. ) * S_hat * Sij_hat[a][b] + pow( filter, 2. ) * SSij_hat[a][b];
			/*
			Mij_cat_i[a][b] = -4. * S_hat_i * Sij_hat_i[a][b] + SSij_hat_i[a][b];
			Mij_cat_j[a][b] = -4. * S_hat_j * Sij_hat_j[a][b] + SSij_hat_j[a][b];
			Mij_cat_k[a][b] = -4. * S_hat_k * Sij_hat_k[a][b] + SSij_hat_k[a][b];*/
		}
		
		Mij[0][0] = Mij_cat[0][0] * csi0 + Mij_cat[0][1] * csi1 + Mij_cat[0][2] * csi2;
		Mij[0][1] = Mij_cat[0][0] * eta0 + Mij_cat[0][1] * eta1 + Mij_cat[0][2] * eta2;
		Mij[0][2] = Mij_cat[0][0] * zet0 + Mij_cat[0][1] * zet1 + Mij_cat[0][2] * zet2;
		Mij[1][0] = Mij_cat[1][0] * csi0 + Mij_cat[1][1] * csi1 + Mij_cat[1][2] * csi2;
		Mij[1][1] = Mij_cat[1][0] * eta0 + Mij_cat[1][1] * eta1 + Mij_cat[1][2] * eta2;
		Mij[1][2] = Mij_cat[1][0] * zet0 + Mij_cat[1][1] * zet1 + Mij_cat[1][2] * zet2;
		Mij[2][0] = Mij_cat[2][0] * csi0 + Mij_cat[2][1] * csi1 + Mij_cat[2][2] * csi2;
		Mij[2][1] = Mij_cat[2][0] * eta0 + Mij_cat[2][1] * eta1 + Mij_cat[2][2] * eta2;
		Mij[2][2] = Mij_cat[2][0] * zet0 + Mij_cat[2][1] * zet1 + Mij_cat[2][2] * zet2;
					
		double num=0, num1=0, denom=0;
		int m, n, l;
		
		/*
			g11 ~ csi*csi ~ dx^4
			G11 ~ dx^-4
		*/
	
	
		for(q=0; q<3; q++)
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			num += Lij[b][a] * Mij[a][q] * G[b][q];
			if(clark) num -= Nij[b][a] * Mij[a][q] * G[b][q];
		}
		
		for(m=0; m<3; m++)
		for(n=0; n<3; n++)
		for(l=0; l<3; l++) {
			denom += Mij[n][m] * Mij[n][l] * G[m][l];
		}
	
		//printf("%f %f\n", num, denom);
		
		LM[k][j][i] = num;
		MM[k][j][i] = denom;
    	}
	
	DAVecRestoreArray(da, user->lLM, &LM);//
	DAVecRestoreArray(da, user->lMM, &MM);//
	//DAVecRestoreArray(da, user->lNM, &NM);//
	
	DALocalToLocalBegin(da, user->lLM, INSERT_VALUES, user->lLM);
	DALocalToLocalEnd(da, user->lLM, INSERT_VALUES, user->lLM);
	DALocalToLocalBegin(da, user->lMM, INSERT_VALUES, user->lMM);
	DALocalToLocalEnd(da, user->lMM, INSERT_VALUES, user->lMM);
	//DALocalToLocalBegin(da, user->lNM, INSERT_VALUES, user->lNM);
	//DALocalToLocalEnd(da, user->lNM, INSERT_VALUES, user->lNM);
	  
	DAVecGetArray(da, user->lLM, &LM);//
	DAVecGetArray(da, user->lMM, &MM);//
	//DAVecGetArray(da, user->lNM, &NM);//
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int c=k, b=j, a=i, flag=0;
		
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
			LM[k][j][i] = LM[c][b][a];
			MM[k][j][i] = MM[c][b][a];
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			Cs[k][j][i] = 0;
			continue;
		}
				
		double weight[3][3][3];
		double LM0[3][3][3], MM0[3][3][3];
		int a, b, c;
			
		for(a=-1; a<=1; a++)
		for(b=-1; b<=1; b++)
		for(c=-1; c<=1; c++) {
			int R=c+1, Q=b+1, P=a+1;
			int K=k+c, J=j+b, I=i+a;
			
			weight[R][Q][P] = 1./aj[K][J][I];
			
			if( nvert[K][J][I]>1.1 ) weight[R][Q][P]=0;
			
			if( i_periodic ) {
				if( I==0 ) I=mx-2;
				else if( I==mx-1 ) I=1;
			}
			else if( ii_periodic ) {
				if( I==0 ) I=-2;
				else if( I==mx-1 ) I=mx+1;
			}
			else if( I==0 || I==mx-1) weight[R][Q][P]=0;
			
			if( j_periodic ) {
				if( J==0 ) J=my-2;
				else if( J==my-1 ) J=1;
			}
			else if( jj_periodic ) {
				if( J==0 ) J=-2;
				else if( J==my-1 ) J=my+1;
			}
			else if( J==0 || j==my-1) weight[R][Q][P]=0;
			
			if( k_periodic ) {
				if( K==0 ) K=mz-2;
				else if( K==mz-1 ) K=1;
			}
			else if( kk_periodic ) {
				if( K==0 ) K=-2;
				else if( K==mz-1 ) K=mz+1;
			}
			else if( K==0 || K==mz-1) weight[R][Q][P]=0;
			
			LM0[R][Q][P] = LM[K][J][I];
			MM0[R][Q][P] = MM[K][J][I];
		}
			
		double C=0;
			
		double LM_avg, MM_avg;//, NM_avg;

		if ( i_homo_filter || j_homo_filter || k_homo_filter || les==3 ) {
			LM_avg = LM[k][j][i];
			MM_avg = MM[k][j][i];
		}
		else {
			LM_avg = integrate_testfilter_simpson(LM0,weight);
			MM_avg = integrate_testfilter_simpson(MM0,weight);
		}
		
		C = 0.5 * LM_avg / (MM_avg + les_eps );
		
		if ( les==3 ) {
			if(ti<100 && tistart==0 && !rstart_flg) { }
			else {
				C = (1.0 - 0.001) * Cs[k][j][i] + 0.001 * C;
			}
		}
		
		if(les==1) Cs[k][j][i] = 0.01;
		else Cs[k][j][i] = PetscMax(C, 0);
	}
	
	if( les==3 ) {}
	else if ( i_homo_filter && k_homo_filter ) {
		std::vector<int> count, total_count;
		std::vector<double> J_LM(my), J_MM(my), LM_tmp(my), MM_tmp(my);
		
		count.resize(my);
		total_count.resize(my);
		
		for(j=0; j<my; j++) {
			LM_tmp[j] = 0;
			MM_tmp[j] = 0;
			count[j] = total_count[j] = 0;
		}
		
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if( nvert[k][j][i]<0.1 ) {
				LM_tmp[j] += LM[k][j][i];
				MM_tmp[j] += MM[k][j][i];
				count[j] ++;
			}
		}
		
		MPI_Allreduce( &LM_tmp[0], &J_LM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &MM_tmp[0], &J_MM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &count[0], &total_count[0], my, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	
		
		for(j=0; j<my; j++) {
			if( total_count[j]>0) {
				J_LM[j] /= (double) (total_count[j]);
				J_MM[j] /= (double) (total_count[j]);
			}
		}
		
		for (j=lys; j<lye; j++)
		for (k=lzs; k<lze; k++)
		for (i=lxs; i<lxe; i++) {
			Cs[k][j][i] = 0.5 * J_LM[j] / ( J_MM[j]+les_eps);
		}
	}
	else if (i_homo_filter || j_homo_filter || k_homo_filter) {
		int plane_size;
		
		if(i_homo_filter) plane_size = my*mz;
		else if(j_homo_filter) plane_size = mx*mz;
		else if(k_homo_filter) plane_size = mx*my;
		
		std::vector<int> count(plane_size), total_count(plane_size);
		std::vector<double> J_LM(plane_size), J_MM(plane_size), LM_tmp(plane_size), MM_tmp(plane_size);
		int pos;
		
		pos=0;
		
		std::fill( LM_tmp.begin(), LM_tmp.end(), 0.);
		std::fill( MM_tmp.begin(), MM_tmp.end(), 0.);
		std::fill( count.begin(), count.end(), 0);
		std::fill( total_count.begin(), total_count.end(), 0);
		
		for(pos=0; pos<plane_size; pos++) {
			LM_tmp[pos] = 0;
			MM_tmp[pos] = 0;
			count[pos] = total_count[pos] = 0;
		}
		
		pos=0;
		
		if(i_homo_filter)  {
			for(k=0; k<mz; k++)
			for(j=0; j<my; j++) {
				if( k>=lzs && k<lze && j>=lys && j<lye) {
					for (i=lxs; i<lxe; i++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		else if(j_homo_filter)  {
			for(k=0; k<mz; k++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && k>=lzs && k<lze) {
					for (j=lys; j<lye; j++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		else if(k_homo_filter)  {
			for(j=0; j<my; j++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && j>=lys && j<lye) {
					for (k=lzs; k<lze; k++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		
		MPI_Allreduce( &LM_tmp[0], &J_LM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &MM_tmp[0], &J_MM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &count[0], &total_count[0], plane_size, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	
		
		pos=0;
		
		for(pos=0; pos<plane_size; pos++) {
			if( total_count[pos]>0) {
				double N = (double) (total_count[pos]);
				J_LM[pos] /= N;
				J_MM[pos] /= N;
			}
		}
		
		pos=0;
		if(i_homo_filter)  {
			for(k=0; k<mz; k++)
			for(j=0; j<my; j++) {
				if( k>=lzs && k<lze && j>=lys && j<lye) {
					for (i=lxs; i<lxe; i++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
		else if(j_homo_filter)  {
			for(k=0; k<mz; k++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && k>=lzs && k<lze) {
					for (j=lys; j<lye; j++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
		else if(k_homo_filter)  {
			for(j=0; j<my; j++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && j>=lys && j<lye) {
					for (k=lzs; k<lze; k++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
	}
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		
		if(nvert[k][j][i]>1.1 || k==0 || k==mz-1 || j==0 || j==my-1 || i==0 || i==mx-1) {
			Cs[k][j][i] = 0;
		}
		else {
			if(nvert[k][j][i]>0.1 && nvert[k][j][i]<1.1) {
				Cs[k][j][i] = PetscMax(wall_cs, Cs[k][j][i]);	// stabilize at high Re, osl 0.005
			}
			Cs[k][j][i] = PetscMin( PetscMax ( Cs[k][j][i], 0 ), max_cs);
		}
	}
	
	DAVecRestoreArray(fda, lSx, &Ax);//
	DAVecRestoreArray(fda, lSy, &Ay);//
	DAVecRestoreArray(fda, lSz, &Az);//
	DAVecRestoreArray(da, lS, &Sabs);//
	DAVecRestoreArray(fda, lUcat_f, &ucat_f);//
	
	DAVecRestoreArray(da, user->lLM, &LM);//
	DAVecRestoreArray(da, user->lMM, &MM);//

	DAVecRestoreArray(fda, Ucont, &ucont);
	DAVecRestoreArray(fda, Ucat,  &ucat);
	DAVecRestoreArray(fda, Csi, &csi);
	DAVecRestoreArray(fda, Eta, &eta);
	DAVecRestoreArray(fda, Zet, &zet);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lAj, &aj); 
	DAVecRestoreArray(da, user->lCs, &Cs);
	
	DAVecRestoreArray(fda, user->lCent, &cent);
	VecDestroy(user->lLM);//
	VecDestroy(user->lMM);//
  
	VecDestroy(lSx);//
	VecDestroy(lSy);//
	VecDestroy(lSz);//
	VecDestroy(lS);//
	VecDestroy(lUcat_f);//

	DAVecRestoreArray(da, user->lIAj, &iaj);  
	DAVecRestoreArray(da, user->lJAj, &jaj);  
	DAVecRestoreArray(da, user->lKAj, &kaj);  
	
	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);

	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);

	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);
	
	DALocalToLocalBegin(da, user->lCs, INSERT_VALUES, user->lCs);
	DALocalToLocalEnd(da, user->lCs, INSERT_VALUES, user->lCs);
	
	DAVecGetArray(da, user->lCs, &Cs);
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


		if(flag) Cs[k][j][i] = Cs[c][b][a];
	}
	DAVecRestoreArray(da, user->lCs, &Cs);
	/*
	if ( ! i_homo_filter && ! j_homo_filter && ! k_homo_filter && les!=3 ) {
		Vec lCs_tmp;
		PetscReal ***Cs_tmp;
		
		VecDuplicate(user->lCs, &lCs_tmp);
		VecCopy(user->lCs, lCs_tmp);
		
		DAVecGetArray(da, user->lAj, &aj);
		DAVecGetArray(da, user->lNvert, &nvert);
		DAVecGetArray(da, user->lCs, &Cs);
		DAVecGetArray(da, lCs_tmp, &Cs_tmp);
		for(k=lzs; k<lze; k++)
		for(j=lys; j<lye; j++)
		for(i=lxs; i<lxe; i++) {
			if(nvert[k][j][i]<0.1) {
				double weight[3][3][3], C[3][3][3];
				int a, b, c;

				for(a=-1; a<=1; a++)
				for(b=-1; b<=1; b++)
				for(c=-1; c<=1; c++) {
					int R=c+1, Q=b+1, P=a+1;
					int K=k+c, J=j+b, I=i+a;
					
					weight[R][Q][P] = 1./aj[K][J][I];
					
					if( nvert[K][J][I]>1.1 ) weight[R][Q][P]=0;
					
					if( i_periodic ) {
						if( I==0 ) I=mx-2;
						else if( I==mx-1 ) I=1;
					}
					else if( ii_periodic ) {
						if( I==0 ) I=-2;
						else if( I==mx-1 ) I=mx+1;
					}
					else if( I==0 || I==mx-1) weight[R][Q][P]=0;
					
					if( j_periodic ) {
						if( J==0 ) J=my-2;
						else if( J==my-1 ) J=1;
					}
					else if( jj_periodic ) {
						if( J==0 ) J=-2;
						else if( J==my-1 ) J=my+1;
					}
					else if( J==0 || j==my-1) weight[R][Q][P]=0;
					
					if( k_periodic ) {
						if( K==0 ) K=mz-2;
						else if( K==mz-1 ) K=1;
					}
					else if( kk_periodic ) {
						if( K==0 ) K=-2;
						else if( K==mz-1 ) K=mz+1;
					}
					else if( K==0 || K==mz-1) weight[R][Q][P]=0;
					
					C[R][Q][P] = Cs_tmp[K][J][I];
				}
				Cs[k][j][i] = integrate_testfilter_simpson(C,weight);
			}
		}
		DAVecRestoreArray(da, user->lAj, &aj);
		DAVecRestoreArray(da, user->lNvert, &nvert);
		DAVecRestoreArray(da, user->lCs, &Cs);
		DAVecRestoreArray(da, lCs_tmp, &Cs_tmp);
		
		VecDestroy(lCs_tmp);
	}
	*/
	double lmax_norm=0, max_norm;
	int p;
	
	if(testfilter_ik) PetscPrintf(PETSC_COMM_WORLD, "\nFilter type : Box filter homogeneous\n");
	else PetscPrintf(PETSC_COMM_WORLD, "Filter type : Box filter 3D\n");
	if(clark) PetscPrintf(PETSC_COMM_WORLD, "Clark model\n");
	
	VecMax(user->lCs, &p, &lmax_norm);
	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Cs  = %e \n", max_norm);
	
};

void Compute_eddy_viscosity_LES(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	
	PetscReal ***Cs, ***lnu_t, ***nvert, ***aj, ***ustar;
	Cmpnts ***csi, ***eta, ***zet, ***ucat;
	
	DAGetLocalInfo(da, &info);
	mx = info.mx, my = info.my, mz = info.mz;
	xs = info.xs, xe = xs + info.xm;
	ys = info.ys, ye = ys + info.ym;
	zs = info.zs, ze = zs + info.zm;

	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	VecSet(user->lNu_t, 0);
	
	DAVecGetArray(fda, user->lUcat,  &ucat);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lNu_t, &lnu_t);
	DAVecGetArray(da, user->lCs, &Cs);
	DAVecGetArray(da, user->lUstar, &ustar);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			lnu_t[k][j][i]=0;
			continue;
		}
		double ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	
		double Sabs = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );
		
		double filter  = pow( 1./aj[k][j][i],1./3.);
		lnu_t[k][j][i] = Cs[k][j][i] * pow ( filter, 2.0 ) * Sabs;

		if(les && wallfunction==2 && nvert[k][j][i]+nvert[k][j][i+1]+nvert[k][j][i-1]+nvert[k][j+1][i]+nvert[k][j-1][i]+nvert[k+1][j][i]+nvert[k-1][j][i]>0.1) lnu_t[k][j][i]=0;
		
		/*
		if( (user->bctype[0]==-1 && i==1) || (user->bctype[1]==-1 && i==mx-2) ){
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double sb = 0.5/aj[k][j][i]/area;
			double yplus = ustar[k][j][i] * sb * user->ren;
			lnu_t[k][j][i] = near_wall_eddy_viscosity(yplus) / user->ren;
		}
		
		if( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 && j==my-2) ){
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double sb = 0.5/aj[k][j][i]/area;
			double yplus = ustar[k][j][i] * sb * user->ren;
			lnu_t[k][j][i] = near_wall_eddy_viscosity(yplus) / user->ren;
		}
		*/
	}
	
	if(immersed && wallfunction==2) {
		DAVecRestoreArray(da, user->lNu_t, &lnu_t);
		
		DALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
		DALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
		
		DAVecGetArray(da, user->lNu_t, &lnu_t);
		
		for(int ibi=0; ibi<NumberOfBodies; ibi++)
		{
			extern IBMNodes *ibm_ptr;
			IBMNodes *ibm = &ibm_ptr[ibi];
			
			IBMListNode *current;
			current = user->ibmlist[ibi].head;
			while (current) {
				IBMInfo *ibminfo = &current->ibm_intp;
				
				current = current->next;
				double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
				int ni = ibminfo->cell;
				int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
				int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
				int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
				i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
				double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
				double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
				double nx = ibm->nf_x[ni], ny = ibm->nf_y[ni], nz = ibm->nf_z[ni];
				
				Cmpnts Ua, Uc;
				if (ni>=0) {
					Ua.x = ibm->u[ibm->nv1[ni]].x * cs1 + ibm->u[ibm->nv2[ni]].x * cs2 + ibm->u[ibm->nv3[ni]].x * cs3;
					Ua.y = ibm->u[ibm->nv1[ni]].y * cs1 + ibm->u[ibm->nv2[ni]].y * cs2 + ibm->u[ibm->nv3[ni]].y * cs3;
					Ua.z = ibm->u[ibm->nv1[ni]].z * cs1 + ibm->u[ibm->nv2[ni]].z * cs2 + ibm->u[ibm->nv3[ni]].z * cs3;
				}
				else {
					Ua.x = Ua.y = Ua.z = 0;
				}
				
				Uc.x = (ucat[kp1][jp1][ip1].x * sk1 + ucat[kp2][jp2][ip2].x * sk2 + ucat[kp3][jp3][ip3].x * sk3);
				Uc.y = (ucat[kp1][jp1][ip1].y * sk1 + ucat[kp2][jp2][ip2].y * sk2 + ucat[kp3][jp3][ip3].y * sk3);
				Uc.z = (ucat[kp1][jp1][ip1].z * sk1 + ucat[kp2][jp2][ip2].z * sk2 + ucat[kp3][jp3][ip3].z * sk3);
				
				//double yplus = ustar[k][j][i] * sb * user->ren;
				double nu = 1./user->ren;
				double nu_t_c = (lnu_t[kp1][jp1][ip1] * sk1 + lnu_t[kp2][jp2][ip2] * sk2 + lnu_t[kp3][jp3][ip3] * sk3);
				double eps=1.e-5;
				double dUc_ds = u_Cabot(nu, sc+eps, ustar[k][j][i], 0) - u_Cabot(nu, sc-eps, ustar[k][j][i], 0);
				dUc_ds /= (2.0 * eps);
				double f1 = ustar[k][j][i] * ustar[k][j][i];
				double f2 = (nu+nu_t_c) *  dUc_ds;
								
				double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
				double un = u_c * nx + v_c * ny + w_c * nz;
				double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
				double ut_mag_c = sqrt( ut*ut + vt*vt + wt*wt );
				double ut_mag_b  = u_Cabot(nu, sb, ustar[k][j][i], 0);
				
				//lnu_t[k][j][i] = nu_t_c * sb / sc;
				//lnu_t[k][j][i] = - nu + (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (ut_mag_c - ut_mag_b);
				//printf("nu_t = %f %f %f\n", ut_mag_c, ut_mag_b, ut_mag_c - (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (0 + nu));
				
				
				//ut_mag_b =  ut_mag_c - (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (lnu_t[k][j][i] + nu);
				//lnu_t[k][j][i] = near_wall_eddy_viscosity(yplus) / user->ren;
				/*
					Uc_t = Uc - (Uc . n) n
					
					
					f1 = ustar*ustar;
					f2 = shear_stress at C;
					(nu+nu_t)*(Uc - Ub)/(sc-sb) = (f1*(sc-sb) + f2*sb) / sc
					nu + nu_t = (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (Uc - Ub)
					nu_t = - nu + (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (Uc - Ub)
				*/
			};
		}
	}
	
	DAVecRestoreArray(fda, user->lUcat,  &ucat);
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	DAVecRestoreArray(da, user->lCs, &Cs);
	DAVecRestoreArray(da, user->lUstar, &ustar);
	
	DALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	DALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	
	DAVecGetArray(da, user->lNu_t, &lnu_t);
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


		if(flag) {
			lnu_t[k][j][i] = lnu_t[c][b][a];
		}
	}
	DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
	double lmax_norm=0, max_norm;
	int p;
	
	VecMax(user->lNu_t, &p, &lmax_norm);
	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Nu_t = %e\n", max_norm);
};


