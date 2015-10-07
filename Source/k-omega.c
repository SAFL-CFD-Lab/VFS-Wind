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

extern PetscInt immersed, NumberOfBodies, ti, tistart, wallfunction;
extern double find_utau_Cabot(double nu,  double u, double y, double guess, double dpdn);
extern PetscErrorCode MyKSPMonitor1(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy);
extern PetscInt inletprofile;

double solid=0.1;

const double C1=10.0, C2=2.0;	// C1=1~10, C2=2~5

double a1=0.31, beta_star = 0.09, sigma = 0.5, sigma_star = 0.5;
double beta1 = 0.075, beta2 = 0.0828;
double alpha1 = 5./9., alpha2 = 0.44;
const double sigma_k1 = 0.85, sigma_k2 = 1.0;
const double sigma_o1 = 0.50, sigma_o2 = 0.856;

#define wall_omega(ren, dist)  ( 6. / beta1 / ren / pow (dist, 2.0 ) )
// See Wilcox p.381

extern int rans; /* 1: Wilcox Low Re; 2: Wilcox High Re, 3: SST Menter */
extern double kappa;

void Get_alpha_beta_star(double ren, double K, double O, double *alpha, double *alpha_star, double *beta_star);

double Upwind(double W, double E, double a)
{
	if(a>0) return W;
	else return E;
};

void Get_alpha_beta_star(double ren, double K, double O, double *alpha, double *alpha_star, double *beta_star)
{
	
	if(rans==1) {
		const double Rb=8.0, Rk = 6.0, Rw=2.7, alpha0_star = beta1/3.0, alpha0 = 0.1;
		double ReT = K/O*ren;
		*alpha_star = (alpha0_star + ReT/Rk)/(1.0+ReT/Rk);
		*alpha = 5./9. * ( alpha0 + ReT/Rw) / ( 1.0 + ReT/Rw ) / (*alpha_star);
		*beta_star = 0.09 * ( 5./18. + pow( ReT/Rb, 4.0 ) ) / ( 1.0 + pow( ReT/Rb, 4.0 ) );
	}
	else *alpha = 5./9., *alpha_star=1., *beta_star = 0.09;
};

void Compute_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
{
	/*if( i_periodic && (i==0 || i==mx-2) ) {
		*dudc = ucat[k][j][1].x - ucat[k][j][mx-2].x;
		*dvdc = ucat[k][j][1].y - ucat[k][j][mx-2].y;
		*dwdc = ucat[k][j][1].z - ucat[k][j][mx-2].z;
	}
	else if( ii_periodic && i==0) {
		*dudc = ucat[k][j][i+1].x - ucat[k][j][-2].x;
		*dvdc = ucat[k][j][i+1].y - ucat[k][j][-2].y;
		*dwdc = ucat[k][j][i+1].z - ucat[k][j][-2].z;
	}
	else if( ii_periodic && i==mx-2) {
		*dudc = ucat[k][j][mx+1].x - ucat[k][j][i].x;
		*dvdc = ucat[k][j][mx+1].y - ucat[k][j][i].y;
		*dwdc = ucat[k][j][mx+1].z - ucat[k][j][i].z;
	}
	else */{
		*dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
		*dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
		*dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
	}
	
	if ((nvert[k][j+1][i])> solid || (nvert[k][j+1][i+1])> solid) {
		*dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
		*dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
		*dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if  ((nvert[k][j-1][i])> solid || (nvert[k][j-1][i+1])> solid) {
		*dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
		*dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
		*dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
		*dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
		*dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
		*dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
	}

	if ((nvert[k+1][j][i])> solid || (nvert[k+1][j][i+1])> solid) {
		*dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
		*dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
		*dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j][i+1])> solid) {
		*dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
		*dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
		*dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
		*dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
		*dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
		*dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
	}
}

void Compute_dscalar_i (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz )
{
	
	*dkdc = K[k][j][i+1] - K[k][j][i];
	
	
	if ((nvert[k][j+1][i])> solid || (nvert[k][j+1][i+1])> solid) {
		*dkde = (K[k][j  ][i+1] + K[k][j  ][i] - K[k][j-1][i+1] - K[k][j-1][i]) * 0.5;
	}
	else if  ((nvert[k][j-1][i])> solid || (nvert[k][j-1][i+1])> solid) {
		*dkde = (K[k][j+1][i+1] + K[k][j+1][i] - K[k][j  ][i+1] - K[k][j  ][i]) * 0.5;
	}
	else {
		*dkde = (K[k][j+1][i+1] + K[k][j+1][i] - K[k][j-1][i+1] - K[k][j-1][i]) * 0.25;
	}

	if ((nvert[k+1][j][i])> solid || (nvert[k+1][j][i+1])> solid) {
		*dkdz = (K[k  ][j][i+1] + K[k  ][j][i] - K[k-1][j][i+1] - K[k-1][j][i]) * 0.5;
	}
	else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j][i+1])> solid) {
		*dkdz = (K[k+1][j][i+1] + K[k+1][j][i] - K[k  ][j][i+1] - K[k  ][j][i]) * 0.5;
	}
	else {
		*dkdz = (K[k+1][j][i+1] + K[k+1][j][i] - K[k-1][j][i+1] - K[k-1][j][i]) * 0.25;
	}
}

void Compute_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
				
{
	if ((nvert[k][j][i+1])> solid || (nvert[k][j+1][i+1])> solid) {
		*dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
		*dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
		*dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k][j][i-1])> solid || (nvert[k][j+1][i-1])> solid) {
		*dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
		*dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
		*dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
		*dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
		*dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
		*dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}
	/*
	if( j_periodic && (j==0 || j==my-2) ) {
		*dude = ucat[k][1][i].x - ucat[k][my-2][i].x;
		*dvde = ucat[k][1][i].y - ucat[k][my-2][i].y;
		*dwde = ucat[k][1][i].z - ucat[k][my-2][i].z;
	}
	else if( jj_periodic && j==0) {
		*dude = ucat[k][j+1][i].x - ucat[k][-2][i].x;
		*dvde = ucat[k][j+1][i].y - ucat[k][-2][i].y;
		*dwde = ucat[k][j+1][i].z - ucat[k][-2][i].z;
	}
	else if( jj_periodic && j==my-2) {
		*dude = ucat[k][my+1][i].x - ucat[k][j][i].x;
		*dvde = ucat[k][my+1][i].y - ucat[k][j][i].y;
		*dwde = ucat[k][my+1][i].z - ucat[k][j][i].z;
	}
	else */{
		*dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
		*dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
		*dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;
	}
	
	if ((nvert[k+1][j][i])> solid || (nvert[k+1][j+1][i])> solid) {
		*dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
		*dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
		*dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j+1][i])> solid) {
		*dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
		*dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
		*dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
		*dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
		*dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
		*dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
	}
};

void Compute_dscalar_j (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz )
				
{
	if ((nvert[k][j][i+1])> solid || (nvert[k][j+1][i+1])> solid) {
		*dkdc = (K[k][j+1][i  ] + K[k][j][i  ] - K[k][j+1][i-1] - K[k][j][i-1]) * 0.5;
	}
	else if ((nvert[k][j][i-1])> solid || (nvert[k][j+1][i-1])> solid) {
		*dkdc = (K[k][j+1][i+1] + K[k][j][i+1] - K[k][j+1][i  ] - K[k][j][i  ]) * 0.5;
	}
	else {
		*dkdc = (K[k][j+1][i+1] + K[k][j][i+1] - K[k][j+1][i-1] - K[k][j][i-1]) * 0.25;
	}

	*dkde = K[k][j+1][i] - K[k][j][i];
	
	if ((nvert[k+1][j][i])> solid || (nvert[k+1][j+1][i])> solid) {
		*dkdz = (K[k  ][j+1][i] + K[k  ][j][i] - K[k-1][j+1][i] - K[k-1][j][i]) * 0.5;
	}
	else if ((nvert[k-1][j][i])> solid || (nvert[k-1][j+1][i])> solid) {
		*dkdz = (K[k+1][j+1][i] + K[k+1][j][i] - K[k  ][j+1][i] - K[k  ][j][i]) * 0.5;
	}
	else {
		*dkdz = (K[k+1][j+1][i] + K[k+1][j][i] - K[k-1][j+1][i] - K[k-1][j][i]) * 0.25;
	}
};

void Compute_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
				
{
	if ((nvert[k][j][i+1])> solid || (nvert[k+1][j][i+1])> solid) {
		*dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
		*dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
		*dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k][j][i-1])> solid || (nvert[k+1][j][i-1])> solid) {
		*dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
		*dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
		*dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
		*dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
		*dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
		*dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}

	if ((nvert[k][j+1][i])> solid || (nvert[k+1][j+1][i])> solid) {
		*dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
		*dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
		*dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if ((nvert[k][j-1][i])> solid || (nvert[k+1][j-1][i])> solid) {
		*dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
		*dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
		*dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
		*dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
		*dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
		*dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
	}
	/*
	if( k_periodic && (k==0 || k==mz-2) ) {
		*dudz = ucat[1][j][i].x - ucat[mz-2][j][i].x;
		*dvdz = ucat[1][j][i].y - ucat[mz-2][j][i].y;
		*dwdz = ucat[1][j][i].z - ucat[mz-2][j][i].z;
	}
	else if( kk_periodic && k==0) {
		*dudz = ucat[1][j][i].x - ucat[-2][j][i].x;
		*dvdz = ucat[1][j][i].y - ucat[-2][j][i].y;
		*dwdz = ucat[1][j][i].z - ucat[-2][j][i].z;
	}
	else if( kk_periodic && k==mz-2) {
		*dudz = ucat[mz+1][j][i].x - ucat[mz-2][j][i].x;
		*dvdz = ucat[mz+1][j][i].y - ucat[mz-2][j][i].y;
		*dwdz = ucat[mz+1][j][i].z - ucat[mz-2][j][i].z;
	}
	else*/ {
		*dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
		*dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
		*dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;
	}
}

void Compute_dscalar_k (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz )
{
	if ((nvert[k][j][i+1])> solid || (nvert[k+1][j][i+1])> solid) {
		*dkdc = (K[k+1][j][i  ] + K[k][j][i  ] - K[k+1][j][i-1] - K[k][j][i-1]) * 0.5;
	}
	else if ((nvert[k][j][i-1])> solid || (nvert[k+1][j][i-1])> solid) {
		*dkdc = (K[k+1][j][i+1] + K[k][j][i+1] - K[k+1][j][i  ] - K[k][j][i  ]) * 0.5;
	}
	else {
		*dkdc = (K[k+1][j][i+1] + K[k][j][i+1] - K[k+1][j][i-1] - K[k][j][i-1]) * 0.25;
	}

	if ((nvert[k][j+1][i])> solid || (nvert[k+1][j+1][i])> solid) {
		*dkde = (K[k+1][j  ][i] + K[k][j  ][i] - K[k+1][j-1][i] - K[k][j-1][i]) * 0.5;
	}
	else if ((nvert[k][j-1][i])> solid || (nvert[k+1][j-1][i])> solid) {
		*dkde = (K[k+1][j+1][i] + K[k][j+1][i] - K[k+1][j  ][i] - K[k][j  ][i]) * 0.5;
	}
	else {
		*dkde = (K[k+1][j+1][i] + K[k][j+1][i] - K[k+1][j-1][i] - K[k][j-1][i]) * 0.25;
	}

	*dkdz = K[k+1][j][i] - K[k][j][i];
}

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
{
	if ((nvert[k][j][i+1])> solid || (!i_periodic &&  !ii_periodic && i==mx-2) ) {
		*dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
		*dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
		*dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
	}
	else if ((nvert[k][j][i-1])> solid || (!i_periodic &&  !ii_periodic && i==1) ) {
		*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
		*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
		*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
	}
	else {
		/*if(i_periodic && i==1) {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][mx-2].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][mx-2].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][mx-2].z ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dudc = ( ucat[k][j][1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][1].z - ucat[k][j][i-1].z ) * 0.5;
		}
		else if(ii_periodic && i==1) {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][-2].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][-2].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][-2].z ) * 0.5;
		}
		else if(ii_periodic && i==mx-2) {
			*dudc = ( ucat[k][j][mx+1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][mx+1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][mx+1].z - ucat[k][j][i-1].z ) * 0.5;
		}
		else*/ {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
		}
	}

	if ((nvert[k][j+1][i])> solid || (!j_periodic &&  !jj_periodic && j==my-2) ) {
		*dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
		*dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
		*dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
	}
	else if ((nvert[k][j-1][i])> solid || (!j_periodic &&  !jj_periodic && j==1) ) {
		*dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
		*dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
		*dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );
	}
	else {
		/*if(j_periodic && j==1) {
			*dude = ( ucat[k][j+1][i].x - ucat[k][my-2][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][my-2][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][my-2][i].z ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dude = ( ucat[k][1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
		else if(jj_periodic && j==1) {
			*dude = ( ucat[k][j+1][i].x - ucat[k][-2][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][-2][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][-2][i].z ) * 0.5;
		}
		else if(jj_periodic && j==my-2) {
			*dude = ( ucat[k][my+1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][my+1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][my+1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
		else */{
			*dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
	}

	if ((nvert[k+1][j][i])> solid || ( !k_periodic &&  !kk_periodic && k==mz-2) ) {
		*dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
		*dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
		*dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
	}
	else if ((nvert[k-1][j][i])> solid || (!k_periodic &&  !kk_periodic && k==0) ) {
		*dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
		*dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
		*dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
	}
	else {
		/*if(k_periodic && k==1) {
			*dudz = ( ucat[k+1][j][i].x - ucat[mz-2][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[mz-2][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[mz-2][j][i].z ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dudz = ( ucat[1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
		else if(kk_periodic && k==1) {
			*dudz = ( ucat[k+1][j][i].x - ucat[-2][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[-2][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[-2][j][i].z ) * 0.5;
		}
		else if(kk_periodic && k==mz-2) {
			*dudz = ( ucat[mz+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[mz+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[mz+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
		else*/ {
			*dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
	}
}

void Compute_dkdo_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts2 ***komega, PetscReal ***nvert, 
				double *dkdc, double *dodc,  double *dkde, double *dode,  double *dkdz, double *dodz )
{
  /*if ((nvert[k][j][i+1])> solid) {
		*dkdc = ( komega[k][j][i].x - komega[k][j][i-1].x );
		*dodc = ( komega[k][j][i].y - komega[k][j][i-1].y );
	}
	else if ((nvert[k][j][i-1])> solid) {
		*dkdc = ( komega[k][j][i+1].x - komega[k][j][i].x );
		*dodc = ( komega[k][j][i+1].y - komega[k][j][i].y );
	}
	else*/ {/*
		if(i_periodic && i==1) {
			*dkdc = ( komega[k][j][i+1].x - komega[k][j][mx-2].x ) * 0.5;
			*dodc = ( komega[k][j][i+1].y - komega[k][j][mx-2].y ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dkdc = ( komega[k][j][1].x - komega[k][j][i-1].x ) * 0.5;
			*dodc = ( komega[k][j][1].y - komega[k][j][i-1].y ) * 0.5;
		}
		else if(ii_periodic && i==1) {
			*dkdc = ( komega[k][j][i+1].x - komega[k][j][-2].x ) * 0.5;
			*dodc = ( komega[k][j][i+1].y - komega[k][j][-2].y ) * 0.5;
		}
		else if(ii_periodic && i==mx-2) {
			*dkdc = ( komega[k][j][mx+1].x - komega[k][j][i-1].x ) * 0.5;
			*dodc = ( komega[k][j][mx+1].y - komega[k][j][i-1].y ) * 0.5;
		}
		else*/ {
			*dkdc = ( komega[k][j][i+1].x - komega[k][j][i-1].x ) * 0.5;
			*dodc = ( komega[k][j][i+1].y - komega[k][j][i-1].y ) * 0.5;
		}
	}

	/*if ((nvert[k][j+1][i])> solid) {
		*dkde = ( komega[k][j][i].x - komega[k][j-1][i].x );
		*dode = ( komega[k][j][i].y - komega[k][j-1][i].y );
	}
	else if ((nvert[k][j-1][i])> solid) {
		*dkde = ( komega[k][j+1][i].x - komega[k][j][i].x );
		*dode = ( komega[k][j+1][i].y - komega[k][j][i].y );
	}
	else*/ {
	  /*if(j_periodic && j==1) {
			*dkde = ( komega[k][j+1][i].x - komega[k][my-2][i].x ) * 0.5;
			*dode = ( komega[k][j+1][i].y - komega[k][my-2][i].y ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dkde = ( komega[k][1][i].x - komega[k][j-1][i].x ) * 0.5;
			*dode = ( komega[k][1][i].y - komega[k][j-1][i].y ) * 0.5;
		}
		else if(jj_periodic && j==1) {
			*dkde = ( komega[k][j+1][i].x - komega[k][-2][i].x ) * 0.5;
			*dode = ( komega[k][j+1][i].y - komega[k][-2][i].y ) * 0.5;
		}
		else if(jj_periodic && j==my-2) {
			*dkde = ( komega[k][my+1][i].x - komega[k][j-1][i].x ) * 0.5;
			*dode = ( komega[k][my+1][i].y - komega[k][j-1][i].y ) * 0.5;
		}
		else*/ {
			*dkde = ( komega[k][j+1][i].x - komega[k][j-1][i].x ) * 0.5;
			*dode = ( komega[k][j+1][i].y - komega[k][j-1][i].y ) * 0.5;
		}
	}

	/*if ((nvert[k+1][j][i])> solid) {
		*dkdz = ( komega[k][j][i].x - komega[k-1][j][i].x );
		*dodz = ( komega[k][j][i].y - komega[k-1][j][i].y );
	}
	else if ((nvert[k-1][j][i])> solid) {
		*dkdz = ( komega[k+1][j][i].x - komega[k][j][i].x );
		*dodz = ( komega[k+1][j][i].y - komega[k][j][i].y );
	}
	else*/ {
	  /*if(k_periodic && k==1) {
			*dkdz = ( komega[k+1][j][i].x - komega[mz-2][j][i].x ) * 0.5;
			*dodz = ( komega[k+1][j][i].y - komega[mz-2][j][i].y ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dkdz = ( komega[1][j][i].x - komega[k-1][j][i].x ) * 0.5;
			*dodz = ( komega[1][j][i].y - komega[k-1][j][i].y ) * 0.5;
		}
		else if(kk_periodic && k==1) {
			*dkdz = ( komega[k+1][j][i].x - komega[-2][j][i].x ) * 0.5;
			*dodz = ( komega[k+1][j][i].y - komega[-2][j][i].y ) * 0.5;
		}
		else if(kk_periodic && k==mz-2) {
			*dkdz = ( komega[mz+1][j][i].x - komega[k-1][j][i].x ) * 0.5;
			*dodz = ( komega[mz+1][j][i].y - komega[k-1][j][i].y ) * 0.5;
		}
		else*/ {
			*dkdz = ( komega[k+1][j][i].x - komega[k-1][j][i].x ) * 0.5;
			*dodz = ( komega[k+1][j][i].y - komega[k-1][j][i].y ) * 0.5;
		}
	}
}

void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz)
{
	if ((nvert[k][j][i+1])> solid) {
		*dkdc = ( K[k][j][i] - K[k][j][i-1] );
	}
	else if ((nvert[k][j][i-1])> solid) {
		*dkdc = ( K[k][j][i+1] - K[k][j][i] );
	}
	else {
		/*if(i_periodic && i==1) {
			*dkdc = ( K[k][j][i+1] - K[k][j][mx-2] ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dkdc = ( K[k][j][1] - K[k][j][i-1] ) * 0.5;
		}
		else if(ii_periodic && i==1) {
			*dkdc = ( K[k][j][i+1] - K[k][j][-2] ) * 0.5;
		}
		else if(ii_periodic && i==mx-2) {
			*dkdc = ( K[k][j][mx+1] - K[k][j][i-1] ) * 0.5;
		}
		else*/ {
			*dkdc = ( K[k][j][i+1] - K[k][j][i-1] ) * 0.5;
		}
	}

	if ((nvert[k][j+1][i])> solid) {
		*dkde = ( K[k][j][i] - K[k][j-1][i] );
	}
	else if ((nvert[k][j-1][i])> solid) {
		*dkde = ( K[k][j+1][i] - K[k][j][i] );
	}
	else {
		/*if(j_periodic && j==1) {
			*dkde = ( K[k][j+1][i] - K[k][my-2][i] ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dkde = ( K[k][1][i] - K[k][j-1][i] ) * 0.5;
		}
		else if(jj_periodic && j==1) {
			*dkde = ( K[k][j+1][i] - K[k][-2][i] ) * 0.5;
		}
		else if(jj_periodic && j==my-2) {
			*dkde = ( K[k][my+1][i] - K[k][j-1][i] ) * 0.5;
		}
		else*/ {
			*dkde = ( K[k][j+1][i] - K[k][j-1][i] ) * 0.5;
		}
	}

	if ((nvert[k+1][j][i])> solid) {
		*dkdz = ( K[k][j][i] - K[k-1][j][i] );
	}
	else if ((nvert[k-1][j][i])> solid) {
		*dkdz = ( K[k+1][j][i] - K[k][j][i] );
	}
	else {
		/*if(k_periodic && k==1) {
			*dkdz = ( K[k+1][j][i] - K[mz-2][j][i] ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dkdz = ( K[1][j][i] - K[k-1][j][i] ) * 0.5;
		}
		else if(kk_periodic && k==1) {
			*dkdz = ( K[k+1][j][i] - K[-2][j][i] ) * 0.5;
		}
		else if(kk_periodic && k==mz-2) {
			*dkdz = ( K[mz+1][j][i] - K[k-1][j][i] ) * 0.5;
		}
		else */{
			*dkdz = ( K[k+1][j][i] - K[k-1][j][i] ) * 0.5;
		}
	}	
}

void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz )
{
	*du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
	*du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
	*du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
	*dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
	*dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
	*dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
	*dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
	*dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;	
	*dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;
};

void Compute_dkdo_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dkdc, double dodc, double dkde, double dode, double dkdz, double dodz,
					double *dk_dx, double *do_dx, double *dk_dy, double *do_dy, double *dk_dz, double *do_dz)
{
	*dk_dx = (dkdc * csi0 + dkde * eta0 + dkdz * zet0) * ajc;
	*dk_dy = (dkdc * csi1 + dkde * eta1 + dkdz * zet1) * ajc;
	*dk_dz = (dkdc * csi2 + dkde * eta2 + dkdz * zet2) * ajc;
	*do_dx = (dodc * csi0 + dode * eta0 + dodz * zet0) * ajc;
	*do_dy = (dodc * csi1 + dode * eta1 + dodz * zet1) * ajc;
	*do_dz = (dodc * csi2 + dode * eta2 + dodz * zet2) * ajc;
};

void Compute_dscalar_dxyz ( double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
							double dkdc, double dkde, double dkdz, double *dk_dx, double *dk_dy, double *dk_dz)
{
	*dk_dx = (dkdc * csi0 + dkde * eta0 + dkdz * zet0) * ajc;
	*dk_dy = (dkdc * csi1 + dkde * eta1 + dkdz * zet1) * ajc;
	*dk_dz = (dkdc * csi2 + dkde * eta2 + dkdz * zet2) * ajc;
};

void RHS_K_Omega(UserCtx *user, Vec KOmega_RHS)
{
	DA		da = user->da, fda = user->fda, fda2 = user->fda2;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	
	PetscReal	***aj;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	Cmpnts2	***K_Omega, ***K_Omega_o, ***komega_rhs;
	Cmpnts	***ucont, ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
	PetscReal	***nvert, ***distance, ***lf1, ***lnu_t;
 
	Vec Fp1, Fp2, Fp3;
	Vec Visc1, Visc2, Visc3;
	Cmpnts2 ***fp1, ***fp2, ***fp3;
	Cmpnts2 ***visc1, ***visc2, ***visc3;
	
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal	***iaj, ***jaj, ***kaj, ***rho, ***mu;

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
	
	VecDuplicate(user->lK_Omega, &Fp1);
	VecDuplicate(user->lK_Omega, &Fp2);
	VecDuplicate(user->lK_Omega, &Fp3);
	VecDuplicate(user->lK_Omega, &Visc1);
	VecDuplicate(user->lK_Omega, &Visc2);
	VecDuplicate(user->lK_Omega, &Visc3);
	
	VecSet(Fp1,0);
	VecSet(Fp2,0);
	VecSet(Fp3,0);
	VecSet(Visc1,0);
	VecSet(Visc2,0);
	VecSet(Visc3,0);
	
	
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	
	//DAVecGetArray(da, user->lSrans, &lS);
	if(rans==3) DAVecGetArray(da, user->lF1, &lf1);
	DAVecGetArray(da, user->lNu_t, &lnu_t);
	
	DAVecGetArray(fda, user->lUcont, &ucont);
	DAVecGetArray(fda, user->lUcat,  &ucat);
	
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
	
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->Distance, &distance);
	
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);
	
	DAVecGetArray(fda2, user->lK_Omega, &K_Omega);
	DAVecGetArray(fda2, user->lK_Omega_o, &K_Omega_o);
	DAVecGetArray(fda2, KOmega_RHS, &komega_rhs);
		
	DAVecGetArray(fda2, Fp1, &fp1);
	DAVecGetArray(fda2, Fp2, &fp2);
	DAVecGetArray(fda2, Fp3, &fp3);
	DAVecGetArray(fda2, Visc1, &visc1);
	DAVecGetArray(fda2, Visc2, &visc2);
	DAVecGetArray(fda2, Visc3, &visc3);		
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs-1; i<lxe; i++) {
		double csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		double eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		double zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
		double ajc = iaj[k][j][i];
		
		double g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		double g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		double g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
		
		double dkdc, dkde, dkdz;
		double dodc, dode, dodz;
		double nu_t;
		
		
		dkdc = K_Omega[k][j][i+1].x - K_Omega[k][j][i].x;
		dodc = K_Omega[k][j][i+1].y - K_Omega[k][j][i].y;

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dkde = (K_Omega[k][j  ][i+1].x + K_Omega[k][j  ][i].x - K_Omega[k][j-1][i+1].x - K_Omega[k][j-1][i].x) * 0.5;
			dode = (K_Omega[k][j  ][i+1].y + K_Omega[k][j  ][i].y - K_Omega[k][j-1][i+1].y - K_Omega[k][j-1][i].y) * 0.5;
		}
		else if  ((nvert[k][j-1][i])> 1.1 || (nvert[k][j-1][i+1])> 1.1) {
			dkde = (K_Omega[k][j+1][i+1].x + K_Omega[k][j+1][i].x - K_Omega[k][j  ][i+1].x - K_Omega[k][j  ][i].x) * 0.5;
			dode = (K_Omega[k][j+1][i+1].y + K_Omega[k][j+1][i].y - K_Omega[k][j  ][i+1].y - K_Omega[k][j  ][i].y) * 0.5;
		}
		else {
			dkde = (K_Omega[k][j+1][i+1].x + K_Omega[k][j+1][i].x - K_Omega[k][j-1][i+1].x - K_Omega[k][j-1][i].x) * 0.25;
			dode = (K_Omega[k][j+1][i+1].y + K_Omega[k][j+1][i].y - K_Omega[k][j-1][i+1].y - K_Omega[k][j-1][i].y) * 0.25;
		}	  

		if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
			dkdz = (K_Omega[k  ][j][i+1].x + K_Omega[k  ][j][i].x - K_Omega[k-1][j][i+1].x - K_Omega[k-1][j][i].x) * 0.5;
			dodz = (K_Omega[k  ][j][i+1].y + K_Omega[k  ][j][i].y - K_Omega[k-1][j][i+1].y - K_Omega[k-1][j][i].y) * 0.5;
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j][i+1])> 1.1) {
			dkdz = (K_Omega[k+1][j][i+1].x + K_Omega[k+1][j][i].x - K_Omega[k  ][j][i+1].x - K_Omega[k  ][j][i].x) * 0.5;
			dodz = (K_Omega[k+1][j][i+1].y + K_Omega[k+1][j][i].y - K_Omega[k  ][j][i+1].y - K_Omega[k  ][j][i].y) * 0.5;
		}
		else {
			dkdz = (K_Omega[k+1][j][i+1].x + K_Omega[k+1][j][i].x - K_Omega[k-1][j][i+1].x - K_Omega[k-1][j][i].x) * 0.25;
			dodz = (K_Omega[k+1][j][i+1].y + K_Omega[k+1][j][i].y - K_Omega[k-1][j][i+1].y - K_Omega[k-1][j][i].y) * 0.25;
		}
		
		double wL=1./aj[k][j][i], wR=1./aj[k][j][i+1];
		double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k][j][i+1].x*wR ) / (wL+wR);	Km_o = PetscMax ( Km_o, 0 );
		double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k][j][i+1].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
		double sigma_k=0, sigma_o=0;
		
		if(rans==3) {
			nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j][i+1]);
			double f1 = 0.5 * ( lf1[k][j][i] + lf1[k][j][i+1] );
			sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
			sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
		}
		else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j][i+1]);
		
		nu_t = PetscMax(nu_t, 0);

		if( nvert[k][j][i]+nvert[k][j][i+1]>0.1 || i==0 || i==mx-2 || periodic ) {
			fp1[k][j][i].x = -ucont[k][j][i].x * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j][i+1].x, ucont[k][j][i].x);
			fp1[k][j][i].y = -ucont[k][j][i].x * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j][i+1].y, ucont[k][j][i].x);
		}
		else {
			fp1[k][j][i].x = -ucont[k][j][i].x * weno3 ( K_Omega[k][j][i-1].x, K_Omega[k][j][i].x, K_Omega[k][j][i+1].x, K_Omega[k][j][i+2].x, ucont[k][j][i].x );
			fp1[k][j][i].y = -ucont[k][j][i].x * weno3 ( K_Omega[k][j][i-1].y, K_Omega[k][j][i].y, K_Omega[k][j][i+1].y, K_Omega[k][j][i+2].y, ucont[k][j][i].x );
		}
		
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || i==0) nu=mu[k][j][i+1];
			else if(nvert[k][j][i+1]>0.1 || i==mx-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j][i+1] );
		}
		
		if(rans==3) {
			visc1[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
			visc1[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
		}
		else  {
			visc1[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
			visc1[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys-1; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		double eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		double zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
		double ajc = jaj[k][j][i];
		
		double g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
		double g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		double g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
		
		double dkdc, dkde, dkdz;
		double dodc, dode, dodz;
		double nu_t;		
		
		if ((nvert[k][j][i+1])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dkdc = (K_Omega[k][j+1][i  ].x + K_Omega[k][j][i  ].x - K_Omega[k][j+1][i-1].x - K_Omega[k][j][i-1].x) * 0.5;
			dodc = (K_Omega[k][j+1][i  ].y + K_Omega[k][j][i  ].y - K_Omega[k][j+1][i-1].y - K_Omega[k][j][i-1].y) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k][j+1][i-1])> 1.1) {
			dkdc = (K_Omega[k][j+1][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k][j+1][i  ].x - K_Omega[k][j][i  ].x) * 0.5;
			dodc = (K_Omega[k][j+1][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k][j+1][i  ].y - K_Omega[k][j][i  ].y) * 0.5;
		}
		else {
			dkdc = (K_Omega[k][j+1][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k][j+1][i-1].x - K_Omega[k][j][i-1].x) * 0.25;
			dodc = (K_Omega[k][j+1][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k][j+1][i-1].y - K_Omega[k][j][i-1].y) * 0.25;
		}

		dkde = K_Omega[k][j+1][i].x - K_Omega[k][j][i].x;
		dode = K_Omega[k][j+1][i].y - K_Omega[k][j][i].y;

		if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dkdz = (K_Omega[k  ][j+1][i].x + K_Omega[k  ][j][i].x - K_Omega[k-1][j+1][i].x - K_Omega[k-1][j][i].x) * 0.5;
			dodz = (K_Omega[k  ][j+1][i].y + K_Omega[k  ][j][i].y - K_Omega[k-1][j+1][i].y - K_Omega[k-1][j][i].y) * 0.5;
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j+1][i])> 1.1) {
			dkdz = (K_Omega[k+1][j+1][i].x + K_Omega[k+1][j][i].x - K_Omega[k  ][j+1][i].x - K_Omega[k  ][j][i].x) * 0.5;
			dodz = (K_Omega[k+1][j+1][i].y + K_Omega[k+1][j][i].y - K_Omega[k  ][j+1][i].y - K_Omega[k  ][j][i].y) * 0.5;
		}
		else {
			dkdz = (K_Omega[k+1][j+1][i].x + K_Omega[k+1][j][i].x - K_Omega[k-1][j+1][i].x - K_Omega[k-1][j][i].x) * 0.25;
			dodz = (K_Omega[k+1][j+1][i].y + K_Omega[k+1][j][i].y - K_Omega[k-1][j+1][i].y - K_Omega[k-1][j][i].y) * 0.25;
		}
		
		double wL=1./aj[k][j][i], wR=1./aj[k][j+1][i];
		double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k][j+1][i].x*wR ) / (wL+wR);	Km_o = PetscMax ( Km_o, 0 );
		double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k][j+1][i].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
		double sigma_k=0, sigma_o=0;
		
		if(rans==3) {
			nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i]);
			double f1 = 0.5 * ( lf1[k][j][i] + lf1[k][j+1][i] );
			
			sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
			sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
		}
		else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i]);
		
		nu_t = PetscMax(nu_t, 0);
		
		if( nvert[k][j][i]+nvert[k][j+1][i]>0.1 || j==0 || j==my-2 || periodic ) {
			fp2[k][j][i].x = -ucont[k][j][i].y * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j+1][i].x, ucont[k][j][i].y);
			fp2[k][j][i].y = -ucont[k][j][i].y * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j+1][i].y, ucont[k][j][i].y);
		}
		else {
			fp2[k][j][i].x = -ucont[k][j][i].y * weno3 ( K_Omega[k][j-1][i].x, K_Omega[k][j][i].x, K_Omega[k][j+1][i].x, K_Omega[k][j+2][i].x, ucont[k][j][i].y );
			fp2[k][j][i].y = -ucont[k][j][i].y * weno3 ( K_Omega[k][j-1][i].y, K_Omega[k][j][i].y, K_Omega[k][j+1][i].y, K_Omega[k][j+2][i].y, ucont[k][j][i].y );
		}
		
		
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || j==0) nu=mu[k][j+1][i];
			else if(nvert[k][j+1][i]>0.1 || j==my-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j+1][i] );
		}
		
		if(rans==3) {
			visc2[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
			visc2[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
		}
		else  {
			visc2[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
			visc2[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
		}
	}
	
	for (k=lzs-1; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		double eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		double zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		double ajc = kaj[k][j][i];
		
		double g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
		double g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
		double g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
		
		double dkdc, dkde, dkdz;
		double dodc, dode, dodz;
		double nu_t;
		
		if ((nvert[k][j][i+1])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
			dkdc = (K_Omega[k+1][j][i  ].x + K_Omega[k][j][i  ].x - K_Omega[k+1][j][i-1].x - K_Omega[k][j][i-1].x) * 0.5;
			dodc = (K_Omega[k+1][j][i  ].y + K_Omega[k][j][i  ].y - K_Omega[k+1][j][i-1].y - K_Omega[k][j][i-1].y) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k+1][j][i-1])> 1.1) {
			dkdc = (K_Omega[k+1][j][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k+1][j][i  ].x - K_Omega[k][j][i  ].x) * 0.5;
			dodc = (K_Omega[k+1][j][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k+1][j][i  ].y - K_Omega[k][j][i  ].y) * 0.5;
		}
		else {
			dkdc = (K_Omega[k+1][j][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k+1][j][i-1].x - K_Omega[k][j][i-1].x) * 0.25;
			dodc = (K_Omega[k+1][j][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k+1][j][i-1].y - K_Omega[k][j][i-1].y) * 0.25;
		}

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dkde = (K_Omega[k+1][j  ][i].x + K_Omega[k][j  ][i].x - K_Omega[k+1][j-1][i].x - K_Omega[k][j-1][i].x) * 0.5;
			dode = (K_Omega[k+1][j  ][i].y + K_Omega[k][j  ][i].y - K_Omega[k+1][j-1][i].y - K_Omega[k][j-1][i].y) * 0.5;
		}
		else if ((nvert[k][j-1][i])> 1.1 || (nvert[k+1][j-1][i])> 1.1) {
			dkde = (K_Omega[k+1][j+1][i].x + K_Omega[k][j+1][i].x - K_Omega[k+1][j  ][i].x - K_Omega[k][j  ][i].x) * 0.5;
			dode = (K_Omega[k+1][j+1][i].y + K_Omega[k][j+1][i].y - K_Omega[k+1][j  ][i].y - K_Omega[k][j  ][i].y) * 0.5;
		}
		else {
			dkde = (K_Omega[k+1][j+1][i].x + K_Omega[k][j+1][i].x - K_Omega[k+1][j-1][i].x - K_Omega[k][j-1][i].x) * 0.25;
			dode = (K_Omega[k+1][j+1][i].y + K_Omega[k][j+1][i].y - K_Omega[k+1][j-1][i].y - K_Omega[k][j-1][i].y) * 0.25;
		}

		dkdz = K_Omega[k+1][j][i].x - K_Omega[k][j][i].x;
		dodz = K_Omega[k+1][j][i].y - K_Omega[k][j][i].y;
		
		double wL=1./aj[k][j][i], wR=1./aj[k+1][j][i];
		double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k+1][j][i].x*wR ) / (wL+wR);	Km_o = PetscMax ( Km_o, 0 );
		double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k+1][j][i].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
		double sigma_k=0, sigma_o=0;
		
		if(rans==3) {
			double f1 = 0.5 * ( lf1[k][j][i] + lf1[k+1][j][i] );
			nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k+1][j][i]);
			
			sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
			sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
		}
		else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k+1][j][i]);

		nu_t = PetscMax(nu_t, 0);
		
		if( nvert[k][j][i]+nvert[k+1][j][i]>0.1 || k==0 || k==mz-2 || periodic ) {
			fp3[k][j][i].x = -ucont[k][j][i].z * Upwind ( K_Omega[k][j][i].x, K_Omega[k+1][j][i].x, ucont[k][j][i].z);
			fp3[k][j][i].y = -ucont[k][j][i].z * Upwind ( K_Omega[k][j][i].y, K_Omega[k+1][j][i].y, ucont[k][j][i].z);
		}
		else {
			fp3[k][j][i].x = -ucont[k][j][i].z * weno3 ( K_Omega[k-1][j][i].x, K_Omega[k][j][i].x, K_Omega[k+1][j][i].x, K_Omega[k+2][j][i].x, ucont[k][j][i].z );
			fp3[k][j][i].y = -ucont[k][j][i].z * weno3 ( K_Omega[k-1][j][i].y, K_Omega[k][j][i].y, K_Omega[k+1][j][i].y, K_Omega[k+2][j][i].y, ucont[k][j][i].z );
		}
		
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || k==0) nu=mu[k+1][j][i];
			else if(nvert[k+1][j][i]>0.1 || k==mz-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k+1][j][i] );
		}
		
		if(rans==3) {
			visc3[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
			visc3[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
		}
		else  {
			visc3[k][j][i].x = (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
			visc3[k][j][i].y = (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
		}
	}
	
	DAVecRestoreArray(fda2, Fp1, &fp1);
	DAVecRestoreArray(fda2, Fp2, &fp2);
	DAVecRestoreArray(fda2, Fp3, &fp3);
	DAVecRestoreArray(fda2, Visc1, &visc1);
	DAVecRestoreArray(fda2, Visc2, &visc2);
	DAVecRestoreArray(fda2, Visc3, &visc3);

	DALocalToLocalBegin(fda2, Fp1, INSERT_VALUES, Fp1);
	DALocalToLocalEnd(fda2, Fp1, INSERT_VALUES, Fp1);
	
	DALocalToLocalBegin(fda2, Fp2, INSERT_VALUES, Fp2);
	DALocalToLocalEnd(fda2, Fp2, INSERT_VALUES, Fp2);
	
	DALocalToLocalBegin(fda2, Fp3, INSERT_VALUES, Fp3);
	DALocalToLocalEnd(fda2, Fp3, INSERT_VALUES, Fp3);
	
	DALocalToLocalBegin(fda2, Visc1, INSERT_VALUES, Visc1);
	DALocalToLocalEnd(fda2, Visc1, INSERT_VALUES, Visc1);
	
	DALocalToLocalBegin(fda2, Visc2, INSERT_VALUES, Visc2);
	DALocalToLocalEnd(fda2, Visc2, INSERT_VALUES, Visc2);
	
	DALocalToLocalBegin(fda2, Visc3, INSERT_VALUES, Visc3);
	DALocalToLocalEnd(fda2, Visc3, INSERT_VALUES, Visc3);
	
	DAVecGetArray(fda2, Fp1, &fp1);
	DAVecGetArray(fda2, Fp2, &fp2);
	DAVecGetArray(fda2, Fp3, &fp3);
	DAVecGetArray(fda2, Visc1, &visc1);
	DAVecGetArray(fda2, Visc2, &visc2);
	DAVecGetArray(fda2, Visc3, &visc3);
	
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
			fp1[k][j][i] = fp1[c][b][a];
			fp2[k][j][i] = fp2[c][b][a];
			fp3[k][j][i] = fp3[c][b][a];
			visc1[k][j][i] = visc1[c][b][a];
			visc2[k][j][i] = visc2[c][b][a];
			visc3[k][j][i] = visc3[c][b][a];
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ( i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1 ) {
			komega_rhs[k][j][i].x = komega_rhs[k][j][i].y = 0;
			continue;
		}
			
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		double Km = K_Omega[k][j][i].x, Om = K_Omega[k][j][i].y;
		double Km_o = K_Omega_o[k][j][i].x, Om_o = K_Omega_o[k][j][i].y;
		
		Km_o = PetscMax ( Km_o, 0 );
		Om_o = PetscMax ( Om_o, 1.e-5 );

		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double ajc = aj[k][j][i];//, d = distance[k][j][i];
		
		double dkdc, dodc, dkde, dode, dkdz, dodz;
		double dk_dx, dk_dy, dk_dz, do_dx,  do_dy, do_dz;
		
		Compute_dkdo_center (i, j, k, mx, my, mz, K_Omega, nvert,  &dkdc,  &dodc,  &dkde,  &dode,  &dkdz,  &dodz );
		Compute_dkdo_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dkdc, dodc, dkde, dode, dkdz, dodz, &dk_dx, &do_dx, &dk_dy, &do_dy, &dk_dz, &do_dz);
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
		double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
		
		double nu_t, inv_nu_t, f1=1;
		
		if(rans==3) {
			nu_t = lnu_t[k][j][i];
			f1 = lf1[k][j][i];
			/*
			if( !i_periodic && !ii_periodic && (i==mx-2 ) ) f1=1;
			if( !j_periodic && !jj_periodic && (j==my-2 ) ) f1=1;
			if( !k_periodic && !kk_periodic && (k==mz-2 ) ) f1=1;
			*/
			if( nvert[k][j][i-1]+nvert[k][j][i+1]>0.1 ) f1=1;
			if( nvert[k][j-1][i]+nvert[k][j+1][i]>0.1 ) f1=1;
			if( nvert[k-1][j][i]+nvert[k+1][j][i]>0.1 ) f1=1;
		}
		else nu_t = lnu_t[k][j][i];
		//else nu_t = Km_o/Om_o;
		
		if( nu_t<1.e-20 ) inv_nu_t = 0;
		else inv_nu_t=1./nu_t;
		
		nu_t = PetscMax(nu_t, 0);
		
		// Wilcox pp. 75
		//double Strain[3][3] = { Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz };
		/*double Rotation[3][3] = { 0, 				0.5*(du_dy-dv_dx), 	0.5*(du_dz-dw_dx),
							-0.5*(du_dy-dv_dx), 	0, 				0.5*(dv_dz-dw_dy), 
							-0.5*(du_dz-dw_dx), 	-0.5*(dv_dz-dw_dy), 	0 };*/
		//double tau[3][3];
		//int m, n;
		
		//lS[k][j][i] = S;
							
		if ( nvert[k][j][i] < 0.1 ) {
			double alpha_star=1.0;
			double r = 1.;
			
			if(levelset) r = rho[k][j][i];
			if(lowRe) {
			  double ren = user->ren;
			  if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
			   
			  Get_alpha_beta_star(ren, Km_o, Om_o, &alpha1, &alpha_star, &beta_star);
			}

			double alpha0 = alpha1, beta0 = beta1;
			double Dissipation = beta_star * Km * Om;
			double Production = nu_t * S * S / r;
			
			/* K RHS */
			// f~ucont*K
			// rhs~f*ajc~U*K*ajc~u*L^2*K/L^3~uK/L
			//komega_rhs[k][j][i].x = ( fp_K1[k][j][i] - fp_K1[k][j][i-1] + fp_K2[k][j][i] - fp_K2[k][j-1][i] + fp_K3[k][j][i] - fp_K3[k-1][j][i] ) * ajc;	// advection, diffusion
			komega_rhs[k][j][i].x = ( fp1[k][j][i].x - fp1[k][j][i-1].x + fp2[k][j][i].x - fp2[k][j-1][i].x + fp3[k][j][i].x - fp3[k-1][j][i].x ) * ajc;	// advection
			komega_rhs[k][j][i].x += ( visc1[k][j][i].x - visc1[k][j][i-1].x + visc2[k][j][i].x - visc2[k][j-1][i].x + visc3[k][j][i].x - visc3[k-1][j][i].x ) * ajc / r;	// diffusion
			
			if(rans==3) {
			  komega_rhs[k][j][i].x += PetscMin ( Production, 10*beta_star*Km_o*Om_o ); //for OSL
			  // komega_rhs[k][j][i].x += PetscMax( PetscMin ( Production, 10*beta_star*Km_o*Om_o ), 0 ); // cross vane indoor levelset; for test
			}
			else komega_rhs[k][j][i].x += Production;
			
			komega_rhs[k][j][i].x -= Dissipation;
			
			if(rans==3) {
				assert( f1>=0 && f1<=1.0 );
				alpha1 = beta1/beta_star - sigma_o1*kappa*kappa/sqrt(beta_star);
				alpha2 = beta2/beta_star - sigma_o2*kappa*kappa/sqrt(beta_star);
				
				alpha0 = alpha1*f1 + alpha2*(1.0-f1);
				beta0 = beta1*f1 + beta2*(1.0-f1);	
			}

			/* Omega RHS */
			//komega_rhs[k][j][i].y = ( fp_Omega1[k][j][i] - fp_Omega1[k][j][i-1] + fp_Omega2[k][j][i] - fp_Omega2[k][j-1][i] + fp_Omega3[k][j][i] - fp_Omega3[k-1][j][i] ) * ajc;	// advection, diffusion
			komega_rhs[k][j][i].y = ( fp1[k][j][i].y - fp1[k][j][i-1].y + fp2[k][j][i].y - fp2[k][j-1][i].y + fp3[k][j][i].y - fp3[k-1][j][i].y ) * ajc;	// advection
			komega_rhs[k][j][i].y += ( visc1[k][j][i].y - visc1[k][j][i-1].y + visc2[k][j][i].y - visc2[k][j-1][i].y + visc3[k][j][i].y - visc3[k-1][j][i].y ) * ajc / r;	// diffusion
			/*
			if(rans==3) komega_rhs[k][j][i].y += alpha0 * S * S;// * inv_nu_t;->not necessary
			else*/
			komega_rhs[k][j][i].y += alpha0 * alpha_star * S * S;// * inv_nu_t;->not necessary
			
			komega_rhs[k][j][i].y -= beta0 * Om * Om;
			
			//if(rans==3 && cross_diffusion) komega_rhs[k][j][i].y += 2.0*(1.0-f1) * sigma_o2 / Om_o *  ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ); // cross-diffusion
			if(rans==3 && cross_diffusion) komega_rhs[k][j][i].y += 2.0*(1.0-f1) * sigma_o2 / Om *  ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ); // cross-diffusion

		}
	}
	
	if(rans==3) DAVecRestoreArray(da, user->lF1, &lf1);
	DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
	DAVecRestoreArray(fda2, Fp1, &fp1);
	DAVecRestoreArray(fda2, Fp2, &fp2);
	DAVecRestoreArray(fda2, Fp3, &fp3);
	DAVecRestoreArray(fda2, Visc1, &visc1);
	DAVecRestoreArray(fda2, Visc2, &visc2);
	DAVecRestoreArray(fda2, Visc3, &visc3);
	
	DAVecRestoreArray(fda, user->lUcont, &ucont);
	DAVecRestoreArray(fda, user->lUcat,  &ucat);
	
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
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->Distance, &distance);
	
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);
	
	DAVecRestoreArray(fda2, user->lK_Omega, &K_Omega);
	DAVecRestoreArray(fda2, user->lK_Omega_o, &K_Omega_o);
	DAVecRestoreArray(fda2, KOmega_RHS, &komega_rhs);
	
	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}
	
	VecDestroy(Fp1);
	VecDestroy(Fp2);
	VecDestroy(Fp3);
	VecDestroy(Visc1);
	VecDestroy(Visc2);
	VecDestroy(Visc3);
};


void K_Omega_IC(UserCtx *user)
{
	DA		da = user->da;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	double nu = 1./user->ren;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	Cmpnts2 ***K_Omega;
	PetscReal	***nvert;

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

	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
	
	// BC for K, Omega
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	// pressure node; cell center
		
		K_Omega[k][j][i].y = C1;
		K_Omega[k][j][i].x = pow ( 10.0, -C2 ) * nu * C1;
		
		if( nvert[k][j][i]>0.1) {
			K_Omega[k][j][i].x = K_Omega[k][j][i].y = 0;
		}
		else if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {
			K_Omega[k][j][i].y = C1;
			K_Omega[k][j][i].x = pow ( 10.0, -C2 ) * nu * C1;
		}
	}
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
	
	DALocalToGlobal(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
	VecCopy(user->K_Omega, user->K_Omega_o);
	
};


void K_Omega_BC(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k, ibi;
	//double nu = 1./user->ren;
	PetscReal	***aj, ***distance, ***rho, ***mu;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	Cmpnts2 ***K_Omega;
	Cmpnts	***csi, ***eta, ***zet, ***ucat;
	PetscReal	***nvert, ***ustar;

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

	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(fda, user->lUcat, &ucat);
	
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lUstar, &ustar);
	DAVecGetArray(da, user->Distance, &distance);
	
	DAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
	
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	
	// BC for K, Omega
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	// pressure node; cell center
		double ren = user->ren;
		//if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
		
		// from saved inflow file
		if(inletprofile==100 && k==1) {
			K_Omega[k-1][j][i] = user->komega_plane[j][i];
			K_Omega[k-1][j][i].x = std::max ( K_Omega[k-1][j][i].x, 0. );
			K_Omega[k-1][j][i].y = std::max ( K_Omega[k-1][j][i].y, 1.e-4);
		}
		else if ( user->bctype[4]==5 && k==1 && nvert[k][j][i]<0.1) {	// inflow
			K_Omega[k-1][j][i].y = 2*C1 - K_Omega[k][j][i].y;
			K_Omega[k-1][j][i].x = 2*pow ( 10.0, -C2 ) / ren * C1 - K_Omega[k][j][i].x;
		}
		
		// slip
		if ( user->bctype[0] == 10 && i==1 ) K_Omega[k][j][i-1] = K_Omega[k][j][i];
		if ( user->bctype[1] == 10 && i==mx-2 ) K_Omega[k][j][i+1] = K_Omega[k][j][i];
		if ( user->bctype[2] == 10 && j==1 ) K_Omega[k][j-1][i] = K_Omega[k][j][i];
		if ( user->bctype[3] == 10 && j==my-2 ) K_Omega[k][j+1][i] = K_Omega[k][j][i];
		if ( user->bctype[4] == 10 && k==1 ) K_Omega[k-1][j][i] = K_Omega[k][j][i];
		if ( user->bctype[5] == 10 && k==mz-2 ) K_Omega[k+1][j][i] = K_Omega[k][j][i];
			
		// outflow
		if ( user->bctype[5] == 4 && k==mz-2 ) {
			K_Omega[k+1][j][i] = K_Omega[k][j][i];
		}
		
		// couette
		if ( user->bctype[3] == 12 && j==my-2 ) {
			double dist = distance[k][j][i];
			K_Omega[k][j][i].y = wall_omega(ren, dist);
			if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x;
		}
		
		
		// wall
		if ( user->bctype[0] == 1 && i<=1 ) {
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dist = distance[k][j][i];
			
			dist = 0.5/aj[k][j][i]/area;
			
			K_Omega[k][j][i].y = wall_omega(ren, dist);
			if(i==1) K_Omega[k][j][i-1].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j][i-1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
		}
		if ( user->bctype[1] == 1 && i>=mx-2 ) {
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
			dist = 0.5/aj[k][j][i]/area;
			
			K_Omega[k][j][i].y = wall_omega(ren, dist);
			if(i==mx-2) K_Omega[k][j][i+1].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j][i+1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
		}
		if ( user->bctype[2] == 1 && j<=1 ) {
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
			dist = 0.5/aj[k][j][i]/area;
			
			K_Omega[k][j][i].y = wall_omega(ren, dist);
			if(j==1) K_Omega[k][j-1][i].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j-1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
		}
		if ( user->bctype[3] == 1 && j>=my-2 ) {
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
			dist = 0.5/aj[k][j][i]/area;
			
			K_Omega[k][j][i].y = wall_omega(ren, dist);
			if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j+1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
		}
				
		// wall function
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[0]==-1 || user->bctype[0]==-2) && i==1) || ( (user->bctype[1]==-1 || user->bctype[1]==-2) &&  i==mx-2) ) && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) {
		  
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double sb, sc; 
			Cmpnts Ua, Uc;
			Ua.x = Ua.y = Ua.z = 0;

			sb = 0.5/aj[k][j][i]/area;
			
			if(i==1) {
				sc = 2*sb + 0.5/aj[k][j][i+1]/area;
				Uc = ucat[k][j][i+1];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j][i-1]/area;
				Uc = ucat[k][j][i-1];
			}
			/*
			double ni[3], nj[3], nk[3];
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;

			wall_function (1./ren, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
			*/
			double utau = ustar[k][j][i];
			
			double Kc = utau*utau/sqrt(0.09);
			//double Oc = utau/sqrt(0.09)/(kappa*sc);
			double Ob = utau/sqrt(0.09)/(kappa*sb);
			
			K_Omega[k][j][i].x = Kc;
			K_Omega[k][j][i].y = Ob;
		}
		
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[2]==-1 || user->bctype[2]==-2) && j==1) || ( (user->bctype[3]==-1 || user->bctype[3]==-2) &&  j==my-2) ) && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1)) {
		  
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double sb, sc; 
			Cmpnts Ua, Uc;
			Ua.x = Ua.y = Ua.z = 0;

			sb = 0.5/aj[k][j][i]/area;
			
			if(j==1) {
				sc = 2*sb + 0.5/aj[k][j+1][i]/area;
				Uc = ucat[k][j+1][i];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j-1][i]/area;
				Uc = ucat[k][j-1][i];
			}
			/*
			double ni[3], nj[3], nk[3];
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;

			wall_function (1./ren, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
			*/
			double utau = ustar[k][j][i];
			
			double Kc = utau*utau/sqrt(0.09);
			//double Oc = utau/sqrt(0.09)/(kappa*sc);
			double Ob = utau/sqrt(0.09)/(kappa*sb);
			
			K_Omega[k][j][i].x = Kc;
			K_Omega[k][j][i].y = Ob;
		}
		
		if ( nvert[k][j][i] > 1.1 ) K_Omega[k][j][i].x = K_Omega[k][j][i].y = 0.;
	}
	DAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);

        DALocalToLocalBegin(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);
        DALocalToLocalEnd(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);

        DAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
	
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
	
		if(flag) K_Omega[k][j][i] = K_Omega[c][b][a];
	}
	
	if(immersed)
	for(ibi=0; ibi<NumberOfBodies; ibi++)
	{
		extern IBMNodes *ibm_ptr;
		IBMNodes *ibm = &ibm_ptr[ibi];
		
		IBMListNode *current;
		current = user->ibmlist[ibi].head;
		while (current) {
			IBMInfo *ibminfo = &current->ibm_intp;
			//int ni = ibminfo->cell;
			current = current->next;
			double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
			int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
			int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
			int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
			i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
			double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
			double Kc = (K_Omega[kp1][jp1][ip1].x * sk1 + K_Omega[kp2][jp2][ip2].x * sk2 + K_Omega[kp3][jp3][ip3].x * sk3);
			
			double ren = user->ren;
			if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
				
			if(wallfunction && ti>0) {
				double utau = ustar[k][j][i];
				const double yplus_min = 0.25;
				sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				double K, Ob;//, Om;
				
				// Wilcox pp.108-109
				K = utau*utau/sqrt(0.09);
				Ob = utau/sqrt(0.09)/(kappa*sb);
				
				K_Omega[k][j][i].x = K;
				K_Omega[k][j][i].y = Ob;
			}
			else {
				const double yplus_min = 0.25;
				
				double utau = ustar[k][j][i];
				
				K_Omega[k][j][i].x = Kc * sb / sc;
				sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				K_Omega[k][j][i].y = wall_omega(ren, sb);	
				
				if ( K_Omega[k][j][i].x < 0 ) K_Omega[k][j][i].x = utau*utau/sqrt(0.09);
			}
			if(user->bctype[4]==5 && k==1) K_Omega[k][j][i]=K_Omega[k-1][j][i];
		};
	}
	
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(fda, user->lUcat, &ucat);
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lUstar, &ustar);
	DAVecRestoreArray(da, user->Distance, &distance);
	
	DAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
	
	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}
	
	DALocalToLocalBegin(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);
	DALocalToLocalEnd(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);
};

PetscErrorCode FormFunction_K_omega(SNES snes, Vec K_Omega, Vec Rhs, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert;
	Cmpnts2 ***rhs;

	DAGetLocalInfo(user->da, &info);
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
	

	DAGlobalToLocalBegin(user->fda2, K_Omega, INSERT_VALUES, user->lK_Omega);
	DAGlobalToLocalEnd(user->fda2, K_Omega, INSERT_VALUES, user->lK_Omega);
	
	K_Omega_BC(user);
	RHS_K_Omega(user, Rhs);
		
	
	VecAXPY(Rhs, -1/user->dt, K_Omega);
	VecAXPY(Rhs, 1/user->dt, user->K_Omega_o);
	
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->fda2, Rhs, &rhs);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1) {
			rhs[k][j][i].x = rhs[k][j][i].y = 0;
		}
		
		// couette
		if ( user->bctype[3] == 12 && j==my-2 ) rhs[k][j][i].y = 0;
	
		//wall_omega
		if( i<=1 && user->bctype[0]==1 ) rhs[k][j][i].y = 0;
		if( i>=mx-2 && user->bctype[1]==1 ) rhs[k][j][i].y = 0;
		
		if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i].y = 0;
		if( j>=my-2 && user->bctype[3]==1 ) rhs[k][j][i].y = 0;
		
		if( k==1 && user->bctype[4]==1 ) rhs[k][j][i].y = 0;
		if( k==mz-2 && user->bctype[5]==1 ) rhs[k][j][i].y = 0;
		
		// wall function k, omega
		if ( 	( i==1 && user->bctype[0] == -1 ) || ( i==mx-2 && user->bctype[1] == -1 ) ||
			( j==1 && user->bctype[2] == -1 ) || ( j==my-2 && user->bctype[3] == -1 ) ||
			( k==1 && user->bctype[4] == -1 ) || ( k==mz-2 && user->bctype[5] == -1 ) ||
			( i==1 && user->bctype[0] == -2 ) || ( i==mx-2 && user->bctype[1] == -2 ) ||
                        ( j==1 && user->bctype[2] == -2 ) || ( j==my-2 && user->bctype[3] == -2 ) ||
                        ( k==1 && user->bctype[4] == -2 ) || ( k==mz-2 && user->bctype[5] == -2 )
			) {
			rhs[k][j][i].x = 0;
			rhs[k][j][i].y = 0;
		}

	}
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->fda2, Rhs, &rhs);
	
	return(0);
}

int snes_ko_created=0;
Vec r_ko;
Mat J_ko;
SNES snes_two_eq;

void Solve_K_Omega(UserCtx *user)
{
	
	KSP ksp;
	PC pc;
	
	
	double norm;
	
	int bi=0;
	double tol=5.e-5;//1.e-6
	
	if(!snes_ko_created) {
		snes_ko_created=1;
		
		VecDuplicate(user[bi].K_Omega, &r_ko);
		SNESCreate(PETSC_COMM_WORLD,&snes_two_eq);
		SNESSetFunction(snes_two_eq,r_ko,FormFunction_K_omega,(void *)&user[bi]);
		MatCreateSNESMF(snes_two_eq, &J_ko);
		SNESSetJacobian(snes_two_eq,J_ko,J_ko,MatMFFDComputeJacobian,(void *)&user[bi]);
		SNESSetType(snes_two_eq, SNESTR);			//SNESTR,SNESLS	
		SNESSetMaxLinearSolveFailures(snes_two_eq,10000);
		SNESSetMaxNonlinearStepFailures(snes_two_eq,10000);		
		SNESKSPSetUseEW(snes_two_eq, PETSC_TRUE);
		SNESKSPSetParametersEW(snes_two_eq,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		SNESSetTolerances(snes_two_eq,PETSC_DEFAULT,tol,PETSC_DEFAULT,50,50000);
			
		SNESGetKSP(snes_two_eq, &ksp);
		KSPSetType(ksp, KSPGMRES);
		//KSPGMRESSetPreAllocateVectors(ksp);
		
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCNONE);
		
		int maxits=20/*10000*/;	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	}
	
	extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy);
	SNESMonitorSet(snes_two_eq,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	
	PetscPrintf(PETSC_COMM_WORLD, "\nSolving K-Omega 0...\n");
	
	SNESSolve(snes_two_eq, PETSC_NULL, user[bi].K_Omega);
	
	SNESGetFunctionNorm(snes_two_eq, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "\nK_Omega SNES residual norm=%.5e\n\n", norm);

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	DAGetLocalInfo(user->da, &info);
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
	
	DAGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	DAGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

	Cmpnts2 ***k_omega, ***lk_omega;
	
	DAVecGetArray(user->fda2, user->K_Omega, &k_omega);
	DAVecGetArray(user->fda2, user->lK_Omega, &lk_omega);

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
		
		/*
                if(k_periodic && k==0) c=mz-2, flag=1;
                else if(k_periodic && k==mz-1) c=1, flag=1;
                if(j_periodic && j==0) b=my-2, flag=1;
                else if(j_periodic && j==my-1) b=1, flag=1;
                if(i_periodic && i==0) a=mx-2, flag=1;
                else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
                else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		*/
		if(flag) k_omega[k][j][i] = lk_omega[c][b][a];
	}
	DAVecRestoreArray(user->fda2, user->K_Omega, &k_omega);
	DAVecRestoreArray(user->fda2, user->lK_Omega, &lk_omega);

	DAGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	DAGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	K_Omega_BC(user);
	DALocalToGlobal(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
};

void K_Omega_Set_Constant(UserCtx *user)
{
	DA		da = user->da;
	DA		fda = user->fda;
	DA		fda2 = user->fda2;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts2 ***K_Omega;
	Cmpnts ***ucat, ***csi, ***eta, ***zet;
	PetscReal	***nvert, ***nu_t, ***f1, ***aj, ***distance, ***rho, ***mu;

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
	
	DAGlobalToLocalBegin(fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	DAGlobalToLocalEnd(fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	
	DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lNu_t, &nu_t);
	if(rans==3) DAVecGetArray(da, user->lF1, &f1);
	DAVecGetArray(da, user->Distance, &distance);
	DAVecGetArray(fda2, user->lK_Omega, &K_Omega);
	DAVecGetArray(fda, user->lUcat, &ucat);
	
	// BC for K, Omega
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	// pressure node; cell center
		double Km = K_Omega[k][j][i].x, Om = K_Omega[k][j][i].y;
		Km = PetscMax ( Km, 0 );
		Om = PetscMax ( Om, 1.e-5 );

		if( nvert[k][j][i]>0.1 || i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {
			nu_t[k][j][i] = 0;
			if(rans==3) f1[k][j][i] = 1;
			continue;
		}
		
		double ren = user->ren, r=1.0;
		if(levelset) {
			ren = rho[k][j][i]/mu[k][j][i];
			r =  rho[k][j][i];
		}

		double tmp, alpha_star=1;
		if(lowRe) {
		  Get_alpha_beta_star(ren, Km, Om, &tmp, &alpha_star, &tmp);
		}

		
		if(rans==3) {
			const double a1=0.31, beta_star=0.09;
			double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
			double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
			double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
			double ajc = aj[k][j][i], d = distance[k][j][i];

			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			
			double dkdc, dodc, dkde, dode, dkdz, dodz;
			double dk_dx, dk_dy, dk_dz, do_dx,  do_dy, do_dz;
			
			Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
					&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			Compute_dkdo_center (i, j, k, mx, my, mz, K_Omega, nvert,  &dkdc,  &dodc,  &dkde,  &dode,  &dkdz,  &dodz );
			Compute_dkdo_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dkdc, dodc, dkde, dode, dkdz, dodz, &dk_dx, &do_dx, &dk_dy, &do_dy, &dk_dz, &do_dz);
			
			double CD = PetscMax ( 2*sigma_o2/Om * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-10 );
			double arg1 = PetscMin ( PetscMax ( sqrt(Km)/beta_star/Om/d, 500./ren/Om/(d*d) ),  4*sigma_o2*Km/CD/(d*d) );
			f1[k][j][i] = tanh ( pow(arg1, 4 ) );
			
			double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
			double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
			double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
			double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
			double vorticity = sqrt( pow( dw_dy - dv_dz, 2. ) + pow ( du_dz - dw_dx, 2. ) +  pow( dv_dx - du_dy, 2. ) );
			double arg2 = PetscMax ( 2*sqrt(Km)/beta_star/Om/d, 500/ren/Om/(d*d) );
			double f2 = tanh ( arg2*arg2 );
			
			nu_t[k][j][i] = a1 * Km * alpha_star / PetscMax ( a1* Om, f2 * S ) * r;//1994
			//nu_t[k][j][i] = a1 * Km * alpha_star / PetscMax ( a1* Om, f2 * vorticity ) * r; //original
		}
		else {
			nu_t[k][j][i] = Km/Om * alpha_star * r;
		}
		
		nu_t[k][j][i] = PetscMax(nu_t[k][j][i], 0);
	}
	
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
			nu_t[k][j][i] = nu_t[c][b][a];
			if(rans==3) f1[k][j][i] = f1[c][b][a];
		}
	}
	
	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lNu_t, &nu_t);
	if(rans==3) DAVecRestoreArray(da, user->lF1, &f1);
	DAVecRestoreArray(da, user->Distance, &distance);
	DAVecRestoreArray(fda2, user->lK_Omega, &K_Omega);
	DAVecRestoreArray(fda, user->lUcat, &ucat);

	DALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	DALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	if(rans==3)  {
		DALocalToLocalBegin(da, user->lF1, INSERT_VALUES, user->lF1);
		DALocalToLocalEnd(da, user->lF1, INSERT_VALUES, user->lF1);
	}
};

