/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

static char help[] = "Testing programming!";

#include <vector>
#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include <stdio.h>
#include <stdlib.h>

#define NEWMETRIC

#ifdef TECIO
#include "TECIO.h"
#endif
PetscInt ti, block_number, Flux_in;
int binary_input=0;
int xyz_input=0;
PetscInt tis, tie, tsteps=5;
PetscReal angle;
int nv_once=0;
int onlyV=0;
int k_average=0;
int j_average=0;
int i_average=0;
int ik_average=0;
int ikc_average=0;	// conditional spatial averaging in ik directions (channel flow)
int reynolds=0;	// 1: contravariant reynolds stress

int i_begin, i_end;
int j_begin, j_end;
int k_begin, k_end;

int pcr=0;
int avg=0, rans=0, rans_output=0, levelset=0;
int vc = 1;

int cs=0;
int i_periodic=0;
int j_periodic=0;
int k_periodic=0;
int kk_periodic=0;
int averaging_option=0;
int pi=-1, pk=-1;
int shear=0;

char prefix[256];

//int l, m, n;
/* Symmetric matrix A -> eigenvectors in columns of V, corresponding eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);


typedef struct {
  PetscReal t, f;
} FlowWave;

typedef struct {
  PassiveScalar u, v, w, p;
} PassiveField;

typedef struct {
  PetscScalar u, v, w;
} Field;

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  PetscScalar x, y;
} Cmpnts2;

typedef struct {
  PassiveScalar csi[3], eta[3], zet[3], aj;
} Metrics;

typedef struct {
  Vec	Ubcs; // An ugly hack, waste of memory
} BCS;


typedef struct {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;	// number of vertices and number of elements
  PetscInt	*nv1, *nv2, *nv3;	// Node index of each triangle
  PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
/*   PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270]; */
  Cmpnts	*u, *uold;
  
} IBMNodes;

typedef struct {
  PetscInt	IM, JM, KM; // dimensions of grid
  DA da;	/* Data structure for scalars (include the grid geometry
		   informaion, to obtain the grid information, 
		   use DAGetCoordinates) */
  DA fda, fda2;	// Data Structure for vectors
  DALocalInfo info;

  Vec	Cent;	// Coordinates of cell centers
  Vec 	Csi, Eta, Zet, Aj;
  Vec 	ICsi, IEta, IZet, IAj;
  Vec 	JCsi, JEta, JZet, JAj;
  Vec 	KCsi, KEta, KZet, KAj;
  Vec 	Ucont;	// Contravariant velocity components
  Vec 	Ucat;	// Cartesian velocity components
  Vec	Ucat_o;
  Vec 	P;
  Vec	Phi;
  Vec	GridSpace;
  Vec	Nvert;
  Vec	Nvert_o;
  BCS	Bcs;

  PetscInt	*nvert;//ody property

  PetscReal	ren;	// Reynolds number
  PetscReal	dt; 	// time step
  PetscReal	st;	// Strouhal number

  PetscReal	r[101], tin[101], uinr[101][1001];

  Vec	lUcont, lUcat, lP, lPhi;
  Vec	lCsi, lEta, lZet, lAj;
  Vec	lICsi, lIEta, lIZet, lIAj;
  Vec	lJCsi, lJEta, lJZet, lJAj;
  Vec	lKCsi, lKEta, lKZet, lKAj;
  Vec	lGridSpace;
  Vec	lNvert, lNvert_o;
  Vec	lCent;
  
  Vec Ucat_sum;		// u, v, w
  Vec Ucat_cross_sum;		// uv, vw, wu
  Vec Ucat_square_sum;	// u^2, v^2, w^2

  PetscInt _this;

  FlowWave *inflow, *kinematics;
  PetscInt number_flowwave, number_kinematics;

} UserCtx;

PetscErrorCode ReadCoordinates(UserCtx *user);
PetscErrorCode QCriteria(UserCtx *user);
PetscErrorCode Velocity_Magnitude(UserCtx *user);
PetscErrorCode Lambda2(UserCtx *user);
PetscErrorCode FormMetrics(UserCtx *user);
void Calc_avg_shear_stress(UserCtx *user);


int file_exist(char *str)
{
  int r=0;

  /*if(!my_rank)*/ {
    FILE *fp=fopen(str, "r");
    if(!fp) {
      r=0;
      printf("\nFILE !!! %s does not exist !!!\n", str);
    }
    else {
      fclose(fp);
      r=1;
    }
  }
  MPI_Bcast(&r, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  return r;
};

void Calculate_Covariant_metrics(double g[3][3], double G[3][3])
{
	/*
		| csi.x  csi.y csi.z |-1		| x.csi  x.eta x.zet | 
		| eta.x eta.y eta.z |	 =	| y.csi   y.eta  y.zet |
		| zet.x zet.y zet.z |		| z.csi  z.eta z.zet |
	
	*/
	const double a11=g[0][0], a12=g[0][1], a13=g[0][2];
	const double a21=g[1][0], a22=g[1][1], a23=g[1][2];
	const double a31=g[2][0], a32=g[2][1], a33=g[2][2];

	double det= a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);
	
	G[0][0] = (a33*a22-a32*a23)/det,	G[0][1] = - (a33*a12-a32*a13)/det, 	G[0][2] = (a23*a12-a22*a13)/det;
	G[1][0] = -(a33*a21-a31*a23)/det, G[1][1] = (a33*a11-a31*a13)/det,	G[1][2] = - (a23*a11-a21*a13)/det;
	G[2][0] = (a32*a21-a31*a22)/det,	G[2][1] = - (a32*a11-a31*a12)/det,	G[2][2] = (a22*a11-a21*a12)/det;
};

void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3])
{
	double g[3][3];
	double G[3][3];
	
	g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
	g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
	g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;
	
	Calculate_Covariant_metrics(g, G);
	double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
	double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
	double xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];
	      
	double nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
	double nx_j = xeta, ny_j = yeta, nz_j = zeta;
	double nx_k = xzet, ny_k = yzet, nz_k = zzet;
	      
	double sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
	double sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
	double sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

	nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
	nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
	nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

	ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
	nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
	nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;
}

double Contravariant_Reynolds_stress(double uu, double uv, double uw, double vv, double vw, double ww,
	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2)
{
	
	double A = uu*csi0*eta0 + vv*csi1*eta1 + ww*csi2*eta2 + uv * (csi0*eta1+csi1*eta0)	+ uw * (csi0*eta2+csi2*eta0) + vw * (csi1*eta2+csi2*eta1);
	double B = sqrt(csi0*csi0+csi1*csi1+csi2*csi2)*sqrt(eta0*eta0+eta1*eta1+eta2*eta2);
	
	return A/B;
}

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
{
	if ((nvert[k][j][i+1])> 0.1) {
		*dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
		*dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
		*dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
	}
	else if ((nvert[k][j][i-1])> 0.1) {
		*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
		*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
		*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
	}
	else {
		if(i_periodic && i==1) {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][mx-2].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][mx-2].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][mx-2].z ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dudc = ( ucat[k][j][1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][1].z - ucat[k][j][i-1].z ) * 0.5;
		}
		else {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
		}
	}

	if ((nvert[k][j+1][i])> 0.1) {
		*dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
		*dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
		*dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
	}
	else if ((nvert[k][j-1][i])> 0.1) {
		*dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
		*dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
		*dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );
	}
	else {
		if(j_periodic && j==1) {
			*dude = ( ucat[k][j+1][i].x - ucat[k][my-2][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][my-2][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][my-2][i].z ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dude = ( ucat[k][1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
		else {
			*dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
	}

	if ((nvert[k+1][j][i])> 0.1) {
		*dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
		*dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
		*dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
	}
	else if ((nvert[k-1][j][i])> 0.1) {
		*dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
		*dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
		*dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
	}
	else {
		if(k_periodic && k==1) {
			*dudz = ( ucat[k+1][j][i].x - ucat[mz-2][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[mz-2][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[mz-2][j][i].z ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dudz = ( ucat[1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
		else {
			*dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
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


void IKavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;

	for (j=ys; j<ye-2; j++) {
		double iksum=0;
		int count=0;
		for (i=xs; i<xe-2; i++)
		for (k=zs; k<ze-2; k++) {
			iksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
			count++;
		}
		double ikavg = iksum/(double)count;
		for (i=xs; i<xe-2; i++)
		for (k=zs; k<ze-2; k++) 	x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ikavg;
	}
}

/*
	pi, pk : # of grid points correcsponding to the period
	conditional averaging
*/
void IKavg_c (float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  //int i, j, k;
	
  if(pi<=0) pi = (xe-xs-2); // no averaging
  if(pk<=0) pk = (ze-zs-2); // no averaging 
	
	int ni = (xe-xs-2) / pi;
	int nk = (ze-zs-2) / pk;

	std::vector< std::vector<float> > iksum (pk);
	
	for(int k=0; k<pk; k++) iksum[k].resize(pi);
	
	for (int j=ys; j<ye-2; j++) {
		
		for(int k=0; k<pk; k++) std::fill( iksum[k].begin(), iksum[k].end(), 0.0 );
		
		for (int i=xs; i<xe-2; i++)
		for (int k=zs; k<ze-2; k++) {
			iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi] += x [k * (mx-2)*(my-2) + j*(mx-2) + i];
		}

		for (int i=xs; i<xe-2; i++)
		for (int k=zs; k<ze-2; k++) x [k * (mx-2)*(my-2) + j*(mx-2) + i] = iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi ] / (ni*nk);
	}
}


void Kavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;

	for (j=ys; j<ye-2; j++) 
	for (i=xs; i<xe-2; i++) {
	  double ksum=0;
	  int count=0;
	  for (k=zs; k<ze-2; k++) {
		ksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
		count++;
	  }
	  double kavg = ksum/(double)count;
	  for (k=zs; k<ze-2; k++) {
		x[k * (mx-2)*(my-2) + j*(mx-2) + i] = kavg;
	  }
	}
      
}

void Javg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++)
	for (i=xs; i<xe-2; i++) {
	  double jsum=0;
	  int count=0;
	  for (j=ys; j<ye-2; j++) {
		jsum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
		count++;
	  }
	  double javg = jsum/(double)count;
	  for (j=ys; j<ye-2; j++) {
		x[k * (mx-2)*(my-2) + j*(mx-2) + i] = javg;
	  }
	}
      
}

void Iavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) {
	  for (j=ys; j<ye-2; j++) {
	    double isum=0;
	    int count=0;
	    for (i=xs; i<xe-2; i++) {
	      isum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
	      count++;
	    }
	    double iavg = isum/(double)count;
	    for (i=xs; i<xe-2; i++) {
	      x[k * (mx-2)*(my-2) + j*(mx-2) + i] = iavg;
	    }
	  }
	}
}

void Iavg(Cmpnts ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) 
	for (j=ys; j<ye-2; j++) {
		Cmpnts isum, iavg;
		isum.x = isum.y = isum.z = 0;
	    
	  int count=0;
		for (i=xs; i<xe-2; i++) {
	     isum.x += x[k+1][j+1][i+1].x;
	     isum.y += x[k+1][j+1][i+1].y;
	     isum.z += x[k+1][j+1][i+1].z;
	     count++;
	  }
	  
	  iavg.x = isum.x /(double)count;
	  iavg.y = isum.y /(double)count;
	  iavg.z = isum.z /(double)count;
	    
		for (i=xs; i<xe-2; i++) {
	  	x[k+1][j+1][i+1] = iavg;
		}
	}
}

void Iavg(PetscReal ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) 
	for (j=ys; j<ye-2; j++) {
		double isum, iavg;
		isum = 0;
	    
	  int count=0;
		for (i=xs; i<xe-2; i++) {
	     isum += x[k+1][j+1][i+1];
	     count++;
	  }
	  iavg = isum /(double)count;
	    
		for (i=xs; i<xe-2; i++) {
	  	x[k+1][j+1][i+1] = iavg;
		}
	}
}

PetscErrorCode TECIOOut_V(UserCtx *user, int only_V)	// seokkoo
{
	PetscInt bi;

	char filen[80];
	sprintf(filen, "%sResult%05d.plt", prefix, ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;
	
	
	
	if(only_V)   {
		if(cs) I = TECINI100((char*)"Flow", (char*)"X Y Z UU Cs", filen, (char*)".", &Debug, &VIsDouble);
		else if(only_V==2) I = TECINI100((char*)"Flow", (char*)"X Y Z UU", filen, (char*)".", &Debug, &VIsDouble);
		else I = TECINI100((char*)"Flow", (char*)"X Y Z UU Nv", filen, (char*)".", &Debug, &VIsDouble);
	}
	else if(rans/* && rans_output*/) {
		if(levelset) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv K Omega Nut Level", filen, (char*)".", &Debug, &VIsDouble);
		else I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv K Omega Nut ", filen, (char*)".", &Debug, &VIsDouble);
	}
	else  {
		if(cs) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Cs", filen, (char*)".", &Debug, &VIsDouble);
		else if(levelset) {
			I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv Level", filen, (char*)".", &Debug, &VIsDouble);
		}
		else I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv", filen, (char*)".", &Debug, &VIsDouble);
		
	}

	for (bi=0; bi<block_number; bi++) {
		DA da = user[bi].da, fda = user[bi].fda;
		DALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o, ***csi, ***eta, ***zet;
		PetscReal	***p, ***nvert, ***level;
		Vec		Coor, zCoor, nCoor;
		VecScatter	ctx;
		Vec K_Omega;
		
		DAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		/*DAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		DAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DAVecGetArray(user[bi].fda, user[bi].Zet, &zet);*/
		DAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[40] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */
		
		
		/*************************/
		printf("mi=%d, mj=%d, mk=%d\n", mx, my, mz);
		printf("xs=%d, xe=%d\n", xs, xe);
		printf("ys=%d, ye=%d\n", ys, ye);
		printf("zs=%d, ze=%d\n", zs, ze);
		//exit(0);
		
		i_begin = 1, i_end = mx-1;	// cross section in tecplot
		j_begin = 1, j_end = my-1;
		k_begin = 1, k_end = mz-1;
		
		PetscOptionsGetInt(PETSC_NULL, "-i_begin", &i_begin, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL, "-i_end", &i_end, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL, "-j_begin", &j_begin, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL, "-j_end", &j_end, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL, "-k_begin", &k_begin, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL, "-k_end", &k_end, PETSC_NULL);
		
		xs = i_begin - 1, xe = i_end+1;
		ys = j_begin - 1, ye = j_end+1;
		zs = k_begin - 1, ze = k_end+1;
		
		printf("xs=%d, xe=%d\n", xs, xe);
		printf("ys=%d, ye=%d\n", ys, ye);
		printf("zs=%d, ze=%d\n", zs, ze);
		//exit(0);
		//xs=0, xe=nsection+1;
		/*************************/
		
		if (vc) {
			LOC[3]=0; LOC[4]=0; LOC[5]=0; LOC[6]=0;
		}
		else if(only_V) {
			LOC[4]=0; LOC[5]=0; LOC[6]=0;
		}
		/*
		IMax = mx-1;
		JMax = my-1;
		KMax = mz-1;
		*/
		IMax = i_end - i_begin + 1;
		JMax = j_end - j_begin + 1;
		KMax = k_end - k_begin + 1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		//III = (mx-1) * (my-1) * (mz-1);
		III = IMax*JMax*KMax;
		
		DAGetCoordinates(da, &Coor);
		DAVecGetArray(fda, Coor, &coor);

		float *x;
		x = new float [III];
		
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].x;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].y;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].z;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);
		
		DAVecRestoreArray(fda, Coor, &coor);
		delete []x;
    		
		if(!vc) {
			x = new float [(mx-1)*(my-1)*(mz-1)];
			DAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
			for (k=0; k<mz-1; k++)
			for (j=0; j<my-1; j++)
			for (i=0; i<mx-1; i++) {
				ucat_o[k][j][i].x = 0.125 *
					(ucat[k][j][i].x + ucat[k][j][i+1].x +
					ucat[k][j+1][i].x + ucat[k][j+1][i+1].x +
					ucat[k+1][j][i].x + ucat[k+1][j][i+1].x +
					ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);
				ucat_o[k][j][i].y = 0.125 *
					(ucat[k][j][i].y + ucat[k][j][i+1].y +
					ucat[k][j+1][i].y + ucat[k][j+1][i+1].y +
					ucat[k+1][j][i].y + ucat[k+1][j][i+1].y +
					ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y);
				ucat_o[k][j][i].z = 0.125 *
					(ucat[k][j][i].z + ucat[k][j][i+1].z +
					ucat[k][j+1][i].z + ucat[k][j+1][i+1].z +
					ucat[k+1][j][i].z + ucat[k+1][j][i+1].z +
					ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z);
			}
	      
			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x;
			}
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);

			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y;
			}
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);

			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z;
			}
		
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);
	      
			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = sqrt( ucat_o[k][j][i].x*ucat_o[k][j][i].x + ucat_o[k][j][i].y*ucat_o[k][j][i].y + ucat_o[k][j][i].z*ucat_o[k][j][i].z );
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);

			DAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
			delete []x;
		}
		else {
		  //x = new float [(mx-2)*(my-2)*(mz-2)];
			//III = (mx-2) * (my-2) * (mz-2);
			III = (IMax-1)*(JMax-1)*(KMax-1);
			x = new float [III];

			if(!only_V)  {
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].x;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].y;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].z;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);
			}
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = sqrt(ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x+ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y+ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z);
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
			delete []x;
		}

		
		III = (IMax-1)*(JMax-1)*(KMax-1);
		//III = (mx-2) * (my-2) * (mz-2);
   		x = new float [III];
		//x.resize (III);
    
		if(!only_V) {
			DAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
			DAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		
		
		if(only_V!=2) {
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = nvert[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
		}
		
		if(only_V==2) {	// Z Vorticity
			/*
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
				double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
				double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
				double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
				double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
				double ajc = aj[k][j][i];
			
				if(i==0 || j==0 || k==0) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0;
				else {
					Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
					Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
										&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = dv_dx - du_dy;
				}
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
			*/
		}
		
		if(!onlyV && rans /*&& rans_output*/) {
			Cmpnts2 ***komega;
			DACreateGlobalVector(user->fda2, &K_Omega);
			PetscViewer	viewer;
			sprintf(filen, "kfield%05d_%1d.dat", ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, K_Omega);
			PetscViewerDestroy(viewer);
			DAVecGetArray(user[bi].fda2, K_Omega, &komega);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x;
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].y;
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x/(komega[k+1][j+1][i+1].y+1.e-20);
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			DAVecRestoreArray(user[bi].fda2, K_Omega, &komega);
			VecDestroy(K_Omega);
		}
		if(levelset) {
			PetscReal ***level;
			Vec Levelset;
			DACreateGlobalVector(user->da, &Levelset);
			PetscViewer	viewer;
			sprintf(filen, "lfield%05d_%1d.dat", ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, Levelset);
			PetscViewerDestroy(viewer);
			DAVecGetArray(user[bi].da, Levelset, &level);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = level[k+1][j+1][i+1];
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			DAVecRestoreArray(user[bi].da, Levelset, &level);
			VecDestroy(Levelset);
		}
		
		delete []x;
		/*
		DAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		DAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);*/
		DAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
		DAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
	}
	I = TECEND100();
	return 0;
}


PetscErrorCode TECIOOut_Averaging(UserCtx *user)	// seokkoo
{
	PetscInt bi;

	char filen[80];
	sprintf(filen, "%sResult%05d-avg.plt", prefix, ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	if(pcr) I = TECINI100((char*)"Averaging", "X Y Z P Velocity_Magnitude Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
	else if(avg==1) {

		if(averaging_option==1) {
			//I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UV_ VW_ UW_ Nv",  filen, (char*)".",  &Debug,  &VIsDouble); //OSL
			//I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UU_ VV_ WW_ Nv",  filen, (char*)".",  &Debug,  &VIsDouble); //OSL
			I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
		}
		else if(averaging_option==2) {
			//I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw  UV_ VW_ UW_ P pp Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
			I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw  P pp Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
		}
		else if(averaging_option==3) {
			//I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UV_ VW_ UW_ P pp Vortx Vorty Vortz vortx2 vorty2 vortz2 Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
			I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw P pp Vortx Vorty Vortz vortx2 vorty2 vortz2 Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
		}
		
	}
	else if(avg==2) I = TECINI100((char*)"Averaging", "X Y Z U V W K Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
	

	for (bi=0; bi<block_number; bi++) {
		DA da = user[bi].da, fda = user[bi].fda;
		DALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o;
		Cmpnts	***u2sum, ***u1sum,  ***usum;
		PetscReal	***p, ***nvert;
		Vec		Coor, zCoor, nCoor;
		//VecScatter ctx;

		Vec X, Y, Z, U, V, W;

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DAGetCoordinates(da, &Coor);
		DAVecGetArray(fda, Coor, &coor);

		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
		DAVecRestoreArray(fda, Coor, &coor);
    
		//delete []x;
		double N=(double)tis+1.0;
		//x = new float [(mx-2)*(my-2)*(mz-2)];

		III = (mx-2) * (my-2) * (mz-2);
	
		if(pcr)  {
			DAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			I = TECDAT100(&III, x, &DIsDouble);
			DAVecRestoreArray(user[bi].da, user[bi].P, &p);

			// Load ucat
			PetscViewer	viewer;
			sprintf(filen, "ufield%05d_%1d.dat", ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user[bi].Ucat));
			PetscViewerDestroy(viewer);
			
			DAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] 
								= sqrt ( ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x + ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y + ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z );
			I = TECDAT100(&III, x, &DIsDouble);
			DAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
		}
		else if(avg==1) {
			PetscViewer viewer;
			char filen[128];
			
			DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
			DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
			DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);

			sprintf(filen, "su0_%05d_%1d.dat", ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(viewer);

			sprintf(filen, "su1_%05d_%1d.dat", ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(viewer);
			
			sprintf(filen, "su2_%05d_%1d.dat", ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(viewer);
			
			DAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			DAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
			DAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
			
			// U
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].x/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// V
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].y/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// W
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].z/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// uu, u rms
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for(i=xs; i<xe-2; i++) {
				double U = usum[k+1][j+1][i+1].x/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// vv, v rms
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double V = usum[k+1][j+1][i+1].y/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].y/N - V*V );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// ww, w rms
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double W = usum[k+1][j+1][i+1].z/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].z/N - W*W );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// uv
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].x/N - UV;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// vw
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].y/N - VW;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);     
	      
			// wu
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].z/N - WU;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);           
			
			/*******************************
			//
			DACreateGlobalVector(user[bi].fda, &user[bi].Csi);
			DACreateGlobalVector(user[bi].fda, &user[bi].Eta);
			DACreateGlobalVector(user[bi].fda, &user[bi].Zet);
			DACreateGlobalVector(user[bi].da, &user[bi].Aj);
			FormMetrics(&(user[bi]));
			
			Cmpnts ***csi, ***eta, ***zet;
			DAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
			DAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
			DAVecGetArray(user[bi].fda, user[bi].Zet, &zet);

			//UV_ or UU_
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double U = usum[k+1][j+1][i+1].x/N;
				double V = usum[k+1][j+1][i+1].y/N;
				double W = usum[k+1][j+1][i+1].z/N;
				
				double uu = ( u2sum[k+1][j+1][i+1].x/N - U*U );
				double vv = ( u2sum[k+1][j+1][i+1].y/N - V*V );
				double ww = ( u2sum[k+1][j+1][i+1].z/N - W*W );
				double uv = u1sum[k+1][j+1][i+1].x/N - U*V;
				double vw = u1sum[k+1][j+1][i+1].y/N - V*W;
				double uw = u1sum[k+1][j+1][i+1].z/N - W*U;
				
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0 = eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				
				// UV_
				//x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	csi0, csi1, csi2, eta0, eta1, eta2);
				// UU_
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	csi0, csi1, csi2, csi0, csi1, csi2);
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
			
			// VW_ or VV_
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double U = usum[k+1][j+1][i+1].x/N;
				double V = usum[k+1][j+1][i+1].y/N;
				double W = usum[k+1][j+1][i+1].z/N;
				
				double uu = ( u2sum[k+1][j+1][i+1].x/N - U*U );
				double vv = ( u2sum[k+1][j+1][i+1].y/N - V*V );
				double ww = ( u2sum[k+1][j+1][i+1].z/N - W*W );
				double uv = u1sum[k+1][j+1][i+1].x/N - U*V;
				double vw = u1sum[k+1][j+1][i+1].y/N - V*W;
				double uw = u1sum[k+1][j+1][i+1].z/N - W*U;
				
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0 = eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				
				// VW_
				//x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	zet0, zet1, zet2, eta0, eta1, eta2);
				// VV_
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	eta0, eta1, eta2, eta0, eta1, eta2);
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
			
			//UW_ or WW_
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double U = usum[k+1][j+1][i+1].x/N;
				double V = usum[k+1][j+1][i+1].y/N;
				double W = usum[k+1][j+1][i+1].z/N;
				
				double uu = ( u2sum[k+1][j+1][i+1].x/N - U*U );
				double vv = ( u2sum[k+1][j+1][i+1].y/N - V*V );
				double ww = ( u2sum[k+1][j+1][i+1].z/N - W*W );
				double uv = u1sum[k+1][j+1][i+1].x/N - U*V;
				double vw = u1sum[k+1][j+1][i+1].y/N - V*W;
				double uw = u1sum[k+1][j+1][i+1].z/N - W*U;
				
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0 = eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				
				// UW_
				//x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	csi0, csi1, csi2, zet0, zet1, zet2);
				// WW_
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Contravariant_Reynolds_stress(uu, uv, uw, vv, vw, ww,	zet0, zet1, zet2, zet0, zet1, zet2);
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
			
			DAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
			DAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
			DAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
			
			VecDestroy(user[bi].Csi);
			VecDestroy(user[bi].Eta);
			VecDestroy(user[bi].Zet);
			VecDestroy(user[bi].Aj);
			
			DAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			DAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
			DAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
			
			VecDestroy(user[bi].Ucat_sum);
			VecDestroy(user[bi].Ucat_cross_sum);
			VecDestroy(user[bi].Ucat_square_sum);
			//
			********************************/
			
			if(averaging_option>=2) {
				Vec P_sum, P_square_sum;
				PetscReal ***psum, ***p2sum;

				DACreateGlobalVector(user[bi].da, &P_sum);
				DACreateGlobalVector(user[bi].da, &P_square_sum);

				sprintf(filen, "sp_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, P_sum);
				PetscViewerDestroy(viewer);

				sprintf(filen, "sp2_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, P_square_sum);
				PetscViewerDestroy(viewer);

				DAVecGetArray(user[bi].da, P_sum, &psum);
				DAVecGetArray(user[bi].da, P_square_sum, &p2sum);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double P = psum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = P;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double P = psum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( p2sum[k+1][j+1][i+1]/N - P*P );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				DAVecRestoreArray(user[bi].da, P_sum, &psum);
				DAVecRestoreArray(user[bi].da, P_square_sum, &p2sum);
				
				VecDestroy(P_sum);
				VecDestroy(P_square_sum);
			}

			if(averaging_option>=3) {
				Vec Vort_sum, Vort_square_sum;
				Cmpnts ***vortsum, ***vort2sum;

				DACreateGlobalVector(user[bi].fda, &Vort_sum);
				DACreateGlobalVector(user[bi].fda, &Vort_square_sum);

				sprintf(filen, "svo_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, Vort_sum);
				PetscViewerDestroy(viewer);

				sprintf(filen, "svo2_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, Vort_square_sum);
				PetscViewerDestroy(viewer);

				DAVecGetArray(user[bi].fda, Vort_sum, &vortsum);
				DAVecGetArray(user[bi].fda, Vort_square_sum, &vort2sum);
                        
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortx = vortsum[k+1][j+1][i+1].x/N;
						x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortx;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vorty = vortsum[k+1][j+1][i+1].y/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vorty;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortz = vortsum[k+1][j+1][i+1].z/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortz;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortx = vortsum[k+1][j+1][i+1].x/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].x/N - vortx*vortx );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vorty = vortsum[k+1][j+1][i+1].y/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].y/N - vorty*vorty );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortz = vortsum[k+1][j+1][i+1].z/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].z/N - vortz*vortz );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);                        

				DAVecRestoreArray(user[bi].fda, Vort_sum, &vortsum);
				DAVecRestoreArray(user[bi].fda, Vort_square_sum, &vort2sum);

				VecDestroy(Vort_sum);
				VecDestroy(Vort_square_sum);
			  
			  //haha
			  /*
			  //TKE budget
			 	Vec Udp_sum, dU2_sum, UUU_sum;
				PetscReal ***udpsum, ***aj;
				Cmpnts ***du2sum, ***uuusum;
				Cmpnts ***csi, ***eta, ***zet;
			
				DACreateGlobalVector(user[bi].da, &Udp_sum);
				DACreateGlobalVector(user[bi].fda, &dU2_sum);
				DACreateGlobalVector(user[bi].fda, &UUU_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Csi);
				DACreateGlobalVector(user[bi].fda, &user[bi].Eta);
				DACreateGlobalVector(user[bi].fda, &user[bi].Zet);
				DACreateGlobalVector(user[bi].da, &user[bi].Aj);
				
	
				sprintf(filen, "su0_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, (user[bi].Ucat_sum));
				PetscViewerDestroy(viewer);
	
				sprintf(filen, "su1_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, (user[bi].Ucat_cross_sum));
				PetscViewerDestroy(viewer);
				
				sprintf(filen, "su2_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, (user[bi].Ucat_square_sum));
				PetscViewerDestroy(viewer);
				
				sprintf(filen, "su3_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, Udp_sum);
				PetscViewerDestroy(viewer);
				
				sprintf(filen, "su4_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, dU2_sum);
				PetscViewerDestroy(viewer);
				
				sprintf(filen, "su5_%05d_%1d.dat", ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, UUU_sum);
				PetscViewerDestroy(viewer);
				
				FormMetrics(&(user[bi]));
				
				DAVecGetArray(user[bi].da, user[bi].Aj, &aj);
				DAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
				DAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
				DAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
				DAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
				DAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
				DAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
				DAVecGetArray(user[bi].da, Udp_sum, &udpsum);
				DAVecGetArray(user[bi].fda, dU2_sum, &du2sum);
				DAVecGetArray(user[bi].fda, UUU_sum, &uuusum);
				
				DAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
				DAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
				DAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
				DAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
				DAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
				DAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
				DAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
				DAVecRestoreArray(user[bi].da, Udp_sum, &udpsum);
				DAVecRestoreArray(user[bi].fda, dU2_sum, &du2sum);
				DAVecRestoreArray(user[bi].fda, UUU_sum, &uuusum);
				
				VecDestroy(user[bi].Csi);
				VecDestroy(user[bi].Eta);
				VecDestroy(user[bi].Zet);
				VecDestroy(user[bi].Aj);
				VecDestroy(user[bi].Ucat_sum);
				VecDestroy(user[bi].Ucat_cross_sum);
				VecDestroy(user[bi].Ucat_square_sum);
				VecDestroy(Udp_sum);
				VecDestroy(dU2_sum);
				VecDestroy(UUU_sum);
				*/
			}
		}
		else if(avg==2) {
			PetscViewer viewer;
			Vec K_sum;
			PetscReal ***ksum;
			char filen[128];
			
			DACreateGlobalVector(user->fda, &user->Ucat_sum);
                        DACreateGlobalVector(user->fda, &user->Ucat_square_sum);
			if(rans) {
				DACreateGlobalVector(user->da, &K_sum);
			}

			sprintf(filen, "su0_%05d_%1d.dat", ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user->Ucat_sum));
			PetscViewerDestroy(viewer);
			
			if(rans) {
				sprintf(filen, "sk_%05d_%1d.dat", ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, K_sum);
				PetscViewerDestroy(viewer);
			}
			else {
				sprintf(filen, "su2_%05d_%1d.dat", ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, (user->Ucat_square_sum));
				PetscViewerDestroy(viewer);
			}
			
			DAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			
			if(rans) DAVecGetArray(user[bi].da, K_sum, &ksum);
			else DAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
	      
			// U
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].x/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// V
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].y/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// W
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].z/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// k
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for(i=xs; i<xe-2; i++) {
				if(rans)  {
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ksum[k+1][j+1][i+1]/N;
				}
				else {
					double U = usum[k+1][j+1][i+1].x/N;
					double V = usum[k+1][j+1][i+1].y/N;
					double W = usum[k+1][j+1][i+1].z/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U ) + ( u2sum[k+1][j+1][i+1].y/N - V*V ) + ( u2sum[k+1][j+1][i+1].z/N - W*W );
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] *= 0.5;
				}
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			DAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			
			if(rans) DAVecRestoreArray(user[bi].da, K_sum, &ksum);
			else DAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
			
			VecDestroy(user->Ucat_sum);
			if(rans) VecDestroy(K_sum);
			else VecDestroy(user->Ucat_square_sum);
		}
	
		DAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
		I = TECDAT100(&III, x, &DIsDouble);
		DAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		
		delete []x;
	}
	I = TECEND100();
	return 0;
}

PetscErrorCode TECIOOutQ(UserCtx *user, int Q)
{
	PetscInt bi;

	char filen[80];
	sprintf(filen, "QCriteria%05d.plt", ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	if(Q==1) {
		printf("qcr=%d, Q-Criterion\n", Q);
		//I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
		I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}
	else if(Q==2) {
		printf("Lambda2-Criterion\n");
		//I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
		I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}
	else if(Q==3) {
		printf("Q-Criterion from saved file\n");
		I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}

	for (bi=0; bi<block_number; bi++) {
		DA da = user[bi].da, fda = user[bi].fda;
		DALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o;
		PetscReal	***p, ***nvert;
		Vec		Coor, zCoor, nCoor;
		VecScatter	ctx;

		Vec X, Y, Z, U, V, W;

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[8] = {1, 1, 1, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		III = (mx-1) * (my-1) * (mz-1);
		float	*x;
		PetscMalloc(mx*my*mz*sizeof(float), &x);	// seokkoo
	
		DAGetCoordinates(da, &Coor);
		DAVecGetArray(fda, Coor, &coor);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		}
		I = TECDAT100(&III, x, &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		}
		I = TECDAT100(&III, x, &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		}
      
		I = TECDAT100(&III, x, &DIsDouble);
		DAVecRestoreArray(fda, Coor, &coor);
    
		III = (mx-2) * (my-2) * (mz-2);

		if(Q==1) {
			QCriteria(user);
			DAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] =p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		else if(Q==2) {
			Lambda2(user);
			DAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		else if(Q==3) {
			char filen2[128];
			PetscViewer	viewer;
			
			sprintf(filen2, "qfield%05d_%1d.dat", ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user[bi].P));
			PetscViewerDestroy(viewer);  
			
			DAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		
		Velocity_Magnitude(user);
		DAVecGetArray(user[bi].da, user[bi].P, &p);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, x, &DIsDouble);
		DAVecRestoreArray(user[bi].da, user[bi].P, &p);
    
		/*
		DAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, x, &DIsDouble);
		DAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		*/
		
		PetscFree(x);
	}
	I = TECEND100();
	
	return 0;
}

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

  PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt	i, j, k, ia, ja, ka, ib, jb, kb;
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
  DAVecGetArray(cda, Csi, &csi);
  DAVecGetArray(cda, Eta, &eta);
  DAVecGetArray(cda, Zet, &zet);
  ierr = DAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DAGetGhostedCoordinates(da, &coords);
  DAVecGetArray(fda, coords, &coor);


  //  VecDuplicate(coords, &Cent);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  /* Calculating transformation metrics in i direction */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<lxe; i++){
	/* csi = de X dz */
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);
				       		    	    	 
	dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
		      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
		      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
		      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);
	  
	csi[k][j][i].x = dyde * dzdz - dzde * dydz;
	csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	csi[k][j][i].z = dxde * dydz - dyde * dxdz;

	
      }
    }
  }

  // Need more work -- lg65
  /* calculating j direction metrics */
  for (k=lzs; k<lze; k++){
    for (j=ys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	/* eta = dz X de */
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);
			    		         	 		   	 
	dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
	dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
	dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);
	  
	eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }

  /* calculating k direction metrics */
  for (k=zs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);
	  
	zet[k][j][i].x = dydc * dzde - dzdc * dyde;
	zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	zet[k][j][i].z = dxdc * dyde - dydc * dxde;

	/*	if ((i==1 || i==mx-2) && j==1 && (k==1 || k==0)) {
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxdc * dyde, dydc * dxde, dzdc);
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxde, dyde, dzde);
	  PetscPrintf(PETSC_COMM_WORLD, "Met %e %e %e\n", zet[k][j][i].x, zet[k][j][i].y, zet[k][j][i].z);
	  }*/
	
      }
    }
  }

  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	if (i>0 && j>0 && k>0) {
	  dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
			 coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
			 coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
			 coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

	  dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x - 
			 coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
	  dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y - 
			 coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
	  dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z - 
			 coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

	  dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
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


  //  PetscPrintf(PETSC_COMM_WORLD, "Local info: %d", info.mx);



  DAVecRestoreArray(cda, Csi, &csi);
  DAVecRestoreArray(cda, Eta, &eta);
  DAVecRestoreArray(cda, Zet, &zet);
  DAVecRestoreArray(da, Aj,  &aj);


  DAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);

  PetscBarrier(PETSC_NULL);
  return 0;
}

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
  PetscViewer	viewer;

  char filen2[128];

  PetscOptionsClearValue("-vecload_block_size");
  sprintf(filen2, "pfield%05d_%1d.dat", ti, user->_this);
  
  PetscViewer	pviewer;
  //Vec temp;
  PetscInt rank;
  PetscReal norm;
	
  if(file_exist(filen2))
  if(!onlyV) {
    //DACreateNaturalVector(user->da, &temp);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
	VecLoadIntoVector(pviewer, (user->P));
	VecNorm(user->P, NORM_INFINITY, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
	PetscViewerDestroy(pviewer);
	//VecDestroy(temp);
  }

  if(nv_once) sprintf(filen2, "nvfield%05d_%1d.dat", 0, user->_this);
  else sprintf(filen2, "nvfield%05d_%1d.dat", ti, user->_this);
  
  if(cs) sprintf(filen2, "cs_%05d_%1d.dat", ti, user->_this);
  
  if( !nv_once || (nv_once && ti==tis) )
  {
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
	VecLoadIntoVector(pviewer, (user->Nvert));
	PetscViewerDestroy(pviewer);  
  }
	
}

PetscErrorCode Ucont_P_Binary_Input1(UserCtx *user)
{
  PetscViewer viewer;
  char filen[128];
  
  sprintf(filen, "ufield%05d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoadIntoVector(viewer, (user->Ucat));
  PetscViewerDestroy(viewer);

  PetscBarrier(PETSC_NULL);
}

PetscErrorCode Ucont_P_Binary_Input_Averaging(UserCtx *user)
{
	PetscViewer viewer;
	char filen[128];
	/*
	sprintf(filen, "su0_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucat_sum));
	PetscViewerDestroy(viewer);

	sprintf(filen, "su1_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucat_cross_sum));
	PetscViewerDestroy(viewer);
	
	sprintf(filen, "su2_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucat_square_sum));
	PetscViewerDestroy(viewer);
	*/
	/*
	sprintf(filen, "sp_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->P));
	PetscViewerDestroy(viewer);
	*/
  
	if(pcr) {
		Vec Ptmp;
		VecDuplicate(user->P, &Ptmp);
		
		sprintf(filen, "pfield%05d_%1d.dat", ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoadIntoVector(viewer, user->P);
		PetscViewerDestroy(viewer);
		
		sprintf(filen, "sp_%05d_%1d.dat", ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoadIntoVector(viewer, Ptmp);
		PetscViewerDestroy(viewer);
		
		
		VecScale(Ptmp, -1./((double)tis+1.0));
		VecAXPY(user->P, 1., Ptmp);
		
		VecDestroy(Ptmp);
	}
	
	if(nv_once) sprintf(filen, "nvfield%05d_%1d.dat", 0, user->_this);
	else sprintf(filen, "nvfield%05d_%1d.dat", ti, user->_this);

	//if(cs) sprintf(filen2, "cs_%05d_%1d.dat", ti, user->_this);

	if( !nv_once || (nv_once && ti==tis) )
	  {
	    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	    VecLoadIntoVector(viewer, (user->Nvert));
	    PetscViewerDestroy(viewer);
	  }
	/*
	if( !nv_once || (nv_once && ti==tis) ) {
		sprintf(filen, "nvfield%05d_%1d.dat", ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoadIntoVector(viewer, (user->Nvert));
		PetscViewerDestroy(viewer);
	}
	*/
	PetscBarrier(PETSC_NULL);
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
	PetscTruth flag;
	
	DA	da, fda;
	Vec	qn, qnm;
	Vec	c;
	UserCtx	*user;

	PetscErrorCode ierr;
		
	IBMNodes	ibm, *ibm0, *ibm1;

	PetscInitialize(&argc, &argv, (char *)0, help);

	
	PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
	
	
	char tmp_str[256];
	PetscOptionsGetString(PETSC_NULL, "-prefix", tmp_str, 256, &flag);
	if(flag)sprintf(prefix, "%s_", tmp_str);
	else sprintf(prefix, "");
		
	PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-xyz", &xyz_input, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ransout", &rans_output, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-levelset", &levelset, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-avg", &avg, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-shear", &shear, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging_option, &flag);	// from control.dat
	
	PetscOptionsGetInt(PETSC_NULL, "-cs", &cs, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-ii_periodic", &i_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-jj_periodic", &j_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-kk_periodic", &k_periodic, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-nv", &nv_once, &flag);
	printf("nv_once=%d\n", nv_once);
	
	int QCR = 0;
	PetscOptionsGetInt(PETSC_NULL, "-qcr", &QCR, PETSC_NULL);
	
	PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
	if (!flag) PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
	
	PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, &flag);
	if (!flag) tie = tis;
    
	PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
	if (!flag) tsteps = 5; /* Default increasement is 5 */
    
	PetscOptionsGetInt(PETSC_NULL, "-onlyV", &onlyV, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-iavg", &i_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-javg", &j_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-kavg", &k_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-ikavg", &ik_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-pcr", &pcr, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-reynolds", &reynolds, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-ikcavg", &ikc_average, &flag);
	if(flag) {
		PetscTruth flag1, flag2;
		PetscOptionsGetInt(PETSC_NULL, "-pi", &pi, &flag1);
		PetscOptionsGetInt(PETSC_NULL, "-pk", &pk, &flag2);
		
		if(!flag1 || !flag2) {
			printf("To use -ikcavg you must set -pi and -pk, which are number of points in i- and k- directions.\n");
			exit(0);
		}
	}
	
  
	if(pcr) avg=1;
	if(i_average) avg=1;
	if(j_average) avg=1;
	if(k_average) avg=1;
	if(ik_average) avg=1;
	if(ikc_average) avg=1;
  
  
	if(i_average + j_average + k_average >1) PetscPrintf(PETSC_COMM_WORLD, "Iavg and Javg cannot be set together !! !\n"), exit(0);
      
	PetscInt rank, bi;

	PetscMalloc(sizeof(IBMNodes), &ibm0);
	PetscMalloc(sizeof(IBMNodes), &ibm1);

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(xyz_input) {block_number=1;}
	else {
		FILE *fd;
		fd = fopen("grid.dat", "r");
		if(binary_input) fread(&block_number, sizeof(int), 1, fd);
		else fscanf(fd, "%i\n", &block_number);
		MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		fclose(fd);
	}
	
	PetscMalloc(block_number*sizeof(UserCtx), &user);
	PetscOptionsGetReal(PETSC_NULL, "-ren", &user->ren, PETSC_NULL);

	ReadCoordinates(user);

	PetscPrintf(PETSC_COMM_WORLD, "read coord!\n");

	for (bi=0; bi<block_number; bi++) {
		DACreateGlobalVector(user[bi].da, &user[bi].Nvert);
		if(shear) {
			DACreateGlobalVector(user[bi].fda, &user[bi].Csi);
			DACreateGlobalVector(user[bi].fda, &user[bi].Eta);
			DACreateGlobalVector(user[bi].fda, &user[bi].Zet);
			DACreateGlobalVector(user[bi].da, &user[bi].Aj);
			FormMetrics(&(user[bi]));
			
			Calc_avg_shear_stress(&(user[bi]));
						
			VecDestroy(user[bi].Csi);
			VecDestroy(user[bi].Eta);
			VecDestroy(user[bi].Zet);
			VecDestroy(user[bi].Aj);
			exit(0);
		}
		else if(!avg) {
			DACreateGlobalVector(user[bi].da, &user[bi].P);
			DACreateGlobalVector(user[bi].fda, &user[bi].Ucat);
			if(!vc) DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_o);
			
			if(QCR)	{
				if(QCR==1 || QCR==2) {
					DACreateGlobalVector(user[bi].fda, &user[bi].Csi);
					DACreateGlobalVector(user[bi].fda, &user[bi].Eta);
					DACreateGlobalVector(user[bi].fda, &user[bi].Zet);
					DACreateGlobalVector(user[bi].da, &user[bi].Aj);
					FormMetrics(&(user[bi]));
				}
			}
		}
		else {
			if(pcr) {
				DACreateGlobalVector(user[bi].da, &user[bi].P);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat);
			}
			else if(avg==1) {
			  /*
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
			  */
			}
			else if(avg==2) {	// just compute k
			  /*
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
				DACreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
			  */
			}
		}
		
	}


  
	if(avg) {
		if(i_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in I direction!\n");
		else if(j_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in J direction!\n");
		else if(k_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in K direction!\n");
		else if(ik_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in IK direction!\n");
		else PetscPrintf(PETSC_COMM_WORLD, "Averaging !\n");
		/*
		DACreateGlobalVector(user[bi].fda, &user->Ucat_sum);
		DACreateGlobalVector(user[bi].fda, &user->Ucat_cross_sum);
		DACreateGlobalVector(user[bi].fda, &user->Ucat_square_sum);
		*/
		
	}
  
  
	for (ti=tis; ti<=tie; ti+=tsteps) {
		for (bi=0; bi<block_number; bi++) {
			if(avg) Ucont_P_Binary_Input_Averaging(&user[bi]);
			else {
				Ucont_P_Binary_Input(&user[bi]);
				Ucont_P_Binary_Input1(&user[bi]);
			}
		}
		
		if (!QCR) {
			if(avg) TECIOOut_Averaging(user);
			else TECIOOut_V(user, onlyV);
			//TECIOOut(user);
		}
		else {
			TECIOOutQ(user, QCR);
		}
	}
	PetscFinalize();
}



PetscErrorCode ReadCoordinates(UserCtx *user)
{
	Cmpnts ***coor;

	Vec Coor;
	PetscInt bi, i, j, k, rank, IM, JM, KM;
	PetscReal *gc;
	FILE *fd;
	PetscReal	d0 = 1.;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscReal	cl = 1.;
	PetscOptionsGetReal(PETSC_NULL, "-cl", &cl, PETSC_NULL);

	char str[256];
	
	if(xyz_input) sprintf(str, "xyz.dat");
	else sprintf(str, "grid.dat");
	
	fd = fopen(str, "r");
	
	if(fd==NULL) printf("Cannot open %s !\n", str),exit(0);

	printf("Begin Reading %s !\n", str);
	  
	if(xyz_input) {i=1;}
	else if(binary_input) {
		fread(&i, sizeof(int), 1, fd);
		if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a text file !\n"),exit(0);
	}
	else {
		fscanf(fd, "%i\n", &i);
		if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a binary file !\n"),exit(0);
	}
	

	for (bi=block_number-1; bi>=0; bi--) {
	  
		std::vector<double> X, Y,Z;
		double tmp;
		
		if(xyz_input) {
			fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
			X.resize(user[bi].IM);
			Y.resize(user[bi].JM);
			Z.resize(user[bi].KM);
			
			for (i=0; i<user[bi].IM; i++) fscanf(fd, "%le %le %le\n", &X[i], &tmp, &tmp);
			for (j=0; j<user[bi].JM; j++) fscanf(fd, "%le %le %le\n", &tmp, &Y[j], &tmp);
			for (k=0; k<user[bi].KM; k++) fscanf(fd, "%le %le %le\n", &tmp, &tmp, &Z[k]);
		}
		else if(binary_input) {
			fread(&(user[bi].IM), sizeof(int), 1, fd);
			fread(&(user[bi].JM), sizeof(int), 1, fd);
			fread(&(user[bi].KM), sizeof(int), 1, fd);
		}
		else fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
		
		IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    
    
		DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
				user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
				PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
				&(user[bi].da));
		if(rans) {
			DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
				user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
				PETSC_DECIDE, 2, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
				&(user[bi].fda2));
		}
		DASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
		DAGetCoordinateDA(user[bi].da, &(user[bi].fda));
	
		DAGetLocalInfo(user[bi].da, &(user[bi].info));

		DALocalInfo	info = user[bi].info;
		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
		
		DAGetGhostedCoordinates(user[bi].da, &Coor);
		DAVecGetArray(user[bi].fda, Coor, &coor);
		
		double buffer;
		
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].x = X[i]/cl;
				else coor[k][j][i].x = buffer/cl;
			}
		}
			
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].y = Y[j]/cl;
				else coor[k][j][i].y = buffer/cl;
			}
		}
	
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].z = Z[k]/cl;
				else coor[k][j][i].z = buffer/cl;
			}
		}

		DAVecRestoreArray(user[bi].fda, Coor, &coor);

		Vec	gCoor;
		DAGetCoordinates(user[bi].da, &gCoor);
		DALocalToGlobal(user[bi].fda, Coor, INSERT_VALUES, gCoor);
		DAGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
		DAGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);

	}
	
	fclose(fd);
	
	printf("Finish Reading %s !\n", str);
  
	for (bi=0; bi<block_number; bi++) {
		user[bi]._this = bi;
	}
	return(0);
}

void Calc_avg_shear_stress(UserCtx *user)
{
	double N=(double)tis+1.0;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***usum, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***psum, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

	char filen[256];
	PetscViewer	viewer;
		
	Vec P_sum;
	DACreateGlobalVector(user->da, &P_sum);
  DACreateGlobalVector(user->fda, &user->Ucat_sum);
  	  
  ti=tis;
  sprintf(filen, "su0_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucat_sum));
	PetscViewerDestroy(viewer);
	
	ti=tis;		
	sprintf(filen, "sp_%05d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, P_sum);
	PetscViewerDestroy(viewer);

  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);
  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->fda, user->Ucat_sum, &usum);
  DAVecGetArray(user->da, P_sum, &psum);


	double force_skin_bottom = 0;
	double force_pressure_bottom = 0;
	double force_bottom = 0;
	double area_bottom = 0;
	
	double force_skin_top = 0;
	double force_pressure_top = 0;
	double force_top = 0;
	double area_top = 0;
			
	j=0;
	for (k=lzs; k<lze; k++)
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j+1][i] < 0.1) {
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			
			dudc=0, dvdc=0, dwdc=0;
			
			dude=usum[k][j+1][i].x * 2.0 / N;
			dvde=usum[k][j+1][i].y * 2.0 / N;
			dwde=usum[k][j+1][i].z * 2.0 / N;
			
			dudz=0, dvdz=0, dwdz=0;
			
			double ajc = aj[k][j+1][i];
			double csi0 = csi[k][j+1][i].x, csi1 = csi[k][j+1][i].y, csi2 = csi[k][j+1][i].z;
			double eta0 = eta[k][j+1][i].x, eta1 = eta[k][j+1][i].y, eta2 = eta[k][j+1][i].z;
			double zet0 = zet[k][j+1][i].x, zet1 = zet[k][j+1][i].y, zet2 = zet[k][j+1][i].z;

			Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, 
					dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
					&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

			double j_area = sqrt( eta[k][j+1][i].x*eta[k][j+1][i].x + eta[k][j+1][i].y*eta[k][j+1][i].y + eta[k][j+1][i].z*eta[k][j+1][i].z );
			double ni[3], nj[3], nk[3];
			double nx, ny, nz;
			Calculate_normal(csi[k][j+1][i], eta[k][j+1][i], zet[k][j+1][i], ni, nj, nk);
			nx = nj[0]; //inward normal
			ny = nj[1]; //inward normal
			nz = nj[2]; //inward normal
			
			
			double Fp = - psum[k][j+1][i] * eta2 / N;
			double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
			//double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;
			
			force_skin_bottom += Fs;
			force_pressure_bottom += Fp;
			force_bottom += Fs + Fp;
			area_bottom += fabs(eta1);	// projected area
		}
	}

	j=my-2;
	for (k=lzs; k<lze; k++)
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j][i] < 0.1) {
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			
			dudc=0, dvdc=0, dwdc=0;
			
			dude = -usum[k][j][i].x * 2.0 / N;
			dvde = -usum[k][j][i].y * 2.0 / N;
			dwde = -usum[k][j][i].z * 2.0 / N;
			
			dudz=0, dvdz=0, dwdz=0;
			
			double ajc = aj[k][j][i];
			double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
			double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
			double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

			Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, 
					dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
					&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

			double j_area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double ni[3], nj[3], nk[3];
			double nx, ny, nz;
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			nx = -nj[0]; //inward normal
			ny = -nj[1]; //inward normal
			nz = -nj[2]; //inward normal
			
			
			double Fp = - psum[k][j][i] * eta2 / N;
			double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
			//double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;
			
			force_skin_top += Fs;
			force_pressure_top += Fp;
			force_top += Fs + Fp;
			area_top += fabs(eta1);	// projected area
		}
	}
	
  DAVecRestoreArray(user->fda, user->Csi, &csi);
  DAVecRestoreArray(user->fda, user->Eta, &eta);
  DAVecRestoreArray(user->fda, user->Zet, &zet);
  DAVecRestoreArray(user->da, user->Aj, &aj);
  DAVecRestoreArray(user->da, user->Nvert, &nvert);
  DAVecRestoreArray(user->fda, user->Ucat_sum, &usum);
  DAVecRestoreArray(user->da, P_sum, &psum);

	VecDestroy(P_sum);
  VecDestroy(user->Ucat_sum);
  
  printf("Top:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
  			area_top, force_top, force_skin_top, force_pressure_top);
  			
	printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
  			force_top/area_top, force_skin_top/area_top, force_pressure_top/area_top,
  			sqrt(fabs(force_top/area_top)), sqrt(fabs(force_top/area_top))*user->ren);
  			
  printf("\n");
  
  printf("Bottom:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
  			area_bottom, force_bottom, force_skin_bottom, force_pressure_bottom);
  			
	printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
  			force_bottom/area_bottom, force_skin_bottom/area_bottom, force_pressure_bottom/area_bottom,
  			sqrt(fabs(force_bottom/area_bottom)), sqrt(fabs(force_bottom/area_bottom))*user->ren);
}

PetscErrorCode Lambda2(UserCtx *user)
{
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);

  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  //PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	      
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double ajc = aj[k][j][i];
	
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
							
	double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);

	
	w11 = 0;
	w12 = 0.5*(du_dy - dv_dx);
	w13 = 0.5*(du_dz - dw_dx);
	w21 = -w12;
	w22 = 0.;
	w23 = 0.5*(dv_dz - dw_dy);
	w31 = -w13;
	w32 = -w23;
	w33 = 0.;
	
	
	double S[3][3], W[3][3], D[3][3];
	
	D[0][0] = du_dx, D[0][1] = du_dy, D[0][2] = du_dz;
	D[1][0] = dv_dx, D[1][1] = dv_dy, D[1][2] = dv_dz;
	D[2][0] = dw_dx, D[2][1] = dw_dy, D[2][2] = dw_dz;
	
	S[0][0] = Sxx;
	S[0][1] = Sxy;
	S[0][2] = Sxz;

	S[1][0] = Syx;
	S[1][1] = Syy;
	S[1][2] = Syz;

	S[2][0] = Szx;
	S[2][1] = Szy;
	S[2][2] = Szz;

	W[0][0] = w11;
	W[0][1] = w12;
	W[0][2] = w13;
	W[1][0] = w21;
	W[1][1] = w22;
	W[1][2] = w23;
	W[2][0] = w31;
	W[2][1] = w32;
	W[2][2] = w33;
	
	// lambda-2
	double A[3][3], V[3][3], d[3];
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) A[row][col]=0;
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		A[row][col] += S[row][0] * S[0][col];
		A[row][col] += S[row][1] * S[1][col];
		A[row][col] += S[row][2] * S[2][col];
	}
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		A[row][col] += W[row][0] * W[0][col];
		A[row][col] += W[row][1] * W[1][col];
		A[row][col] += W[row][2] * W[2][col];
	}
	
	if(nvert[k][j][i]<0.1) {
		eigen_decomposition(A, V, d);
		q[k][j][i] = d[1];
	}
	else q[k][j][i] = 1000.0;
/*	
	// delta criterion
	double DD[3][3];
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) DD[row][col]=0;
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		DD[row][col] += D[row][0] * D[0][col];
		DD[row][col] += D[row][1] * D[1][col];
		DD[row][col] += D[row][2] * D[2][col];
	}
	double tr_DD = DD[0][0] + DD[1][1] + DD[2][2];
	double det_D = D[0][0]*(D[2][2]*D[1][1]-D[2][1]*D[1][2])-D[1][0]*(D[2][2]*D[0][1]-D[2][1]*D[0][2])+D[2][0]*(D[1][2]*D[0][1]-D[1][1]*D[0][2]);
	
	//double Q = -0.5*tr_DD;
	
	double SS=0, WW=0;
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		SS+=S[row][col]*S[row][col];
		WW+=W[row][col]*W[row][col];
	}
	double Q = 0.5*(WW - SS);
	
	double R = - det_D;
	if(nvert[k][j][i]<0.1) {
		q[k][j][i] = pow( 0.5*R, 2. )  + pow( Q/3., 3.);
	}
	else q[k][j][i] = -10;
	if(q[k][j][i]<0) q[k][j][i]=-10;
*/	
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->fda, user->Csi, &csi);
  DAVecRestoreArray(user->fda, user->Eta, &eta);
  DAVecRestoreArray(user->fda, user->Zet, &zet);

  DAVecRestoreArray(user->da, user->Aj, &aj);
  DAVecRestoreArray(user->da, user->Nvert, &nvert);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode QCriteria(UserCtx *user)
{

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;
	
  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);

  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double ajc = aj[k][j][i];
	
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
							
	double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	so = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz;
	
	w11 = 0;
	w12 = 0.5*(du_dy - dv_dx);
	w13 = 0.5*(du_dz - dw_dx);
	w21 = -w12;
	w22 = 0.;
	w23 = 0.5*(dv_dz - dw_dy);
	w31 = -w13;
	w32 = -w23;
	w33 = 0.;
	
	wo = w11*w11 + w12*w12 + w13*w13 + w21*w21 + w22*w22 + w23*w23 + w31*w31 + w32*w32 + w33*w33;

/*
	so = ( d11 *  d11 + d22 * d22 + d33 * d33) + 0.5* ( (d12 + d21) * (d12 + d21) + (d13 + d31) * (d13 + d31) + (d23 + d32) * (d23 + d32) );
	wo = 0.5 * ( (d12 - d21)*(d12 - d21) + (d13 - d31) * (d13 - d31) + (d23 - d32) * (d23 - d32) );
	V19=0.5 * ( (V13 - V11)*(V13 - V11) + (V16 - V12) * (V16 - V12) + (V17 - V15) * (V17 - V15) ) - 0.5 * ( V10 *  V10 + V14 * V14 + V18 * V18) - 0.25* ( (V13 + V11) * (V13 + V11) + (V16 + V12) * (V16 + V12) + (V17 + V15) * (V17 + V15) )
*/
	
	if( nvert[k][j][i]>0.1 ) q[k][j][i] = -100;
	else q[k][j][i] = (wo - so) / 2.;
      }
    }
  }

	DAVecRestoreArray(user->fda, user->Ucat, &ucat);
	DAVecRestoreArray(user->fda, user->Csi, &csi);
	DAVecRestoreArray(user->fda, user->Eta, &eta);
	DAVecRestoreArray(user->fda, user->Zet, &zet);

	DAVecRestoreArray(user->da, user->Aj, &aj);
	DAVecRestoreArray(user->da, user->Nvert, &nvert);
	DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode Velocity_Magnitude(UserCtx *user)	// store at P
{
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	q[k][j][i] = sqrt( ucat[k][j][i].x*ucat[k][j][i].x + ucat[k][j][i].y*ucat[k][j][i].y + ucat[k][j][i].z*ucat[k][j][i].z );
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}


PetscErrorCode ibm_read(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal	t, dr;
  double xt;
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata0", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
    n_v =0;
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%le", &xt);
    ibm->n_v = n_v;

    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));

    for (i=0; i<n_v; i++) {
      fscanf(fd, "%le %le %le %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i], &t, &t, &t);
      x_bp[i] = x_bp[i] / 28.;
      y_bp[i] = y_bp[i] / 28.;
      z_bp[i] = z_bp[i] / 28.;
    }
    ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

    for (i=0; i<n_v; i++) {
      PetscReal temp;
      temp = ibm->y_bp0[i];
      ibm->y_bp0[i] = ibm->z_bp0[i];
      ibm->z_bp0[i] = -temp;
    }


    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fscanf(fd, "%i\n", &n_elmt);
    ibm->n_elmt = n_elmt;
    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }
    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

    fclose(fd);

    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

      
    }
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

/*     for (i=0; i<n_elmt; i++) { */
/*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
/*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode ibm_read_ucd(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal	t, dr;
  PetscInt 	temp;
  double xt;
  char string[128];
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata", "r");
    if (!fd) SETERRQ(1, "Cannot open IBM node file")
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);

      fscanf(fd, "%i %i %i %i %i\n", &n_v, &n_elmt, &temp, &temp, &temp);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;

      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

      
      PetscReal cl = 1.;

      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
      
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &temp, &x_bp[i], &y_bp[i], &z_bp[i]);
	x_bp[i] = x_bp[i] / cl;
	y_bp[i] = y_bp[i] / cl;
	z_bp[i] = z_bp[i] / cl;
	
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;
      }
      ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

	

      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      char str[20];
      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &temp, &temp, str, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      fclose(fd);
    }
    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

      
    }
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    /*     for (i=0; i<n_elmt; i++) { */
    /*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
    /*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    for (i=0; i<ibm->n_v; i++) {
      ibm->x_bp[i] = ibm->x_bp0[i];
      ibm->y_bp[i] = ibm->y_bp0[i];
      ibm->z_bp[i] = ibm->z_bp0[i];

      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode Combine_Elmt(IBMNodes *ibm, IBMNodes *ibm0, IBMNodes *ibm1)
{

  PetscInt i;

  ibm->n_v = ibm0->n_v + ibm1->n_v;
  ibm->n_elmt = ibm0->n_elmt + ibm1->n_elmt;

  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  for (i=0; i<ibm0->n_v; i++) {
    ibm->x_bp[i] = ibm0->x_bp[i];
    ibm->y_bp[i] = ibm0->y_bp[i];
    ibm->z_bp[i] = ibm0->z_bp[i];

    ibm->u[i] = ibm0->u[i];
    ibm->uold[i] = ibm0->uold[i];
    //    ibm->u[i].x = 0.;
/*     PetscPrintf(PETSC_COMM_WORLD, "Vel %e %e %e\n", ibm->u[i].x, ibm->u[i].y, ibm->u[i].z); */
  }
  for (i=0; i<ibm0->n_elmt; i++) {
    ibm->nv1[i] = ibm0->nv1[i];
    ibm->nv2[i] = ibm0->nv2[i];
    ibm->nv3[i] = ibm0->nv3[i];

    ibm->nf_x[i] = ibm0->nf_x[i];
    ibm->nf_y[i] = ibm0->nf_y[i];
    ibm->nf_z[i] = ibm0->nf_z[i];
  }

  for (i=ibm0->n_v; i<n_v; i++) {
    ibm->x_bp[i] = ibm1->x_bp[i-ibm0->n_v];
    ibm->y_bp[i] = ibm1->y_bp[i-ibm0->n_v];
    ibm->z_bp[i] = ibm1->z_bp[i-ibm0->n_v];
    ibm->u[i].x = 0.;
    ibm->u[i].y = 0.;
    ibm->u[i].z = 0.;
  }

  for (i=ibm0->n_elmt; i<n_elmt; i++) {
    ibm->nv1[i] = ibm1->nv1[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv2[i] = ibm1->nv2[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv3[i] = ibm1->nv3[i-ibm0->n_elmt] + ibm0->n_v;

    ibm->nf_x[i] = ibm1->nf_x[i-ibm0->n_elmt];
    ibm->nf_y[i] = ibm1->nf_y[i-ibm0->n_elmt];
    ibm->nf_z[i] = ibm1->nf_z[i-ibm0->n_elmt];
  }

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%05d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 1670-96);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=96; i<n_elmt; i++) {
	if (fabs(ibm->nf_z[i]) > 0.5 ||
	    (fabs(ibm->nf_z[i]) < 0.5 &&
	     (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
	      ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	}
      }
      fclose(f);

      sprintf(filen, "leaflet%05d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 96);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<96; i++) {
	if (fabs(ibm->nf_z[i]) > 0.5 ||
	    (fabs(ibm->nf_z[i]) < 0.5 &&
	     (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
	      ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	}
      }
      fclose(f);

    }
  }

  return 0;
}

PetscErrorCode Elmt_Move(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
  for (i=0; i<n_elmt; i++) {
      n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
      dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
      dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
      dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];

      dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e];
      dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e];
      dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e];

      ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
      ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] +
		ibm->nf_z[i]*ibm->nf_z[i]);

      ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

/*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
      //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);

  }
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      //      ibm->uold[i] = ibm->u[i];

      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}

PetscErrorCode Elmt_Move1(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
  for (i=0; i<n_elmt; i++) {
      n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
      dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
      dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
      dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];

      dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e];
      dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e];
      dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e];

      ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
      ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] +
		ibm->nf_z[i]*ibm->nf_z[i]);

      ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

/*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
      //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);

  }
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      ibm->uold[i] = ibm->u[i];

      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}


/*****************************************************************/
#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (int i = 0; i < n-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
  double e[n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);
}
