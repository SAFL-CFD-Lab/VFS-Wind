/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

Vec	LevelSet,	LevelSet0, LevelSet_o;
double dtau;
//double thickness=0.02;
extern double	sign(double	a);
extern double	M(double a,	double b);
extern int immersed, NumberOfBodies;
extern int i_periodic, j_periodic, k_periodic;

void Compute_dlevel_center_levelset	(int i,	int	j, int k,	 int mx, int my, int mz, double	sgn, int wall_distance,	PetscReal	***level,	PetscReal	***nvert,	double *dldc,	double *dlde,	double *dldz);

void Init_LevelSet_Vectors(UserCtx *user)
{
	VecDuplicate(user->P,	&LevelSet);
	VecDuplicate(user->P,	&LevelSet0);
	VecDuplicate(user->P,	&LevelSet_o);
	//VecDuplicate(user->lP, &lLevelSet);
};

void Destroy_LevelSet_Vectors(UserCtx	*user)
{
	VecDestroy (LevelSet);
	VecDestroy (LevelSet0);
	VecDestroy (LevelSet_o);
	//VecDestroy (lLevelSet);
};

void Distance_Function_RHS (UserCtx *user, Vec Levelset_RHS, int wall_distance)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs,	xe,	ys,	ye,	zs,	ze;
	PetscInt	mx,	my,	mz;
	PetscInt	i, j,	k;
	PetscReal	dpdc,	dpde,	dpdz;
	
	Vec		Csi	=	user->lCsi,	Eta	=	user->lEta,	Zet	=	user->lZet;
	Vec		Aj	=	user->lAj;

	Cmpnts	***csi,	***eta,	***zet;
	PetscReal	***aj, ***level, ***level0,	***rhs,	***grad_level;

	PetscInt	lxs, lys,	lzs, lxe,	lye, lze;
	PetscReal	***nvert;
	
	Vec	L, lLevelset0;	// grad	level, level0
	
	xs = info.xs;	xe = xs	+	info.xm;
	ys = info.ys;	ye = ys	+	info.ym;
	zs = info.zs;	ze = zs	+	info.zm;
	
	lxs	=	xs;	lxe	=	xe;
	lys	=	ys;	lye	=	ye;
	lzs	=	zs;	lze	=	ze;
	
	mx = info.mx;	my = info.my;	mz = info.mz;
	
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx)	lxe	=	xe-1;
	if (ye==my)	lye	=	ye-1;
	if (ze==mz)	lze	=	ze-1;
	
	VecSet(Levelset_RHS, 0);
	
	VecDuplicate(user->lP, &L);
	VecDuplicate(user->lP, &lLevelset0);
	
	DAGlobalToLocalBegin(user->da, LevelSet0, INSERT_VALUES, lLevelset0);
	DAGlobalToLocalEnd(user->da, LevelSet0,	INSERT_VALUES, lLevelset0);
	
	DAVecGetArray(fda, Csi,	&csi);
	DAVecGetArray(fda, Eta,	&eta);
	DAVecGetArray(fda, Zet,	&zet);
	DAVecGetArray(da,	 Aj,	&aj);
	DAVecGetArray(da,	user->lNvert,	&nvert);
	DAVecGetArray(da,	user->lLevelset, &level);
	DAVecGetArray(da,	lLevelset0,	&level0);
	DAVecGetArray(da,	Levelset_RHS,	&rhs);
	
	DAVecGetArray(da,  L,  &grad_level);
	for	(k=lzs;	k<lze; k++)
	for	(j=lys;	j<lye; j++)
	for	(i=lxs;	i<lxe; i++) {
		double dldc, dlde, dldz;
		double dl_dx, dl_dy, dl_dz;
		
		double csi0=csi[k][j][i].x,csi1=csi[k][j][i].y, csi2=csi[k][j][i].z;
		double eta0=eta[k][j][i].x,eta1=eta[k][j][i].y, eta2=eta[k][j][i].z;
		double zet0=zet[k][j][i].x,zet1=zet[k][j][i].y, zet2=zet[k][j][i].z;
		double ajc = aj[k][j][i];
		
		double dx=pow(1./aj[k][j][i],1./3.);
		if(dthick_set) dx = dthick;

		double sgn = sign1(level0[k][j][i], dx);
		//if(wall_distance)	sgn	=	sign(level0[k][j][i]);
		
		//Compute_dlevel_center_levelset (i, j, k, mx, my, mz, sgn, wall_distance, level, nvert, &dldc, &dlde, &dldz);
		Compute_dlevel_center_levelset (i, j, k, mx, my, mz, sign(level0[k][j][i]), wall_distance, level, nvert, &dldc, &dlde, &dldz); //100521
		
		Compute_dscalar_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		
		grad_level[k][j][i] = sqrt( dl_dx*dl_dx	+ dl_dy*dl_dy +	dl_dz*dl_dz );
		if(nvert[k][j][i]>0.1) grad_level[k][j][i]=0;
	}
	DAVecRestoreArray(da,  L,  &grad_level);

	DALocalToLocalBegin(user->da, L, INSERT_VALUES,	L);
	DALocalToLocalEnd(user->da, L, INSERT_VALUES, L);
	
	DAVecGetArray(da,  L,  &grad_level);
	
	// Neumann,	periodic conditions
	if(xs==0 ||	xe==mx)	{
		int	from,	to;
		for	(k=lzs;	k<lze; k++)
		for	(j=lys;	j<lye; j++)	{
			if(xs==0)	{
				i	=	1, from	=	i, to	=	0;
				
				if(i_periodic) from	=	mx-2;
				else if(ii_periodic) from	=	-2;
				
				grad_level[k][j][to] = grad_level[k][j][from];
			}
			
			if(xe==mx) {
				i	=	mx-2,	from = i,	to = mx-1;
				
				if(i_periodic) from	=	1;
				else if(ii_periodic) from	=	mx+1;
				
				grad_level[k][j][to] = grad_level[k][j][from];
			}
		}
	}
	
	if(ys==0 ||	ye==my)	{
		int	from,	to;
				
		for	(k=lzs;	k<lze; k++)
		for	(i=lxs;	i<lxe; i++)	{
			if(ys==0)	{
				j	=	1, from	=	j, to	=	0;
				
				if(j_periodic) from	=	my-2;
				else if(jj_periodic) from	=	-2;
				
				grad_level[k][to][i] = grad_level[k][from][i];
			}
			
			if(ye==my) {
				j	=	my-2,	from = j,	to = my-1;
				
				if(j_periodic) from	=	1;
				else if(jj_periodic) from	=	my+1;
				
				grad_level[k][to][i] = grad_level[k][from][i];
			}
		}
	}
	
	if(zs==0 ||	ze==mz)	{
		int	from,	to;
		
		for	(j=lys;	j<lye; j++)
		for	(i=lxs;	i<lxe; i++)	{
			if(zs==0)	{
				k	=	1, from	=	k, to	=	0;
				
				if(k_periodic) from	=	mz-2;
				else if(kk_periodic) from	=	-2;
				
				grad_level[to][j][i] = grad_level[from][j][i];
			}
			
			if(ze==mz) {
				k	=	mz-2,	from = k,	to = mz-1;
				
				if(k_periodic) from	=	1;
				else if(kk_periodic) from	=	mz+1;
				
				grad_level[to][j][i] = grad_level[from][j][i];
			}
		}
	}
	
	DAVecRestoreArray(da,  L,  &grad_level);
	DALocalToLocalBegin(user->da, L, INSERT_VALUES,	L);
	DALocalToLocalEnd(user->da, L, INSERT_VALUES, L);
	DAVecGetArray(da,  L,  &grad_level);

	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++)	{
		if (i<= 0 || i>= mx-1 || j<=0 || j>=my-1 || k<=0 || k>=mz-1 || nvert[k][j][i]>1.1){
			rhs[k][j][i]=0.;
			continue;
		}
		
		if(nvert[k][j][i]>0.1) {
		  rhs[k][j][i]=0.;
		  continue;
		}
		
		if(wall_distance) {
			if(nvert[k][j][i]>0.1) { rhs[k][j][i]=0.; continue; }
			if(i <= 1 && (user->bctype[0]==1 || user->bctype[0]==-1 || user->bctype[0]==-2)) { rhs[k][j][i]=0.; continue; }
			if(i >=	mx-2 &&	(user->bctype[1]==1 || user->bctype[1]==-1 || user->bctype[1]==-2)) { rhs[k][j][i]=0.; continue; }
			if(j <=1 && (user->bctype[2]==1 || user->bctype[2]==-1 || user->bctype[2]==-2)) { rhs[k][j][i]=0.; continue; }
			if(j >=my-2 && (user->bctype[3]==1 || user->bctype[3]==-1 || user->bctype[3]==-2 || user->bctype[3]==12)) { rhs[k][j][i]=0.; continue; }
			if(k<=1	&& (user->bctype[4]==1 || user->bctype[4]==-1 || user->bctype[4]==-2)) { rhs[k][j][i]=0.; continue; }
			if(k>=mz-2 && (user->bctype[5]==1 || user->bctype[5]==-1 || user->bctype[5]==-2)){ rhs[k][j][i]=0.; continue; }
		}		
		else if( !wall_distance	&& user->bctype[4]==5 && user->bctype[5]==4) {
		  //if ( fix_inlet && k==1 ) { rhs[k][j][i] = 0; continue; }
		  //haha if ( fix_outlet && k==mz-2 ) { rhs[k][j][i] = 0; continue; } // important to stabilize outlet
		}
		
		double dx=pow(1./aj[k][j][i],1./3.);
		if(dthick_set) dx = dthick;
		double sgn = sign1(level0[k][j][i],dx);
		//if(wall_distance)	sgn	=	sign(level0[k][j][i]);
		
		double denom[3][3][3], num[3][3][3], weight[3][3][3];
		
		for(int	p=-1;	p<=1;	p++)
		for(int	q=-1;	q<=1;	q++)
		for(int	r=-1;	r<=1;	r++) {
			int	R=r+1, Q=q+1,	P=p+1;
			int	K=k+r, J=j+q,	I=i+p;
			double phi = level[K][J][I], grad	=	grad_level[K][J][I], dx=pow(1./aj[K][J][I],1./3.);
			if(dthick_set) dx	=	dthick;

			double f = dH(phi,dx)	*	grad;
			
			double _sgn	=	sign1( level0[K][J][I],	dx );
			//if(wall_distance)	_sgn = sign(level0[K][J][I]);
			
			num[R][Q][P] = dH(phi,dx) * _sgn * ( 1.	- grad );
			denom[R][Q][P] = dH(phi,dx) * f;
		}
		
		for(int	p=-1;	p<=1;	p++)
		for(int	q=-1;	q<=1;	q++)
		for(int	r=-1;	r<=1;	r++) {
		  int	R=r+1, Q=q+1,	P=p+1;
		  int	K=k+r, J=j+q,	I=i+p;
		  if( (!i_periodic && !ii_periodic && (I==0 || I==mx-1 ) ) ||
			(!j_periodic && !jj_periodic && (J==0 || J==my-1 ) ) || 
			(!k_periodic && !kk_periodic && (K==0 || K==mz-1) ) ||	
			nvert[K][J][I]>0.1) {
		    num[R][Q][P] = num[1][1][1];
		    denom[R][Q][P] = denom[1][1][1];
		  }
		}
		
		get_weight (i, j, k, mx, my, mz, aj, nvert, 0.1, weight);
		
		double numerator = integrate_testfilter(num, weight);
		double denominator = integrate_testfilter(denom, weight);
		
		double correction;
		
		if(	fabs(denominator)<1.e-10 ) correction=0;
		else {
			double grad	=	grad_level[k][j][i];
			double phi = level[k][j][i];
			double dx=pow(1./aj[k][j][i],1./3.);
			if(dthick_set) dx	=	dthick;

			double f = dH(phi,dx)	*	grad;
			correction = - numerator / denominator;
			correction *=	dH(phi,dx) * grad;
		}
		
		double dlevel_dx, dlevel_dy, dlevel_dz;
		double dldc, dlde, dldz;
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double ajc = aj[k][j][i];
		
		//Compute_dlevel_center_levelset (i, j, k, mx, my, mz, sgn, wall_distance, level, nvert, &dldc, &dlde, &dldz);
		Compute_dlevel_center_levelset (i, j, k, mx, my, mz, sign(level0[k][j][i]), wall_distance, level, nvert, &dldc, &dlde, &dldz); //100521

		Compute_dscalar_dxyz (csi0, csi1, csi2,	eta0, eta1, eta2, zet0,	zet1, zet2, ajc, dldc, dlde, dldz, &dlevel_dx, &dlevel_dy, &dlevel_dz);
		rhs[k][j][i] = sgn * ( 1. - sqrt( dlevel_dx*dlevel_dx +	dlevel_dy*dlevel_dy + dlevel_dz*dlevel_dz ) );
			
		if(nvert[k][j][i+1]+nvert[k][j][i-1]+nvert[k][j+1][i]+nvert[k][j-1][i]+nvert[k+1][j][i]+nvert[k-1][j][i]>0.1) {
		  // correction = 0;
		}

		if(!wall_distance) rhs[k][j][i] += correction;	// Sussman Fetami
	}

	DAVecRestoreArray(da,	 L,	 &grad_level);
	DAVecRestoreArray(fda, Csi,	&csi);
	DAVecRestoreArray(fda, Eta,	&eta);
	DAVecRestoreArray(fda, Zet,	&zet);
	DAVecRestoreArray(da,	 Aj,	&aj);
	DAVecRestoreArray(da,	user->lNvert,	&nvert);
	DAVecRestoreArray(da,	user->lLevelset, &level);
	DAVecRestoreArray(da,	lLevelset0,	&level0);
	DAVecRestoreArray(da,	Levelset_RHS,	&rhs);
	
	VecDestroy(L);
	VecDestroy(lLevelset0);
}

void Distance_Function_IC(UserCtx	*user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs,	xe,	ys,	ye,	zs,	ze;	// Local grid	information
	PetscInt	mx,	my,	mz;	// Dimensions	in three directions
	PetscInt	i, j,	k, ibi;
	PetscInt	lxs, lxe,	lys, lye,	lzs, lze;
	Cmpnts	***csi,	***eta,	***zet;
	PetscReal	***nvert,	***L0, ***aj;

	DAGetLocalInfo(da, &info);
	mx = info.mx;	my = info.my;	mz = info.mz;
	xs = info.xs;	xe = xs	+	info.xm;
	ys = info.ys;	ye = ys	+	info.ym;
	zs = info.zs;	ze = zs	+	info.zm;

	lxs	=	xs;	lxe	=	xe;
	lys	=	ys;	lye	=	ye;
	lzs	=	zs;	lze	=	ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx)	lxe	=	xe-1;
	if (ye==my)	lye	=	ye-1;
	if (ze==mz)	lze	=	ze-1;

	DAVecGetArray(da,	user->lNvert,	&nvert);
	DAVecGetArray(da,	LevelSet0, &L0);
	
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da,	user->lAj, &aj);
	
	
	if(immersed)
	for(ibi=0; ibi<NumberOfBodies; ibi++)
	{
		//IBMNodes *ibm	=	&ibm_ptr[ibi];
		IBMListNode	*current;
		current	=	user->ibmlist[ibi].head;
		while	(current) {
			IBMInfo	*ibminfo = &current->ibm_intp;
			current	= current->next;
			double sb = ibminfo->d_s;
			i=ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
			L0[k][j][i] = sb;
		};
	}	
	
	for	(k=zs; k<ze; k++)
	for	(j=ys; j<ye; j++)
	for	(i=xs; i<xe; i++) {
		int	flag=0;
		double dist=1.e10;
		
		if(nvert[k][j][i]>0.1) continue;
		
		if((user->bctype[0]==1|| user->bctype[0]==-1 || user->bctype[0]==-2) && i==1){
			double area=sqrt(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k][j][i-1] = dist;
		}
		
		
		if((user->bctype[1]==1|| user->bctype[1]==-1 || user->bctype[1]==-2) && i==mx-2 ) {
			double area=sqrt(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k][j][i+1] = dist;
		}
		if((user->bctype[2]==1|| user->bctype[2]==-1 || user->bctype[2]==-2) && j==1){
			double area=sqrt(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k][j-1][i] = dist;
		}
		if((user->bctype[3]==1|| user->bctype[3]==-1 || user->bctype[3]==-2 || user->bctype[3]==12) &&j==my-2){
			double area=sqrt(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k][j+1][i] = dist;
		}
		if((user->bctype[4]==1|| user->bctype[4]==-1 || user->bctype[4]==-2) && k==1){
			double area=sqrt(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k-1][j][i] = dist;
		}
		if((user->bctype[5]==1 || user->bctype[5]==-1 || user->bctype[5]==-2) && k==mz-2 ) {
			double area=sqrt(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z);
			dist = PetscMin(dist,0.5/aj[k][j][i]/area);
			flag=1;
			L0[k+1][j][i] = dist;
		}

		if(flag) L0[k][j][i] = dist;
		else L0[k][j][i] = 0.05;
		
		//if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) L0[k][j][i]=0;
	}
	

	
	
	DAVecRestoreArray(da,	user->lNvert,	&nvert);
	DAVecRestoreArray(da,	LevelSet0, &L0);
	
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(da,	user->lAj, &aj);
	
	//DALocalToGlobal(da,	lLevelSet0,	INSERT_VALUES, LevelSet);	// copy
	//VecCopy	(LevelSet, LevelSet_o);	
	VecCopy	(LevelSet0,	LevelSet); 
	VecCopy	(LevelSet0,	LevelSet_o); 
};

PetscErrorCode FormFunction_Distance(SNES	snes,	Vec	L, Vec Rhs,	void *ptr)
{
	UserCtx	*user	=	(UserCtx*)ptr;

	DALocalInfo	info;
	PetscInt	xs,	xe,	ys,	ye,	zs,	ze;	// Local grid	information
	PetscInt	mx,	my,	mz;	// Dimensions	in three directions
	PetscInt	i, j,	k;
	PetscInt	lxs, lxe,	lys, lye,	lzs, lze;
	PetscReal	***level;

	DAGetLocalInfo(user->da, &info);
	mx = info.mx;	my = info.my;	mz = info.mz;
	xs = info.xs;	xe = xs	+	info.xm;
	ys = info.ys;	ye = ys	+	info.ym;
	zs = info.zs;	ze = zs	+	info.zm;

	lxs	=	xs;	lxe	=	xe;
	lys	=	ys;	lye	=	ye;
	lzs	=	zs;	lze	=	ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx)	lxe	=	xe-1;
	if (ye==my)	lye	=	ye-1;
	if (ze==mz)	lze	=	ze-1;
	
	DAGlobalToLocalBegin(user->da, L,	INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, L,	INSERT_VALUES, user->lLevelset);
	
	DAVecGetArray(user->da,	user->lLevelset, &level);
	
	if(xs==0 ||	xe==mx)	{
		int	from,	to;
		for	(k=lzs;	k<lze; k++)
		for	(j=lys;	j<lye; j++)	{
			if(xs==0)	{
				i	=	1, from	=	i, to	=	0;
				
				if(i_periodic) from	=	mx-2;
				else if(ii_periodic) from	=	-2;
				
				level[k][j][to]	=	level[k][j][from];
			}
			
			if(xe==mx) {
				i	=	mx-2,	from = i,	to = mx-1;
				
				if(i_periodic) from	=	1;
				else if(ii_periodic) from	=	mx+1;
				
				level[k][j][to]	=	level[k][j][from];
			}
		}
	}
	
	if(ys==0 ||	ye==my)	{
		int	from,	to;
				
		for	(k=lzs;	k<lze; k++)
		for	(i=lxs;	i<lxe; i++)	{
			if(ys==0)	{
				j	=	1, from	=	j, to	=	0;
				
				if(j_periodic) from	=	my-2;
				else if(jj_periodic) from	=	-2;
				
				level[k][to][i]	=	level[k][from][i];
			}
			
			if(ye==my) {
				j	=	my-2,	from = j,	to = my-1;
				
				if(j_periodic) from	=	1;
				else if(jj_periodic) from	=	my+1;
				
				level[k][to][i]	=	level[k][from][i];
			}
		}
	}
	
	if(zs==0 ||	ze==mz)	{
		int	from,	to;
		
		for	(j=lys;	j<lye; j++)
		for	(i=lxs;	i<lxe; i++)	{
			if(zs==0)	{
				k	=	1, from	=	k, to	=	0;
				
				if(k_periodic) from	=	mz-2;
				else if(kk_periodic) from	=	-2;
				
				level[to][j][i]	=	level[from][j][i];
			}
			
			if(ze==mz) {
				k	=	mz-2,	from = k,	to = mz-1;
				
				if(k_periodic) from	=	1;
				else if(kk_periodic) from	=	mz+1;
				
				level[to][j][i]	=	level[from][j][i];
			}
		}
	}
	
	DAVecRestoreArray(user->da,	user->lLevelset, &level);
	
	Distance_Function_RHS(user,	Rhs, 1);
	VecAXPY(Rhs, -1/dtau,	L);
	VecAXPY(Rhs, 1/dtau, LevelSet_o);
	
	return(0);
}

void Solve_Distance_Explicit(UserCtx *user)
{
	DALocalInfo	info;
	PetscInt	xs,	xe,	ys,	ye,	zs,	ze;	// Local grid	information
	PetscInt	mx,	my,	mz;	// Dimensions	in three directions
	PetscInt	i, j,	k;
	PetscInt	lxs, lxe,	lys, lye,	lzs, lze;
	PetscReal	***level;

	DAGetLocalInfo(user->da, &info);
	mx = info.mx;	my = info.my;	mz = info.mz;
	xs = info.xs;	xe = xs	+	info.xm;
	ys = info.ys;	ye = ys	+	info.ym;
	zs = info.zs;	ze = zs	+	info.zm;

	lxs	=	xs;	lxe	=	xe;
	lys	=	ys;	lye	=	ye;
	lzs	=	zs;	lze	=	ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx)	lxe	=	xe-1;
	if (ye==my)	lye	=	ye-1;
	if (ze==mz)	lze	=	ze-1;
	
	DAGlobalToLocalBegin(user->da, LevelSet, INSERT_VALUES,	user->lLevelset);
	DAGlobalToLocalEnd(user->da, LevelSet, INSERT_VALUES,	user->lLevelset);
	
	DAVecGetArray(user->da,	user->lLevelset, &level);
	
	if(xs==0 ||	xe==mx)	{
		int	from,	to;
		for	(k=lzs;	k<lze; k++)
		for	(j=lys;	j<lye; j++)	{
			if(xs==0)	{
				i	=	1, from	=	i, to	=	0;
				
				if(i_periodic) from	=	mx-2;
				else if(ii_periodic) from	=	-2;
				
				level[k][j][to]	=	level[k][j][from];
			}
			
			if(xe==mx) {
				i	=	mx-2,	from = i,	to = mx-1;
				
				if(i_periodic) from	=	1;
				else if(ii_periodic) from	=	mx+1;
				
				level[k][j][to]	=	level[k][j][from];
			}
		}
	}
	
	if(ys==0 ||	ye==my)	{
		int	from,	to;
				
		for	(k=lzs;	k<lze; k++)
		for	(i=lxs;	i<lxe; i++)	{
			if(ys==0)	{
				j	=	1, from	=	j, to	=	0;
				
				if(j_periodic) from	=	my-2;
				else if(jj_periodic) from	=	-2;
				
				level[k][to][i]	=	level[k][from][i];
			}
			
			if(ye==my) {
				j	=	my-2,	from = j,	to = my-1;
				
				if(j_periodic) from	=	1;
				else if(jj_periodic) from	=	my+1;
				
				level[k][to][i]	=	level[k][from][i];
			}
		}
	}
	
	if(zs==0 ||	ze==mz)	{
		int	from,	to;
		
		for	(j=lys;	j<lye; j++)
		for	(i=lxs;	i<lxe; i++)	{
			if(zs==0)	{
				k	=	1, from	=	k, to	=	0;
				
				if(k_periodic) from = mz-2;
				else if(kk_periodic) from = -2;
				
				level[to][j][i]	= level[from][j][i];
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(k_periodic) from = 1;
				else if(kk_periodic) from = mz+1;
				
				level[to][j][i] = level[from][j][i];
			}
		}
	}
	
	DAVecRestoreArray(user->da, user->lLevelset, &level);
	
	Vec Rhs;
	VecDuplicate(user->P, &Rhs);
	Distance_Function_RHS(user, Rhs, 1);
	VecAXPY(LevelSet, dtau,	Rhs);
	VecDestroy(Rhs);
	return;
}

void Solve_Distance(UserCtx	*user, int iter)
{
	SNES snes_distance;
	KSP ksp;
	PC pc;
	Vec	r;
	Mat	J;
	double norm;
	
	int	bi=0;
		
	VecDuplicate(LevelSet, &r);
	
	
	SNESCreate(PETSC_COMM_WORLD,&snes_distance);
	SNESSetFunction(snes_distance,r,FormFunction_Distance,(void	*)&user[bi]);
	MatCreateSNESMF(snes_distance, &J);
	SNESSetJacobian(snes_distance,J,J,MatMFFDComputeJacobian,(void *)&user[bi]);
		
	
	SNESSetType(snes_distance, SNESTR);			//SNESTR,SNESLS	
	double tol=1.e-2;
	SNESSetMaxLinearSolveFailures(snes_distance,10000);
	SNESSetMaxNonlinearStepFailures(snes_distance,10000);		
	SNESKSPSetUseEW(snes_distance, PETSC_TRUE);
	SNESKSPSetParametersEW(snes_distance,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	SNESSetTolerances(snes_distance,PETSC_DEFAULT,tol,PETSC_DEFAULT,5,50000);	// snes	iter
		
	SNESGetKSP(snes_distance, &ksp);
	KSPSetType(ksp,	KSPGMRES);
	//KSPGMRESSetPreAllocateVectors(ksp);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCNONE);
	
	//int	maxits=10;	
	int	maxits=4;	// ksp iter
	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
	KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	
	extern PetscErrorCode	MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void	*dummy);
	SNESMonitorSet(snes_distance,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	SNESSolve(snes_distance, PETSC_NULL, LevelSet);
	
	SNESGetFunctionNorm(snes_distance, &norm);
	//PetscPrintf(PETSC_COMM_WORLD,	"\nDistance	SNES residual	norm=%.5e\n\n",	norm);
	VecDestroy(r);
	MatDestroy(J);
	SNESDestroy(snes_distance);
};

void Compute_Distance_Function(UserCtx *user)
{
	DALocalInfo	info;
	PetscInt	i, j,	k;
	DAGetLocalInfo(user->da, &info);
	PetscInt	mx = info.mx,	my = info.my,	mz = info.mz;
	PetscInt	xs = info.xs,	xe = xs	+	info.xm;
	PetscInt	ys = info.ys,	ye = ys	+	info.ym;
	PetscInt	zs = info.zs,	ze = zs	+	info.zm;
	PetscInt	lxs	=	xs,	lxe	=	xe;
	PetscInt	lys	=	ys,	lye	=	ye;
	PetscInt	lzs	=	zs,	lze	=	ze;
	PetscReal	***level;
	
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx)	lxe	=	xe-1;
	if (ye==my)	lye	=	ye-1;
	if (ze==mz)	lze	=	ze-1;
	
	int	iter=0;
	Init_LevelSet_Vectors(user);
	Distance_Function_IC(user);
		
	Vec	D;
	double norm=100, norm_old;
		
	dtau = dx_min*0.25;
	//dtau = dx_min	*	0.5;
	
	VecDuplicate(LevelSet, &D);
	
	do {
		norm_old = norm;
		
		PetscPrintf(PETSC_COMM_WORLD,	"\nSolving Distance %d... (dtau=%f)\n",	iter, dtau);
		//Solve_Distance_Explicit(user);
		Solve_Distance(user, iter);
		VecWAXPY(D, -1., LevelSet_o, LevelSet);
		VecNorm(D, NORM_INFINITY, &norm);
		PetscPrintf(PETSC_COMM_WORLD, "\nNorm=%.5e\n\n", norm);
		
		PetscReal ***level;
		DAVecGetArray(user->da, LevelSet, &level);
		for	(k=zs; k<ze; k++)
		for	(j=ys; j<ye; j++)
		for	(i=xs; i<xe; i++)	if(level[k][j][i]<0) level[k][j][i]=0.01;
		DAVecRestoreArray(user->da, LevelSet, &level);
	
		if(iter > 10 && (fabs(norm)<1.e-5 || iter>1000) ) break;
		//if(iter>10 && iter%50==0) dtau *=1.04;
		VecCopy(LevelSet, LevelSet_o);
		
		iter++;
	}	while(1);
	VecDestroy(D);
	
	Vec	lLevelset;
	PetscReal	***llevel;
	VecDuplicate(user->lP, &lLevelset);

	DAGlobalToLocalBegin(user->da, LevelSet, INSERT_VALUES,	lLevelset);
	DAGlobalToLocalEnd(user->da, LevelSet, INSERT_VALUES,	lLevelset);

	DAVecGetArray(user->da,	LevelSet,	&level);
	DAVecGetArray(user->da,	lLevelset, &llevel);

	for	(k=zs; k<ze; k++)
	for	(j=ys; j<ye; j++)
	for	(i=xs; i<xe; i++)	{
		int	c=k, b=j,	a=i, flag=0;
		
		// Neumann conditions
		if(i==0) a=1,	flag=1;
		if(i==mx-1)	a=mx-2,	flag=1;
		if(j==0) b=1,	flag=1;
		if(j==my-1)	b=my-2,	flag=1;
		if(k==0) c=1,	flag=1;
		if(k==mz-1)	c=mz-2,	flag=1;
		
		if(i_periodic	&& i==0) a=mx-2, flag=1;
		else if(i_periodic &&	i==mx-1) a=1,	flag=1;
		
		if(j_periodic	&& j==0) b=my-2, flag=1;
		else if(j_periodic &&	j==my-1) b=1,	flag=1;
		
		if(k_periodic	&& k==0) c=mz-2, flag=1;
		else if(k_periodic &&	k==mz-1) c=1,	flag=1;
		
		
		if(ii_periodic &&	i==0)	a=-2,	flag=1;
		else if(ii_periodic	&& i==mx-1)	a=mx+1,	flag=1;
		
		if(jj_periodic &&	j==0)	b=-2,	flag=1;
		else if(jj_periodic	&& j==my-1)	b=my+1,	flag=1;
		
		if(kk_periodic &&	k==0)	c=-2,	flag=1;
		else if(kk_periodic	&& k==mz-1)	c=mz+1,	flag=1;
		
		if(flag) level[k][j][i]	=	llevel[c][b][a];
	}
	DAVecRestoreArray(user->da,	LevelSet,	&level);
	DAVecRestoreArray(user->da,	lLevelset, &llevel);

	VecDestroy(lLevelset);
	VecCopy(LevelSet,	user->Distance);	
	Destroy_LevelSet_Vectors(user);
};

/*
void Compute_dlevel_center_levelset(int i, int j, int k, int mx, int my, int mz, double sgn, int wall_distance, PetscReal ***level, PetscReal ***nvert, double *dldc, double *dlde, double *dldz)
{
	double dplus0, dminus0,	dplusL,	dminusL, dplusR, dminusR;
	double _dplus, _dminus;
	double solid=0.1;

	// i-direction
	if ( i==1 || (i!=mx-2 && !wall_distance && nvert[k][j][i-1]>solid ) )	{
		dplus0=level[k][j][i+1]-level[k][j][i],	dminus0=level[k][j][i]-level[k][j][i-1];
		dplusL=level[k][j][i]-level[k][j][i-1],	dminusL=0;
		dplusR=level[k][j][i+2]-level[k][j][i+1], dminusR=level[k][j][i+1]-level[k][j][i];
		if(nvert[k][j][i+2]>1.1) dplusR=dplus0*0;
	}
	else if	( i==mx-2 || (i!=1 && !wall_distance && nvert[k][j][i+1]>solid ) ) {
		dplus0=level[k][j][i+1]-level[k][j][i],	dminus0=level[k][j][i]-level[k][j][i-1];
		dplusL=level[k][j][i]-level[k][j][i-1],	dminusL=level[k][j][i-1]-level[k][j][i-2];
		dplusR=0, dminusR=level[k][j][i+1]-level[k][j][i];
		if(nvert[k][j][i-2]>1.1) dminusL=dminus0*0;
	}
	else {
		dplus0=level[k][j][i+1]-level[k][j][i],	dminus0=level[k][j][i]-level[k][j][i-1];
		dplusL=level[k][j][i]-level[k][j][i-1],	dminusL=level[k][j][i-1]-level[k][j][i-2];
		dplusR=level[k][j][i+2]-level[k][j][i+1], dminusR=level[k][j][i+1]-level[k][j][i];
		if(nvert[k][j][i+2]>1.1) dplusR=dplus0*0;
		if(nvert[k][j][i-2]>1.1) dminusL=dminus0*0;
	}
	_dplus	= dplus0 -  0.5 * M( dplus0*dminus0,dplusR*dminusR );
	_dminus	= dminus0 + 0.5 * M( dplus0*dminus0, dplusL*dminusL );
	if(	sgn*dplus0<0 &&	sgn*(dminus0+dplus0)<0)	*dldc=_dplus;
	else if( sgn*dminus0>0 &&	sgn*(dplus0+dminus0)>0)	*dldc=_dminus;
	else if( sgn*dminus0<0 &&	sgn*dplus0>0)	*dldc=0;
	else *dldc=0.5*(_dplus+_dminus);

	// j-direction
	if ( j==1 || (j!=my-2 && !wall_distance && nvert[k][j-1][i]>solid) )	{
	        dplus0=level[k][j+1][i]-level[k][j][i], dminus0=level[k][j][i]-level[k][j-1][i];
	        dplusL=level[k][j][i]-level[k][j-1][i], dminusL=0;
	        dplusR=level[k][j+2][i]-level[k][j+1][i], dminusR=level[k][j+1][i]-level[k][j][i];
	        if(nvert[k][j+2][i]>1.1) dplusR=dplus0*0;
	}
	else if	( j==my-2 || (j!=1 && !wall_distance && nvert[k][j+1][i]>solid) ) {
		dplus0=level[k][j+1][i]-level[k][j][i],	dminus0=level[k][j][i]-level[k][j-1][i];
		dplusL=level[k][j][i]-level[k][j-1][i],	dminusL=level[k][j-1][i]-level[k][j-2][i];
		dplusR=0, dminusR=level[k][j+1][i]-level[k][j][i];
		if(nvert[k][j-2][i]>1.1) dminusL=dminus0*0;
	}
	else {
		dplus0=level[k][j+1][i]-level[k][j][i],	dminus0=level[k][j][i]-level[k][j-1][i];
		dplusL=level[k][j][i]-level[k][j-1][i],	dminusL=level[k][j-1][i]-level[k][j-2][i];
		dplusR=level[k][j+2][i]-level[k][j+1][i],	dminusR=level[k][j+1][i]-level[k][j][i];

		if(nvert[k][j+2][i]>1.1) dplusR=dplus0*0;
		if(nvert[k][j-2][i]>1.1) dminusL=dminus0*0;
	}

	_dplus = dplus0 - 0.5 * M( dplus0*dminus0, dplusR*dminusR );
	_dminus = dminus0 + 0.5 * M( dplus0*dminus0, dplusL*dminusL );

	if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dlde=_dplus;
	else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dlde=_dminus;
	else if( sgn*dminus0<0 && sgn*dplus0>0) *dlde=0;
	else *dlde=0.5*(_dplus+_dminus);
	
	// k-direction
	if ( k==1 || (k!=mz-2 && !wall_distance && nvert[k-1][j][i]>solid) ) {
		dplus0=level[k+1][j][i]-level[k][j][i],	dminus0=level[k][j][i]-level[k-1][j][i];
		dplusL=level[k][j][i]-level[k-1][j][i],	dminusL=0;
		dplusR=level[k+2][j][i]-level[k+1][j][i], dminusR=level[k+1][j][i]-level[k][j][i];
		if(nvert[k+2][j][i]>1.1) dplusR=dplus0*0;
	}
	else if	( k==mz-2 || (k!=1 && !wall_distance && nvert[k+1][j][i]>solid) ) {
		dplus0=level[k+1][j][i]-level[k][j][i],	dminus0=level[k][j][i]-level[k-1][j][i];
		dplusL=level[k][j][i]-level[k-1][j][i],	dminusL=level[k-1][j][i]-level[k-2][j][i];
		dplusR=0, dminusR=level[k+1][j][i]-level[k][j][i];
		if(nvert[k-2][j][i]>1.1) dminusL=dminus0*0;
	}
	else {
		dplus0=level[k+1][j][i]-level[k][j][i],	dminus0=level[k][j][i]-level[k-1][j][i];
		dplusL=level[k][j][i]-level[k-1][j][i],	dminusL=level[k-1][j][i]-level[k-2][j][i];
		dplusR=level[k+2][j][i]-level[k+1][j][i], dminusR=level[k+1][j][i]-level[k][j][i];

		if(nvert[k+2][j][i]>1.1) dplusR=dplus0*0;
		if(nvert[k-2][j][i]>1.1) dminusL=dminus0*0;
	}

	_dplus	= dplus0 -  0.5 * M( dplus0*dminus0, dplusR*dminusR );
	_dminus	= dminus0 + 0.5	* M( dplus0*dminus0, dplusL*dminusL );

	if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dldz=_dplus;
	else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dldz=_dminus;
	else if( sgn*dminus0<0 && sgn*dplus0>0)	*dldz=0;
	else *dldz=0.5*(_dplus+_dminus);
}
*/

void Compute_dlevel_center_levelset_100303 (int i, int j, int k,  int mx, int my, int mz, double sgn, int wall_distance, PetscReal ***level, PetscReal ***nvert, double *dldc, double *dlde, double *dldz)
{
        double dplus0, dminus0, dplusL, dminusL, dplusR, dminusR;
	double _dplus, _dminus;

	// i-direction                                                                                                                                                           
	if ( (i!=mx-2 && !wall_distance && (int)(nvert[k][j][i-1]+0.1)==1) || i==1) {
	  dplus0=level[k][j][i+1]-level[k][j][i], dminus0=level[k][j][i]-level[k][j][i-1];
	  dplusL=level[k][j][i]-level[k][j][i-1], dminusL=0;
	  dplusR=level[k][j][i+2]-level[k][j][i+1], dminusR=level[k][j][i+1]-level[k][j][i];
        }
        else if ( (i!=1 && !wall_distance && (int)(nvert[k][j][i+1]+0.1)==1) || i==mx-2) {
	  dplus0=level[k][j][i+1]-level[k][j][i], dminus0=level[k][j][i]-level[k][j][i-1];
	  dplusL=level[k][j][i]-level[k][j][i-1], dminusL=level[k][j][i-1]-level[k][j][i-2];
	  dplusR=0, dminusR=level[k][j][i+1]-level[k][j][i];
        }
	else {
                dplus0=level[k][j][i+1]-level[k][j][i], dminus0=level[k][j][i]-level[k][j][i-1];
                dplusL=level[k][j][i]-level[k][j][i-1], dminusL=level[k][j][i-1]-level[k][j][i-2];
                dplusR=level[k][j][i+2]-level[k][j][i+1], dminusR=level[k][j][i+1]-level[k][j][i];
        }
        _dplus    = dplus0    -  0.5 * M( dplus0*dminus0, dplusR*dminusR );
        _dminus = dminus0 + 0.5 * M( dplus0*dminus0, dplusL*dminusL );

        if( sgn*dplus0<0 && sgn*(dminus0+dplus0)<0) *dldc=_dplus;
        else if( sgn*dminus0>0 && sgn*(dplus0+dminus0)>0) *dldc=_dminus;
        else if( sgn*dminus0<0 && sgn*dplus0>0) *dldc=0;
        else *dldc=0.5*(_dplus+_dminus);

	// j-direction                                                                                                                                                           
        if ( (j!=my-2 && !wall_distance && (int)(nvert[k][j-1][i]+0.1)==1) || j==1) {
	  dplus0=level[k][j+1][i]-level[k][j][i], dminus0=level[k][j][i]-level[k][j-1][i];
	  dplusL=level[k][j][i]-level[k][j-1][i], dminusL=0;
	  dplusR=level[k][j+2][i]-level[k][j+1][i], dminusR=level[k][j+1][i]-level[k][j][i];
	  
        }
        else if ( (j!=1 && !wall_distance && (int)(nvert[k][j+1][i]+0.1)==1) || j==my-2) {
	  dplus0=level[k][j+1][i]-level[k][j][i], dminus0=level[k][j][i]-level[k][j-1][i];
	  dplusL=level[k][j][i]-level[k][j-1][i], dminusL=level[k][j-1][i]-level[k][j-2][i];
	  dplusR=0, dminusR=level[k][j+1][i]-level[k][j][i];
        }
        else {
                dplus0=level[k][j+1][i]-level[k][j][i], dminus0=level[k][j][i]-level[k][j-1][i];
                dplusL=level[k][j][i]-level[k][j-1][i], dminusL=level[k][j-1][i]-level[k][j-2][i];
                dplusR=level[k][j+2][i]-level[k][j+1][i], dminusR=level[k][j+1][i]-level[k][j][i];
        }

        _dplus    = dplus0    -  0.5 * M( dplus0*dminus0, dplusR*dminusR );
        _dminus = dminus0 + 0.5 * M( dplus0*dminus0, dplusL*dminusL );

        if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dlde=_dplus;
        else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dlde=_dminus;
        else if( sgn*dminus0<0 && sgn*dplus0>0) *dlde=0;
        else *dlde=0.5*(_dplus+_dminus);
	
	// k-direction
	if ( (k!=mz-2 && !wall_distance && (int)(nvert[k-1][j][i]+0.1)==1) || k==1) {
	  dplus0=level[k+1][j][i]-level[k][j][i], dminus0=level[k][j][i]-level[k-1][j][i];
	  dplusL=level[k][j][i]-level[k-1][j][i], dminusL=0;
	  dplusR=level[k+2][j][i]-level[k+1][j][i], dminusR=level[k+1][j][i]-level[k][j][i];
        }
        else if ( (k!=1 && !wall_distance && (int)(nvert[k+1][j][i]+0.1)==1) || k==mz-2) {
	  dplus0=level[k+1][j][i]-level[k][j][i], dminus0=level[k][j][i]-level[k-1][j][i];
	  dplusL=level[k][j][i]-level[k-1][j][i], dminusL=level[k-1][j][i]-level[k-2][j][i];
	  dplusR=0, dminusR=level[k+1][j][i]-level[k][j][i];
        }
        else {
                dplus0=level[k+1][j][i]-level[k][j][i], dminus0=level[k][j][i]-level[k-1][j][i];
                dplusL=level[k][j][i]-level[k-1][j][i], dminusL=level[k-1][j][i]-level[k-2][j][i];
                dplusR=level[k+2][j][i]-level[k+1][j][i], dminusR=level[k+1][j][i]-level[k][j][i];
        }

        _dplus    = dplus0    -  0.5 * M( dplus0*dminus0, dplusR*dminusR );
        _dminus = dminus0 + 0.5 * M( dplus0*dminus0, dplusL*dminusL );

        if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dldz=_dplus;
        else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dldz=_dminus;
        else if( sgn*dminus0<0 && sgn*dplus0>0) *dldz=0;
        else *dldz=0.5*(_dplus+_dminus);
}


void Compute_dlevel_center_levelset /*100202*/(int i, int j, int k,  int mx, int my, int mz, double sgn, int wall_distance, PetscReal ***level, PetscReal ***nvert, double *dldc, double *dlde, double *dldz)
{
	double dplus0, dminus0;
	double _dplus, _dminus;
	
	// i-direction
	int iLL=i-2, iL=i-1, iR=i+1, iRR=i+2;
	
	if ( /*!wall_distance &&*/ nvert[k][j][i-1]>0.1 ) {
	  //iL = i;
		iLL = iL;
	}
	if ( i!=1 && nvert[k][j][i-2]>1.1 ) {
		iLL = iL;
	}
	if ( /*!wall_distance &&*/ nvert[k][j][i+1]>0.1 ) {
	  //iR = i;
		iRR = iR;
	}
	if ( i!=mx-2 && nvert[k][j][i+2]>1.1 ) {
		iRR = iR;
	}
	if ( i==1 ) iLL = iL;
	if ( i==mx-2 ) iRR = iR;
	
	dplus0=level[k][j][iR]-level[k][j][i], dminus0=level[k][j][i]-level[k][j][iL];
	
	_dplus  =  dplus0 - 0.5 * M(level[k][j][iR]-2.*level[k][j][i]+level[k][j][iL], level[k][j][iRR]-2.*level[k][j][iR]+level[k][j][i]);
	_dminus = dminus0 + 0.5 * M(level[k][j][iR]-2.*level[k][j][i]+level[k][j][iL], level[k][j][i]-2.*level[k][j][iL]+level[k][j][iLL]);
		
	if( sgn*dplus0<0 && sgn*(dminus0+dplus0)<0) *dldc=_dplus;
	else if( sgn*dminus0>0 && sgn*(dplus0+dminus0)>0) *dldc=_dminus;
	else if( sgn*dminus0<0 && sgn*dplus0>0) *dldc=0;
	else *dldc=0.5*(_dplus+_dminus);
		
	// j-direction
	int jLL=j-2, jL=j-1, jR=j+1, jRR=j+2;
	
	if ( /*!wall_distance &&*/ nvert[k][j-1][i]>0.1 ) {
	  //jL = j;
		jLL = jL;
	}
	if ( j!=1 && nvert[k][j-2][i]>1.1 ) {
		jLL = jL;
	}
	if ( /*!wall_distance && */nvert[k][j+1][i]>0.1 ) {
	  //jR = j;
		jRR = jR;
	}
	if ( j!=my-2 && nvert[k][j+2][i]>1.1 ) {
		jRR = jR;
	}
	if ( j==1 ) jLL = jL;
        if ( j==my-2 ) jRR = jR;
	
	dplus0=level[k][jR][i]-level[k][j][i], dminus0=level[k][j][i]-level[k][jL][i];
	
	_dplus  =  dplus0 - 0.5 * M(level[k][jR][i]-2.*level[k][j][i]+level[k][jL][i], level[k][jRR][i]-2.*level[k][jR][i]+level[k][j][i]);
	_dminus = dminus0 + 0.5 * M(level[k][jR][i]-2.*level[k][j][i]+level[k][jL][i], level[k][j][i]-2.*level[k][jL][i]+level[k][jLL][i]);
			
	if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dlde=_dplus;
	else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dlde=_dminus;
	else if( sgn*dminus0<0 && sgn*dplus0>0) *dlde=0;
	else *dlde=0.5*(_dplus+_dminus);
	
	// k-direction
	int kLL=k-2, kL=k-1, kR=k+1, kRR=k+2;
	
	if ( /*!wall_distance &&*/ nvert[k-1][j][i]>0.1 ) {
	  //kL = k;
		kLL = kL;
	}
	if ( k!=1 && nvert[k-2][j][i]>1.1 ) {
		kLL = kL;
	}
	if ( /*!wall_distance &&*/ nvert[k+1][j][i]>0.1 ) {
	  //kR = k;
		kRR = kR;
	}
	if ( k!=mz-2 && nvert[k+2][j][i]>1.1 ) {
		kRR = kR;
	}
	if ( k==1 ) kLL = kL;
        if ( k==mz-2 ) kRR = kR;
	
	dplus0=level[kR][j][i]-level[k][j][i], dminus0=level[k][j][i]-level[kL][j][i];
	
	_dplus  =  dplus0 - 0.5 * M(level[kR][j][i]-2.*level[k][j][i]+level[kL][j][i], level[kRR][j][i]-2.*level[kR][j][i]+level[k][j][i]);
	_dminus = dminus0 + 0.5 * M(level[kR][j][i]-2.*level[k][j][i]+level[kL][j][i], level[k][j][i]-2.*level[kL][j][i]+level[kLL][j][i]);
		
	if( sgn*dplus0<0 && sgn*dminus0<-sgn*dplus0) *dldz=_dplus;
	else if( sgn*dminus0>0 && sgn*dplus0>-sgn*dminus0) *dldz=_dminus;
	else if( sgn*dminus0<0 && sgn*dplus0>0) *dldz=0;
	else *dldz=0.5*(_dplus+_dminus);
}
