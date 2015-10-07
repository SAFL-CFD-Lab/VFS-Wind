/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

// llevel for periodicity !!
#include "variables.h"

double dtau_levelset;
extern Vec LevelSet, LevelSet0, LevelSet_o;
extern double M(double a, double b);
extern int immersed, NumberOfBodies;
extern int i_periodic, j_periodic, k_periodic;
extern PetscInt   tiout;

void Init_Levelset_Vectors(UserCtx *user)
{
	VecDuplicate(user->P, &LevelSet);
	VecDuplicate(user->P, &LevelSet0);
	VecDuplicate(user->P, &LevelSet_o);
//	VecDuplicate(user->lP, &lLevelSet);
};

void Destroy_Levelset_Vectors(UserCtx *user)
{
	VecDestroy (LevelSet);
	VecDestroy (LevelSet0);
	VecDestroy (LevelSet_o);
//	VecDestroy (lLevelSet);
};

void Initialize_free_surface_location_vector(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	mx, my, mz;

	Cmpnts ***cent;
	DAVecGetArray(fda, user->lCent, &cent);	
	DAGetLocalInfo(da, &info);
	mx = info.mx; my = info.my; mz = info.mz;
	
	int xs = info.xs, xe = xs + info.xm;
	int ys = info.ys, ye = ys + info.ym;
	int zs = info.zs, ze = zs + info.zm;

	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	//temporarily allocate
	double *buffer_x, *buffer_z;
	buffer_x = (double *) malloc( sizeof(double *) * (mx-2));
	buffer_z = (double *) malloc( sizeof(double *) * (mz-2));
	for(int i=0; i<mx-2; i++) {
		buffer_x[i]=-10000000;
	}
	for(int k=0; k<mz-2; k++) {
		buffer_z[k]=-10000000;
	}
	

	user->free_surface_location = (double *) malloc( sizeof(double *) * (mx-2) * (mz-2));
	user->free_surface_x = (double *) malloc( sizeof(double *) * (mx-2));
	user->free_surface_z = (double *) malloc( sizeof(double *) * (mz-2));
	for(int i=0; i<(mx-2); i++) 
	for(int k=0; k<(mz-2); k++) {
		user->free_surface_location[i*(mz-2)+k]=0;
	}
	for(int i=0; i<(mx-2); i++) {
		user->free_surface_x[i]=0.;
	}
	for(int k=0; k<(mz-2); k++) {
		user->free_surface_z[k]=0.;
	}
	for (int i=lxs; i<lxe; i++) {
		int j=1,k=1;
		if(j<=lye && j>=lys && k<=lze && k>=lzs){
			buffer_x[i-1]=cent[k][j][i].x;//for cartesian meshes only. i-1 in buffer_x is i in cent
		}
	}
	for (int k=lzs; k<lze; k++) {
		int j=1,i=1;
		if(j<=lye && j>=lys && i<=lxe && i>=lxs){
			buffer_z[k-1]=cent[k][j][i].z;//for cartesian meshes only. j-1 in buffer_x is j in cent
		}
	}
	MPI_Allreduce ( &buffer_x[0], &user->free_surface_x[0], mx-2, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
	MPI_Allreduce ( &buffer_z[0], &user->free_surface_z[0], mz-2, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
	// free memeory
	free ( buffer_x );
	free ( buffer_z );
	DAVecRestoreArray(fda, user->lCent, &cent);
};

void Calc_free_surface_location(UserCtx *user)
{
	DA da = user->da, fda = user->fda;
	DALocalInfo info;
	PetscInt mx, my, mz;
	PetscReal ***level, ***nvert;
	Cmpnts ***cent;

	DAGetLocalInfo(da, &info);
	mx = info.mx; my = info.my; mz = info.mz;

	int xs = info.xs, xe = xs + info.xm;
	int ys = info.ys, ye = ys + info.ym;
	int zs = info.zs, ze = zs + info.zm;
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	double *buffer;

	//i,j,k=0 in cent corresponds to ghost
	//i,j,k=0 in buffer and free_surface_location it is the first node
	//so cent[k][j][i].x    ====   free_surface_x[(i-1)]
	
	//temporarily allocate
	buffer = (double *) malloc( sizeof(double *) * (mx-2) * (mz-2));
	for(int i=0; i<(mx-2); i++)
	for(int k=0; k<(mz-2); k++) {
		buffer[i*(mz-2)+k]=-10000;
	}

	DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, user->lCent, &cent);

	for (int k=lzs; k<lze; k++)
	for (int j=lys; j<lye; j++)
	for (int i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>0.1) continue;

		if( level[k][j][i]>=0 && level[k][j+1][i]<0 ) {	// water surface is above my cell center
			buffer[(i-1)*(mz-2)+(k-1)] = cent[k][j][i].y + level[k][j][i];
		}
		else if( level[k][j][i]<0 && level[k][j-1][i]>=0 ) {	// water surface is below my cell center
			buffer[(i-1)*(mz-2)+(k-1)] = cent[k][j][i].y + level[k][j][i];
		}
	}
	DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, user->lCent, &cent);

	MPI_Reduce ( &buffer[0], &user->free_surface_location[0], (mx-2)*(mz-2), MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);

	for(int i=0; i<mx-2; i++)
	for(int k=0; k<mz-2; k++) {
		if( (int) (user->free_surface_location[i*(mz-2)+k])==-10000 ) {
			user->free_surface_location[i*(mz-2)+k] = 0;
		}
	}
	if (ti == (ti/tiout)*tiout) {
		int rank=0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		if (!rank) {
			FILE *f;
			char filen[80];
			sprintf(filen, "FreeSurfaceElev_%06d.dat",ti);
			f = fopen(filen, "w");
			PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"X\", \"Y\", \"Z\"\n");
			PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE I=%d J=%d K=%d F=POINT \n", mx-2,1,mz-2);
			for(int k=0; k<mz-2; k++)
			for(int i=0; i<mx-2; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le\n",user->free_surface_x[i],user->free_surface_location[i*(mz-2)+k],user->free_surface_z[k]);
			}
			fclose(f);
		}
	}

	// free memeory
	free ( buffer );
	PetscPrintf(PETSC_COMM_WORLD,  "free_surface_location computed\n");
};


double dist(Cmpnts &a, Cmpnts &b)
{
  return sqrt( pow(a.x-b.x,2.) + pow(a.y-b.y,2.) + pow(a.z-b.z,2.) );
};


void Levelset_BC(UserCtx *user)
{
	DA		da = user->da;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal	***nvert, ***level, ***llevel, ***aj;
	Cmpnts ***cent;

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
	/*		
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
        DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	*/
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->da, user->Levelset, &level);
	DAVecGetArray(user->da, user->lLevelset, &llevel);
	DAVecGetArray(user->fda, user->lCent, &cent);
	DAVecGetArray(user->da, user->lAj, &aj);

	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if((int)(nvert[k][j][i]+0.1)==1
		   ) {
			int count=0;
			double sum=0;
			
			if (i<mx-3 && nvert[k][j][i+1]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k][j][i+1]);
			  double AC = dist(cent[k][j][i], cent[k][j][i+2]);
			  double lB = llevel[k][j][i+1];
			  double lC = llevel[k][j][i+2];
			  double val = lC - AC/AB * (lC - lB);

			  //if(nvert[k][j][i+2]<0.1) sum += val, count++;
			  sum += llevel[k][j][i+1], count++;
			}
			if (i>2 && nvert[k][j][i-1]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k][j][i-1]);
                          double AC = dist(cent[k][j][i], cent[k][j][i-2]);
                          double lB = level[k][j][i-1];
                          double lC = level[k][j][i-2];
                          double val = lC - AC/AB * (lC - lB);
			  //if(nvert[k][j][i-2]<0.1) sum += val, count++;
			  sum += llevel[k][j][i-1], count++;
			}
			/*
			if (j<my-3 && nvert[k][j+1][i]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k][j+1][i]);
                          double AC = dist(cent[k][j][i], cent[k][j+2][i]);
                          double lB = llevel[k][j+1][i];
                          double lC = llevel[k][j+2][i];
                          double val = lC - AC/AB * (lC - lB);
			  if(nvert[k][j+2][i]<0.1) sum += val, count++;
			}
			if (j>2 && nvert[k][j-1][i]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k][j-1][i]);
                          double AC = dist(cent[k][j][i], cent[k][j-2][i]);
                          double lB = llevel[k][j-1][i];
                          double lC = llevel[k][j-2][i];
                          double val = lC - AC/AB * (lC - lB);
			  if(nvert[k][j-2][i]<0.1) sum += val, count++;
			}
			*/
			if (k<mz-3 && nvert[k+1][j][i]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k+1][j][i]);
                          double AC = dist(cent[k][j][i], cent[k+2][j][i]);
                          double lB = llevel[k+1][j][i];
                          double lC = llevel[k+2][j][i];
                          double val = lC - AC/AB * (lC - lB);
			  // if(nvert[k+2][j][i]<0.1) sum += val, count++;
			  sum += llevel[k+1][j][i], count++;
			}
			if (k>2 && nvert[k-1][j][i]<0.1) {
			  double AB = dist(cent[k][j][i], cent[k-1][j][i]);
                          double AC = dist(cent[k][j][i], cent[k-2][j][i]);
                          double lB = llevel[k-1][j][i];
                          double lC = llevel[k-2][j][i];
                          double val = lC - AC/AB * (lC - lB);
			  //if(nvert[k-2][j][i]<0.1) sum += val, count++;
			  sum += llevel[k-1][j][i], count++;
			}
			
			if(count) {	// prevent NaN
				level[k][j][i] = sum / (double)count;
			}
		}
	}
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				level[k][j][to] = llevel[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				level[k][j][to] = llevel[k][j][from];
			}
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
				
				level[k][to][i] = llevel[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = llevel[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				
				if(user->bctype[4]==5 && fix_inlet) {
				  double y=cent[k][j][i].y, z=cent[k][j][i].z;
				  if(inlet_y_flag) {
				    level[k][j][i] = (inlet_y - y);
				    level[to][j][i] = (inlet_y - y);
				  }
				  else if(inlet_z_flag) {
				    level[k][j][i] = (inlet_z - z);
				    level[to][j][i] = (inlet_z - z);
				  }
				}
				else {
				  if(k_periodic) from = mz-2;
				  else if(kk_periodic) from = -2;
				  level[to][j][i] = llevel[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(user->bctype[5]==4 && fix_outlet) {
				  double y=cent[k][j][i].y, z=cent[k][j][i].z;
				  if(inlet_y_flag) {
				    level[k][j][i] = (outlet_y - y);
				    level[to][j][i] = (outlet_y - y);
				  }
				  else if(inlet_z_flag) {
				    level[k][j][i] = (outlet_z - z);
				    level[to][j][i] = (outlet_z - z);
				  }
				}
				else {
				  if(k_periodic) from = 1;
				  else if(kk_periodic) from = mz+1;
				  level[to][j][i] = llevel[from][j][i];
				}
			}
		}
	}
	
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, user->Levelset, &level);
	DAVecRestoreArray(user->da, user->lLevelset, &llevel);
	DAVecRestoreArray(user->fda, user->lCent, &cent);
	DAVecRestoreArray(user->da, user->lAj, &aj);
	
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
}

void Levelset_Function_IC(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts	***csi, ***eta, ***zet, ***cent;
	PetscReal	***nvert, ***level,***llevel,  ***aj, ***p;

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
	DAVecGetArray(da, user->Levelset, &level);
	DAVecGetArray(da, user->lLevelset, &llevel);
		
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(fda, user->lCent, &cent);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->P, &p);
		
	std::vector<double> elevation_k ( mz );	// for inlet-outlet interpolation
	if( (user->bctype[4]==5 && user->bctype[5]==4) ) {
		if(inlet_y_flag) {
			elevation_k[0] = inlet_y;
			elevation_k[mz-1] = outlet_y;
			
			for (k=1; k<=mz-2; k++) {
				elevation_k[k] = ( outlet_y - inlet_y ) / (double) (mz -2 - 1) * (double) ( k - 1 ) + inlet_y;
			}
		}
		else if(inlet_z_flag) {
			elevation_k[0] = inlet_z;
			elevation_k[mz-1] = outlet_z;
			
			for (k=1; k<=mz-2; k++) {
				elevation_k[k] = ( outlet_z - inlet_z ) / (double) (mz -2 - 1) * (double) ( k - 1 ) + inlet_z;
			}
		}
		
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double x0, y0, z0, R=0.15;
		double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
		
		if( (user->bctype[4]==5 && user->bctype[5]==4) ) {
			if(inlet_y_flag) level[k][j][i] = elevation_k[k] - y;
			else if(inlet_z_flag) level[k][j][i] = elevation_k[k] - z;
			p[k][j][i] = 0;
			if(level[k][j][i]>0) {
			  if(inlet_y_flag)  p[k][j][i] = - rho_water * gravity_y * level[k][j][i];
			  else if(inlet_z_flag) p[k][j][i] = - rho_water * gravity_z * level[k][j][i]; 
			}
			continue;
		}
		else{ 
			if(level_in==1) {
				level[k][j][i] = level_in_height - z;
				p[k][j][i] = 0;
				if(level[k][j][i]>0) {
					p[k][j][i] = - rho_water * gravity_z * level[k][j][i];
				}			
			}
			if(level_in==2) {
				level[k][j][i] = level_in_height - y;
				p[k][j][i] = 0;
				if(level[k][j][i]>0) {
					p[k][j][i] = - rho_water * gravity_y * level[k][j][i];
				}			
			}				
			if(sloshing==1) // 1d sloshing
			{
				double a=sloshing_a, b=sloshing_b, d=sloshing_d;	// a=0.05 (nonlinear), a=0.001 (linear)
				double k2 = 2*M_PI/b;
				double xi = d + a * cos ( 2.0 * M_PI * x / b );
				if(!inviscid) xi = d + a * cos ( k2*x );
				level[k][j][i] = xi - z;
			}
			if(sloshing==2) { // 2d sloshing
				double L = 20., d=1., Beta=0.25, a=0.1; // L = width of the 3D tank
				double eta0 = a * exp ( -Beta * ( pow(x-L/2, 2) + pow(y-L/2, 2) ) );
				level[k][j][i] = eta0 + d - z;
				p[k][j][i] = 0;
				if(level[k][j][i]>0) {
					p[k][j][i] = - rho_water * gravity_z * level[k][j][i];
				}	
			}		
		}
		if( nvert[k][j][i]>1.1 || i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) level[k][j][i] = 1.*sign(level[k][j][i]);
	}
		
	// Neumann , periodic conditions
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				level[k][j][to] = llevel[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				level[k][j][to] = llevel[k][j][from];
			}
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
				
				level[k][to][i] = llevel[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = llevel[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
				
				if(user->bctype[4]==5) {
					if(inlet_y_flag) level[to][j][i] = 2.*(inlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(inlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = mz-2;
					else if(kk_periodic) from = -2;
					level[to][j][i] = llevel[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;

				if(user->bctype[5]==4 && fix_outlet) {
					/*if(inlet_y_flag) level[k][j][i] = outlet_y - y;
					else if(inlet_z_flag) level[k][j][i] = outlet_z - z;
					*/
					if(inlet_y_flag) level[to][j][i] = 2.*(outlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(outlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = 1;
					else if(kk_periodic) from = mz+1;
					level[to][j][i] = llevel[from][j][i];
				}
			}
		}
	}
	
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->Levelset, &level);
	DAVecRestoreArray(da, user->lLevelset, &llevel);
		
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(fda, user->lCent, &cent);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->P, &p);
	
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
        DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);

	DAGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);
};

PetscErrorCode FormFunction_Levelset (SNES snes, Vec L, Vec Rhs, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal ***level;
	Cmpnts ***cent;

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
	
	DAVecGetArray(user->fda, user->lCent, &cent);
	
	
	DAGlobalToLocalBegin(user->da, L, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, L, INSERT_VALUES, user->lLevelset);
	
	DAVecGetArray(user->da, user->lLevelset, &level);
		
		
	// Neumann , periodic conditions
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				level[k][j][to] = level[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				level[k][j][to] = level[k][j][from];
			}
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
				
				level[k][to][i] = level[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = level[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
				
				if(user->bctype[4]==5) {
					if(inlet_y_flag) level[to][j][i] = 2.*(inlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(inlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = mz-2;
					else if(kk_periodic) from = -2;
					level[to][j][i] = level[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;

				if(user->bctype[5]==4 && fix_outlet) {
					/*if(inlet_y_flag) level[k][j][i] = outlet_y - y;
					else if(inlet_z_flag) level[k][j][i] = outlet_z - z;
					*/
					if(inlet_y_flag) level[to][j][i] = 2.*(outlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(outlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = 1;
					else if(kk_periodic) from = mz+1;
					level[to][j][i] = level[from][j][i];
				}
			}
		}
	}
	
	DAVecRestoreArray(user->fda, user->lCent, &cent);
	DAVecRestoreArray(user->da, user->lLevelset, &level);
	
	Distance_Function_RHS(user, Rhs, 0);
	VecAXPY(Rhs, -1./dtau_levelset, L);
	VecAXPY(Rhs, 1./dtau_levelset, LevelSet_o);
	
	return(0);
}

void Solve_Reinit_explicit(UserCtx *user, int iter)
{
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal ***level;
	Cmpnts ***cent;

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
	
	DAVecGetArray(user->fda, user->lCent, &cent);
	
	DAGlobalToLocalBegin(user->da, LevelSet, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, LevelSet, INSERT_VALUES, user->lLevelset);
	
	DAVecGetArray(user->da, user->lLevelset, &level);
		
		
	// Neumann , periodic conditions
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				level[k][j][to] = level[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				level[k][j][to] = level[k][j][from];
			}
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
				
				level[k][to][i] = level[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = level[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
				
				if(user->bctype[4]==5) {
					if(inlet_y_flag) level[to][j][i] = 2.*(inlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(inlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = mz-2;
					else if(kk_periodic) from = -2;
					level[to][j][i] = level[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;

				if(user->bctype[5]==4 && fix_outlet) {
					if(inlet_y_flag) level[to][j][i] = 2.*(outlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(outlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = 1;
					else if(kk_periodic) from = mz+1;
					level[to][j][i] = level[from][j][i];
				}
			}
		}
	}
	
	DAVecRestoreArray(user->fda, user->lCent, &cent);
	DAVecRestoreArray(user->da, user->lLevelset, &level);
	
	Vec Rhs;
	
	VecDuplicate(user->P, &Rhs);
		
	Distance_Function_RHS(user, Rhs, 0);
	
	VecAXPY(LevelSet, dtau_levelset, Rhs);
	VecDestroy(Rhs);
	//VecAXPY(Rhs, -1./dtau_levelset, L);
	//VecAXPY(Rhs, 1./dtau_levelset, LevelSet_o);
	return;
}

void Solve_Reinit_implicit(UserCtx *user, int iter)
{
	SNES snes_distance;
	KSP ksp;
	PC pc;
	Vec r;
	Mat J;
	double norm;
	
	int bi=0;
	VecDuplicate(LevelSet, &r);
	SNESCreate(PETSC_COMM_WORLD,&snes_distance);
	SNESSetFunction(snes_distance,r,FormFunction_Levelset,(void *)&user[bi]);
	MatCreateSNESMF(snes_distance, &J);
	SNESSetJacobian(snes_distance,J,J,MatMFFDComputeJacobian,(void *)&user[bi]);
		
	
	SNESSetType(snes_distance, SNESTR);			//SNESTR,SNESLS	
	double tol=1.e-2;
	SNESSetMaxLinearSolveFailures(snes_distance,10000);
	SNESSetMaxNonlinearStepFailures(snes_distance,10000);		
	SNESKSPSetUseEW(snes_distance, PETSC_TRUE);
	SNESKSPSetParametersEW(snes_distance,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	SNESSetTolerances(snes_distance,PETSC_DEFAULT,tol,PETSC_DEFAULT,5,50000);
		
	SNESGetKSP(snes_distance, &ksp);
	KSPSetType(ksp, KSPGMRES);
	//KSPGMRESSetPreAllocateVectors(ksp);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCNONE);
	
	int maxits=4;	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
	KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	
	extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy);
	SNESMonitorSet(snes_distance,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	
	PetscPrintf(PETSC_COMM_WORLD, "\nSolving Levelset %d...\n", iter);
	SNESSolve(snes_distance, PETSC_NULL, LevelSet);
	VecCopy(LevelSet, user->Levelset);
	
	SNESGetFunctionNorm(snes_distance, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "\nDistance SNES residual norm=%.5e\n\n", norm);
	VecDestroy(r);
	MatDestroy(J);
	SNESDestroy(snes_distance);
};

void Reinit_Levelset(UserCtx *user)
{
	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	
	int iter=0, maxit=levelset_it;
	
	//if(ti<tistart+3) maxit=15;
	
	Init_Levelset_Vectors(user);
	
	Vec D;
	
	VecCopy (user->Levelset, LevelSet0);
	VecCopy (LevelSet0, LevelSet); 
	VecCopy (LevelSet0, LevelSet_o); 
	
	double norm, norm_old;
	
	//dtau_levelset = dx_min * 0.10;//0.15
	//dtau_levelset = dx_min * 0.05;
	dtau_levelset = dx_min * levelset_tau;

	VecDuplicate(LevelSet, &D);
	do {
		norm_old = norm;
		//Solve_Reinit_implicit(user, iter);
		Solve_Reinit_explicit(user, iter);
		VecWAXPY(D, -1., LevelSet_o, LevelSet);
		VecNorm(D, NORM_INFINITY, &norm);
			
		if(/*fabs(norm)<1.e-5 ||*/ iter++>= maxit ) break;
		VecCopy(LevelSet, LevelSet_o);
		
		//if(iter>500 && iter%500==0) dtau_levelset *= 2;
	} while(1);
	VecDestroy(D);
		
	DALocalInfo	info;
	PetscInt	i, j, k;
	DAGetLocalInfo(user->da, &info);
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	xs = info.xs, xe = xs + info.xm;
	PetscInt	ys = info.ys, ye = ys + info.ym;
	PetscInt	zs = info.zs, ze = zs + info.zm;
	PetscInt	lxs = xs, lxe = xe;
	PetscInt	lys = ys, lye = ye;
	PetscInt	lzs = zs, lze = ze;
	PetscReal ***level, ***llevel, ***nvert;
	Cmpnts ***cent;
	
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	VecCopy(LevelSet, user->Levelset);
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->da, user->Levelset, &level);
	DAVecGetArray(user->da, user->lLevelset, &llevel);
	DAVecGetArray(user->fda, user->lCent, &cent);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>0.1) {
			if ( nvert[k][j][i+1]+nvert[k][j][i-1]+nvert[k+1][j][i]+nvert[k-1][j][i] > 3.9 ) continue;
			double count=0;
			double sum=0;
			
			if (nvert[k][j][i+1]>0.1) sum += llevel[k][j][i+1], count+=1.;
			if (nvert[k][j][i- 1]>0.1) sum += llevel[k][j][i -1], count+=1.;
			if (nvert[k+1][j][i]>0.1) sum += llevel[k+1][j][i], count+=1.;
			if (nvert[k- 1][j][i]>0.1) sum += llevel[k -1][j][i], count+=1.;
			
			level[k][j][i] = sum / count;
		}
	}
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				level[k][j][to] = llevel[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				level[k][j][to] = llevel[k][j][from];
			}
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
				
				level[k][to][i] = llevel[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = llevel[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				
				if(user->bctype[4]==5) {
				  double y=cent[k][j][i].y, z=cent[k][j][i].z;
				  if(inlet_y_flag) level[to][j][i] = 2.*(inlet_y - y) - llevel[k][j][i];
				  else if(inlet_z_flag) level[to][j][i] = 2.*(inlet_z - z) - llevel[k][j][i];
                                }
				else {
				  if(k_periodic) from = mz-2;
				  else if(kk_periodic) from = -2;
				  level[to][j][i] = llevel[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(user->bctype[5]==4 && fix_outlet) {
				  double y=cent[k][j][i].y, z=cent[k][j][i].z;
				  if(inlet_y_flag) level[to][j][i] = 2.*(outlet_y - y) - llevel[k][j][i];
				  else if(inlet_z_flag) level[to][j][i] = 2.*(outlet_z - z) - llevel[k][j][i];
                                }
				else {
				  if(k_periodic) from = 1;
				  else if(kk_periodic) from = mz+1;
				  level[to][j][i] = llevel[from][j][i];
				}
			}
		}
	}
	
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, user->Levelset, &level);
	DAVecRestoreArray(user->da, user->lLevelset, &llevel);
	DAVecRestoreArray(user->fda, user->lCent, &cent);

	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
        DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	
	Destroy_Levelset_Vectors(user);
	
	PetscGetTime(&te);
	cput=te-ts;
	if (!my_rank) {
		FILE *f;
		char filen[80];
  		sprintf(filen, "%s/Converge_dU", path);
		f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d(levelset) %.2e(s) %le\n", ti, cput, norm);
		fclose(f);
	}
};

void Compute_Density(UserCtx *user)
{
	DA		da = user->da;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal	***nvert, ***level, ***rho, ***mu, ***aj;

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
	DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, user->lDensity, &rho);
	DAVecGetArray(da, user->lMu, &mu);
	DAVecGetArray(da, user->lAj, &aj);

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		
		if(ti==tistart) {
			if( 	levelset && user->bctype[4]==5 && k==mz-2 && 
				(level[k][j][i]<0 || level[k][j][i]*level[k][j-1][i]<0 || level[k][j][i]*level[k][j+1][i]<0 )
			) {
				//nvert[k][j][i] = 1.0;
			}
		}
		
		double dx=pow(1./aj[k][j][i],1./3.);
		if(dthick_set) dx = dthick;

		rho[k][j][i] = rho_air + (rho_water - rho_air) * H ( level[k][j][i], dx );
		mu[k][j][i] = mu_air + (mu_water - mu_air) * H ( level[k][j][i], dx );
		
		if(rho[k][j][i]<0) printf("Negative density !\n");
		if(mu[k][j][i]<0) printf("Negative viscosity!\n");

		if(i==1) rho[k][j][i-1]=rho[k][j][i];
		if(i==mx-2) rho[k][j][i+1]=rho[k][j][i];
		if(j==1) rho[k][j-1][i]=rho[k][j][i];
		if(j==my-2) rho[k][j+1][i]=rho[k][j][i];
		if(k==1) rho[k-1][j][i]=rho[k][j][i];
		if(k==mz-2) rho[k+1][j][i]=rho[k][j][i];
		
		
		if(i==1) mu[k][j][i-1]=mu[k][j][i];
		if(i==mx-2) mu[k][j][i+1]=mu[k][j][i];
		if(j==1) mu[k][j-1][i]=mu[k][j][i];
		if(j==my-2) mu[k][j+1][i]=mu[k][j][i];
		if(k==1) mu[k-1][j][i]=mu[k][j][i];
		if(k==mz-2) mu[k+1][j][i]=mu[k][j][i];
		
		if( nvert[k][j][i]>0.1) rho[k][j][i] = 0;
		if( nvert[k][j][i]>0.1) mu[k][j][i] = 0;
	}
		
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, user->lDensity, &rho);
	DAVecRestoreArray(da, user->lMu, &mu);
	DAVecRestoreArray(da, user->lAj, &aj);
	
	DALocalToLocalBegin(da, user->lDensity, INSERT_VALUES, user->lDensity);
	DALocalToLocalEnd(da, user->lDensity, INSERT_VALUES, user->lDensity);
	
	DALocalToLocalBegin(da, user->lMu, INSERT_VALUES, user->lMu);
	DALocalToLocalEnd(da, user->lMu, INSERT_VALUES, user->lMu);
	
	if(ti==tistart) {
		DALocalToLocalBegin(da, user->lNvert, INSERT_VALUES, user->lNvert);
		DALocalToLocalEnd(da, user->lNvert, INSERT_VALUES, user->lNvert);
	}
	
	DAVecGetArray(da, user->lDensity, &rho);
	DAVecGetArray(da, user->lMu, &mu);
	
	// Neumann , periodic conditions
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				rho[k][j][to] = rho[k][j][from];
				mu[k][j][to] = mu[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				rho[k][j][to] = rho[k][j][from];
				mu[k][j][to] = mu[k][j][from];
			}
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
				
				rho[k][to][i] = rho[k][from][i];
				mu[k][to][i] = mu[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				rho[k][to][i] = rho[k][from][i];
				mu[k][to][i] = mu[k][from][i];
			}
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
				
				rho[to][j][i] = rho[from][j][i];
				mu[to][j][i] = mu[from][j][i];
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(k_periodic) from = 1;
				else if(kk_periodic) from = mz+1;
				
				rho[to][j][i] = rho[from][j][i];
				mu[to][j][i] = mu[from][j][i];
			}
		}
	}

	
	DAVecRestoreArray(da, user->lDensity, &rho);
	DAVecRestoreArray(da, user->lMu, &mu);
};

double minmod(double m1, double m2)	// for UNO
{
	return 0.5 * ( sign(m1)+sign(m2) ) * std::min ( fabs(m1), fabs(m2) );
}

double eno2(double f0, double f1, double f2, double f3, double a)
{
	if( a > 0 ) return ( f1 + 0.5*M(f2-f1, f1-f0) );
	else  return ( f2 - 0.5*M(f3-f2, f2-f1) );
};


double weno3(double f0, double f1, double f2, double f3, double wavespeed)
{
	double fL, fC, fR;
	
	if(wavespeed>0)  {
		fL = f0; 
		fC = f1; 
		fR = f2;
	}
	else {
		// mirror
		fL = f3; 
		fC = f2; 
		fR = f1; 
	}
	
	double d0=2./3., d1=1./3.;	// weno3
	//if(wavespeed<=0) d0=1./3., d1=2./3.;
	
	const double eps=1.e-6;
		
	double beta0 = pow( fC - fR,  2. );
	double beta1 = pow( fL - fC,  2. );
	
	double alpha0 = d0 / pow( eps + beta0, 2. );
	double alpha1 = d1 / pow( eps + beta1, 2. );

	double sumalpha = alpha0 + alpha1;
	
	double w0 = alpha0 / sumalpha;
	double w1 = alpha1 / sumalpha;
	
	double u0 = ( fC*0.5  + fR*0.5 );
	double u1 = ( -fL*0.5 + fC*1.5 );
	
	return w0*u0 + w1*u1;
};

/*
High Order Finite Difference and Finite Volume
WENO Schemes and Discontinuous Galerkin Methods
for CFD, by Shu, ICASE
*/
double weno5(double f0, double f1, double f2, double f3, double f4, double f5, double wavespeed)
{
	double A, B, C, D, E;
	
	if(wavespeed>0) {
		A = f0;
		B = f1;
		C = f2;
		D = f3;
		E = f4;
	}
	else {
		// mirror
		A = f5;
		B = f4;
		C = f3;
		D = f2;
		E = f1;
	}
	
	double eps = 1.e-6;// * std::max ( A*A, std::max ( B*B, std::max ( C*C, std::max ( D*D, E*E ) ) ) ) + 1.e-99;
	double d0, d1, d2;
	
	//if(wavespeed<=0) d0=0.1, d1=0.6, d2=0.3;
	//else 
	  d0=0.3, d1=0.6, d2=0.1;
	
	double beta0 = 13./12. * pow( A - 2. * B + C, 2. ) + 1./4. * pow ( A - 4. * B  + 3. * C, 2. );
	double beta1 = 13./12. * pow( B - 2. * C + D, 2. ) + 1./4. * pow ( B - D, 2. );
	double beta2 = 13./12. * pow( C - 2. * D + E, 2. ) + 1./4. * pow ( 3. * C - 4. * D  + E, 2. );
	
	double alpha0 = d0 / pow( eps + beta0, 2.);
	double alpha1 = d1 / pow( eps + beta1, 2.);
	double alpha2 = d2 / pow( eps + beta2, 2.);
	
	double sumalpha = alpha0 + alpha1 + alpha2;
		
	double w0 = alpha0 / sumalpha;
	double w1 = alpha1 / sumalpha;
	double w2 = alpha2 / sumalpha;
	
	double u0 = 2./6. * A - 7./6. * B + 11./6. * C;
	double u1 = -1./6. * B + 5./6. * C + 2./6. * D;
	double u2 = 2./6. * C + 5./6. * E - 1./6. * E;
	
	return w0*u0 + w1*u1 + w2*u2;
	
};

void Levelset_Advect_RHS(UserCtx *user, Vec DRHS)
{
	DA 		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs, xe, ys, ye, zs, ze;
	PetscInt	mx, my, mz;
	PetscInt	i, j, k;
	
	Vec	Aj  = user->lAj;

	Cmpnts	***ucont, ***kzet, ***cent;
	PetscReal	***aj;
	PetscReal	***level, ***rhs, ***nvert;
	PetscInt	lxs, lys, lzs, lxe, lye, lze;

	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;
  	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  	mx = info.mx; my = info.my; mz = info.mz;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	VecSet(DRHS,0);
	
	DAVecGetArray(fda, user->lCent, &cent);
	
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
		
	DAVecGetArray(user->da, user->lLevelset, &level);
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				level[k][j][to] = level[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				level[k][j][to] = level[k][j][from];
			}
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
				
				level[k][to][i] = level[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				level[k][to][i] = level[k][from][i];
			}
		}
	}
	
	if(zs==0 || ze==mz) {
		int from, to;
		
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if(zs==0) {
				k = 1, from = k, to = 0;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
				
				if(user->bctype[4]==5) {
					if(inlet_y_flag) level[to][j][i] = 2.*(inlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(inlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = mz-2;
					else if(kk_periodic) from = -2;
					level[to][j][i] = level[from][j][i];
				}
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				double x = cent[k][j][i].x, y = cent[k][j][i].y, z = cent[k][j][i].z;
				
				if(user->bctype[5]==4 && fix_outlet) {
					if(inlet_y_flag) level[to][j][i] = 2.*(outlet_y - y) - level[k][j][i];
					else if(inlet_z_flag) level[to][j][i] = 2.*(outlet_z - z) - level[k][j][i];
				}
				else {
					if(k_periodic) from = 1;
					else if(kk_periodic) from = mz+1;
					level[to][j][i] = level[from][j][i];
				}
			}
		}
	}
	
	DAVecRestoreArray(user->da, user->lLevelset, &level);
  
	DAVecGetArray(da,  Aj,  &aj);
	DAVecGetArray(fda,  user->lKZet,  &kzet);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, user->lUcont, &ucont);
	DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, DRHS, &rhs);
  
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double U_dpdc=0, U_dpde=0, U_dpdz=0;
		
		if (nvert[k][j][i]>0.1) {
			rhs[k][j][i]=0.;
			continue;
		}
		
		double densityL, densityR;
		double aL=ucont[k][j][i-1].x, aR=ucont[k][j][i].x;	// wave speed
		
		if ( i==mx-2 ) {
			densityL = level[k][j][i-1] + 0.5*M(level[k][j][i]-level[k][j][i-1], level[k][j][i-1]-level[k][j][i-2]);
			
			if( !i_periodic && !ii_periodic ) densityR = level[k][j][i];
			else densityR = level[k][j][i] + 0.5*M(level[k][j][i+1]-level[k][j][i], level[k][j][i]-level[k][j][i-1]);
		}
		else if ( i==1 ) {
			if( !i_periodic && !ii_periodic ) densityL = level[k][j][i];// - 0.5*M(level[k][j][i+1]-level[k][j][i], level[k][j][i]-level[k][j][i-1]);
			else densityL = level[k][j][i] - 0.5*M(level[k][j][i+1]-level[k][j][i], level[k][j][i]-level[k][j][i-1]);
			
			densityR = level[k][j][i+1] - 0.5*M(level[k][j][i+2]-level[k][j][i+1], level[k][j][i+1]-level[k][j][i]);
		}
		else {
			int iR=i+1, iRR=i+2;
			int iL=i-1, iLL=i-2;
			
			//if(nvert[k][j][i-1]>0.1) iL = i;
			if(nvert[k][j][i-2]>0.1) iLL = iL;
			//if(nvert[k][j][i+1]>0.1) iR = i;
			if(nvert[k][j][i+2]>0.1) iRR = iR;
			
			densityL = weno3(level[k][j][iLL],level[k][j][iL],level[k][j][i],level[k][j][iR],aL);
			densityR = weno3(level[k][j][iL],level[k][j][i],level[k][j][iR],level[k][j][iRR],aR);
		}
		//U_dpdc = densityR*ucont[k][j][i].x - densityL*ucont[k][j][i-1].x;
		U_dpdc = 0.5*(ucont[k][j][i].x+ucont[k][j][i-1].x) * (densityR - densityL);
		
		aL=ucont[k][j-1][i].y, aR=ucont[k][j][i].y;
		if ( j==my-2 ) {
			densityL = level[k][j-1][i] + 0.5*M(level[k][j][i]-level[k][j-1][i], level[k][j-1][i]-level[k][j-2][i]);
			
			if( !j_periodic && !jj_periodic ) densityR = level[k][j][i];// + 0.5*M(level[k][j+1][i]-level[k][j][i], level[k][j][i]-level[k][j-1][i]);
			else densityR = level[k][j][i] + 0.5*M(level[k][j+1][i]-level[k][j][i], level[k][j][i]-level[k][j-1][i]);
		}
		else if ( j==1 ) {
			if( !j_periodic && !jj_periodic ) densityL = level[k][j][i];// - 0.5*M(level[k][j+1][i]-level[k][j][i], level[k][j][i]-level[k][j-1][i]);
			else densityL = level[k][j][i] - 0.5*M(level[k][j+1][i]-level[k][j][i], level[k][j][i]-level[k][j-1][i]);
			
			densityR = level[k][j+1][i] - 0.5*M(level[k][j+2][i]-level[k][j+1][i], level[k][j+1][i]-level[k][j][i]);
		}
		else {
			int jR=j+1, jRR=j+2;
			int jL=j-1, jLL=j-2;
			
			//if(nvert[k][j-1][i]>0.1) jL = j;
			if(nvert[k][j-2][i]>0.1) jLL = jL;
			//if(nvert[k][j+1][i]>0.1) jR = j;
			if(nvert[k][j+2][i]>0.1) jRR = jR;
			
			densityL = weno3(level[k][jLL][i],level[k][jL][i],level[k][j][i],level[k][jR][i],aL);
			densityR = weno3(level[k][jL][i],level[k][j][i],level[k][jR][i],level[k][jRR][i],aR);
		}
		//U_dpde = densityR*ucont[k][j][i].y - densityL*ucont[k][j-1][i].y;
		U_dpde = 0.5*(ucont[k][j][i].y+ucont[k][j-1][i].y) * (densityR - densityL);
		
		aL=ucont[k-1][j][i].z, aR=ucont[k][j][i].z;
		if ( k==mz-2 ) {
			densityL = level[k-1][j][i] + 0.5*M(level[k][j][i]-level[k-1][j][i], level[k-1][j][i]-level[k-2][j][i]);
			
			if( user->bctype[5]!=4 &&  (!k_periodic && !kk_periodic) ) densityR = level[k][j][i];// + 0.5*M(level[k+1][j][i]-level[k][j][i], level[k][j][i]-level[k-1][j][i]);
			else if( user->bctype[5]==4) densityR = 0.5 * ( level[k][j][i]+level[k+1][j][i] );// outlet condition
			else densityR = level[k][j][i] + 0.5*M(level[k+1][j][i]-level[k][j][i], level[k][j][i]-level[k-1][j][i]);
		}
		else if ( k==1 ) {
			if( user->bctype[4]!=5 &&  (!k_periodic && !kk_periodic)  ) densityL = level[k][j][i];// - 0.5*M(level[k+1][j][i]-level[k][j][i], level[k][j][i]-level[k-1][j][i]);
			else if( user->bctype[4]==5) densityL = 0.5 * ( level[k][j][i]+level[k-1][j][i] );// intlet condition
			else densityL = level[k][j][i] - 0.5*M(level[k+1][j][i]-level[k][j][i], level[k][j][i]-level[k-1][j][i]);
			
			densityR = level[k+1][j][i] - 0.5*M(level[k+2][j][i]-level[k+1][j][i], level[k+1][j][i]-level[k][j][i]);
		}
		else {
			int kR=k+1, kRR=k+2;
			int kL=k-1, kLL=k-2;
			
			//if(nvert[k-1][j][i]>0.1) kL = k;
			if(nvert[k-2][j][i]>0.1) kLL = kL;
			//if(nvert[k+1][j][i]>0.1) kR = k;
			if(nvert[k+2][j][i]>0.1) kRR = kR;
			
			densityL = weno3(level[kLL][j][i],level[kL][j][i],level[k][j][i],level[kR][j][i],aL);
			densityR = weno3(level[kL][j][i],level[k][j][i],level[kR][j][i],level[kRR][j][i],aR);
		}
		//U_dpdz = densityR*ucont[k][j][i].z - densityL*ucont[k-1][j][i].z;
		U_dpdz = 0.5*(ucont[k][j][i].z+ucont[k-1][j][i].z) * (densityR - densityL);
		
		if( user->bctype[4]==5 &&  k<=inlet_buffer_k ) {
			// continuously feeds inlet levelset in the horizontal direction
			/*
			U_dpde = 0; 
			U_dpdc = 0; 
			*/
		}
			
		rhs[k][j][i] = - ( U_dpdc + U_dpde + U_dpdz ) * aj[k][j][i];	// advection
		
		if( user->bctype[4]==5 && user->bctype[5]==4 ) {
		  //if ( k<=inlet_buffer_k ) { rhs[k][j][i] = 0; }
			//if ( fix_outlet && k==mz-2 ) { rhs[k][j][i] = 0; }
		}
		
	}
	
	DAVecRestoreArray(fda, user->lCent, &cent);
	DAVecRestoreArray(da,  Aj,  &aj);
	DAVecRestoreArray(fda,  user->lKZet,  &kzet);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, user->lUcont, &ucont);
	DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, DRHS, &rhs);
	
	//DALocalToGlobal(user->da, user->lLevelset, INSERT_VALUES, user->Levelset);
}

void Advect_Levelset(UserCtx *user, double dt)
{

	//VecCopy(user->Levelset, user->Levelset_o);
	Vec R0, R1;
	VecDuplicate(user->Levelset, &R0);	// allocation
	VecDuplicate(user->Levelset, &R1);	// allocation
	
	Levelset_Advect_RHS(user, R0);
	VecAXPY(user->Levelset, dt, R0);        /* U(1) = U(n) + dt * RHS(n) */
	
	VecWAXPY(user->Levelset, dt, R0, user->Levelset_o);	
	Levelset_Advect_RHS(user, R1);
	VecWAXPY(user->Levelset, 0.5*dt, R0, user->Levelset_o);
	VecAXPY(user->Levelset, 0.5*dt, R1);
	
	VecDestroy(R0);		// free
	VecDestroy(R1);		// free
	
	DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
	
	
}
void Compute_Surface_Tension(UserCtx *user)
{
	DA 		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs, xe, ys, ye, zs, ze;
	PetscInt	mx, my, mz;
	PetscInt	i, j, k;
	PetscInt	lxs, lys, lzs, lxe, lye, lze;

	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;
  	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  	mx = info.mx; my = info.my; mz = info.mz;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	Vec Curv, Grad_abs, Heaviside;
	Cmpnts ***curv, ***stension;
	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal ***aj, ***iaj, ***jaj, ***kaj;
	PetscReal ***level, ***nvert, ***density, ***grad, ***h;
	
	
	VecDuplicate(user->lUcont, &Curv);
	VecDuplicate(user->lP, &Grad_abs);
	VecDuplicate(user->lP, &Heaviside);
	
	DAVecGetArray(fda, Curv, &curv);
	DAVecGetArray(da, Grad_abs, &grad);
	DAVecGetArray(da, Heaviside, &h);
	DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, user->lDensity, &density);
	DAVecGetArray(da, user->lNvert, &nvert);
	
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
	
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);
	
	DAVecGetArray(fda, user->lST, &stension);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double dldc, dlde, dldz;
		double dl_dx, dl_dy, dl_dz;
		double ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		
		Compute_dscalar_center (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		grad[k][j][i] = sqrt ( dl_dx*dl_dx + dl_dy*dl_dy + dl_dz*dl_dz );
		double dx=pow(1./aj[k][j][i],1./3.);
		if(dthick_set) dx = dthick;

		h[k][j][i] = H ( level[k][j][i], dx );
		/*
		int c=k, b=j, a=i, flag=0;
		
		// Neumann conditions
		if(i==1) a=0, flag=1;
		if(i==mx-2) a=mx-1, flag=1;
		if(j==1) b=0, flag=1;
		if(j==my-2) b=my-1, flag=1;
		if(k==1) c=0, flag=1;
		if(k==mz-2) c=mz-1, flag=1;
		
		if(flag) {
			grad[c][b][a] = grad[k][j][i];
			h[c][b][a] = h[k][j][i];
		}
		*/
	}
	
	DAVecRestoreArray(da, Grad_abs, &grad);
	DAVecRestoreArray(da, Heaviside, &h);
	
	DALocalToLocalBegin(da, Grad_abs, INSERT_VALUES, Grad_abs);
	DALocalToLocalEnd(da, Grad_abs, INSERT_VALUES, Grad_abs);
	
	DALocalToLocalBegin(da, Heaviside, INSERT_VALUES, Heaviside);
	DALocalToLocalEnd(da, Heaviside, INSERT_VALUES, Heaviside);
	
	DAVecGetArray(da, Grad_abs, &grad);
	DAVecGetArray(da, Heaviside, &h);
	
		
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				grad[k][j][to] = grad[k][j][from];
				h[k][j][to] = h[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				grad[k][j][to] = grad[k][j][from];
				h[k][j][to] = h[k][j][from];
			}
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
				
				grad[k][to][i] = grad[k][from][i];
				h[k][to][i] = h[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				grad[k][to][i] = grad[k][from][i];
				h[k][to][i] = h[k][from][i];
			}
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
				
				grad[to][j][i] = grad[from][j][i];
				h[to][j][i] = h[from][j][i];
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(k_periodic) from = 1;
				else if(kk_periodic) from = mz+1;
				
				grad[to][j][i] = grad[from][j][i];
				h[to][j][i] = h[from][j][i];
			}
		}
	}
	
	
	// i direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(j==0 || k==0) continue;
		
		double ajc = iaj[k][j][i];
		double csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		double eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		double zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
		
		double dldc, dlde, dldz;
		Compute_dscalar_i (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
		
		double dl_dx, dl_dy, dl_dz;
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		double abs_grad = sqrt ( dl_dx*dl_dx + dl_dy*dl_dy + dl_dz*dl_dz );// + 1.e-20;
		
		//double abs_grad = 0.5 * ( grad[k][j][i] + grad[k][j][i+1] );
		//if(tistart==0 && ti<25) abs_grad+=1.e-20;

		double g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		double g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		double g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
		
		if(abs_grad<1.e-10) curv[k][j][i].x=0.;
		else curv[k][j][i].x = (g11 * dldc + g21 * dlde + g31 * dldz) * ajc / abs_grad;
		if (!i_periodic && !ii_periodic && (i==0 || i==mx-2)) curv[k][j][i].x = 0;
		if ( nvert[k][j][i] + nvert[k][j][i+1] > 0.1 ) curv[k][j][i].x = 0;
	}
	
	// j direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(i==0 || k==0) continue;

		double ajc = jaj[k][j][i];
		double csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		double eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		double zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
		
		double dldc, dlde, dldz;
		Compute_dscalar_j (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
		
		
		double dl_dx, dl_dy, dl_dz;
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		double abs_grad = sqrt ( dl_dx*dl_dx + dl_dy*dl_dy + dl_dz*dl_dz );// + 1.e-20;
		
		//double abs_grad = 0.5 * ( grad[k][j][i] + grad[k][j+1][i] );
		//if(tistart==0 && ti<25) abs_grad+=1.e-20;

		double g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
		double g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		double g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
		
		if(abs_grad<1.e-10) curv[k][j][i].y=0.;
		else curv[k][j][i].y = (g11 * dldc + g21 * dlde + g31 * dldz) * ajc / abs_grad;
		if(!j_periodic && !jj_periodic && (j==0 || j==my-2)) curv[k][j][i].y = 0;
		if ( nvert[k][j][i] + nvert[k][j+1][i] > 0.1 ) curv[k][j][i].y = 0;
	}
	
	// k direction
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==mx-1 || j==my-1 || k==mz-1) continue;
		if(i==0 || j==0) continue;
		
		double ajc = kaj[k][j][i];
		double csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		double eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		double zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		
		double dldc, dlde, dldz;
		Compute_dscalar_k (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
		
		double dl_dx, dl_dy, dl_dz;
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		double abs_grad = sqrt ( dl_dx*dl_dx + dl_dy*dl_dy + dl_dz*dl_dz );// + 1.e-20;
		
		//double abs_grad = 0.5 * ( grad[k][j][i] + grad[k+1][j][i] );
		//if(tistart==0 && ti<25) abs_grad+=1.e-20;

		double g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
		double g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
		double g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
		
		if(abs_grad<1.e-10) curv[k][j][i].z=0.;
		else curv[k][j][i].z = (g11 * dldc + g21 * dlde + g31 * dldz) * ajc / abs_grad;
		if(!k_periodic && !kk_periodic && (k==0 || k==mz-2)) curv[k][j][i].z = 0;
		if ( nvert[k][j][i] + nvert[k+1][j][i] > 0.1 ) curv[k][j][i].z = 0;
	}
		
	DAVecRestoreArray(fda, Curv, &curv);
		
	DALocalToLocalBegin(fda, Curv, INSERT_VALUES, Curv);
	DALocalToLocalEnd(fda, Curv, INSERT_VALUES, Curv);
		
	DAVecGetArray(fda, Curv, &curv);
	

	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				curv[k][j][to] = curv[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				curv[k][j][to] = curv[k][j][from];
			}
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
				
				curv[k][to][i] = curv[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				curv[k][to][i] = curv[k][from][i];
			}
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
				
				curv[to][j][i] = curv[from][j][i];
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(k_periodic) from = 1;
				else if(kk_periodic) from = mz+1;
				
				curv[to][j][i] = curv[from][j][i];
			}
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double phi = level[k][j][i];
		double ajc = aj[k][j][i];
		double kappa = (curv[k][j][i].x - curv[k][j][i-1].x + curv[k][j][i].y - curv[k][j-1][i].y + curv[k][j][i].z - curv[k-1][j][i].z) * ajc;
		double dx=pow(1./aj[k][j][i],1./3.);
		if(dthick_set) dx = dthick;

		double gradH = dH (phi, dx); 
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double dldc, dlde, dldz, dl_dx, dl_dy, dl_dz;
		double dhdc, dhde, dhdz, dh_dx, dh_dy, dh_dz;
		
		double sigma = 7.28e-2;
		
		Compute_dscalar_center (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);
		
		Compute_dscalar_center (i, j, k, mx, my, mz, h, nvert, &dhdc, &dhde, &dhdz );
		Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dhdc, dhde, dhdz, &dh_dx, &dh_dy, &dh_dz);
		
		double A = sigma * kappa / density[k][j][i];// * gradH;// / density[k][j][i];
		
		
		stension[k][j][i].x = - A * dl_dx * gradH;
		stension[k][j][i].y = - A * dl_dy * gradH;
		stension[k][j][i].z = - A * dl_dz * gradH;
		
		if ( nvert[k][j][i]> 0.1 ) Set (&stension[k][j][i], 0);
	
	}
	
	DAVecRestoreArray(fda, Curv, &curv);
	DAVecRestoreArray(da, Grad_abs, &grad);
	DAVecRestoreArray(da, Heaviside, &h);
	DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, user->lDensity, &density);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	
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
	
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);
	
	DAVecRestoreArray(fda, user->lST, &stension);
	
	DALocalToLocalBegin(fda, user->lST, INSERT_VALUES, user->lST);
	DALocalToLocalEnd(fda, user->lST, INSERT_VALUES, user->lST);
	
	DAVecGetArray(fda, user->lST, &stension);
	
	if(xs==0 || xe==mx) {
		int from, to;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			if(xs==0) {
				i = 1, from = i, to = 0;
				
				if(i_periodic) from = mx-2;
				else if(ii_periodic) from = -2;
				
				stension[k][j][to] = stension[k][j][from];
			}
			
			if(xe==mx) {
				i = mx-2, from = i, to = mx-1;
				
				if(i_periodic) from = 1;
				else if(ii_periodic) from = mx+1;
				
				stension[k][j][to] = stension[k][j][from];
			}
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
				
				stension[k][to][i] = stension[k][from][i];
			}
			
			if(ye==my) {
				j = my-2, from = j, to = my-1;
				
				if(j_periodic) from = 1;
				else if(jj_periodic) from = my+1;
				
				stension[k][to][i] = stension[k][from][i];
			}
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
				
				stension[to][j][i] = stension[from][j][i];
			}
			
			if(ze==mz) {
				k = mz-2, from = k, to = mz-1;
				
				if(k_periodic) from = 1;
				else if(kk_periodic) from = mz+1;
				
				stension[to][j][i] = stension[from][j][i];
			}
		}
	}
	DAVecRestoreArray(fda, user->lST, &stension);
	
	
	VecDestroy(Curv);
	VecDestroy(Grad_abs);
	VecDestroy(Heaviside);
}

double H (double p, double dx)
{
	double eps;
	
	if(dthick_set) eps = dthick * 1.5;
	else eps = dx * 1.5;
 
	if ( p < -eps ) return 0;
	else if ( fabs(p) <= eps ) return 0.5 * ( 1.0 + p/eps + 1./M_PI * sin ( M_PI * p / eps ) );
	else return 1;
};

double dH (double p, double dx)
{
	double eps;
	
	if(dthick_set) eps = dthick * 1.5;
	else eps = dx * 1.5;
	
	if ( p < -eps ) return 0;
	else if ( fabs(p) <= eps ) return 0.5 / eps * ( 1.0 + cos ( M_PI * p / eps ) );
	else return 0;
};

// Puckett
double mean ( double A, double B )
{
	return 2.0 * A * B / ( A + B );
	//return 0.5 * ( A + B );
}

double sign1(double a, double dx)
{
  //return sign(a);
  	return 2.0 * ( H (a, dx) - 0.5 );
};

// Yue & Patel
double mod_sign(double d0, double grad0, double e)
{
	return d0 / sqrt ( d0*d0 + pow(grad0*e, 2.0) );
}


