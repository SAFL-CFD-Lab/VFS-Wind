/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

double threshold=0.1;

extern PetscReal Flux_in, angle;
extern PetscInt block_number;
extern int pseudo_periodic;
extern void Calculate_Covariant_metrics(double g[3][3], double G[3][3]);
extern int initial_perturbation;

extern void WIND_vel_interpolate(UserCtx *user, double *u1, double *v1, double *w1, double y1, double z1, int j, int k);

#define INLET 5
#define OUTLET 4
#define SOLIDWALL 1
#define SYMMETRIC 3
#define FARFIELD 6

int ucont_plane_allocated=0;

void pseudo_periodic_BC(UserCtx *user)
{
	DALocalInfo info = user->info;
	int	xs = info.xs, xe = info.xs + info.xm;
	int 	ys = info.ys, ye = info.ys + info.ym;
	int	zs = info.zs, ze = info.zs + info.zm;
	int	mx = info.mx, my = info.my, mz = info.mz;

	int i, j, k;
	int kplane=(mz*6)/10;
	
	if ( !ucont_plane_allocated ) {
		ucont_plane_allocated = 1;
		
		user->ucont_plane = (Cmpnts **)malloc( sizeof(Cmpnts *) * my );
		if(rans) user->komega_plane = (Cmpnts2 **)malloc( sizeof(Cmpnts2 *) * my );
		
		for(j=0; j<my; j++) {
			user->ucont_plane[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * mx );
			if(rans) user->komega_plane[j] = (Cmpnts2 *)malloc( sizeof(Cmpnts2) * mx );
		}
	}
	
	Cmpnts ***ucat;
	Cmpnts2 ***komega;
	
	std::vector< std::vector<Cmpnts> > ucont_plane_tmp (my);
	std::vector< std::vector<Cmpnts2> > komega_plane_tmp (my);
	
	for( j=0; j<my; j++) ucont_plane_tmp[j].resize(mx);
	for( j=0; j<my; j++) komega_plane_tmp[j].resize(mx);
	
	
	for(j=0; j<my; j++)
	for(i=0; i<mx; i++) {
		ucont_plane_tmp[j][i].x = ucont_plane_tmp[j][i].y = ucont_plane_tmp[j][i].z = 0;
		komega_plane_tmp[j][i].x = komega_plane_tmp[j][i].y = 0;
	}
	
	DAVecGetArray(user->fda, user->Ucat, &ucat);
	if(rans) DAVecGetArray(user->fda2, user->K_Omega, &komega);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(k==kplane) {
			ucont_plane_tmp[j][i] = ucat[k][j][i];
			if(rans) komega_plane_tmp[j][i] = komega[k][j][i];
		}
	}
	
	DAVecRestoreArray(user->fda, user->Ucat, &ucat);
	if(rans) DAVecRestoreArray(user->fda2, user->K_Omega, &komega);
	
	for(j=0; j<my; j++) {
		MPI_Allreduce( &ucont_plane_tmp[j][0], &user->ucont_plane[j][0], mx*3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		if(rans) MPI_Allreduce( &komega_plane_tmp[j][0], &user->komega_plane[j][0], mx*2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	}

};

PetscReal InletInterpolation(PetscReal r, UserCtx *user)
{
  PetscInt i;
  PetscReal temp;
  PetscInt tstep;
  tstep = ti - 1000 * (ti/1000);

/*   return(0); */

/*   return (2*(1-(r*r)*4)); */
/*   if (ti>0)  */
/*   return(1.); */
/*   else */
/*    return 1.; */
  return(Flux_in);

/*   if (r>1.) { */
/*     temp = 0.; */
/*     return (temp); */
/*   } */
/*   else { */
/*     temp = 2 * (1 - (r*2) * (r*2)); */
/*     return(temp); */
/*   } */
  
/*   if (ti==0) { */
/*     temp = 0.; */
/*     return(temp); */
/*   } */
  if (r>1.) {
    temp = user->uinr[99][tstep];
    //    PetscPrintf(PETSC_COMM_WORLD, "Inlet %e", temp);
    return(temp);
  }
  for (i=0; i<100; i++) {
    if (r>= (user->r[i]) && r< (user->r[i+1])) {
      temp = user->uinr[i][tstep] + (user->uinr[i+1][tstep] - user->uinr[i][tstep]) *
	(r-user->r[i]) / (user->r[i+1]-user->r[i]);
      //      PetscPrintf(PETSC_COMM_WORLD, "Inlet %e", temp);
      return (temp);
    }
  }   
}

double inletArea=0;
int k_area_allocated=0;

void Calc_Inlet_Area(UserCtx *user)
{
	inletArea=0;
		
	PetscInt i, j, k;

	Cmpnts	***csi, ***eta, ***zet;
	
	DA da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
 
	PetscReal	***nvert, ***level, ***aj;
	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	if(levelset)  DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
	DAVecGetArray(da, user->lAj, &aj);

	if(!k_area_allocated) {
	  k_area_allocated=1;
		user->k_area = new double [mz];
	}
        
	std::vector<double> lArea(mz);
	
	std::fill ( lArea.begin(), lArea.end(), 0 );

        for(k=lzs; k<lze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if (nvert[k-1][j][i]+nvert[k][j][i] < 0.1) {
			if(j>=1 && j<=my-2 && i>=1 && i<=mx-2) {
				double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );			
				if (air_flow_levelset_periodic && level[k][j][i]>-dthick) area=0.;
				lArea[k] += area;
			}
		}
	}

        MPI_Allreduce( &lArea[0], &user->k_area[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	user->k_area[0] = user->k_area[1];
	for (k=1; k<=mz-2; k++) user->mean_k_area += user->k_area[k];
        user->mean_k_area/= (double)(mz-2);

	if (user->bctype[4] == INLET || k_periodic || kk_periodic) {
		double lArea=0;
		if (zs==0) {
			k = 0;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {

				double k_area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				double j_area = sqrt( eta[k+1][j][i].x*eta[k+1][j][i].x + eta[k+1][j][i].y*eta[k+1][j][i].y + eta[k+1][j][i].z*eta[k+1][j][i].z );

				if (nvert[k+1][j][i] < threshold) {
				  /*if( levelset && level[k][j][i]<0) { 
				    }*/
				  
				  if( levelset ) {
				    double dx = pow(1./aj[k][j][i], 1./3.);
				    if(dthick_set) dx=dthick;
				    double vf=0;
						if (!air_flow_levelset) vf= H(level[k+1][j][i], dx);
						if (air_flow_levelset && level[k][j][i]<=-dthick)vf = 1.0;// H(-level[k+1][j][i], dx);
						if (air_flow_levelset && level[k][j][i]>-dthick)vf = 0.;// H(-level[k+1][j][i], dx);
						if (air_flow_levelset_periodic && level[k][j][i]<=-dthick)vf = 1.0;// H(-level[k+1][j][i], dx);
						if (air_flow_levelset_periodic && level[k][j][i]>-dthick)vf = 0.;// H(-level[k+1][j][i], dx);						
				    lArea += k_area * vf;
				  }
				  
				  else {
				    lArea += k_area;
				  }
				}
			}
		}
		
		PetscGlobalSum(&lArea, &inletArea, PETSC_COMM_WORLD);    
	}
       
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	DAVecRestoreArray(da, user->lAj, &aj);
	if(levelset)  DAVecRestoreArray (da, user->lLevelset, &level);
	
	PetscPrintf(PETSC_COMM_WORLD, "\n*** (Fluid) Inlet Area:%f, (Fluid) Inlet Flux: %f\n\n", inletArea, inlet_flux);
};

PetscErrorCode InflowFlux(UserCtx *user) 
{
	
  
	PetscInt i, j, k;
	PetscReal r, uin, xc, yc;
	Vec Coor;
	Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet;
	Cmpnts ***icsi;
	
	PetscReal  /*H=4.1,*/ Umax=1.5;

	DA da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert, ***level, ***aj;	//seokkoo
	
	DAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
	if(levelset) DAVecGetArray(da, user->lLevelset, &level);
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	DAGetGhostedCoordinates(da, &Coor);
	DAVecGetArray(fda, Coor, &coor);
	DAVecGetArray(fda, user->Ucont, &ucont);
	DAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
	DAVecGetArray(fda, user->Ucat,  &ucat);
  
	DAVecGetArray(fda, user->lCsi,  &csi);
	DAVecGetArray(fda, user->lEta,  &eta);
	DAVecGetArray(fda, user->lZet,  &zet);
	DAVecGetArray(da, user->lAj,  &aj);

	DAVecGetArray(fda, user->lICsi,  &icsi);
  
	PetscReal FluxIn=0., FluxIn_gas=0.;
	
	double lFluxIn0=0, sumFluxIn0=0;
	double lFluxIn1=0, sumFluxIn1=0;
  
	srand( time(NULL)) ;	// seokkoo
	int fluct=0;
      
	if (user->bctype[4] == INLET) {
		if (zs==0) {
			k = 0;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				xc = (coor[k+1][j][i].x + coor[k+1][j-1][i].x + coor[k+1][j][i-1].x + coor[k+1][j-1][i-1].x) * 0.25;
				yc = (coor[k+1][j][i].y + coor[k+1][j-1][i].y + coor[k+1][j][i-1].y + coor[k+1][j-1][i-1].y) * 0.25;
				r = sqrt(xc * xc + yc * yc);
				double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				
				if (inletprofile == 0) uin = InletInterpolation(r, user);
				else if (inletprofile==1) {
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
				}
				else if (inletprofile==2) {	// uniform flow with noise
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
				}
				else if (inletprofile == 3) uin = InletInterpolation(r, user);
				
				
				else if (inletprofile==10) {	// Power law for hemisphere case
					//fluct=1;
					double delta = 0.45263, a=5.99;
					if( yc>=delta ) uin=1.0;
					else if(yc<=0) uin=0.0;
					else uin = pow( yc/delta, 1./a );
				} 
				else if (inletprofile==11) {	// backward facing step
					fluct=1;
					yc -= 3;
					double delta=1.2;
					if( 2-fabs(yc)<delta ) uin = pow( (2-fabs(yc))/delta, 1./7 );
					else uin=1;
				} 
				else if (inletprofile==12) {	// pipe shear stress test
					//		fluct=1;
					uin = 2*(1 - pow(r/0.5,2.0) );
				}
				else if (inletprofile==13) {	// periodic channel flow
				}
				else if (inletprofile==101) {	// periodic channel flow
				}
				else if(inletprofile==14) {	// curved pipe flow (Anwer)
					double R[11] = {0.000,0.062,0.125,0.188,0.250,0.312,0.354,0.399,0.438,0.469,0.500};
					double W[11] = {1.120,1.115,1.110,1.093,1.050,1.010,0.963,0.915,0.825,0.730,0.000};
					
					int ii;
					for(ii=1; ii<11; ii++) {
						if( R[ii]>=r && R[ii-1]<r ) break;
					}
					uin = ( W[ii] - W[ii-1] ) / ( R[ii] - R[ii-1] ) * ( r - R[ii-1] ) + W[ii-1];
					
					if(r>0.5) uin = 0;
				}
				else if(inletprofile==15) {	// round jet (Longmire)
					if(r>0.5) uin = 0;
					else uin = 1;
				}
				else if (inletprofile==16) {	// periodic pipe flow
				}
				else if(inletprofile==20){}// enright test
				else if (inletprofile==21) {    // laminar poiseulle flow profile
				  double w_bulk=1.0;//inlet_flux/SumArea;
				  uin = 1.5 * w_bulk * yc * ( 2 - yc );
				}
				else if (inletprofile == 100) {}	// saved data for LES
				
				else {
					PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
					uin = 0.;
				}
				
				if(pseudo_periodic || inletprofile==100) {	// pseudo-periodic BC in k-direction
					fluct=0;
				}
				
				if (nvert[k+1][j][i] < threshold) {	//seokkoo
					if(pseudo_periodic || inletprofile==100) {
						double u = user->ucont_plane[j][i].x;
						double v = user->ucont_plane[j][i].y;
						double w = user->ucont_plane[j][i].z;
						
						ucat[k][j][i].x = u;
						ucat[k][j][i].y = v;
						ucat[k][j][i].z = w;
						ucont[k][j][i].z = 0.5*(ucat[k][j][i].x+ucat[k+1][j][i].x) * zet[k][j][i].x + 0.5*(ucat[k][j][i].y+ucat[k+1][j][i].y) * zet[k][j][i].y + 0.5*(ucat[k][j][i].z+ucat[k+1][j][i].z) * zet[k][j][i].z;
					}
					else if(air_flow_levelset==2){
						double t=user->dt*(double)ti;
						double x = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
						double y = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
						//printf("x coordinate: %f, y coordinate: %f \n",x,y);
						double u=0.;//streamwise
						double v=0.;
						double w=0.;//vertical;
						//printf("inflow flix\n");
						if(level[k][j][i]<=-dthick){
							if(x>user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1])x=2*user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1]-x;
							if(x<user[0].wave_inf[0].WIND_Y[0])x=2*user[0].wave_inf[0].WIND_Y[0]-x;
							WIND_vel_interpolate(&user[0], &u,  &v,  &w, x, y, user[0].wave_inf[0].WIND_id_x[i-lxs][j-lys], user[0].wave_inf[0].WIND_id_y[i-lxs][j-lys]);	
							if(x>user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1] || x<user[0].wave_inf[0].WIND_Y[0])u=-u;
							//printf("%f %f %f \n",u,v,w);
						}
						else{
							u=0.,v=0.,w=0.;
						}
						ucat[k][j][i].x = - ucat[k+1][j][i].x + 2*u;
						ucat[k][j][i].y = - ucat[k+1][j][i].y + 2*v;
						ucat[k][j][i].z = - ucat[k+1][j][i].z + 2*w;							
						ucont[k][j][i].x = 0.5*(ucat[k][j][i].x+ucat[k+1][j][i].x) * csi[k][j][i].x + 0.5*(ucat[k][j][i].y+ucat[k+1][j][i].y) * csi[k][j][i].y + 0.5*(ucat[k][j][i].z+ucat[k+1][j][i].z) * csi[k][j][i].z;
						ucont[k][j][i].y = 0.5*(ucat[k][j][i].x+ucat[k+1][j][i].x) * eta[k][j][i].x + 0.5*(ucat[k][j][i].y+ucat[k+1][j][i].y) * eta[k][j][i].y + 0.5*(ucat[k][j][i].z+ucat[k+1][j][i].z) * eta[k][j][i].z;
						ucont[k][j][i].z = 0.5*(ucat[k][j][i].x+ucat[k+1][j][i].x) * zet[k][j][i].x + 0.5*(ucat[k][j][i].y+ucat[k+1][j][i].y) * zet[k][j][i].y + 0.5*(ucat[k][j][i].z+ucat[k+1][j][i].z) * zet[k][j][i].z;
						fluct=0;
					}							
					else {
						if(levelset && !air_flow_levelset) {
						  if( level[k][j][i]<0 ) { 
						    ucont[k][j][i].z = uin * area * (1.);//haha 
						  }
						  else ucont[k][j][i].z = uin * area;
						}
						else if(levelset && air_flow_levelset) {
						  if( level[k][j][i]<=-dthick ) {
						    ucont[k][j][i].z = uin * area * (1.);
						  }
						  else ucont[k][j][i].z = 0.;//uin * area;
						}	
						else ucont[k][j][i].z = uin * area;
						
						ucont[k][j][i].x = 0;
						ucont[k][j][i].y = 0;
						Cmpnts u;
						Contra2Cart_single(csi[k][j][i], eta[k][j][i], zet[k][j][i], ucont[k][j][i], &u);
						
						ucat[k][j][i].x = - ucat[k+1][j][i].x + 2*u.x;
						ucat[k][j][i].y = - ucat[k+1][j][i].y + 2*u.y;
						ucat[k][j][i].z = - ucat[k+1][j][i].z + 2*u.z;
						
					}
					ubcs[k][j][i] = ucat[k][j][i];
				}
				else {
					ucat[k][j][i].z = 0;	//seokkoo
					ubcs[k][j][i].z = 0;
					ucont[k][j][i].z = 0;
				}
				if(j>=0 && j<=my-2 && i>=0 && i<=mx-2) {
					lFluxIn0 += ucont[k][j][i].z;

					if(fluct) {
						double F = 0.03;			// 3% noise
						int n = rand() % 20000;		// RAND_MAX = 65535
						n -= 10000;
						ucont[k][j][i].z *= ( 1 + ((double)n)/10000. * F );	// uin * (1+-0.xx)
					}
					if( levelset && !air_flow_levelset) {
					  double dx = pow(1./aj[k][j][i], 1./3.);
					  if(dthick_set) dx=dthick;
					  double vf = H(level[k+1][j][i], dx);
					  lFluxIn1 += ucont[k][j][i].z * vf;
					}
					else if(levelset && air_flow_levelset){
						if (level[k][j][i]<=-dthick)lFluxIn1 += ucont[k][j][i].z;			
					}
					else lFluxIn1 += ucont[k][j][i].z;
					
				}
			}
		}
		
		PetscGlobalSum(&lFluxIn0, &sumFluxIn0, PETSC_COMM_WORLD);
		PetscGlobalSum(&lFluxIn1, &sumFluxIn1, PETSC_COMM_WORLD);
			
		extern int ti;
		if(pseudo_periodic || inletprofile==100 || inletprofile==14) {
		  //if(inlet_flux<0) inlet_flux=1.0*inletArea;	// bulk velocity is unity.
			PetscPrintf(PETSC_COMM_WORLD,  "\nConstant Flux is %f !\n\n", inlet_flux);
			sumFluxIn0=inlet_flux;
		}
		
		PetscPrintf(PETSC_COMM_WORLD, "\n***** Fluxin0:%f, Fluxin1:%f, Area:%f\n\n", sumFluxIn0, sumFluxIn1, inletArea);
		
		if (zs==0 && (fluct || pseudo_periodic || inletprofile==100 || inletprofile==14) && inlet_flux>0 ) {
			k = 0;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				//double A = zet[k][j][i].z;
				double A = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				if (nvert[k+1][j][i] < threshold) ucont[k][j][i].z += (sumFluxIn0 - sumFluxIn1) * A / inletArea;
			}
		}
		
		/* initial perturbation */
		/*
		if(ti==0 && (pseudo_periodic || initial_perturbation)) {
			PetscPrintf(PETSC_COMM_WORLD, "\nGenerating initial perturbation\n");
			for(k=0; k<mz; k++) {	// for all CPU
				lFluxIn1=0;
				if( k>zs && k<lze && k!=0) {
					for (j=lys; j<lye; j++)
					for (i=lxs; i<lxe; i++) {
						//double area  = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
						if (nvert[k][j][i]+nvert[k+1][j][i] < threshold) {
							double F = 0.25;			// 30% noise
							int n = rand() % 20000;		// RAND_MAX = 65535
							n -= 10000;
							ucont[k][j][i].z *= ( 1 + ((double)n)/10000. * F );	// uin * (1+-0.xx)
							ucont[k][j][i].x += (  ((double)n)/10000. * F );	// uin * (1+-0.xx)
							ucont[k][j][i].y += (  ((double)n)/10000. * F );	// uin * (1+-0.xx)
							//printf("%f %d\n", ucont[k][j][i].z/area, k);
							lFluxIn1 += ucont[k][j][i].z;
						}
					}
				}
				PetscGlobalSum(&lFluxIn1, &sumFluxIn1, PETSC_COMM_WORLD);
				
				if( k>zs && k<lze) {
					for (j=lys; j<lye; j++)
					for (i=lxs; i<lxe; i++) {

						if (nvert[k][j][i]+nvert[k+1][j][i] < threshold) {
							double A = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
							ucont[k][j][i].z += (inlet_flux - sumFluxIn1) * A / inletArea;
							//printf("%f %f \n", ucont[k][j][i].z/A, AreaSum);	// here
						}
					}
				}
			}
		}
		*/			

		FluxIn = 0;
		if (zs==0) {
			k = 0;
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {
				if (nvert[k+1][j][i] < threshold) {
					if( levelset && air_flow_levelset ) {
						//if (level[k][j][i]<=-dthick) 
						FluxIn += ucont[k][j][i].z;
					}
					else FluxIn += ucont[k][j][i].z;
				}
			}
		}
		
       }
	else if (user->bctype[5] == INLET) {
		if (ze==mz) {
			k = mz-1;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
	  
				xc = (coor[k-1][j][i].x + coor[k-1][j-1][i].x + coor[k-1][j][i-1].x + coor[k-1][j-1][i-1].x) * 0.25;
				yc = (coor[k-1][j][i].y + coor[k-1][j-1][i].y + coor[k-1][j][i-1].y + coor[k-1][j-1][i-1].y) * 0.25;
				r = sqrt(xc * xc + yc * yc);
				
				if (inletprofile == 0) {
					uin = InletInterpolation(r, user);
				} 
				else if (inletprofile == 1) {
					uin=1.;
				} 
				else if (inletprofile == -1) {
					uin = -1.;
				} 
				else if (inletprofile == 2) {
				  //uin = 4.*Umax*yc*(H-yc)/(H*H);//InletInterpolation(r, user);
				} 
				else if (inletprofile == 3) {
					uin = InletInterpolation(r, user);
				} 
				else {
					PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
					uin = 0.;
				}
				ucat[k][j][i].z = uin;
				ubcs[k][j][i].z = uin;
				ucont[k-1][j][i].z = uin * zet[k][j][i].z;
			}
		}
		
		FluxIn = 0;
		if (ze==mz) {
			k = ze-2;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i] < threshold) FluxIn += ucont[k-1][j][i].z;
			}
		}
		else {
			FluxIn = 0;
		}
	}
        
	else if(user->bctype[0]==11) {
		FluxIn = 0;
		if (xs==0) {
			i = 0;
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				if (nvert[k][j][i+1] < threshold) {
					double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
					
					if( zc <= 0 ) {
						double u=0, v=0, w=1.;
						ucont[k][j][i].x = u * icsi[k][j][i].x + v * icsi[k][j][i].y + w * icsi[k][j][i].z;
						
						FluxIn += ucont[k][j][i].x;
					}
					else { }	// outflow
				}
			}
		}
	}
	//exit(0);
	PetscGlobalSum(&FluxIn, &FluxInSum, PETSC_COMM_WORLD);
	PetscGlobalSum(&FluxIn_gas, &FluxInSum_gas, PETSC_COMM_WORLD);	// 100203
	
	user->FluxInSum = FluxInSum;
  
	DAVecRestoreArray(fda, Coor, &coor);
	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
	DAVecRestoreArray(fda, user->Ucat,  &ucat);
  
	DAVecRestoreArray(fda, user->lCsi,  &csi);
	DAVecRestoreArray(fda, user->lEta,  &eta);
	DAVecRestoreArray(fda, user->lZet,  &zet);
	DAVecRestoreArray(da, user->lAj,  &aj);

	DAVecRestoreArray(fda, user->lICsi,  &icsi);
	
	DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo
	if(levelset) DAVecRestoreArray(da, user->lLevelset, &level);
	
	return 0;
}

PetscErrorCode OutflowFlux(UserCtx *user) {
  
  PetscInt i, j, k;
  PetscReal FluxOut;
  Vec Coor;
  Cmpnts	***ucont, ***ubcs, ***ucat, ***coor;
  


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
  
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
  DAVecGetArray(fda, user->Ucont, &ucont);
  DAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DAVecGetArray(fda, user->Ucat,  &ucat);
 
	PetscReal	***nvert;	//seokkoo
	DAVecGetArray(da, user->lNvert, &nvert);	//seokkoo  
	
	PetscReal	***level;	//seokkoo
	if(levelset) DAVecGetArray(da, user->lLevelset, &level);
	
/*     DAVecGetArray(fda, user->Csi,  &csi); */
/*     DAVecGetArray(fda, user->Eta,  &eta); */
/*     DAVecGetArray(fda, user->Zet,  &zet); */

  FluxOut = 0;
  
  if (user->bctype[5] == 4) {    
    if (ze==mz) {
      k = mz-2;
      for (j=lys; j<lye; j++) {
				for (i=lxs; i<lxe; i++) {
					if (nvert[k][j][i] < threshold){ //seokkoo
						if(air_flow_levelset){
							if(level[k][j][i]<=-dthick) FluxOut += ucont[k][j][i].z;
						}
						else{
							FluxOut += ucont[k][j][i].z;
						}
					}
				}
			}
    }
    else {
      FluxOut = 0;
    }
  } 
	else if (user->bctype[4] == 4) {    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
		if (nvert[k+1][j][i] < threshold) //seokkoo
			FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }
  }

  PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
  user->FluxOutSum = FluxOutSum;

  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DAVecRestoreArray(fda, user->Ucat,  &ucat);
  if(levelset) DAVecRestoreArray(da, user->lLevelset, &level); 
  
  DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo
  //  VecDestroy(Coor);
/*     DAVecRestoreArray(fda, user->Csi,  &csi); */
/*     DAVecRestoreArray(fda, user->Eta,  &eta); */
/*     DAVecRestoreArray(fda, user->Zet,  &zet); */
	return 0;
}

PetscErrorCode Blank_Interface(UserCtx *user) {
  PetscInt counter=0;
  PetscInt ci, cj, ck;
  
  PetscInt i;

  PetscReal ***nvert;

/*   for (bi=0; bi<block_number; bi++) { */
    DA		da =user->da;
    DALocalInfo	info = user->info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+2;
    if (ys==0) lys = ys+2;
    if (zs==0) lzs = zs+2;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-3;
    if (ze==mz) lze = ze-3;
    
    DAVecGetArray(user->da, user->Nvert, &nvert);

    for (i=0; i<user->itfcptsnumber; i++) {
      ci = user->itfcI[i];
      cj = user->itfcJ[i];
      ck = user->itfcK[i];

      if (ci>xs && ci<lxe &&
	  cj>lys && cj<lye &&
	  ck>lzs && ck<lze) {
	counter++;
	nvert[ck][cj][ci]=1.5;

      }

    }

    PetscPrintf(PETSC_COMM_WORLD, "Interface pts blanked!!!! %i,\n", counter);

    DAVecRestoreArray(user->da, user->Nvert, &nvert);

    DAGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DAGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

    return(0);
}

PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;

  Vec	hostU, nhostU;
  Cmpnts ***itfc, ***ucat;
  PetscReal *hostu;

/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCreateSeq(PETSC_COMM_SELF, */
/* 		 3*user[bi].info.mx*user[bi].info.my*user[bi].info.mz, */
/* 		 &(user[bi].nhostU)); */
/*   } */
  // First calculate Phi components at grid nodes
  for (bi=0; bi<block_number; bi++) {
    DALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    DAVecGetArray(user[bi].fda, user[bi].Itfc, &itfc);
    DAVecGetArray(user[bi].fda, user[bi].lUcat, &ucat);

    DACreateNaturalVector(user[bi].fda, &nhostU);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (k<mz-1 && j<my-1 && i<mx-1) {
	    itfc[k][j][i].x = 0.125 * (ucat[k  ][j  ][i  ].x +
				       ucat[k  ][j  ][i+1].x +
				       ucat[k  ][j+1][i  ].x +
				       ucat[k  ][j+1][i+1].x +
				       ucat[k+1][j  ][i  ].x +
				       ucat[k+1][j  ][i+1].x +
				       ucat[k+1][j+1][i  ].x +
				       ucat[k+1][j+1][i+1].x);

	    itfc[k][j][i].y = 0.125 * (ucat[k  ][j  ][i  ].y +
				       ucat[k  ][j  ][i+1].y +
				       ucat[k  ][j+1][i  ].y +
				       ucat[k  ][j+1][i+1].y +
				       ucat[k+1][j  ][i  ].y +
				       ucat[k+1][j  ][i+1].y +
				       ucat[k+1][j+1][i  ].y +
				       ucat[k+1][j+1][i+1].y);
	    itfc[k][j][i].z = 0.125 * (ucat[k  ][j  ][i  ].z +
				       ucat[k  ][j  ][i+1].z +
				       ucat[k  ][j+1][i  ].z +
				       ucat[k  ][j+1][i+1].z +
				       ucat[k+1][j  ][i  ].z +
				       ucat[k+1][j  ][i+1].z +
				       ucat[k+1][j+1][i  ].z +
				       ucat[k+1][j+1][i+1].z);
	  }
	}
      }
    }
    DAVecRestoreArray(user[bi].fda, user[bi].Itfc, &itfc);
    DAVecRestoreArray(user[bi].fda, user[bi].lUcat, &ucat);

    DAGlobalToNaturalBegin(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    DAGlobalToNaturalEnd(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    //    DALocalToGlobal(user[bi].fda, user[bi].lItfc, INSERT_VALUES, user[bi].Itfc);
    VecScatter tolocalall;
    VecScatterCreateToAll(nhostU, &tolocalall, &(user[bi].nhostU));
    VecScatterBegin(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterEnd(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(tolocalall);
    VecDestroy(nhostU);
  }

  VecCreateSeq(PETSC_COMM_SELF, 24, &hostU);
  for (bi=0; bi<block_number; bi++) {
    DALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    

    PetscInt    lmx, lmy, lmz;
    DAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];


      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);
      }
      // Is the point a local point?
      // Get the host cell information from the host CPU
      // Update
    }
    PetscBarrier(PETSC_NULL);
    DAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DALocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DALocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  }
  VecDestroy(hostU);
  
  for (bi=0; bi<block_number; bi++) {

    DALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts ***ucont, ***kzet, ***jeta;

    DAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DAVecGetArray(user[bi].fda, user[bi].lKZet, &kzet);
    DAVecGetArray(user[bi].fda, user[bi].lJEta, &jeta);
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if (k>0 && j>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if (k>0 && j>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=1;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }




    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (j>0 && i>0) {
	    ucont[k][j][i].z = (0.25 * (itfc[k  ][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j-1][i  ].x +
					itfc[k  ][j-1][i-1].x) *
				kzet[k][j][i].x +
				0.25 * (itfc[k  ][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j-1][i  ].y +
					itfc[k  ][j-1][i-1].y) *
				kzet[k][j][i].y +
				0.25 * (itfc[k  ][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j-1][i  ].z +
					itfc[k  ][j-1][i-1].z) *
				kzet[k][j][i].z);
				    
/* 	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x + */
/* 				       itfc[k  ][j  ][i-1].x + */
/* 				       itfc[k  ][j-1][i  ].x + */
/* 				       itfc[k-1][j  ][i  ].x + */
/* 				       itfc[k  ][j-1][i-1].x + */
/* 				       itfc[k-1][j  ][i-1].x + */
/* 				       itfc[k-1][j-1][i  ].x + */
/* 				       itfc[k-1][j-1][i-1].x); */
	
/* 	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y + */
/* 				       itfc[k  ][j  ][i-1].y + */
/* 				       itfc[k  ][j-1][i  ].y + */
/* 				       itfc[k-1][j  ][i  ].y + */
/* 				       itfc[k  ][j-1][i-1].y + */
/* 				       itfc[k-1][j  ][i-1].y + */
/* 				       itfc[k-1][j-1][i  ].y + */
/* 				       itfc[k-1][j-1][i-1].y); */
	
/* 	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z + */
/* 				       itfc[k  ][j  ][i-1].z + */
/* 				       itfc[k  ][j-1][i  ].z + */
/* 				       itfc[k-1][j  ][i  ].z + */
/* 				       itfc[k  ][j-1][i-1].z + */
/* 				       itfc[k-1][j  ][i-1].z + */
/* 				       itfc[k-1][j-1][i  ].z + */
/* 				       itfc[k-1][j-1][i-1].z); */

/* 	    if (i==1 || j==1 || i==mx-2 || j==my-2) { */
/* 	      ucat[k][j][i].x *= 2.; */
/* 	      ucat[k][j][i].y *= 2.; */
/* 	      ucat[k][j][i].z *= 2.; */
/* 	    } */

/* 	    ucx = (ucat[k][j][i].x + ucat[k+1][j][i].x) * 0.5; */
/* 	    ucy = (ucat[k][j][i].y + ucat[k+1][j][i].y) * 0.5; */
/* 	    ucz = (ucat[k][j][i].z + ucat[k+1][j][i].z) * 0.5; */

/* 	    ucont[k][j][i].z = ucx * kzet[k][j][i].x + */
/* 	      ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z; */
	  }
	}
      }

    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (j>0 && i>0) {
	    ucont[k][j][i].z = (0.25 * (itfc[k  ][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j-1][i  ].x +
					itfc[k  ][j-1][i-1].x) *
				kzet[k][j][i].x +
				0.25 * (itfc[k  ][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j-1][i  ].y +
					itfc[k  ][j-1][i-1].y) *
				kzet[k][j][i].y +
				0.25 * (itfc[k  ][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j-1][i  ].z +
					itfc[k  ][j-1][i-1].z) *
				kzet[k][j][i].z);
/* 	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x + */
/* 				       itfc[k  ][j  ][i-1].x + */
/* 				       itfc[k  ][j-1][i  ].x + */
/* 				       itfc[k-1][j  ][i  ].x + */
/* 				       itfc[k  ][j-1][i-1].x + */
/* 				       itfc[k-1][j  ][i-1].x + */
/* 				       itfc[k-1][j-1][i  ].x + */
/* 				       itfc[k-1][j-1][i-1].x); */
	
/* 	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y + */
/* 				       itfc[k  ][j  ][i-1].y + */
/* 				       itfc[k  ][j-1][i  ].y + */
/* 				       itfc[k-1][j  ][i  ].y + */
/* 				       itfc[k  ][j-1][i-1].y + */
/* 				       itfc[k-1][j  ][i-1].y + */
/* 				       itfc[k-1][j-1][i  ].y + */
/* 				       itfc[k-1][j-1][i-1].y); */
	
/* 	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z + */
/* 				       itfc[k  ][j  ][i-1].z + */
/* 				       itfc[k  ][j-1][i  ].z + */
/* 				       itfc[k-1][j  ][i  ].z + */
/* 				       itfc[k  ][j-1][i-1].z + */
/* 				       itfc[k-1][j  ][i-1].z + */
/* 				       itfc[k-1][j-1][i  ].z + */
/* 				       itfc[k-1][j-1][i-1].z); */

/* 	    if (i==1 || j==1 || i==mx-2 || j==my-2) { */
/* 	      ucat[k][j][i].x *= 2.; */
/* 	      ucat[k][j][i].y *= 2.; */
/* 	      ucat[k][j][i].z *= 2.; */
/* 	    } */


/* 	    ucx = (ucat[k][j][i].x + ucat[k-1][j][i].x) * 0.5; */
/* 	    ucy = (ucat[k][j][i].y + ucat[k-1][j][i].y) * 0.5; */
/* 	    ucz = (ucat[k][j][i].z + ucat[k-1][j][i].z) * 0.5; */

/* 	    ucont[k-1][j][i].z = ucx * kzet[k-1][j][i].x + */
/* 	      ucy * kzet[k-1][j][i].y + ucz * kzet[k-1][j][i].z; */
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	
	    ucont[k][j][i].y = (0.25 * (itfc[k-1][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j  ][i  ].x +
					itfc[k-1][j  ][i-1].x) *
				jeta[k][j][i].x +
				0.25 * (itfc[k-1][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j  ][i  ].y +
					itfc[k-1][j  ][i-1].y) *
				jeta[k][j][i].y +
				0.25 * (itfc[k-1][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j  ][i  ].z +
					itfc[k-1][j  ][i-1].z) *
				jeta[k][j][i].z);
				
	
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	
	    ucont[k][j][i].y = (0.25 * (itfc[k-1][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j  ][i  ].x +
					itfc[k-1][j  ][i-1].x) *
				jeta[k][j][i].x +
				0.25 * (itfc[k-1][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j  ][i  ].y +
					itfc[k-1][j  ][i-1].y) *
				jeta[k][j][i].y +
				0.25 * (itfc[k-1][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j  ][i  ].z +
					itfc[k-1][j  ][i-1].z) *
				jeta[k][j][i].z);

	
	  }
	}
      }
    }



    DAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
    DAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
    DAVecRestoreArray(user[bi].fda, user[bi].lKZet, &kzet);
    DAVecRestoreArray(user[bi].fda, user[bi].lJEta, &jeta);
    
  }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(user[bi].nhostU);
    Contra2Cart(&(user[bi]));
  }

  
	return 0;
}

PetscErrorCode Block_Interface_P(UserCtx *user) {
  PetscInt bi;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  //  VecScatter tolocal;
  Vec	hostP;
  Vec	nhostP;
  PetscReal ***itfcp, *hostp, ***phi;

  VecCreateSeq(PETSC_COMM_SELF, 8, &hostP);

  
  for (bi=0; bi<block_number; bi++) {
    VecCreateSeq(PETSC_COMM_SELF,
		 user[bi].info.mx*user[bi].info.my*user[bi].info.mz,
		 &user[bi].nhostP);
/*     PetscInt N; */
/*     VecGetSize(user[bi].nhostP, &N); */
/*     PetscPrintf(PETSC_COMM_SELF, "Number %i\n", N); */
  }

  // First calculate Phi components at grid nodes
  for (bi=0; bi<block_number; bi++) {
    DALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    DAVecGetArray(user[bi].da, user[bi].ItfcP, &itfcp);
    DAVecGetArray(user[bi].da, user[bi].lPhi, &phi);

    DACreateNaturalVector(user[bi].da, &nhostP);

    PetscBarrier(PETSC_NULL);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (k<mz-1 && j<my-1 && i<mx-1) {
	    itfcp[k][j][i] = 0.125 * (phi[k  ][j  ][i  ] +
				      phi[k  ][j  ][i+1] +
				      phi[k  ][j+1][i  ] +
				      phi[k  ][j+1][i+1] +
				      phi[k+1][j  ][i  ] +
				      phi[k+1][j  ][i+1] +
				      phi[k+1][j+1][i  ] +
				      phi[k+1][j+1][i+1]);
	  }
	}
      }
    }
    if (bi==0 && zs!=0) {
      PetscPrintf(PETSC_COMM_SELF, "%le PPP\n", itfcp[28][21][21]);
    }
    DAVecRestoreArray(user[bi].da, user[bi].ItfcP, &itfcp);
    DAVecRestoreArray(user[bi].da, user[bi].lPhi, &phi);

    
    PetscBarrier(PETSC_NULL);

    DAGlobalToNaturalBegin(user[bi].da, user[bi].ItfcP, INSERT_VALUES, nhostP);
    DAGlobalToNaturalEnd(user[bi].da, user[bi].ItfcP, INSERT_VALUES, nhostP);
    
/*     DALocalToGlobal(user[bi].da, user[bi].lItfcP, INSERT_VALUES, user[bi].ItfcP); */
    VecScatter tolocalall;
    VecScatterCreateToAll(nhostP, &tolocalall, &(user[bi].nhostP));
/*     DAGlobalToNaturalAllCreate(user[bi].da, &tolocalall); */
    VecScatterBegin(tolocalall, nhostP, user[bi].nhostP, INSERT_VALUES, SCATTER_FORWARD);
    
    VecScatterEnd(tolocalall, nhostP, user[bi].nhostP, INSERT_VALUES,  SCATTER_FORWARD);
    VecScatterDestroy(tolocalall);

    VecGetArray(user[bi].nhostP, &hostp);
    if (bi==0) {
      PetscPrintf(PETSC_COMM_SELF, "%i %le PPP\n", zs, hostp[28*42*42+21*42+21]);
    }

    VecRestoreArray(user[bi].nhostP, &hostp);
    VecDestroy(nhostP);
  }

  for (bi=0; bi<block_number; bi++) {
    DALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    

    PetscInt    lmx, lmy, lmz;

    DAVecGetArray(user[bi].da, user[bi].lItfcP, &itfcp);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];
      
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze) {
	
	VecGetArray(user[hb].nhostP, &hostp);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfcp[ck][cj][ci] = (hostp[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))]
			     * (1-x) * (1-y) * (1-z) + //i,j,k
			     hostp[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))]
			     * x     * (1-y) * (1-z) + //i+1,j,k
			     hostp[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))]
			     * (1-x) * y     * (1-z) + //i,j+1,k
			     hostp[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))]
			     * x     * y     * (1-z) + //i+1,j+1,k
			     hostp[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))]
			     * (1-x) * (1-y) * z     + //i,j,k+1
			     hostp[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))]
			     * x     * (1-y) * z     + //i+1,j,k+1
			     hostp[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))]
			     * (1-x) * y     * z     + //i,j+1,k+1
			     hostp[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))]
			     * x     * y     * z); //i+1,j+1,k+1


	VecRestoreArray(user[hb].nhostP, &hostp);
      }
    }
    DAVecRestoreArray(user[bi].da, user[bi].lItfcP, &itfcp);
    DALocalToLocalBegin(user[bi].da, user[bi].lItfcP, INSERT_VALUES,
			user[bi].lItfcP);
    DALocalToLocalEnd(user[bi].da, user[bi].lItfcP, INSERT_VALUES,
		      user[bi].lItfcP);

  }
  VecDestroy(hostP);
  
  PetscReal ***p;

  for (bi=0; bi<block_number; bi++) {
    DALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    Cmpnts ***kcsi, ***keta, ***kzet, ***ucont;
    PetscReal ***kaj;
    DAVecGetArray(user[bi].da, user[bi].Phi, &p);
    DAVecGetArray(user[bi].da, user[bi].lItfcP, &itfcp);
    
    DAVecGetArray(user[bi].fda, user[bi].lKCsi, &kcsi);
    DAVecGetArray(user[bi].fda, user[bi].lKEta, &keta);
    DAVecGetArray(user[bi].fda, user[bi].lKZet, &kzet);
    DAVecGetArray(user[bi].da, user[bi].lKAj, &kaj);
    DAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	
	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=1;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[4] == 0 && zs==0) {
      k=1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
      if (ti>0) {
	PetscPrintf(PETSC_COMM_WORLD, "PP %le\n", p[k][21][21]);
      }
/*       k=1; */
/*       for (j=ys+1; j<ye-1; j++) { */
/* 	for (i=xs+1; i<xe-1; i++) { */
/* 	  ucont[k][j][i].z -= ((p[k+1][j][i]-p[k][j][i]) * */
/* 	    		 (kzet[k][j][i].x * kzet[k][j][i].x + */
/* 			  kzet[k][j][i].y * kzet[k][j][i].y + */
/* 			  kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */

/* 	  ucont[k][j][i].z -= */
/* 	    (0.25 * (p[k][j][i+1] + p[k+1][j][i+1] - */
/* 		     p[k][j][i-1] - p[k+1][j][i-1]) * */
/* 	     (kcsi[k][j][i].x * kzet[k][j][i].x + */
/* 	      kcsi[k][j][i].y * kzet[k][j][i].y + */
/* 	      kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	  ucont[k][j][i].z -=	       */
/* 	    (0.25 * (p[k][j+1][i] + p[k+1][j+1][i] - */
/* 		     p[k][j-1][i] - p[k+1][j-1][i]) * */
/* 	     (keta[k][j][i].x * kzet[k][j][i].x + */
/* 	      keta[k][j][i].y * kzet[k][j][i].y + */
/* 	      keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	} */
/*       } */

    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
      if (ti>0) {
	PetscPrintf(PETSC_COMM_SELF, "PP0 %le\n", p[k][21][21]);
      }
/*       k=mz-3; */
/*       for (j=ys+1; j<ye-1; j++) { */
/* 	for (i=xs+1; i<xe-1; i++) { */
/* 	  ucont[k][j][i].z -= ((p[k+1][j][i]-p[k][j][i]) * */
/* 	    		 (kzet[k][j][i].x * kzet[k][j][i].x + */
/* 			  kzet[k][j][i].y * kzet[k][j][i].y + */
/* 			  kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */

/* 	  ucont[k][j][i].z -= */
/* 	    (0.25 * (p[k][j][i+1] + p[k+1][j][i+1] - */
/* 		     p[k][j][i-1] - p[k+1][j][i-1]) * */
/* 	     (kcsi[k][j][i].x * kzet[k][j][i].x + */
/* 	      kcsi[k][j][i].y * kzet[k][j][i].y + */
/* 	      kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	  ucont[k][j][i].z -=	       */
/* 	    (0.25 * (p[k][j+1][i] + p[k+1][j+1][i] - */
/* 		     p[k][j-1][i] - p[k+1][j-1][i]) * */
/* 	     (keta[k][j][i].x * kzet[k][j][i].x + */
/* 	      keta[k][j][i].y * kzet[k][j][i].y + */
/* 	      keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	} */
/*       } */
    }
    DAVecRestoreArray(user[bi].da, user[bi].lItfcP, &itfcp);
    DAVecRestoreArray(user[bi].da, user[bi].Phi, &p);
    DAVecRestoreArray(user[bi].fda, user[bi].lKCsi, &kcsi);
    DAVecRestoreArray(user[bi].fda, user[bi].lKEta, &keta);
    DAVecRestoreArray(user[bi].fda, user[bi].lKZet, &kzet);
    DAVecRestoreArray(user[bi].da, user[bi].lKAj, &kaj);
    DAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
  }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(user[bi].nhostP);
  }
	return 0;
}


/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   5:	Inlet
   4:	Outlet
   8:   Characteristic BC
*/
int outflow_scale=1;

PetscErrorCode FormBCS(UserCtx *user, FSInfo *fsi)
{
  DA da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  Vec		Coor;
  Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet;
  PetscScalar	FluxIn, FluxOut, FluxOut2, ratio;
  PetscScalar   lArea, lArea2, AreaSum,AreaSum2, ***level;
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
  PetscInt      moveframe=1;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


// seokkoo
	extern PetscErrorCode CalcVolumeFlux(UserCtx *user, Vec lUcor, PetscReal *ibm_Flux, PetscReal *ibm_Area);
	double ibm_Flux=0, ibm_Area=0;
	//CalcVolumeFlux(user, user->Ucont, &ibm_Flux, &ibm_Area);
	ibm_Flux=0, ibm_Area=0;
    
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
/*   DAVecGetArray(fda, user->Ucont, &ucont); */
  DAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
	if(levelset) DAVecGetArray(da, user->lLevelset,  &level);

  DAVecGetArray(fda, user->lCsi,  &csi);
  DAVecGetArray(fda, user->lEta,  &eta);
  DAVecGetArray(fda, user->lZet,  &zet);


  //PetscInt ttemp;
  /*for (ttemp=0; ttemp<5; ttemp++) */{	// 5 times? why. seokkoo
  Contra2Cart(user);
  DAVecGetArray(fda, user->Ucat,  &ucat);

/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
  DAVecGetArray(fda, user->lUcont, &ucont);


  if (user->bctype[5]==6) {
    if (moveframe) {
      V_frame.x= -(fsi->S_new[1]-fsi->S_old[1]);
      V_frame.y= -(fsi->S_new[3]-fsi->S_old[3]);
      V_frame.z= -(fsi->S_new[5]-fsi->S_old[5]);
    } else {
      V_frame.x=0.;
      V_frame.y=0.;
      V_frame.z=0.;
    }

    if (moveframe) {
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (k>1 && k<ze-1)
	    ucont[k-1][j][i].z += V_frame.z * 0.5*(zet[k-1][j  ][i  ].z + zet[k][j][i].z) +
	                          V_frame.y * 0.5*(zet[k-1][j  ][i  ].y + zet[k][j][i].y) +
	                          V_frame.x * 0.5*(zet[k-1][j  ][i  ].x + zet[k][j][i].x)  ;
	    if (j>1 && j<ye-1)
	    ucont[k][j-1][i].y += V_frame.z * eta[k  ][j-1][i  ].z +
	                          V_frame.y * eta[k  ][j-1][i  ].y +
	                          V_frame.x * eta[k  ][j-1][i  ].x;
	    if (i>1 && i<xe-1)
	    ucont[k][j][i-1].x += V_frame.z * csi[k  ][j  ][i-1].z +
	                          V_frame.y * csi[k  ][j  ][i-1].y +
	                          V_frame.x * csi[k  ][j  ][i-1].x;
	  }
	}
      }
      
    }
        
  }

  DAVecRestoreArray(fda, user->lUcont, &ucont);

  DAVecGetArray(fda, user->Ucont, &ucont);

  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;	
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i+1].x;
	  FarFluxIn += ucont[k][j][i].x;
	  FarAreaIn += csi[k][j][i].x;
	}
      }
    }
  }
  
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
/* 	  FarFluxIn -= ucont[k][j][i-1].x; */
	  FarFluxOut += ucont[k][j][i-1].x;
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j+1][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
/* 	  FarFluxIn -= ucont[k][j-1][i].y; */
	  FarFluxOut += ucont[k][j-1][i].y;
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k+1][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
/* 	  FarFluxIn -= ucont[k-1][j][i].z; */
	  FarFluxOut += ucont[k-1][j][i].z;
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }

  PetscGlobalSum(&FarFluxIn, &FarFluxInSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&FarFluxOut, &FarFluxOutSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&FarAreaIn, &FarAreaInSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&FarAreaOut, &FarAreaOutSum, PETSC_COMM_WORLD);

  if (user->bctype[5]==6) {
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ;
    VelDiffIn  = FluxDiff / FarAreaInSum ;
    if (fabs(FluxDiff) < 1.e-6) VelDiffIn = 0.;
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = FluxDiff / FarAreaOutSum ;
    if (fabs(FluxDiff) < 1.e-6) VelDiffOut = 0.;
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;
    if (moveframe) {
      V_frame.x= -(fsi->S_new[1]-fsi->S_old[1]);
      V_frame.y= -(fsi->S_new[3]-fsi->S_old[3]);
      V_frame.z= -(fsi->S_new[5]-fsi->S_old[5]);
    } else {
      V_frame.x=0.;
      V_frame.y=0.;
      V_frame.z=0.;
    }

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
    PetscPrintf(PETSC_COMM_WORLD, "Cop Vel  Diff %d %le %le %le \n", ti, V_frame.x,V_frame.y,V_frame.z);    
        
  }


  // scale global mass conservation

  if (user->bctype[5]==6) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut + V_frame.z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut + V_frame.y;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut + V_frame.x;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn + V_frame.x;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i+1].x;
	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn + V_frame.y;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j+1][i].y;
	}
      }
    }
  }
  

  if (user->bctype[4]==6) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn + V_frame.z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k+1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = FluxInSum + FarFluxInSum;
    PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / FluxOutSum;
    if (fabs(FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    PetscPrintf(PETSC_COMM_WORLD, "Char Ratio %d %le %le %le %le %d %d\n", ti, ratio, FluxIn, FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (ti==0 || ti==1) 
	    if (inletprofile<0) 
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6) 
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else 
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
/*     OUTLET BC :4 */
/* ==================================================================================             */
/*
// temp for plate !!
   for (k=lzs; k<lze; k++) { 
     for (j=lys; j<lye; j++) { 
       for (i=lxs; i<lxe; i++) { 
 	if (k==121 && j>96 && j<145)
 	  ucont[k][j][i].z = 0.;
       } 
     } 
   } 
*/
	double FluxOut_gas=0, lArea_gas=0, AreaSum_gas=0, ratio_gas;
  
	if(user->bctype[0]==11) {
		PetscReal	***nvert;
		DAVecGetArray(da, user->lNvert, &nvert);
		
		lArea=0.;
		if (xs==0) {
			i = 0;
			/*
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
				if( zc > 0 ) ubcs[k][j][i] = ucat[k][j][i+1];
			}
			*/
			
			FluxOut = 0;
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
				if( zc > 0 && nvert[k][j][i+1] < threshold) {
					
					double u=ucat[k][j][i+1].x, v=ucat[k][j][i+1].y, w=ucat[k][j][i+1].z;
					ucat[k][j][i].x=u;
					ucat[k][j][i].y=v;
					ucat[k][j][i].z=w;
					
					ucont[k][j][i].x = u*csi[k][j][i].x + v*csi[k][j][i].y + w*csi[k][j][i].z;
					
					FluxOut +=  ucont[k][j][i].x;
					//lArea += sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
					lArea += fabs(csi[k][j][i].z);
				}
			}
		}
		else FluxOut = 0.;
		
		FluxIn = FluxInSum + FarFluxInSum;
		PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
		PetscGlobalSum(&lArea, &AreaSum, PETSC_COMM_WORLD);
		 
		//here
		FluxOutSum *= -1;
		ratio = (FluxInSum - FluxOutSum) / AreaSum;
		
		double FluxOut_new=0, FluxOut_new_sum;
		
		if(outflow_scale) {
			PetscPrintf(PETSC_COMM_WORLD, "Time %d, Vel correction=%e, FluxIn=%e, FluxOut=%e, Area=%f\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum);
		
			if (xs==0) {
				i=0;
				for (j=lys; j<lye; j++) 
				for (k=lzs; k<lze; k++) {
					double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
					if( zc > 0 && nvert[k][j][i+1] < threshold) {
						double Area = sqrt( csi[k][j][i+1].x*csi[k][j][i+1].x + csi[k][j][i+1].y*csi[k][j][i+1].y + csi[k][j][i+1].z*csi[k][j][i+1].z );
						Area = csi[k][j][i+1].z;
						ucont[k][j][i].x += (FluxInSum - FluxOutSum) * Area / AreaSum;
						FluxOut_new += ucont[k][j][i].x;
					}
				}
			}
		}		
		PetscGlobalSum(&FluxOut_new, &FluxOut_new_sum, PETSC_COMM_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "Corrected FluxOut=%e\n", FluxOut_new_sum);
		DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	}
		
	if (user->bctype[5]==4) {
		PetscReal	***nvert;	//seokkoo
		DAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
		
		lArea=0.;
		if (ze==mz) {
			k = ze-1;
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {  
				ubcs[k][j][i].x = ucat[k-1][j][i].x;
				ubcs[k][j][i].y = ucat[k-1][j][i].y;
				if (nvert[k-1][j][i] < threshold) ubcs[k][j][i].z = ucat[k-1][j][i].z;
				else ubcs[k][j][i].z = ucat[k-1][j][i].z * 0;
			}
						
			FluxOut = 0;
			FluxOut2= 0;			
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {
				if (nvert[k-1][j][i] < threshold) //seokkoo
				{/*
					if ( levelset && ! (level[k-1][j][i] > 0 || ( level[k-1][j][i] <0 && level[k-1][j-1][i] > 0 ) ) ) { 
						FluxOut_gas +=  ucont[k-1][j][i].z;
						lArea_gas += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
					}
					else */{
					  //if( levelset && (level[k-1][j][i] > 0 || ( level[k-1][j][i] <0 && level[k-1][j-1][i] > 0 ) ) ) 
					  {//water part
					    //if(ucont[k-1][j][i].z<0) ucont[k-1][j][i].z=0; // 100216 fix levelset outlet problem
					  }          
					  //haha
					  // ucont[k-1][j][i].z = ucat[k-1][j][i].x*zet[k-1][j][i].x + ucat[k-1][j][i].y*zet[k-1][j][i].y + ucat[k-1][j][i].z*zet[k-1][j][i].z;
						if (levelset && air_flow_levelset){
							if(level[k-1][j][i]<=-dthick){
								FluxOut += ucont[k-1][j][i].z ;							
								lArea += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
							}
							else{
								FluxOut2 +=  ucont[k-1][j][i].z;
								lArea2 += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
							}
						}
						else{
							FluxOut +=  ucont[k-1][j][i].z;					
							lArea += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
						}
					}
				}
			}
		}
		else{	
			FluxOut = 0.;
			FluxOut2 = 0.;		
		}
		
		FluxIn = FluxInSum + FarFluxInSum;
		PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
		PetscGlobalSum(&FluxOut_gas, &FluxOutSum_gas, PETSC_COMM_WORLD);
		
		PetscGlobalSum(&lArea, &AreaSum, PETSC_COMM_WORLD);
		PetscGlobalSum(&lArea_gas, &AreaSum_gas, PETSC_COMM_WORLD);
		 
		ratio = (FluxInSum - FluxOutSum) / AreaSum;
		ratio_gas = (FluxInSum_gas - FluxOutSum_gas) / AreaSum_gas;
		
		/**
		FluxInSum += FluxInSum_gas;
		FluxOutSum += FluxOutSum_gas;
		AreaSum += AreaSum_gas;
		
		FluxInSum_gas = FluxInSum;
		FluxOutSum_gas = FluxOutSum;
		AreaSum_gas = AreaSum;
		
		ratio = ratio_gas = (FluxInSum - FluxOutSum) / AreaSum;
		**/
		
		if(outflow_scale) {
			PetscPrintf(PETSC_COMM_WORLD, "Time %d, Vel correction=%e, FluxIn=%e, FluxOut=%e, Area=%f\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum);
			//PetscPrintf(PETSC_COMM_WORLD, "Vel correction=%e, FluxIn=%e, FluxOut=%e, Area=%f\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum);
		
			if (ze==mz) {
				k = ze-1;
				for (j=lys; j<lye; j++) 
				for (i=lxs; i<lxe; i++) {
					double Area = sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
					
					if (nvert[k-1][j][i] < threshold) {
						/*if( levelset && ! (level[k-1][j][i] > 0 || ( level[k-1][j][i] <0 && level[k-1][j-1][i] > 0 ) ) ) {
							ucont[k-1][j][i].z += (FluxInSum_gas - FluxOutSum_gas) * Area / AreaSum_gas;
							ucont[k-1][j][i].z = 0;
							//ucat[k][j][i] = ucat[k-1][j][i];
						}
						else */{
							if (levelset && user->bctype[5]==4 && air_flow_levelset) {
								if(level[k-1][j][i]<=-dthick){
									ucont[k-1][j][i].z += (FluxInSum - FluxOutSum) * Area / (AreaSum);
								}
								else{
									//ucont[k-1][j][i].z += 0.;//(FluxInSum2 - FluxOutSum2) * Area / AreaSum2;
								}
							}
							else ucont[k-1][j][i].z += (FluxInSum - FluxOutSum) * Area / AreaSum;
							//ucat[k][j][i] = ucat[k-1][j][i];
							//ubcs[k][j][i].z += ratio;
							/*
							Cmpnts u, ucon=ucont[k-1][j][i];
							ucon.x = ucon.y = 0;
							Contra2Cart_single(csi[k-1][j][i], eta[k-1][j][i], zet[k-1][j][i], ucon, &u);
							
							ucat[k][j][i].x = - ucat[k-1][j][i].x + 2*u.x;
							ucat[k][j][i].y = - ucat[k-1][j][i].y + 2*u.y;
							ucat[k][j][i].z = - ucat[k-1][j][i].z + 2*u.z;
							*/
						}
					}
				}
			}
		}		
		
		DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	}
	
  else if (user->bctype[5]==0) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	}
      }
    }
  } 
  else if (user->bctype[5]==2) {
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 1.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 0.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
/*
  if (user->bctype[3]==2) {
    // Designed for driven cavity problem 
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
        for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = 1.;// - ucat[k-1][j][i].x;
          ubcs[k][j][i].y = 0.;//- ucat[k-1][j][i].y;
          ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
        }
      }
    }
  }
*/
  
/*   OUTLET at k==0 */
  if (user->bctype[4]==4) {
    lArea=0.;
	PetscReal	***nvert;	//seokkoo
	DAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
    if (zs == 0) {
      k = zs;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j][i] < threshold) //seokkoo 
		{
			FluxOut += (ucat[k+1][j][i].x * (zet[k][j][i].x + zet[k+1][j][i].x) +
				ucat[k+1][j][i].y * (zet[k][j][i].y + zet[k+1][j][i].y) +
				ucat[k+1][j][i].z * (zet[k][j][i].z + zet[k+1][j][i].z)) * 0.5;

			lArea += zet[k][j][i].z;
		}
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    DAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
    
    FluxIn = FluxInSum + FarFluxInSum;
    PetscGlobalSum(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lArea, &AreaSum, PETSC_COMM_WORLD);

    ratio = (FluxInSum - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Ratio b %d  %le %le %le %d %d\n", ti, ratio, FluxInSum, FluxOutSum, zs, ze);

    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  //  PetscPrintf(PETSC_COMM_WORLD, "Ratio b %d %d %d\n", i,j,k);

	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z + ratio;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;
	}
      }
    }
  }
  
  // slip
	Cmpnts ***lucont;	// for use of ucont[k-1] etc..
	DAVecGetArray(fda, user->lUcont, &lucont);
  
  
	if ( (user->bctype[0]==1 || user->bctype[0]==10) && xs==0) {
		i= 0;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (j=lys; j<lye; j++) {
			ucont[k][j][i].x=0;
		}
	}
      
	if ( (user->bctype[1]==1 || user->bctype[1]==10) && xe==mx) {
		i= xe-1;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			ucont[k][j][i-1].x=0;
		}
	}
	
	if ( (user->bctype[2]==1 || user->bctype[2]==10) && ys==0) {
		j= 0;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (i=lxs; i<lxe; i++) {
			ucont[k][j][i].y=0;
		}
	}
	
	if ( (user->bctype[3]==1 || user->bctype[3]==10) && ye==my) {
		j= ye-1;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (i=lxs; i<lxe; i++) {
			ucont[k][j-1][i].y=0;
		}
	}
	/*
	if ( (user->bctype[4]==1 || user->bctype[4]==10) && zs==0) {
		k= 0;
		for (i=lxs; i<lxe; i++)
		for (j=lys; j<lye; j++) {
			ucont[k][j][i].z=0;
		}
	}
	
	if ( (user->bctype[5]==1 || user->bctype[5]==10) && ze==mz) {
		k= ze-1;
		for (i=lxs; i<lxe; i++)
		for (j=lys; j<lye; j++) {
			ucont[k-1][j][i].z=0;
		}
	}*/
    
  DAVecRestoreArray(fda, user->lUcont, &lucont);
//  end slip
  
  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    
  DAVecRestoreArray(fda, user->Ucat, &ucat);
  
  Contra2Cart(user);
  DAVecGetArray(fda, user->Ucat, &ucat);
  



/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */
  if (user->bctype[0]==3) {
	  
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }

/* ==================================================================================             */
/*   INTERFACE BC */
/* ==================================================================================             */
  if (user->bctype[0]==0) {
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i+1].x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==0) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i-1].y;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }

  if (user->bctype[2]==0) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = ucat[k][j+1][i].y;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==0) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = ucat[k][j-1][i].y;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }

  if (user->bctype[4]==0) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==0) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
  // boundary conditions on ghost nodes
  
  //removed by seokkoo
	  /*
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }
*/

  
  // 0 velocity on the corner point
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }
  DAVecRestoreArray(fda, user->Ucat,  &ucat);
  }

  DAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DAVecRestoreArray(fda, Coor, &coor);

  DAVecRestoreArray(fda, user->lCsi,  &csi);
  DAVecRestoreArray(fda, user->lEta,  &eta);
  DAVecRestoreArray(fda, user->lZet,  &zet);
  
  if(levelset) DAVecRestoreArray(da, user->lLevelset,  &level);

  //  DAVecRestoreArray(fda, user->Ucont_o, &ucont_o);

  DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  return(0);
}

PetscErrorCode fluxin(UserCtx *user)
{
  PetscInt  iRotate;

  PetscInt ts_p_cycle;
  PetscInt opening, closing;
  PetscInt open_steps, close_steps;

/*   ts_p_cycle = 500; */
/*   opening = 0; */
/*   open_steps = 50; */

/*   closing = 225; */
/*   close_steps = 50; */

  //ts_p_cycle = 1290;
  ts_p_cycle = 2500;//10000;//5000;

  opening = 10;
  open_steps = 100;

  closing = 580;
  close_steps = 80;

  PetscReal t_rel;

  iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);

  // if (angle>.0 && iRotate>1058) iRotate-=angle;

  t_rel = iRotate * (1. / ts_p_cycle) * 860 + 6.8  ; //+7.15;
/*   t_rel = (iRotate-940) * (1. / ts_p_cycle) * 860/2 + */
/*                    940. * (1. / 2500.     ) * 860 + 6.8;     */

  PetscInt i;
  PetscTruth interpolated = PETSC_FALSE;
  //PetscPrintf(PETSC_COMM_WORLD, "Inflow00 Rate %d %e %e %i\n",ti, Flux_in, t_rel, user->number_flowwave);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    for (i=0; i<user->number_flowwave-1; i++) {
/*       PetscPrintf(PETSC_COMM_WORLD, "Inflow Rate %e %e\n", Flux_in, user->inflow[i].t); */
      if (t_rel >= user->inflow[i].t && t_rel <= user->inflow[i+1].t) {
	Flux_in = user->inflow[i].f + (user->inflow[i+1].f - user->inflow[i].f) /
	  (user->inflow[i+1].t - user->inflow[i].t) *
	  (t_rel - user->inflow[i].t);
/* 	PetscPrintf(PETSC_COMM_SELF, "Inflow Rate %i %e %e %e %e %e %e\n", i, Flux_in, t_rel, user->inflow[i].f, user->inflow[i].t, user->inflow[i+1].f, user->inflow[i+1].t); */
	interpolated = PETSC_TRUE;
      }
      if (interpolated) break;
    }  
    //if (t_rel > 350) Flux_in = 0.;

    MPI_Bcast(&Flux_in, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else {
    MPI_Bcast(&Flux_in, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

/*   if (Flux_in<0.) { */
/*     user->bctype[5]= 5; */
/*     user->bctype[4]= 4; */
/*     PetscPrintf(PETSC_COMM_WORLD, "IINNFLOW change inlet!!!!!!!%e\n", Flux_in); */
/*   } else { */
/*     user->bctype[5]= 4; */
/*     user->bctype[4]= 5; */
/*   } */
    

/*   PetscPrintf(PETSC_COMM_SELF, "IINNFLOW %e\n", Flux_in); */
/*   Flux_in = PetscMax(0.000, sin( (iRotate) * (1./ts_p_cycle) * 2 * 3.1415926)); */

  /*
  if (iRotate >= opening && iRotate<=opening + open_steps) {
    angle = -rg + rg/(PetscReal)open_steps * (iRotate);
  }
  else if (iRotate>(closing - close_steps) && iRotate <=closing) {
    angle = -rg/(PetscReal)close_steps  * (iRotate-(closing-close_steps));
  }
  else if (iRotate>closing) {
    angle = -rg;
  }
  else {
    angle = 0;
  }
  */

/*   angle = 0.; */
  //Flux_in = 1.;
  PetscPrintf(PETSC_COMM_WORLD, "Angle %d %le %le flux-in %le intp%d\n",ti, t_rel, angle, Flux_in, interpolated);

  return 0;
}

PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont)
{
  DA fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	lFluxOut = 0., ratio;
  Cmpnts ***ucont, ***zet, ***ucat, ***ubcs;
  
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  if (user->bctype[5] == 4) {

    //    OutflowFlux(user);

    Contra2Cart(user);
    DAVecGetArray(fda, Ucont, &ucont);
    DAVecGetArray(fda, user->lZet, &zet);
    DAVecGetArray(fda, user->Ucat, &ucat);
    DAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
    /* Inflow flux at previous time step is 0, under this condition, it's assumed
       the flux difference at two time steps is uniformly distributed 
       to all outflow boundary nodes.*/
    if (ti==1) {//Flux_in_old < 1.e-6) {
      PetscReal lArea = 0., AreaSum;
      
      if (ze == mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lArea += zet[k][j][i].z;
	  }
	}
      }
      PetscGlobalSum(&lArea, &AreaSum, PETSC_COMM_WORLD);

      PetscReal vd;
      vd = (user->FluxInSum) / AreaSum;
      PetscPrintf(PETSC_COMM_SELF, "FluxOOOO %e %e %e\n", vd, Flux_in, Flux_in);
      
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].z = vd;
	    ucont[k-1][j][i].z = vd * zet[k-1][j][i].z;
	  }
	}
      }
    }
    /* Scale the outflow flux to ensure global flux conservation */
    else {
      lFluxOut = 0.;
      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lFluxOut += (ucat[k][j][i].x * (zet[k][j][i].x + zet[k-1][j][i].x) +
			 ucat[k][j][i].y * (zet[k][j][i].y + zet[k-1][j][i].y) +
			 ucat[k][j][i].z * (zet[k][j][i].z + zet[k-1][j][i].z))
	      * 0.5;
	  }
	}
      }
      PetscGlobalSum(&lFluxOut, &FluxOutSum, PETSC_COMM_WORLD);
      PetscBarrier(PETSC_NULL);
      ratio = user->FluxInSum / FluxOutSum;

      PetscPrintf(PETSC_COMM_WORLD, "Ratio %e %e\n", ratio, FluxOutSum);
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z * ratio;
	    ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  }
	}
      }
      
    }
    DAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
    DAVecRestoreArray(fda, Ucont, &ucont);
    DAVecRestoreArray(fda, user->lZet, &zet);
    DAVecRestoreArray(fda, user->Ucat, &ucat);

/*     DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*     DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

/*     Contra2Cart(user, user->lUcont, user->Ucat); */
  }
  return 0;
}

PetscErrorCode SetInitialGuessToOne(UserCtx *user)
{
	DA da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	Cmpnts ***ucont, ***cent;
	Cmpnts ***icsi, ***jeta, ***kzet, ***zet;
  
	PetscInt i, j, k;

	PetscReal	***nvert, ***p, ***level;	//seokkoo
	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Vec Coor;
	Cmpnts	***coor;
	DAGetGhostedCoordinates(da, &Coor);
	
	if(levelset) DAVecGetArray(da, user->lLevelset, &level);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->P, &p);
	DAVecGetArray(fda, Coor, &coor);
	DAVecGetArray(fda, user->Ucont, &ucont);
	DAVecGetArray(fda, user->lICsi,  &icsi);
	DAVecGetArray(fda, user->lJEta,  &jeta);
	DAVecGetArray(fda, user->lKZet,  &kzet);
	DAVecGetArray(fda, user->lZet,  &zet);
	DAVecGetArray(fda, user->lCent,  &cent);
  
	double lArea=0, SumArea;
	if(zs==0) {
		k=0;
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) if(nvert[k][j][i]+nvert[k+1][j][i]<0.1) lArea+=sqrt( kzet[k][j][i].x*kzet[k][j][i].x + kzet[k][j][i].y*kzet[k][j][i].y + kzet[k][j][i].z*kzet[k][j][i].z );
	}
	PetscGlobalSum(&lArea, &SumArea, PETSC_COMM_WORLD);
		      
	
	  
	
	double a = 2.*M_PI;
	double lambda = user->ren/2. - sqrt ( pow(user->ren/2., 2.0) + pow(a, 2.0) );
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=xs; i<lxe; i++) {
		double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	// centx[].x
		double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	// centx[].y
		double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

		if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].x = - cos (xi*a) * sin (yi*a) * icsi[k][j][i].x;
		}
		else if(inletprofile==18) {	// 2D Kovasznay flow
			ucont[k][j][i].x = ( 1.0 - exp ( lambda*xi )*cos( a*yi ) ) * icsi[k][j][i].x;
		}
		else if(inletprofile==20) {// Enright test
		  ucont[k][j][i].x = 2 * pow( sin (M_PI*xi), 2.) * sin (2.*M_PI*yi) * sin (2.*M_PI*zi) * icsi[k][j][i].x;
		}
		else if( inletprofile==0 ) ucont[k][j][i].x = 0;
	}
	
	for (k=lzs; k<lze; k++)
	for (j=ys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	
		double xj = (coor[k  ][j][i  ].x + coor[k-1][j][i  ].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;
		double yj = (coor[k  ][j][i  ].y + coor[k-1][j][i  ].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;
		double zj = (coor[k  ][j][i  ].z + coor[k-1][j][i  ].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25;

		if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].y = sin (xj*a) * cos (yj*a) * jeta[k][j][i].y;
		}
		else if(inletprofile==18) {	// 2D Kovasznay flow
			ucont[k][j][i].y = ( lambda/a * exp ( lambda*xj )*sin( a*yj ) ) * jeta[k][j][i].y;
		}
		else if(inletprofile==20) {// Enright test
		  ucont[k][j][i].y = - sin (2.*M_PI*xj) * pow(sin (M_PI*yj), 2.) * sin (2.*M_PI*zj) * jeta[k][j][i].y;
                }
		else if( inletprofile==0 ) ucont[k][j][i].y = 0;

	}
	
	for (k=zs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double xk = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
		double yk = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
		double zk = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;
		if(inletprofile==20) {// Enright test
			ucont[k][j][i].z = - sin (2.*M_PI*xk) * sin (2.*M_PI*yk) * pow( sin (M_PI*zk), 2.) * kzet[k][j][i].z;
		}
		else if( inletprofile==0 ) ucont[k][j][i].z = 0;
		else if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].z = 0;
			if(k) {
				p[k][j][i] = - 0.25 * ( cos(2.*a*cent[k][j][i].x) + cos(2.*a*cent[k][j][i].y) );
			}
		}
		else if(inletprofile==18) {	// 18: 2D Kovasznay flow
			ucont[k][j][i].z = 0;
		}
		else if(inletprofile==20) {}
		else if(nvert[k][j][i]+nvert[k+1][j][i]<0.1) {
			double w;
		      
			//double area = sqrt( kzet[k][j][i].x*kzet[k][j][i].x + kzet[k][j][i].y*kzet[k][j][i].y + kzet[k][j][i].z*kzet[k][j][i].z );
			double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double xc = (coor[k+1][j][i].x + coor[k+1][j-1][i].x + coor[k+1][j][i-1].x + coor[k+1][j-1][i-1].x) * 0.25;
			double yc = (coor[k+1][j][i].y + coor[k+1][j-1][i].y + coor[k+1][j][i-1].y + coor[k+1][j-1][i-1].y) * 0.25;
			double zc = (coor[k+1][j][i].z + coor[k+1][j-1][i].z + coor[k+1][j][i-1].z + coor[k+1][j-1][i-1].z) * 0.25;
		      
			if(inletprofile==10) {
				double delta = 0.45263, a=5.99;
				if( yc>=delta ) w=1.0;
				else if(yc<=0) w=0.0;
				else w = pow( yc/delta, 1./a );
			}
			else if(inletprofile==12) {	// pipe shear stress test
				double r = sqrt(xc * xc + yc * yc);
				w = 2*(1 - pow(r/0.5,2.0) );
			}
			else if(inletprofile==13) {	// periodic channel flow
				double w_bulk=inlet_flux/SumArea;
				w = 1.5 * w_bulk * (yc/channel_height*2.) * ( 2. - yc /channel_height*2.);
				if(air_flow_levelset_periodic && level[k][j][i]>=-dthick)w=0;
			}
			else if(inletprofile==101) {	// periodic channel flow
				double w_bulk=inlet_flux/SumArea;
				w = 1.5 * w_bulk * yc/1000. * ( 2 - yc/1000. );
				if(air_flow_levelset_periodic && level[k][j][i]>=-dthick)w=0;
			}			
			else if(inletprofile==15) w = 0;	// jet
			else if(inletprofile==16) {		// periodic pipe flow
				double r = sqrt(xc * xc + yc * yc);
				double w_bulk=inlet_flux/SumArea;
				w = w_bulk*2.*(1. - pow(r/0.5,2.0) );
			}
			else if(inletprofile==-9) w=0;
			else if(inlet_flux>0) {
				/*
				w=inlet_flux/user->k_area[k];*/
				w=inlet_flux/inletArea;
			}
			else w=1;

			if(levelset && air_flow_levelset && level[k][j][i]>=-dthick) w=0.;
			if(levelset && air_flow_levelset_periodic && level[k][j][i]>=-dthick) w=0.;
			
			ucont[k][j][i].z = w * area;
			//	printf("%f %f %f\n", w, inlet_flux,inletArea);
		       
			if(inletprofile==19) {
				double u=0, v=0, w=1;
				ucont[k][j][i].x = u * icsi[k][j][i].x  + v  * icsi[k][j][i].y + w * icsi[k][j][i].z;
				ucont[k][j][i].y = u * jeta[k][j][i].x  + v  * jeta[k][j][i].y + w * jeta[k][j][i].z;
				ucont[k][j][i].z = u * kzet[k][j][i].x  + v  * kzet[k][j][i].y + w * kzet[k][j][i].z;
			}
			if(air_flow_levelset==2){
				double x = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
				double y = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
				//printf("x coordinate: %f, y coordinate: %f\n",x,y);
				double u=0.;
				double v=0.;
				w=0.;
				if(level[k][j][i]<=-dthick){
					if(x>user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1])x=2*user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1]-x;
					if(x<user[0].wave_inf[0].WIND_Y[0])x=2*user[0].wave_inf[0].WIND_Y[0]-x;
					WIND_vel_interpolate(&user[0], &u,  &v,  &w, x, y, user[0].wave_inf[0].WIND_id_x[i-lxs][j-lys], user[0].wave_inf[0].WIND_id_y[i-lxs][j-lys]);	
					if(x>user[0].wave_inf[0].WIND_Y[user[0].wave_inf[0].NXMOD-1] || x<user[0].wave_inf[0].WIND_Y[0])u=-u;
				}
				else{
					u=0.,v=0.,w=0.;
				}
				u=0.,v=0.;
				//ucont[k][j][i].x = u * icsi[k][j][i].x  + v  * icsi[k][j][i].y + w * icsi[k][j][i].z;
				//ucont[k][j][i].y = u * jeta[k][j][i].x  + v  * jeta[k][j][i].y + w * jeta[k][j][i].z;
				ucont[k][j][i].z = u * kzet[k][j][i].x  + v  * kzet[k][j][i].y + w * kzet[k][j][i].z;				
			}				
			
			
		}
	}
    
	srand( time(NULL)) ;	// seokkoo
	for (i = 0; i < (rand() % 3000); i++) (rand() % 3000);	//seokkoo
  
	if( initial_perturbation) {
		PetscPrintf(PETSC_COMM_WORLD, "\nGenerating initial perturbation\n");
		for(k=lzs; k<lze; k++) {	// for all CPU
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i]+nvert[k+1][j][i] < threshold) {
				  int n1, n2, n3;
					double F;
					
					F  = 1.00; // 100%
					n1 = rand() % 20000 - 10000;
					n2 = rand() % 20000 - 10000;
					n3 = rand() % 20000 - 10000;
					/*
					if(inletprofile==13) {
						ucont[k][j][i].z += 10*((double)n)/10000.* F * kzet[k][j][i].z;
						
						F  = 0.1;
						n = rand() % 20000; n -= 10000;
						ucont[k][j][i].x = ((double)n)/10000. * F * icsi[k][j][i].x;
					
						F  = 0.1;
						n = rand() % 20000; n -= 10000;
						ucont[k][j][i].y = ((double)n)/10000. * F * jeta[k][j][i].y;
					}
					else */
					if(!levelset) {
						ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*F );
						ucont[k][j][i].x = ((double)n3)/10000. * 0.1 * ucont[k][j][i].z;
						ucont[k][j][i].y = ((double)n2)/10000. * 0.1 * ucont[k][j][i].z;
						ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*F );
					}
					else if (levelset && air_flow_levelset_periodic && level[k][j][i]<-dthick){
						//ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*F );
						double coeff=.8;
						ucont[k][j][i].x = ((double)n3)/10000.*coeff * 0.1 * ucont[k][j][i].z;
						ucont[k][j][i].y = ((double)n2)/10000.*coeff * 0.1 * ucont[k][j][i].z;
						ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*coeff*F );
					}
				}
			}
		}
	}
	
	if(levelset) DAVecRestoreArray(da, user->lLevelset, &level);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->P, &p);
	DAVecRestoreArray(fda, Coor, &coor);
	DAVecRestoreArray(fda, user->Ucont, &ucont);
	DAVecRestoreArray(fda, user->lICsi,  &icsi);
	DAVecRestoreArray(fda, user->lJEta,  &jeta);
	DAVecRestoreArray(fda, user->lKZet,  &kzet);
	DAVecRestoreArray(fda, user->lZet,  &zet);
	DAVecRestoreArray(fda, user->lCent,  &cent);
	
	DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	
	DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
	
	Contra2Cart(user);
	
	VecCopy(user->Ucont, user->Ucont_o);
        DAGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
        DAGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);

        DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
        DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);

        DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);
        DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);

	
	return 0;
}
