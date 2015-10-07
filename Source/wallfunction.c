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

double integrate_F(double nu, double utau, double yb);
double integrate_F(double nu, double utau, double yb, double ks);
double integrate_Fy(double nu, double utau, double ya, double yb);
double u_loglaw_roughness(double nu, double y, double utau, double ks);
double f_loglaw_roughness(double nu, double u, double y, double utau, double ks);
double df_loglaw_roughness(double nu, double u, double y, double utau, double ks);
double find_utau_loglaw_roughness(double nu, double u, double y, double guess, double ks);
double sign(double a);

double kappa=0.41, B=5.5;

void wall_function_freesurface (double nu, double sb, Cmpnts Ua, Cmpnts Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Ub.x - Ua.x, v_c = Ub.y - Ua.y, w_c = Ub.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	*ustar = find_utau_Cabot(nu, ut_mag, sb, 0.01, 0);
}

void wall_function (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	*ustar = find_utau_Cabot(nu, ut_mag, sc, 0.01, 0);
	//double ut_mag_modeled = u_Werner(1./user->ren, sb, utau);
	double ut_mag_modeled = u_Cabot(nu, sb, *ustar, 0);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void wall_function_roughness (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	*ustar = find_utau_Cabot_roughness(nu, ut_mag, sc, 0.01, 0, ks);
	double ut_mag_modeled = u_Cabot_roughness(nu, sb, *ustar, 0, ks);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void wall_function_roughness_loglaw (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	*ustar = find_utau_loglaw_roughness(nu, ut_mag, sc, 0.01, ks);
	double ut_mag_modeled = u_loglaw_roughness(nu, sb, *ustar, ks);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut;// + sb/sc * un * nx;
	(*Ub).y = vt;// + sb/sc * un * ny;
	(*Ub).z = wt;// + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void wall_function_slip (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
  double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
  double un = u_c * nx + v_c * ny + w_c * nz;
  double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
  double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
  *ustar = find_utau_Cabot(nu, ut_mag, sc, 0.01, 0);
  //double ut_mag_modeled = u_Werner(1./user->ren, sb, utau);
  double ut_mag_modeled = u_Cabot(nu, sb, *ustar, 0);

  if(ut_mag>1.e-10) {
    ut *= ut_mag_modeled/ut_mag;
    vt *= ut_mag_modeled/ut_mag;
    wt *= ut_mag_modeled/ut_mag;
  }
  else ut=vt=wt=0;

  // u = ut + (u.n)n
  (*Ub).x = ut;// + sb/sc * un * nx;
  (*Ub).y = vt;// + sb/sc * un * ny;
  (*Ub).z = wt;// + sb/sc * un * nz;

  (*Ub).x += Ua.x;
  (*Ub).y += Ua.y;
  (*Ub).z += Ua.z;
}

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	
	*ustar = sqrt ( ut_mag / sc / user->ren );
	
	(*Ub).x = sb/sc * u_c;
	(*Ub).y = sb/sc * v_c;
	(*Ub).z = sb/sc * w_c;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, double nx, double ny, double nz)

{
  
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
	
	/*
	//Iman
	double Uan = Ua.x * nx + Ua.y * ny + Ua.z * nz;
	double Ucn = Uc.x * nx + Uc.y * ny + Uc.z * nz;
	double Ubn = Uan + (Ucn - Uan) * sb/sc;
	//   in the Tangential direction Ub_t = Uc_t                                                                                                                                          
	double ut = Uc.x - Ucn * nx, vt = Uc.y - Ucn * ny, wt = Uc.z - Ucn * nz;

        (*Ub).x = ut;
        (*Ub).y = vt;
        (*Ub).z = wt;

        (*Ub).x += Ubn*nx;
        (*Ub).y += Ubn*ny;
        (*Ub).z += Ubn*nz;
	*/
	//printf("%f %f %f = %f %f %f = %f %f %f\n", Ua.x, Ua.y, Ua.z, Uc.x, Uc.y, Uc.z, (*Ub).x, (*Ub).y, (*Ub).z);
}

double u_Werner(double nu, double y, double utau)
{
	double yplus = utau * y / nu; 
	double A=8.3, B=1./7.;
	
	if(yplus<=11.81) return yplus*utau;
	else return A * pow(yplus, B) * utau;
};

double f_Werner(double nu, double u, double y, double utau)
{
	double ym=11.81*nu/utau;
	double A=8.3, B=1./7.;
	
	if( fabs(u) <= nu/(2*ym) * pow(A, 2./(1.-B) ) ) return utau*utau - u/y*nu;
	else return utau*utau -  u/fabs(u) * pow( 0.5*(1-B) * pow(A, (1+B)/(1-B)) * pow(nu/y, 1+B) + (1+B)/A * pow(nu/y, B) * fabs(u), 2/(1+B) );
}

double df_Werner(double nu, double u, double y, double utau)
{
	double eps=1.e-7;
	return ( f_Werner(nu, u, y, utau+eps) - f_Werner(nu, u, y, utau-eps) ) / ( 2*eps ) ;
}

double u_Cabot(double nu, double y, double utau, double dpdn)
{
	return utau * utau * integrate_F( nu, utau, y );// + dpdn * integrate_Fy( nu, utau, 0, y );
};

double u_Cabot_roughness(double nu, double y, double utau, double dpdn, double ks)
{
  return utau * utau * integrate_F( nu, utau, y, ks );// + dpdn * integrate_Fy( nu, utau, 0, y );
};

double u_loglaw_roughness(double nu, double y, double utau, double ks)
{
  // return utau / kappa * log ( (y+ks)/ks );
  
  double yp =  utau * y / nu, kp =  utau * ks / nu;
  return utau * ( 1. / kappa * log(y/ks) + 8.5 );
  //return utau * ( 1. / kappa * log(yp) + 5.2 - 
  //		  (5.2 - 8.5 +  1./kappa * log(kp) ) * sin(0.4258 * ( log(kp) - 0.811 ) ) );
  // wallfunction of Cebeci and Bradshaw, 1977 
  // See p.859 of An analytical wall-function for turbulent flows and heat transfer over rough walls (2006)

};

double f_Cabot(double nu, double u, double y, double utau, double dpdn)
{
	return utau * utau * integrate_F( nu, utau, y ) - u;
}

double f_Cabot_roughness(double nu, double u, double y, double utau, double dpdn, double ks)
{
	return utau * utau * integrate_F( nu, utau, y, ks ) - u;
}

double f_loglaw_roughness(double nu, double u, double y, double utau, double ks)
{
  return u_loglaw_roughness(nu, y, utau, ks) - u;//utau / kappa * log ( (y+ks)/ks ) - u;
  //double yp =  utau * y / nu, kp =  utau * ks / nu;
  //return utau * ( 1. / kappa * log(yp) + 5.2 -
  //                (5.2 - 8.5 +  1./kappa * log(kp) ) * sin(0.4258 * ( log(kp) - 0.811 ) ) ) - u;
}

double df_Cabot(double nu, double u, double y, double utau, double dpdn)
{
	double eps=1.e-7;
	return ( f_Cabot(nu, u, y, utau+eps, dpdn) - f_Cabot(nu, u, y, utau-eps, dpdn) ) / ( 2*eps ) ;
}

double df_Cabot_roughness(double nu, double u, double y, double utau, double dpdn, double ks)
{
	double eps=1.e-7;
	return ( f_Cabot_roughness(nu, u, y, utau+eps, dpdn, ks) - f_Cabot_roughness (nu, u, y, utau-eps, dpdn, ks) ) / ( 2*eps ) ;
}

double df_loglaw_roughness(double nu, double u, double y, double utau, double ks)
{
  //return 1. / kappa * log ( (y+ks)/ks );
  double eps=1.e-7;
  return ( f_loglaw_roughness(nu, u, y, utau+eps, ks) - f_loglaw_roughness(nu, u, y, utau-eps, ks) ) / ( 2*eps ) ;
}

double nu_t(double yplus)// in fact, nu_t/nu
{
	return kappa * yplus * pow ( 1. - exp( - yplus / 19. ) , 2.0 );
};
/*
	integrate F dy 	= integrate [ 1/ (nut + nu) ] dy
				= integrate [ 1/ (nu) / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] dy
				= integrate [ 1/ utau / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] d(yp)
				= 1/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)

	integrate Fy dy	= y/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
				= nu/utau^2 * integrate [  yp / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
*/

int pre_integrate_flag=0;
int n_yp=0;
int interval_yp=2;
int max_yp=1e6;
double *integration_buffer;

void pre_integrate ()
{
	if(pre_integrate_flag) return;
	else pre_integrate_flag=1;
	
	n_yp = ( max_yp / interval_yp ) ;
	
	integration_buffer = new double [ n_yp + 1 ];
	
	integration_buffer[0] = 0.;
	
	for(int i=1; i<=n_yp; i++) {
		int N=24;
		double ya_p = (double)(i-1) * interval_yp;
		double yb_p = (double)(i) * interval_yp;
		double ydiff = yb_p - ya_p, dy_p = ydiff / (double)N;
		std::vector<double> E(N+1);
		double val=0, ybegin_p=ya_p;
		
		for(int k=0; k<=N; k++) E[k] =  1. / ( 1. + nu_t(ybegin_p + dy_p*k ) );
		
		for(int k=0; k<N; k++) {
			double F[4];
			F[0] = E[k];
			F[1] = 1./ ( 1. + nu_t (ybegin_p+dy_p*1./3.) );
			F[2] = 1./ ( 1. + nu_t (ybegin_p+dy_p*2./3.) );
			F[3] = E [k+1];
			val += dy_p  / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.;
			ybegin_p += dy_p;
		}
		integration_buffer[i] = integration_buffer[i-1] + val;
	}
};

double integrate_F (double nu, double utau, double yb)
{
	double val=0;
	
	pre_integrate ();
	
	double ya_plus = 0 * utau / nu;
	double yb_plus = yb * utau / nu;
	
	
	int ib = (int) ( yb_plus / (double) interval_yp );
	int N=4;
	
	if ( yb_plus <= (double) max_yp ) {
		double int_a = 0;
		double int_b = ( integration_buffer[ib+1] - integration_buffer[ib] ) / (double) (interval_yp) * ( yb_plus - (double)(ib) * interval_yp ) + integration_buffer[ib];
		
		return ( int_b - int_a ) / utau;
	}
	else {
		val = integration_buffer [ n_yp ];
		ya_plus = max_yp;
	}
	
	
	double ydiff = yb_plus - ya_plus, dy = ydiff / (double)N;
	std::vector<double> E(N+1);
	double ybegin=ya_plus;
	
	for(int i=0; i<=N; i++) {
		E[i] =  1./ ( 1. + nu_t(ybegin + dy*i ) );
	}
	
	for(int i=0; i<N; i++) {
	        double F[4];
	        F[0] = E[i];
		F[1] = 1./ ( 1. + nu_t(ybegin+dy*1./3.) );
		F[2] = 1./ ( 1. + nu_t(ybegin+dy*2./3.) );
		F[3] = E[i+1];
		val += dy / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.;
		ybegin += dy;
	}
	
	val /= utau;
	return val;
};

// ks: roughness size
double integrate_F(double nu, double utau, double yb, double ks)
{
  double ks_plus = /*fabs*/(utau * ks / nu);
  
  //  if(nu<1.e-40) printf("Zero viscosity !!!\n");

	double yp_shift = 0.9*(sqrt(ks_plus)-(ks_plus)*exp(-ks_plus/6.));
	if(yp_shift<0) return 0;
	
	if( fabs(utau)>1.e-9 ) {
	  double y_shift = yp_shift * nu / utau;
	  return integrate_F (nu, utau, yb + y_shift) - integrate_F (nu, utau, y_shift);
	}
	else return 0;//integrate_F (nu, utau, yb);
};

double find_utau_Cabot(double nu, double u, double y, double guess, double dpdn)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x = x0 - f_Cabot(nu, u, y, x0, dpdn)/df_Cabot(nu, u, y, x0, dpdn);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double find_utau_Cabot_roughness(double nu, double u, double y, double guess, double dpdn, double ks)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x = x0 - f_Cabot_roughness(nu, u, y, x0, dpdn, ks)/df_Cabot_roughness(nu, u, y, x0, dpdn, ks);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double find_utau_loglaw_roughness(double nu, double u, double y, double guess, double ks)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x = x0 - f_loglaw_roughness(nu, u, y, x0, ks)/df_loglaw_roughness(nu, u, y, x0, ks);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double find_utau_Werner(double nu, double u, double y, double guess)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<20; i++) {
		x = x0 - f_Werner(nu, u, y, x0)/df_Werner(nu, u, y, x0);
		if( fabs(x0 - x) < 1.e-7 ) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=19 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};


double sign(double a)
{
	if(a>0) return 1;
	else if(a<0) return -1;
	else return 0;
}
