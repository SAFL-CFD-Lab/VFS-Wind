/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

/*Modeling the rotor blades */

#include "variables.h"

#include <algorithm>
using namespace std;
extern PetscInt ti, tiout;
extern	double dfunc_2h(double r);
extern	double dfunc_h(double r);
extern	double dfunc_2hs1(double r);
extern	double dfunc_2hs2(double r);
extern	double dfunc_2hs3(double r);
extern	double dfunc_2hs4(double r);
extern	double dfunc_2hs5(double r);
extern	double dfunc_2hs6(double r);
extern	double dfunc_4hs1(double r);

extern  double dfunc_nhs2(double r, double n); // with plateau, n half width

extern double dfunc_4h2peak(double r);

extern double dfunc_nhs1(double r, double n);

extern	double dfunc_4h(double r);
extern	double dfunc_4htail(double r);
extern  double dfunc_s3h(double r);
extern  double dfunc_s4h(double r);
extern  double dfunc_sc4h(double r);
extern double dfunc_nh(double r, double n);
extern double dfunc_exp(double r, double n);

extern  double dfunc_4uniform(double r);
extern  double dfunc_6uniform(double r);
extern	PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects); 

extern Cmpnts ArbitraryRotate(Cmpnts p,double theta,Cmpnts r);
//double  coef_cr1 = 2.2, coef_cr2 = 1.0, coef_cr3 = 4.0;

int Itpwidth=(int)halfwidth_dfunc+2;
//int Itpwidth=20;



PetscErrorCode surface_read_xpatch(IBMNodes *ibm, int ibi, FSInfo *fsi, char fname[80], double reflength)
{

  	int	rank;
  	int	n_v , n_elmt ;
  	PetscReal	*x_bp , *y_bp , *z_bp ;
  	int	*nv1 , *nv2 , *nv3 ;
  	PetscReal	*nf_x, *nf_y, *nf_z;
  	int	i,ii;
  	int	n1e, n2e, n3e;
  	PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  	PetscReal     dr;
  	PetscReal     *dA ;//area
  	PetscReal	*nt_x, *nt_y, *nt_z;
  	PetscReal	*ns_x, *ns_y, *ns_z;

  	char   ss[20];
  	char string[128];

      	PetscPrintf(PETSC_COMM_WORLD, "xc=%le, yc=%le, zc=%le for %i th turbine\n", fsi->x_c, fsi->y_c, fsi->z_c, ibi);
      	PetscPrintf(PETSC_COMM_WORLD, "nx=%le, ny=%le, nz=%le for %i th turbine\n", fsi->nx_tb, fsi->ny_tb, fsi->nz_tb, ibi);

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  	if(!rank) { // root processor read in the data
    		FILE *fd;
    		PetscPrintf(PETSC_COMM_SELF, "READ %s\n", fname);
    		char filen[160];  
    		sprintf(filen,"%s/%s" ,path, fname);
 
    		fd = fopen(filen, "r"); 
    		if (!fd) printf("Cannot open %s !!", filen),exit(0);
    		else printf("Opened %s !\n", filen);
    		n_v =0;

    		if (fd) {

      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      
      			fscanf(fd, "%i",&n_v);
      			PetscPrintf(PETSC_COMM_SELF, "number of nodes %d\n",n_v);
 
      			ibm->n_v = n_v;
      
      			MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      			x_bp = ibm->x_bp;	// seokkoo
      			y_bp = ibm->y_bp;
      			z_bp = ibm->z_bp;

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
 
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

      			for (i=0; i<n_v; i++) {
				fscanf(fd, "%le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
		
				x_bp[i]=x_bp[i]/reflength;
				y_bp[i]=y_bp[i]/reflength;
				z_bp[i]=z_bp[i]/reflength;


			        ibm->x_bp_i[i] = x_bp[i];
			        ibm->y_bp_i[i] = y_bp[i];
			        ibm->z_bp_i[i] = z_bp[i];

			        x_bp[i] += fsi->x_c;
			        y_bp[i] += fsi->y_c;
			        z_bp[i] += fsi->z_c;

			        ibm->x_bp0[i] = x_bp[i];
			        ibm->y_bp0[i] = y_bp[i];
			        ibm->z_bp0[i] = z_bp[i];

	
				ibm->x_bp[i] = x_bp[i];
				ibm->y_bp[i] = y_bp[i];
				ibm->z_bp[i] = z_bp[i];

				ibm->x_bp_o[i] = x_bp[i];
				ibm->y_bp_o[i] = y_bp[i];
				ibm->z_bp_o[i] = z_bp[i];

				ibm->tmprt[i] = 0.;

				ibm->u[i].x = 0.;
				ibm->u[i].y = 0.;
				ibm->u[i].z = 0.;

				ibm->uold[i].x = 0.;
				ibm->uold[i].y = 0.;
				ibm->uold[i].z = 0.;

				ibm->urm1[i].x = 0.;
				ibm->urm1[i].y = 0.;
				ibm->urm1[i].z = 0.;
			}
      			i=0;
		      	PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

			// no
      			MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
		      	MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

			fscanf(fd, "%i", &ii);
      			PetscPrintf(PETSC_COMM_SELF, "number of -- %d\n", ii);

		
  			char string1[128];
			fscanf(fd, "%s", &string1);
//      			fgets(string1, 64, fd);
//      			PetscPrintf(PETSC_COMM_SELF, "the string %s\n",string1);
			fscanf(fd, "%i %i", &n_elmt, &ii);
			ibm->n_elmt=n_elmt;

      			PetscPrintf(PETSC_COMM_SELF, "number of elemts %d\n",n_elmt);
      			MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      			PetscMalloc(n_elmt*sizeof(int), &nv1);
      			PetscMalloc(n_elmt*sizeof(int), &nv2);
      			PetscMalloc(n_elmt*sizeof(int), &nv3);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

			// rotor model
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));


        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));


		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
       		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
                	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));


        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

			PetscMalloc(n_elmt*sizeof(int), &(ibm->color));
			PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

      			for (i=0; i<n_elmt; i++) {
				fscanf(fd, "%i %i %i %i %i %i\n", &nv1[i], &nv2[i], &nv3[i], &ii, &(ibm->color[i]), &ii);
				nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      			}
      			ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

			int color_translate = ibm->color[0]-1;
      			for (i=0; i<n_elmt; i++) {
				ibm->color[i]-=color_translate;
      			}

		        i=0;
		        PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

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
      
		      	if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
			  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
				ns_x[i] = 1.;     
				ns_y[i] = 0.;     
				ns_z[i] = 0. ;
	
				nt_x[i] = 0.;
				nt_y[i] = 1.;
				nt_z[i] = 0.;
		        } else {
				ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
				ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
				ns_z[i] = 0. ;
	
				nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
				nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
				nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
		        }
      
		        dA[i] = dr/2.; 
      
      			// Calc the center of the element
      			ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
		        ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
		        ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
		}
    
    
    		ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    		ibm->dA = dA;
    		ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
	        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    
    
	        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
    		MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    		int ti=0;
    		FILE *f;
    		sprintf(filen, "%s/%s_%2.2d_nf.dat",path,fname,ibi);
    		f = fopen(filen, "w");
    		PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z, color\n");
    		PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-13]=CELLCENTERED)\n", n_v, n_elmt);
    		for (i=0; i<n_v; i++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    		}
    		for (i=0; i<n_v; i++) {
   	   		PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
	    	}
    		for (i=0; i<n_v; i++) {	
	  	      	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    		}
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	   	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->color[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	        }
    
    		fclose(f);

  	}
  	else if (rank) {
    		MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	        ibm->n_v = n_v;
    
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
	  
		x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
		y_bp = ibm->y_bp;
		z_bp = ibm->z_bp;
	  
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
    
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
 
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

	        for (i=0; i<n_v; i++) {
			ibm->tmprt[i] = 0.;

	  		ibm->u[i].x = 0.;
  		        ibm->u[i].y = 0.;
	    	        ibm->u[i].z = 0.;

		        ibm->uold[i].x = 0.;
		        ibm->uold[i].y = 0.;
	   	        ibm->uold[i].z = 0.;
      
	  	        ibm->urm1[i].x = 0.;
		        ibm->urm1[i].y = 0.;
	 	        ibm->urm1[i].z = 0.;      
	        }

	        MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        
	        MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
	        MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

 	        MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	       	MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	       	MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    

	       	MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	       	n_elmt = ibm->n_elmt;

	       	PetscMalloc(n_elmt*sizeof(int), &nv1);
	       	PetscMalloc(n_elmt*sizeof(int), &nv2);
	       	PetscMalloc(n_elmt*sizeof(int), &nv3);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

   	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

	       	ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
	       	ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
	       	ibm->dA = dA;
	       	ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
	       	ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

	       	// rotor model
	       	//


        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));



        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));


	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
          	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
              	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

		PetscMalloc(n_elmt*sizeof(int), &(ibm->color));

		PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

	    	MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	}
	PetscPrintf(PETSC_COMM_WORLD, "Read %s !\n", fname);

    	for (i=0; i<ibm->n_elmt; i++) {
		ibm->Urel_mean[i]=0.0;
		ibm->Fr_mean[i]=0.0;
		ibm->Ft_mean[i]=0.0;
		ibm->Fa_mean[i]=0.0;
		ibm->Ur_mean[i]=0.0;
		ibm->Ut_mean[i]=0.0;
		ibm->Ua_mean[i]=0.0;
		ibm->AOA_mean[i]=0.0;

	}

	count_AL=0.0;
  	return(0);
}






PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ, char fname[80], double reflength)
{
  	int	rank;
  	int	n_v , n_elmt ;
  	PetscReal	*x_bp , *y_bp , *z_bp ;
  	int	*nv1 , *nv2 , *nv3 ;
  	PetscReal	*nf_x, *nf_y, *nf_z;
  	int	i,ii;
  	int	n1e, n2e, n3e;
  	PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  	PetscReal     dr;
  	PetscReal     *dA ;//area
  	PetscReal	*nt_x, *nt_y, *nt_z;
  	PetscReal	*ns_x, *ns_y, *ns_z;

  	char   ss[20];
  	char string[128];

      	PetscPrintf(PETSC_COMM_WORLD, "xc=%le, yc=%le, zc=%le for %i th turbine\n", fsi->x_c, fsi->y_c, fsi->z_c, ibi);
      	PetscPrintf(PETSC_COMM_WORLD, "nx=%le, ny=%le, nz=%le for %i th turbine\n", fsi->nx_tb, fsi->ny_tb, fsi->nz_tb, ibi);

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  	if(!rank) { // root processor read in the data
    		FILE *fd;
    		PetscPrintf(PETSC_COMM_SELF, "READ %s\n", fname);
    		char filen[160];  
    		sprintf(filen,"%s/%s" ,path, fname);
 
    		fd = fopen(filen, "r"); 
    		if (!fd) printf("Cannot open %s !!", filen),exit(0);
    		else printf("Opened %s !\n", filen);
    		n_v =0;

    		if (fd) {
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      
      			fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);
      			PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
      
      			ibm->n_v = n_v;
      			ibm->n_elmt = n_elmt;      
      
      			MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      			x_bp = ibm->x_bp;	// seokkoo
      			y_bp = ibm->y_bp;
      			z_bp = ibm->z_bp;

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
 
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));
			
      			for (i=0; i<n_v; i++) {
				fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
				x_bp[i]=x_bp[i]/reflength;
				y_bp[i]=y_bp[i]/reflength;
				z_bp[i]=z_bp[i]/reflength;

				// rotate the axis. Assume the one from gridgen is (0,0,1)
				/*
				double n_rot[3];
				n_rot[0]=-fsi->ny_tb; n_rot[1]=fsi->nx_tb; n_rot[2]=0.0;
				double ux=n_rot[0], uy=n_rot[1], uz=n_rot[2];

				double ang_rot;
				double rrr=fsi->nz_tb;
				ang_rot=acos(rrr);
				double cos_=cos(ang_rot), sin_=sin(ang_rot);
				*/

				// Problem for some cases    .........................
				/*
				Cmpnts p, r, q;
				p.x=x_bp[i]; p.y=y_bp[i]; p.z=z_bp[i];
				r.x=-fsi->ny_tb; r.y=fsi->nx_tb; r.z=0.0;
				double theta=acos(fsi->nz_tb);
				q=ArbitraryRotate(p,theta,r);
				x_bp[i]=q.x;	
				y_bp[i]=q.y;	
				z_bp[i]=q.z;	
				*/
				// 				.................................
				
				/*
				double mat_rot[3][3];

				mat_rot[0][0]=cos_+pow(ux,2)*(1.0-cos_);
				mat_rot[0][1]=ux*uy*(1.0-cos_)-uz*sin_;
				mat_rot[0][2]=ux*uz*(1.0-cos_)+uy*sin_;

				mat_rot[1][0]=uy*ux*(1.0-cos_)+uz*sin_;
				mat_rot[1][1]=cos_+pow(uy,2)*(1.0-cos_);
				mat_rot[1][2]=uy*uz*(1.0-cos_)-ux*sin_;

				mat_rot[2][0]=uz*ux*(1.0-cos_)-uy*sin_;
				mat_rot[2][1]=uz*uy*(1.0-cos_)+ux*sin_;
				mat_rot[2][2]=cos_+pow(uz,2)*(1.0-cos_);
				*/

				// double xb=x_bp[i], yb=y_bp[i], zb=z_bp[i];

				/*
				x_bp[i]=mat_rot[0][0]*xb+mat_rot[0][1]*yb+mat_rot[0][2]*zb;
				y_bp[i]=mat_rot[1][0]*xb+mat_rot[1][1]*yb+mat_rot[1][2]*zb;
				z_bp[i]=mat_rot[2][0]*xb+mat_rot[2][1]*yb+mat_rot[2][2]*zb;
				*/

				/*
				double wCv[3];
				wCv[0]=uy*zb-uz*yb; wCv[1]=uz*xb-ux*zb; wCv[2]=ux*yb-uy*xb;
				double wPv=ux*xb+uy*yb+uz*zb;
				x_bp[i]=xb*cos_+wCv[0]*sin_+wPv*ux*(1.0-cos_);
				y_bp[i]=yb*cos_+wCv[1]*sin_+wPv*uy*(1.0-cos_);
				z_bp[i]=zb*cos_+wCv[2]*sin_+wPv*uz*(1.0-cos_);
				*/

			        ibm->x_bp_i[i] = x_bp[i];
			        ibm->y_bp_i[i] = y_bp[i];
			        ibm->z_bp_i[i] = z_bp[i];

			        x_bp[i] += fsi->x_c;
			        y_bp[i] += fsi->y_c;
			        z_bp[i] += fsi->z_c;

			        ibm->x_bp0[i] = x_bp[i];
			        ibm->y_bp0[i] = y_bp[i];
			        ibm->z_bp0[i] = z_bp[i];


				if (OneDmZ) {
					double R=fsi->r_rotor/reflength_wt;
					double rr = loc_refvel*2.0*R;
					x_bp[i]-=rr*fsi->nx_tb; 
					y_bp[i]-=rr*fsi->ny_tb; 
					z_bp[i]-=rr*fsi->nz_tb; 
				}
	
				ibm->x_bp[i] = x_bp[i];
				ibm->y_bp[i] = y_bp[i];
				ibm->z_bp[i] = z_bp[i];

				ibm->x_bp_o[i] = x_bp[i];
				ibm->y_bp_o[i] = y_bp[i];
				ibm->z_bp_o[i] = z_bp[i];

				ibm->tmprt[i] = 0.;

				ibm->u[i].x = 0.;
				ibm->u[i].y = 0.;
				ibm->u[i].z = 0.;

				ibm->uold[i].x = 0.;
				ibm->uold[i].y = 0.;
				ibm->uold[i].z = 0.;

				ibm->urm1[i].x = 0.;
				ibm->urm1[i].y = 0.;
				ibm->urm1[i].z = 0.;
			}
      			i=0;
		      	PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

      			MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
		      	MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      			PetscMalloc(n_elmt*sizeof(int), &nv1);
      			PetscMalloc(n_elmt*sizeof(int), &nv2);
      			PetscMalloc(n_elmt*sizeof(int), &nv3);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

			// rotor model
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));


        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


 		       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));





		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
       		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
                	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));


        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));


      			for (i=0; i<n_elmt; i++) {
				fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, &ss, nv1+i, nv2+i, nv3+i);
				nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
		      	}
		        ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

		        i=0;
		        PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

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
      
		      	if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
			  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
				ns_x[i] = 1.;     
				ns_y[i] = 0.;     
				ns_z[i] = 0. ;
	
				nt_x[i] = 0.;
				nt_y[i] = 1.;
				nt_z[i] = 0.;
		        } else {
				ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
				ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
				ns_z[i] = 0. ;
	
				nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
				nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
				nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
		        }
      
		        dA[i] = dr/2.; 
      
      			// Calc the center of the element
      			ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
		        ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
		        ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
		}
    
    
    		ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    		ibm->dA = dA;
    		ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
	        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    
    
	        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
    		MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		int ti=0;
    		FILE *f;
    		sprintf(filen, "%s/%s_%2.2d_nf.dat",path,fname,ibi);
    		f = fopen(filen, "w");
    		PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    		PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    		for (i=0; i<n_v; i++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    		}
    		for (i=0; i<n_v; i++) {
   	   		PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
	    	}
    		for (i=0; i<n_v; i++) {	
	  	      	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    		}
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
		        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	   	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]);
	        }
	        for (i=0; i<n_elmt; i++) {
	 	        PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	        }
    
    		fclose(f);

  	}
  	else if (rank) {
    		MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	        ibm->n_v = n_v;
    
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
	  
		x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
		y_bp = ibm->y_bp;
		z_bp = ibm->z_bp;
	  
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
    
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
 
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
	        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
	        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

	        for (i=0; i<n_v; i++) {

			ibm->tmprt[i] = 0.;

	  		ibm->u[i].x = 0.;
  		        ibm->u[i].y = 0.;
	    	        ibm->u[i].z = 0.;

		        ibm->uold[i].x = 0.;
		        ibm->uold[i].y = 0.;
	   	        ibm->uold[i].z = 0.;
      
	  	        ibm->urm1[i].x = 0.;
		        ibm->urm1[i].y = 0.;
	 	        ibm->urm1[i].z = 0.;      
	        }

	        MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        
	        MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
	        MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

 	        MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	       	MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
	       	MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    

	       	MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	       	n_elmt = ibm->n_elmt;

	       	PetscMalloc(n_elmt*sizeof(int), &nv1);
	       	PetscMalloc(n_elmt*sizeof(int), &nv2);
	       	PetscMalloc(n_elmt*sizeof(int), &nv3);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

   	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

	       	ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
	       	ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
	       	ibm->dA = dA;
	       	ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
	       	ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

	       	// rotor model
	       	//


        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));



        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

	       	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));




	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
          	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
              	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));


	    	MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

	        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	}
	PetscPrintf(PETSC_COMM_WORLD, "Read %s !\n", fname);

    	for (i=0; i<ibm->n_elmt; i++) {
		ibm->Urel_mean[i]=0.0;
		ibm->Fr_mean[i]=0.0;
		ibm->Ft_mean[i]=0.0;
		ibm->Fa_mean[i]=0.0;
		ibm->Ur_mean[i]=0.0;
		ibm->Ut_mean[i]=0.0;
		ibm->Ua_mean[i]=0.0;
		ibm->AOA_mean[i]=0.0;

	}

	count_AL=0.0;

  	return(0);
}

PetscErrorCode ACL_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, double reflength)
{
  	int	rank;
  	int	n_v , n_elmt ;
  	PetscReal	*x_bp , *y_bp , *z_bp ;
  	int	*nv1 , *nv2 , *nv3 ;
  	PetscReal	*nf_x, *nf_y, *nf_z;
  	int	i,ii;
  	int	n1e, n2e, n3e;
  	PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  	PetscReal     dr;
  	//Added 4/1/06 iman
  	PetscReal     *dA ;//area
  	PetscReal	*nt_x, *nt_y, *nt_z;
  	PetscReal	*ns_x, *ns_y, *ns_z;

	int nv_blade, nelmt_blade;

  	char   ss[20];
  	//double xt;
  	char string[128];

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  	if(!rank) { // root processor read in the data
    		FILE *fd;
    		PetscPrintf(PETSC_COMM_SELF, "READ acldata\n");
    		char filen[80];  
    		sprintf(filen,"%s/acldata%3.3d" , path, 0);

    		fd = fopen(filen, "r"); 
    		if (!fd) printf("Cannot open %s !!", filen),exit(0);
    		else printf("Opened %s !\n", filen);

    		n_v =0;

    		if (fd) {

      			fscanf(fd, "%i ",&n_v);
      			n_elmt = n_v - 1;
      			PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements & blades %d %d %d \n",n_v, n_elmt, num_blade);
        		nv_blade=n_v;
			nelmt_blade=n_elmt; 

      			ibm->n_v = n_v * num_blade;
      			ibm->n_elmt = n_elmt * num_blade;      
 
      			n_v=ibm->n_v;
      			n_elmt=ibm->n_elmt;      
 

      			MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	    
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      			x_bp = ibm->x_bp; // seokkoo
      			y_bp = ibm->y_bp;
      			z_bp = ibm->z_bp;

        		//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhx));
        		//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhy));
        		//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhz));

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      			PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

			int nb,j;
      			for (nb=0; nb<num_blade; nb++) {
        			if (nb != 0) fscanf(fd, "%i ", &ii);
        			for (j=0; j<nv_blade; j++) {
          				i = nb*nv_blade + j;
          
	  				fscanf(fd, "%le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
				
					x_bp[i]=x_bp[i]/reflength;
					y_bp[i]=y_bp[i]/reflength;
					z_bp[i]=z_bp[i]/reflength;


					// rotate the axis. Assume the one from gridgen is (0,0,1)
	

					// not work for some cases .................
					/*

					Cmpnts p, r, q;
					p.x=x_bp[i]; p.y=y_bp[i]; p.z=z_bp[i];
					r.x=-fsi->ny_tb; r.y=fsi->nx_tb; r.z=0.0;
					double theta=acos(fsi->nz_tb);
					q=ArbitraryRotate(p,theta,r);
					x_bp[i]=q.x;	
					y_bp[i]=q.y;	
					z_bp[i]=q.z;	
					*/
					// ....................
	
					//double n_rot[3];
					//n_rot[0]=-fsi->ny_tb; n_rot[1]=fsi->nx_tb; n_rot[2]=0.0;
					//double ux=n_rot[0], uy=n_rot[1], uz=n_rot[2];

					//double ang_rot;
					//double rrr=fsi->nz_tb;
					//ang_rot=acos(rrr);
					//double cos_=cos(ang_rot), sin_=sin(ang_rot);

					/*
					double mat_rot[3][3];

					mat_rot[0][0]=cos_+pow(ux,2)*(1.0-cos_);
					mat_rot[0][1]=ux*uy*(1.0-cos_)-uz*sin_;
					mat_rot[0][2]=ux*uz*(1.0-cos_)+uy*sin_;

					mat_rot[1][0]=uy*ux*(1.0-cos_)+uz*sin_;
					mat_rot[1][1]=cos_+pow(uy,2)*(1.0-cos_);
					mat_rot[1][2]=uy*uz*(1.0-cos_)-ux*sin_;

					mat_rot[2][0]=uz*ux*(1.0-cos_)-uy*sin_;
					mat_rot[2][1]=uz*uy*(1.0-cos_)+ux*sin_;
					mat_rot[2][2]=cos_+pow(uz,2)*(1.0-cos_);
					*/

					//double xb=x_bp[i], yb=y_bp[i], zb=z_bp[i];

					/*
					x_bp[i]=mat_rot[0][0]*xb+mat_rot[0][1]*yb+mat_rot[0][2]*zb;
					y_bp[i]=mat_rot[1][0]*xb+mat_rot[1][1]*yb+mat_rot[1][2]*zb;
					z_bp[i]=mat_rot[2][0]*xb+mat_rot[2][1]*yb+mat_rot[2][2]*zb;
					*/

					//double wCv[3];
					//wCv[0]=uy*zb-uz*yb; wCv[1]=uz*xb-ux*zb; wCv[2]=ux*yb-uy*xb;
					//double wPv=ux*xb+uy*yb+uz*zb;
					//x_bp[i]=xb*cos_+wCv[0]*sin_+wPv*ux*(1.0-cos_);
					//y_bp[i]=yb*cos_+wCv[1]*sin_+wPv*uy*(1.0-cos_);
					//z_bp[i]=zb*cos_+wCv[2]*sin_+wPv*uz*(1.0-cos_);

	        			ibm->x_bp_i[i] = x_bp[i];
        				ibm->y_bp_i[i] = y_bp[i];
        				ibm->z_bp_i[i] = z_bp[i];

	        			x_bp[i] += fsi->x_c;
        				y_bp[i] += fsi->y_c;
        				z_bp[i] += fsi->z_c;

	        			ibm->x_bp0[i] = x_bp[i];
        				ibm->y_bp0[i] = y_bp[i];
        				ibm->z_bp0[i] = z_bp[i];

	
					ibm->x_bp[i] = x_bp[i];
					ibm->y_bp[i] = y_bp[i];
					ibm->z_bp[i] = z_bp[i];

					ibm->x_bp_o[i] = x_bp[i];
					ibm->y_bp_o[i] = y_bp[i];
					ibm->z_bp_o[i] = z_bp[i];


					ibm->tmprt[i] = 0.;

					ibm->u[i].x = 0.;
					ibm->u[i].y = 0.;
					ibm->u[i].z = 0.;

					ibm->uold[i].x = 0.;
					ibm->uold[i].y = 0.;
					ibm->uold[i].z = 0.;

					ibm->urm1[i].x = 0.;
					ibm->urm1[i].y = 0.;
					ibm->urm1[i].z = 0.;
	      			}
			}
      			i=0;
      			PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

      			MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

     	 		MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      			MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      			MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      			PetscMalloc(n_elmt*sizeof(int), &nv1);
      			PetscMalloc(n_elmt*sizeof(int), &nv2);
      			PetscMalloc(n_elmt*sizeof(int), &nv3);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      			// Added 4/1/06 iman
      			PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      			// Added 6/4/06 iman
      			//PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      			// end added


      			// added 12-7-2010 xyang
      			//

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));



        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


	        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));




	        	PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
        		PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
        		PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

        		PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
        		PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
        		PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

       			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
              		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
              		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

      			if (surface_p_out) {
	        		PetscMalloc(n_elmt*sizeof(int), &(ibm->ib_elmt));
	        		PetscMalloc(n_elmt*sizeof(int), &(ibm->jb_elmt));
        			PetscMalloc(n_elmt*sizeof(int), &(ibm->kb_elmt));

        			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->xb_elmt));
        			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->yb_elmt));
	        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->zb_elmt));

        			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->p_elmt));
        			PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->tau_elmt));
	      		}

        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                       	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);


			PetscMalloc(n_elmt*sizeof(int), &ibm->color);

			// colored as blade index 
      			for (nb=0; nb<num_blade; nb++) {
        		for (j=0; j<nelmt_blade; j++) {
	  			i = nb*nelmt_blade + j;
	  			ibm->color[i] = nb+1;
        		}
      			}


      			for (nb=0; nb<num_blade; nb++) {
        		for (j=0; j<nelmt_blade; j++) {
	  			i = nb*nelmt_blade + j;
				ii = nb*nv_blade + j;
	  			nv1[i] = ii; nv2[i] = ii + 1; nv3[i] = ii + 1;
        		}
      			}

	      		ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      			i=0;
	      		PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);

      			fclose(fd);
	    	}
     
		int nb, j;
		// for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line. and dA denotes the length of each element.
      		for (nb=0; nb<num_blade; nb++) {
      
    		for (j=0; j<nelmt_blade; j++) {
      			i = nb*nelmt_blade + j;

      			n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; 
      			dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
      			dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
      			dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
      
      
      			nf_x[i] = dx12;
      			nf_y[i] = dy12;
      			nf_z[i] = dz12;
      
      			dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      			nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
			ns_x[i] = 0.;     
			ns_y[i] = 0.;     
			ns_z[i] = 0. ;
	
			nt_x[i] = 0.;
			nt_y[i] = 0.;
			nt_z[i] = 0.;
     
			dA[i] = dr;

      			// Calc the center of the element
      			ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e])/2.;
      			ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e])/2.;
      			ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e])/2.;
	
      		}
      
    		}
     
    		ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
	    	//Added 4/1/06 iman
	    	ibm->dA = dA;
	    	ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
	    	ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    
    
    
	    	MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
	    	MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
	    	// Added 4/1/06 iman
	    	MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	    	// Added 4/2/06 iman
	    	MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	    	MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
	    	// Added 6/4/06 iman
	    	MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	    	MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);


	    	MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	 	/*    PetscFree(dA); */
		/*     PetscFree(nf_x);PetscFree(nf_y);PetscFree(nf_z); */
		/*     PetscFree(nt_x);PetscFree(nt_y);PetscFree(nt_z); */
		/*     PetscFree(ns_x);PetscFree(ns_y);PetscFree(ns_z); */
		/*     PetscFree(nv1);PetscFree(nv2);PetscFree(nv3); */
		/*     PetscFree(x_bp);PetscFree(y_bp);PetscFree(z_bp); */


    		int ti=0;
    		FILE *f;
    		sprintf(filen, "%s/line%3.3d_%2.2d_nf.dat",path,ti,ibi);
    		f = fopen(filen, "w");
    		for (i=0; i<ibm->n_elmt; i++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%e %le %le \n", ibm->cent_x[i], ibm->cent_y[i], ibm->cent_z[i]);
    		}
    		fclose(f);



  	}
  	else if (rank) {
    		MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    		ibm->n_v = n_v;
	
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
	  
		x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
		y_bp = ibm->y_bp;
		z_bp = ibm->z_bp;

        	//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhx));
        	//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhy));
        	//PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dhz));
	  
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
 
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    		PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    		PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    		PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));


    		for (i=0; i<n_v; i++) {
			ibm->tmprt[i] = 0.;
      			ibm->u[i].x = 0.;
      			ibm->u[i].y = 0.;
      			ibm->u[i].z = 0.;

      			ibm->uold[i].x = 0.;
      			ibm->uold[i].y = 0.;
      			ibm->uold[i].z = 0.;
      
      			ibm->urm1[i].x = 0.;
      			ibm->urm1[i].y = 0.;
      			ibm->urm1[i].z = 0.;      
    		}
 
      		MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      		MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      		MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
       
    		MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    		MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    		ibm->n_elmt = n_elmt;

    		PetscMalloc(n_elmt*sizeof(int), &nv1);
    		PetscMalloc(n_elmt*sizeof(int), &nv2);
    		PetscMalloc(n_elmt*sizeof(int), &nv3);

    		PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    		//Added 4/1/06 iman
    		PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    		//Added 4/2/06 iman
    		PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    		PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    		PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

    		ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    		ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
    		// Added 4/2/06 iman
    		ibm->dA = dA;
    		ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    		ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    		// Added 6/4/06
    		//PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
    		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
    		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
    		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));


		// added 12-7-2010 xyang
		//

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urel_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));



        		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_rel));

      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


        	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));




      		PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
      		PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
      		PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

      		PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
      		PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
      		PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

       		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
              	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
              	PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));
                         

    		if (surface_p_out) {
 		 	PetscMalloc(n_elmt*sizeof(int), &(ibm->ib_elmt));
      			PetscMalloc(n_elmt*sizeof(int), &(ibm->jb_elmt));
      			PetscMalloc(n_elmt*sizeof(int), &(ibm->kb_elmt));

      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->xb_elmt));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->yb_elmt));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->zb_elmt));

      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->p_elmt));
      			PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->tau_elmt));
    		}


        	PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
       		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
       		PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
                PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);

		PetscMalloc(n_elmt*sizeof(int), &ibm->color);

    		MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		//Added 4/2/06 iman
    		MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
  	}



    	for (i=0; i<ibm->n_elmt; i++) {
		ibm->Urel_mean[i]=0.0;
		ibm->Fr_mean[i]=0.0;
		ibm->Ft_mean[i]=0.0;
		ibm->Fa_mean[i]=0.0;
		ibm->Ur_mean[i]=0.0;
		ibm->Ut_mean[i]=0.0;
		ibm->Ua_mean[i]=0.0;
		ibm->AOA_mean[i]=0.0;


	}

	count_AL=0.0;
	PetscPrintf(PETSC_COMM_WORLD, "Finish Reading acldata file !\n");
  	return(0);
}


// Call in main.c, after reading the grid for turbines
PetscErrorCode Pre_process_save3(UserCtx *user, IBMNodes *ibm, int NumberOfObjects)

{

	PetscReal ts1, te1;
  PetscGetTime(&ts1);  

  	PetscBarrier(PETSC_NULL);

  	DA              da = user->da, fda = user->fda;
  	DALocalInfo  	info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj;

  	Vec		Coor;


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


	int *iclose, *jclose, *kclose;

	int *sum_iclose, *sum_jclose, *sum_kclose;


  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);

  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);


  	int n1e, n2e, n3e;

  	for (ibi=0; ibi<NumberOfObjects; ibi++) {


		iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

		sum_iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			iclose[l]=0;
			jclose[l]=0;
			kclose[l]=0;

			sum_iclose[l]=0;
			sum_jclose[l]=0;
			sum_kclose[l]=0;

		}

	  	PetscBarrier(PETSC_NULL);

	PetscReal ts2, te2;
        PetscGetTime(&ts2);  

		int count[ibm[ibi].n_elmt], sum_count[ibm[ibi].n_elmt];
	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			int imark=0;
			double dmin=1.e6;
			int ic,jc,kc;
	                for (i=lxs; i<lxe; i++)
        	        for (j=lys; j<lye; j++) 
                	for (k=lzs; k<lze; k++){ 
			
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					ic=i; jc=j; kc=k;
				}
                	}
	
			count[l]=1;	

			i=ic; j=jc; k=kc;
		        double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dhx=1.0/aj[k][j][i]/area;

		     	area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dhy=1.0/aj[k][j][i]/area;

		       	area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double dhz=1.0/aj[k][j][i]/area;	

			double dh_grid = 2.0*sqrt(dhx*dhx+dhy*dhy+dhz*dhz);

			//double dmin_global;
			//MPI_Allreduce (&dmin, &dmin_global, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);

			double diff=dmin-dh_grid;
			if (diff>1.e-6) {
				count[l]=0;	
				ic=0; jc=0; kc=0;
				iclose[l]=0; jclose[l]=0; kclose[l]=0;
			} else {
				iclose[l]=ic; jclose[l]=jc; kclose[l]=kc;
				i=ic; j=jc; k=kc;
			}

	        }
	  	PetscBarrier(PETSC_NULL);
        PetscGetTime(&te2);  
      	PetscPrintf(PETSC_COMM_WORLD, "Time for finding the closest  %le\n", te2-ts2);




	  	PetscBarrier(PETSC_NULL);

		MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&count[0], &sum_count[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);


	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			iclose[l]=sum_iclose[l]/sum_count[l];
			jclose[l]=sum_jclose[l]/sum_count[l];
			kclose[l]=sum_kclose[l]/sum_count[l];
		}


	        int iii1 = 0;
		int iii2 = 0;
	        int iii = 0;

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

                        int n1e = ibm[ibi].nv1[i], n2e = ibm[ibi].nv2[i], n3e = ibm[ibi].nv3[i];

			int ic=iclose[l];
			int jc=jclose[l];
			int kc=kclose[l];

			int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
			if (ic1>lxe||ic2<lxs) {
				ibm[ibi].i_min[l]=lxs;
				ibm[ibi].i_max[l]=lxs;
			} else {
				ibm[ibi].i_min[l]=PetscMax(ic1, lxs);
				ibm[ibi].i_max[l]=PetscMin(ic2, lxe);
			}


			int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 
	       		if (jc1>lye||jc2<lys) {
				ibm[ibi].j_min[l]=lys;
				ibm[ibi].j_max[l]=lys;
			} else {
				ibm[ibi].j_min[l]=PetscMax(jc1, lys);
				ibm[ibi].j_max[l]=PetscMin(jc2, lye);
			}


			int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 
	       		if (kc1>lze||kc2<lzs) {
				ibm[ibi].k_min[l]=lzs;
				ibm[ibi].k_max[l]=lzs;
			} else {
				ibm[ibi].k_min[l]=PetscMax(kc1, lzs);
				ibm[ibi].k_max[l]=PetscMin(kc2, lze);
			}


		}


	        free(iclose);
		free(jclose);
		free(kclose);

	        free(sum_iclose);
		free(sum_jclose);
		free(sum_kclose);

	}

  	DAVecRestoreArray(fda, Coor, &coor);

	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

  	PetscBarrier(PETSC_NULL);

        PetscGetTime(&te1);  
      	PetscPrintf(PETSC_COMM_WORLD, "Time for rotor model search  %le\n", te1-ts1);

  	return(0);

}







// Call in main.c, after reading the grid for turbines
PetscErrorCode Pre_process_save2(UserCtx *user, IBMNodes *ibm, int NumberOfObjects)
{

	PetscReal ts1, te1;
	PetscGetTime(&ts1);  
	PetscBarrier(PETSC_NULL);
	DA              da = user->da, fda = user->fda;
	DALocalInfo  	info;
	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt        mx, my, mz; // Dimensions in three directions
	PetscInt        i, j, k, l, ibi;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts		***coor, ***csi, ***eta, ***zet;
	PetscReal 	***aj;
	Vec		Coor;
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
	double *dhx_, *dhy_, *dhz_;
	int *iclose, *jclose, *kclose;
	double *sum_dhx_, *sum_dhy_, *sum_dhz_;
	int *sum_iclose, *sum_jclose, *sum_kclose;
	DAGetGhostedCoordinates(da, &Coor);
	DAVecGetArray(fda, Coor, &coor);
	DAVecGetArray(fda, user->lCsi,  &csi);
	DAVecGetArray(fda, user->lEta,  &eta);
	DAVecGetArray(fda, user->lZet,  &zet);
	DAVecGetArray(da,  user->lAj,  &aj);
	int n1e, n2e, n3e;
	for (ibi=0; ibi<NumberOfObjects; ibi++) {

		dhx_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		dhy_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		dhz_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));

		iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

		sum_dhx_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		sum_dhy_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		sum_dhz_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));

		sum_iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

		for (l=0; l<ibm[ibi].n_elmt; l++) {
			dhx_[l] = 0.0;
			dhy_[l] = 0.0;
			dhz_[l] = 0.0;
			sum_dhx_[l] = 0.0;
			sum_dhy_[l] = 0.0;
			sum_dhz_[l] = 0.0;
			iclose[l]=0;
			jclose[l]=0;
			kclose[l]=0;
			sum_iclose[l]=0;
			sum_jclose[l]=0;
			sum_kclose[l]=0;
		}


		int count[ibm[ibi].n_elmt], sum_count[ibm[ibi].n_elmt];
	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			int imark=0;
			double dmin=1.e6;
			int ic,jc,kc;
	                for (i=lxs; i<lxe; i++)
        	        for (j=lys; j<lye; j++) 
                	for (k=lzs; k<lze; k++){ 
			
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					ic=i; jc=j; kc=k;
				}
                	}
	
			count[l]=1;	
			double dmin_global;
			MPI_Allreduce (&dmin, &dmin_global, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
			double diff=fabs(dmin-dmin_global);
			if (diff>1.e-6) {
				count[l]=0;	
				ic=0; jc=0; kc=0;
				iclose[l]=0; jclose[l]=0; kclose[l]=0;
				dhx_[l]=0.0;
				dhy_[l]=0.0;
				dhz_[l]=0.0;
			} else {

				iclose[l]=ic; jclose[l]=jc; kclose[l]=kc;
				i=ic; j=jc; k=kc;
		                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
				dhx_[l]=1.0/aj[k][j][i]/area;

			       	area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
				dhy_[l]=1.0/aj[k][j][i]/area;

			       	area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				dhz_[l]=1.0/aj[k][j][i]/area;	
			}

	        }



	  	PetscBarrier(PETSC_NULL);

		MPI_Allreduce (&dhx_[0], &sum_dhx_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&dhy_[0], &sum_dhy_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&dhz_[0], &sum_dhz_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&count[0], &sum_count[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

        	       	ibm[ibi].dhx[l]=sum_dhx_[l]/sum_count[l];
	               	ibm[ibi].dhy[l]=sum_dhy_[l]/sum_count[l];
        	       	ibm[ibi].dhz[l]=sum_dhz_[l]/sum_count[l];

	               	dhx_[l]=sum_dhx_[l];
	               	dhy_[l]=sum_dhy_[l];
	               	dhz_[l]=sum_dhz_[l];

			iclose[l]=sum_iclose[l]/sum_count[l];
			jclose[l]=sum_jclose[l]/sum_count[l];
			kclose[l]=sum_kclose[l]/sum_count[l];
//			if (l==0) {
//		      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d iclose %d width %d \n",ibi, iclose[0], Itpwidth);
//		      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d jclose %d width %d \n",ibi, jclose[0], Itpwidth);
//		      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d kclose %d width %d \n",ibi, kclose[0], Itpwidth);
//			}
		}


	        int iii1 = 0;
		int iii2 = 0;
	        int iii = 0;

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

                        int n1e = ibm[ibi].nv1[i], n2e = ibm[ibi].nv2[i], n3e = ibm[ibi].nv3[i];

			/*
			double d12=sqrt(pow(ibm[ibi].x_bp[n1e]-ibm[ibi].x_bp[n2e],2)+pow(ibm[ibi].y_bp[n1e]-ibm[ibi].y_bp[n2e],2)+pow(ibm[ibi].z_bp[n1e]-ibm[ibi].z_bp[n2e],2));
			double d23=sqrt(pow(ibm[ibi].x_bp[n3e]-ibm[ibi].x_bp[n2e],2)+pow(ibm[ibi].y_bp[n3e]-ibm[ibi].y_bp[n2e],2)+pow(ibm[ibi].z_bp[n3e]-ibm[ibi].z_bp[n2e],2));
			double d31=sqrt(pow(ibm[ibi].x_bp[n1e]-ibm[ibi].x_bp[n3e],2)+pow(ibm[ibi].y_bp[n1e]-ibm[ibi].y_bp[n3e],2)+pow(ibm[ibi].z_bp[n1e]-ibm[ibi].z_bp[n3e],2));

			double d_max=max(d12,d23);
			d_max=max(d_max,d31);
			
			Itpwidth=(int)d_max/ibm[ibi].dhx[l]+4;
			*/

			int ic=iclose[l];
			int jc=jclose[l];
			int kc=kclose[l];

			int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
		
			if (ic1>lxe||ic2<lxs) {
				ibm[ibi].i_min[l]=lxs;
				ibm[ibi].i_max[l]=lxs;
			} else {
				ibm[ibi].i_min[l]=PetscMax(ic1, lxs);
				ibm[ibi].i_max[l]=PetscMin(ic2, lxe);
			}

			// Itpwidth=(int)d_max/ibm[ibi].dhy[l]+4;
			int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 

	       		if (jc1>lye||jc2<lys) {
				ibm[ibi].j_min[l]=lys;
				ibm[ibi].j_max[l]=lys;
			} else {
				ibm[ibi].j_min[l]=PetscMax(jc1, lys);
				ibm[ibi].j_max[l]=PetscMin(jc2, lye);
			}

			// Itpwidth=(int)d_max/ibm[ibi].dhz[l]+4;
			int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 

	       		if (kc1>lze||kc2<lzs) {
				ibm[ibi].k_min[l]=lzs;
				ibm[ibi].k_max[l]=lzs;
			} else {
				ibm[ibi].k_min[l]=PetscMax(kc1, lzs);
				ibm[ibi].k_max[l]=PetscMin(kc2, lze);
			}


		}

		// testttttt RTD

		/*
	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			ibm[ibi].i_min[l]=lxs;
			ibm[ibi].i_max[l]=lxe;

			ibm[ibi].j_min[l]=lys;
			ibm[ibi].j_max[l]=lye;


			ibm[ibi].k_min[l]=lzs;
			ibm[ibi].k_max[l]=lze;

		}
		*/


		/*
		if (rotor_model == 2) {

		int imin_g=1000000, imax_g=-1000000, jmin_g=1000000, jmax_g=-1000000, kmin_g=1000000, kmax_g=-1000000;

	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			if (imin_g>ibm[ibi].i_min[l]) imin_g=ibm[ibi].i_min[l];	
			if (imax_g<ibm[ibi].i_max[l]) imax_g=ibm[ibi].i_max[l];	

			if (jmin_g>ibm[ibi].j_min[l]) jmin_g=ibm[ibi].j_min[l];	
			if (jmax_g<ibm[ibi].j_max[l]) jmax_g=ibm[ibi].j_max[l];	

			if (kmin_g>ibm[ibi].k_min[l]) kmin_g=ibm[ibi].k_min[l];	
			if (kmax_g<ibm[ibi].k_max[l]) kmax_g=ibm[ibi].k_max[l];	
		}

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			ibm[ibi].i_min[l]=imin_g;	
			ibm[ibi].j_min[l]=jmin_g;	
			ibm[ibi].k_min[l]=kmin_g;	

			ibm[ibi].i_max[l]=imax_g;	
			ibm[ibi].j_max[l]=jmax_g;	
			ibm[ibi].k_max[l]=kmax_g;	

		}	

		}
		*/

		free(dhx_);
		free(dhy_); 
		free(dhz_);
	        free(iclose);
		free(jclose);
		free(kclose);

		free(sum_dhx_);
		free(sum_dhy_); 
		free(sum_dhz_);
	        free(sum_iclose);
		free(sum_jclose);
		free(sum_kclose);

	}

  	DAVecRestoreArray(fda, Coor, &coor);

	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

  	PetscBarrier(PETSC_NULL);

        PetscGetTime(&te1);  
      	PetscPrintf(PETSC_COMM_WORLD, "Time for rotor model search  %le\n", te1-ts1);

  	return(0);

}




// Call in main.c, after reading the grid for turbines
PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm, int NumberOfObjects)
{
	DA              da = user->da, fda = user->fda;
	DALocalInfo  	info;
	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt        mx, my, mz; // Dimensions in three directions
	PetscInt        i, j, k, l, ibi;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts		***coor, ***csi, ***eta, ***zet;
	PetscReal 	***aj;
	Vec		Coor;
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
	DAGetGhostedCoordinates(da, &Coor);
	DAVecGetArray(fda, Coor, &coor);

	DAVecGetArray(fda, user->lCsi,  &csi);
	DAVecGetArray(fda, user->lEta,  &eta);
	DAVecGetArray(fda, user->lZet,  &zet);
	DAVecGetArray(da,  user->lAj,  &aj);

	int n1e, n2e, n3e;

	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		int iclose[ibm[ibi].n_elmt], jclose[ibm[ibi].n_elmt], kclose[ibm[ibi].n_elmt]; 
		int sum_iclose[ibm[ibi].n_elmt], sum_jclose[ibm[ibi].n_elmt], sum_kclose[ibm[ibi].n_elmt]; 
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			iclose[l]=0;
			jclose[l]=0;
			kclose[l]=0;
			sum_iclose[l]=0;
			sum_jclose[l]=0;
			sum_kclose[l]=0;
		}
		int count[ibm[ibi].n_elmt], sum_count[ibm[ibi].n_elmt];
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			int imark=0;
			double dmin=1.e6;
			int ic,jc,kc;
			dmin=1.e6;
			for (i=lxs; i<lxe; i++){
				j=lys+1;
				k=lzs+1;
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					ic=i;
				}
			}
			dmin=1.e6;
			for (j=lys; j<lye; j++){
				k=lzs+1; 
				i=lxs+1;
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					jc=j;
				}
			}
			dmin=1.e6;
			for (k=lzs; k<lze; k++){
				i=lxs+1;
				j=lys+1;
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					kc=k;
				}
			}
			
			double r1=ibm[ibi].cent_x[l]-coor[kc][jc][ic].x, r2=ibm[ibi].cent_y[l]-coor[kc][jc][ic].y, r3=ibm[ibi].cent_z[l]-coor[kc][jc][ic].z; 
			dmin=sqrt(r1*r1+r2*r2+r3*r3);
			count[l]=1;	
			i=ic; j=jc; k=kc;
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dhx=1.0/aj[k][j][i]/area;
			area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dhy=1.0/aj[k][j][i]/area;
			area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double dhz=1.0/aj[k][j][i]/area;	
			double dh_grid = 2.0*sqrt(dhx*dhx+dhy*dhy+dhz*dhz);
			double diff=dmin-dh_grid;
			if (diff>1.e-6) {
				count[l]=0;	
				ic=0; jc=0; kc=0;
				iclose[l]=0; jclose[l]=0; kclose[l]=0;
			} 
			else {
				iclose[l]=ic; jclose[l]=jc; kclose[l]=kc;
				i=ic; j=jc; k=kc;
			}
		}
		PetscBarrier(PETSC_NULL);
		MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&count[0], &sum_count[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		for (l=0; l<ibm[ibi].n_elmt; l++) {
			iclose[l]=sum_iclose[l]/(sum_count[l]);
			jclose[l]=sum_jclose[l]/(sum_count[l]);
			kclose[l]=sum_kclose[l]/(sum_count[l]);
		}
		int iii1 = 0;
		int iii2 = 0;
		int iii = 0;
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			int n1e = ibm[ibi].nv1[i], n2e = ibm[ibi].nv2[i], n3e = ibm[ibi].nv3[i];
			int ic=iclose[l];
			int jc=jclose[l];
			int kc=kclose[l];
			int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
			if (ic1>lxe||ic2<lxs) {
				ibm[ibi].i_min[l]=lxs;
				ibm[ibi].i_max[l]=lxs;
			} 
			else {
				ibm[ibi].i_min[l]=PetscMax(ic1, lxs);
				ibm[ibi].i_max[l]=PetscMin(ic2, lxe);
			}
			int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 

	    if (jc1>lye||jc2<lys) {
				ibm[ibi].j_min[l]=lys;
				ibm[ibi].j_max[l]=lys;
			} 
			else {
				ibm[ibi].j_min[l]=PetscMax(jc1, lys);
				ibm[ibi].j_max[l]=PetscMin(jc2, lye);
			}
			int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 
			if (kc1>lze||kc2<lzs) {
				ibm[ibi].k_min[l]=lzs;
				ibm[ibi].k_max[l]=lzs;
			} 
			else {
				ibm[ibi].k_min[l]=PetscMax(kc1, lzs);
				ibm[ibi].k_max[l]=PetscMin(kc2, lze);
			}
		}
	}
	DAVecRestoreArray(fda, Coor, &coor);
	DAVecRestoreArray(fda, user->lCsi,  &csi);
	DAVecRestoreArray(fda, user->lEta,  &eta);
	DAVecRestoreArray(fda, user->lZet,  &zet);
	DAVecRestoreArray(da,  user->lAj,  &aj);
	return(0);
}



/* Read information of ACL, like chord, twist angle, lift and drag coefficients */
// call in the main.c for actuator line simulation only
PetscErrorCode airfoil_ACL(ACL *acl, IBMNodes *ibm,  FSInfo *fsi)
{
	int	rank;
	int	n_CD, n_CL;
	int	ifoil, i;
	char string[128];
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if(!rank) {
		FILE *fd;
		PetscPrintf(PETSC_COMM_SELF, "READ airfoil data for ACL\n");
		char filen[80]; 
		for(ifoil=0; ifoil<num_foiltype; ifoil++) {
			sprintf(filen,"%s/FOIL%2.2d" , path, ifoil);
			fd = fopen(filen, "r"); 
			if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open airfoil data file"), exit(0);
			if (fd) {
				fgets(string, 128, fd);
				fgets(string, 128, fd);
				fscanf(fd, "%i ",&(acl[ifoil].num_AC));
				MPI_Bcast(&(acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);
				PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].r_AC));  
				PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].angle_twistInp));  
				PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].chord_bladeInp)); 
			 	for(i=0; i<acl[ifoil].num_AC; i++) {
					fscanf(fd, "%le %le %le", &(acl[ifoil].r_AC[i]), &(acl[ifoil].angle_twistInp[i]), &(acl[ifoil].chord_bladeInp[i]));
					acl[ifoil].r_AC[i] = acl[ifoil].r_AC[i] / reflength_wt;
					acl[ifoil].chord_bladeInp[i] = acl[ifoil].chord_bladeInp[i] / reflength_wt;
				} 
				acl[ifoil].r_beg = acl[ifoil].r_AC[0];
				acl[ifoil].r_end = acl[ifoil].r_AC[acl[ifoil].num_AC-1];
				MPI_Bcast(acl[ifoil].r_AC, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(acl[ifoil].angle_twistInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(acl[ifoil].chord_bladeInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(&(acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(&(acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			}
			fclose(fd);
			sprintf(filen,"%s/CD%2.2d" , path, ifoil);
			fd = fopen(filen, "r"); 
			if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open CD file"), exit(0);
			if (fd) {
				fgets(string, 128, fd);
				fgets(string, 128, fd);
				fscanf(fd, "%i ",&(acl[ifoil].num_CD));
				MPI_Bcast(&(acl[ifoil].num_CD), 1, MPI_INT, 0, PETSC_COMM_WORLD);
				PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].ang_CD));
				PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].CDInp));
				for(i=0; i<acl[ifoil].num_CD; i++) {
					fscanf(fd, "%le %le", &(acl[ifoil].ang_CD[i]), &(acl[ifoil].CDInp[i]));
				}
				MPI_Bcast(acl[ifoil].ang_CD, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(acl[ifoil].CDInp, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
			}
			fclose(fd);
			sprintf(filen,"%s/CL%2.2d" , path, ifoil);
			fd = fopen(filen, "r");
			if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open CL file"), exit(0);
			if (fd) {
				fgets(string, 128, fd);
				fgets(string, 128, fd);
				fscanf(fd, "%i ",&(acl[ifoil].num_CL));
				MPI_Bcast(&(acl[ifoil].num_CL), 1, MPI_INT, 0, PETSC_COMM_WORLD);
				PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].ang_CL));
				PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].CLInp));
				for(i=0; i<acl[ifoil].num_CL; i++) {
					fscanf(fd, "%le %le", &(acl[ifoil].ang_CL[i]), &(acl[ifoil].CLInp[i]));
				}
				MPI_Bcast(acl[ifoil].ang_CL, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(acl[ifoil].CLInp, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
			}
			fclose(fd);
		}
	} 
	else if(rank) {
		for(ifoil=0; ifoil<num_foiltype; ifoil++) {
			// angle of twist , chord length of blade
			MPI_Bcast(&(acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);
			PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].r_AC));  
			PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].angle_twistInp));  
			PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].chord_bladeInp)); 
			MPI_Bcast(acl[ifoil].r_AC, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(acl[ifoil].angle_twistInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(acl[ifoil].chord_bladeInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(&(acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(&(acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			// drag coefficients
			MPI_Bcast(&(acl[ifoil].num_CD), 1, MPI_INT, 0, PETSC_COMM_WORLD);
			PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].ang_CD));
			PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].CDInp));
			MPI_Bcast(acl[ifoil].ang_CD, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(acl[ifoil].CDInp, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
			// left coefficient
			MPI_Bcast(&(acl[ifoil].num_CL), 1, MPI_INT, 0, PETSC_COMM_WORLD);
			PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].ang_CL));
			PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].CLInp));
			MPI_Bcast(acl[ifoil].ang_CL, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(acl[ifoil].CLInp, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
		}
	}
	int ibi, j;
	double r, fac1, fac2;
	for (ibi=0;ibi<NumberOfTurbines;ibi++) {
		for (i=0;i<ibm[ibi].n_elmt;i++) {
			r = sqrt( pow((ibm[ibi].cent_x[i]-fsi[ibi].x_c),2.0) + pow((ibm[ibi].cent_y[i]-fsi[ibi].y_c),2.0) + pow((ibm[ibi].cent_z[i]-fsi[ibi].z_c),2.0));
			for(ifoil=0; ifoil<num_foiltype; ifoil++) {
				if( r>=acl[ifoil].r_beg && r<=acl[ifoil].r_end ) {
					for(j=0; j<acl[ifoil].num_AC-1; j++) {
						if (r>=acl[ifoil].r_AC[j] && r<=acl[ifoil].r_AC[j+1]) {
							fac1 = (acl[ifoil].r_AC[j+1]-r)/(acl[ifoil].r_AC[j+1]-acl[ifoil].r_AC[j]);
							fac2 = (-acl[ifoil].r_AC[j]+r)/(acl[ifoil].r_AC[j+1]-acl[ifoil].r_AC[j]);
							ibm[ibi].angle_twist[i]	= fac1*acl[ifoil].angle_twistInp[j] + fac2*acl[ifoil].angle_twistInp[j+1];
							ibm[ibi].chord_blade[i]	= fac1*acl[ifoil].chord_bladeInp[j] + fac2*acl[ifoil].chord_bladeInp[j+1];
						}
					}
				}
			}
		}
	}
	PetscPrintf(PETSC_COMM_WORLD, "Finish READing airfoil data for ACL\n");
	return(0);
}

// calculate the force on the actuator disk, called in the solvers.c beforce solving the momentum equation
PetscErrorCode Calc_F_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

  	PetscInt	l, ibi;
  	double	pi = 3.141592653589793, a = 0.25;
  	double 	C_T;  // = 4.0 / 3.0;
  	double 	U_ref, A_sum;
  	PetscReal	sign;
  	double 	Uref_x, Uref_y, Uref_z;


  	Calc_U_lagr(user, ibm, fsi, NumberOfObjects);
  
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
	

		if (rotor_model==4) {
			C_T=ibm[ibi].CT;
		} else {

			double indf_ax=ibm[ibi].indf_axis;
		  	C_T = 4.0 * indf_ax * (1-indf_ax);
	  		C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );
		}


		if (rotor_model==4) U_ref=ibm[ibi].U_ref;
		else {
	    		U_ref = 0.0;
    			A_sum = 0.0;
    			Uref_x = 0.0;
	    		Uref_y = 0.0;
    			Uref_z = 0.0;

	    		for (l=0; l<ibm[ibi].n_elmt; l++) {
      				Uref_x += ibm[ibi].U_lagr_x[l] * ibm[ibi].dA[l];
      				Uref_y += ibm[ibi].U_lagr_y[l] * ibm[ibi].dA[l];
	      			Uref_z += ibm[ibi].U_lagr_z[l] * ibm[ibi].dA[l];
      				A_sum += ibm[ibi].dA[l];
	    		}
    			Uref_x /= A_sum;
    			Uref_y /= A_sum;
	    		Uref_z /= A_sum;

    			U_ref = Uref_x*fsi[ibi].nx_tb+Uref_y*fsi[ibi].ny_tb+Uref_z*fsi[ibi].nz_tb;
		}


    		PetscPrintf(PETSC_COMM_WORLD, "**** The U_ref at %i th body: %le\n", ibi, U_ref);

    		sign = U_ref / fabs(U_ref+1.e-9);
		//U_ref = 1.0; //RTD
    		//sign = 1.0;
    		for (l=0; l<ibm[ibi].n_elmt; l++) {
      			ibm[ibi].F_lagr_x[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].nx_tb;
      			ibm[ibi].F_lagr_y[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].ny_tb;
      			ibm[ibi].F_lagr_z[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].nz_tb;  
    		}


 	}

 	return(0);

}

// still working on 
// calculate the force on the actuator disk, called in the solvers.c beforce solving the momentum equation
PetscErrorCode Calc_F_lagr_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

  	PetscInt	l, ibi;
  	double	pi = 3.141592653589793, a = 0.25;
  	double 	C_T;  // = 4.0 / 3.0;
  	double 	U_ref, A_sum;
  	PetscReal	sign;
  	double 	Uref_x, Uref_y, Uref_z;
  	int 		nv1, nv2, nv3;
	double 		r1, r2, r3, dh, rx, ry, rz;

	double nfx, nfy, nfz;


	Cmpnts nr, na, nt; 
	Cmpnts Ut, Un, Ubt, Ubn, Ub;

  	Calc_U_lagr(user, ibm, fsi, NumberOfObjects);

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
	
		double Force_axial = 0.5*(ibm[ibi].axialprojectedarea/(reflength_wt*reflength_wt))*ibm[ibi].axialforcecoefficient*ibm[ibi].U_ref*ibm[ibi].U_ref;
		double Force_axial_fromfriction = 0.5*(ibm[ibi].tangentialprojectedarea/(reflength_wt*reflength_wt))*ibm[ibi].tangentialforcecoefficient*ibm[ibi].U_ref*ibm[ibi].U_ref;


		double sumforce_pressure=0.0, sumforce_friction=0.0;

		// calculate the friction and pressure  coefficients
    		for (l=0; l<ibm[ibi].n_elmt; l++) {

			nv1=ibm[ibi].nv1[l];
			nv2=ibm[ibi].nv2[l];
			nv3=ibm[ibi].nv3[l];


		       	rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1];
		       	ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1];
	          	rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

	              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
	                rx = ibm[ibi].x_bp[nv3] - ibm[ibi].x_bp[nv2];
	                ry = ibm[ibi].y_bp[nv3] - ibm[ibi].y_bp[nv2];
	                rz = ibm[ibi].z_bp[nv3] - ibm[ibi].z_bp[nv2];

	                r2 = sqrt(rx*rx + ry*ry + rz*rz );

	                rx = ibm[ibi].x_bp[nv1] - ibm[ibi].x_bp[nv3];
	                ry = ibm[ibi].y_bp[nv1] - ibm[ibi].y_bp[nv3];
	                rz = ibm[ibi].z_bp[nv1] - ibm[ibi].z_bp[nv3];

	                r3 = sqrt(rx*rx + ry*ry + rz*rz );

			dh = 3.0*(r1 + r2 + r3)/3.0;
	
			//dh = ibm[ibi].dh;



			nfx=ibm[ibi].nf_x[l];
			nfy=ibm[ibi].nf_y[l];
			nfz=ibm[ibi].nf_z[l];

			double UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 

			Un.x = UUn*nfx; 
			Un.y = UUn*nfy; 
			Un.z = UUn*nfz; 

			Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
			Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
			Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;


			Ub.x = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
			Ub.y = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
			Ub.z = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
			UUn = Ub.x*nfx +  Ub.y*nfy +  Ub.z*nfz; 

			Ubn.x = UUn*nfx; 
			Ubn.y = UUn*nfy; 
			Ubn.z = UUn*nfz; 

			Ubt.x = Ub.x - Ubn.x;
			Ubt.y = Ub.y - Ubn.y;
			Ubt.z = Ub.z - Ubn.z;

			double Urefx = Ubn.x-Un.x;
			double Urefy = Ubn.y-Un.y;
			double Urefz = Ubn.z-Un.z;

			double signx = Urefx/(fabs(Urefx)+1.e-11);
			double signy = Urefy/(fabs(Urefy)+1.e-11);
			double signz = Urefz/(fabs(Urefz)+1.e-11);

			//double fx, fy, fz;
	      		//fx =  ibm[ibi].dA[l]*dh*Urefx/user->dt;
	      		//fy =  ibm[ibi].dA[l]*dh*Urefy/user->dt;
	      		//fz =  ibm[ibi].dA[l]*dh*Urefz/user->dt;  


			double fx, fy, fz;
	      		fx =  0.5 * ibm[ibi].dA[l] * (Urefx * Urefx) * signx;
	      		fy =  0.5 * ibm[ibi].dA[l] * (Urefy * Urefy) * signy;
	      		fz =  0.5 * ibm[ibi].dA[l] * (Urefz * Urefz) * signz;  

			sumforce_pressure += fx*fsi[ibi].nx_tb + fy*fsi[ibi].ny_tb + fz*fsi[ibi].nz_tb;

			Urefx = Ubt.x-Ut.x;
			Urefy = Ubt.y-Ut.y;
			Urefz = Ubt.z-Ut.z;

			signx = Urefx/(fabs(Urefx)+1.e-11);
			signy = Urefy/(fabs(Urefy)+1.e-11);
			signz = Urefz/(fabs(Urefz)+1.e-11);
			
	      		fx = 0.5 * ibm[ibi].dA[l] * (Urefx * Urefx) * signx;
      			fy = 0.5 * ibm[ibi].dA[l] * (Urefy * Urefy) * signy;
      			fz = 0.5 * ibm[ibi].dA[l] * (Urefz * Urefz) * signz;  

			sumforce_friction += fx*fsi[ibi].nx_tb + fy*fsi[ibi].ny_tb + fz*fsi[ibi].nz_tb;
    		}


		if (sumforce_friction>1.e-19) ibm[ibi].friction_factor = 0.0; 
		ibm[ibi].friction_factor = fmax(0.0, Force_axial_fromfriction / (-sumforce_friction-1.e-19)); 

		if (sumforce_pressure>1.e-19) ibm[ibi].pressure_factor = 0.0; 
		else ibm[ibi].pressure_factor = fmax(0.0,(Force_axial - Force_axial_fromfriction) / (-sumforce_pressure-1.e-19)); 

		double sum_dh=0.0;
		sumforce_pressure=0.0; sumforce_friction=0.0;
    		for (l=0; l<ibm[ibi].n_elmt; l++) {

			nv1=ibm[ibi].nv1[l];
			nv2=ibm[ibi].nv2[l];
			nv3=ibm[ibi].nv3[l];

		       	rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1];
		       	ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1];
	          	rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

	              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
	                rx = ibm[ibi].x_bp[nv3] - ibm[ibi].x_bp[nv2];
	                ry = ibm[ibi].y_bp[nv3] - ibm[ibi].y_bp[nv2];
	                rz = ibm[ibi].z_bp[nv3] - ibm[ibi].z_bp[nv2];

	                r2 = sqrt(rx*rx + ry*ry + rz*rz );

	                rx = ibm[ibi].x_bp[nv1] - ibm[ibi].x_bp[nv3];
	                ry = ibm[ibi].y_bp[nv1] - ibm[ibi].y_bp[nv3];
	                rz = ibm[ibi].z_bp[nv1] - ibm[ibi].z_bp[nv3];

	                r3 = sqrt(rx*rx + ry*ry + rz*rz );

			dh = 3.0*(r1 + r2 + r3)/3.0;

			if (nacelle_model==2) dh = ibm[ibi].dh;

			sum_dh += dh;

			nfx=ibm[ibi].nf_x[l];
			nfy=ibm[ibi].nf_y[l];
			nfz=ibm[ibi].nf_z[l];


			double UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 

			Un.x = UUn*nfx; 
			Un.y = UUn*nfy; 
			Un.z = UUn*nfz; 

			Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
			Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
			Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;


			Ub.x = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
			Ub.y = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
			Ub.z = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
			UUn = Ub.x*nfx +  Ub.y*nfy +  Ub.z*nfz; 

			Ubn.x = UUn*nfx; 
			Ubn.y = UUn*nfy; 
			Ubn.z = UUn*nfz; 

			Ubt.x = Ub.x - Ubn.x;
			Ubt.y = Ub.y - Ubn.y;
			Ubt.z = Ub.z - Ubn.z;


			double Urefx = Ubn.x-Un.x;
			double Urefy = Ubn.y-Un.y;
			double Urefz = Ubn.z-Un.z;

			double signx = Urefx/(fabs(Urefx)+1.e-11);
			double signy = Urefy/(fabs(Urefy)+1.e-11);
			double signz = Urefz/(fabs(Urefz)+1.e-11);

			double Uref=ibm[ibi].U_ref;	

			
	      		//ibm[ibi].F_lagr_x[l] =  dh*Urefx/user->dt;
	      		//ibm[ibi].F_lagr_y[l] =  dh*Urefy/user->dt;
	      		//ibm[ibi].F_lagr_z[l] =  dh*Urefz/user->dt;  
			
			double fx, fy, fz;

			if (nacelle_model == 1 || nacelle_model == 2) {
		      		fx =  dh*Urefx/user->dt;
		      		fy =  dh*Urefy/user->dt;
	      			fz =  dh*Urefz/user->dt;  
			} else {
	      			fx =  0.5 * ibm[ibi].pressure_factor * (Urefx * Urefx) * signx;
	      			fy =  0.5 * ibm[ibi].pressure_factor * (Urefy * Urefy) * signy;
	      			fz =  0.5 * ibm[ibi].pressure_factor * (Urefz * Urefz) * signz;  
			}	

	      		ibm[ibi].F_lagr_x[l] = fx;
	      		ibm[ibi].F_lagr_y[l] = fy;
	      		ibm[ibi].F_lagr_z[l] = fz;  


			sumforce_pressure += (fx*fsi[ibi].nx_tb + fy*fsi[ibi].ny_tb + fz*fsi[ibi].nz_tb)*ibm[ibi].dA[l];


			Urefx = Ubt.x-Ut.x;
			Urefy = Ubt.y-Ut.y;
			Urefz = Ubt.z-Ut.z;

			signx = Urefx/(fabs(Urefx)+1.e-11);
			signy = Urefy/(fabs(Urefy)+1.e-11);
			signz = Urefz/(fabs(Urefz)+1.e-11);
			
	      		fx = 0.5 * ibm[ibi].friction_factor * (Urefx * Urefx) * signx;
      			fy = 0.5 * ibm[ibi].friction_factor * (Urefy * Urefy) * signy;
      			fz = 0.5 * ibm[ibi].friction_factor * (Urefz * Urefz) * signz;  

	     		ibm[ibi].F_lagr_x[l] += fx;
      			ibm[ibi].F_lagr_y[l] += fy;
      			ibm[ibi].F_lagr_z[l] += fz;  

			sumforce_friction += (fx*fsi[ibi].nx_tb + fy*fsi[ibi].ny_tb + fz*fsi[ibi].nz_tb)*ibm[ibi].dA[l];
    		}

	        if (!rank) {
		        FILE *f;
                	char filen[80];
	       	        sprintf(filen, "nacelleforcecoefficients%2.2d_%2.2d", nacelle_model, ibi);
   			if (ti==1) {
				f = fopen(filen, "w");
    				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=time axialforcecoefficient axialprojectedarea tangentialforcecoefficient tangentialprojectedarea frictioncoefficient pressurecoefficient Uref Force_axial_prescrbied sumaxialforce_pressure sumaxialforce_friction dh\n");
			} else f = fopen(filen, "a");

        	        PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le %le \n",ti, ibm[ibi].axialforcecoefficient, ibm[ibi].axialprojectedarea, ibm[ibi].tangentialforcecoefficient, ibm[ibi].tangentialprojectedarea, ibm[ibi].friction_factor, ibm[ibi].pressure_factor, ibm[ibi].U_ref, Force_axial, sumforce_pressure, sumforce_friction, sum_dh/(double)ibm[ibi].n_elmt );
        		fclose(f);

        	}
 	}

 	return(0);

}



// not for turbine simulations
PetscErrorCode Calc_Ftmprt_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

  	int	l, ibi;
  	int 		nv1, nv2, nv3;
	double 		r1, r2, r3, rr, rx, ry, rz;


        int rank=0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {
             	nv1 = ibm[ibi].nv1[l];
            	nv2 = ibm[ibi].nv2[l];
              	nv3 = ibm[ibi].nv3[l];

               	rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1];
             	ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1];
             	rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
                rx = ibm[ibi].x_bp[nv3] - ibm[ibi].x_bp[nv2];
                ry = ibm[ibi].y_bp[nv3] - ibm[ibi].y_bp[nv2];
                rz = ibm[ibi].z_bp[nv3] - ibm[ibi].z_bp[nv2];

                r2 = sqrt(rx*rx + ry*ry + rz*rz );

                rx = ibm[ibi].x_bp[nv1] - ibm[ibi].x_bp[nv3];
                ry = ibm[ibi].y_bp[nv1] - ibm[ibi].y_bp[nv3];
                rz = ibm[ibi].z_bp[nv1] - ibm[ibi].z_bp[nv3];

                r3 = sqrt(rx*rx + ry*ry + rz*rz );

		double tmprtb;
		if (rotor_model == 3 || rotor_model == 6) {
			rr = r1*r1;
			tmprtb = (ibm[ibi].tmprt[nv1]+ibm[ibi].tmprt[nv2])/2.0;
		} else {
			rr = (r1 + r2 + r3)/3.0;
			tmprtb = (ibm[ibi].tmprt[nv1]+ibm[ibi].tmprt[nv2]+ibm[ibi].tmprt[nv3])/3.0;
		}
	
      		ibm[ibi].Ftmprt_lagr[l] = rr * (tmprtb - ibm[ibi].Tmprt_lagr[l])/user->dt;    

    	}

  	}
	
          
 return(0);

}





// not for turbine simulations
PetscErrorCode Calc_F_lagr_noslip(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

  	int	l, ibi;
  	int 		nv1, nv2, nv3;
	double 		r1, r2, r3, rr, rx, ry, rz;


        int rank=0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  	Calc_U_lagr(user, ibm, fsi, NumberOfObjects);

	double U_ref;
 
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {
             	nv1 = ibm[ibi].nv1[l];
            	nv2 = ibm[ibi].nv2[l];
              	nv3 = ibm[ibi].nv3[l];

               	rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1];
             	ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1];
             	rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
                rx = ibm[ibi].x_bp[nv3] - ibm[ibi].x_bp[nv2];
                ry = ibm[ibi].y_bp[nv3] - ibm[ibi].y_bp[nv2];
                rz = ibm[ibi].z_bp[nv3] - ibm[ibi].z_bp[nv2];

                r2 = sqrt(rx*rx + ry*ry + rz*rz );

                rx = ibm[ibi].x_bp[nv1] - ibm[ibi].x_bp[nv3];
                ry = ibm[ibi].y_bp[nv1] - ibm[ibi].y_bp[nv3];
                rz = ibm[ibi].z_bp[nv1] - ibm[ibi].z_bp[nv3];

                r3 = sqrt(rx*rx + ry*ry + rz*rz );

		rr = (r1 + r2 + r3)/3.0;
	
//		if (rotor_model) {
//			if (rotor_model==2 || rotor_model==3 || rotor_model==4 || rotor_model == 6) U_ref=ibm[ibi].U_ref;
//			if (rotor_model==1) U_ref=refvel_cfd;

			/*
			double indf_ax=1.0;
		  	double C_T = 4.0 * indf_ax * (1-indf_ax);
	  		double C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );
			*/



			// U_ref=ibm[ibi].U_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].U_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].U_lagr_z[l]*fsi[ibi].nz_tb;
			//double U_ref=0.6;

//			double CD=ibm[ibi].CD_bluff;
  //  			double sign = U_ref / (fabs(U_ref)+1.e-20);
      			//ibm[ibi].F_lagr_x[l]=-0.5*CD*(U_ref*U_ref)*sign*fsi[ibi].nx_tb/(L_nacelle/reflength_wt);
      			//ibm[ibi].F_lagr_y[l]=-0.5*CD*(U_ref*U_ref)*sign*fsi[ibi].ny_tb/(L_nacelle/reflength_wt);
      			//ibm[ibi].F_lagr_z[l]=-0.5*CD*(U_ref*U_ref)*sign*fsi[ibi].nz_tb/(L_nacelle/reflength_wt); 

			
//			ibm[ibi].F_lagr_x[l] = rr * (ibm[ibi].u[l].x - ibm[ibi].U_lagr_x[l])/user->dt;
//	      		ibm[ibi].F_lagr_y[l] = rr * (ibm[ibi].u[l].y - ibm[ibi].U_lagr_y[l])/user->dt;
  //    			ibm[ibi].F_lagr_z[l] = rr * (ibm[ibi].u[l].z - ibm[ibi].U_lagr_z[l])/user->dt;    
			

//		} else 
		{

			double ux = (ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x+ibm[ibi].u[nv3].x)/3.0;
			double uy = (ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y+ibm[ibi].u[nv3].y)/3.0;
			double uz = (ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z+ibm[ibi].u[nv3].z)/3.0;
      			ibm[ibi].F_lagr_x[l] = rr * (ux - ibm[ibi].U_lagr_x[l])/user->dt;
	      		ibm[ibi].F_lagr_y[l] = rr * (uy - ibm[ibi].U_lagr_y[l])/user->dt;
      			ibm[ibi].F_lagr_z[l] = rr * (uz - ibm[ibi].U_lagr_z[l])/user->dt;    

      			//ibm[ibi].F_lagr_x[l] = ux;
	      		//ibm[ibi].F_lagr_y[l] = uy;
      			//ibm[ibi].F_lagr_z[l] = uz;    

		}

    	}

  	}
	
          	for(ibi=0;ibi<NumberOfObjects;ibi++) PetscPrintf(PETSC_COMM_WORLD, "IBDelta %d Uref %le CD %le Lag_Force Fx %le Fy %le Fz %le\n", ibi, U_ref, ibm[ibi].CD_bluff, ibm[ibi].F_lagr_x[0], ibm[ibi].F_lagr_y[0], ibm[ibi].F_lagr_z[0]);
//          	for(ibi=0;ibi<NumberOfObjects;ibi++) PetscPrintf(PETSC_COMM_WORLD, "IBDelta %d Lag_Ux %le\n", ibi, ibm[ibi].U_lagr_x[0]);
 //         	for(ibi=0;ibi<NumberOfObjects;ibi++) PetscPrintf(PETSC_COMM_WORLD, "IBDelta %d Lag_Uy %le\n", ibi, ibm[ibi].U_lagr_y[0]);
         	for(ibi=0;ibi<NumberOfObjects;ibi++) PetscPrintf(PETSC_COMM_WORLD, "IBDelta %d Lag_Uz %le\n", ibi, ibm[ibi].U_lagr_z[0]);
 return(0);

}



// calculate force at largrangian points using actuator line model
//
// called in solvers.c before solving the momentum equation
PetscErrorCode Calc_F_lagr_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, ACL *acl, int NumberOfObjects)
{
	PetscInt      	l, ibi, j;
	double		pi = 3.141592653589793;
	PetscReal		A_sum, U_ref, rr, rx, ry, rz, tmp, u_tangent;
	Cmpnts	      	U_b, n_relV, n_L, n_blade, n_rot;	
	int		nv1, nv2, nv3, ifoil;
	PetscReal 	fac1, fac2, r;
	int istall;
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//Calc_U_lagr(user, ibm, fsi, NumberOfObjects);
	double sign_L;
	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			nv1 = ibm[ibi].nv1[l];
			nv2 = ibm[ibi].nv2[l];
			nv3 = ibm[ibi].nv3[l];
			rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
			ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
			rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;
			double r = sqrt(rx*rx+ry*ry+rz*rz)+1.e-20;
			n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;
			double ux, uy, uz;
			if (rotor_model == 3 || rotor_model == 5 || rotor_model == 6) {
				ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
				uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
				uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 
			} 
			else if (rotor_model == 2) {
				double fac=1.0/3.0;
				ux=fac*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x+ibm[ibi].u[nv3].x); 
				uy=fac*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y+ibm[ibi].u[nv3].y); 
				uz=fac*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z+ibm[ibi].u[nv3].z); 
			}			
			double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));
			n_rot.x=ux/(Ublade+1.e-20);
			n_rot.y=uy/(Ublade+1.e-20);
			n_rot.z=uz/(Ublade+1.e-20);
			U_b.x=ibm[ibi].U_lagr_x[l]-ux;
			U_b.y=ibm[ibi].U_lagr_y[l]-uy;
			U_b.z=ibm[ibi].U_lagr_z[l]-uz;
			double Ur=U_b.x*n_blade.x+U_b.y*n_blade.y+U_b.z*n_blade.z;
			U_b.x-=Ur*n_blade.x;
			U_b.y-=Ur*n_blade.y;
			U_b.z-=Ur*n_blade.z;
			U_ref=sqrt(pow(U_b.x,2)+pow(U_b.y,2)+pow(U_b.z,2));
			n_relV.x = U_b.x/(U_ref+1.e-20); n_relV.y = U_b.y/(U_ref+1.e-20); n_relV.z = U_b.z/(U_ref+1.e-20); 
			tmp = -n_relV.x*n_rot.x-n_relV.y*n_rot.y-n_relV.z*n_rot.z;
			tmp /= (1+1e-9);
			ibm[ibi].angle_attack[l] = acos(tmp) * 180.0 / pi - ibm[ibi].angle_twist[l] - ibm[ibi].pitch[0];
			n_L.x = n_relV.y*n_blade.z - n_relV.z*n_blade.y; 
			n_L.y = n_relV.z*n_blade.x - n_relV.x*n_blade.z; 
			n_L.z = n_relV.x*n_blade.y - n_relV.y*n_blade.x;
			tmp = fsi[ibi].nx_tb*n_L.x + fsi[ibi].ny_tb*n_L.y + fsi[ibi].nz_tb*n_L.z;
			if (tmp < 0.0) n_L.x = -n_L.x, n_L.y = -n_L.y, n_L.z = -n_L.z;   
			for(ifoil=0; ifoil<num_foiltype; ifoil++) {
				if( r>=acl[ifoil].r_beg && r<=acl[ifoil].r_end ) {
					int inrange=0;
					for(j=0; j<acl[ifoil].num_CD-1; j++) {
						if (ibm[ibi].angle_attack[l]>=acl[ifoil].ang_CD[j] && ibm[ibi].angle_attack[l]<=acl[ifoil].ang_CD[j+1]) {
							fac1 = (acl[ifoil].ang_CD[j+1]-ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CD[j+1]-acl[ifoil].ang_CD[j]);
							fac2 = (-acl[ifoil].ang_CD[j]+ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CD[j+1]-acl[ifoil].ang_CD[j]);
							ibm[ibi].CD[l] = fac1*acl[ifoil].CDInp[j] + fac2*acl[ifoil].CDInp[j+1];
							inrange=1;
						}
					}
					if (!inrange) {
						if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CD[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CD[l] = (-1.2 - acl[ifoil].CDInp[0]) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CD[0]) / (-45.0-acl[ifoil].ang_CD[0]) + acl[ifoil].CDInp[0];
						if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CD[l] = -1.2;
						if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CD[l] = -1.2;
						if (ibm[ibi].angle_attack[l]>acl[ifoil].ang_CD[acl[ifoil].num_CD-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CD[l] = (1.2 - acl[ifoil].CDInp[acl[ifoil].num_CD-1]) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CD[acl[ifoil].num_CD-1]) / (45.0-acl[ifoil].ang_CD[acl[ifoil].num_CD-1]) + acl[ifoil].CDInp[acl[ifoil].num_CD-1];
						if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CD[l] = 1.2;
					}
					inrange=0;
					for(j=0; j<acl[ifoil].num_CL-1; j++) {
						if (ibm[ibi].angle_attack[l]>=acl[ifoil].ang_CL[j] && ibm[ibi].angle_attack[l]<=acl[ifoil].ang_CL[j+1]) {
							fac1 = (acl[ifoil].ang_CL[j+1]-ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CL[j+1]-acl[ifoil].ang_CL[j]);
							fac2 = (-acl[ifoil].ang_CL[j]+ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CL[j+1]-acl[ifoil].ang_CL[j]);
							ibm[ibi].CL[l] = fac1*acl[ifoil].CLInp[j] + fac2*acl[ifoil].CLInp[j+1];
							inrange=1;
						}
					}
					if (!inrange) {
						if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CL[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CL[l] = (-1.05 - acl[ifoil].CLInp[0] ) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CL[0] ) / (-45.0-acl[ifoil].ang_CL[0] ) + acl[ifoil].CLInp[0] ;
						if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);
						if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CL[l] = 0.0;
						if (ibm[ibi].angle_attack[l]>acl[ifoil].ang_CL[acl[ifoil].num_CL-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CL[l] = (1.05 - acl[ifoil].CLInp[acl[ifoil].num_CL-1] ) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CL[acl[ifoil].num_CL-1] ) / (45.0-acl[ifoil].ang_CL[acl[ifoil].num_CL-1] ) + acl[ifoil].CLInp[acl[ifoil].num_CL-1] ;
						if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);
					}
				}
			}
			// add 3D and rotational effects  // Chariaropoulos PK and Hansen MOL, JFE 2000:122:330-6
			if (correction3D_CH) {
				double  coef_cr1 = c1_CH, coef_cr2 = c2_CH, coef_cr3 = c3_CH;
				double angle_twist = ibm[ibi].angle_twist[l]*pi/180;
				double angle_AOA = ibm[ibi].angle_attack[l]*pi/180;
				double angle_inflow = angle_twist + angle_AOA;
				ibm[ibi].CD[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * pow(cos(angle_twist), coef_cr3) * ( ibm[ibi].CD[l] - 0.001 );
				ibm[ibi].CL[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * pow(cos(angle_twist), coef_cr3) * ( 2.0*pi*pi*fabs(ibm[ibi].angle_attack[l])/180.0 - ibm[ibi].CL[l]);
			} 
			double factor;
			if (rotor_model==2) {
				factor=1.0;
			} 
			else if (rotor_model == 3 || rotor_model == 5 || rotor_model == 6) {
				factor=ibm[ibi].chord_blade[l];
			}
			if (ti>tistart+1) ibm[ibi].Urel_mean[l] += U_ref;
			ibm[ibi].U_rel[l] = U_ref;
			// force from C_D
			ibm[ibi].F_lagr_x[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * n_relV.x * factor;
			ibm[ibi].F_lagr_y[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * n_relV.y * factor;
			ibm[ibi].F_lagr_z[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * n_relV.z * factor; 
			// add force from C_L
			ibm[ibi].F_lagr_x[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * n_L.x * factor;
			ibm[ibi].F_lagr_y[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * n_L.y * factor;
			ibm[ibi].F_lagr_z[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * n_L.z * factor;
			// Modification to 2D force
			double F1=1.0;
			// Shen Model 
			if (Shen_AL) {
				double R=fsi[ibi].r_rotor/reflength_wt;
				double pi = 3.1415926;
				double phi1 = (ibm[ibi].angle_attack[l] + ibm[ibi].angle_twist[l] + ibm[ibi].pitch[0])*pi/180.0;
				double Coef_a = a_shen; //0.125;
				double Coef_b = b_shen; //21;
				double Coef_c = c_shen; //0.1;
				double g1 = exp(-Coef_a*(num_blade*fabs(ibm[ibi].Tipspeedratio)-Coef_b)) + Coef_c; 
				//g1=g1_AL;
				F1 = 2.0*acos(exp(-g1*num_blade*(max(R-r,0.0))/2.0/r/fabs(sin(phi1))))/pi;
			}
			// Prandtl model 
			if (Prandtl_AL) {
				double R=fsi[ibi].r_rotor/reflength_wt;
				double pi = 3.1415926;
				double g1 = sqrt(1+pow(ibm[ibi].Tipspeedratio,2)); 
				F1 = 2.0*acos(exp(-g1*num_blade*(max(R-r,0.0))/2.0/R))/pi;
			}
			// Shen1_AL is calibrated for 2.5 MW clipper turbine 
			if (Shen1_AL) {
				double R=fsi[ibi].r_rotor/reflength_wt;
				double pi = 3.1415926;
				double phi1 = (ibm[ibi].angle_attack[l] + ibm[ibi].angle_twist[l] + ibm[ibi].pitch[0])*pi/180.0;
				//double Coef_a = 0.0879; //0.00373;
				//double Coef_b = -2.3738; //-0.61076;
				//double Coef_c = 12.1076; //3.2461;
				double Lambda = fabs(ibm[ibi].Tipspeedratio);
				double g1;

				if (Lambda<7.0636) {
					g1 = -1.7636*Lambda+11.6112;
				} 
				else if (Lambda>7.0636 && Lambda<9.8146) {
					g1 = -0.5805*Lambda+3.2542;
				} 
				else {
					g1 = -0.5077*Lambda+2.5397;
				}
				g1*=correction_ALShen1;
				g1 = exp(g1);
				F1 = 2.0*acos(exp(-g1*num_blade*(max(R-r,0.0))/2.0/r/fabs(sin(phi1))))/pi;
			}
			ibm[ibi].F_lagr_x[l] *= F1;
			ibm[ibi].F_lagr_y[l] *= F1;
			ibm[ibi].F_lagr_z[l] *= F1;
			if (correction3D_CL) {
				if (r<0.75*fsi[ibi].r_rotor/reflength_wt) {
					double angle_twist = ibm[ibi].angle_twist[l]*pi/180;
					double angle_AOA = ibm[ibi].angle_attack[l]*pi/180;
					double angle_inflow = angle_twist + angle_AOA;
					double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;
					double F_correc = -0.5*c0_CL * pow( ibm[ibi].chord_blade[l]/r, 1) * pow(cos(angle_inflow), 2) * U_ref * U_ref * factor;
					ibm[ibi].F_lagr_x[l]+=F_correc*nx;	
					ibm[ibi].F_lagr_y[l]+=F_correc*ny;	
					ibm[ibi].F_lagr_z[l]+=F_correc*nz;	
				}
			}
		}
	}
	// smooth the force in the blended region	
	if (smoothforce_AL) {
		int i,j,k;
		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			int Elmt_blade=ibm[ibi].n_elmt/num_blade;
			double tmp_forcex[ibm[ibi].n_elmt];
			double tmp_forcey[ibm[ibi].n_elmt];
			double tmp_forcez[ibm[ibi].n_elmt];
			int ii=0;
			do {
				ii++;
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					tmp_forcex[l]=ibm[ibi].F_lagr_x[l];
					tmp_forcey[l]=ibm[ibi].F_lagr_y[l];
					tmp_forcez[l]=ibm[ibi].F_lagr_z[l];
				}
				for (k=0; k<num_blade; k++) {
					int j1=k*Elmt_blade+4;
					int j2=(k+1)*Elmt_blade-4;
					for (j=j1; j<j2; j++) {
						double sum_forcex=0.0;
						double sum_forcey=0.0;
						double sum_forcez=0.0;
						for (i=j-4; i<j+4; i++) {
							double weight=(double)(i-j);
							sum_forcex+=tmp_forcex[i]*dfunc_s4h(weight);
							sum_forcey+=tmp_forcey[i]*dfunc_s4h(weight);
							sum_forcez+=tmp_forcez[i]*dfunc_s4h(weight);
						}
						ibm[ibi].F_lagr_x[j]=sum_forcex;
						ibm[ibi].F_lagr_y[j]=sum_forcey;
						ibm[ibi].F_lagr_z[j]=sum_forcez;
					}	
				}
			} while(ii<=smoothforce_AL);
		}
	}
	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			r = sqrt( pow((ibm[ibi].cent_x[l]-fsi[ibi].x_c),2) + pow((ibm[ibi].cent_y[l]-fsi[ibi].y_c),2) + pow((ibm[ibi].cent_z[l]-fsi[ibi].z_c),2));
			if (r<r_nacelle/reflength_wt) {
				ibm[ibi].F_lagr_x[l]=0.0;	
				ibm[ibi].F_lagr_y[l]=0.0;	
				ibm[ibi].F_lagr_z[l]=0.0;	
			}
		}
	}
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (ti>tistart+1) {
		PetscPrintf(PETSC_COMM_WORLD, "AL: export\n");
		count_AL += 1.0;
		FILE *f;
		char filen[80];
		double fac = 1.0 / count_AL;
		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			if (!rank) {
				sprintf(filen, "ForceAOA@ACL%2.2d", ibi);
				f = fopen(filen, "w");
				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = \"r\", \"Ft\", \"Fa\", \"AOA\", \"C\", \"Urel\"  \n");
			}
			for (l=0; l<ibm[ibi].n_elmt; l++) {
				nv1 = ibm[ibi].nv1[l];
				nv2 = ibm[ibi].nv2[l];
				nv3 = ibm[ibi].nv3[l];
				rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
				ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
				rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;
				double r = sqrt(rx*rx+ry*ry+rz*rz);
				n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;
				double ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
				double uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
				double uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 
				double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));
				n_rot.x=ux/(Ublade+1.e-20);
				n_rot.y=uy/(Ublade+1.e-20);
				n_rot.z=uz/(Ublade+1.e-20);
				//ibm[ibi].Fr_mean[l] += ibm[ibi].F_lagr_x[l]*n_blade.x+ibm[ibi].F_lagr_y[l]*n_blade.y+ibm[ibi].F_lagr_z[l]*n_blade.z;
				ibm[ibi].Ft_mean[l] += ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z;
				ibm[ibi].Fa_mean[l] += ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb;
				//ibm[ibi].Ur_mean[l] += ibm[ibi].U_lagr_x[l]*n_blade.x+ibm[ibi].U_lagr_y[l]*n_blade.y+ibm[ibi].U_lagr_z[l]*n_blade.z;
				//ibm[ibi].Ut_mean[l] += ibm[ibi].U_lagr_x[l]*n_rot.x+ibm[ibi].U_lagr_y[l]*n_rot.y+ibm[ibi].U_lagr_z[l]*n_rot.z;
				//ibm[ibi].Ua_mean[l] += ibm[ibi].U_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].U_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].U_lagr_z[l]*fsi[ibi].nz_tb;
				ibm[ibi].AOA_mean[l] += ibm[ibi].angle_attack[l];
				//double Fr = -ibm[ibi].Fr_mean[l]*fac;
				double Ft = -ibm[ibi].Ft_mean[l]*fac;
				double Fa = -ibm[ibi].Fa_mean[l]*fac;
				//double Ur = ibm[ibi].Ur_mean[l]*fac;
				//double Ut = ibm[ibi].Ut_mean[l]*fac;
				//double Ua = ibm[ibi].Ua_mean[l]*fac;
				double AOA = ibm[ibi].AOA_mean[l]*fac;
				if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le \n", r, Ft, Fa, AOA, ibm[ibi].chord_blade[l], ibm[ibi].Urel_mean[l]*fac);
			}
			if (!rank) fclose(f);
		}
	}

	if (ti==tistart+1 || ti == (ti/tiout) * tiout) {
		FILE *f;
		char filen[80];
		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			if (!rank) {
				sprintf(filen, "ForceAOA@ACL%06d_%2.2d", ti, ibi);
				f = fopen(filen, "w");
				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = \"r\", \"Ft\", \"Fa\", \"AOA\", \"C\", \"Urel\"  \n");
			}
			for (l=0; l<ibm[ibi].n_elmt; l++) {
				nv1 = ibm[ibi].nv1[l];
				nv2 = ibm[ibi].nv2[l];
				nv3 = ibm[ibi].nv3[l];
				rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
				ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
				rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;
				double r = sqrt(rx*rx+ry*ry+rz*rz);
				n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;
				double ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
				double uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
				double uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 
				double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));
				n_rot.x=ux/(Ublade+1.e-20);
				n_rot.y=uy/(Ublade+1.e-20);
				n_rot.z=uz/(Ublade+1.e-20);
				double Ft = -(ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z);
				double Fa = -(ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb);
				double AOA = ibm[ibi].angle_attack[l];
				if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le \n", r, Ft, Fa, AOA, ibm[ibi].chord_blade[l], ibm[ibi].U_rel[l]);
			}
			if (!rank) fclose(f);
		}
	}
	PetscPrintf(PETSC_COMM_WORLD, "AL: lagrange force finished \n");
	return(0);
}




/* Interpolate the velocity at the Lagrangian points */
// subroutine for Calc_F_lagr_ACL and Calc_F_lagr
PetscErrorCode Calc_Tmprt_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

  	DA              da = user->da, fda = user->fda;
  	DALocalInfo     info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts	***ucat, ***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj;

	PetscReal 	***tmprt;

  	Vec		Coor;

  	PetscReal	dfunc;

  	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

  	double r1, r2, r3;

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


  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);
  	DAVecGetArray(da, user->lTmprt, &tmprt);
  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);

  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
    	for (l=0; l<ibm[ibi].n_elmt; l++) {
      		ibm[ibi].Tmprt_lagr[l] = 0.0;
    	}
  	}


        clock_t start, end;
        double elapsed;

        start = clock();
	double ni[3], nj[3], nk[3];

	for (ibi=0; ibi<NumberOfObjects; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {

	    	for (k = ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
	      	for (j = ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
	        for (i = ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {


        		double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dhx_=1.0/aj[k][j][i]/area;

			area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dhy_=1.0/aj[k][j][i]/area;

			area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double dhz_=1.0/aj[k][j][i]/area;	



			xc=coor[k][j][i].x;
			yc=coor[k][j][i].y;
			zc=coor[k][j][i].z;

	            	double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);

			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

			r1=(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
			r2=(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
			r3=(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 


	                dfunc =  dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
	            	ibm[ibi].Tmprt_lagr[l] += tmprt[k][j][i] * dfunc;

      		}

  	}
  	}

        end = clock();
        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

        PetscPrintf(PETSC_COMM_WORLD, "Time for Tmprt_larg local %le \n", elapsed);



  	PetscBarrier(PETSC_NULL);

        start = clock();
  	for (ibi=0; ibi<NumberOfObjects; ibi++)  {

		std::vector<double> tmprt_local (ibm[ibi].n_elmt); 
		std::vector<double> tmprt_sum (ibm[ibi].n_elmt); 
	
  		for (i=0; i<ibm[ibi].n_elmt; i++ ) {
			tmprt_local[i] = ibm[ibi].Tmprt_lagr[i];
  		}

	  	MPI_Allreduce( &(tmprt_local[0]), &(tmprt_sum[0]), ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

	  	for (i=0; i<ibm[ibi].n_elmt; i++ ) {
        		ibm[ibi].Tmprt_lagr[i] = tmprt_sum[i];
	  	}

	  	PetscBarrier(PETSC_NULL);
	}

        end = clock();
        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

        PetscPrintf(PETSC_COMM_WORLD, "Time for U_larg global sum %le \n", elapsed);
  	PetscBarrier(PETSC_NULL);

        PetscPrintf(PETSC_COMM_WORLD, "Tmprt_lagr at ibi %d, i %d = %le \n", 0,0,ibm[0].Tmprt_lagr[0]);

  	if (ii_periodicWT || jj_periodicWT || kk_periodicWT) {

    		double xc_min = 1.0e6;        
 	   	double yc_min = 1.0e6;
	    	double zc_min = 1.0e6;

  		for (ibi=0; ibi<NumberOfObjects; ibi++) {
      			xc_min = PetscMin(xc_min, fsi[ibi].x_c);
	      		yc_min = PetscMin(yc_min, fsi[ibi].y_c);
      			zc_min = PetscMin(zc_min, fsi[ibi].z_c);
	    	}

//        	PetscPrintf(PETSC_COMM_WORLD, "HERE 1 \n");
		int IndWT[Nz_WT][Ny_WT][Nx_WT];

		double fac_x=1.0/Sx_WT, fac_y=1.0/Sy_WT, fac_z=1.0/Sz_WT;
        	for (ibi=0; ibi<NumberOfObjects; ibi++) {
			int ii=(int)((fsi[ibi].x_c-xc_min+1.e-9)*fac_x);
			int jj=(int)((fsi[ibi].y_c-yc_min+1.e-9)*fac_y);
			int kk=(int)((fsi[ibi].z_c-zc_min+1.e-9)*fac_z);

        		PetscPrintf(PETSC_COMM_WORLD, "ibi ii jj kk %d %d %d %d \n", ibi, ii, jj, kk);
			IndWT[kk][jj][ii]=ibi;
		}

		int i, j, k;

//        	PetscPrintf(PETSC_COMM_WORLD, "HERE 2 \n");
		for (k=0;k<Nz_WT;k++)
		for (j=0;j<Ny_WT;j++)
		for (i=0;i<Nx_WT;i++){

			if (ii_periodicWT && i==0) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[k][j][Nx_WT-1];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].Tmprt_lagr[l]=ibm[ibi].Tmprt_lagr[l]+ibm[ibi1].Tmprt_lagr[l];
				}
			}

                        if (jj_periodicWT && j==0) {
                                ibi=IndWT[k][j][i];
                                int ibi1=IndWT[k][Ny_WT-1][i];
                                for (l=0; l<ibm[ibi].n_elmt; l++) {
                                        ibm[ibi].Tmprt_lagr[l]=ibm[ibi].Tmprt_lagr[l]+ibm[ibi1].Tmprt_lagr[l];
                                }
                        }

                        if (kk_periodicWT && k==0) {
                                ibi=IndWT[k][j][i];
                                int ibi1=IndWT[Nz_WT-1][j][i];
                                for (l=0; l<ibm[ibi].n_elmt; l++) {
                                        ibm[ibi].Tmprt_lagr[l]=ibm[ibi].Tmprt_lagr[l]+ibm[ibi1].Tmprt_lagr[l];
                                }
                        }
		}


//        	PetscPrintf(PETSC_COMM_WORLD, "HERE 3 \n");
                for (k=0;k<Nz_WT;k++)
                for (j=0;j<Ny_WT;j++)
                for (i=0;i<Nx_WT;i++){

                        if (ii_periodicWT && i==Nx_WT-1) {
                                ibi=IndWT[k][j][i];
                                int ibi1=IndWT[k][j][0];
                                for (l=0; l<ibm[ibi].n_elmt; l++) {
                                        ibm[ibi].Tmprt_lagr[l]=ibm[ibi1].Tmprt_lagr[l];
                                }
                        }

                        if (jj_periodicWT && j==Ny_WT-1) {
                                ibi=IndWT[k][j][i];
                                int ibi1=IndWT[k][0][i];
                                for (l=0; l<ibm[ibi].n_elmt; l++) {
                                        ibm[ibi].Tmprt_lagr[l]=ibm[ibi1].Tmprt_lagr[l];
                                }
                        }

                        if (kk_periodicWT && k==Nz_WT-1) {
                                ibi=IndWT[k][j][i];
                                int ibi1=IndWT[0][j][i];
                                for (l=0; l<ibm[ibi].n_elmt; l++) {
                                        ibm[ibi].Tmprt_lagr[l]=ibm[ibi1].Tmprt_lagr[l];
                                }
                        }
                }


	}


  	DAVecRestoreArray(fda, Coor, &coor);
  	DAVecRestoreArray(da, user->lTmprt, &tmprt);
  	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

	return(0);
}






/* Interpolate the velocity at the Lagrangian points */
// subroutine for Calc_F_lagr_ACL and Calc_F_lagr
PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{
	DA              da = user->da, fda = user->fda;
	DALocalInfo     info;
	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt        mx, my, mz; // Dimensions in three directions
	PetscInt        i, j, k, l, ibi;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts	***ucat, ***coor, ***csi, ***eta, ***zet;
	PetscReal 	***aj;
	Vec		Coor;
	PetscReal	dfunc;
	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
	double r1, r2, r3;
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
	DAGetGhostedCoordinates(da, &Coor);
	DAVecGetArray(fda, Coor, &coor);
	DAVecGetArray(fda, user->lUcat, &ucat);
	DAVecGetArray(fda, user->lCsi,  &csi);
	DAVecGetArray(fda, user->lEta,  &eta);
	DAVecGetArray(fda, user->lZet,  &zet);
	DAVecGetArray(da,  user->lAj,  &aj);

	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			ibm[ibi].U_lagr_x[l] = 0.0;
			ibm[ibi].U_lagr_y[l] = 0.0;
			ibm[ibi].U_lagr_z[l] = 0.0;
		}
	}
	clock_t start, end;
	double elapsed;
	start = clock();
	double ni[3], nj[3], nk[3];
	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			for (k = ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
			for (j = ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
			for (i = ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {
				xc = 0.125 *
				(coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
				coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
				coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
				coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
				yc = 0.125 *
				(coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
				coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
				coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
				coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
				zc = 0.125 *
				(coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
				coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
				coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
				coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);
				double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
				double dhx_=1.0/aj[k][j][i]/area;
				area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
				double dhy_=1.0/aj[k][j][i]/area;
				area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				double dhz_=1.0/aj[k][j][i]/area;	
				double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);
				Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
				r1=(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
				r2=(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
				r3=(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 
				dfunc =  dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
				ibm[ibi].U_lagr_x[l] += ucat[k][j][i].x * dfunc;
				ibm[ibi].U_lagr_y[l] += ucat[k][j][i].y * dfunc;
				ibm[ibi].U_lagr_z[l] += ucat[k][j][i].z * dfunc;          
			}
		}
	}

	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	PetscPrintf(PETSC_COMM_WORLD, "Time for U_larg local %le \n", elapsed);
	PetscBarrier(PETSC_NULL);
	start = clock();
	for (ibi=0; ibi<NumberOfObjects; ibi++)  {
		std::vector<Cmpnts> u_local (ibm[ibi].n_elmt); 
		std::vector<Cmpnts> u_sum (ibm[ibi].n_elmt); 
		for (i=0; i<ibm[ibi].n_elmt; i++ ) {
			u_local[i].x = ibm[ibi].U_lagr_x[i];
			u_local[i].y = ibm[ibi].U_lagr_y[i];
			u_local[i].z = ibm[ibi].U_lagr_z[i];
		}
		MPI_Allreduce( &(u_local[0]), &(u_sum[0]), 3*ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		for (i=0; i<ibm[ibi].n_elmt; i++ ) {
			ibm[ibi].U_lagr_x[i] = u_sum[i].x;
			ibm[ibi].U_lagr_y[i] = u_sum[i].y;
			ibm[ibi].U_lagr_z[i] = u_sum[i].z;
		}
		PetscBarrier(PETSC_NULL);
	}
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	PetscPrintf(PETSC_COMM_WORLD, "Time for U_larg global sum %le \n", elapsed);
	PetscBarrier(PETSC_NULL);

	if (ii_periodicWT || jj_periodicWT || kk_periodicWT) {
		double xc_min = 1.0e6;        
		double yc_min = 1.0e6;
		double zc_min = 1.0e6;
		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			xc_min = PetscMin(xc_min, fsi[ibi].x_c);
			yc_min = PetscMin(yc_min, fsi[ibi].y_c);
			zc_min = PetscMin(zc_min, fsi[ibi].z_c);
		}
		int IndWT[Nz_WT][Ny_WT][Nx_WT];
		double fac_x=1.0/Sx_WT, fac_y=1.0/Sy_WT, fac_z=1.0/Sz_WT;
		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			int ii=(int)((fsi[ibi].x_c-xc_min+1.e-9)*fac_x);
			int jj=(int)((fsi[ibi].y_c-yc_min+1.e-9)*fac_y);
			int kk=(int)((fsi[ibi].z_c-zc_min+1.e-9)*fac_z);
			PetscPrintf(PETSC_COMM_WORLD, "ibi ii jj kk %d %d %d %d \n", ibi, ii, jj, kk);
			IndWT[kk][jj][ii]=ibi;
		}
		int i, j, k;
		for (k=0;k<Nz_WT;k++)
		for (j=0;j<Ny_WT;j++)
		for (i=0;i<Nx_WT;i++){
			if (ii_periodicWT && i==0) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[k][j][Nx_WT-1];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
				ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
				ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
				ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
				}
			}
			if (jj_periodicWT && j==0) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[k][Ny_WT-1][i];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
					ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
					ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
				}
			}
			if (kk_periodicWT && k==0) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[Nz_WT-1][j][i];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
					ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
					ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
				}
			}
		}
		for (k=0;k<Nz_WT;k++)
		for (j=0;j<Ny_WT;j++)
		for (i=0;i<Nx_WT;i++){
			if (ii_periodicWT && i==Nx_WT-1) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[k][j][0];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
					ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
					ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
				}
			}
			if (jj_periodicWT && j==Ny_WT-1) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[k][0][i];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
					ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
					ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
				}
			}
			if (kk_periodicWT && k==Nz_WT-1) {
				ibi=IndWT[k][j][i];
				int ibi1=IndWT[0][j][i];
				for (l=0; l<ibm[ibi].n_elmt; l++) {
					ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
					ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
					ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
				}
			}
		}
	}
	if (MoveFrame) {
		for (ibi=0; ibi<NumberOfObjects; ibi++)
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			ibm[ibi].U_lagr_x[l] += u_frame;
			ibm[ibi].U_lagr_y[l] += v_frame;
			ibm[ibi].U_lagr_z[l] += w_frame;
		}
	}
	DAVecRestoreArray(fda, Coor, &coor);
	DAVecRestoreArray(fda, user->lUcat, &ucat);
	DAVecRestoreArray(fda, user->lCsi,  &csi);
	DAVecRestoreArray(fda, user->lEta,  &eta);
	DAVecRestoreArray(fda, user->lZet,  &zet);
	DAVecRestoreArray(da,  user->lAj,  &aj);
	return(0);
}

// distribute the force on the turbines to background grid.
//
// called after Calc_F_lagr_ACL or Calc_F_lagr before solving the momentum equatiion
PetscErrorCode Calc_F_eul_old(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects, double dh, int df)
{

  	DA              da = user->da, fda = user->fda;
  	DALocalInfo     info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***lf_eul, ***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj, ***nvert;

  	Vec		Coor;

  	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
  
  	double		dfunc;

  	double r1, r2, r3;

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


  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);
  	DAVecGetArray(fda, user->lF_eul, &lf_eul);
  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);

	DAVecGetArray(da, user->lNvert, &nvert);

  	double ni[3], nj[3], nk[3];
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		/*
      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d i_min %d i_max %d \n",ibi, ibm[ibi].i_min[0], ibm[ibi].i_max[0]);
      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d j_min %d j_max %d \n",ibi, ibm[ibi].j_min[0], ibm[ibi].j_max[0]);
      		PetscPrintf(PETSC_COMM_WORLD, "ibi %d k_min %d k_max %d \n",ibi, ibm[ibi].k_min[0], ibm[ibi].k_max[0]);
		*/
  	for (l=0; l<ibm[ibi].n_elmt; l++) {

	    	for (k=ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
	      	for (j=ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
	        for (i=ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {


//                xc = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;
//                yc = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;
//                zc = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25;

			double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	
			double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	
			double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

			double xj = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;	
			double yj = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;	
			double zj = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25; 

			double xk = (coor[k  ][j  ][i].x + coor[k][j-1][i].x + coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;	
			double yk = (coor[k  ][j  ][i].y + coor[k][j-1][i].y + coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;	
			double zk = (coor[k  ][j  ][i].z + coor[k][j-1][i].z + coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25; 

			
			xc=coor[k][j][i].x;
			yc=coor[k][j][i].y;
			zc=coor[k][j][i].z;

			
			//xi=xc; xj=xc; xk=xc;
			//yi=yc; yj=yc; yk=yc;
			//zi=zc; zj=zc; zk=zc;
		

	            	double rxi= (xi - ibm[ibi].cent_x[l]), ryi = (yi - ibm[ibi].cent_y[l]), rzi = (zi - ibm[ibi].cent_z[l]);
	            	double rxj= (xj - ibm[ibi].cent_x[l]), ryj = (yj - ibm[ibi].cent_y[l]), rzj = (zj - ibm[ibi].cent_z[l]);
	            	double rxk= (xk - ibm[ibi].cent_x[l]), ryk = (yk - ibm[ibi].cent_y[l]), rzk = (zk - ibm[ibi].cent_z[l]);

			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

			//double dhx_=ibm[ibi].dhx[l];
			//double dhy_=ibm[ibi].dhy[l];
			//double dhz_=ibm[ibi].dhz[l];


        		double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dhx_=1.0/aj[k][j][i]/area;

			area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dhy_=1.0/aj[k][j][i]/area;

			area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double dhz_=1.0/aj[k][j][i]/area;	





			double r1i=fabs(rxi*ni[0]+ryi*ni[1]+rzi*ni[2])/dhx_; 
			double r2i=fabs(rxi*nj[0]+ryi*nj[1]+rzi*nj[2])/dhy_; 
			double r3i=fabs(rxi*nk[0]+ryi*nk[1]+rzi*nk[2])/dhz_; 

			double r1j=fabs(rxj*ni[0]+ryj*ni[1]+rzj*ni[2])/dhx_; 
			double r2j=fabs(rxj*nj[0]+ryj*nj[1]+rzj*nj[2])/dhy_; 
			double r3j=fabs(rxj*nk[0]+ryj*nk[1]+rzj*nk[2])/dhz_; 

			double r1k=fabs(rxk*ni[0]+ryk*ni[1]+rzk*ni[2])/dhx_; 
			double r2k=fabs(rxk*nj[0]+ryk*nj[1]+rzk*nj[2])/dhy_; 
			double r3k=fabs(rxk*nk[0]+ryk*nk[1]+rzk*nk[2])/dhz_; 



	          	vol_eul = aj[k][j][i];

			/*
			if (deltafunc == 1) dfunc = vol_eul * dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);
			if (deltafunc == 2) dfunc = vol_eul * dfunc_4h(r1) * dfunc_4h(r2) * dfunc_4h(r3);
			if (deltafunc == 3) dfunc = vol_eul * dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
			if (deltafunc == 4) dfunc = vol_eul * dfunc_s4h(r1) * dfunc_s4h(r2) * dfunc_s4h(r3);
			*/

			//dfunc = vol_eul * dfunc_nh(r1,3) * dfunc_nh(r2,3) * dfunc_nh(r3,3);
			// dfunc = vol_eul * dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
			//
			double dfunci, dfuncj, dfunck;
			if (df == -2) {
				dfunci = vol_eul * dfunc_4h2peak(r1i) * dfunc_4h2peak(r2i) * dfunc_4h2peak(r3i);
				dfuncj = vol_eul * dfunc_4h2peak(r1j) * dfunc_4h2peak(r2j) * dfunc_4h2peak(r3j);
				dfunck = vol_eul * dfunc_4h2peak(r1k) * dfunc_4h2peak(r2k) * dfunc_4h2peak(r3k);
			} else if (df == -1) {
				dfunci = vol_eul * dfunc_h(r1i) * dfunc_h(r2i) * dfunc_h(r3i);
				dfuncj = vol_eul * dfunc_h(r1j) * dfunc_h(r2j) * dfunc_h(r3j);
				dfunck = vol_eul * dfunc_h(r1k) * dfunc_h(r2k) * dfunc_h(r3k);
			} else if (df == 0) {
				dfunci = vol_eul * dfunc_2h(r1i) * dfunc_2h(r2i) * dfunc_2h(r3i);
				dfuncj = vol_eul * dfunc_2h(r1j) * dfunc_2h(r2j) * dfunc_2h(r3j);
				dfunck = vol_eul * dfunc_2h(r1k) * dfunc_2h(r2k) * dfunc_2h(r3k);
			} else if (df == 1) {
				dfunci = vol_eul * dfunc_2hs1(r1i) * dfunc_2hs1(r2i) * dfunc_2hs1(r3i);
				dfuncj = vol_eul * dfunc_2hs1(r1j) * dfunc_2hs1(r2j) * dfunc_2hs1(r3j);
				dfunck = vol_eul * dfunc_2hs1(r1k) * dfunc_2hs1(r2k) * dfunc_2hs1(r3k);
			} else if (df == 2) {			
				dfunci = vol_eul * dfunc_2hs2(r1i) * dfunc_2hs2(r2i) * dfunc_2hs2(r3i);
				dfuncj = vol_eul * dfunc_2hs2(r1j) * dfunc_2hs2(r2j) * dfunc_2hs2(r3j);
				dfunck = vol_eul * dfunc_2hs2(r1k) * dfunc_2hs2(r2k) * dfunc_2hs2(r3k);
			} else if (df == 3) {			
				dfunci = vol_eul * dfunc_2hs3(r1i) * dfunc_2hs3(r2i) * dfunc_2hs3(r3i);
				dfuncj = vol_eul * dfunc_2hs3(r1j) * dfunc_2hs3(r2j) * dfunc_2hs3(r3j);
				dfunck = vol_eul * dfunc_2hs3(r1k) * dfunc_2hs3(r2k) * dfunc_2hs3(r3k);
			} else if (df == 4) {			
				dfunci = vol_eul * dfunc_2hs4(r1i) * dfunc_2hs4(r2i) * dfunc_2hs4(r3i);
				dfuncj = vol_eul * dfunc_2hs4(r1j) * dfunc_2hs4(r2j) * dfunc_2hs4(r3j);
				dfunck = vol_eul * dfunc_2hs4(r1k) * dfunc_2hs4(r2k) * dfunc_2hs4(r3k);
			} else if (df == 5) {			
				double n = halfwidth_dfunc;
				dfunci = vol_eul * dfunc_nhs1(r1i,n) * dfunc_nhs1(r2i,n) * dfunc_nhs1(r3i,n);
				dfuncj = vol_eul * dfunc_nhs1(r1j,n) * dfunc_nhs1(r2j,n) * dfunc_nhs1(r3j,n);
				dfunck = vol_eul * dfunc_nhs1(r1k,n) * dfunc_nhs1(r2k,n) * dfunc_nhs1(r3k,n);
			} else if (df == 6) {		
				double n = halfwidth_dfunc;
				dfunci = vol_eul * dfunc_nhs2(r1i,n) * dfunc_nhs2(r2i,n) * dfunc_nhs2(r3i,n);
				dfuncj = vol_eul * dfunc_nhs2(r1j,n) * dfunc_nhs2(r2j,n) * dfunc_nhs2(r3j,n);
				dfunck = vol_eul * dfunc_nhs2(r1k,n) * dfunc_nhs2(r2k,n) * dfunc_nhs2(r3k,n);
			} else if (df == 7) {		
				double n = halfwidth_dfunc;
				dfunci = vol_eul * dfunc_exp(r1i,n) * dfunc_exp(r2i,n) * dfunc_exp(r3i,n);
				dfuncj = vol_eul * dfunc_exp(r1j,n) * dfunc_exp(r2j,n) * dfunc_exp(r3j,n);
				dfunck = vol_eul * dfunc_exp(r1k,n) * dfunc_exp(r2k,n) * dfunc_exp(r3k,n);
			} else { 
				dfunci = vol_eul * dfunc_s3h(r1i) * dfunc_s3h(r2i) * dfunc_s3h(r3i);
				dfuncj = vol_eul * dfunc_s3h(r1j) * dfunc_s3h(r2j) * dfunc_s3h(r3j);
				dfunck = vol_eul * dfunc_s3h(r1k) * dfunc_s3h(r2k) * dfunc_s3h(r3k);
			}

			//dfunc = vol_eul * dfunc_nh(r1, dfunc_wd) * dfunc_nh(r2, dfunc_wd) * dfunc_nh(r3, dfunc_wd);

			//dfunc = vol_eul * dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);

	            	lf_eul[k][j][i].x += ibm[ibi].F_lagr_x[l] * dfunci * dh * ibm[ibi].dA[l] * csi[k][j][i].x +
	                                     ibm[ibi].F_lagr_y[l] * dfunci * dh * ibm[ibi].dA[l] * csi[k][j][i].y +
                                 	     ibm[ibi].F_lagr_z[l] * dfunci * dh * ibm[ibi].dA[l] * csi[k][j][i].z;

	            	lf_eul[k][j][i].y += ibm[ibi].F_lagr_x[l] * dfuncj * dh * ibm[ibi].dA[l] *  eta[k][j][i].x +
        	                             ibm[ibi].F_lagr_y[l] * dfuncj * dh * ibm[ibi].dA[l] *  eta[k][j][i].y +
                	                     ibm[ibi].F_lagr_z[l] * dfuncj * dh * ibm[ibi].dA[l] *  eta[k][j][i].z;

	              	lf_eul[k][j][i].z += ibm[ibi].F_lagr_x[l] * dfunck * dh * ibm[ibi].dA[l] *  zet[k][j][i].x +
        	                             ibm[ibi].F_lagr_y[l] * dfunck * dh * ibm[ibi].dA[l] *  zet[k][j][i].y +
                	                     ibm[ibi].F_lagr_z[l] * dfunck * dh * ibm[ibi].dA[l] *  zet[k][j][i].z;
       
//			ibm[ibi].F_lagr_x[l]=1.00022323;
//			lf_eul[k][j][i].x+=ibm[ibi].F_lagr_x[l] ; 
//			lf_eul[k][j][i].y+=ibm[ibi].F_lagr_x[l] ; 
//			lf_eul[k][j][i].z+=ibm[ibi].F_lagr_x[l] ; 
	      	}
	}
	}

//          	for(ibi=0;ibi<NumberOfObjects;ibi++) PetscPrintf(PETSC_COMM_WORLD, "eul IBDelta %d Lag_Force %le\n", ibi, ibm[ibi].F_lagr_z[0]);


	for (k=lzs;k<lze;k++) 
	for (j=lys;j<lye;j++) 
	for (i=lxs;i<lxe;i++) {

		int ii, jj, kk;
		double _nvert;
		_nvert = 0.0;

		for (kk=k-1;kk<k+2;kk++) 
		for (jj=j-1;jj<j+2;jj++) 
		for (ii=i-1;ii<i+2;ii++) {
			_nvert += nvert[kk][jj][ii];
		}

		//_nvert = nvert[k][j][i];
		if ( _nvert >0.1 ) { 
                        lf_eul[k][j][i].x=0.0;
                        lf_eul[k][j][i].y=0.0;
                        lf_eul[k][j][i].z=0.0;
		}


		if (i==1) {
			lf_eul[k][j][i-1].x=0.0;
			lf_eul[k][j][i-1].y=0.0;
			lf_eul[k][j][i-1].z=0.0;
		}

		if (j==1) {
			lf_eul[k][j-1][i].x=0.0;
			lf_eul[k][j-1][i].y=0.0;
			lf_eul[k][j-1][i].z=0.0;
		}

		if (k==1) {
			lf_eul[k-1][j][i].x=0.0;
			lf_eul[k-1][j][i].y=0.0;
			lf_eul[k-1][j][i].z=0.0;
		}

		if (i==mx-2) {
			lf_eul[k][j][i+1].x=0.0;
			lf_eul[k][j][i+1].y=0.0;
			lf_eul[k][j][i+1].z=0.0;
		}

		if (j==my-2) {
			lf_eul[k][j+1][i].x=0.0;
			lf_eul[k][j+1][i].y=0.0;
			lf_eul[k][j+1][i].z=0.0;
		}

		if (k==mz-2) {
			lf_eul[k+1][j][i].x=0.0;
			lf_eul[k+1][j][i].y=0.0;
			lf_eul[k+1][j][i].z=0.0;
		}



	}
  
 	DAVecRestoreArray(fda, Coor, &coor);
  	DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
  	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

	DAVecRestoreArray(da, user->lNvert, &nvert);

	DALocalToGlobal(fda, user->lF_eul, INSERT_VALUES, user->F_eul);


	/*
	DAVecGetArray(fda, user->lF_eul, &lf_eul);
	DAVecGetArray(da, user->lNvert, &nvert);
	for (k=lzs;k<lze;k++) 
	for (j=lys;j<lye;j++) 
	for (i=lxs;i<lxe;i++) {
		if (fabs(lf_eul[k][j][i].x)>1.e-9 || fabs(lf_eul[k][j][i].y)>1.e-9 || fabs(lf_eul[k][j][i].z)>1.e-9 ) {
			nvert[k][j][i]:
		}
	}
  	DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	*/

	
	//if (ti%tiout==0) {
	//	char fname[80];
	//	sprintf(fname,"F_eul%06d_", ti);
	//	TECIOOut_rhs(user, user->lF_eul, fname);
	//}

  return(0);
}



// distribute the force on the turbines to background grid.
//
// called after Calc_F_lagr_ACL or Calc_F_lagr before solving the momentum equatiion
PetscErrorCode Calc_Ftmprt_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df)
{


	//PetscReal ts, te;
	//int my_rank;
	//MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);


        //PetscGetTime(&ts);  // xiaolei


  	DA              da = user->da, fda = user->fda;
  	DALocalInfo     info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj, ***nvert;

	PetscReal 	***lftmprt_eul;

  	Vec		Coor;

  	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
  
  	double		dfunc;

  	double r1, r2, r3;

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


  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);
  	DAVecGetArray(da, user->lFtmprt_eul, &lftmprt_eul);
  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);

	DAVecGetArray(da, user->lNvert, &nvert);


  	double ni[3], nj[3], nk[3];
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {

	    	for (k=ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
	      	for (j=ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
	        for (i=ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {

			xc=coor[k][j][i].x;
			yc=coor[k][j][i].y;
			zc=coor[k][j][i].z;

	            	double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);

			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

        		double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double dhx_=1.0/aj[k][j][i]/area;

			area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double dhy_=1.0/aj[k][j][i]/area;

			area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			double dhz_=1.0/aj[k][j][i]/area;	

			double r1=fabs(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
			double r2=fabs(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
			double r3=fabs(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 

	          	vol_eul = aj[k][j][i];

			//
			double dfunc;
			if (df == 0) {
				dfunc = vol_eul * dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);
			} else if (df == 7) {		
				double n = halfwidth_dfunc;
				dfunc = vol_eul * dfunc_exp(r1,n) * dfunc_exp(r2,n) * dfunc_exp(r3,n);
			} else { 
				dfunc = vol_eul * dfunc_sc4h(r1) * dfunc_sc4h(r2) * dfunc_sc4h(r3);
			}


	            	lftmprt_eul[k][j][i] += ibm[ibi].Ftmprt_lagr[l] * dfunc * ibm[ibi].dA[l];

	      	}
	}
	}



	for (k=lzs;k<lze;k++) 
	for (j=lys;j<lye;j++) 
	for (i=lxs;i<lxe;i++) {

		int ii, jj, kk;
		double _nvert;
		_nvert = 0.0;

		for (kk=k-1;kk<k+2;kk++) 
		for (jj=j-1;jj<j+2;jj++) 
		for (ii=i-1;ii<i+2;ii++) {
			_nvert += nvert[kk][jj][ii];
		}

		//_nvert = nvert[k][j][i];
		if ( _nvert >2.9 ) { 
                        lftmprt_eul[k][j][i]=0.0;
		}



		if (i==1) {
			lftmprt_eul[k][j][i-1]=0.0;
		}

		if (j==1) {
			lftmprt_eul[k][j-1][i]=0.0;
		}

		if (k==1) {
			lftmprt_eul[k-1][j][i]=0.0;
		}

		if (i==mx-2) {
			lftmprt_eul[k][j][i+1]=0.0;
		}

		if (j==my-2) {
			lftmprt_eul[k][j+1][i]=0.0;
		}

		if (k==mz-2) {
			lftmprt_eul[k+1][j][i]=0.0;
		}



	}
  
 	DAVecRestoreArray(fda, Coor, &coor);
  	DAVecRestoreArray(da, user->lFtmprt_eul, &lftmprt_eul);
  	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

	DAVecRestoreArray(da, user->lNvert, &nvert);

	DALocalToLocalBegin(da, user->lFtmprt_eul, INSERT_VALUES, user->lFtmprt_eul);
	DALocalToLocalEnd(da, user->lFtmprt_eul, INSERT_VALUES, user->lFtmprt_eul);
	DALocalToGlobal(da, user->lFtmprt_eul, INSERT_VALUES, user->lFtmprt_eul);	

	//if (ti==tistart || ti%tiout==0) {
	//	char fname[80];
	//	sprintf(fname,"F_eul%06d_", ti);
	//	TECIOOut_rhs(user, user->lF_eul, fname);
	//}

  return(0);
}






// distribute the force on the turbines to background grid.
//
// called after Calc_F_lagr_ACL or Calc_F_lagr before solving the momentum equatiion
PetscErrorCode Calc_F_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df)
{


	//PetscReal ts, te;
	//int my_rank;
	//MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);


        //PetscGetTime(&ts);  // xiaolei


  	DA              da = user->da, fda = user->fda;
  	DALocalInfo     info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***lf_eul, ***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj, ***nvert;

  	Vec		Coor;

  	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
  
  	double		dfunc;

  	double r1, r2, r3;

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


  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);
  	DAVecGetArray(fda, user->lF_eul, &lf_eul);
  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);

	DAVecGetArray(da, user->lNvert, &nvert);


  	double ni[3], nj[3], nk[3];
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {

	    	for (k=ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
	      	for (j=ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
	        for (i=ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {


			double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	
			double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	
			double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

			double xj = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;	
			double yj = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;	
			double zj = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25; 

			double xk = (coor[k  ][j  ][i].x + coor[k  ][j-1][i].x + coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;	
			double yk = (coor[k  ][j  ][i].y + coor[k  ][j-1][i].y + coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;	
			double zk = (coor[k  ][j  ][i].z + coor[k  ][j-1][i].z + coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25; 
			
			xc=coor[k][j][i].x;
			yc=coor[k][j][i].y;
			zc=coor[k][j][i].z;


			//xi=xc; xj=xc; xk=xc;
			//yi=yc; yj=yc; yk=yc;
			//zi=zc; zj=zc; zk=zc;
		

	            	double rxi= (xi - ibm[ibi].cent_x[l]), ryi = (yi - ibm[ibi].cent_y[l]), rzi = (zi - ibm[ibi].cent_z[l]);
	            	double rxj= (xj - ibm[ibi].cent_x[l]), ryj = (yj - ibm[ibi].cent_y[l]), rzj = (zj - ibm[ibi].cent_z[l]);
	            	double rxk= (xk - ibm[ibi].cent_x[l]), ryk = (yk - ibm[ibi].cent_y[l]), rzk = (zk - ibm[ibi].cent_z[l]);

			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			
			double dhx_, dhy_, dhz_;

			if (forcewidthfixed) {
				dhx_=dhi_fixed;
				dhy_=dhj_fixed;
				dhz_=dhk_fixed;
			} else {
	        		double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
				dhx_=1.0/aj[k][j][i]/area;

				area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
				dhy_=1.0/aj[k][j][i]/area;

				area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				dhz_=1.0/aj[k][j][i]/area;	
			}


			double r1i=fabs(rxi*ni[0]+ryi*ni[1]+rzi*ni[2])/dhx_; 
			double r2i=fabs(rxi*nj[0]+ryi*nj[1]+rzi*nj[2])/dhy_; 
			double r3i=fabs(rxi*nk[0]+ryi*nk[1]+rzi*nk[2])/dhz_; 

			double r1j=fabs(rxj*ni[0]+ryj*ni[1]+rzj*ni[2])/dhx_; 
			double r2j=fabs(rxj*nj[0]+ryj*nj[1]+rzj*nj[2])/dhy_; 
			double r3j=fabs(rxj*nk[0]+ryj*nk[1]+rzj*nk[2])/dhz_; 

			double r1k=fabs(rxk*ni[0]+ryk*ni[1]+rzk*ni[2])/dhx_; 
			double r2k=fabs(rxk*nj[0]+ryk*nj[1]+rzk*nj[2])/dhy_; 
			double r3k=fabs(rxk*nk[0]+ryk*nk[1]+rzk*nk[2])/dhz_; 


	          	//vol_eul = aj[k][j][i];
	          	vol_eul = 1.0/(dhx_*dhy_*dhz_);

			//
			double dfunci, dfuncj, dfunck;
			if (df == 0) {
				dfunci = vol_eul * dfunc_2h(r1i) * dfunc_2h(r2i) * dfunc_2h(r3i);
				dfuncj = vol_eul * dfunc_2h(r1j) * dfunc_2h(r2j) * dfunc_2h(r3j);
				dfunck = vol_eul * dfunc_2h(r1k) * dfunc_2h(r2k) * dfunc_2h(r3k);
			} else if (df == 7) {		
				double n = halfwidth_dfunc;
				dfunci = vol_eul * dfunc_exp(r1i,n) * dfunc_exp(r2i,n) * dfunc_exp(r3i,n);
				dfuncj = vol_eul * dfunc_exp(r1j,n) * dfunc_exp(r2j,n) * dfunc_exp(r3j,n);
				dfunck = vol_eul * dfunc_exp(r1k,n) * dfunc_exp(r2k,n) * dfunc_exp(r3k,n);
			} else { 
				dfunci = vol_eul * dfunc_s4h(r1i) * dfunc_s4h(r2i) * dfunc_s4h(r3i);
				dfuncj = vol_eul * dfunc_s4h(r1j) * dfunc_s4h(r2j) * dfunc_s4h(r3j);
				dfunck = vol_eul * dfunc_s4h(r1k) * dfunc_s4h(r2k) * dfunc_s4h(r3k);
			}


	            	lf_eul[k][j][i].x += ibm[ibi].F_lagr_x[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].x +
	                                     ibm[ibi].F_lagr_y[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].y +
                                 	     ibm[ibi].F_lagr_z[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].z;

	            	lf_eul[k][j][i].y += ibm[ibi].F_lagr_x[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].x +
        	                             ibm[ibi].F_lagr_y[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].y +
                	                     ibm[ibi].F_lagr_z[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].z;

	              	lf_eul[k][j][i].z += ibm[ibi].F_lagr_x[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].x +
        	                             ibm[ibi].F_lagr_y[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].y +
                	                     ibm[ibi].F_lagr_z[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].z;
	      	}
	}
	}



	for (k=lzs;k<lze;k++) 
	for (j=lys;j<lye;j++) 
	for (i=lxs;i<lxe;i++) {

		int ii, jj, kk;
		double _nvert;
		_nvert = 0.0;

		for (kk=k-1;kk<k+2;kk++) 
		for (jj=j-1;jj<j+2;jj++) 
		for (ii=i-1;ii<i+2;ii++) {
			_nvert += nvert[kk][jj][ii];
		}

		//_nvert = nvert[k][j][i];
		if ( _nvert >2.9 ) { 
                        lf_eul[k][j][i].x=0.0;
                        lf_eul[k][j][i].y=0.0;
                        lf_eul[k][j][i].z=0.0;
		}



		if (i==1) {
			lf_eul[k][j][i-1].x=0.0;
			lf_eul[k][j][i-1].y=0.0;
			lf_eul[k][j][i-1].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;

		}

		if (j==1) {
			lf_eul[k][j-1][i].x=0.0;
			lf_eul[k][j-1][i].y=0.0;
			lf_eul[k][j-1][i].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;

			// tower
			//lf_eul[k][j+1][i].x=0.0;
			//lf_eul[k][j+1][i].y=0.0;
			//lf_eul[k][j+1][i].z=0.0;

		}

		if (k==1) {
			lf_eul[k-1][j][i].x=0.0;
			lf_eul[k-1][j][i].y=0.0;
			lf_eul[k-1][j][i].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;
		}

		if (i==mx-2) {
			lf_eul[k][j][i+1].x=0.0;
			lf_eul[k][j][i+1].y=0.0;
			lf_eul[k][j][i+1].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;

		}

		if (j==my-2) {
			lf_eul[k][j+1][i].x=0.0;
			lf_eul[k][j+1][i].y=0.0;
			lf_eul[k][j+1][i].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;
		}

		if (k==mz-2) {
			lf_eul[k+1][j][i].x=0.0;
			lf_eul[k+1][j][i].y=0.0;
			lf_eul[k+1][j][i].z=0.0;

			lf_eul[k][j][i].x=0.0;
			lf_eul[k][j][i].y=0.0;
			lf_eul[k][j][i].z=0.0;
		}



	}
  
 	DAVecRestoreArray(fda, Coor, &coor);
  	DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
  	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

	DAVecRestoreArray(da, user->lNvert, &nvert);


        //PetscGetTime(&te);  // xiaolei
      	//printf( "rank %d Time for F_eul  %le\n", my_rank, te-ts);

	DALocalToLocalBegin(fda, user->lF_eul, INSERT_VALUES, user->lF_eul);
	DALocalToLocalEnd(fda, user->lF_eul, INSERT_VALUES, user->lF_eul);
	DALocalToGlobal(fda, user->lF_eul, INSERT_VALUES, user->F_eul);
	//PetscReal ts, te;
	//PetscGetTime(&ts);
	//PetscGetTime(&te);


	
	//if (ti==tistart || ti%tiout==0) {
	//	char fname[80];
	//	sprintf(fname,"F_eul%06d_", ti);
	//	TECIOOut_rhs(user, user->lF_eul, fname);
	//}

  return(0);
}




// exporting the forces on the turbines actuator disk model
PetscErrorCode Calc_forces_rotor(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi, char fname[80], int NumberOfObjects)
{
  	PetscInt      	l, ibi;
  	double        	pi = 3.141592653589793;
  	PetscReal	A_Sum, F_Sum, P_Sum, U_Sum, rx, ry, rz, M_Sum;


  	for (ibi=0; ibi<NumberOfObjects; ibi++) {


		double C_T;	
		if (rotor_model==1) {	
		double indf_ax=ibm[ibi].indf_axis;
  
	  	C_T = 4.0 * indf_ax * (1-indf_ax);

	  	C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );
		}

		if (rotor_model==4) {
			C_T=ibm[ibi].CT;
		}



    		A_Sum = 0.0; P_Sum = 0.0; U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0;
		double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;

    		for (l=0; l<ibm[ibi].n_elmt; l++) {

      			A_Sum += ibm[ibi].dA[l] ;

			double F_axis=ibm[ibi].F_lagr_x[l]*nx+ibm[ibi].F_lagr_y[l]*ny+ibm[ibi].F_lagr_z[l]*nz;
			double U_axis=ibm[ibi].U_lagr_x[l]*nx+ibm[ibi].U_lagr_y[l]*ny+ibm[ibi].U_lagr_z[l]*nz;

      			F_Sum += F_axis*ibm[ibi].dA[l] ;
      			P_Sum += F_axis*ibm[ibi].dA[l]*U_axis;

	      		U_Sum += U_axis*ibm[ibi].dA[l] ;

      			rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;	
      			ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;	
      			rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;	


      			double M_x= ry*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] - rz*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l];  
      			double M_y= rz*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] - rx*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l];  
      			double M_z= rx*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] - ry*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l];  

			M_Sum+=M_x*nx+M_y*ny+M_z*nz;


    		}
   
		U_Sum=U_Sum/A_Sum;
		
	if(floating_turbine_case){
		double HubtoCG=152.34; //Distance from the hub to the CG
		fsi[0].Force_rotor_z=-F_Sum;
		fsi[0].Mom_rotor_x=-F_Sum*HubtoCG;
	}	
		
	    	int rank=0;
    		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	    	if (!rank) {
      			FILE *f;
      			char filen[80];
	      		sprintf(filen, "%s/%s_%2.2d",path, fname,ibi);
      			if (ti==1) {
				f = fopen(filen, "w");
      				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time\", \"CT\", \"F\", \"Ud\", \"P\", \"U_ref\", \"M_axis\" \n");
			} else f = fopen(filen, "a");

      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le \n", ti*user[0].dt, C_T,-F_Sum,U_Sum,-P_Sum, ibm[ibi].U_ref, M_Sum);
	      		fclose(f);
	    	}


        	if ( ti==tistart+1 || ti % tiout == 0 /*ti == (ti/tiout) * tiout */) {
		        int rank=0;
        		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
			int i;

			double dt = user->dt;
		        if (!rank) {
        		        FILE *f;
	                	char filen[80];
        		        sprintf(filen, "%s/%s%06d_%03d_nf.dat",path,fname, ti,ibi);

		                f = fopen(filen, "w");
	
				int n_v=ibm[ibi].n_v; 	
				int n_elmt=ibm[ibi].n_elmt; 	

		                PetscFPrintf(PETSC_COMM_WORLD, f, "TITLE = \"Actuator disk mesh\" \n ");
			    	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,fx,fy,fz,u,v,w,ubx,uby,ubz\n");
		    		PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-15]=CELLCENTERED)\n", n_v, n_elmt);
                		PetscFPrintf(PETSC_COMM_WORLD, f, "STRANDID=0.1 SOLUTIONTIME=%le \n", ((double)ti)*dt);

		    		for (i=0; i<n_v; i++) {
	      				PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].x_bp[i]);
		    		}
    				for (i=0; i<n_v; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].y_bp[i]);
		    		}
		    		for (i=0; i<n_v; i++) {	
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].z_bp[i]);
		    		}

		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].nf_x[i]);
		    		}
    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].nf_y[i]);
    				}
		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].nf_z[i]);
		    		}

		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].F_lagr_x[i]);
		    		}
    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].F_lagr_y[i]);
    				}
		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].F_lagr_z[i]);
		    		}

		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].U_lagr_x[i]);
		    		}
    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].U_lagr_y[i]);
    				}
		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm[ibi].U_lagr_z[i]);
		    		}

		    		for (i=0; i<n_elmt; i++) {
					int nv1 = ibm[ibi].nv1[i];
					int nv2 = ibm[ibi].nv2[i];
					int nv3 = ibm[ibi].nv3[i];

					double Ubx = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", Ubx);
		    		}

		    		for (i=0; i<n_elmt; i++) {
					int nv1 = ibm[ibi].nv1[i];
					int nv2 = ibm[ibi].nv2[i];
					int nv3 = ibm[ibi].nv3[i];

					double Uby = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", Uby);
		    		}


		    		for (i=0; i<n_elmt; i++) {
					int nv1 = ibm[ibi].nv1[i];
					int nv2 = ibm[ibi].nv2[i];
					int nv3 = ibm[ibi].nv3[i];

					double Ubz = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", Ubz);
		    		}

   				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm[ibi].nv1[i]+1, ibm[ibi].nv2[i]+1, ibm[ibi].nv3[i]+1);
		    		}
 

                		fclose(f);

		        }
        	}
	}

    	return(0);
}

// exporting the forces on the turbines, actuator line model
PetscErrorCode Calc_forces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi)
{
	PetscInt      l, ibi;
	double        pi = 3.141592653589793;
	PetscReal A_Sum = 0.0, P_Sum = 0.0, U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0, Fx_Sum=0.0, Fy_Sum=0.0, Fz_Sum=0.0;
	PetscReal F_z, P, r, rx, ry, rz; 
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	for (ibi=0; ibi<NumberOfTurbines; ibi++) {
		A_Sum = 0.0; P_Sum = 0.0; U_Sum = 0.0; F_Sum=0.0; M_Sum=0.0; Fx_Sum=0.0; Fy_Sum=0.0; Fz_Sum=0.0;
		double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;
		for (l=0; l<ibm[ibi].n_elmt; l++) {
			A_Sum += ibm[ibi].dA[l] ;
			double F_axis=ibm[ibi].F_lagr_x[l]*nx+ibm[ibi].F_lagr_y[l]*ny+ibm[ibi].F_lagr_z[l]*nz;
			double U_axis=ibm[ibi].U_lagr_x[l]*nx+ibm[ibi].U_lagr_y[l]*ny+ibm[ibi].U_lagr_z[l]*nz;
			F_Sum += F_axis*ibm[ibi].dA[l] ;
			P_Sum += F_axis*ibm[ibi].dA[l]*U_axis;
			U_Sum += U_axis*ibm[ibi].dA[l] ;		
			Fx_Sum += ibm[ibi].F_lagr_x[l]*ibm[ibi].dA[l] ;
			Fy_Sum += ibm[ibi].F_lagr_y[l]*ibm[ibi].dA[l] ;
			Fz_Sum += ibm[ibi].F_lagr_z[l]*ibm[ibi].dA[l] ;
			rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;	
			ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;	
			rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;	
			double M_x= ry*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] - rz*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l];  
			double M_y= rz*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] - rx*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l];  
			double M_z= rx*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] - ry*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l];  
			M_Sum+=M_x*nx+M_y*ny+M_z*nz;
		}
		U_Sum=U_Sum/A_Sum;
		fsi[ibi].Torque_aero=-M_Sum;   // 20140124 
		fsi[ibi].Force_axis=-F_Sum;   // 20140124 
		if(floating_turbine_case){
			double HubtoCG=152.34; //Distance from the hub to the CG
			fsi[ibi].Force_rotor_z=-F_Sum;
			fsi[ibi].Mom_rotor_x=-F_Sum*HubtoCG;
			PetscPrintf(PETSC_COMM_WORLD,"AD moment x:%f, AD Force z:%f \n",fsi[ibi].Mom_rotor_x,fsi[ibi].Force_rotor_z);
		}			
		double P_moment=M_Sum*fsi[ibi].angvel_axis;
		double R0=fsi[0].r_rotor/reflength_wt;
		double rho_cfd=1.0; // no levelset 	
		double T_wt=fsi[0].r_rotor/refvel_wt; // the 0 turbine cannot be in the other turbines' wakes.  
		double T_cfd=R0/refvel_cfd; // the 0 turbine cannot be in the other turbines' wakes.  
		double ratio_rho=rho_air/rho_cfd;
		double ratio_L=reflength_wt;	
		double ratio_T=T_wt/T_cfd;	
		double ratio_V=refvel_wt/refvel_cfd;
		double Time_real=ti*user->dt*ratio_T;	
		double ang_real=fsi[ibi].ang_axis;
		double angvel_real=fsi[ibi].angvel_axis/ratio_T;
		double Force_real=ratio_rho*fsi[ibi].Force_axis*pow(ratio_V,2)*pow(ratio_L,2);
		double Torquea_real=ratio_rho*fsi[ibi].Torque_aero*pow(ratio_V,2)*pow(ratio_L,3);
		double Torquec_real=fsi[ibi].Torque_generator;
		double Torquec_cfd=fsi[ibi].Torque_generator/(ratio_rho*pow(ratio_V,2)*pow(ratio_L,3));
		double Uref_real=ibm[ibi].U_ref*ratio_V;
		double Ud_real=U_Sum*ratio_V;
		if (!rank) {
			FILE *f;
			char filen[80];
			sprintf(filen, "Turbine_AL%2.2d_%2.2d",rotor_model,ibi);
			if (ti==1) {
				f = fopen(filen, "w");
				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time\", \"angle\", \"angvel_axis\", \"Force_axis\", \"Torque_fluid\", \"Torque_generator\", \"Uref\", \"Ud\" ,\"TSR\" \n");
			} 
			else f = fopen(filen, "a");
			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le \n",ti*user->dt, fsi[ibi].ang_axis, fsi[ibi].angvel_axis,fsi[ibi].Force_axis, fsi[ibi].Torque_aero, Torquec_cfd, ibm[ibi].U_ref, U_Sum, ibm[ibi].Tipspeedratio);
			fclose(f);
			sprintf(filen, "Turbine_AL_real%2.2d_%2.2d",rotor_model,ibi);
			if (ti==1) {
				f = fopen(filen, "w");
				PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time (s)\", \"angle (rad)\", \"angvel_axis (s<sup>-1</sup>)\", \"Force_axis (N)\", \"Torque_fluid (Nm)\", \"Torque_generator (Nm)\", \"Uref (m/s)\", \"Ud (m/s)\", \"TSR \"\n");
			} 
			else f = fopen(filen, "a");
			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le \n",Time_real, ang_real, angvel_real, Force_real, Torquea_real, Torquec_real, Uref_real, Ud_real, ibm[ibi].Tipspeedratio);
			fclose(f);
		}
	}
	return(0);
}

// exporting the forces on the turbines, actuator line model
PetscErrorCode Calc_eulforces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi)
{
  	double        pi = 3.141592653589793;
	double A_Sum = 0.0, P_Sum = 0.0, U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0, Fx_Sum=0.0, Fy_Sum=0.0, Fz_Sum=0.0;
  	PetscReal F_z, P, r, rx, ry, rz; 

  	int rank=0;
  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);



  	DA              da = user->da, fda = user->fda;
  	DALocalInfo     info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***lf_eul, ***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj, ***nvert;

  	Vec		Coor;

  	PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
  
  	double		dfunc;

  	double r1, r2, r3;

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




  	DAVecGetArray(fda, user->lF_eul, &lf_eul);
  	DAVecGetArray(fda, user->lCsi,  &csi);
  	DAVecGetArray(fda, user->lEta,  &eta);
  	DAVecGetArray(fda, user->lZet,  &zet);
  	DAVecGetArray(da,  user->lAj,  &aj);

  	DAGetGhostedCoordinates(da, &Coor);
  	DAVecGetArray(fda, Coor, &coor);

	A_Sum = 0.0; P_Sum = 0.0; U_Sum = 0.0; F_Sum=0.0; M_Sum=0.0; Fx_Sum=0.0; Fy_Sum=0.0; Fz_Sum=0.0;
	for (k=lzs;k<lze;k++) 
	for (j=lys;j<lye;j++) 
	for (i=lxs;i<lxe;i++) {

	        double area_x = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	        double area_y = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	        double area_z = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );

		double vol = 1.0/aj[k][j][i];
		double nx=fsi[0].nx_tb, ny=fsi[0].ny_tb, nz=fsi[0].nz_tb;

		double F_axis=lf_eul[k][j][i].x*vol*nx/area_x+lf_eul[k][j][i].y*vol*ny/area_y+lf_eul[k][j][i].z*vol*nz/area_z;

      		F_Sum += F_axis;


		double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	
		double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	
		double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

		double xj = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;	
		double yj = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;	
		double zj = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25; 

		double xk = (coor[k  ][j  ][i].x + coor[k][j-1][i].x + coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;	
		double yk = (coor[k  ][j  ][i].y + coor[k][j-1][i].y + coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;	
		double zk = (coor[k  ][j  ][i].z + coor[k][j-1][i].z + coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25; 


      		double rxi = xi-fsi[0].x_c;	
      		double ryi = yi-fsi[0].y_c;	
      		double rzi = zi-fsi[0].z_c;	


      		double rxj = xj-fsi[0].x_c;	
      		double ryj = yj-fsi[0].y_c;	
      		double rzj = zj-fsi[0].z_c;	

      		double rxk = xk-fsi[0].x_c;	
      		double ryk = yk-fsi[0].y_c;	
      		double rzk = zk-fsi[0].z_c;	

      		double M_x= ryk*lf_eul[k][j][i].z*vol/area_z - rzj*lf_eul[k][j][i].y*vol/area_y ;  
      		double M_y= rzi*lf_eul[k][j][i].x*vol/area_x - rxk*lf_eul[k][j][i].z*vol/area_z ;  
      		double M_z= rxj*lf_eul[k][j][i].y*vol/area_y - ryi*lf_eul[k][j][i].x*vol/area_x ;  

		M_Sum+=M_x*nx+M_y*ny+M_z*nz;

    	}

	double Fsum_global, Msum_global;

	MPI_Allreduce (&F_Sum, &Fsum_global, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	MPI_Allreduce (&M_Sum, &Msum_global, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);



	if (!rank) {
      		FILE *f;
      		char filen[80];
	     	sprintf(filen, "Turbine_eulAL%2.2d_%2.2d",rotor_model,0);
      		if (ti==1) {
			f = fopen(filen, "w");
      			PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time\", \"angle\", \"angvel_axis\", \"Force_axis\", \"Torque_fluid\", \"TSR\" \n");
		} else f = fopen(filen, "a");

      		PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le \n",ti*user->dt, fsi[ibi].ang_axis, fsi[ibi].angvel_axis, -Fsum_global, -Msum_global, ibm[ibi].Tipspeedratio);
	      	fclose(f);	
	}


  	DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
  	DAVecRestoreArray(fda, user->lCsi,  &csi);
  	DAVecRestoreArray(fda, user->lEta,  &eta);
  	DAVecRestoreArray(fda, user->lZet,  &zet);
  	DAVecRestoreArray(da,  user->lAj,  &aj);

  	DAVecRestoreArray(fda, Coor, &coor);

	return(0);
}

PetscErrorCode Uref_ACL(UserCtx *user, IBMNodes *ibm_ACL, IBMNodes *ibm_ACD, FSInfo *fsi_wt, int NumberOfObjects)
{
/*
Calculates the reference velocity (U\_ref) for actuator line model. 
This value corresponds to the space averaged velocity along a disk of same diameter as the 
rotor and located some distance upstream of the turbine. The value is multiplied by the 
disk normal which points downstream.
*/
	PetscInt      l, ibi;
	PetscReal     A_Sum, U_Sum;
	Calc_U_lagr(user, ibm_ACD, fsi_wt, NumberOfObjects);
	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		double nx=fsi_wt[ibi].nx_tb, ny=fsi_wt[ibi].ny_tb, nz=fsi_wt[ibi].nz_tb;
		U_Sum = 0.0; A_Sum = 0.0; 
		for (l=0; l<ibm_ACD[ibi].n_elmt; l++) {
			double U_axis=ibm_ACD[ibi].U_lagr_x[l]*nx+ibm_ACD[ibi].U_lagr_y[l]*ny+ibm_ACD[ibi].U_lagr_z[l]*nz;
			U_Sum += U_axis*ibm_ACD[ibi].dA[l];
			A_Sum += ibm_ACD[ibi].dA[l];
		}
		U_Sum /= A_Sum;
		ibm_ACL[ibi].U_ref = U_Sum; // / (1.0 - indf_a);
	}
	return(0);
}

// Translate the coordinates of the disk 1D upstream. Assume z is the streamwise direction from negtive to positive
PetscErrorCode Trans1DUp_XYZ(IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

	int ibi, i;
	
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		
		double R=fsi[ibi].r_rotor/reflength_wt;

                for (i=0; i<ibm[ibi].n_v; i++){
                        ibm[ibi].z_bp[i]-=2.0*R;
                }

                for (i=0; i<ibm[ibi].n_elmt; i++){
                        ibm[ibi].cent_z[i]-=2.0*R;
                }
	}

    	return(0);
}


// update the coordinates for grid nodes of actuator disks/lines for moving frame
PetscErrorCode UpdateXYZ_MoveFrame(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

	int ibi, i;

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {
		fsi[ibi].x_c=fsi[ibi].x_c0+ti*user->dt*u_frame;
		fsi[ibi].y_c=fsi[ibi].y_c0+ti*user->dt*v_frame;
		fsi[ibi].z_c=fsi[ibi].z_c0+ti*user->dt*w_frame;


		if (ii_periodicWT) {
                	double Lp_x=Nx_WT*Sx_WT;
			if ((fsi[ibi].x_c-xmin)<1.e-9) {
                                int np_x = (int)((xmax-fsi[ibi].x_c+1.e-9)/Lp_x);
                                fsi[ibi].x_c+=Lp_x*np_x;
			}
			if ((fsi[ibi].x_c-xmax)>-1.e-9) {
                                int np_x = (int)((fsi[ibi].x_c-xmin+1.e-9)/Lp_x);
				fsi[ibi].x_c-=Lp_x*np_x;
			}
		}

		if (jj_periodicWT) {
                	double Lp_y=Ny_WT*Sy_WT;
                        if ((fsi[ibi].y_c-ymin)<1.e-9) {
				int np_y = (int)((ymax-fsi[ibi].y_c+1.e-9)/Lp_y);
				fsi[ibi].y_c+=Lp_y*np_y;
			}
                        if ((fsi[ibi].y_c-ymax)>-1.e-9) {
                                int np_y = (int)((fsi[ibi].y_c-ymin+1.e-9)/Lp_y);
				fsi[ibi].y_c-=Lp_y*np_y;
			}
		}

		if (kk_periodicWT) {
                	double Lp_z=Nz_WT*Sz_WT;
                        if ((fsi[ibi].z_c-zmin)<1.e-9) {
				int np_z = (int)((zmax-fsi[ibi].z_c+1.e-9)/Lp_z);
				fsi[ibi].z_c+=Lp_z*np_z;
			}
                        if ((fsi[ibi].z_c-zmax)>-1.e-9) {
                                int np_z = (int)((fsi[ibi].z_c-zmin+1.e-9)/Lp_z);
				fsi[ibi].z_c-=Lp_z*np_z;
			}
		}

//			if (!rank) printf("z_c: %d %le %le \n", ibi, fsi[ibi].z_c0, fsi[ibi].z_c);
                for (i=0; i<ibm[ibi].n_v; i++){
//			printf("z: %d %le %le \n", rank, ibm[ibi].z_bp[i], ibm[ibi].z_bp0[i]+fsi[ibi].z_c0);
//			int a;
//			cout << "Hi \n"; 
//			cin >> a;
                        ibm[ibi].x_bp[i]=ibm[ibi].x_bp_i[i]+fsi[ibi].x_c;
                        ibm[ibi].y_bp[i]=ibm[ibi].y_bp_i[i]+fsi[ibi].y_c;
                        ibm[ibi].z_bp[i]=ibm[ibi].z_bp_i[i]+fsi[ibi].z_c;
                }

                for (i=0; i<ibm[ibi].n_elmt; i++){
                        int n1e = ibm[ibi].nv1[i], n2e = ibm[ibi].nv2[i], n3e = ibm[ibi].nv3[i];
                        if (rotor_model == 3 || rotor_model == 6) {
                                ibm[ibi].cent_x[i]= (ibm[ibi].x_bp[n1e]+ibm[ibi].x_bp[n2e])/2.;
                                ibm[ibi].cent_y[i]= (ibm[ibi].y_bp[n1e]+ibm[ibi].y_bp[n2e])/2.;
                                ibm[ibi].cent_z[i]= (ibm[ibi].z_bp[n1e]+ibm[ibi].z_bp[n2e])/2.;
                        }
                        if (rotor_model == 1 || rotor_model == 2) {
                                ibm[ibi].cent_x[i]= (ibm[ibi].x_bp[n1e]+ibm[ibi].x_bp[n2e]+ibm[ibi].x_bp[n3e])/3.;
                                ibm[ibi].cent_y[i]= (ibm[ibi].y_bp[n1e]+ibm[ibi].y_bp[n2e]+ibm[ibi].y_bp[n3e])/3.;
                                ibm[ibi].cent_z[i]= (ibm[ibi].z_bp[n1e]+ibm[ibi].z_bp[n2e]+ibm[ibi].z_bp[n3e])/3.;
                        }

                }


	}

	if (!rank) {
		FILE *fd;
              	char str[256];
               	sprintf(str, "./CenterWT.dat");
               	fd = fopen(str, "w");

      		for (ibi=0;ibi<NumberOfObjects;ibi++) {
                       	fprintf(fd, "%le %le %le \n", fsi[ibi].x_c, fsi[ibi].y_c, fsi[ibi].z_c);
               	}
        		fclose(fd);
	}




    	return(0);
}

/* ==================================================================================             */
PetscErrorCode Export_lagrdata(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension)
{

  	int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  	PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  	int i,j,k;
  	int n1e, n2e, n3e;
  	PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  	PetscReal rx,ry,rz;
  	int n1;

	if (dimension == 2) {
        	if ( ti==tistart+1 || ti % tiout == 0 /*ti == (ti/tiout) * tiout */) {
		        int rank=0;
        		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
			int i;

		        if (!rank) {
        		        FILE *f;
	                	char filen[80];
        		        sprintf(filen, "%s/%s%06d_%03d_nf.dat",path,fname,ti,ibi);

		                f = fopen(filen, "w");
	
				int n_v=ibm->n_v; 	
				int n_elmt=ibm->n_elmt; 	

		                PetscFPrintf(PETSC_COMM_WORLD, f, "TITLE = \"Actuator disk mesh\" \n ");
			    	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x, y, z, ub_x, ub_y, ub_z, tmprtb, n_x,n_y,n_z, F_lagr_x, F_lagr_y, F_lagr_z, s2l, color, U_lagr_x, U_lagr_y, U_lagr_z, Tmprt_lagr, Ftmprt_lagr\n");
		    		PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-7]=NODAL,[8-20]=CELLCENTERED)\n", n_v, n_elmt);
                		PetscFPrintf(PETSC_COMM_WORLD, f, "STRANDID=0.1 SOLUTIONTIME=%le \n", ((double)ti)*dt);

		    		for (i=0; i<n_v; i++) {
	      				PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->x_bp[i]);
		    		}
    				for (i=0; i<n_v; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->y_bp[i]);
		    		}
		    		for (i=0; i<n_v; i++) {	
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->z_bp[i]);
		    		}

		    		for (i=0; i<n_v; i++) {
	      				PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].x);
		    		}
    				for (i=0; i<n_v; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].y);
		    		}
		    		for (i=0; i<n_v; i++) {	
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].z);
		    		}
		    		for (i=0; i<n_v; i++) {	
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->tmprt[i]);
		    		}



		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_x[i]);
		    		}
    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_y[i]);
    				}
		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_z[i]);
		    		}

		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_x[i]);
		    		}
    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_y[i]);
    				}
		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_z[i]);
		    		}

 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->s2l[i]);
		    		}


 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->color[i]);
		    		}
 


 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_x[i]);
		    		}

 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_y[i]);
		    		}

 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_z[i]);
		    		}



 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Tmprt_lagr[i]);
		    		}

 		    		for (i=0; i<n_elmt; i++) {
      					PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Ftmprt_lagr[i]);
		    		}


    				for (i=0; i<n_elmt; i++) {
		      			PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
		    		}
 

                		fclose(f);

		        }
        	}

	}

	if (dimension == 1) {
		if ( ti==tistart+1 || ti % tiout == 0 /*ti == (ti/tiout) * tiout */) {
	    	int rank=0;
	    	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

		if (!rank) {
        	        FILE *f;
	    		char filen[80];  
        	        sprintf(filen, "%s/%s%06d_%03d_nf.dat", path,fname, ti,ibi);
	                f = fopen(filen, "w");

			PetscFPrintf(PETSC_COMM_WORLD, f, "TITLE = \"Actuator line mesh\" \n ");
			PetscFPrintf(PETSC_COMM_WORLD, f, "VARIABLES = \"X\", \"Y\", \"Z\", \"ub_x\", \"ub_y\", \"ub_z\", \"tmprtb\", \"color\", \"U_lagr_x\", \"U_lagr_y\", \"U_lagr_z\", \"F_lagr_x\", \"F_lagr_y\", \"F_lagr_z\", \"Tmprt_lagr\", \"Ftmprt_lagr\"\n");
			PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"P_1\", DATAPACKING=BLOCK, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG, VARLOCATION=([1-7]=NODAL,[8-16]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
			PetscFPrintf(PETSC_COMM_WORLD, f, "STRANDID=0.1 SOLUTIONTIME=%le \n", ((double)ti)*dt);

	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->x_bp[i]);
	                }
	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->y_bp[i]);
	                }
	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->z_bp[i]);
	                }

	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].x);
	                }
	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].y);
	                }
	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].z);
	                }

	                for (i=0; i<ibm->n_v; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->tmprt[i]);
	                }


	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%d \n", ibm->color[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_x[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_y[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_z[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_x[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_y[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_z[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->Tmprt_lagr[i]);
	                }

	                for (i=0; i<ibm->n_elmt; i++) {
        	                PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->Ftmprt_lagr[i]);
	                }




			for (i=0; i<ibm->n_elmt; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d \n", ibm->nv1[i]+1, ibm->nv2[i]+1);
			}


	                fclose(f);
			
	  	}
	 	} 
	}

  	return(0);
}




/* ==================================================================================             */
PetscErrorCode rotor_Rot(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension)
{
	int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
	PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
	int i,j,k;
	int n1e, n2e, n3e;
	PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
	PetscReal rx,ry,rz;
	int n1;
	double rot_angle;
	Cmpnts p, q, nr, na, nt; 
	if (ti==tistart && rstart_turbinerotation) {
		for (i=0; i<n_v; i++) {
			p.x=ibm->x_bp0[i]-FSinfo->x_c;
			p.y=ibm->y_bp0[i]-FSinfo->y_c;
			p.z=ibm->z_bp0[i]-FSinfo->z_c;
			na.x=FSinfo->nx_tb;	
			na.y=FSinfo->ny_tb;	
			na.z=FSinfo->nz_tb;	
			double theta=FSinfo->ang_axis;
			q=ArbitraryRotate(p,theta,na);
			ibm->x_bp[i]=q.x+FSinfo->x_c;
			ibm->y_bp[i]=q.y+FSinfo->y_c;
			ibm->z_bp[i]=q.z+FSinfo->z_c;
		}
	}
	PetscPrintf(PETSC_COMM_WORLD, "angvel of turbine %le \n", FSinfo->angvel_axis);
	for (i=0; i<n_v; i++) {
		p.x=ibm->x_bp[i]-FSinfo->x_c;
		p.y=ibm->y_bp[i]-FSinfo->y_c;
		p.z=ibm->z_bp[i]-FSinfo->z_c;
		na.x=FSinfo->nx_tb;	//rotor normal direction pointing downstream and imported in tur
		na.y=FSinfo->ny_tb;	
		na.z=FSinfo->nz_tb;	
		double theta=FSinfo->angvel_axis*dt;
		q=ArbitraryRotate(p,theta,na);
		ibm->x_bp[i]=q.x+FSinfo->x_c;
		ibm->y_bp[i]=q.y+FSinfo->y_c;
		ibm->z_bp[i]=q.z+FSinfo->z_c;
		double rx = ibm->x_bp[i]-FSinfo->x_c;
		double ry = ibm->y_bp[i]-FSinfo->y_c;
		double rz = ibm->z_bp[i]-FSinfo->z_c;
		double rr = sqrt(rx*rx+ry*ry+rz*rz)+1.e-19;
		nr.x = rx/rr; 
		nr.y = ry/rr; 
		nr.z = rz/rr;
		nt.x=na.y*nr.z-na.z*nr.y;
		nt.y=na.z*nr.x-na.x*nr.z;
		nt.z=na.x*nr.y-na.y*nr.x;
		double Ut=FSinfo->angvel_axis*rr;
		ibm->u[i].x = Ut*nt.x;
		ibm->u[i].y = Ut*nt.y;
		ibm->u[i].z = Ut*nt.z;
	}
	if (dimension == 2) {
		//This corresponds to the actuator disk
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
			if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
			(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
				ibm->ns_x[i] = 1.;
				ibm->ns_y[i] = 0.;
				ibm->ns_z[i] = 0 ;
				// nt = ns x nf
				ibm->nt_x[i] = 0.;
				ibm->nt_y[i] = 1.;
				ibm->nt_z[i] = 0.;
			} 
			else {
				ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->ns_z[i] = 0 ;
				// nt = ns x nf
				ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
			}
			ibm->dA[i] = dr/2.;
			ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e]+ibm->x_bp[n3e])/3.;
			ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e]+ibm->y_bp[n3e])/3.;
			ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e]+ibm->z_bp[n3e])/3.;
		}
	}
	int nb;
	int n_elmt_1 = (ibm->n_elmt)/num_blade;
	if (dimension == 1) {
		//This corresponds to the actuator line
		//for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line. and dA denotes the length of each element.
		for (nb=0; nb<num_blade; nb++) {
			for (j=0; j<n_elmt_1; j++) {
				i = nb * n_elmt_1 + j;
				n1e = ibm->nv1[i]; n2e = ibm->nv2[i]; 
				dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
				dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
				dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
				ibm->nf_x[i] = dx12;
				ibm->nf_y[i] = dy12;
				ibm->nf_z[i] = dz12;
				dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + ibm->nf_z[i]*ibm->nf_z[i]);
				ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
				ibm->ns_x[i] = 0.;     
				ibm->ns_y[i] = 0.;     
				ibm->ns_z[i] = 0. ;
				ibm->nt_x[i] = 0.;
				ibm->nt_y[i] = 0.;
				ibm->nt_z[i] = 0.;
				ibm->dA[i] = dr;
				ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e])/2.;
				ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e])/2.;
				ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e])/2.;
			}
		}
	}
	return(0);
}

PetscErrorCode rotor_Rot_6dof_fsi(FSInfo *FSinfoPlat, FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension)
{
/*
This functions is an extension of rotor_Rot for the case that the turbine rotor 
is moving (such as a floating turbine) in adition to the rotor rotation
*/
	int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
	PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
	int i,j,k;
	int n1e, n2e, n3e;
	PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
	PetscReal rx,ry,rz;
	int n1;
	double rot_angle;
	Cmpnts p, q, nr, na, nt; 
	for (i=0; i<n_v; i++) {
		//get turbine coordinates at the initial time and centered at the origin
		p.x=ibm->x_bp0[i]-FSinfo->x_c0;
		p.y=ibm->y_bp0[i]-FSinfo->y_c0;
		p.z=ibm->z_bp0[i]-FSinfo->z_c0;
		//rotor normal direction pointing downstream (value specified at turbine control file turbine.imp)
		na.x=FSinfo->nx_tb;	
		na.y=FSinfo->ny_tb;	
		na.z=FSinfo->nz_tb;
		double theta=FSinfo->ang_axis;//Rotor rotation accumulated since beginning of the simulation
		q=ArbitraryRotate(p,theta,na);//apply rotor rotation
		//Assign translation and rotation
		q.x=q.x+FSinfo->x_c0;
		q.y=q.y+FSinfo->y_c0;
		q.z=q.z+FSinfo->z_c0;
		if(turbine_prescrived_motion){
			double oscilating_freq=18.849555921538759430775860299677;//3Hz
			double pitch_amplitude=0.06981317007977318307694763073954;//4deg
			double heave_amplitude=0.0035;
			if(turbine_prescrived_motion_pitch){
				FSinfo[ibi].S_ang_r[0]=pitch_amplitude*sin((double)dt*(double)ti*oscilating_freq);
				FSinfo[ibi].S_ang_r[2]=0.;
				FSinfo[ibi].S_ang_r[4]=0.;
				FSinfo->S_new[0]=0.;
				FSinfo->S_new[2]=0.;
				FSinfo->S_new[4]=0.;				
			}
			if(turbine_prescrived_motion_heave){
				FSinfo->S_new[0]=0.;
				FSinfo->S_new[2]=heave_amplitude*sin((double)dt*(double)ti*oscilating_freq);
				FSinfo->S_new[4]=0.;
				FSinfo[ibi].S_ang_r[0]=0.;
				FSinfo[ibi].S_ang_r[2]=0.;
				FSinfo[ibi].S_ang_r[4]=0.;				
			}
		}
		else if(turbine_6dof_fsi_motion){
			FSinfo->S_new[0]=FSinfoPlat->S_new[0];
			FSinfo->S_new[2]=FSinfoPlat->S_new[2];
			FSinfo->S_new[4]=FSinfoPlat->S_new[4];
			FSinfo[ibi].S_ang_r[0]=FSinfoPlat[ibi].S_ang_r[0];
			FSinfo[ibi].S_ang_r[2]=FSinfoPlat[ibi].S_ang_r[2];
			FSinfo[ibi].S_ang_r[4]=FSinfoPlat[ibi].S_ang_r[4];			
		}
		else{
			FSinfo->S_new[0]=0.;
			FSinfo->S_new[2]=0;
			FSinfo->S_new[4]=0.;
			FSinfo[ibi].S_ang_r[0]=0.;
			FSinfo[ibi].S_ang_r[2]=0.;
			FSinfo[ibi].S_ang_r[4]=0.;	
		}

		//add to the rotor rotation the turbine rigid body motion
		//rotations
		rotate_xyz6dof (ti, dt, FSinfo[ibi].S_ang_r[0],FSinfo[ibi].S_ang_r[2],FSinfo[ibi].S_ang_r[4], x_r, y_r, z_r, q.x, q.y, q.z, &ibm->x_bp[i], &ibm->y_bp[i], &ibm->z_bp[i], &rot_angle);
		//x_r, y_r, z_r is the center of rotation of the rigid body motion specified at control.dat
		//translate the body
    ibm->x_bp[i] = ibm->x_bp[i]+(FSinfo->S_new[0]);
    ibm->y_bp[i] = ibm->y_bp[i]+(FSinfo->S_new[2]);
    ibm->z_bp[i] = ibm->z_bp[i]+(FSinfo->S_new[4]);		
		//apply the velocity at the body
		ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / dt ;
		ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / dt ;
		ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / dt ;
		//store rotor nodes to old values		
		ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];		
	}
	//displace the center of rotation of the rotor
	rotate_xyz6dof (ti, dt, FSinfo[ibi].S_ang_r[0],FSinfo[ibi].S_ang_r[2],FSinfo[ibi].S_ang_r[4], x_r, y_r, z_r, FSinfo->x_c0, FSinfo->y_c0, FSinfo->z_c0, &FSinfo->x_c, &FSinfo->y_c, &FSinfo->z_c, &rot_angle);
	FSinfo->x_c+=FSinfo->S_new[0];
	FSinfo->y_c+=FSinfo->S_new[2];
	FSinfo->z_c+=FSinfo->S_new[4];	
	if (dimension == 2) {
		//This corresponds to the actuator disk
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
			if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
			(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
				ibm->ns_x[i] = 1.;
				ibm->ns_y[i] = 0.;
				ibm->ns_z[i] = 0 ;
				// nt = ns x nf
				ibm->nt_x[i] = 0.;
				ibm->nt_y[i] = 1.;
				ibm->nt_z[i] = 0.;
			} 
			else {
				ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->ns_z[i] = 0 ;
				// nt = ns x nf
				ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
				ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
			}
			ibm->dA[i] = dr/2.;
			ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e]+ibm->x_bp[n3e])/3.;
			ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e]+ibm->y_bp[n3e])/3.;
			ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e]+ibm->z_bp[n3e])/3.;
		}
	}
	int nb;
	int n_elmt_1 = (ibm->n_elmt)/num_blade;
	if (dimension == 1) {
		//This corresponds to the actuator line
		//for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line. and dA denotes the length of each element.
		for (nb=0; nb<num_blade; nb++) {
			for (j=0; j<n_elmt_1; j++) {
				i = nb * n_elmt_1 + j;
				n1e = ibm->nv1[i]; n2e = ibm->nv2[i]; 
				dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
				dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
				dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
				ibm->nf_x[i] = dx12;
				ibm->nf_y[i] = dy12;
				ibm->nf_z[i] = dz12;
				dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + ibm->nf_z[i]*ibm->nf_z[i]);
				ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
				ibm->ns_x[i] = 0.;     
				ibm->ns_y[i] = 0.;     
				ibm->ns_z[i] = 0. ;
				ibm->nt_x[i] = 0.;
				ibm->nt_y[i] = 0.;
				ibm->nt_z[i] = 0.;
				ibm->dA[i] = dr;
				ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e])/2.;
				ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e])/2.;
				ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e])/2.;
			}
		}
	}
	return(0);
}


PetscErrorCode refAL_Rot(FSInfo *FSinfo, IBMNodes *ibm, IBMNodes *ibm_ref, int ibi)
{
	int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
	PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
	int i,j,k;
	int n1e, n2e, n3e;
	PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
	PetscReal rx,ry,rz;
	int n1;
	double rot_angle;
	Cmpnts p, q, nr, na, nt; 
	for (i=0; i<n_v; i++) {
		p.x=ibm->x_bp[i]-FSinfo->x_c;
		p.y=ibm->y_bp[i]-FSinfo->y_c;
		p.z=ibm->z_bp[i]-FSinfo->z_c;
		na.x=FSinfo->nx_tb;	
		na.y=FSinfo->ny_tb;	
		na.z=FSinfo->nz_tb;	
		double isign = (1.e-19+FSinfo->angvel_axis)/(fabs(FSinfo->angvel_axis)+1.e-19);
		if (isign>0.0) isign = 1.0;
		else isign = -1.0;
		double theta = isign * refangle_AL*M_PI/180.0;
		q=ArbitraryRotate(p,theta,na);
		ibm_ref->x_bp[i]=q.x+FSinfo->x_c;
		ibm_ref->y_bp[i]=q.y+FSinfo->y_c;
		ibm_ref->z_bp[i]=q.z+FSinfo->z_c;
	}
	int nb;
	int n_elmt_1 = (ibm->n_elmt)/num_blade;
	for (nb=0; nb<num_blade; nb++) {
		for (j=0; j<n_elmt_1; j++) {
			i = nb * n_elmt_1 + j;
			n1e = ibm_ref->nv1[i]; n2e = ibm_ref->nv2[i]; 
			dx12 = ibm_ref->x_bp[n2e] - ibm_ref->x_bp[n1e];
			dy12 = ibm_ref->y_bp[n2e] - ibm_ref->y_bp[n1e];
			dz12 = ibm_ref->z_bp[n2e] - ibm_ref->z_bp[n1e];
			ibm_ref->nf_x[i] = dx12;
			ibm_ref->nf_y[i] = dy12;
			ibm_ref->nf_z[i] = dz12;
			dr = sqrt(ibm_ref->nf_x[i]*ibm_ref->nf_x[i] + ibm_ref->nf_y[i]*ibm_ref->nf_y[i] + ibm_ref->nf_z[i]*ibm_ref->nf_z[i]);
			ibm_ref->nf_x[i] /=dr; ibm_ref->nf_y[i]/=dr; ibm_ref->nf_z[i]/=dr;
			ibm_ref->ns_x[i] = 0.;     
			ibm_ref->ns_y[i] = 0.;     
			ibm_ref->ns_z[i] = 0. ;
			ibm_ref->nt_x[i] = 0.;
			ibm_ref->nt_y[i] = 0.;
			ibm_ref->nt_z[i] = 0.;
			ibm_ref->dA[i] = dr;
			ibm_ref->cent_x[i]= (ibm_ref->x_bp[n1e]+ibm_ref->x_bp[n2e])/2.;
			ibm_ref->cent_y[i]= (ibm_ref->y_bp[n1e]+ibm_ref->y_bp[n2e])/2.;
			ibm_ref->cent_z[i]= (ibm_ref->z_bp[n1e]+ibm_ref->z_bp[n2e])/2.;
		}
	}
	return(0);
}





PetscErrorCode calc_s2l(IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

	int ibi, elmt_s, elmt_l;
	double rs, rl, dd_min, dd;

	double xs, ys, zs, xl, yl, zl;	
	int colors, colorl;
	double nx, ny, nz;

  	for (ibi=0; ibi<NumberOfObjects; ibi++) {

		for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) {

			dd_min = 100000;
			xs = ibm_surface[ibi].cent_x[elmt_s] - fsi_surface[ibi].x_c;
			ys = ibm_surface[ibi].cent_y[elmt_s] - fsi_surface[ibi].y_c;
			zs = ibm_surface[ibi].cent_z[elmt_s] - fsi_surface[ibi].z_c;

			colors = ibm_surface[ibi].color[elmt_s];
			for (elmt_l=0; elmt_l<ibm_line[ibi].n_elmt; elmt_l++) {
				xl = ibm_line[ibi].cent_x[elmt_l] - fsi_line[ibi].x_c;
				yl = ibm_line[ibi].cent_y[elmt_l] - fsi_line[ibi].y_c;
				zl = ibm_line[ibi].cent_z[elmt_l] - fsi_line[ibi].z_c;

				rl = sqrt(xl*xl+yl*yl+zl*zl);

				nx = xl/rl; ny = yl/rl; nz = zl/rl;

				rs = fabs(xs*nx+ys*ny+zs*nz);

				colorl = ibm_line[ibi].color[elmt_l];
		
      				//PetscPrintf(PETSC_COMM_WORLD, "rs=%le, rl=%le, colors=%d, colorl=%d for %d th line segment\n", rs, rl, colors, colorl, elmt_l);
				double icolor = double(colors-colorl);
				dd=fabs(rs-rl);
				if (dd-dd_min<1.e-9 && colors==colorl) {
					dd_min=dd;
					ibm_surface[ibi].s2l[elmt_s] = elmt_l;
				}
			}
		}
  	}

    	return(0);
}


PetscErrorCode ForceProjection_l2s(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

	int ibi, elmt_s, elmt_l, l;
	
  	for (ibi=0; ibi<NumberOfObjects; ibi++) {

		for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) {
			elmt_l=ibm_surface[ibi].s2l[elmt_s];
                        ibm_surface[ibi].F_lagr_x[elmt_s] = ibm_line[ibi].F_lagr_x[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
                        ibm_surface[ibi].F_lagr_y[elmt_s] = ibm_line[ibi].F_lagr_y[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
                        ibm_surface[ibi].F_lagr_z[elmt_s] = ibm_line[ibi].F_lagr_z[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
		}
  	}

	/*
	if (radialforce_AL) {
		int nv1, nv2, nv3;
		double rx, ry, rz, rr, r1, r2, r3;
		Cmpnts n_blade;

		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			for (l=0; l<ibm_surface[ibi].n_elmt; l++) {
				nv1 = ibm_surface[ibi].nv1[l];
				nv2 = ibm_surface[ibi].nv2[l];
				nv3 = ibm_surface[ibi].nv3[l];

		               	rx = ibm_surface[ibi].x_bp[nv2] - ibm_surface[ibi].x_bp[nv1];
		             	ry = ibm_surface[ibi].y_bp[nv2] - ibm_surface[ibi].y_bp[nv1];
		             	rz = ibm_surface[ibi].z_bp[nv2] - ibm_surface[ibi].z_bp[nv1];

		              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
		                rx = ibm_surface[ibi].x_bp[nv3] - ibm_surface[ibi].x_bp[nv2];
		                ry = ibm_surface[ibi].y_bp[nv3] - ibm_surface[ibi].y_bp[nv2];
		                rz = ibm_surface[ibi].z_bp[nv3] - ibm_surface[ibi].z_bp[nv2];

		                r2 = sqrt(rx*rx + ry*ry + rz*rz );

		                rx = ibm_surface[ibi].x_bp[nv1] - ibm_surface[ibi].x_bp[nv3];
		                ry = ibm_surface[ibi].y_bp[nv1] - ibm_surface[ibi].y_bp[nv3];
		                rz = ibm_surface[ibi].z_bp[nv1] - ibm_surface[ibi].z_bp[nv3];

		                r3 = sqrt(rx*rx + ry*ry + rz*rz );

				rr = (r1 + r2 + r3)/3.0;



				rx = ibm_surface[ibi].cent_x[l]-fsi_surface[ibi].x_c;
				ry = ibm_surface[ibi].cent_y[l]-fsi_surface[ibi].y_c;
				rz = ibm_surface[ibi].cent_z[l]-fsi_surface[ibi].z_c;
	
				double r = sqrt(rx*rx+ry*ry+rz*rz);
        	                n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;


				double Ur = ibm_surface[ibi].U_lagr_x[l]*n_blade.x+ibm_surface[ibi].U_lagr_y[l]*n_blade.y+ibm_surface[ibi].U_lagr_z[l]*n_blade.z;

				double Fr =  10 * rr * (0 - Ur)/user->dt;

                        	ibm_surface[ibi].F_lagr_x[l]+=Fr*n_blade.x;
	                        ibm_surface[ibi].F_lagr_y[l]+=Fr*n_blade.y;
        	                ibm_surface[ibi].F_lagr_z[l]+=Fr*n_blade.z;
				

			}
		}
	}

	*/

	/*
	if (AL_Noslip) {
		int nv1, nv2, nv3;
		double rx, ry, rz, rr, r1, r2, r3;
		Cmpnts n_blade;

		for (ibi=0; ibi<NumberOfObjects; ibi++) {
			for (l=0; l<ibm_surface[ibi].n_elmt; l++) {
				nv1 = ibm_surface[ibi].nv1[l];
				nv2 = ibm_surface[ibi].nv2[l];
				nv3 = ibm_surface[ibi].nv3[l];

		               	rx = ibm_surface[ibi].x_bp[nv2] - ibm_surface[ibi].x_bp[nv1];
		             	ry = ibm_surface[ibi].y_bp[nv2] - ibm_surface[ibi].y_bp[nv1];
		             	rz = ibm_surface[ibi].z_bp[nv2] - ibm_surface[ibi].z_bp[nv1];

		              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
		                rx = ibm_surface[ibi].x_bp[nv3] - ibm_surface[ibi].x_bp[nv2];
		                ry = ibm_surface[ibi].y_bp[nv3] - ibm_surface[ibi].y_bp[nv2];
		                rz = ibm_surface[ibi].z_bp[nv3] - ibm_surface[ibi].z_bp[nv2];

		                r2 = sqrt(rx*rx + ry*ry + rz*rz );

		                rx = ibm_surface[ibi].x_bp[nv1] - ibm_surface[ibi].x_bp[nv3];
		                ry = ibm_surface[ibi].y_bp[nv1] - ibm_surface[ibi].y_bp[nv3];
		                rz = ibm_surface[ibi].z_bp[nv1] - ibm_surface[ibi].z_bp[nv3];

		                r3 = sqrt(rx*rx + ry*ry + rz*rz );

				rr = (r1 + r2 + r3)/3.0;



				rx = ibm_surface[ibi].cent_x[l]-fsi_surface[ibi].x_c;
				ry = ibm_surface[ibi].cent_y[l]-fsi_surface[ibi].y_c;
				rz = ibm_surface[ibi].cent_z[l]-fsi_surface[ibi].z_c;
	
				double r = sqrt(rx*rx+ry*ry+rz*rz);
        	                n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;


				double Ur = ibm_surface[ibi].U_lagr_x[l]*n_blade.x+ibm_surface[ibi].U_lagr_y[l]*n_blade.y+ibm_surface[ibi].U_lagr_z[l]*n_blade.z;

				double Fr =  rr * (0 - Ur)/user->dt;

                        	ibm_surface[ibi].F_lagr_x[l]=rr * (0 - ibm_surface[ibi].U_lagr_x[l])/user->dt;
	                        ibm_surface[ibi].F_lagr_y[l]=rr * (0 - ibm_surface[ibi].U_lagr_y[l])/user->dt;
        	                ibm_surface[ibi].F_lagr_z[l]=rr * (0 - ibm_surface[ibi].U_lagr_z[l])/user->dt;
				


			}
		}
	}

	*/


    	return(0);
}



double dfunc_4htail(double r)
{
  if (fabs(r) <= 1.0) {
    return 1-r*r;
  } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) {
    return 2-3*fabs(r)+r*r;
  } else {
    return 0.0;
  }
}

double dfunc_4h(double r)
{
  if (fabs(r) <= 1.0) {
    return (3.0-2.0*fabs(r)+sqrt(1.0+4.0*fabs(r)-4.0*r*r))/8.0;
  } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) {
    return (5.0-2.0*fabs(r)-sqrt(-7.0+12.0*fabs(r)-4.0*r*r))/8.0;
  } else {
    return 0.0;
  }
}

double dfunc_s3h(double r)
{
  if (fabs(r) <= 1.0) {
    return 17.0/48.0 + sqrt(3.0)*3.14159265/108.0 + fabs(r)/4.0 - r*r/4.0 + (1.0-2.0*fabs(r))*sqrt(-12.0*r*r+12.0*fabs(r)+1.0)/16.0 - sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-1.0)/2.0)/12.0;
  } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) {
    return 55.0/48.0 - sqrt(3.0)*3.14159265/108.0 - 13.0*fabs(r)/12.0 + r*r/4.0 + (2.0*fabs(r)-3.0)*sqrt(-12.0*r*r+36.0*fabs(r)-23.0)/48.0 + sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-3.0)/2.0)/36.0;
  } else { 
    return 0.0;
  }

}


double dfunc_2h(double r)
{
  if (fabs(r) < 1.0) {
    return 1.0-fabs(r);
  } else {
    return 0.0;
  }
}


double dfunc_s4h(double r)
{

  if (fabs(r) <= 0.5) {
    return 3.0/8.0+3.14159265/32.0-pow(r,2)/4.0;
  } else if (fabs(r) >= 0.5 && fabs(r) <= 1.5) {
    return 1.0/4.0+(1.0-fabs(r))*sqrt(-2.0+8.0*fabs(r)-4.0*pow(r,2))/8.0-asin(sqrt(2.0)*(fabs(r)-1.0))/8.0;
  } else if (fabs(r) >= 1.5 && fabs(r) <= 2.5) {
    return
17.0/16.0-3.14159265/64.0-3.0*fabs(r)/4.0+pow(r,2)/8.0+(fabs(r)-2.0)*sqrt(-14.0+16.0*fabs(r)-4.0*pow(r,2))/16.0+asin(sqrt(2.0)*(fabs(r)-2.0))/16.0;
  } else {
    return 0.0;
  }

}



double dfunc_sc4h(double r)
{

  if (fabs(r) <= 1.5) {
    return 0.25*(M_PI+2*sin(0.25*M_PI*(2*fabs(r)+1))-2*sin(0.25*M_PI*(2*fabs(r)-1)))/M_PI;
  } else if (fabs(r) >= 1.5 && fabs(r) <= 2.5) {
    return -0.125*(-5*M_PI+2*M_PI*fabs(r)+4*sin(0.25*M_PI*(2*fabs(r)-1)))/M_PI;
  } else {
    return 0.0;
  }

}




Cmpnts ArbitraryRotate(Cmpnts p, double theta, Cmpnts r) 
{
	Cmpnts q = {0.0,0.0,0.0};
	double costheta,sintheta;
	double rr=sqrt(r.x*r.x+r.y*r.y+r.z*r.z)+1.0e-11;
	r.x=r.x/rr; r.y=r.y/rr; r.z=r.z/rr;
	costheta = cos(theta);
	sintheta = sin(theta);
	q.x += (costheta + (1 - costheta) * r.x * r.x) * p.x;
	q.x += ((1 - costheta) * r.x * r.y - r.z * sintheta) * p.y;
	q.x += ((1 - costheta) * r.x * r.z + r.y * sintheta) * p.z;
	q.y += ((1 - costheta) * r.x * r.y + r.z * sintheta) * p.x;
	q.y += (costheta + (1 - costheta) * r.y * r.y) * p.y;
	q.y += ((1 - costheta) * r.y * r.z - r.x * sintheta) * p.z;
	q.z += ((1 - costheta) * r.x * r.z - r.y * sintheta) * p.x;
	q.z += ((1 - costheta) * r.y * r.z + r.x * sintheta) * p.y;
	q.z += (costheta + (1 - costheta) * r.z * r.z) * p.z;
	return(q);
}


double dfunc_4uniform(double r)
{

  if (fabs(r) < 2.1) {
    return 0.25;
  } else {
    return 0.0;
  }
}


double dfunc_6uniform(double r)
{

  if (fabs(r) < 3.00001) {
    return 0.166666667;
  } else {
    return 0.0;
  }
}





double dfunc_nh(double r, double n)
{
  if (fabs(r) < n) {
    return (n-fabs(r))/pow(n,2);
  } else {
    return 0.0;
  }
}


double dfunc_exp(double r, double n)
{
    return exp(-(r/n)*(r/n))/(pow(n,1)*pow(3.1415926,0.5));
    //return exp(-r*r)/pow(3.1415926,1);
}



double dfunc_2hs1(double r)
{
  	if (fabs(r) <=0.5){
    		return 3.0/4.0-pow(r,2.0);
  	} else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
		return pow((2*fabs(r) - 3.0),2.0)/8.0;
  	} else {
    		return 0.0;
  	}
}


double dfunc_2hs2(double r)
{
  	if (fabs(r) <=1.0){
    		return pow(fabs(r),3.0)/2.0 - pow(r,2) + 2.0/3.0;
  	} else if (fabs(r) >= 1.0 && fabs(r) <=2.0) {
		return -pow((fabs(r) - 2.0),3.0)/6.0;
  	} else {
    		return 0.0;
  	}
}



double dfunc_2hs3(double r)
{
  	if (fabs(r) <=0.5){
    		return pow(r,4.0)/4. - (5.*pow(r,2.))/8. + 115./192.;
  	} else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
		return - pow(r,4.)/6. + (5.*pow(fabs(r),3.))/6. - (5.*pow(r,2.))/4. + (5.*fabs(r))/24. + 55./96.;
  	} else if (fabs(r) >= 1.5 && fabs(r) <=2.5) {
		return pow((2.*fabs(r) - 5.),4.)/384.;
  	} else {
    		return 0.0;
  	}
}



double dfunc_2hs4(double r)
{
  	if (fabs(r) <=1.){
    		return - pow(fabs(r),5.)/12. + pow(r,4.)/4. - pow(r,2.)/2. + 11./20.;
  	} else if (fabs(r) >= 1. && fabs(r) <=2.) {
		return pow(fabs(r),5.)/24. - (3.*pow(r,4.))/8. + (5.*pow(fabs(r),3.))/4. - (7.*pow(r,2.))/4. + (5.*fabs(r))/8. + 17./40.;
  	} else if (fabs(r) >= 2. && fabs(r) <=3.) {
		return -pow((fabs(r) - 3.),5.)/120.;
  	} else {
    		return 0.0;
  	}
}



double dfunc_2hs5(double r)
{
  	if (fabs(r) <=0.5){
    		return - pow(r,6.)/36. + (7.*pow(r,4.))/48. - (77.*pow(r,2.))/192. + 5887./11520.;
  	} else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
		return pow(r,6.)/48. - (7.*pow(fabs(r),5.))/48. + (21.*pow(r,4.))/64. - (35.*pow(fabs(r),3.))/288. - (91.*pow(r,2.))/256. - (7.*fabs(r))/768. + 7861./15360.;
  	} else if (fabs(r) >= 1.5 && fabs(r) <=2.5) {
		return - pow(r,6.)/120. + (7.*pow(fabs(r),5.))/60. - (21.*pow(r,4.))/32. + (133.*pow(fabs(r),3.))/72. - (329.*pow(r,2.))/128. + (1267.*fabs(r))/960. + 1379./7680.;
  	} else if (fabs(r) >= 2.5 && fabs(r) <=3.5) {
		return pow((2.*fabs(r) - 7.),6.)/46080.;
  	} else {
    		return 0.0;
  	}
}



double dfunc_2hs6(double r)
{
  	if (fabs(r) <=1.){
    		return pow(fabs(r),7.)/144. - pow(r,6.)/36. + pow(r,4.)/9. - pow(r,2.)/3. + 151./315.;
  	} else if (fabs(r) >= 1. && fabs(r) <=2.) {
		return - pow(fabs(r),7.)/240. + pow(r,6.)/20. - (7.*pow(fabs(r),5.))/30. + pow(r,4.)/2. - (7.*pow(fabs(r),3.))/18. - pow(r,2.)/10. - (7.*fabs(r))/90. + 103./210.;
  	} else if (fabs(r) >= 2. && fabs(r) <=3.) {
		return pow(fabs(r),7.)/720. - pow(r,6.)/36. + (7.*pow(fabs(r),5.))/30. - (19.*pow(r,4.))/18. + (49.*pow(fabs(r),3.))/18. - (23.*pow(r,2.))/6. + (217.*fabs(r))/90. - 139./630.;
  	} else if (fabs(r) >= 3. && fabs(r) <=4.) {
		return -pow((fabs(r) - 4.),7.)/5040.;
  	} else {
    		return 0.0;
  	}
}


double dfunc_4hs1(double r)
{
  	if (fabs(r) <=0.5){
    		return 15./64. - pow(r,2)/16.;
  	} else if (fabs(r) >= 0.5 && fabs(r) <=3.5) {
		return 1.0/4.0 - fabs(r)/16.0;
  	} else if (fabs(r) >= 3.5 && fabs(r) <=4.5) {
		return pow((2.*fabs(r) - 9.),2.)/128.;
  	} else {
    		return 0.0;
  	}
}


double dfunc_nhs2(double r, double n)
{
  	if (fabs(r) <= n-1.){
    		return 1./(2.*n);
  	} else if (fabs(r) >= n-1. && fabs(r) <= n) {
		return fabs(r) - n/2. - (2.*pow(r,2.) + 3.*fabs(r) - 1.)/(4.*n) + 3./4.;
  	} else if (fabs(r) >= n && fabs(r) <= n+1.) {
		return ((n - fabs(r) + 1.)*(2.*n - 2.*fabs(r) + 1.))/(4.*n);
  	} else {
    		return 0.0;
  	}

}



double dfunc_nhs1(double r, double n)
{
  	if (fabs(r) <= n-0.5){
    		return 1./(2.*n);
  	} else if (fabs(r) >= n-0.5 && fabs(r) <= n+0.5) {
		return 1./2. - (fabs(r) - 1./2.)/(2.*n);
  	} else {
    		return 0.0;
  	}

}


double dfunc_h(double r)
{
  if (fabs(r) < 0.5) {
    return 1.0;
  } else {
    return 0.0;
  }
}


double dfunc_4h2peak(double r)
{
	if (fabs(r)>2.0) {
		return 0.0;
	} else if (r>=-2.0 && r<=0.0) {
		return 0.5*(1.0-fabs(r+1.0));
	} else if (r>=0.0 && r<=2.0) {
		return 0.5*(1.0-fabs(r-1.0));
	}
}




PetscErrorCode deltafunc_test( ) 
{

	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (!rank) {
    		FILE *f;
	    	char filen[80];  
        	sprintf(filen, "Deltafunction.dat");
	        f = fopen(filen, "w");

		int i;
		double dx, r;

		dx=100.0/1000.0;
		double n = halfwidth_dfunc;
	        for (i=0; i<1001; i++) {
			double r=-50.0+dx*(double)i;
        	     	PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le \n", r, dfunc_s3h(r), dfunc_exp(r,n),  dfunc_sc4h(r),  dfunc_2hs3(r),  dfunc_2hs4(r),  dfunc_nhs1(r,n),  dfunc_nhs2(r,n));
	       	}

	       	fclose(f);



		double m0, m1, m0_s1, m1_s1, m0_s2, m1_s2;
		m0_s1=0.0;
		m1_s1=0.0;
		m0_s2=0.0;
		m1_s2=0.0;
		m0=0.0;
		m1=0.0;


		int NN = 1000;

		int N=NN/2;

		srand( time(NULL)) ;



		int rn = rand() % 200 - 100;

		double rrn = (double)rn/100.;
	        for (i=0; i<NN; i++) {

			/*
			double r = i-N+rrn;
			if (r<3.0) {
				m0_s1 += dfunc_sc4h(r);
				m1_s1 += r*dfunc_sc4h(r);
			}
			*/

			double nwitdh=2.0;
			double r = (i-N+rrn)/nwitdh;
			if (r<3.0) {
				m0_s1 += dfunc_s4h(r)/nwitdh;
				m1_s1 += r*dfunc_s4h(r)/nwitdh;
			}


			r = i-N+rrn;
			m0_s2 += dfunc_exp(r,n);
			m1_s2 += r*dfunc_exp(r,n);

			if (r<n) {
				m0 += dfunc_exp(r,n);
				m1 += r*dfunc_exp(r,n);
			}

		}
		PetscPrintf(PETSC_COMM_WORLD, "dfunc_s4h (2.5 width distributed): halfwidth = %le, 0 moment = %le, random shift = %le \n", 2.0, m0_s1, rrn);
		PetscPrintf(PETSC_COMM_WORLD, "dfunc_s4h (2.5 width distributed): halfwidth = %le, 1 moment = %le, random shift = %le \n", 2.0, m1_s1, rrn);

		PetscPrintf(PETSC_COMM_WORLD, "dfunc_exp (%le width distributed): halfwidth = %le 0 moment = %le, random shift = %le \n", n, n, m0, rrn);
		PetscPrintf(PETSC_COMM_WORLD, "dfunc_exp (%le width distributed): halfwidth = %le 1 moment = %le, random shift = %le \n", n, n, m1, rrn);

		PetscPrintf(PETSC_COMM_WORLD, "dfunc_exp (500 width distributed): halfwidth = %le, 0 moment = %le, random shift = %le \n", n, m0_s2, rrn);
		PetscPrintf(PETSC_COMM_WORLD, "dfunc_exp (500 width distributed): halfwidth = %le, 1 moment = %le, random shift = %le \n", n, m1_s2, rrn);

	}

	return(0);
}

PetscErrorCode TurbineTorqueControl_Input(FSInfo *fsi, IBMNodes *ibm)
{
	int  i, ibi;
	for (ibi=0; ibi<NumberOfTurbines; ibi++) {
		FILE *f;
		char filen[80];  
		sprintf(filen, "TurbineTorqueControl%6.6d_%3.3d.dat", (int)ti, ibi);
		if (f) printf("Read %s for turbine restart file!! \n", filen);
		f = fopen(filen, "r");
		if (!f) printf("Cannot open %s !!", filen);
		double angvel, TSR;
		fscanf(f, "%le %le %le %le %le \n", &(fsi[ibi].Torque_aero), &(fsi[ibi].Torque_generator), &(fsi[ibi].ang_axis), &(angvel), &(TSR)); 
		if (!FixTipSpeedRatio) ibm[ibi].Tipspeedratio=TSR;
		if (!fixturbineangvel) fsi[ibi].angvel_axis=0.0;
		fclose(f);
	}
  return(0);
}

PetscErrorCode TurbineTorqueControl_Output(FSInfo *fsi, IBMNodes *ibm)
{
	int rank, i, ibi;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	PetscBarrier(PETSC_NULL);
	if (!rank) {
		for (ibi=0; ibi<NumberOfTurbines; ibi++) {
			FILE *f;
			char filen[80];
			sprintf(filen, "TurbineTorqueControl%6.6d_%3.3d.dat", (int)ti+1, ibi);
			f = fopen(filen, "w");
			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le \n", fsi[ibi].Torque_aero, fsi[ibi].Torque_generator, fsi[ibi].ang_axis, fsi[ibi].angvel_axis, ibm[ibi].Tipspeedratio); 
			fclose(f);
		}
	}	
	return(0);
}

// Calculating rotational velocity of turbine 
PetscErrorCode Calc_turbineangvel(PetscReal dt, IBMNodes *ibm, FSInfo *fsi)
{
	PetscInt      l, ibi;
	double        pi = 3.141592653589793;
	PetscReal A_Sum = 0.0, P_Sum = 0.0, U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0, Fx_Sum=0.0, Fy_Sum=0.0, Fz_Sum=0.0;
	PetscReal F_z, P, r, rx, ry, rz; 
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	double R0=fsi[0].r_rotor/reflength_wt;
	double rho_cfd=1.0; // no levelset 	
	double T_wt=fsi[0].r_rotor/refvel_wt; // the 0 turbine cannot be in the other turbines' wakes.  
	double T_cfd=R0/refvel_cfd; // the 0 turbine cannot be in the other turbines' wakes.  
	double ratio_rho=rho_air/rho_cfd;
	double ratio_L=reflength_wt;	
	double ratio_T=T_wt/T_cfd;	
	double ratio_V=refvel_wt/refvel_cfd;

	if(ti==tistart) {
		if (rstart_turbinerotation) {
			PetscPrintf(PETSC_COMM_WORLD, "Read turbine rotation restart file \n");
			TurbineTorqueControl_Input(fsi, ibm);
			for (ibi=0; ibi<NumberOfTurbines; ibi++) PetscPrintf(PETSC_COMM_WORLD, "angvel %f ang %f \n", fsi[ibi].angvel_axis, fsi[ibi].ang_axis);
		} 
		else {
			for (ibi=0; ibi<NumberOfTurbines; ibi++) {
				fsi[ibi].ang_axis=0.0;
				double R=fsi[ibi].r_rotor/reflength_wt;
				double TSR;
				if (FixTipSpeedRatio) {
					TSR=ibm[ibi].Tipspeedratio;
					fsi[ibi].angvel_axis=TSR*ibm[ibi].U_ref/R;
				} 
				else if (turbinetorquecontrol) {
					TSR=fsi[ibi].TSR_max;
					fsi[ibi].angvel_axis=TSR*ibm[ibi].U_ref/R;
					PetscPrintf(PETSC_COMM_WORLD, "TSR %f \n", fsi[ibi].TSR_max);
				} 
				else if (fixturbineangvel) {
					fsi[ibi].angvel_axis=fsi[ibi].angvel_fixed;
					ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;
				}
				PetscPrintf(PETSC_COMM_WORLD, "angvel of turbine at first time step %f \n", fsi[ibi].angvel_axis);
			}
		}
	}
	else {
		for (ibi=0; ibi<NumberOfTurbines; ibi++) {
			double R=fsi[ibi].r_rotor/reflength_wt;
			if (FixTipSpeedRatio) fsi[ibi].angvel_axis=ibm[ibi].Tipspeedratio*ibm[ibi].U_ref/R;
			else if (fixturbineangvel) {
				fsi[ibi].angvel_axis=fsi[ibi].angvel_fixed;
				ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;
			} 
			else if (turbinetorquecontrol){
				// J, K and dt needs normalization 	
				double J=fsi[ibi].J_rotation;
				//double mu_s=1.0/J;
				double K=0.5*M_PI*pow(fsi[ibi].r_rotor,2+X_control)*fsi[ibi].CP_max/pow(fsi[ibi].TSR_max,X_control);
				double angvel_wt = fsi[ibi].angvel_axis / ratio_T; 
				double UU=ibm[ibi].U_ref*ratio_V;
				double Torque_c, Torque_a;
				if (!fixturbinegeneratortorque) fsi[ibi].Torque_generator=K*pow(angvel_wt,X_control-1)*pow(UU,3-X_control); 
				Torque_c=fsi[ibi].Torque_generator;
				Torque_a=fsi[ibi].Torque_aero*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
				double Torque=Torque_a-Torque_c;
				double dt_wt=dt * ratio_T;
				double angvel_0=fsi[ibi].angvel_axis / ratio_T;
				double angvel_1=angvel_0 + dt_wt*Torque/J; 
				fsi[ibi].angvel_axis=angvel_1 * ratio_T;
				ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;
				double err = (angvel_1-angvel_0)*J/dt_wt-Torque;
				PetscPrintf(PETSC_COMM_WORLD, "err of solving Torque equation  %le inertial %le Torque %le \n", err, (angvel_1-angvel_0)*J/dt_wt, Torque);
			}
		}
	}

	for (ibi=0; ibi<NumberOfTurbines; ibi++) {
		double angvel_1=fsi[ibi].angvel_axis;
		fsi[ibi].ang_axis+=angvel_1*dt;
	}
	if ( ti % tiout == 0) TurbineTorqueControl_Output(fsi, ibm);
	return(0);
}