/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c, L_dim;
extern char path[256];

PetscErrorCode ibm_surface_out(IBMNodes *ibm, PetscInt ti,
			       PetscInt ibi)
{
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "%s/surface%3.3d_%2.2d.dat",path,ti,ibi);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", ibm->n_v, ibm->n_elmt);
      for (i=0; i<ibm->n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
  return(0);
}

PetscErrorCode ibm_read(IBMNodes *ibm1, IBMNodes *ibm2)
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
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscReal     cl=30.;

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
    ibm1->n_v = n_v/2;
    ibm2->n_v = n_v/2;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v; i++) {
      fscanf(fd, "%le %le %le %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i], &t, &t, &t);
      x_bp[i] = x_bp[i] / cl;// 28.
      y_bp[i] = y_bp[i] / cl;
      z_bp[i] = z_bp[i] / cl;
    }
/*     ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */   

    PetscReal temp;
    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i];
      ibm1->x_bp0[i] = z_bp[i];
      ibm1->y_bp0[i] = x_bp[i]-CMy_c;
      ibm1->z_bp0[i] = -temp+CMz_c;

      ibm1->x_bp[i] = z_bp[i];
      ibm1->y_bp[i] = x_bp[i]-CMy_c;
      ibm1->z_bp[i] = -temp+CMz_c;

      ibm1->x_bp_o[i] = z_bp[i];
      ibm1->y_bp_o[i] = x_bp[i]-CMy_c;
      ibm1->z_bp_o[i] = -temp+CMz_c;   

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   
    }

    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i+n_v/2];
      ibm2->x_bp0[i] = z_bp[i+n_v/2];
      ibm2->y_bp0[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp0[i] = -temp+CMz_c;

      ibm2->x_bp[i] = z_bp[i+n_v/2];
      ibm2->y_bp[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp[i] = -temp+CMz_c;

      ibm2->x_bp_o[i] = z_bp[i+n_v/2];
      ibm2->y_bp_o[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp_o[i] = -temp+CMz_c;

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   
    }

/*     ibm1->x_bp=ibm1->x_bp0;ibm1->y_bp=ibm1->y_bp0;ibm1->z_bp=ibm1->z_bp0; */
/*     ibm2->x_bp=ibm2->x_bp0;ibm2->y_bp=ibm2->y_bp0;ibm2->z_bp=ibm2->z_bp0; */
/*     ibm1->x_bp_o=ibm1->x_bp0;ibm1->y_bp_o=ibm1->y_bp0;ibm1->z_bp_o=ibm1->z_bp0; */
/*     ibm2->x_bp_o=ibm2->x_bp0;ibm2->y_bp_o=ibm2->y_bp0;ibm2->z_bp_o=ibm2->z_bp0; */

    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fscanf(fd, "%i\n", &n_elmt);
    ibm1->n_elmt = n_elmt/2;
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm2->n_elmt = n_elmt/2;
    //MPI_Bcast(&(ibm2->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "elmts , %d \n", n_elmt);

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

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
    
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }

    //ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

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

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) 
      {
	if (nf_z[i]>0) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
	} else {
	ns_x[i] = -1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = -1.;
	nt_z[i] = 0.;
	}	  
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
/*       ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.; */
/*       ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.; */
/*       ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	 */
    }
   
    //ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    for (i=0; i<n_elmt/2; i++) {
      ibm1->nv1[i]=nv1[i];
      ibm1->nv2[i]=nv2[i];
      ibm1->nv3[i]=nv3[i];

      ibm1->nf_x[i]=  nf_z[i];
      ibm1->nf_y[i]=  nf_x[i];
      ibm1->nf_z[i]= -nf_y[i];

      ibm1->ns_x[i]=  ns_z[i];
      ibm1->ns_y[i]=  ns_x[i];
      ibm1->ns_z[i]= -ns_y[i];

      ibm1->nt_x[i]=  nt_z[i];
      ibm1->nt_y[i]=  nt_x[i];
      ibm1->nt_z[i]= -nt_y[i];
    }

    for (i=0; i<n_elmt/2; i++) {
      ibm2->nv1[i]=nv1[i+n_elmt/2]-n_v/2;
      ibm2->nv2[i]=nv2[i+n_elmt/2]-n_v/2;
      ibm2->nv3[i]=nv3[i+n_elmt/2]-n_v/2;

      ibm2->nf_x[i]=  nf_z[i+n_elmt/2];
      ibm2->nf_y[i]=  nf_x[i+n_elmt/2];
      ibm2->nf_z[i]= -nf_y[i+n_elmt/2];

      ibm2->ns_x[i]=  ns_z[i+n_elmt/2];
      ibm2->ns_y[i]=  ns_x[i+n_elmt/2];
      ibm2->ns_z[i]= -ns_y[i+n_elmt/2];

      ibm2->nt_x[i]=  nt_z[i+n_elmt/2];
      ibm2->nt_y[i]=  nt_x[i+n_elmt/2];
      ibm2->nt_z[i]= -nt_y[i+n_elmt/2];

    }

/*     for (i=0; i<n_elmt; i++) { */
/*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
/*     } */

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_v = n_v/2;

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v/2; i++) {

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   

    }

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     ibm1->x_bp=ibm1->x_bp0;ibm1->y_bp=ibm1->y_bp0;ibm1->z_bp=ibm1->z_bp0; */
/*     ibm2->x_bp=ibm2->x_bp0;ibm2->y_bp=ibm2->y_bp0;ibm2->z_bp=ibm2->z_bp0; */
/*     ibm1->x_bp_o=ibm1->x_bp0;ibm1->y_bp_o=ibm1->y_bp0;ibm1->z_bp_o=ibm1->z_bp0; */
/*     ibm2->x_bp_o=ibm2->x_bp0;ibm2->y_bp_o=ibm2->y_bp0;ibm2->z_bp_o=ibm2->z_bp0; */


    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_elmt = n_elmt/2;
    ibm2->n_elmt = n_elmt/2;
    
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

/*     ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3; */
/*     ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z; */


/*     for (i=0; i<n_elmt/2; i++) { */
/*       ibm1->nv1[i]=nv1[i]; */
/*       ibm1->nv2[i]=nv2[i]; */
/*       ibm1->nv3[i]=nv3[i]; */

/*       ibm1->nf_x[i]= nf_x[i]; */
/*       ibm1->nf_y[i]= nf_y[i]; *//*       ibm1->nf_z[i]= nf_z[i]; */
/*     } */

/*     for (i=0; i<n_elmt/2; i++) { */
/*       ibm2->nv1[i]=nv1[i+n_elmt/2]; */
/*       ibm2->nv2[i]=nv2[i+n_elmt/2]; */
/*       ibm2->nv3[i]=nv3[i+n_elmt/2]; */

/*       ibm2->nf_x[i]= nf_x[i+n_elmt/2]; */
/*       ibm2->nf_y[i]= nf_y[i+n_elmt/2]; */
/*       ibm2->nf_z[i]= nf_z[i+n_elmt/2]; */
/*     } */

/*     for (i=0; i<n_elmt; i++) { */
    i=10;
      PetscPrintf(PETSC_COMM_SELF, " ibm1 xbp %d %le %le %le\n", i, ibm1->z_bp0[i], ibm1->z_bp_o[i], ibm1->z_bp[i]);
/*     } */

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

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


//removed by seokkoo
/*
PetscErrorCode ibm_read_thin(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
//   PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270]; 
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
    fd = fopen("NodeData_A1_m", "r");
    
    n_v = 3270;
    ibm->n_v = n_v;

    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscInt t1, j;
    for (j=0; j<101; j++)
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%d %le %le %le %le", &t1, &(ibm->x_bp_in[j][i]),
	       &(ibm->y_bp_in[j][i]), &(ibm->z_bp_in[j][i]), &t);
	ibm->x_bp_in[j][i] = ibm->x_bp_in[j][i] / 24.;
	ibm->y_bp_in[j][i] = ibm->y_bp_in[j][i] / 24.;
	ibm->z_bp_in[j][i] = ibm->z_bp_in[j][i] / 24.;
    }

    for (j=45; j<50; j++) {
      for (i=0; i<n_v; i++) {
	ibm->x_bp_in[j][i] = ibm->x_bp_in[50-j][i];
	ibm->y_bp_in[j][i] = ibm->y_bp_in[50-j][i];
	ibm->z_bp_in[j][i] = ibm->z_bp_in[50-j][i];
      }
    }
	
    MPI_Bcast(ibm->x_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fclose(fd);

    n_elmt = 6228;
    ibm->n_elmt = n_elmt;
    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_z);
    fd = fopen("Elements_id_A1", "r");
    
    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i, %i, %i, %i\n", &t1, nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }
    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

    fclose(fd);

    for (i=0; i<n_elmt; i++) {
      nv1[i+n_elmt] = nv1[i];
      nv2[i+n_elmt] = nv3[i];
      nv3[i+n_elmt] = nv2[i];
    }

    n_elmt = n_elmt*2;
    ibm->n_elmt = n_elmt;

    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    n_v = ibm->n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

    MPI_Bcast(ibm->x_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    n_elmt = ibm->n_elmt;


    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    n_elmt = 2 * n_elmt;
    ibm->n_elmt = n_elmt;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }


  return(0);
}
*/
PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;

  char   ss[20];
  //double xt;
  char string[128];

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ ibmdata\n");
    char filen[80];  
    sprintf(filen,"%s/ibmdata%2.2d" , path, ibi);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
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
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
	    
	    /*
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);	// removed by seokkoo 03.04.2009
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
      */
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      x_bp = ibm->x_bp;	// seokkoo
      y_bp = ibm->y_bp;
      z_bp = ibm->z_bp;
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
	x_bp[i] = x_bp[i]*L_dim + CMx_c;//0.25 ;// 24.;	
	y_bp[i] = y_bp[i]*L_dim + CMy_c ;//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
	z_bp[i] = z_bp[i]*L_dim + CMz_c ;//+ ibi*1.5;//2.;//8.;//15.;//2.   ;// 24.;

	//Beginning of initial rotation
	double z_bp_,y_bp_,x_bp_;//Temporary variables for initial rotation
	//Center of rotation (_r) and initial angle (angle_ 0) sfecified at control file
	// initial rotation around x-axis
	y_bp_ = y_r + (y_bp[i]-y_r)*cos(angle_x0) - (z_bp[i]-z_r)*sin(angle_x0);
	z_bp_ = z_r + (y_bp[i]-y_r)*sin(angle_x0) + (z_bp[i]-z_r)*cos(angle_x0);
	y_bp[i] = y_bp_;
	z_bp[i] = z_bp_;	
	// initial rotation around y-axis
	x_bp_ = x_r + (z_bp[i]-z_r)*sin(angle_y0) + (x_bp[i]-x_r)*cos(angle_y0);		
	z_bp_ = z_r + (z_bp[i]-z_r)*cos(angle_y0) - (x_bp[i]-x_r)*sin(angle_y0);
	x_bp[i] = x_bp_;
	z_bp[i] = z_bp_;
	// initial rotation around z-axis
	x_bp_ = x_r + (x_bp[i]-x_r)*cos(angle_z0) - (y_bp[i]-y_r)*sin(angle_z0);
	y_bp_ = y_r + (x_bp[i]-x_r)*sin(angle_z0) + (y_bp[i]-y_r)*cos(angle_z0);		
	x_bp[i] = x_bp_;
	y_bp[i] = y_bp_;
	//End of initial rotation
	

		
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->x_bp0[i] = x_bp[i];
	ibm->y_bp0[i] = y_bp[i];
	ibm->z_bp0[i] = z_bp[i];

	ibm->x_bp_o[i] = x_bp[i];
	ibm->y_bp_o[i] = y_bp[i];
	ibm->z_bp_o[i] = z_bp[i];

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

/*       ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */
/*       ibm->x_bp_o = x_bp; ibm->y_bp_o = y_bp; ibm->z_bp_o = z_bp; */

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

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
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
      
      
	//seokkoo begin
	{	//only for rank 0
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv1);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv2);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv3);
		
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_x_bp);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_y_bp);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_z_bp);
	}
	
	PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->count);
	PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->local2global_elmt);
	
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);
	//seokkoo end

      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, &ss, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
	      
		// seokkoo
	      ibm->_nv1[i] = nv1[i];
	      ibm->_nv2[i] = nv2[i];
	      ibm->_nv3[i] = nv3[i];
	      // seokkoo
	      ibm->_x_bp[i] = ibm->x_bp[i];
	      ibm->_y_bp[i] = ibm->y_bp[i];
	      ibm->_z_bp[i] = ibm->z_bp[i];

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
      
      // Temp sol. 2D
/*       if (fabs(nf_x[i])<.5) */
/* 	nf_x[i]=0.; */

      // Addedd 4/2/06 iman
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
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
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

 /*    PetscFree(dA); */
/*     PetscFree(nf_x);PetscFree(nf_y);PetscFree(nf_z); */
/*     PetscFree(nt_x);PetscFree(nt_y);PetscFree(nt_z); */
/*     PetscFree(ns_x);PetscFree(ns_y);PetscFree(ns_z); */
/*     PetscFree(nv1);PetscFree(nv2);PetscFree(nv3); */
/*     PetscFree(x_bp);PetscFree(y_bp);PetscFree(z_bp); */
    PetscInt ti=0;
    FILE *f;
    //char filen[80];
    sprintf(filen, "%s/surface%3.3d_%2.2d_nf.dat",path,ti,ibi);
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
    
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt); */
/*     for (i=0; i<n_v; i++) { */
      
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*     } */
    fclose(f);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    /*
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);	// removed by seokkoo 03.04.2009
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    */
	
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
	  
	x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
	y_bp = ibm->y_bp;
	z_bp = ibm->z_bp;
	  
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

    for (i=0; i<n_v; i++) {
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

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

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
    
		//seokkoo
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->count);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->local2global_elmt);
		//seokkoo
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);

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
  }
	PetscPrintf(PETSC_COMM_WORLD, "Read ucd file !\n");
  return(0);
}

PetscErrorCode ibm_read_ucd_old(IBMNodes *ibm)
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
	x_bp[i] = (x_bp[i]    ) / cl;
	y_bp[i] = (y_bp[i]+6. ) / cl ;
	z_bp[i] = (z_bp[i]+15.) / cl;
	
/* 	ibm->x_bp[i] = x_bp[i]; */
/* 	ibm->y_bp[i] = y_bp[i]; */
/* 	ibm->z_bp[i] = z_bp[i]; */

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
