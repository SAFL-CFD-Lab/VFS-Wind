/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
extern PetscInt ti;
extern PetscReal angle;
extern PetscInt tiout;

PetscErrorCode Elmt_Move(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.52;
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

PetscErrorCode Elmt_Init_Cyl(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal rcx = -0.122, rcz = -0.32, z0 = 2.7;
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  //angle = 0;
  angle =-angle * 3.1415926/180.;
  //  angle = 0;
  for (i=0; i<n_v/2; i++) {
    // change for stat case 4/4/06 iman
    /*
    ibm->x_bp[i] = (ibm->x_bp0[i] - rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] - rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;
    */
    ibm->x_bp[i] = (ibm->x_bp0[i] );//- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = ibm->z_bp0[i];//(ibm->x_bp0[i] - rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;

  }
  rcx = -rcx;
  // change for stat case 4/3/06 iman
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] );//- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = ibm->z_bp0[i];//(ibm->x_bp0[i] - rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }
  /*
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
  */
/*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
      //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);

  //  }

  if (ti>0) {
    for (i=0; i<n_v; i++) {
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

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
}

PetscErrorCode Init_Elmt(IBMNodes *ibm, IBMNodes *ibm0, IBMNodes *ibm1)
{
  ibm->n_v = ibm0->n_v + ibm1->n_v;
  ibm->n_elmt = ibm0->n_elmt + ibm1->n_elmt;

  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
  PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
  PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

  PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->nv1));
  PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->nv2));
  PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->nv3));

  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->nf_x));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->nf_y));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->nf_z));

  PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
  PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
  
  PetscInt i;
  for (i=0; i<n_v; i++) {
    ibm->uold[i].x = 0.;
    ibm->uold[i].y = 0.;
    ibm->uold[i].z = 0.;

    ibm->u[i].x = 0.;
    ibm->u[i].y = 0.;
    ibm->u[i].z = 0.;
  }
  return 0;
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
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }

  return 0;
}

PetscErrorCode InflowWaveFormRead(UserCtx *user)
{
  PetscInt rank;
  PetscInt number=0;
  PetscReal num_beats = 70.;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *fd, *fd1;
/*     fd = fopen("inflowwave", "r"); */
    fd = fopen("inflow.dat", "r");
    FlowWave temp[3000];
    PetscInt condition = 1;
    float t, f;
    while (condition != EOF) {
      condition = fscanf(fd, "%f %f\n", &t, &f);
      if (condition != EOF) {
	temp[number].t = t; temp[number].f = f;
	number++;
      }
    }
    fclose(fd);

    MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    user->number_flowwave = number;
    PetscPrintf(PETSC_COMM_SELF, " flowwave number %d\n",user->number_flowwave );

    PetscMalloc(number*sizeof(FlowWave), &user->inflow);
    PetscInt i;
    PetscReal maxin = 0;
    for (i=0; i<number; i++) {
      maxin = PetscMax(maxin, temp[i].f);
      user->inflow[i].t = temp[i].t;// * num_beats / 60.;
/*       if (temp[i].t >= 291.5) user->inflow[i].f = 0.; */
/*       else user->inflow[i].f = temp[i].f; */

      user->inflow[i].f = temp[i].f;// > 0 ? temp[i].f : 0;
    }
    for (i=0; i<number; i++) {
      user->inflow[i].f /= maxin;
    }
    fd1 = fopen("newinflow", "w");
    for (i=0; i<number; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd1, "%e %e\n", user->inflow[i].t,
		   user->inflow[i].f);
    }
    fclose(fd1);
/*     MPI_Bcast(user->inflow, 2*number, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }
  else {
    MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    user->number_flowwave = number;
    PetscMalloc(number*sizeof(FlowWave), &user->inflow);
/*     MPI_Bcast(user->inflow, 2*number, MPIU_REAL, 0, PETSC_COMM_WORLD);     */
  }
  return 0;
}

/* PetscErrorCode InflowWaveFormRead(UserCtx *user) */
/* { */
/*   PetscInt rank; */
/*   PetscInt number=0; */
/*   PetscReal num_beats = 70.; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   if (!rank) { */
/*     FILE *fd, *fd1; */
/*     fd = fopen("inflow.dat", "r"); */
/*     FlowWave temp[3000]; */
/*     PetscInt condition = 1; */
/*     float t, f; */
/*     while (condition != EOF) { */
/*       condition = fscanf(fd, "%f %f\n", &t, &f); */
/*       if (condition != EOF) { */
/* 	temp[number].t = t; temp[number].f = f; */
/* 	number++; */
/*       } */
/*     } */
/*     fclose(fd); */

/*     MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD); */
/*     user->number_flowwave = number; */
/*     PetscMalloc(number*sizeof(FlowWave), &user->inflow); */
/*     PetscInt i; */
/*     PetscReal maxin = 0; */
/*     for (i=0; i<number; i++) { */
/*       maxin = PetscMax(maxin, temp[i].f); */
/*       user->inflow[i].t = temp[i].t;// * num_beats / 60. / 1000.; */
/*       if (temp[i].t >= 291.5) user->inflow[i].f = 0.; */
/*       else user->inflow[i].f = temp[i].f; */
/*     } */
/*     for (i=0; i<number; i++) { */
/*       user->inflow[i].f /= maxin; */
/*     } */
/*     fd1 = fopen("newinflow", "w"); */
/*     for (i=0; i<number; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd1, "%e %e\n", user->inflow[i].t, */
/* 		   user->inflow[i].f); */
/*     } */
/*     fclose(fd1); */

/*     fd = fopen("leaflet.dat", "r"); */
/*     number = 0; */
/*     condition = 1; */
/*     while (condition != EOF) { */
/*       condition = fscanf(fd, "%f %f\n", &t, &f); */
/*       if (condition != EOF) { */
/* 	temp[number].t = t; temp[number].f = f; */
/* 	number++; */
/*       } */
/*     } */
    
/*     fclose(fd); */
/*     MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD); */
/*     user->number_kinematics = number; */

/*     PetscMalloc(number*sizeof(FlowWave), &user->kinematics); */

/*     for (i=0; i<number; i++) { */
/*       user->kinematics[i].t = temp[i].t + 10; */
/*       user->kinematics[i].f = temp[i].f; */
/*     } */
/*     PetscPrintf(PETSC_COMM_WORLD, "Read Kinematics %i\n", user->number_kinematics); */

/*     MPI_Bcast(user->inflow, 2*number, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*   } */
/*   else { */
/*     MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD); */
/*     user->number_flowwave = number; */
/*     PetscMalloc(number*sizeof(FlowWave), &user->inflow); */

/*     MPI_Bcast(&number, 1, MPI_INT, 0, PETSC_COMM_WORLD); */
/*     user->number_kinematics = number; */

/*     MPI_Bcast(user->inflow, 2*number, MPIU_REAL, 0, PETSC_COMM_WORLD);     */
/*   } */
/*   PetscPrintf(PETSC_COMM_SELF, "Flow Wave %i\n", rank); */

/*   return 0; */
/* } */
