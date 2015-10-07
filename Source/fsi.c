/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
extern PetscInt radi,fish,TwoD;
extern char path[256];

PetscErrorCode linear_intp(Cpt2D p, Cpt2D p1, Cpt2D p2, 
			   IBMInfo *ibminfo, PetscInt number,
			   PetscInt nvert)
{
  PetscReal  x12, y12, xp2, yp2, Cr;
  x12 = p1.x - p2.x; y12 = p1.y - p2.y;
  xp2 = p.x - p2.x; yp2 = p.y - p2.y;

  if (fabs(x12)>1e-7) {
    Cr=xp2/x12;
  }  else if (fabs(y12)>1e-7) {
    Cr=yp2/y12;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "%Something Wrong!!! Linear intp two points are the same!!!\n");
  }
  
  if (nvert==1) {
    ibminfo[number].cr1 = 0.;
    ibminfo[number].cr2 = Cr;    
    ibminfo[number].cr3 = 1-Cr;
  } else if (nvert==2) {
    ibminfo[number].cr1 = Cr;
    ibminfo[number].cr2 = 0.;    
    ibminfo[number].cr3 = 1-Cr;
  } else if (nvert==3) {
    ibminfo[number].cr1 = Cr;
    ibminfo[number].cr2 = 1-Cr;    
    ibminfo[number].cr3 = 0.;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "%Wrong Nvert in Linear intp!!!\n");
  }
  return(0);  
}

PetscErrorCode triangle_intpp2_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].cs11 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].cs22 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].cs33 = 1. - ibminfo[number].cs11 - ibminfo[number].cs22;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].cs33<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Cs %le %le %le\n", ibminfo[number].cs33, ibminfo[number].cs22, ibminfo[number].cs11);

}

PetscErrorCode FsiInitialize(PetscInt n_elmt, FSInfo *fsi,PetscInt ibi)
{
  /*  Note: This routine shouldn't be called before ibm_read */
  PetscInt  i;

  PetscMalloc(n_elmt*sizeof(IBMInfo), &(fsi->fsi_intp));
  PetscMalloc(n_elmt*sizeof(SurfElmtInfo), &(fsi->elmtinfo));

  for (i=0;i<6;i++) {
    fsi->S_old[i]=0.;
    fsi->S_new[i]=0.;
    fsi->S_realm1[i]=0.;
    fsi->S_real[i]=0.;

    fsi->S_ang_n[i]=0.;
    fsi->S_ang_o[i]=0.;
    fsi->S_ang_r[i]=0.;
    fsi->S_ang_rm1[i]=0.;
  }

/*   fsi->S_new[4]= */
  
  /*  
  fsi->S_new[4]=-4.721795e-03;
  fsi->S_old[4]=0.;
  fsi->S_real[4]=-4.646872e-03;
  fsi->S_realm1[4]=-4.572523e-03;
  
  fsi->S_new[5]=-7.523306e-03;
  fsi->S_old[5]=-7.523306e-03;
  fsi->S_real[5]=-7.465149e-03;
  fsi->S_realm1[5]=-7.465149e-03;
  */
  fsi->F_x_old = 0.;
  fsi->F_y_old = 0.;
  fsi->F_z_old = 0.;

  fsi->F_x_real = 0.;
  fsi->F_y_real = 0.;
  fsi->F_z_real = 0.;

  fsi->M_x_old = 0.;
  fsi->M_y_old = 0.;
  fsi->M_z_old = 0.;

  fsi->M_x_real = 0.;
  fsi->M_y_real = 0.;
  fsi->M_z_real = 0.;

  fsi->x_c=0.0;fsi->y_c=0.;fsi->z_c=0.;

  fsi->red_vel=0.52;//1.5;//0.3921;
  fsi->damp=.02;
  fsi->mu_s=500.;//0.025;//0.00568;

  fsi->Max_xbc= 1e23;fsi->Max_ybc= 1e23;fsi->Max_zbc= 1e23;
  fsi->Min_xbc=-1e23;fsi->Min_ybc=-1e23;fsi->Min_zbc=-1e23;
  

  PetscOptionsGetReal(PETSC_NULL, "-red_vel", &(fsi->red_vel), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-damp", &(fsi->damp), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-mu_s", &(fsi->mu_s), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-x_c", &(fsi->x_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-y_c", &(fsi->y_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-z_c", &(fsi->z_c), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-Max_xbc", &(fsi->Max_xbc), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-Min_xbd", &(fsi->Min_xbc), PETSC_NULL);

/*   Check Later!!!!!!!!!!!!!!!! */
  //fsi->y_c += 0.5;
  //fsi->x_c += 0.5;
 // fsi->z_c +=ibi*1.5;
/*   fsi->S_new[4]= fsi->z_c; */
/*   fsi->S_new[2]= fsi->y_c; */

  return(0);
}

PetscErrorCode SetPressure(UserCtx *user) 
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt      gxs, gxe, gys, gye, gzs, gze; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  PetscInt      i,j,k;
  Cmpnts        ***coor, ***ucat;
  PetscReal     ***p;

  DAVecGetArray(fda, user->lCent,&coor);
  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(da, user->lP, &p);

  for (k=gzs; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	p[k][j][i]=4.*coor[k][j][i].y;
	//p[k][j][i]=16.*coor[k][j][i].y*coor[k][j][i].y;
	ucat[k][j][i].y=1.;
	//if (j==23)  PetscPrintf(PETSC_COMM_SELF, "%le %le\n",coor[k][j][i].y, p[k][j][i]);
      }
    }
  }
  
  DAVecRestoreArray(fda, user->lCent,&coor);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lP, &p);
  return(0);
}

PetscErrorCode Closest_NearBndryPt_ToSurfElmt(UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo, FSInfo *fsi, PetscInt ibi)
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  //PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

/*   if (xs==0) lxs = xs+1; */
/*   if (ys==0) lys = ys+1; */
/*   if (zs==0) lzs = zs+1; */

/*   if (xe==mx) lxe = xe-1; */
/*   if (ye==my) lye = ye-1; */
/*   if (ze==mz) lze = ze-1; */

  PetscInt	i, j, k;

  //PetscInt      nbnumber = user->ibmnumber;
  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt, nbn;
  Cmpnts        ***coor, *nbncoor;
  PetscReal     xc,yc,zc; // tri shape surface elmt center coords
  PetscReal     nfx,nfy,nfz;
  PetscReal     x,y,z;    // near bndry pt coords
  PetscReal     dist,*distMin, dmin;     // distance between surf elmt to near bndry pt

  IBMInfo       *ibminfo;
  IBMListNode   *current;

/*   PetscMalloc(nbnumber*sizeof(PetscReal), &dist); */
/*   PetscMalloc(nbnumber*sizeof(PetscReal), &distMin); */
/*   PetscMalloc(nbnumber*sizeof(Cmpnts), &nbncoor); */
  DAVecGetArray(fda, user->Cent,&coor);

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  //PetscPrintf(PETSC_COMM_SELF, "n_elmt, nbnumber %d %d\n",n_elmt, nbnumber);
 
  for (elmt=0; elmt<n_elmt; elmt++) {
    xc=ibm->cent_x[elmt]; 
    yc=ibm->cent_y[elmt]; 
    zc=ibm->cent_z[elmt]; 

    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    dmin=1e10;

    if (xc>=fsi->Min_xbc && xc<=fsi->Max_xbc) { //elmt inside the fluid domain
    //if (xc>=0.022 && xc<=0.078) { //elmt inside the fluid domain
    //if (xc>=0.1 && xc<=0.3) { //elmt inside the fluid domain
    //if (zc>4.7 && zc<4.7515 || yc>.1 && yc<.25005 || xc>.1 && xc<.25005) { //elmt inside the fluid domain
    //if (zc>0.0000 && yc>0.0000  && xc>0.0000) {
    //if (nfy>=0.0000  && nfz>=0.0000) {
      elmtinfo[elmt].n_P=1;
      //for (nbn=0; nbn<nbnumber; nbn ++) {
      current = user->ibmlist[ibi].head;
      while (current) {
	ibminfo = &current->ibm_intp;
	current = current->next;	
	//i=ibminfo[nbn].ni; j=ibminfo[nbn].nj; k=ibminfo[nbn].nk;
	i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      
	if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze) {	        
	  x=coor[k][j][i].x;
	  y=coor[k][j][i].y;
	  z=coor[k][j][i].z;
	
	  dist=(x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc);

	if (dmin>dist) {
	  dmin=dist;
	  elmtinfo[elmt].Clsnbpt_i=i;
	  elmtinfo[elmt].Clsnbpt_j=j;
	  elmtinfo[elmt].Clsnbpt_k=k;

	}
	}
      }      
/*       PetscPrintf(PETSC_COMM_SELF, "%d %d %d %d %d %le %le %le\n",elmt,elmtinfo[elmt].n_P,elmtinfo[elmt].Clsnbpt_i, */
/* 		elmtinfo[elmt].Clsnbpt_j,elmtinfo[elmt].Clsnbpt_k,xc,yc,zc ); */

    }
  }
    
  DAVecRestoreArray(fda, user->Cent,&coor);
  return(0);
}

PetscErrorCode distance(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;
  
  PetscReal dx1, dy1, dz1, dx2, dy2, dz2, r;

  dx1 = p3.x - p1.x;
  dy1 = p3.y - p1.y;
  dz1 = p3.z - p1.z;
 
  dx2 = p4.x - p2.x;
  dy2 = p4.y - p2.y;
  dz2 = p4.z - p2.z;

  xn1 = dy1 * dz2 - dz1 * dy2;
  yn1 = - (dx1 * dz2 - dz1 * dx2);
  zn1 = dx1 * dy2 - dy1 * dx2;

  r = sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1);
  xn1 /= r; yn1 /= r; zn1 /= r;

  xc = 0.25 * (p1.x + p2.x + p3.x + p4.x);
  yc = 0.25 * (p1.y + p2.y + p3.y + p4.y);
  zc = 0.25 * (p1.z + p2.z + p3.z + p4.z);

  *d = (p.x - xc) * xn1 + (p.y - yc) * yn1 + (p.z - zc) * zn1;
  if (PetscAbsReal(*d)<1.e-6) *d=0.;
  return (0);
}

PetscTruth ISInsideCell(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{
  // k direction
  distance(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  if (d[4]<0) return(PETSC_FALSE);
  distance(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));
  if (d[5]<0) return(PETSC_FALSE);

  // j direction
  distance(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  if (d[2]<0) return(PETSC_FALSE);

  distance(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));
  if (d[3]<0) return(PETSC_FALSE);

  // i direction
  distance(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  if (d[0]<0) return(PETSC_FALSE);
  
  distance(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));
  if (d[1]<0) return(PETSC_FALSE);
  return(PETSC_TRUE);
}

PetscErrorCode GridCellaroundSurElmt(UserCtx *user, IBMNodes *ibm,SurfElmtInfo *elmtinfo)
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-2;
  if (ye==my) lye = ye-2;
  if (ze==mz) lze = ze-2;

  PetscInt      zstart,zend,ystart,yend,xstart,xend;
  PetscInt	i, j, k, inbn,jnbn,knbn, notfound;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt, nradius=20;
  Cmpnts        ***coor,pc,cell[8];
  PetscReal     d[6];
  PetscReal     AroundCellSum,Aroundcellnum;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  Aroundcellnum=0;
  for (elmt=0; elmt<n_elmt; elmt++) {
  if (elmtinfo[elmt].n_P>0){ 
    pc.x=ibm->cent_x[elmt]; 
    pc.y=ibm->cent_y[elmt]; 
    pc.z=ibm->cent_z[elmt]; 

    inbn=elmtinfo[elmt].Clsnbpt_i;
    jnbn=elmtinfo[elmt].Clsnbpt_j;
    knbn=elmtinfo[elmt].Clsnbpt_k;

    zstart=knbn-nradius; zend=knbn+nradius;
    ystart=jnbn-nradius; yend=jnbn+nradius;
    xstart=inbn-nradius; xend=inbn+nradius;

    if (zstart<lzs) zstart=lzs;
    if (ystart<lys) ystart=lys;
    if (xstart<lxs) xstart=lxs;
    if (zend>lze) zend=lze;
    if (yend>lye) yend=lye;
    if (xend>lxe) xend=lxe;

    notfound=0;
    elmtinfo[elmt].FoundAroundcell=0;

    for (k=zstart; k<zend; k++) {
      for (j=ystart; j<yend; j++) {
	for (i=xstart; i<xend; i++) {
	  cell[0] = coor[k  ][j  ][i  ];
	  cell[1] = coor[k  ][j  ][i+1];
	  cell[2] = coor[k  ][j+1][i+1];
	  cell[3] = coor[k  ][j+1][i  ];
	  
	  cell[4] = coor[k+1][j  ][i  ];
	  cell[5] = coor[k+1][j  ][i+1];
	  cell[6] = coor[k+1][j+1][i+1];
	  cell[7] = coor[k+1][j+1][i  ];
	
	  if(ISInsideCell(pc, cell, d)){
	    elmtinfo[elmt].icell=i;
	    elmtinfo[elmt].jcell=j;
	    elmtinfo[elmt].kcell=k;
	    elmtinfo[elmt].FoundAroundcell=1;
	    //if (j==22)  PetscPrintf(PETSC_COMM_SELF, "%le %le\n",coor[k][j][i].y, coor[k][j+1][i].y);
	    // correction if pt exactly on one side of the cell
	    if (fabs(d[0])<1e-8 && 
		ibm->nf_x[elmt]*(cell[1].x-cell[0].x)<0.) elmtinfo[elmt].icell=i-1;
	    if (fabs(d[1])<1e-8 && 
		ibm->nf_x[elmt]*(cell[0].x-cell[1].x)<0.) elmtinfo[elmt].icell=i+1;
	    if (fabs(d[2])<1e-8 && 
		ibm->nf_y[elmt]*(cell[3].y-cell[0].y)<0.) elmtinfo[elmt].jcell=j-1;
	    if (fabs(d[3])<1e-8 && 
		ibm->nf_y[elmt]*(cell[0].y-cell[3].y)<0.) elmtinfo[elmt].jcell=j+1;
	    if (fabs(d[4])<1e-8 && 
		ibm->nf_z[elmt]*(cell[4].z-cell[0].z)<0.) elmtinfo[elmt].kcell=k-1;
	    if (fabs(d[5])<1e-8 && 
		ibm->nf_z[elmt]*(cell[0].z-cell[4].z)<0.) elmtinfo[elmt].kcell=k+1;

	    Aroundcellnum+=1.;
	    notfound=1;
	  }
	}
      }	
    }

/*      PetscPrintf(PETSC_COMM_SELF, "%d %d %d %d %d %le %le %le\n",elmt,elmtinfo[elmt].icell,elmtinfo[elmt].jcell,  */
/* 		 elmtinfo[elmt].kcell,rank,pc.x,pc.y,pc.z);  */

  }
  }
  PetscGlobalSum(&Aroundcellnum, &AroundCellSum, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "n_elmt & number of cells %d %le\n",n_elmt,AroundCellSum);

  DAVecRestoreArray(fda, user->lCent,&coor);

  return(0);
}

PetscErrorCode GridCellaround2ndElmt(UserCtx *user, IBMNodes *ibm,
				     Cmpnts pc,PetscInt elmt,
				     PetscInt knbn,PetscInt jnbn,
				     PetscInt inbn, PetscInt *kin,
				     PetscInt *jin, PetscInt *iin)
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-2;
  if (ye==my) lye = ye-2;
  if (ze==mz) lze = ze-2;

  PetscInt      zstart,zend,ystart,yend,xstart,xend;
  PetscInt	i, j, k, notfound;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      nradius=3;
  Cmpnts        ***coor,cell[8];
  PetscReal     d[6];
  PetscReal     AroundCellSum,Aroundcellnum;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  Aroundcellnum=0;
/*   for (elmt=0; elmt<n_elmt; elmt++) { */
/*   if (elmtinfo[elmt].n_P>0){  */
/*     pc.x=ibm->cent_x[elmt];  */
/*     pc.y=ibm->cent_y[elmt];  */
/*     pc.z=ibm->cent_z[elmt];  */

/*     inbn=elmtinfo[elmt].Clsnbpt_i; */
/*     jnbn=elmtinfo[elmt].Clsnbpt_j; */
/*     knbn=elmtinfo[elmt].Clsnbpt_k; */

    zstart=knbn-nradius; zend=knbn+nradius;
    ystart=jnbn-nradius; yend=jnbn+nradius;
    xstart=inbn-nradius; xend=inbn+nradius;

    if (zstart<lzs) zstart=lzs;
    if (ystart<lys) ystart=lys;
    if (xstart<lxs) xstart=lxs;
    if (zend>lze) zend=lze;
    if (yend>lye) yend=lye;
    if (xend>lxe) xend=lxe;

    notfound=0;
/*     elmtinfo[elmt].FoundAroundcell=0; */

    for (k=zstart; k<zend; k++) {
      for (j=ystart; j<yend; j++) {
	for (i=xstart; i<xend; i++) {
	  cell[0] = coor[k  ][j  ][i  ];
	  cell[1] = coor[k  ][j  ][i+1];
	  cell[2] = coor[k  ][j+1][i+1];
	  cell[3] = coor[k  ][j+1][i  ];
	  
	  cell[4] = coor[k+1][j  ][i  ];
	  cell[5] = coor[k+1][j  ][i+1];
	  cell[6] = coor[k+1][j+1][i+1];
	  cell[7] = coor[k+1][j+1][i  ];
	
	  if(ISInsideCell(pc, cell, d)){
	    *kin=k;
	    *jin=j;
	    *iin=i;
/* 	    elmtinfo[elmt].icell=i; */
/* 	    elmtinfo[elmt].jcell=j; */
/* 	    elmtinfo[elmt].kcell=k; */
/* 	    elmtinfo[elmt].FoundAroundcell=1; */
	    //if (j==22)  PetscPrintf(PETSC_COMM_SELF, "%le %le\n",coor[k][j][i].y, coor[k][j+1][i].y);
	    // correction if pt exactly on one side of the cell
	    if (fabs(d[0])<1e-8 && 
		ibm->nf_x[elmt]*(cell[1].x-cell[0].x)<0.) *iin=i-1;//elmtinfo[elmt].icell=i-1;
	    if (fabs(d[1])<1e-8 && 
		ibm->nf_x[elmt]*(cell[0].x-cell[1].x)<0.) *iin=i+1;//elmtinfo[elmt].icell=i+1;
	    if (fabs(d[2])<1e-8 && 
		ibm->nf_y[elmt]*(cell[3].y-cell[0].y)<0.) *jin=j-1;//elmtinfo[elmt].jcell=j-1;
	    if (fabs(d[3])<1e-8 && 
		ibm->nf_y[elmt]*(cell[0].y-cell[3].y)<0.) *jin=j+1;//elmtinfo[elmt].jcell=j+1;
	    if (fabs(d[4])<1e-8 && 
		ibm->nf_z[elmt]*(cell[4].z-cell[0].z)<0.) *kin=k-1;//elmtinfo[elmt].kcell=k-1;
	    if (fabs(d[5])<1e-8 && 
		ibm->nf_z[elmt]*(cell[0].z-cell[4].z)<0.) *kin=k+1;//elmtinfo[elmt].kcell=k+1;

	    Aroundcellnum+=1.;
	    notfound=1;
	    break;
	  }
	}
      }	
    }

    if (!notfound) PetscPrintf(PETSC_COMM_SELF, "2nd Around Cell WAS NOT FOUND!!!!!!!!!!!! %d %d %d %d\n", elmt,inbn,jnbn,knbn);
/*      PetscPrintf(PETSC_COMM_SELF, "%d %d %d %d %d %le %le %le\n",elmt,elmtinfo[elmt].icell,elmtinfo[elmt].jcell,  */
/* 		 elmtinfo[elmt].kcell,rank,pc.x,pc.y,pc.z);  */

/*   } */
/*   } */
/*   PetscGlobalSum(&Aroundcellnum, &AroundCellSum, PETSC_COMM_WORLD); */
/*   PetscPrintf(PETSC_COMM_WORLD, "n_elmt & number of cells %d %le\n",n_elmt,AroundCellSum); */

  DAVecRestoreArray(fda, user->lCent,&coor);

  return(0);
}


PetscErrorCode triangle_intp_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].cr1 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].cr2 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].cr3 = 1. - ibminfo[number].cr1 - ibminfo[number].cr2;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].cr3<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Cr %d %le %le %le\n", number,ibminfo[number].cr3, ibminfo[number].cr2, ibminfo[number].cr1);
  
}

PetscErrorCode triangle_intp2_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].cs1 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].cs2 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].cs3 = 1. - ibminfo[number].cs1 - ibminfo[number].cs2;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].cs3<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Cs %d %le %le %le\n",number, ibminfo[number].cs3, ibminfo[number].cs2, ibminfo[number].cs1);

}



PetscErrorCode fsi_InterceptionPoint(Cmpnts p, Cmpnts pc[8], 
	       PetscReal nvertpc[8], PetscReal nfx, PetscReal nfy,
	       PetscReal nfz, IBMInfo *ibminfo, PetscInt number, 
	       Cmpnts *intp, PetscInt *Need3rdPoint)
{
  PetscInt 	triangles[3][12];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;//, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo[number].imode = -100;
  // k-plane
  triangles[0][0]  = 0; triangles[1][0]  = 1; triangles[2][0]  = 2;
  triangles[0][1]  = 0; triangles[1][1]  = 2; triangles[2][1]  = 3;
  triangles[0][2]  = 4; triangles[1][2]  = 5; triangles[2][2]  = 6;
  triangles[0][3]  = 4; triangles[1][3]  = 6; triangles[2][3]  = 7;
  // i-plane
  triangles[0][4]  = 0; triangles[1][4]  = 4; triangles[2][4]  = 3;
  triangles[0][5]  = 3; triangles[1][5]  = 4; triangles[2][5]  = 7;
  triangles[0][6]  = 1; triangles[1][6]  = 5; triangles[2][6]  = 2;
  triangles[0][7]  = 2; triangles[1][7]  = 6; triangles[2][7]  = 5;
  // j-plane
  triangles[0][8]  = 0; triangles[1][8]  = 4; triangles[2][8]  = 1;
  triangles[0][9]  = 1; triangles[1][9]  = 4; triangles[2][9]  = 5;
  triangles[0][10] = 3; triangles[1][10] = 7; triangles[2][10] = 2;
  triangles[0][11] = 2; triangles[1][11] = 7; triangles[2][11] = 6;

  for (i=0; i<12; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;   //a1=p -p1
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;//a2=p2-p1
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;//a3=p3-p1

    // area of the parralelogram since h=1 and V=ah=nf.(a2xa3)
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      // the distance of the point from the triangle plane
      // d = Vol/area = a1.(a2xa3)/area
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/

	  // Calculate the interpolatin Coefficients
	  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
	      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
	      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    *Need3rdPoint = 0;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp_fsi(pjp, pj1, pj2, pj3, ibminfo,number);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    *Need3rdPoint = 1;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intp(pjp, pj2, pj3, ibminfo, number,1);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intp(pjp, pj2, pj3, ibminfo, number,1);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intp(pjp, pj2, pj3, ibminfo,number,1);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    *Need3rdPoint = 1;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intp(pjp, pj1, pj3, ibminfo, number,2);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intp(pjp, pj1, pj3, ibminfo, number,2);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intp(pjp, pj1, pj3, ibminfo,number,2);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }

	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    *Need3rdPoint = 1;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intp(pjp, pj1, pj2, ibminfo, number,3);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intp(pjp, pj1, pj2, ibminfo, number,3);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intp(pjp, pj1, pj2, ibminfo,number,3);
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    *Need3rdPoint = 1;
	    ibminfo[number].cr1 = 1.;
	    ibminfo[number].cr2 = 0.;
	    ibminfo[number].cr3 = 0.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    *Need3rdPoint = 1;
	    ibminfo[number].cr1 = 0.;
	    ibminfo[number].cr2 = 1.;
	    ibminfo[number].cr3 = 0.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    *Need3rdPoint = 1;
	    ibminfo[number].cr1 = 0.;
	    ibminfo[number].cr2 = 0.;
	    ibminfo[number].cr3 = 1.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else {
	    PetscPrintf(PETSC_COMM_WORLD, "%Something Wrong! All host nodes are blanked!!!!!\n");
	    return(1);
	  }

	  *intp = pint;

	  ibminfo[number].d_i = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo[number].imode = cell;

	  if (ibminfo[number].cr1<1e-6 &&
	      ibminfo[number].cr2<1e-6 &&
	      ibminfo[number].cr3<1e-6)
	    PetscPrintf(PETSC_COMM_SELF, "0 fsi Coeff!!!! %d  %d %le %le %le %le \n", number, ibminfo[number].imode,ibminfo[number].d_i, nvertpc[triangles[0][i]],nvertpc[triangles[1][i]],nvertpc[triangles[2][i]]);

	  return (0);
	}
      }
    }
  }
  return(0);
}

PetscErrorCode Find_fsi_interp_Coeff(IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo)
{

/* Note:  ibminfo returns the interpolation info 
   for the fsi (fsi_intp) */

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	i, j, k;
  PetscInt	i2, j2, k2;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt, ip[8],jp[8],kp[8];
  Cmpnts        ***coor,pc[8],p, intp;
  PetscReal	***nvert,nvertpc[8];
  PetscReal     nfx,nfy,nfz;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  DAVecGetArray(da, user->lNvert, &nvert);

  for (elmt=0; elmt<n_elmt; elmt++) {
  if (elmtinfo[elmt].n_P>0 && elmtinfo[elmt].FoundAroundcell>0) {
    //p=ibm->cent[elmt];
    p.x=ibm->cent_x[elmt]; 
    p.y=ibm->cent_y[elmt]; 
    p.z=ibm->cent_z[elmt]; 

    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    i=elmtinfo[elmt].icell;
    j=elmtinfo[elmt].jcell;
    k=elmtinfo[elmt].kcell;

    // normal correction for near domain bndry pts
/*     if (i==1) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i+1].x); */
/*     } */
    
/*     if (i==mx-3) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i-1].x); */
/*     } */
    
/*     if (j==1) { */
/*       nfy=-0.0001*PetscSign(coor[k][j][i].y-coor[k][j+1][i].y); */
/*     } */
    
/*     if (j==my-2) { */
/*       nfy=-0.0001*PetscSign(coor[k][j][i].y-coor[k][j-1][i].y); */
/*     } */

/*     if (k==1) { */
/*       nfz=-0.0001*PetscSign(coor[k][j][i].z-coor[k+1][j][i].z); */
/*     } */

/*     if (k==mz-2) { */
/*       nfz=-0.0001*PetscSign(coor[k][j][i].z-coor[k-1][j][i].z); */
/*     } */

    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {

      pc[0] = coor[k  ][j  ][i  ];
      pc[1] = coor[k  ][j  ][i+1];
      pc[2] = coor[k  ][j+1][i+1];
      pc[3] = coor[k  ][j+1][i  ];
      
      pc[4] = coor[k+1][j  ][i  ];
      pc[5] = coor[k+1][j  ][i+1];
      pc[6] = coor[k+1][j+1][i+1];
      pc[7] = coor[k+1][j+1][i  ];

      nvertpc[0] = nvert[k  ][j  ][i  ];
      nvertpc[1] = nvert[k  ][j  ][i+1];
      nvertpc[2] = nvert[k  ][j+1][i+1];
      nvertpc[3] = nvert[k  ][j+1][i  ];
      
      nvertpc[4] = nvert[k+1][j  ][i  ];
      nvertpc[5] = nvert[k+1][j  ][i+1];
      nvertpc[6] = nvert[k+1][j+1][i+1];
      nvertpc[7] = nvert[k+1][j+1][i  ];

      kp[0]=k  ;jp[0]=j  ;ip[0]=i  ;
      kp[1]=k  ;jp[1]=j  ;ip[1]=i+1;
      kp[2]=k  ;jp[2]=j+1;ip[2]=i+1;
      kp[3]=k  ;jp[3]=j+1;ip[3]=i  ;
      kp[4]=k+1;jp[4]=j  ;ip[4]=i  ;
      kp[5]=k+1;jp[5]=j  ;ip[5]=i+1;
      kp[6]=k+1;jp[6]=j+1;ip[6]=i+1;
      kp[7]=k+1;jp[7]=j+1;ip[7]=i  ;

      fsi_InterceptionPoint(p,pc,nvertpc, nfx, nfy,
	      nfz, ibminfo, elmt, &intp, 
	      &(elmtinfo[elmt].Need3rdPoint));

      switch (ibminfo[elmt].imode) {
      case(0): {
	ibminfo[elmt].i1=ip[0]; ibminfo[elmt].j1 = jp[0]; ibminfo[elmt].k1 = kp[0];
	ibminfo[elmt].i2=ip[1]; ibminfo[elmt].j2 = jp[1]; ibminfo[elmt].k2 = kp[1];
	ibminfo[elmt].i3=ip[2]; ibminfo[elmt].j3 = jp[2]; ibminfo[elmt].k3 = kp[2];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (1): {
	ibminfo[elmt].i1=ip[0]; ibminfo[elmt].j1 = jp[0]; ibminfo[elmt].k1 = kp[0];
	ibminfo[elmt].i2=ip[2]; ibminfo[elmt].j2 = jp[2]; ibminfo[elmt].k2 = kp[2];
	ibminfo[elmt].i3=ip[3]; ibminfo[elmt].j3 = jp[3]; ibminfo[elmt].k3 = kp[3];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j,k-1,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j,k-1,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (2): {
	ibminfo[elmt].i1=ip[4]; ibminfo[elmt].j1 = jp[4]; ibminfo[elmt].k1 = kp[4];
	ibminfo[elmt].i2=ip[5]; ibminfo[elmt].j2 = jp[5]; ibminfo[elmt].k2 = kp[5];
	ibminfo[elmt].i3=ip[6]; ibminfo[elmt].j3 = jp[6]; ibminfo[elmt].k3 = kp[6];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j,k+1,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j,k+1,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (3): {
	ibminfo[elmt].i1=ip[4]; ibminfo[elmt].j1 = jp[4]; ibminfo[elmt].k1 = kp[4];
	ibminfo[elmt].i2=ip[6]; ibminfo[elmt].j2 = jp[6]; ibminfo[elmt].k2 = kp[6];
	ibminfo[elmt].i3=ip[7]; ibminfo[elmt].j3 = jp[7]; ibminfo[elmt].k3 = kp[7];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j,k+1,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j,k+1,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (4): {
	ibminfo[elmt].i1=ip[0]; ibminfo[elmt].j1 = jp[0]; ibminfo[elmt].k1 = kp[0];
	ibminfo[elmt].i2=ip[4]; ibminfo[elmt].j2 = jp[4]; ibminfo[elmt].k2 = kp[4];
	ibminfo[elmt].i3=ip[3]; ibminfo[elmt].j3 = jp[3]; ibminfo[elmt].k3 = kp[3];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i-1,j,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i-1,j,k,elmt,intp,ibminfo,user,ibm); */       
      }
      case (5): {
	ibminfo[elmt].i1=ip[3]; ibminfo[elmt].j1 = jp[3]; ibminfo[elmt].k1 = kp[3];
	ibminfo[elmt].i2=ip[4]; ibminfo[elmt].j2 = jp[4]; ibminfo[elmt].k2 = kp[4];
	ibminfo[elmt].i3=ip[7]; ibminfo[elmt].j3 = jp[7]; ibminfo[elmt].k3 = kp[7];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i-1,j,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i-1,j,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (6): {
	ibminfo[elmt].i1=ip[1]; ibminfo[elmt].j1 = jp[1]; ibminfo[elmt].k1 = kp[1];
	ibminfo[elmt].i2=ip[5]; ibminfo[elmt].j2 = jp[5]; ibminfo[elmt].k2 = kp[5];
	ibminfo[elmt].i3=ip[2]; ibminfo[elmt].j3 = jp[2]; ibminfo[elmt].k3 = kp[2];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i+1,j,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i+1,j,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (7): {
	ibminfo[elmt].i1=ip[2]; ibminfo[elmt].j1 = jp[2]; ibminfo[elmt].k1 = kp[2];
	ibminfo[elmt].i2=ip[6]; ibminfo[elmt].j2 = jp[6]; ibminfo[elmt].k2 = kp[6];
	ibminfo[elmt].i3=ip[5]; ibminfo[elmt].j3 = jp[5]; ibminfo[elmt].k3 = kp[5];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i+1,j,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i+1,j,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (8): {
	ibminfo[elmt].i1=ip[0]; ibminfo[elmt].j1 = jp[0]; ibminfo[elmt].k1 = kp[0];
	ibminfo[elmt].i2=ip[4]; ibminfo[elmt].j2 = jp[4]; ibminfo[elmt].k2 = kp[4];
	ibminfo[elmt].i3=ip[1]; ibminfo[elmt].j3 = jp[1]; ibminfo[elmt].k3 = kp[1];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j-1,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j-1,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (9): {
	ibminfo[elmt].i1=ip[1]; ibminfo[elmt].j1 = jp[1]; ibminfo[elmt].k1 = kp[1];
	ibminfo[elmt].i2=ip[4]; ibminfo[elmt].j2 = jp[4]; ibminfo[elmt].k2 = kp[4];
	ibminfo[elmt].i3=ip[5]; ibminfo[elmt].j3 = jp[5]; ibminfo[elmt].k3 = kp[5];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j-1,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j-1,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (10): {
	ibminfo[elmt].i1=ip[3]; ibminfo[elmt].j1 = jp[3]; ibminfo[elmt].k1 = kp[3];
	ibminfo[elmt].i2=ip[7]; ibminfo[elmt].j2 = jp[7]; ibminfo[elmt].k2 = kp[7];
	ibminfo[elmt].i3=ip[2]; ibminfo[elmt].j3 = jp[2]; ibminfo[elmt].k3 = kp[2];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j+1,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j+1,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      case (11): {
	ibminfo[elmt].i1=ip[2]; ibminfo[elmt].j1 = jp[2]; ibminfo[elmt].k1 = kp[2];
	ibminfo[elmt].i2=ip[7]; ibminfo[elmt].j2 = jp[7]; ibminfo[elmt].k2 = kp[7];
	ibminfo[elmt].i3=ip[6]; ibminfo[elmt].j3 = jp[6]; ibminfo[elmt].k3 = kp[6];
	GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2);
	Find_fsi_2nd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm);
/* 	Find_fsi_2nd_interp_Coeff2(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff(i,j+1,k,elmt,intp,ibminfo,user,ibm); */
/* 	Find_fsi_2nd_interp_Coeff2(i,j+1,k,elmt,intp,ibminfo,user,ibm); */
	break;
      }
      }
      if (ibminfo[elmt].imode<0) 
	PetscPrintf(PETSC_COMM_SELF, "FSI Interpolation Coeffients Were not Found!!!! %d  %d %le \n", elmt, ibminfo[elmt].imode,ibminfo[elmt].d_i);
      //PetscPrintf(PETSC_COMM_SELF, "FSI Interpolatoion host %d  %d  %d \n", ibminfo[elmt].i1,ibminfo[elmt].j1,ibminfo[elmt].k1);

    }
  }
  }
  DAVecRestoreArray(fda, user->lCent,&coor);  
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode linear_intpp(Cpt2D p, Cpt2D p1, Cpt2D p2, 
			   IBMInfo *ibminfo, PetscInt number,
			   PetscInt nvert)
{
  PetscReal  x12, y12, xp2, yp2, Cr;
  x12 = p1.x - p2.x; y12 = p1.y - p2.y;
  xp2 = p.x - p2.x; yp2 = p.y - p2.y;

  if (fabs(x12)>1e-7) {
    Cr=xp2/x12;
  }  else if (fabs(y12)>1e-7) {
    Cr=yp2/y12;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "%Something Wrong!!! Linear intp two points are the same!!!\n");
  }
  
  if (nvert==1) {
    ibminfo[number].cr11 = 0.;
    ibminfo[number].cr22 = Cr;    
    ibminfo[number].cr33 = 1-Cr;
  } else if (nvert==2) {
    ibminfo[number].cr11 = Cr;
    ibminfo[number].cr22 = 0.;    
    ibminfo[number].cr33 = 1-Cr;
  } else if (nvert==3) {
    ibminfo[number].cr11 = Cr;
    ibminfo[number].cr22 = 1-Cr;    
    ibminfo[number].cr33 = 0.;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "%Wrong Nvert in Linear intp!!!\n");
  }
  return(0);  
}

PetscErrorCode triangle_intpp_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].cr11 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].cr22 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].cr33 = 1. - ibminfo[number].cr11 - ibminfo[number].cr22;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].cr33<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Crii %le %le %le\n", ibminfo[number].cr33, ibminfo[number].cr22, ibminfo[number].cr11);
  
}

PetscErrorCode fsi_InterceptionPoint2(Cmpnts p, Cmpnts pc[8], 
	       PetscReal nvertpc[8], PetscReal nfx, PetscReal nfy, PetscReal nfz, 
	       IBMInfo *ibminfo, PetscInt number)
{
  PetscInt 	triangles[3][12];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;//, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo[number].iimode = -100;
  // k-plane
  triangles[0][0]  = 0; triangles[1][0]  = 1; triangles[2][0]  = 3;
  triangles[0][1]  = 1; triangles[1][1]  = 2; triangles[2][1]  = 3;
  triangles[0][2]  = 4; triangles[1][2]  = 5; triangles[2][2]  = 7;
  triangles[0][3]  = 5; triangles[1][3]  = 6; triangles[2][3]  = 7;
  // i-plane
  triangles[0][4]  = 0; triangles[1][4]  = 7; triangles[2][4]  = 3;
  triangles[0][5]  = 0; triangles[1][5]  = 4; triangles[2][5]  = 7;
  triangles[0][6]  = 1; triangles[1][6]  = 5; triangles[2][6]  = 6;
  triangles[0][7]  = 2; triangles[1][7]  = 6; triangles[2][7]  = 1;
  // j-plane
  triangles[0][8]  = 0; triangles[1][8]  = 5; triangles[2][8]  = 1;
  triangles[0][9]  = 0; triangles[1][9]  = 4; triangles[2][9]  = 5;
  triangles[0][10] = 3; triangles[1][10] = 7; triangles[2][10] = 6;
  triangles[0][11] = 2; triangles[1][11] = 3; triangles[2][11] = 6;

  for (i=0; i<12; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;   //a1=p -p1
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;//a2=p2-p1
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;//a3=p3-p1

    // area of the parralelogram since h=1 and V=ah=nf.(a2xa3)
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      // the distance of the point from the triangle plane
      // d = Vol/area = a1.(a2xa3)/area
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/

	  // Calculate the interpolatin Coefficients
	  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
	      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
	      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {

	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intpp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intpp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intpp_fsi(pjp, pj1, pj2, pj3, ibminfo,number);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intpp(pjp, pj2, pj3, ibminfo, number,1);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intpp(pjp, pj2, pj3, ibminfo, number,1);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intpp(pjp, pj2, pj3, ibminfo,number,1);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intpp(pjp, pj1, pj3, ibminfo, number,2);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intpp(pjp, pj1, pj3, ibminfo, number,2);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intpp(pjp, pj1, pj3, ibminfo,number,2);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }

	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      linear_intpp(pjp, pj1, pj2, ibminfo, number,3);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      linear_intpp(pjp, pj1, pj2, ibminfo, number,3);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      linear_intpp(pjp, pj1, pj2, ibminfo,number,3);
	      triangle_intpp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    ibminfo[number].cr11 = 1.;
	    ibminfo[number].cr22 = 0.;
	    ibminfo[number].cr33 = 0.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) != 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) == 3) {
	    ibminfo[number].cr11 = 0.;
	    ibminfo[number].cr22 = 1.;
	    ibminfo[number].cr33 = 0.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else  if ((int)(nvertpc[triangles[0][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[1][i]] +0.5) == 3 &&
		      (int)(nvertpc[triangles[2][i]] +0.5) != 3) {
	    ibminfo[number].cr11 = 0.;
	    ibminfo[number].cr22 = 0.;
	    ibminfo[number].cr33 = 1.;
	    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	      pjp.x = pint.y; pjp.y = pint.z;
	      pj1.x = p1.y;   pj1.y = p1.z;
	      pj2.x = p2.y;   pj2.y = p2.z;
	      pj3.x = p3.y;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	      pjp.x = pint.x; pjp.y = pint.z;
	      pj1.x = p1.x;   pj1.y = p1.z;
	      pj2.x = p2.x;   pj2.y = p2.z;
	      pj3.x = p3.x;   pj3.y = p3.z;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	      pjp.x = pint.y; pjp.y = pint.x;
	      pj1.x = p1.y;   pj1.y = p1.x;
	      pj2.x = p2.y;   pj2.y = p2.x;
	      pj3.x = p3.y;   pj3.y = p3.x;
	      triangle_intp2_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	    }
	  } else {
	    PetscPrintf(PETSC_COMM_WORLD, "%Something Wrong! All host nodes are blanked!!!!!\n");
	    return(1);
	  }

	  ibminfo[number].d_ii = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo[number].iimode = cell;

	  if (ibminfo[number].cr11<1e-6 &&
	      ibminfo[number].cr22<1e-6 &&
	      ibminfo[number].cr33<1e-6)
	    PetscPrintf(PETSC_COMM_SELF, "0 fsi Coeff!!!! %d  %d %le %le %le %le \n", number, ibminfo[number].iimode,ibminfo[number].d_ii, nvertpc[triangles[0][i]],nvertpc[triangles[1][i]],nvertpc[triangles[2][i]]);

	  return (0);
	}
      }
    }
  }
  return(0);
}

PetscErrorCode Find_fsi_interp_Coeff2(IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo)
{

/* Note:  ibminfo returns the interpolation info 
   for the fsi (fsi_intp) */

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	i, j, k;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt, ip[8],jp[8],kp[8];
  Cmpnts        ***coor,pc[8],p;
  PetscReal	***nvert,nvertpc[8];
  PetscReal     nfx,nfy,nfz;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  DAVecGetArray(da, user->lNvert, &nvert);

  for (elmt=0; elmt<n_elmt; elmt++) {
  if (elmtinfo[elmt].n_P>0 && elmtinfo[elmt].FoundAroundcell>0) {
    //p=ibm->cent[elmt];
    p.x=ibm->cent_x[elmt]; 
    p.y=ibm->cent_y[elmt]; 
    p.z=ibm->cent_z[elmt]; 

    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    i=elmtinfo[elmt].icell;
    j=elmtinfo[elmt].jcell;
    k=elmtinfo[elmt].kcell;

    // normal correction for near domain bndry pts
/*     if (i==1) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i+1].x); */
/*     } */
    
/*     if (i==mx-3) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i-1].x); */
/*     } */

/*     if (j==1) { */
/*       nfy=-0.0001*PetscSign(coor[k][j][i].y-coor[k][j+1][i].y); */
/*     } */
    
/*     if (j==my-2) { */
/*       nfy=-0.0001*PetscSign(coor[k][j][i].y-coor[k][j-1][i].y); */
/*     } */

/*     if (k==1) { */
/*       nfz=-0.0001*PetscSign(coor[k][j][i].z-coor[k+1][j][i].z); */
/*     } */

/*     if (k==mz-2) { */
/*       nfz=-0.0001*PetscSign(coor[k][j][i].z-coor[k-1][j][i].z); */
/*     } */

    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {

      pc[0] = coor[k  ][j  ][i  ];
      pc[1] = coor[k  ][j  ][i+1];
      pc[2] = coor[k  ][j+1][i+1];
      pc[3] = coor[k  ][j+1][i  ];
      
      pc[4] = coor[k+1][j  ][i  ];
      pc[5] = coor[k+1][j  ][i+1];
      pc[6] = coor[k+1][j+1][i+1];
      pc[7] = coor[k+1][j+1][i  ];

      nvertpc[0] = nvert[k  ][j  ][i  ];
      nvertpc[1] = nvert[k  ][j  ][i+1];
      nvertpc[2] = nvert[k  ][j+1][i+1];
      nvertpc[3] = nvert[k  ][j+1][i  ];
      
      nvertpc[4] = nvert[k+1][j  ][i  ];
      nvertpc[5] = nvert[k+1][j  ][i+1];
      nvertpc[6] = nvert[k+1][j+1][i+1];
      nvertpc[7] = nvert[k+1][j+1][i  ];

      kp[0]=k  ;jp[0]=j  ;ip[0]=i  ;
      kp[1]=k  ;jp[1]=j  ;ip[1]=i+1;
      kp[2]=k  ;jp[2]=j+1;ip[2]=i+1;
      kp[3]=k  ;jp[3]=j+1;ip[3]=i  ;
      kp[4]=k+1;jp[4]=j  ;ip[4]=i  ;
      kp[5]=k+1;jp[5]=j  ;ip[5]=i+1;
      kp[6]=k+1;jp[6]=j+1;ip[6]=i+1;
      kp[7]=k+1;jp[7]=j+1;ip[7]=i  ;

      fsi_InterceptionPoint2(p,pc,nvertpc, nfx, nfy,
	      nfz, ibminfo, elmt);

      switch (ibminfo[elmt].iimode) {
      case(0): {
	ibminfo[elmt].i11=ip[0]; ibminfo[elmt].j11 = jp[0]; ibminfo[elmt].k11 = kp[0];
	ibminfo[elmt].i22=ip[1]; ibminfo[elmt].j22 = jp[1]; ibminfo[elmt].k22 = kp[1];
	ibminfo[elmt].i33=ip[3]; ibminfo[elmt].j33 = jp[3]; ibminfo[elmt].k33 = kp[3];
	break;
      }
      case (1): {
	ibminfo[elmt].i11=ip[1]; ibminfo[elmt].j11 = jp[1]; ibminfo[elmt].k11 = kp[1];
	ibminfo[elmt].i22=ip[2]; ibminfo[elmt].j22 = jp[2]; ibminfo[elmt].k22 = kp[2];
	ibminfo[elmt].i33=ip[3]; ibminfo[elmt].j33 = jp[3]; ibminfo[elmt].k33 = kp[3];
	break;
      }
      case (2): {
	ibminfo[elmt].i11=ip[4]; ibminfo[elmt].j11 = jp[4]; ibminfo[elmt].k11 = kp[4];
	ibminfo[elmt].i22=ip[5]; ibminfo[elmt].j22 = jp[5]; ibminfo[elmt].k22 = kp[5];
	ibminfo[elmt].i33=ip[7]; ibminfo[elmt].j33 = jp[7]; ibminfo[elmt].k33 = kp[7];
	break;
      }
      case (3): {
	ibminfo[elmt].i11=ip[5]; ibminfo[elmt].j11 = jp[5]; ibminfo[elmt].k11 = kp[5];
	ibminfo[elmt].i22=ip[6]; ibminfo[elmt].j22 = jp[6]; ibminfo[elmt].k22 = kp[6];
	ibminfo[elmt].i33=ip[7]; ibminfo[elmt].j33 = jp[7]; ibminfo[elmt].k33 = kp[7];
	break;
      }
      case (4): {
	ibminfo[elmt].i11=ip[0]; ibminfo[elmt].j11 = jp[0]; ibminfo[elmt].k11 = kp[0];
	ibminfo[elmt].i22=ip[7]; ibminfo[elmt].j22 = jp[7]; ibminfo[elmt].k22 = kp[7];
	ibminfo[elmt].i33=ip[3]; ibminfo[elmt].j33 = jp[3]; ibminfo[elmt].k33 = kp[3];
	break;
      }
      case (5): {
	ibminfo[elmt].i11=ip[0]; ibminfo[elmt].j11 = jp[0]; ibminfo[elmt].k11 = kp[0];
	ibminfo[elmt].i22=ip[4]; ibminfo[elmt].j22 = jp[4]; ibminfo[elmt].k22 = kp[4];
	ibminfo[elmt].i33=ip[7]; ibminfo[elmt].j33 = jp[7]; ibminfo[elmt].k33 = kp[7];
	break;
      }
      case (6): {
	ibminfo[elmt].i11=ip[1]; ibminfo[elmt].j11 = jp[1]; ibminfo[elmt].k11 = kp[1];
	ibminfo[elmt].i22=ip[5]; ibminfo[elmt].j22 = jp[5]; ibminfo[elmt].k22 = kp[5];
	ibminfo[elmt].i33=ip[6]; ibminfo[elmt].j33 = jp[6]; ibminfo[elmt].k33 = kp[6];
	break;
      }
      case (7): {
	ibminfo[elmt].i11=ip[2]; ibminfo[elmt].j11 = jp[2]; ibminfo[elmt].k11 = kp[2];
	ibminfo[elmt].i22=ip[6]; ibminfo[elmt].j22 = jp[6]; ibminfo[elmt].k22 = kp[6];
	ibminfo[elmt].i33=ip[1]; ibminfo[elmt].j33 = jp[1]; ibminfo[elmt].k33 = kp[1];
	break;
      }
      case (8): {
	ibminfo[elmt].i11=ip[0]; ibminfo[elmt].j11 = jp[0]; ibminfo[elmt].k11 = kp[0];
	ibminfo[elmt].i22=ip[5]; ibminfo[elmt].j22 = jp[5]; ibminfo[elmt].k22 = kp[5];
	ibminfo[elmt].i33=ip[1]; ibminfo[elmt].j33 = jp[1]; ibminfo[elmt].k33 = kp[1];
	break;
      }
      case (9): {
	ibminfo[elmt].i11=ip[0]; ibminfo[elmt].j11 = jp[0]; ibminfo[elmt].k11 = kp[0];
	ibminfo[elmt].i22=ip[4]; ibminfo[elmt].j22 = jp[4]; ibminfo[elmt].k22 = kp[4];
	ibminfo[elmt].i33=ip[5]; ibminfo[elmt].j33 = jp[5]; ibminfo[elmt].k33 = kp[5];
	break;
      }
      case (10): {
	ibminfo[elmt].i11=ip[3]; ibminfo[elmt].j11 = jp[3]; ibminfo[elmt].k11 = kp[3];
	ibminfo[elmt].i22=ip[7]; ibminfo[elmt].j22 = jp[7]; ibminfo[elmt].k22 = kp[7];
	ibminfo[elmt].i33=ip[6]; ibminfo[elmt].j33 = jp[6]; ibminfo[elmt].k33 = kp[6];
	break;
      }
      case (11): {
	ibminfo[elmt].i11=ip[2]; ibminfo[elmt].j11 = jp[2]; ibminfo[elmt].k11 = kp[2];
	ibminfo[elmt].i22=ip[3]; ibminfo[elmt].j22 = jp[3]; ibminfo[elmt].k22 = kp[3];
	ibminfo[elmt].i33=ip[6]; ibminfo[elmt].j33 = jp[6]; ibminfo[elmt].k33 = kp[6];
	break;
      }
      }
      if (ibminfo[elmt].iimode<0) 
	PetscPrintf(PETSC_COMM_SELF, "FSI Interpolation Coeffients Were not Found!!!! %d  %d %le \n", elmt, ibminfo[elmt].imode,ibminfo[elmt].d_i);
      //PetscPrintf(PETSC_COMM_SELF, "FSI Interpolatoion host %d  %d  %d \n", ibminfo[elmt].i1,ibminfo[elmt].j1,ibminfo[elmt].k1);

    }
  }
  }
  DAVecRestoreArray(fda, user->lCent,&coor);  
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode triangle_2nd_intp_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].ct1 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].ct2 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].ct3 = 1. - ibminfo[number].ct1 - ibminfo[number].ct2;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].ct3<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Ct %d %le %le %le\n",number,ibminfo[number].ct3, ibminfo[number].ct2, ibminfo[number].ct1);
  return(0);
}

PetscErrorCode fsi_2nd_InterceptionPoint(Cmpnts p, Cmpnts pc[8], 
	       PetscReal nfx, PetscReal nfy, PetscReal nfz, 
	       IBMInfo *ibminfo, PetscInt number, Cmpnts *intp)
{
  PetscInt 	triangles[3][12];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;//, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo[number].smode = -100;
  // k-plane
  triangles[0][0]  = 0; triangles[1][0]  = 1; triangles[2][0]  = 2;
  triangles[0][1]  = 0; triangles[1][1]  = 2; triangles[2][1]  = 3;
  triangles[0][2]  = 4; triangles[1][2]  = 5; triangles[2][2]  = 6;
  triangles[0][3]  = 4; triangles[1][3]  = 6; triangles[2][3]  = 7;
  // i-plane
  triangles[0][4]  = 0; triangles[1][4]  = 4; triangles[2][4]  = 3;
  triangles[0][5]  = 3; triangles[1][5]  = 4; triangles[2][5]  = 7;
  triangles[0][6]  = 1; triangles[1][6]  = 5; triangles[2][6]  = 2;
  triangles[0][7]  = 2; triangles[1][7]  = 6; triangles[2][7]  = 5;
  // j-plane
  triangles[0][8]  = 0; triangles[1][8]  = 4; triangles[2][8]  = 1;
  triangles[0][9]  = 1; triangles[1][9]  = 4; triangles[2][9]  = 5;
  triangles[0][10] = 3; triangles[1][10] = 7; triangles[2][10] = 2;
  triangles[0][11] = 2; triangles[1][11] = 7; triangles[2][11] = 6;

  for (i=0; i<12; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;   //a1=p -p1
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;//a2=p2-p1
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;//a3=p3-p1

    // area of the parralelogram since h=1 and V=ah=nf.(a2xa3)
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      // the distance of the point from the triangle plane
      // d = Vol/area = a1.(a2xa3)/area
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/

	  // Calculate the interpolatin Coefficients
	  
	  if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	    pjp.x = pint.y; pjp.y = pint.z;
	    pj1.x = p1.y;   pj1.y = p1.z;
	    pj2.x = p2.y;   pj2.y = p2.z;
	    pj3.x = p3.y;   pj3.y = p3.z;
	    triangle_2nd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	    pjp.x = pint.x; pjp.y = pint.z;
	    pj1.x = p1.x;   pj1.y = p1.z;
	    pj2.x = p2.x;   pj2.y = p2.z;
	    pj3.x = p3.x;   pj3.y = p3.z;
	    triangle_2nd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	    pjp.x = pint.y; pjp.y = pint.x;
	    pj1.x = p1.y;   pj1.y = p1.x;
	    pj2.x = p2.y;   pj2.y = p2.x;
	    pj3.x = p3.y;   pj3.y = p3.x;
	    triangle_2nd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo,number);
	  }
	  
	  ibminfo[number].d_s = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo[number].smode = cell;

	  *intp = pint;

	  if (ibminfo[number].ct1<1e-6 &&
	      ibminfo[number].ct2<1e-6 &&
	      ibminfo[number].ct3<1e-6)
	    PetscPrintf(PETSC_COMM_SELF, "0 fsi Coeff 2nd fsi!!!! %d  %d %le %le %le %le \n", number, ibminfo[number].smode,ibminfo[number].d_s);
	  
	  return (0);
	}
      }
    }
  }
  return(0);
}

PetscErrorCode Find_fsi_2nd_interp_Coeff(PetscInt i, PetscInt j,
					 PetscInt k, PetscInt elmt,
					 Cmpnts   p,
					 IBMInfo *ibminfo,
					 UserCtx *user, 
					 IBMNodes *ibm)
{

/* Note:  ibminfo returns the interpolation info 
   for the fsi (fsi_intp) */

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  //PetscInt	i, j, k;
  PetscInt	i2, j2, k2;

  //PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      ip[8],jp[8],kp[8];
  Cmpnts        ***coor,pc[8];
  Cmpnts        intp;
  //PetscReal	***nvert,nvertpc[8];
  PetscReal     nfx,nfy,nfz;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  //DAVecGetArray(da, user->lNvert, &nvert);

  nfx=ibm->nf_x[elmt];
  nfy=ibm->nf_y[elmt];
  nfz=ibm->nf_z[elmt];

    // normal correction for near domain bndry pts
/*     if (i==1) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i+1].x); */
/*     } */
    
/*     if (i==mx-3) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i-1].x); */
/*     } */


  // normal correction for near domain bndry pts
  
  if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
    
    pc[0] = coor[k  ][j  ][i  ];
    pc[1] = coor[k  ][j  ][i+1];
    pc[2] = coor[k  ][j+1][i+1];
    pc[3] = coor[k  ][j+1][i  ];
    
    pc[4] = coor[k+1][j  ][i  ];
    pc[5] = coor[k+1][j  ][i+1];
    pc[6] = coor[k+1][j+1][i+1];
    pc[7] = coor[k+1][j+1][i  ];
        
    kp[0]=k  ;jp[0]=j  ;ip[0]=i  ;
    kp[1]=k  ;jp[1]=j  ;ip[1]=i+1;
    kp[2]=k  ;jp[2]=j+1;ip[2]=i+1;
    kp[3]=k  ;jp[3]=j+1;ip[3]=i  ;
    kp[4]=k+1;jp[4]=j  ;ip[4]=i  ;
    kp[5]=k+1;jp[5]=j  ;ip[5]=i+1;
    kp[6]=k+1;jp[6]=j+1;ip[6]=i+1;
    kp[7]=k+1;jp[7]=j+1;ip[7]=i  ;
    
    fsi_2nd_InterceptionPoint(p,pc, nfx, nfy,
			      nfz, ibminfo,elmt, &intp);
    
    switch (ibminfo[elmt].smode) {
    case(0): {
      ibminfo[elmt].ii1=ip[0]; ibminfo[elmt].jj1 = jp[0]; ibminfo[elmt].kk1 = kp[0];
      ibminfo[elmt].ii2=ip[1]; ibminfo[elmt].jj2 = jp[1]; ibminfo[elmt].kk2 = kp[1];
      ibminfo[elmt].ii3=ip[2]; ibminfo[elmt].jj3 = jp[2]; ibminfo[elmt].kk3 = kp[2];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (1): {
      ibminfo[elmt].ii1=ip[0]; ibminfo[elmt].jj1 = jp[0]; ibminfo[elmt].kk1 = kp[0];
      ibminfo[elmt].ii2=ip[2]; ibminfo[elmt].jj2 = jp[2]; ibminfo[elmt].kk2 = kp[2];
      ibminfo[elmt].ii3=ip[3]; ibminfo[elmt].jj3 = jp[3]; ibminfo[elmt].kk3 = kp[3];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (2): {
      ibminfo[elmt].ii1=ip[4]; ibminfo[elmt].jj1 = jp[4]; ibminfo[elmt].kk1 = kp[4];
      ibminfo[elmt].ii2=ip[5]; ibminfo[elmt].jj2 = jp[5]; ibminfo[elmt].kk2 = kp[5];
      ibminfo[elmt].ii3=ip[6]; ibminfo[elmt].jj3 = jp[6]; ibminfo[elmt].kk3 = kp[6];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (3): {
      ibminfo[elmt].ii1=ip[4]; ibminfo[elmt].jj1 = jp[4]; ibminfo[elmt].kk1 = kp[4];
      ibminfo[elmt].ii2=ip[6]; ibminfo[elmt].jj2 = jp[6]; ibminfo[elmt].kk2 = kp[6];
      ibminfo[elmt].ii3=ip[7]; ibminfo[elmt].jj3 = jp[7]; ibminfo[elmt].kk3 = kp[7];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (4): {
      ibminfo[elmt].ii1=ip[0]; ibminfo[elmt].jj1 = jp[0]; ibminfo[elmt].kk1 = kp[0];
      ibminfo[elmt].ii2=ip[4]; ibminfo[elmt].jj2 = jp[4]; ibminfo[elmt].kk2 = kp[4];
      ibminfo[elmt].ii3=ip[3]; ibminfo[elmt].jj3 = jp[3]; ibminfo[elmt].kk3 = kp[3];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (5): {
      ibminfo[elmt].ii1=ip[3]; ibminfo[elmt].jj1 = jp[3]; ibminfo[elmt].kk1 = kp[3];
      ibminfo[elmt].ii2=ip[4]; ibminfo[elmt].jj2 = jp[4]; ibminfo[elmt].kk2 = kp[4];
      ibminfo[elmt].ii3=ip[7]; ibminfo[elmt].jj3 = jp[7]; ibminfo[elmt].kk3 = kp[7];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (6): {
      ibminfo[elmt].ii1=ip[1]; ibminfo[elmt].jj1 = jp[1]; ibminfo[elmt].kk1 = kp[1];
      ibminfo[elmt].ii2=ip[5]; ibminfo[elmt].jj2 = jp[5]; ibminfo[elmt].kk2 = kp[5];
      ibminfo[elmt].ii3=ip[2]; ibminfo[elmt].jj3 = jp[2]; ibminfo[elmt].kk3 = kp[2];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (7): {
      ibminfo[elmt].ii1=ip[2]; ibminfo[elmt].jj1 = jp[2]; ibminfo[elmt].kk1 = kp[2];
      ibminfo[elmt].ii2=ip[6]; ibminfo[elmt].jj2 = jp[6]; ibminfo[elmt].kk2 = kp[6];
      ibminfo[elmt].ii3=ip[5]; ibminfo[elmt].jj3 = jp[5]; ibminfo[elmt].kk3 = kp[5];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (8): {
      ibminfo[elmt].ii1=ip[0]; ibminfo[elmt].jj1 = jp[0]; ibminfo[elmt].kk1 = kp[0];
      ibminfo[elmt].ii2=ip[4]; ibminfo[elmt].jj2 = jp[4]; ibminfo[elmt].kk2 = kp[4];
      ibminfo[elmt].ii3=ip[1]; ibminfo[elmt].jj3 = jp[1]; ibminfo[elmt].kk3 = kp[1];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (9): {
      ibminfo[elmt].ii1=ip[1]; ibminfo[elmt].jj1 = jp[1]; ibminfo[elmt].kk1 = kp[1];
      ibminfo[elmt].ii2=ip[4]; ibminfo[elmt].jj2 = jp[4]; ibminfo[elmt].kk2 = kp[4];
      ibminfo[elmt].ii3=ip[5]; ibminfo[elmt].jj3 = jp[5]; ibminfo[elmt].kk3 = kp[5];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (10): {
      ibminfo[elmt].ii1=ip[3]; ibminfo[elmt].jj1 = jp[3]; ibminfo[elmt].kk1 = kp[3];
      ibminfo[elmt].ii2=ip[7]; ibminfo[elmt].jj2 = jp[7]; ibminfo[elmt].kk2 = kp[7];
      ibminfo[elmt].ii3=ip[2]; ibminfo[elmt].jj3 = jp[2]; ibminfo[elmt].kk3 = kp[2];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    case (11): {
      ibminfo[elmt].ii1=ip[2]; ibminfo[elmt].jj1 = jp[2]; ibminfo[elmt].kk1 = kp[2];
      ibminfo[elmt].ii2=ip[7]; ibminfo[elmt].jj2 = jp[7]; ibminfo[elmt].kk2 = kp[7];
      ibminfo[elmt].ii3=ip[6]; ibminfo[elmt].jj3 = jp[6]; ibminfo[elmt].kk3 = kp[6];
/*       GridCellaround2ndElmt(user,ibm,intp,elmt,k,j,i,&k2,&j2,&i2); */
/*       Find_fsi_3rd_interp_Coeff(i2,j2,k2,elmt,intp,ibminfo,user,ibm); */
      break;
    }
    }
    if (ibminfo[elmt].smode<0) 
      PetscPrintf(PETSC_COMM_SELF, "FSI 2nd Interpolation Coeffients Were not Found!!!! %d  %d %le %d %d %d\n", elmt, ibminfo[elmt].smode,ibminfo[elmt].d_s,i,j,k);
    //PetscPrintf(PETSC_COMM_SELF, "FSI Interpolatoion host %d  %d  %d \n", ibminfo[elmt].i1,ibminfo[elmt].j1,ibminfo[elmt].k1);
    
  }

  DAVecRestoreArray(fda, user->lCent,&coor);  
  //DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode triangle_3rd_intp_fsi(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].ct11 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].ct22 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].ct33 = 1. - ibminfo[number].ct11 - ibminfo[number].ct22;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  if (ibminfo[number].ct33<0.)
    PetscPrintf(PETSC_COMM_WORLD, "SOMETHING WRONG!!!! fsi_intp Ct ii %d %le %le %le\n",number,ibminfo[number].ct33, ibminfo[number].ct22, ibminfo[number].ct11);
  
}

PetscErrorCode fsi_3rd_InterceptionPoint(Cmpnts p, Cmpnts pc[8], 
	       PetscReal nfx, PetscReal nfy, PetscReal nfz, 
	       IBMInfo *ibminfo, PetscInt number, Cmpnts *intp)
{
  PetscInt 	triangles[3][12];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;//, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo[number].ssmode = -100;
  // k-plane
  triangles[0][0]  = 0; triangles[1][0]  = 1; triangles[2][0]  = 2;
  triangles[0][1]  = 0; triangles[1][1]  = 2; triangles[2][1]  = 3;
  triangles[0][2]  = 4; triangles[1][2]  = 5; triangles[2][2]  = 6;
  triangles[0][3]  = 4; triangles[1][3]  = 6; triangles[2][3]  = 7;
  // i-plane
  triangles[0][4]  = 0; triangles[1][4]  = 4; triangles[2][4]  = 3;
  triangles[0][5]  = 3; triangles[1][5]  = 4; triangles[2][5]  = 7;
  triangles[0][6]  = 1; triangles[1][6]  = 5; triangles[2][6]  = 2;
  triangles[0][7]  = 2; triangles[1][7]  = 6; triangles[2][7]  = 5;
  // j-plane
  triangles[0][8]  = 0; triangles[1][8]  = 4; triangles[2][8]  = 1;
  triangles[0][9]  = 1; triangles[1][9]  = 4; triangles[2][9]  = 5;
  triangles[0][10] = 3; triangles[1][10] = 7; triangles[2][10] = 2;
  triangles[0][11] = 2; triangles[1][11] = 7; triangles[2][11] = 6;

  for (i=0; i<12; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;   //a1=p -p1
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;//a2=p2-p1
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;//a3=p3-p1

    // area of the parralelogram since h=1 and V=ah=nf.(a2xa3)
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      // the distance of the point from the triangle plane
      // d = Vol/area = a1.(a2xa3)/area
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/

	  // Calculate the interpolatin Coefficients
	  
	  if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	    pjp.x = pint.y; pjp.y = pint.z;
	    pj1.x = p1.y;   pj1.y = p1.z;
	    pj2.x = p2.y;   pj2.y = p2.z;
	    pj3.x = p3.y;   pj3.y = p3.z;
	    triangle_3rd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	    pjp.x = pint.x; pjp.y = pint.z;
	    pj1.x = p1.x;   pj1.y = p1.z;
	    pj2.x = p2.x;   pj2.y = p2.z;
	    pj3.x = p3.x;   pj3.y = p3.z;
	    triangle_3rd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	    pjp.x = pint.y; pjp.y = pint.x;
	    pj1.x = p1.y;   pj1.y = p1.x;
	    pj2.x = p2.y;   pj2.y = p2.x;
	    pj3.x = p3.y;   pj3.y = p3.x;
	    triangle_3rd_intp_fsi(pjp, pj1, pj2, pj3, ibminfo,number);
	  }
	  
	  ibminfo[number].d_ss = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo[number].ssmode = cell;

	  *intp = pint;

	  if (ibminfo[number].ct11<1e-6 &&
	      ibminfo[number].ct22<1e-6 &&
	      ibminfo[number].ct33<1e-6)
	    PetscPrintf(PETSC_COMM_SELF, "0 fsi Coeff 3rd fsi!!!! %d  %d %le %le %le %le \n", number, ibminfo[number].smode,ibminfo[number].d_s);
	  
	  return (0);
	}
      }
    }
  }
  return(0);
}

PetscErrorCode Find_fsi_3rd_interp_Coeff(PetscInt i, PetscInt j,
					 PetscInt k, PetscInt elmt,
					 Cmpnts   p,
					 IBMInfo *ibminfo,
					 UserCtx *user, 
					 IBMNodes *ibm)
{

/* Note:  ibminfo returns the interpolation info 
   for the fsi (fsi_intp) */

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  //PetscInt	i, j, k;
  PetscInt	i2, j2, k2;

  //PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      ip[8],jp[8],kp[8];
  Cmpnts        ***coor,pc[8];
  Cmpnts        intp;
  //PetscReal	***nvert,nvertpc[8];
  PetscReal     nfx,nfy,nfz;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  //DAVecGetArray(da, user->lNvert, &nvert);

  nfx=ibm->nf_x[elmt];
  nfy=ibm->nf_y[elmt];
  nfz=ibm->nf_z[elmt];

    // normal correction for near domain bndry pts
/*     if (i==1) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i+1].x); */
/*     } */
    
/*     if (i==mx-3) { */
/*       nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i-1].x); */
/*     } */


  // normal correction for near domain bndry pts
  
  if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
    
    pc[0] = coor[k  ][j  ][i  ];
    pc[1] = coor[k  ][j  ][i+1];
    pc[2] = coor[k  ][j+1][i+1];
    pc[3] = coor[k  ][j+1][i  ];
    
    pc[4] = coor[k+1][j  ][i  ];
    pc[5] = coor[k+1][j  ][i+1];
    pc[6] = coor[k+1][j+1][i+1];
    pc[7] = coor[k+1][j+1][i  ];
        
    kp[0]=k  ;jp[0]=j  ;ip[0]=i  ;
    kp[1]=k  ;jp[1]=j  ;ip[1]=i+1;
    kp[2]=k  ;jp[2]=j+1;ip[2]=i+1;
    kp[3]=k  ;jp[3]=j+1;ip[3]=i  ;
    kp[4]=k+1;jp[4]=j  ;ip[4]=i  ;
    kp[5]=k+1;jp[5]=j  ;ip[5]=i+1;
    kp[6]=k+1;jp[6]=j+1;ip[6]=i+1;
    kp[7]=k+1;jp[7]=j+1;ip[7]=i  ;
    
    fsi_3rd_InterceptionPoint(p,pc, nfx, nfy,
			      nfz, ibminfo,elmt, &intp);
    
    switch (ibminfo[elmt].ssmode) {
    case(0): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[1]; ibminfo[elmt].jj22 = jp[1]; ibminfo[elmt].kk22 = kp[1];
      ibminfo[elmt].ii33=ip[2]; ibminfo[elmt].jj33 = jp[2]; ibminfo[elmt].kk33 = kp[2];
      break;
    }
    case (1): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[2]; ibminfo[elmt].jj22 = jp[2]; ibminfo[elmt].kk22 = kp[2];
      ibminfo[elmt].ii33=ip[3]; ibminfo[elmt].jj33 = jp[3]; ibminfo[elmt].kk33 = kp[3];
      break;
    }
    case (2): {
      ibminfo[elmt].ii11=ip[4]; ibminfo[elmt].jj11 = jp[4]; ibminfo[elmt].kk11 = kp[4];
      ibminfo[elmt].ii22=ip[5]; ibminfo[elmt].jj22 = jp[5]; ibminfo[elmt].kk22 = kp[5];
      ibminfo[elmt].ii33=ip[6]; ibminfo[elmt].jj33 = jp[6]; ibminfo[elmt].kk33 = kp[6];
      break;
    }
    case (3): {
      ibminfo[elmt].ii11=ip[4]; ibminfo[elmt].jj11 = jp[4]; ibminfo[elmt].kk11 = kp[4];
      ibminfo[elmt].ii22=ip[6]; ibminfo[elmt].jj22 = jp[6]; ibminfo[elmt].kk22 = kp[6];
      ibminfo[elmt].ii33=ip[7]; ibminfo[elmt].jj33 = jp[7]; ibminfo[elmt].kk33 = kp[7];
      break;
    }
    case (4): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[3]; ibminfo[elmt].jj33 = jp[3]; ibminfo[elmt].kk33 = kp[3];
      break;
    }
    case (5): {
      ibminfo[elmt].ii11=ip[3]; ibminfo[elmt].jj11 = jp[3]; ibminfo[elmt].kk11 = kp[3];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[7]; ibminfo[elmt].jj33 = jp[7]; ibminfo[elmt].kk33 = kp[7];
      break;
    }
    case (6): {
      ibminfo[elmt].ii11=ip[1]; ibminfo[elmt].jj11 = jp[1]; ibminfo[elmt].kk11 = kp[1];
      ibminfo[elmt].ii22=ip[5]; ibminfo[elmt].jj22 = jp[5]; ibminfo[elmt].kk22 = kp[5];
      ibminfo[elmt].ii33=ip[2]; ibminfo[elmt].jj33 = jp[2]; ibminfo[elmt].kk33 = kp[2];
      break;
    }
    case (7): {
      ibminfo[elmt].ii11=ip[2]; ibminfo[elmt].jj11 = jp[2]; ibminfo[elmt].kk11 = kp[2];
      ibminfo[elmt].ii22=ip[6]; ibminfo[elmt].jj22 = jp[6]; ibminfo[elmt].kk22 = kp[6];
      ibminfo[elmt].ii33=ip[5]; ibminfo[elmt].jj33 = jp[5]; ibminfo[elmt].kk33 = kp[5];
      break;
    }
    case (8): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[1]; ibminfo[elmt].jj33 = jp[1]; ibminfo[elmt].kk33 = kp[1];
      break;
    }
    case (9): {
      ibminfo[elmt].ii11=ip[1]; ibminfo[elmt].jj11 = jp[1]; ibminfo[elmt].kk11 = kp[1];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[5]; ibminfo[elmt].jj33 = jp[5]; ibminfo[elmt].kk33 = kp[5];
      break;
    }
    case (10): {
      ibminfo[elmt].ii11=ip[3]; ibminfo[elmt].jj11 = jp[3]; ibminfo[elmt].kk11 = kp[3];
      ibminfo[elmt].ii22=ip[7]; ibminfo[elmt].jj22 = jp[7]; ibminfo[elmt].kk22 = kp[7];
      ibminfo[elmt].ii33=ip[2]; ibminfo[elmt].jj33 = jp[2]; ibminfo[elmt].kk33 = kp[2];
      break;
    }
    case (11): {
      ibminfo[elmt].ii11=ip[2]; ibminfo[elmt].jj11 = jp[2]; ibminfo[elmt].kk11 = kp[2];
      ibminfo[elmt].ii22=ip[7]; ibminfo[elmt].jj22 = jp[7]; ibminfo[elmt].kk22 = kp[7];
      ibminfo[elmt].ii33=ip[6]; ibminfo[elmt].jj33 = jp[6]; ibminfo[elmt].kk33 = kp[6];
      break;
    }
    }
    if (ibminfo[elmt].ssmode<0) 
      PetscPrintf(PETSC_COMM_SELF, "FSI 2nd Interpolation Coeffients Were not Found!!!! %d  %d %le %d %d %d\n", elmt, ibminfo[elmt].smode,ibminfo[elmt].d_s,i,j,k);
    //PetscPrintf(PETSC_COMM_SELF, "FSI Interpolatoion host %d  %d  %d \n", ibminfo[elmt].i1,ibminfo[elmt].j1,ibminfo[elmt].k1);
    
  }

  DAVecRestoreArray(fda, user->lCent,&coor);  
  //DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode triangle_2nd_intp_fsi2(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo, PetscInt number)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo[number].ct11 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo[number].ct22 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo[number].ct33 = 1. - ibminfo[number].ct11 - ibminfo[number].ct22;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  return(0);
}

PetscErrorCode fsi_2nd_InterceptionPoint2(Cmpnts p, Cmpnts pc[8], 
	       PetscReal nfx, PetscReal nfy, PetscReal nfz, 
	       IBMInfo *ibminfo, PetscInt number)
{
  PetscInt 	triangles[3][12];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;//, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo[number].ssmode = -100;
  // k-plane
  triangles[0][0]  = 0; triangles[1][0]  = 1; triangles[2][0]  = 3;
  triangles[0][1]  = 1; triangles[1][1]  = 2; triangles[2][1]  = 3;
  triangles[0][2]  = 4; triangles[1][2]  = 5; triangles[2][2]  = 7;
  triangles[0][3]  = 5; triangles[1][3]  = 6; triangles[2][3]  = 7;
  // i-plane
  triangles[0][4]  = 0; triangles[1][4]  = 7; triangles[2][4]  = 3;
  triangles[0][5]  = 0; triangles[1][5]  = 4; triangles[2][5]  = 7;
  triangles[0][6]  = 1; triangles[1][6]  = 5; triangles[2][6]  = 6;
  triangles[0][7]  = 2; triangles[1][7]  = 6; triangles[2][7]  = 1;
  // j-plane
  triangles[0][8]  = 0; triangles[1][8]  = 5; triangles[2][8]  = 1;
  triangles[0][9]  = 0; triangles[1][9]  = 4; triangles[2][9]  = 5;
  triangles[0][10] = 3; triangles[1][10] = 7; triangles[2][10] = 6;
  triangles[0][11] = 2; triangles[1][11] = 3; triangles[2][11] = 6;

  for (i=0; i<12; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;   //a1=p -p1
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;//a2=p2-p1
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;//a3=p3-p1

    // area of the parralelogram since h=1 and V=ah=nf.(a2xa3)
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      // the distance of the point from the triangle plane
      // d = Vol/area = a1.(a2xa3)/area
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/

	  // Calculate the interpolatin Coefficients
	  
	  if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	    pjp.x = pint.y; pjp.y = pint.z;
	    pj1.x = p1.y;   pj1.y = p1.z;
	    pj2.x = p2.y;   pj2.y = p2.z;
	    pj3.x = p3.y;   pj3.y = p3.z;
	    triangle_2nd_intp_fsi2(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	    pjp.x = pint.x; pjp.y = pint.z;
	    pj1.x = p1.x;   pj1.y = p1.z;
	    pj2.x = p2.x;   pj2.y = p2.z;
	    pj3.x = p3.x;   pj3.y = p3.z;
	    triangle_2nd_intp_fsi2(pjp, pj1, pj2, pj3, ibminfo, number);
	  }
	  else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	    pjp.x = pint.y; pjp.y = pint.x;
	    pj1.x = p1.y;   pj1.y = p1.x;
	    pj2.x = p2.y;   pj2.y = p2.x;
	    pj3.x = p3.y;   pj3.y = p3.x;
	    triangle_2nd_intp_fsi2(pjp, pj1, pj2, pj3, ibminfo,number);
	  }
	  
	  ibminfo[number].d_ss = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo[number].ssmode = cell;
	  
	  if (ibminfo[number].ct11<1e-6 &&
	      ibminfo[number].ct22<1e-6 &&
	      ibminfo[number].ct33<1e-6)
	    PetscPrintf(PETSC_COMM_SELF, "0 fsi Coeff 2nd fsi!!!! %d  %d %le %le %le %le \n", number, ibminfo[number].imode,ibminfo[number].d_i);
	  
	  return (0);
	}
      }
    }
  }
  return(0);
}

PetscErrorCode Find_fsi_2nd_interp_Coeff2(PetscInt i, PetscInt j,
					 PetscInt k, PetscInt elmt,
					 Cmpnts   p,
					 IBMInfo *ibminfo,
					 UserCtx *user, 
					 IBMNodes *ibm)
{

/* Note:  ibminfo returns the interpolation info 
   for the fsi (fsi_intp) */

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  //PetscInt	i, j, k;

  //PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      ip[8],jp[8],kp[8];
  Cmpnts        ***coor,pc[8];
  //PetscReal	***nvert,nvertpc[8];
  PetscReal     nfx,nfy,nfz;

  PetscInt	rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  DAVecGetArray(fda, user->lCent,&coor);
  //DAVecGetArray(da, user->lNvert, &nvert);

  nfx=ibm->nf_x[elmt];
  nfy=ibm->nf_y[elmt];
  nfz=ibm->nf_z[elmt];

  // normal correction for near domain bndry pts
/*   if (i==1) { */
/*     nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i+1].x); */
/*   } */
  
/*   if (i==mx-3) { */
/*     nfx=0.;//-0.0001*PetscSign(coor[k][j][i].x-coor[k][j][i-1].x); */
/*   } */
  
  if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
    
    pc[0] = coor[k  ][j  ][i  ];
    pc[1] = coor[k  ][j  ][i+1];
    pc[2] = coor[k  ][j+1][i+1];
    pc[3] = coor[k  ][j+1][i  ];
    
    pc[4] = coor[k+1][j  ][i  ];
    pc[5] = coor[k+1][j  ][i+1];
    pc[6] = coor[k+1][j+1][i+1];
    pc[7] = coor[k+1][j+1][i  ];
        
    kp[0]=k  ;jp[0]=j  ;ip[0]=i  ;
    kp[1]=k  ;jp[1]=j  ;ip[1]=i+1;
    kp[2]=k  ;jp[2]=j+1;ip[2]=i+1;
    kp[3]=k  ;jp[3]=j+1;ip[3]=i  ;
    kp[4]=k+1;jp[4]=j  ;ip[4]=i  ;
    kp[5]=k+1;jp[5]=j  ;ip[5]=i+1;
    kp[6]=k+1;jp[6]=j+1;ip[6]=i+1;
    kp[7]=k+1;jp[7]=j+1;ip[7]=i  ;
    
    fsi_2nd_InterceptionPoint2(p,pc, nfx, nfy,
			      nfz, ibminfo, elmt);
    
    switch (ibminfo[elmt].ssmode) {
    case(0): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[1]; ibminfo[elmt].jj22 = jp[1]; ibminfo[elmt].kk22 = kp[1];
      ibminfo[elmt].ii33=ip[3]; ibminfo[elmt].jj33 = jp[3]; ibminfo[elmt].kk33 = kp[3];
      break;
    }
    case (1): {
      ibminfo[elmt].ii11=ip[1]; ibminfo[elmt].jj11 = jp[1]; ibminfo[elmt].kk11 = kp[1];
      ibminfo[elmt].ii22=ip[2]; ibminfo[elmt].jj22 = jp[2]; ibminfo[elmt].kk22 = kp[2];
      ibminfo[elmt].ii33=ip[3]; ibminfo[elmt].jj33 = jp[3]; ibminfo[elmt].kk33 = kp[3];
      break;
    }
    case (2): {
      ibminfo[elmt].ii11=ip[4]; ibminfo[elmt].jj11 = jp[4]; ibminfo[elmt].kk11 = kp[4];
      ibminfo[elmt].ii22=ip[5]; ibminfo[elmt].jj22 = jp[5]; ibminfo[elmt].kk22 = kp[5];
      ibminfo[elmt].ii33=ip[7]; ibminfo[elmt].jj33 = jp[7]; ibminfo[elmt].kk33 = kp[7];
      break;
    }
    case (3): {
      ibminfo[elmt].ii11=ip[5]; ibminfo[elmt].jj11 = jp[5]; ibminfo[elmt].kk11 = kp[5];
      ibminfo[elmt].ii22=ip[6]; ibminfo[elmt].jj22 = jp[6]; ibminfo[elmt].kk22 = kp[6];
      ibminfo[elmt].ii33=ip[7]; ibminfo[elmt].jj33 = jp[7]; ibminfo[elmt].kk33 = kp[7];
      break;
    }
    case (4): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[7]; ibminfo[elmt].jj22 = jp[7]; ibminfo[elmt].kk22 = kp[7];
      ibminfo[elmt].ii33=ip[3]; ibminfo[elmt].jj33 = jp[3]; ibminfo[elmt].kk33 = kp[3];
      break;
    }
    case (5): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[7]; ibminfo[elmt].jj33 = jp[7]; ibminfo[elmt].kk33 = kp[7];
      break;
    }
    case (6): {
      ibminfo[elmt].ii11=ip[1]; ibminfo[elmt].jj11 = jp[1]; ibminfo[elmt].kk11 = kp[1];
      ibminfo[elmt].ii22=ip[5]; ibminfo[elmt].jj22 = jp[5]; ibminfo[elmt].kk22 = kp[5];
      ibminfo[elmt].ii33=ip[6]; ibminfo[elmt].jj33 = jp[6]; ibminfo[elmt].kk33 = kp[6];
      break;
    }
    case (7): {
      ibminfo[elmt].ii11=ip[2]; ibminfo[elmt].jj11 = jp[2]; ibminfo[elmt].kk11 = kp[2];
      ibminfo[elmt].ii22=ip[6]; ibminfo[elmt].jj22 = jp[6]; ibminfo[elmt].kk22 = kp[6];
      ibminfo[elmt].ii33=ip[1]; ibminfo[elmt].jj33 = jp[1]; ibminfo[elmt].kk33 = kp[1];
      break;
    }
    case (8): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[5]; ibminfo[elmt].jj22 = jp[5]; ibminfo[elmt].kk22 = kp[5];
      ibminfo[elmt].ii33=ip[1]; ibminfo[elmt].jj33 = jp[1]; ibminfo[elmt].kk33 = kp[1];
      break;
    }
    case (9): {
      ibminfo[elmt].ii11=ip[0]; ibminfo[elmt].jj11 = jp[0]; ibminfo[elmt].kk11 = kp[0];
      ibminfo[elmt].ii22=ip[4]; ibminfo[elmt].jj22 = jp[4]; ibminfo[elmt].kk22 = kp[4];
      ibminfo[elmt].ii33=ip[5]; ibminfo[elmt].jj33 = jp[5]; ibminfo[elmt].kk33 = kp[5];
      break;
    }
    case (10): {
      ibminfo[elmt].ii11=ip[3]; ibminfo[elmt].jj11 = jp[3]; ibminfo[elmt].kk11 = kp[3];
      ibminfo[elmt].ii22=ip[7]; ibminfo[elmt].jj22 = jp[7]; ibminfo[elmt].kk22 = kp[7];
      ibminfo[elmt].ii33=ip[6]; ibminfo[elmt].jj33 = jp[6]; ibminfo[elmt].kk33 = kp[6];
      break;
    }
    case (11): {
      ibminfo[elmt].ii11=ip[2]; ibminfo[elmt].jj11 = jp[2]; ibminfo[elmt].kk11 = kp[2];
      ibminfo[elmt].ii22=ip[3]; ibminfo[elmt].jj22 = jp[3]; ibminfo[elmt].kk22 = kp[3];
      ibminfo[elmt].ii33=ip[6]; ibminfo[elmt].jj33 = jp[6]; ibminfo[elmt].kk33 = kp[6];
      break;
    }
    }
    if (ibminfo[elmt].ssmode<0) 
      PetscPrintf(PETSC_COMM_SELF, "FSI 2nd Interpolation Coeffients 2 Were not Found!!!! %d  %d %le \n", elmt, ibminfo[elmt].ssmode,ibminfo[elmt].d_ss);
    //PetscPrintf(PETSC_COMM_SELF, "FSI Interpolatoion host %d  %d  %d \n", ibminfo[elmt].i1,ibminfo[elmt].j1,ibminfo[elmt].k1);
    
  } //if 

  DAVecRestoreArray(fda, user->lCent,&coor);  
  //DAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode fsi_interpolation_coeff(UserCtx *user, IBMNodes *ibm, IBMInfo *fsi_intp,SurfElmtInfo *elmtinfo, FSInfo *fsi)
{
/*  Note: this subroutine needs the information of ibmnodes, 
    Therefore shouldn't be called before ibm_search and 
    FsiInitialize */

  Closest_NearBndryPt_ToSurfElmt(user, ibm, elmtinfo, fsi, 0);
  PetscPrintf(PETSC_COMM_WORLD, "Closest nbn\n"); 
  //PetscBarrier(PETSC_NULL);

  GridCellaroundSurElmt(user, ibm, elmtinfo);
  //PetscBarrier(PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Closest grid cell\n" ); 
  Find_fsi_interp_Coeff(fsi_intp, user, ibm, elmtinfo);
  PetscBarrier(PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "fsi interp coeff\n" ); 

  Find_fsi_interp_Coeff2(fsi_intp, user, ibm, elmtinfo);
  PetscPrintf(PETSC_COMM_WORLD, "fsi interp coeff 2\n" ); 

  return(0);
}

PetscErrorCode Calc_fsi_surf_stress(IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo)
{
  DA	        da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt;
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscInt	i, j, k;
  PetscReal     cr1, cr2, cr3;
  PetscReal     cv1, cv2, cv3;
  PetscReal     di;
  PetscReal     nt_x, nt_y, nt_z;
  PetscReal     ns_x, ns_y, ns_z;
  PetscReal     ***p;
  Cmpnts        ***ucat, uinp;

  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(da, user->lP, &p);

  for (elmt=0; elmt<n_elmt; elmt++) {
    //PetscPrintf(PETSC_COMM_SELF, "n_P %d\n",ibm->n_P[elmt]);
    
  if (elmtinfo[elmt].n_P>0  && elmtinfo[elmt].FoundAroundcell>0) {
    ip1 = ibminfo[elmt].i1; jp1 = ibminfo[elmt].j1; kp1 = ibminfo[elmt].k1;
    ip2 = ibminfo[elmt].i2; jp2 = ibminfo[elmt].j2; kp2 = ibminfo[elmt].k2;
    ip3 = ibminfo[elmt].i3; jp3 = ibminfo[elmt].j3; kp3 = ibminfo[elmt].k3;
    
    cr1 = ibminfo[elmt].cr1; cr2 = ibminfo[elmt].cr2; cr3 = ibminfo[elmt].cr3;
    
    di  = ibminfo[elmt].d_i; 

    i=elmtinfo[elmt].icell;
    j=elmtinfo[elmt].jcell;
    k=elmtinfo[elmt].kcell;

    PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3);
    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
      //if (ip1>=xs && ip1<xe && jp1>=ys && jp1<ye && kp1>=zs && kp1<ze) {    
      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];
      elmtinfo[elmt].P= (cv1 * cr1 + cv2 * cr2 + cv3 * cr3);

/*       if (i==1) */
      	PetscPrintf(PETSC_COMM_SELF, "intp P %d %le %le %le %le %le %le %le %le\n",elmt,elmtinfo[elmt].P,cr1,cr2,cr3,cv1,cv2,cv3,di);

/*       if (fabs(elmtinfo[elmt].P)<1e-5) */
/*       	PetscPrintf(PETSC_COMM_SELF, "intp P %le %le %le %le %d %d %d %d \n",elmtinfo[elmt].P,cr1,cr2,cr3,ip1,jp1,kp1,elmt); */
 
      cv1 = ucat[kp1][jp1][ip1].x;
      cv2 = ucat[kp2][jp2][ip2].x;
      cv3 = ucat[kp3][jp3][ip3].x;
      uinp.x = (cv1 * cr1 + cv2 * cr2 + cv3 * cr3);

      cv1 = ucat[kp1][jp1][ip1].y;
      cv2 = ucat[kp2][jp2][ip2].y;
      cv3 = ucat[kp3][jp3][ip3].y;
      uinp.y = (cv1 * cr1 + cv2 * cr2 + cv3 * cr3);

      PetscPrintf(PETSC_COMM_SELF, "intp u_y %d %le %le %le %le\n",elmt,uinp.y,cv1,cv2,cv3);

      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;
      uinp.z = (cv1 * cr1 + cv2 * cr2 + cv3 * cr3);
      
      PetscPrintf(PETSC_COMM_SELF, "intp u_z %d %le %le %le %le\n",elmt,uinp.z,cv1,cv2,cv3);
	
      cv1 = ( ibm->u[ibm->nv1[elmt]].x
	     +ibm->u[ibm->nv2[elmt]].x 
	     +ibm->u[ibm->nv3[elmt]].x)/3.;
      cv2 = ( ibm->u[ibm->nv1[elmt]].y
	     +ibm->u[ibm->nv2[elmt]].y
	     +ibm->u[ibm->nv3[elmt]].y)/3.;
      cv3 = ( ibm->u[ibm->nv1[elmt]].z
	     +ibm->u[ibm->nv2[elmt]].z
	     +ibm->u[ibm->nv3[elmt]].z)/3.;
      
      ns_x= ibm->ns_x[elmt];	
      ns_y= ibm->ns_y[elmt];	
      ns_z= ibm->ns_z[elmt];
      
      nt_x= ibm->nt_x[elmt];	
      nt_y= ibm->nt_y[elmt];
      nt_z= ibm->nt_z[elmt];
      
      if (di>1e-10) {
      elmtinfo[elmt].Tow_ws=((uinp.x-cv1)*ns_x + (uinp.y-cv2)*ns_y +
			     (uinp.z-cv3)*ns_z)/di;

      elmtinfo[elmt].Tow_wt=((uinp.x-cv1)*nt_x + (uinp.y-cv2)*nt_y +
			     (uinp.z-cv3)*nt_z)/di;
      } else {
	PetscPrintf(PETSC_COMM_WORLD, "zero di %le ",di);
      }      

      PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_wt,uinp.x,uinp.y,uinp.z,di);
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_wt,nt_x,nt_y,nt_z,di);
     /*  PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,uinp.x,uinp.y,uinp.z,di); */
/*       PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,ns_x,ns_y,ns_z,di); */

    }
  }
  }
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lP, &p);
  return(0);

}

PetscErrorCode Calc_fsi_surf_stress2(IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo)
{
  DA	        da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt;
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscInt	ip11, ip22, ip33, jp11, jp22, jp33, kp11, kp22, kp33;
  PetscInt	i, j, k;

  PetscReal     cr1, cr2, cr3;
  PetscReal     cv1, cv2, cv3;
  PetscReal     cr11, cr22, cr33;
  PetscReal     cv11, cv22, cv33;
  PetscReal     cs1, cs2, cs3;
  PetscReal     cs11, cs22, cs33;

  PetscInt	iip11, iip22, iip33, jjp11, jjp22, jjp33, kkp11, kkp22, kkp33;
  PetscInt	iip1, iip2, iip3, jjp1, jjp2, jjp3, kkp1, kkp2, kkp3;
  PetscReal     ct1, ct2, ct3;
  PetscReal     ct11, ct22, ct33;
  PetscReal     ds,sd;

  PetscReal     di;
  PetscReal     nf_x, nf_y, nf_z;
  PetscReal     nt_x, nt_y, nt_z;
  PetscReal     ns_x, ns_y, ns_z;
  PetscReal     ***p;
  Cmpnts        ***ucat, uinp;
  PetscReal	***nvert;

  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(da, user->lP, &p);
  DAVecGetArray(da, user->lNvert, &nvert);

  for (elmt=0; elmt<n_elmt; elmt++) {
    //PetscPrintf(PETSC_COMM_SELF, "n_P %d\n",ibm->n_P[elmt]);
    
  if (elmtinfo[elmt].n_P>0  && elmtinfo[elmt].FoundAroundcell>0) {
    ip1 = ibminfo[elmt].i1; jp1 = ibminfo[elmt].j1; kp1 = ibminfo[elmt].k1;
    ip2 = ibminfo[elmt].i2; jp2 = ibminfo[elmt].j2; kp2 = ibminfo[elmt].k2;
    ip3 = ibminfo[elmt].i3; jp3 = ibminfo[elmt].j3; kp3 = ibminfo[elmt].k3;
    
    cr1 = ibminfo[elmt].cr1; cr2 = ibminfo[elmt].cr2; cr3 = ibminfo[elmt].cr3;
    cs1 = ibminfo[elmt].cs1; cs2 = ibminfo[elmt].cs2; cs3 = ibminfo[elmt].cs3;

    ip11 = ibminfo[elmt].i11; jp11 = ibminfo[elmt].j11; kp11 = ibminfo[elmt].k11;
    ip22 = ibminfo[elmt].i22; jp22 = ibminfo[elmt].j22; kp22 = ibminfo[elmt].k22;
    ip33 = ibminfo[elmt].i33; jp33 = ibminfo[elmt].j33; kp33 = ibminfo[elmt].k33;
    
    cr11 = ibminfo[elmt].cr11; cr22 = ibminfo[elmt].cr22; cr33 = ibminfo[elmt].cr33;
    cs11 = ibminfo[elmt].cs11; cs22 = ibminfo[elmt].cs22; cs33 = ibminfo[elmt].cs33;

    iip1 = ibminfo[elmt].ii1; jjp1 = ibminfo[elmt].jj1; kkp1 = ibminfo[elmt].kk1;
    iip2 = ibminfo[elmt].ii2; jjp2 = ibminfo[elmt].jj2; kkp2 = ibminfo[elmt].kk2;
    iip3 = ibminfo[elmt].ii3; jjp3 = ibminfo[elmt].jj3; kkp3 = ibminfo[elmt].kk3;
    
    iip11 = ibminfo[elmt].ii11; jjp11 = ibminfo[elmt].jj11; kkp11 = ibminfo[elmt].kk11;
    iip22 = ibminfo[elmt].ii22; jjp22 = ibminfo[elmt].jj22; kkp22 = ibminfo[elmt].kk22;
    iip33 = ibminfo[elmt].ii33; jjp33 = ibminfo[elmt].jj33; kkp33 = ibminfo[elmt].kk33;

    ct1 = ibminfo[elmt].ct1; ct2 = ibminfo[elmt].ct2; ct3 = ibminfo[elmt].ct3;
    ct11 = ibminfo[elmt].ct11; ct22 = ibminfo[elmt].ct22; ct33 = ibminfo[elmt].ct33;


/*     di  = 0.5*(ibminfo[elmt].d_i+ibminfo[elmt].d_ii);  */
    di  = (ibminfo[elmt].d_i); 
    if (elmtinfo[elmt].Need3rdPoint) {
      ds  = (ibminfo[elmt].d_s);
      di  = di + ds;
    }

    i=elmtinfo[elmt].icell;
    j=elmtinfo[elmt].jcell;
    k=elmtinfo[elmt].kcell;

/*     PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3); */
    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
      //if (ip1>=xs && ip1<xe && jp1>=ys && jp1<ye && kp1>=zs && kp1<ze) {    
/*       cv1 = p[kp1][jp1][ip1]; */
/*       cv2 = p[kp2][jp2][ip2]; */
/*       cv3 = p[kp3][jp3][ip3]; */
/*       cv11 = p[kp11][jp11][ip11]; */
/*       cv22 = p[kp22][jp22][ip22]; */
/*       cv33 = p[kp33][jp33][ip33]; */

/*       elmtinfo[elmt].P= 0.5*(cv1 * cr1 + cv2 * cr2 + cv3 * cr3 + */
/* 			     cv11*cr11 + cv22*cr22 + cv33*cr33); */

      if (elmtinfo[elmt].Need3rdPoint) {
      cv1 = p[kkp1][jjp1][iip1];
      cv2 = p[kkp2][jjp2][iip2];
      cv3 = p[kkp3][jjp3][iip3];
      elmtinfo[elmt].P = cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      } else {

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];
/*       cv11 = p[kp11][jp11][ip11]; */
/*       cv22 = p[kp22][jp22][ip22]; */
/*       cv33 = p[kp33][jp33][ip33]; */

/*       elmtinfo[elmt].P= 0.5*(cv1 * cr1 + cv2 * cr2 + cv3 * cr3 + */
/* 			     cv11*cr11 + cv22*cr22 + cv33*cr33); */
      elmtinfo[elmt].P= cv1 * cr1 + cv2 * cr2 + cv3 * cr3 ;
      }

/*       if ((int)(nvert[kp1][jp1][ip1]+0.5)>1) { */
/* 	ucat[kp1][jp1][ip1].z=0.; */
/* 	ucat[kp1][jp1][ip1].y=0.; */
/* 	ucat[kp1][jp1][ip1].x=0.; */
/*       } */
/*       if ((int)(nvert[kp2][jp2][ip2]+0.5)>1) { */
/* 	ucat[kp2][jp2][ip2].z=0.; */
/* 	ucat[kp2][jp2][ip2].y=0.; */
/* 	ucat[kp2][jp2][ip2].x=0.; */
/*       } */
/*       if ((int)(nvert[kp3][jp3][ip3]+0.5)>1) { */
/* 	ucat[kp3][jp3][ip3].z=0.; */
/* 	ucat[kp3][jp3][ip3].y=0.; */
/* 	ucat[kp3][jp3][ip3].x=0.; */
/*       } */

/*      if ((int)(nvert[kkp1][jjp1][iip1]+0.5)>1) { */
/* 	ucat[kkp1][jjp1][iip1].z=0.; */
/* 	ucat[kkp1][jjp1][iip1].y=0.; */
/* 	ucat[kkp1][jjp1][iip1].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp2][jjp2][iip2]+0.5)>1) { */
/* 	ucat[kkp2][jjp2][iip2].z=0.; */
/* 	ucat[kkp2][jjp2][iip2].y=0.; */
/* 	ucat[kkp2][jjp2][iip2].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp3][jjp3][iip3]+0.5)>1) { */
/* 	ucat[kkp3][jjp3][iip3].z=0.; */
/* 	ucat[kkp3][jjp3][iip3].y=0.; */
/* 	ucat[kkp3][jjp3][iip3].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp11][jjp11][iip11]+0.5)>1) { */
/* 	ucat[kkp11][jjp11][iip11].z=0.; */
/* 	ucat[kkp11][jjp11][iip11].y=0.; */
/* 	ucat[kkp11][jjp11][iip11].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp22][jjp22][iip22]+0.5)>1) { */
/* 	ucat[kkp22][jjp22][iip22].z=0.; */
/* 	ucat[kkp22][jjp22][iip22].y=0.; */
/* 	ucat[kkp22][jjp22][iip22].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp33][jjp33][iip33]+0.5)>1) { */
/* 	ucat[kkp33][jjp33][iip33].z=0.; */
/* 	ucat[kkp33][jjp33][iip33].y=0.; */
/* 	ucat[kkp33][jjp33][iip33].x=0.; */
/*       } */

      if (elmt==184 || elmt==231) 
      	PetscPrintf(PETSC_COMM_SELF, "intp P %d %le %le %le %le %le %le %le %le\n",elmt,elmtinfo[elmt].P,cr1,cr2,cr3,cv1,cv2,cv3,di);

/*       if (fabs(elmtinfo[elmt].P)<1e-5) */
/*       	PetscPrintf(PETSC_COMM_SELF, "intp P %le %le %le %le %d %d %d %d \n",elmtinfo[elmt].P,cr1,cr2,cr3,ip1,jp1,kp1,elmt); */
 
      if (elmtinfo[elmt].Need3rdPoint) {

      cv1 = ucat[kkp1][jjp1][iip1].x;
      cv2 = ucat[kkp2][jjp2][iip2].x;
      cv3 = ucat[kkp3][jjp3][iip3].x;
/*       cv11 = ucat[kkp11][jjp11][iip11].x; */
/*       cv22 = ucat[kkp22][jjp22][iip22].x; */
/*       cv33 = ucat[kkp33][jjp33][iip33].x; */

/*       uinp.x = 0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		    cv11*ct11 + cv22*ct22 + cv33*ct33); */
/*       uinp.x = 0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		    cv11*cs11 + cv22*cs22 + cv33*cs33); */
      uinp.x = cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;

      } else {
      cv1 = ucat[kp1][jp1][ip1].x;
      cv2 = ucat[kp2][jp2][ip2].x;
      cv3 = ucat[kp3][jp3][ip3].x;
/*       cv11 = ucat[kp11][jp11][ip11].x; */
/*       cv22 = ucat[kp22][jp22][ip22].x; */
/*       cv33 = ucat[kp33][jp33][ip33].x; */
      uinp.x = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3 );
      }

      if (elmt==184 || elmt==231){
      PetscPrintf(PETSC_COMM_SELF, "intp u_x 1  %d %le %le %le %le %le %le %le\n",elmt,uinp.x,cv1,cv2,cv3,ct1,ct2,ct3);
      PetscPrintf(PETSC_COMM_SELF, "intp u_x 11 %d %le %le %le %le %le %le %le\n",elmt,uinp.x,cv11,cv22,cv33,ct11,ct22,ct33);
      }

      if (elmtinfo[elmt].Need3rdPoint) {
      cv1 = ucat[kkp1][jjp1][iip1].y;
      cv2 = ucat[kkp2][jjp2][iip2].y;
      cv3 = ucat[kkp3][jjp3][iip3].y;
/*       cv11 = ucat[kkp11][jjp11][iip11].y; */
/*       cv22 = ucat[kkp22][jjp22][iip22].y; */
/*       cv33 = ucat[kkp33][jjp33][iip33].y; */

/*       uinp.y =  0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		     cv11*ct11 + cv22*ct22 + cv33*ct33); */
/*       uinp.y =  0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		     cv11*cs11 + cv22*cs22 + cv33*cs33); */
      uinp.y =  cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      } else {	
      cv1 = ucat[kp1][jp1][ip1].y;
      cv2 = ucat[kp2][jp2][ip2].y;
      cv3 = ucat[kp3][jp3][ip3].y;
/*       cv11 = ucat[kp11][jp11][ip11].y; */
/*       cv22 = ucat[kp22][jp22][ip22].y; */
/*       cv33 = ucat[kp33][jp33][ip33].y; */
      uinp.y =  (cv1 * cs1 + cv2 * cs2 + cv3 * cs3 );
      }

      if (elmt==184 || elmt==231){
      PetscPrintf(PETSC_COMM_SELF, "intp u_y 1  %d %le %le %le %le %le %le %le\n",elmt,uinp.y,cv1,cv2,cv3,ct1,ct2,ct3);
      PetscPrintf(PETSC_COMM_SELF, "intp u_y 11 %d %le %le %le %le %le %le %le\n",elmt,uinp.y,cv11,cv22,cv33,ct11,ct22,ct33);
      }

      if (elmtinfo[elmt].Need3rdPoint) {
      cv1 = ucat[kkp1][jjp1][iip1].z;
      cv2 = ucat[kkp2][jjp2][iip2].z;
      cv3 = ucat[kkp3][jjp3][iip3].z;
/*       cv11 = ucat[kkp11][jjp11][iip11].z; */
/*       cv22 = ucat[kkp22][jjp22][iip22].z; */
/*       cv33 = ucat[kkp33][jjp33][iip33].z; */

/*       uinp.z =  0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		     cv11*cs11 + cv22*cs22 + cv33*cs33); */

      uinp.z =  cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      } else {
      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;
/*       cv11 = ucat[kp11][jp11][ip11].z; */
/*       cv22 = ucat[kp22][jp22][ip22].z; */
/*       cv33 = ucat[kp33][jp33][ip33].z; */
      uinp.z =  (cv1 * cs1 + cv2 * cs2 + cv3 * cs3 );
      }
/*       uinp.z =  0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		     cv11*ct11 + cv22*ct22 + cv33*ct33); */

      if (elmt==184 || elmt==231){
      PetscPrintf(PETSC_COMM_SELF, "intp u_z 1  %d %le %le %le %le %le %le %le\n",elmt,uinp.z,cv1,cv2,cv3,ct1,ct2,ct3);
      PetscPrintf(PETSC_COMM_SELF, "intp u_z 11 %d %le %le %le %le %le %le %le\n",elmt,uinp.z,cv11,cv22,cv33,ct11,ct22,ct33);
      }
      cv1 = ( ibm->u[ibm->nv1[elmt]].x
	     +ibm->u[ibm->nv2[elmt]].x 
	     +ibm->u[ibm->nv3[elmt]].x)/3.;
      cv2 = ( ibm->u[ibm->nv1[elmt]].y
	     +ibm->u[ibm->nv2[elmt]].y
	     +ibm->u[ibm->nv3[elmt]].y)/3.;
      cv3 = ( ibm->u[ibm->nv1[elmt]].z
	     +ibm->u[ibm->nv2[elmt]].z
	     +ibm->u[ibm->nv3[elmt]].z)/3.;
      
      ns_x= ibm->ns_x[elmt];	
      ns_y= ibm->ns_y[elmt];	
      ns_z= ibm->ns_z[elmt];
      
      nt_x= ibm->nt_x[elmt];	
      nt_y= ibm->nt_y[elmt];
      nt_z= ibm->nt_z[elmt];
      
      nf_x= ibm->nf_x[elmt];	
      nf_y= ibm->nf_y[elmt];
      nf_z= ibm->nf_z[elmt];
      
      if (di>1e-10) {
      elmtinfo[elmt].Tow_ws=    ((uinp.x-cv1)*ns_x + 
				 (uinp.y-cv2)*ns_y +
				 (uinp.z-cv3)*ns_z)/di;

      elmtinfo[elmt].Tow_wt=    ((uinp.x-cv1)*nt_x + 
				 (uinp.y-cv2)*nt_y +
				 (uinp.z-cv3)*nt_z)/di;
      
      elmtinfo[elmt].Tow_wn=    ((uinp.x-cv1)*nf_x + 
				 (uinp.y-cv2)*nf_y +
				 (uinp.z-cv3)*nf_z)/di;
      
      } else {
	PetscPrintf(PETSC_COMM_WORLD, "zero di %le ",di);
      }      

      if (elmt==184 || elmt==231) {
/*       PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3); */
      PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",iip1,jjp1,kkp1,iip2,jjp2,kkp2,iip3,jjp3,kkp3);
      PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",iip11,jjp11,kkp11,iip22,jjp22,kkp22,iip33,jjp33,kkp33);
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_wt,uinp.x,uinp.y,uinp.z,di);
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_wt,nt_x,nt_y,nt_z,di);
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,uinp.x,uinp.y,uinp.z,di);
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,ns_x,ns_y,ns_z,di);
      }
    }
  }
  }

  PetscInt rank;
  PetscInt n_v=ibm->n_v;
  PetscInt ti=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    //if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "Stress%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,tow_t,tow_s,p,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z,di\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-16]=CELLCENTERED)\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
      }
      for (i=0; i<n_v; i++) {	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].Tow_wt);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].Tow_ws);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].P);
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
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibminfo[i].d_i);
      }

      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
      //}
  }

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lP, &p);
  return(0);

}

PetscErrorCode Calc_fsi_surf_stress_advanced(IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo)
{
  DA	        da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt      n_elmt = ibm->n_elmt;
  PetscInt      elmt;
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscInt	ip11, ip22, ip33, jp11, jp22, jp33, kp11, kp22, kp33;
  PetscInt	iip1, iip2, iip3, jjp1, jjp2, jjp3, kkp1, kkp2, kkp3;
  PetscInt	iip11, iip22, iip33, jjp11, jjp22, jjp33, kkp11, kkp22, kkp33;

  PetscInt	i, j, k;

  PetscReal     cr1, cr2, cr3;
  PetscReal     cv1, cv2, cv3;
  PetscReal     cr11, cr22, cr33;
  PetscReal     cv11, cv22, cv33;
  PetscReal     cs1, cs2, cs3;
  PetscReal     cs11, cs22, cs33;
  PetscReal     ct1, ct2, ct3;
  PetscReal     ct11, ct22, ct33;
  PetscReal     phia,phib,phic,phid;
  PetscReal     di,sd,sc,sb;
  PetscReal     nt_x, nt_y, nt_z;
  PetscReal     ns_x, ns_y, ns_z;
  PetscReal     ***p;
  Cmpnts        ***ucat, uinp, uinp2, uinp3;
  PetscReal	***nvert;
  PetscReal     lhs[3][3], rhs_l[3][3];
  PetscReal     dtm;

  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(da, user->lP, &p);
  DAVecGetArray(da, user->lNvert, &nvert);

  for (elmt=0; elmt<n_elmt; elmt++) {
    elmtinfo[elmt].P=0.;
    elmtinfo[elmt].Tow_wt=0.;
    elmtinfo[elmt].Tow_ws=0.;
  }

  for (elmt=0; elmt<n_elmt; elmt++) {
    //PetscPrintf(PETSC_COMM_SELF, "n_P %d\n",ibm->n_P[elmt]);
    
  if (elmtinfo[elmt].n_P>0  && elmtinfo[elmt].FoundAroundcell>0) {
    ip1 = ibminfo[elmt].i1; jp1 = ibminfo[elmt].j1; kp1 = ibminfo[elmt].k1;
    ip2 = ibminfo[elmt].i2; jp2 = ibminfo[elmt].j2; kp2 = ibminfo[elmt].k2;
    ip3 = ibminfo[elmt].i3; jp3 = ibminfo[elmt].j3; kp3 = ibminfo[elmt].k3;
    
    cr1 = ibminfo[elmt].cr1; cr2 = ibminfo[elmt].cr2; cr3 = ibminfo[elmt].cr3;
    cs1 = ibminfo[elmt].cs1; cs2 = ibminfo[elmt].cs2; cs3 = ibminfo[elmt].cs3;

    ip11 = ibminfo[elmt].i11; jp11 = ibminfo[elmt].j11; kp11 = ibminfo[elmt].k11;
    ip22 = ibminfo[elmt].i22; jp22 = ibminfo[elmt].j22; kp22 = ibminfo[elmt].k22;
    ip33 = ibminfo[elmt].i33; jp33 = ibminfo[elmt].j33; kp33 = ibminfo[elmt].k33;
    
    cr11 = ibminfo[elmt].cr11; cr22 = ibminfo[elmt].cr22; cr33 = ibminfo[elmt].cr33;
    cs11 = ibminfo[elmt].cs11; cs22 = ibminfo[elmt].cs22; cs33 = ibminfo[elmt].cs33;

    iip11 = ibminfo[elmt].ii11; jjp11 = ibminfo[elmt].jj11; kkp11 = ibminfo[elmt].kk11;
    iip22 = ibminfo[elmt].ii22; jjp22 = ibminfo[elmt].jj22; kkp22 = ibminfo[elmt].kk22;
    iip33 = ibminfo[elmt].ii33; jjp33 = ibminfo[elmt].jj33; kkp33 = ibminfo[elmt].kk33;

    iip1 = ibminfo[elmt].ii1; jjp1 = ibminfo[elmt].jj1; kkp1 = ibminfo[elmt].kk1;
    iip2 = ibminfo[elmt].ii2; jjp2 = ibminfo[elmt].jj2; kkp2 = ibminfo[elmt].kk2;
    iip3 = ibminfo[elmt].ii3; jjp3 = ibminfo[elmt].jj3; kkp3 = ibminfo[elmt].kk3;

    ct1 = ibminfo[elmt].ct1; ct2 = ibminfo[elmt].ct2; ct3 = ibminfo[elmt].ct3;
    ct11 = ibminfo[elmt].ct11; ct22 = ibminfo[elmt].ct22; ct33 = ibminfo[elmt].ct33;
    
    di  = 0.5*(ibminfo[elmt].d_i+ibminfo[elmt].d_ii); 
    
/*   // just check ??? */
/*     elmtinfo[elmt].Need3rdPoint=1; */
    
    if (elmtinfo[elmt].Need3rdPoint) {
      di = di + ibminfo[elmt].d_s; // distance to 2nd pt
      sd = di + ibminfo[elmt].d_ss; // distance to 3rd pt
    } else {
     /*  ds  = 0.5*(ibminfo[elmt].d_s+ibminfo[elmt].d_ss);*/
      sd = ibminfo[elmt].d_s + di;
/*       sb = di; */
/*       sc = ibminfo[elmt].d_s  + sb; */
/*       sd = ibminfo[elmt].d_ss + sc; */
    }

    i=elmtinfo[elmt].icell;
    j=elmtinfo[elmt].jcell;
    k=elmtinfo[elmt].kcell;

/*     PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %le %le %le %le %le %le\n",elmt,ct1,ct2,ct3,ct11,ct22,ct33); */
/*     PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",iip1,jjp1,kkp1,iip2,jjp2,kkp2,iip3,jjp3,kkp3); */
/*     PetscPrintf(PETSC_COMM_SELF, "intp hosts di %d %d %d %d %d %d %d %d %d\n",ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3); */
    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {
      //if (ip1>=xs && ip1<xe && jp1>=ys && jp1<ye && kp1>=zs && kp1<ze) {    
/*       if ((int)(nvert[kp1][jp1][ip1]+0.5)>1) { */
/* 	ucat[kp1][jp1][ip1].z=0.; */
/* 	ucat[kp1][jp1][ip1].y=0.; */
/* 	ucat[kp1][jp1][ip1].x=0.; */
/*       } */
/*       if ((int)(nvert[kp2][jp2][ip2]+0.5)>1) { */
/* 	ucat[kp2][jp2][ip2].z=0.; */
/* 	ucat[kp2][jp2][ip2].y=0.; */
/* 	ucat[kp2][jp2][ip2].x=0.; */
/*       } */
/*       if ((int)(nvert[kp3][jp3][ip3]+0.5)>1) { */
/* 	ucat[kp3][jp3][ip3].z=0.; */
/* 	ucat[kp3][jp3][ip3].y=0.; */
/* 	ucat[kp3][jp3][ip3].x=0.; */
/*       } */
/*       if ((int)(nvert[kp11][jp11][ip11]+0.5)>1) { */
/* 	ucat[kp11][jp11][ip11].z=0.; */
/* 	ucat[kp11][jp11][ip11].y=0.; */
/* 	ucat[kp11][jp11][ip11].x=0.; */
/*       } */
/*       if ((int)(nvert[kp22][jp22][ip22]+0.5)>1) { */
/* 	ucat[kp22][jp22][ip22].z=0.; */
/* 	ucat[kp22][jp22][ip22].y=0.; */
/* 	ucat[kp22][jp22][ip22].x=0.; */
/*       } */
/*       if ((int)(nvert[kp33][jp33][ip33]+0.5)>1) { */
/* 	ucat[kp33][jp33][ip33].z=0.; */
/* 	ucat[kp33][jp33][ip33].y=0.; */
/* 	ucat[kp33][jp33][ip33].x=0.; */
/*       } */

/*      if ((int)(nvert[kkp1][jjp1][iip1]+0.5)>1) { */
/* 	ucat[kkp1][jjp1][iip1].z=0.; */
/* 	ucat[kkp1][jjp1][iip1].y=0.; */
/* 	ucat[kkp1][jjp1][iip1].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp2][jjp2][iip2]+0.5)>1) { */
/* 	ucat[kkp2][jjp2][iip2].z=0.; */
/* 	ucat[kkp2][jjp2][iip2].y=0.; */
/* 	ucat[kkp2][jjp2][iip2].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp3][jjp3][iip3]+0.5)>1) { */
/* 	ucat[kkp3][jjp3][iip3].z=0.; */
/* 	ucat[kkp3][jjp3][iip3].y=0.; */
/* 	ucat[kkp3][jjp3][iip3].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp11][jjp11][iip11]+0.5)>1) { */
/* 	ucat[kkp11][jjp11][iip11].z=0.; */
/* 	ucat[kkp11][jjp11][iip11].y=0.; */
/* 	ucat[kkp11][jjp11][iip11].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp22][jjp22][iip22]+0.5)>1) { */
/* 	ucat[kkp22][jjp22][iip22].z=0.; */
/* 	ucat[kkp22][jjp22][iip22].y=0.; */
/* 	ucat[kkp22][jjp22][iip22].x=0.; */
/*       } */
/*       if ((int)(nvert[kkp33][jjp33][iip33]+0.5)>1) { */
/* 	ucat[kkp33][jjp33][iip33].z=0.; */
/* 	ucat[kkp33][jjp33][iip33].y=0.; */
/* 	ucat[kkp33][jjp33][iip33].x=0.; */
/*       } */
      
      if (elmtinfo[elmt].Need3rdPoint) {
      cv1 = p[kkp1][jjp1][iip1];
      cv2 = p[kkp2][jjp2][iip2];
      cv3 = p[kkp3][jjp3][iip3];
      elmtinfo[elmt].P = cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      } else {

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];
      cv11 = p[kp11][jp11][ip11];
      cv22 = p[kp22][jp22][ip22];
      cv33 = p[kp33][jp33][ip33];

      elmtinfo[elmt].P= 0.5*(cv1 * cr1 + cv2 * cr2 + cv3 * cr3 +
			     cv11*cr11 + cv22*cr22 + cv33*cr33);
      }
/*       if (i==1) */
/*       	PetscPrintf(PETSC_COMM_SELF, "intp P %d %le %le %le %le %le %le %le %le\n",elmt,elmtinfo[elmt].P,cr1,cr2,cr3,cv1,cv2,cv3,di); */

/*       if (fabs(elmtinfo[elmt].P)<1e-5) */
/*       	PetscPrintf(PETSC_COMM_SELF, "intp P %le %le %le %le %d %d %d %d \n",elmtinfo[elmt].P,cr1,cr2,cr3,ip1,jp1,kp1,elmt); */
 
      cv1 = ucat[kp1][jp1][ip1].x;
      cv2 = ucat[kp2][jp2][ip2].x;
      cv3 = ucat[kp3][jp3][ip3].x;
      cv11 = ucat[kp11][jp11][ip11].x;
      cv22 = ucat[kp22][jp22][ip22].x;
      cv33 = ucat[kp33][jp33][ip33].x;

      uinp.x = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
/*       uinp.x = 0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		    cv11*cs11 + cv22*cs22 + cv33*cs33); */
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
	PetscPrintf(PETSC_COMM_SELF, "intp u_x %d %le %le %le %le\n",elmt,uinp.x,cv1,cv2,cv3);

      cv1 = ucat[kp1][jp1][ip1].y;
      cv2 = ucat[kp2][jp2][ip2].y;
      cv3 = ucat[kp3][jp3][ip3].y;
      cv11 = ucat[kp11][jp11][ip11].y;
      cv22 = ucat[kp22][jp22][ip22].y;
      cv33 = ucat[kp33][jp33][ip33].y;

      uinp.y = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
/*       uinp.y =  0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		     cv11*cs11 + cv22*cs22 + cv33*cs33); */
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
	PetscPrintf(PETSC_COMM_SELF, "intp u_y %d %le %le %le %le\n",elmt,uinp.y,cv1,cv2,cv3);

      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;
      cv11 = ucat[kp11][jp11][ip11].z;
      cv22 = ucat[kp22][jp22][ip22].z;
      cv33 = ucat[kp33][jp33][ip33].z;

      uinp.z = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
/*       uinp.z =  0.5*(cv1 * cs1 + cv2 * cs2 + cv3 * cs3 + */
/* 		     cv11*cs11 + cv22*cs22 + cv33*cs33); */
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
	PetscPrintf(PETSC_COMM_SELF, "intp u_z %d %le %le %le %le\n",elmt,uinp.z,cv1,cv2,cv3);

      // 2nd pt
      cv1 = ucat[kkp1][jjp1][iip1].x;
      cv2 = ucat[kkp2][jjp2][iip2].x;
      cv3 = ucat[kkp3][jjp3][iip3].x;
/*       cv11 = ucat[kkp11][jjp11][iip11].x; */
/*       cv22 = ucat[kkp22][jjp22][iip22].x; */
/*       cv33 = ucat[kkp33][jjp33][iip33].x; */

/*       uinp2.x = 0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		    cv11*ct11 + cv22*ct22 + cv33*ct33); */
      uinp2.x = cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 2 u_x %d %le %le %le %le\n",elmt,uinp2.x,cv1,cv2,cv3);


      cv1 = ucat[kkp1][jjp1][iip1].y;
      cv2 = ucat[kkp2][jjp2][iip2].y;
      cv3 = ucat[kkp3][jjp3][iip3].y;
/*       cv11 = ucat[kkp11][jjp11][iip11].y; */
/*       cv22 = ucat[kkp22][jjp22][iip22].y; */
/*       cv33 = ucat[kkp33][jjp33][iip33].y; */

/*       uinp2.y =  0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		     cv11*ct11 + cv22*ct22 + cv33*ct33); */
      uinp2.y =  cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 2 u_y %d %le %le %le %le\n",elmt,uinp2.y,cv1,cv2,cv3);

      cv1 = ucat[kkp1][jjp1][iip1].z;
      cv2 = ucat[kkp2][jjp2][iip2].z;
      cv3 = ucat[kkp3][jjp3][iip3].z;
/*       cv11 = ucat[kkp11][jjp11][iip11].z; */
/*       cv22 = ucat[kkp22][jjp22][iip22].z; */
/*       cv33 = ucat[kkp33][jjp33][iip33].z; */

/*       uinp2.z =  0.5*(cv1 * ct1 + cv2 * ct2 + cv3 * ct3 + */
/* 		     cv11*ct11 + cv22*ct22 + cv33*ct33); */
      uinp2.z =  cv1 * ct1 + cv2 * ct2 + cv3 * ct3 ;
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 2 u_z %d %le %le %le %le\n",elmt,uinp2.z,cv1,cv2,cv3);

      // 3rd pt
/*       cv1 = ucat[kkp1][jjp1][iip1].x; */
/*       cv2 = ucat[kkp2][jjp2][iip2].x; */
/*       cv3 = ucat[kkp3][jjp3][iip3].x; */
      cv11 = ucat[kkp11][jjp11][iip11].x;
      cv22 = ucat[kkp22][jjp22][iip22].x;
      cv33 = ucat[kkp33][jjp33][iip33].x;

      uinp3.x = (cv11*ct11 + cv22*ct22 + cv33*ct33);
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 3 u_x %d %le %le %le %le\n",elmt,uinp3.x,cv1,cv2,cv3);

/*       cv1 = ucat[kkp1][jjp1][iip1].y; */
/*       cv2 = ucat[kkp2][jjp2][iip2].y; */
/*       cv3 = ucat[kkp3][jjp3][iip3].y; */
      cv11 = ucat[kkp11][jjp11][iip11].y;
      cv22 = ucat[kkp22][jjp22][iip22].y;
      cv33 = ucat[kkp33][jjp33][iip33].y;

      uinp3.y = (cv11*ct11 + cv22*ct22 + cv33*ct33);
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 3 u_y %d %le %le %le %le\n",elmt,uinp3.y,cv1,cv2,cv3);

/*       cv1 = ucat[kkp1][jjp1][iip1].z; */
/*       cv2 = ucat[kkp2][jjp2][iip2].z; */
/*       cv3 = ucat[kkp3][jjp3][iip3].z; */
      cv11 = ucat[kkp11][jjp11][iip11].z;
      cv22 = ucat[kkp22][jjp22][iip22].z;
      cv33 = ucat[kkp33][jjp33][iip33].z;

      uinp3.z = (cv11*ct11 + cv22*ct22 + cv33*ct33);
      if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp 3 u_z %d %le %le %le %le\n",elmt,uinp3.z,cv1,cv2,cv3);
      
/*       PetscPrintf(PETSC_COMM_SELF, "intp u_z %d %le %le %le %le\n",elmt,uinp.z,cv1,cv2,cv3); */
	
      cv1 = ( ibm->u[ibm->nv1[elmt]].x
	     +ibm->u[ibm->nv2[elmt]].x 
	     +ibm->u[ibm->nv3[elmt]].x)/3.;
      cv2 = ( ibm->u[ibm->nv1[elmt]].y
	     +ibm->u[ibm->nv2[elmt]].y
	     +ibm->u[ibm->nv3[elmt]].y)/3.;
      cv3 = ( ibm->u[ibm->nv1[elmt]].z
	     +ibm->u[ibm->nv2[elmt]].z
	     +ibm->u[ibm->nv3[elmt]].z)/3.;
     
      ns_x= ibm->ns_x[elmt];	
      ns_y= ibm->ns_y[elmt];	
      ns_z= ibm->ns_z[elmt];

      lhs[0][0] = sb;
      lhs[0][1] = sb*sb;
      lhs[0][2] = sb*sb*sb;
      lhs[1][0] = sc;
      lhs[1][1] = sc*sc;
      lhs[1][2] = sc*sc*sc;
      lhs[2][0] = sd;
      lhs[2][1] = sd*sd;
      lhs[2][2] = sd*sd*sd;

      rhs_l[0][1] = lhs[0][1];
      rhs_l[0][2] = lhs[0][2];
      rhs_l[1][1] = lhs[1][1];
      rhs_l[1][2] = lhs[1][2];
      rhs_l[2][1] = lhs[2][1];
      rhs_l[2][2] = lhs[2][2];
      dtm = detmnt(lhs);
      
      phia= cv1*ns_x + cv2*ns_y + cv3*ns_z;
      if (elmtinfo[elmt].Need3rdPoint) {
	phib= uinp2.x*ns_x + uinp2.y*ns_y + uinp2.z*ns_z;
	phic= uinp3.x*ns_x + uinp3.y*ns_y + uinp3.z*ns_z;
	elmtinfo[elmt].Tow_ws = ( -(phib-phia)*sd*sd + (phic-phia)*di*di )
	                        /( di*di* sd - sd*sd * di);

      } else {
	phib= uinp.x*ns_x + uinp.y*ns_y + uinp.z*ns_z;
	phic= uinp2.x*ns_x + uinp2.y*ns_y + uinp2.z*ns_z;
	phid= uinp3.x*nt_x + uinp3.y*nt_y + uinp3.z*nt_z;
	rhs_l[0][0] = phib-phia;
	rhs_l[1][0] = phic-phia;
	rhs_l[2][0] = phid-phia;
/* 	elmtinfo[elmt].Tow_ws = detmnt(rhs_l)/dtm; */
	elmtinfo[elmt].Tow_ws = ( -(phib-phia)*sd*sd + (phic-phia)*di*di )
	                        /( di*di* sd - sd*sd * di);
      }

      nt_x= ibm->nt_x[elmt];	
      nt_y= ibm->nt_y[elmt];
      nt_z= ibm->nt_z[elmt];

      phia= cv1*nt_x + cv2*nt_y + cv3*nt_z;
      if (elmtinfo[elmt].Need3rdPoint) {
	phib= uinp2.x*nt_x + uinp2.y*nt_y + uinp2.z*nt_z;
	phic= uinp3.x*nt_x + uinp3.y*nt_y + uinp3.z*nt_z;
	elmtinfo[elmt].Tow_wt = ( -(phib-phia)*sd*sd + (phic-phia)*di*di )
	                          /( di*di* sd - sd*sd * di);
      } else {
	phib= uinp.x*nt_x + uinp.y*nt_y + uinp.z*nt_z;
	phic= uinp2.x*nt_x + uinp2.y*nt_y + uinp2.z*nt_z;
	phid= uinp3.x*nt_x + uinp3.y*nt_y + uinp3.z*nt_z;
	rhs_l[0][0] = phib-phia;
	rhs_l[1][0] = phic-phia;
	rhs_l[2][0] = phid-phia;
/* 	elmtinfo[elmt].Tow_wt = detmnt(rhs_l)/dtm; */
	elmtinfo[elmt].Tow_wt = ( -(phib-phia)*sd*sd + (phic-phia)*di*di )
	                          /( di*di* sd - sd*sd * di);
      }
      
/*       if (di>1e-10) { */
/*       elmtinfo[elmt].Tow_ws =   ((uinp.x-cv1)*ns_x +  */
/* 				 (uinp.y-cv2)*ns_y + */
/* 				 (uinp.z-cv3)*ns_z)/di; */

/*       elmtinfo[elmt].Tow_wt =   ((uinp.x-cv1)*nt_x +  */
/* 				 (uinp.y-cv2)*nt_y + */
/* 				 (uinp.z-cv3)*nt_z)/di; */
      
/*       } else { */
/* 	PetscPrintf(PETSC_COMM_WORLD, "zero di %le ",di); */
/*       }       */

     if (elmt==1124 || elmt==1110 || elmt==1395 || elmt==1179) 
      PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %d %le %le %le %le %le %le \n",elmt,elmtinfo[elmt].Need3rdPoint,elmtinfo[elmt].Tow_wt,phia,phib,phic,di,sd);
/*       PetscPrintf(PETSC_COMM_SELF, "intp Tow_t %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_wt,nt_x,nt_y,nt_z,di); */
     /*  PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,uinp.x,uinp.y,uinp.z,di); */
/*       PetscPrintf(PETSC_COMM_SELF, "intp Tow_s %d %le %le %le %le %le \n",elmt,elmtinfo[elmt].Tow_ws,ns_x,ns_y,ns_z,di); */

    }
  }
  }

  PetscInt rank;
  PetscInt n_v=ibm->n_v;
  PetscInt ti=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    //if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "Stress%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,tow_t,tow_s,p,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z,di\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-16]=CELLCENTERED)\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
      }
      for (i=0; i<n_v; i++) {	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].Tow_wt);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].Tow_ws);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", elmtinfo[i].P);
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
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibminfo[i].d_i);
      }

      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
      //}
  }

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lP, &p);
  return(0);

}

PetscErrorCode ibm_Surf_stress(UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo, PetscInt ibi)
{
  DA		da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      nbn,nbn2, nbn2D[1000];
  //PetscInt      nbnumber = user->ibmnumber, nbnumber2D;
  PetscInt      i,j,k;
  PetscReal     sb, sc ;//, lhs[3][3], rhs_l[3][3];
  PetscReal     cv1, cv2, cv3;
  Cmpnts	***ucat;
  Cmpnts	***coor;
  PetscReal     ***p;
  // PetscReal cs1, cs2, cs3;
  PetscInt	ni;
  //PetscReal	ucx, ucy, ucz;
  // Added 4/3/06 iman
  PetscInt      n_elmt=ibm->n_elmt;
  PetscInt      n_P[n_elmt];
  PetscReal     Ps[n_elmt], Tow_ws[n_elmt],Tow_wt[n_elmt],n_Psum[n_elmt],n_Pr[n_elmt];
  PetscReal     PsSum[n_elmt], Tow_wsSum[n_elmt],Tow_wtSum[n_elmt];
  PetscReal     nt_x, nt_y, nt_z;
  PetscReal     ns_x, ns_y, ns_z;
  PetscReal     x,y,z,x_c,y_c=6.0,z_c=15.;
  PetscReal     x2,y2,z2;

  IBMInfo       *ibminfo;
  IBMListNode   *current;

/*   PetscInt rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   if (!rank) { */

    //Added 4/27/06 iman
    //reset the values which will be interpolated
    for (i=0;i<ibm->n_elmt;i++) {
      Ps[i]  = 0.;
      n_P[i] = 0;
      Tow_ws[i]=0.; Tow_wt[i]=0.;
    }
    //nbnumber2D=0;

    DAVecGetArray(fda, user->lUcat, &ucat);
    DAVecGetArray(da, user->lP, &p);
    
    
    current = user->ibmlist[ibi].head;
    while (current) {
/*     for (nbn=0; nbn<nbnumber; nbn++) { */
      ibminfo = &current->ibm_intp;
      current = current->next;	
      //i=ibminfo[nbn].ni; j=ibminfo[nbn].nj; k=ibminfo[nbn].nk;
/*       i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk; */
/*       sb = ibminfo[nbn].d_s; sc = sb + ibminfo[nbn].d_i; */
/*       ni = ibminfo[nbn].cell; */

      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;      
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;
      ni = ibminfo->cell;

/*       // get all the adj pts in an i-plane */
/*       if (i==1) { */
/* 	nbnumber2D ++; */
/* 	nbn2D[nbnumber2D-1]=nbn; */
/*       } */

      //PetscPrintf(PETSC_COMM_WORLD, "IBM_intp sk: %le %le %le %le %le %le\n",sk1,sk2,sk3,cv1,cv2,cv3);
    
      //Added 4/3/06 iman
      if (ni>=0 && i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze) {	       
	n_P[ni]++;
	n_Pr[ni]=(double)(n_P[ni]);
	//Ps=ibm->P[ni];      
	Ps[ni]+=( p[k][j][i]);
	     //+ (n_P-1)*Ps )/n_P;      
/* 	- user->st/user->dt* */
/* 	           (ibm->nf_z[ni]*(ibm->u[ni].z-ibm->u_o[ni].z) */
/* 		   +ibm->nf_y[ni]*(ibm->u[ni].y-ibm->u_o[ni].y) */
/* 		   +ibm->nf_x[ni]*(ibm->u[ni].x-ibm->u_o[ni].x)) */
	  //ibm->P[ni]=Ps;
	  //ibm->n_P[ni]=n_P;
/* 	if (i==4) { */
/* 	  PetscPrintf(PETSC_COMM_WORLD, "Pijk=%le, Pelm=%le, n_P=%d, %d %d %d %le\n",p[k][j][i],ibm->P[ni],ibm->n_P[ni],k,j,i,ibm->nf_x[ni]); */
/* 	} */
	
	
	cv1 = ( ibm->u[ibm->nv1[ni]].x
		+ibm->u[ibm->nv2[ni]].x 
		+ibm->u[ibm->nv3[ni]].x)/3.;
	cv2 = ( ibm->u[ibm->nv1[ni]].y
		+ibm->u[ibm->nv2[ni]].y
		+ibm->u[ibm->nv3[ni]].y)/3.;
	cv3 = ( ibm->u[ibm->nv1[ni]].z
		+ibm->u[ibm->nv2[ni]].z
		+ibm->u[ibm->nv3[ni]].z)/3.;
	
	ns_x= ibm->ns_x[ni];	
	ns_y= ibm->ns_y[ni];	
	ns_z= ibm->ns_z[ni];
	
	nt_x= ibm->nt_x[ni];	
	nt_y= ibm->nt_y[ni];
	nt_z= ibm->nt_z[ni];
	
	if (fabs(sb)>1e-10) {
	Tow_ws[ni] +=((ucat[k][j][i].x-cv1)*ns_x + (ucat[k][j][i].y-cv2)*ns_y +
		       (ucat[k][j][i].z-cv3)*ns_z)/sb;  
	//ibm->Tow_ws[ni]=(ibm->Tow_ws[ni]*(n_P-1)+Tow_ws)/n_P;
	
	Tow_wt[ni] +=((ucat[k][j][i].x-cv1)*nt_x + (ucat[k][j][i].y-cv2)*nt_y +
		(ucat[k][j][i].z-cv3)*nt_z)/sb;  
	//ibm->Tow_wt[ni]=(ibm->Tow_wt[ni]*(n_P-1)+Tow_wt)/n_P;

	//MPI_Bcast(&(ibm->n_[elmt]), 1, MPI_INTEGER, 0, PETSC_COMM_WORLD);

	} else {
	  PetscPrintf(PETSC_COMM_WORLD, " sb < 1e-10 !!!!!!!!!!!!!!!\n");
	}
      } else {
	
      }
      
    }

    PetscBarrier(PETSC_NULL);
    
    for (nbn=0; nbn<n_elmt; nbn++) {
      PetscGlobalSum(&n_Pr[nbn], &n_Psum[nbn], PETSC_COMM_WORLD);
      PetscGlobalSum(&Ps[nbn], &PsSum[nbn], PETSC_COMM_WORLD);
      PetscGlobalSum(&Tow_ws[nbn], &Tow_wsSum[nbn], PETSC_COMM_WORLD);
      PetscGlobalSum(&Tow_wt[nbn], &Tow_wtSum[nbn], PETSC_COMM_WORLD);
    }
    
    for (nbn=0; nbn<n_elmt; nbn++) {
      elmtinfo[nbn].n_P=(int)(n_Psum[nbn]+.1);
      if (n_Psum[nbn]>1e-6){
	elmtinfo[nbn].P=PsSum[nbn]/n_Psum[nbn];
	elmtinfo[nbn].Tow_ws=Tow_wsSum[nbn]/n_Psum[nbn];
	elmtinfo[nbn].Tow_wt=Tow_wtSum[nbn]/n_Psum[nbn];
	elmtinfo[nbn].FoundAroundcell = 1;
	//PetscPrintf(PETSC_COMM_WORLD, " n_P %d %le %le %le\n", n_Psum[nbn],PsSum[nbn],Tow_ws[nbn],Tow_wt[nbn]);
      } else {
	elmtinfo[nbn].P=0.;
	elmtinfo[nbn].Tow_ws=0.;
	elmtinfo[nbn].Tow_wt=0.;
      }
    }

    DAVecGetArray(fda, user->Cent, &coor);

    
/*     PetscInt  sort1[1000],sort2[1000],sort[2000],sort1no=0,sort2no=0; */
/*     PetscInt  nbn3; */
/*     Cmpnts    nf[1000],coor_pt[1000]; */
/*     PetscReal r,p_pt[1000]; */
/*     // sort nbn2D into two parts */
/*     for (nbn2=0; nbn2<nbnumber2D; nbn2++) { */
/*       nbn=nbn2D[nbn2]; */
/*       i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk;       */

/*       if (i>=xs && i<=xe && j>=ys && j<=ye && k>=zs && k<=ze) {	 */
/*       x=coor[k][j][i].x; y=coor[k][j][i].y; z=coor[k][j][i].z; */

/*       // calculate the normal(nf) and coord of pt on cyl(coor_pt) */
/*       nf[nbn2].x= 0; */
/*       nf[nbn2].y= y-y_c; */
/*       nf[nbn2].z= z-z_c; */
/*       r = sqrt(nf[nbn2].y*nf[nbn2].y+nf[nbn2].z*nf[nbn2].z); */
/*       nf[nbn2].y= nf[nbn2].y/r ; */
/*       nf[nbn2].z= nf[nbn2].z/r ; */
      
/*       coor_pt[nbn2].x=x;  */
/*       coor_pt[nbn2].y=y_c+(y-y_c)/(2.*r); */
/*       coor_pt[nbn2].z=z_c+(z-z_c)/(2.*r); */

/*       p_pt[nbn2]=p[k][j][i]; */
/*       if (coor_pt[nbn2].z < z_c) { */
/* 	sort1no++; */
/* 	sort1[sort1no-1]=nbn2; */
/*       } else { */
/* 	sort2no++; */
/* 	sort2[sort2no-1]=nbn2; */
/*       } */
/*       } */
/*     } */

/*     PetscBarrier(PETSC_NULL); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "%d %d  sortno\n",sort1no,sort2no); */

/*     // sort the sort1 based on increasing y */
/*     for (nbn2=0; nbn2<sort1no; nbn2++) { */
/*       for (nbn3= nbn2+1; nbn3< sort1no; nbn3++) { */
/* 	nbn=sort1[nbn2]; */
/* 	//i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk;       */
/* 	x=coor_pt[nbn].x; y=coor_pt[nbn].y; z=coor_pt[nbn].z; */
	
/* 	nbn=sort1[nbn3]; */
/* 	x2=coor_pt[nbn].x; y2=coor_pt[nbn].y; z2=coor_pt[nbn].z; */

/* 	//i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk;       */
/* 	//x2=coor[k][j][i].x; y2=coor[k][j][i].y; z2=coor[k][j][i].z; */
	
/* 	if (y>y2) { */
/* 	  nbn=sort1[nbn2]; */
/* 	  sort1[nbn2]=sort1[nbn3]; */
/* 	  sort1[nbn3]=nbn; */
/* 	} */
/*       } */
/*     } */
/*     //PetscPrintf(PETSC_COMM_WORLD, " nbnumber2D==nbn3\n"); */
      
/*     // sort the sort2 based on decreasing y */
/*     for (nbn2=0; nbn2<sort2no; nbn2++) { */
/*       for (nbn3= nbn2+1; nbn3< sort2no; nbn3++) { */
/* 	nbn=sort2[nbn2]; */
/* 	x=coor_pt[nbn].x; y=coor_pt[nbn].y; z=coor_pt[nbn].z; */

/* 	//i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk;       */
/* 	//x=coor[k][j][i].x; y=coor[k][j][i].y; z=coor[k][j][i].z; */
	
/* 	nbn=sort2[nbn3]; */
/* 	x2=coor_pt[nbn].x; y2=coor_pt[nbn].y; z2=coor_pt[nbn].z; */
/* 	//i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk;       */
/* 	//x2=coor[k][j][i].x; y2=coor[k][j][i].y; z2=coor[k][j][i].z; */
	
/* 	if (y<y2) { */
/* 	  nbn=sort2[nbn2]; */
/* 	  sort2[nbn2]=sort2[nbn3]; */
/* 	  sort2[nbn3]=nbn; */
/* 	} */
/*       } */
/*     } */
/*     //PetscPrintf(PETSC_COMM_WORLD, " nbnumber2D==nbn3\n"); */

/*     // assembel sort1 and sort2 into sort */
/*     nbn3=-1; */
/*     for (nbn2=0; nbn2<sort1no; nbn2++) { */
/*       nbn3++; */
/*       sort[nbn3]=sort1[nbn2]; */
/*     } */
/*     for (nbn2=0; nbn2<sort2no; nbn2++) { */
/*       nbn3++; */
/*       sort[nbn3]=sort2[nbn2]; */
/*     } */

/*     PetscReal P_y,P_z,dl,F_z=0.,F_y=0.; */
/*     PetscPrintf(PETSC_COMM_SELF, "%d %d nbnumber2D==nbn3\n", nbnumber2D, nbn3); */
/*     for (nbn3=0; nbn3<nbnumber2D-1; nbn3++) { */
/*       nbn2=sort[nbn3]; */
/*       nbn=nbn2D[nbn2]; */
/*       i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk; */
/*       P_z=p_pt[nbn2]*nf[nbn2].z; */
/*       P_y=p_pt[nbn2]*nf[nbn2].y; */
/*       x=coor_pt[nbn2].x; y=coor_pt[nbn2].y; z=coor_pt[nbn2].z; */

/*       nbn2=sort[nbn3+1]; */
/*       nbn=nbn2D[nbn2]; */
/*       i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk; */
/*       P_z+=p_pt[nbn2]*nf[nbn2].z; */
/*       P_z=P_z/2.; */
/*       P_y+=p_pt[nbn2]*nf[nbn2].y; */
/*       P_y=P_y/2.; */
/*       x2=coor_pt[nbn2].x; y2=coor_pt[nbn2].y; z2=coor_pt[nbn2].z; */

/*       dl=sqrt((y-y2)*(y-y2)+(z-z2)*(z-z2)); */
      
/*       F_z +=P_z*dl; */
/*       F_y +=P_y*dl; */
/*       //tscPrintf(PETSC_COMM_WORLD, "%d %d %le %le %le %le %le %le %le %le\n", j,k,z-z_c,y,y2,dl,nf[nbn2].z,nf[nbn2].y,P_z,P_y); */
/*     } */
/*     if (nbnumber2D>0) { */
/*       // start of the chain */
/*       nbn3=0; */
/*       nbn2=sort[nbn3]; */
/*       nbn=nbn2D[nbn2];       */
/*       P_z=p_pt[nbn2]*nf[nbn2].z; */
/*       P_y=p_pt[nbn2]*nf[nbn2].y; */
/*       x=coor_pt[nbn2].x; y=coor_pt[nbn2].y; z=coor_pt[nbn2].z; */

/*       nbn2=sort[nbn3+1]; */
/*       nbn=nbn2D[nbn2]; */
/*       x2=coor_pt[nbn2].x; y2=coor_pt[nbn2].y; z2=coor_pt[nbn2].z; */

/*       dl=sqrt((y-y2)*(y-y2)+(z-z2)*(z-z2)); */
      
/*       F_z +=P_z*dl/2.; */
/*       F_y +=P_y*dl/2.; */

/*       // end of chain */
/*       nbn3=nbnumber2D-2; */
/*       nbn2=sort[nbn3];       */
/*       x=coor_pt[nbn2].x; y=coor_pt[nbn2].y; z=coor_pt[nbn2].z; */

/*       nbn2=sort[nbn3+1]; */
/*       nbn=nbn2D[nbn2]; */
/*       P_z+=p_pt[nbn2]*nf[nbn2].z;       */
/*       P_y+=p_pt[nbn2]*nf[nbn2].y; */
/*       x2=coor_pt[nbn2].x; y2=coor_pt[nbn2].y; z2=coor_pt[nbn2].z; */

/*       dl=sqrt((y-y2)*(y-y2)+(z-z2)*(z-z2)); */
      
/*       F_z +=P_z*dl/2.; */
/*       F_y +=P_y*dl/2.; */
        
/* /\*     nbn2=sort[nbnumber2D-1]; *\/ */
/* /\*     nbn3=sort[0]; *\/ */
/* /\*     nbn=nbn2D[nbn2]; *\/ */
/* /\*     i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk; *\/ */
/* /\*     P_z=p_pt[nbn2]*nf[nbn2].z; *\/ */
/* /\*     P_y=p_pt[nbn2]*nf[nbn2].y; *\/ */
    
/* /\*     nbn=nbn2D[nbn3]; *\/ */
/* /\*     i = ibminfo[nbn].ni; j= ibminfo[nbn].nj; k = ibminfo[nbn].nk; *\/ */
/* /\*     P_z+=p_pt[nbn3]*nf[nbn3].z; *\/ */
/* /\*     P_z=P_z/2.; *\/ */
/* /\*     P_y+=p_pt[nbn3]*nf[nbn3].y; *\/ */
/* /\*     P_y=P_y/2.; *\/ */
    
/* /\*     x=coor_pt[nbn2].x; y=coor_pt[nbn2].y; z=coor_pt[nbn2].z; *\/ */
/* /\*     x2=coor_pt[nbn3].x; y2=coor_pt[nbn3].y; z2=coor_pt[nbn3].z; *\/ */
/* /\*     dl=sqrt((y-y2)*(y-y2)+(z-z2)*(z-z2)); *\/ */
    
/* /\*     F_z +=P_z*dl; *\/ */
/* /\*     F_y +=P_y*dl; *\/ */

/*     F_z=2*F_z; */
/*     F_y=2*F_y; */
/*     } */
/*     PetscReal F_ztotal, F_ytotal; */
/*     PetscGlobalSum(&F_z , &F_ztotal, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&F_y , &F_ytotal, PETSC_COMM_WORLD); */

/*     //PetscPrintf(PETSC_COMM_WORLD, "%d %le %le %le %le %le %le\n", nbn2,z-z_c,y,y2,dl,P_z,P_y); */
/*     PetscPrintf(PETSC_COMM_WORLD, "F_z, F_y 2D %d %le %le %le %le\n", nbn2,F_z,F_y,F_ztotal,F_ytotal); */
/*     PetscInt rank; */
/*     MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*     if (!rank) { */
/*       FILE *f; */
/*       char filen[80]; */
/*       sprintf(filen, "Force_Coeff2D"); */
/*       f = fopen(filen, "a"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le\n",ti,F_y,F_z, F_ytotal,F_ztotal); */
/*       fclose(f); */
/*     } */
    
    DAVecRestoreArray(fda, user->Cent,&coor);
    DAVecRestoreArray(fda, user->lUcat, &ucat);
    DAVecRestoreArray(da, user->lP, &p);

    return(0);
}


PetscErrorCode Calc_forces(FSInfo *FSinfo, IBMNodes *ibm, SurfElmtInfo *elmtinfo,PetscReal Re,PetscInt ti)   
{
  PetscReal      pi=3.141592654;
  PetscInt       i,n_elmt,elmt;
  PetscReal     *dA ;         // area of an element
  PetscReal     *P;    //Press on the surface elmt
  PetscInt      *n_P; //number of Press Pts on the elmt
  PetscReal     *Tow_ws, *Tow_wt, *Tow_wn; //wall shear stress of the elmt
  PetscReal     *nf_x, *nf_y, *nf_z; //normal dir
  PetscReal     *nt_x, *nt_y, *nt_z; //tangent dir
  PetscReal     *ns_x, *ns_y, *ns_z; //azimuthal dir
  PetscReal      F_x,F_y,F_z, A_tot; //Forces and Area
  PetscReal      Cp_x,Cp_y,Cp_z; //Pressure Forces
  PetscReal      Cs_x,Cs_y,Cs_z; //Surface Forces
  PetscReal      F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force
  PetscReal      Cp_xSum,Cp_ySum,Cp_zSum; //Pressure Force

  n_elmt=ibm->n_elmt;  
  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force  %le %d %d\n",Re, n_elmt, ibm->n_elmt);

  // Allocate memory
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z); */
    
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area */
  PetscMalloc(n_elmt*sizeof(PetscReal), &(P));  //Press
  PetscMalloc(n_elmt*sizeof(PetscReal), &(n_P)); //no. of Press pts
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_ws)); //wall shear stress
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_wt)); //wall shear stress
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_wn)); //wall shear stress

/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z); */
  
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z); */

  // Init values
  nf_x=ibm->nf_x; nf_y=ibm->nf_y; nf_z=ibm->nf_z;
  nt_x=ibm->nt_x; nt_y=ibm->nt_y; nt_z=ibm->nt_z;
  ns_x=ibm->ns_x; ns_y=ibm->ns_y; ns_z=ibm->ns_z;
  dA  =ibm->dA;
  
  for (elmt=0; elmt<n_elmt; elmt++) {
    n_P[elmt]   =elmtinfo[elmt].n_P; 
    P[elmt]     =elmtinfo[elmt].P;
    Tow_ws[elmt]=elmtinfo[elmt].Tow_ws; 
    Tow_wt[elmt]=elmtinfo[elmt].Tow_wt;
    Tow_wn[elmt]=elmtinfo[elmt].Tow_wn;
  }

  // Calc forces
  F_x=0.;F_y=0.;F_z=0.;A_tot=0.;
  Cp_x=0.;Cp_y=0.;Cp_z=0.;
  Cs_x=0.;Cs_y=0.;Cs_z=0.;
  for (i=0; i<n_elmt; i++) {
    if (n_P[i]>0  && elmtinfo[i].FoundAroundcell>0) {
      Cp_x+=(-P[i]*nf_x[i])*dA[i]; 
      Cp_y+=(-P[i]*nf_y[i])*dA[i]; 
      Cp_z+=(-P[i]*nf_z[i])*dA[i]; 

      Cs_x+=(Tow_wn[i]*nf_x[i] + Tow_ws[i]*ns_x[i] + Tow_wt[i]*nt_x[i])/Re*dA[i]; 
      Cs_y+=(Tow_wn[i]*nf_y[i] + Tow_ws[i]*ns_y[i] + Tow_wt[i]*nt_y[i])/Re*dA[i]; 
      Cs_z+=(Tow_wn[i]*nf_z[i] + Tow_ws[i]*ns_z[i] + Tow_wt[i]*nt_z[i])/Re*dA[i]; 

      //F_x+=(-P[i]*nf_x[i] + Tow_ws[i]*ns_x[i]/Re + Tow_wt[i]*nt_x[i]/Re)*dA[i]; 
      //F_y+=(-P[i]*nf_y[i] + Tow_ws[i]*ns_y[i]/Re + Tow_wt[i]*nt_y[i]/Re)*dA[i]; 
      //F_z+=(-P[i]*nf_z[i] + Tow_ws[i]*ns_z[i]/Re + Tow_wt[i]*nt_z[i]/Re)*dA[i];      
      A_tot +=dA[i];
    }
  }

  // Total Forces on each processor
  F_x=Cp_x + Cs_x; 
  F_y=Cp_y + Cs_y;
  F_z=Cp_z + Cs_z;
  
  // Global Sum
  PetscGlobalSum(&F_x, &F_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_y, &F_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_z, &F_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&A_tot, &A_totSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Cp_x, &Cp_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_y, &Cp_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_z, &Cp_zSum, PETSC_COMM_WORLD);

  // Scale Check later !!!!!
  F_xSum=F_xSum/A_totSum*pi*2.;
  F_ySum=F_ySum/A_totSum*pi*2.;
  F_zSum=F_zSum/A_totSum*pi*2.;

  Cp_xSum=Cp_xSum/A_totSum*pi*2.;
  Cp_ySum=Cp_ySum/A_totSum*pi*2.;
  Cp_zSum=Cp_zSum/A_totSum*pi*2.;    

  // store results in FSinfo
  FSinfo->F_x = F_xSum; FSinfo->F_y = F_ySum; FSinfo->F_z = F_zSum;
  FSinfo->A_tot = A_totSum;

  // output values
  PetscPrintf(PETSC_COMM_SELF, "F_x,F_y,F_z, %le %le %le %le\n",F_x,F_y,F_z,A_tot);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum,Cs_x,Cs_y,Cs_z, A_totSum);
    fclose(f);
  }

  // free memory
  PetscFree(P);
  PetscFree(n_P);
  PetscFree(Tow_ws);
  PetscFree(Tow_wt);
  PetscFree(Tow_wn);

  return(0);
}

PetscErrorCode Calc_forces2(FSInfo *FSinfo, IBMNodes *ibm, SurfElmtInfo *elmtinfo,PetscReal Re,PetscInt ti)   
{
  PetscReal      pi=3.141592654;
  PetscInt       i,n_elmt,elmt;
  PetscReal     *dA ;         // area of an element
  PetscReal     *P;    //Press on the surface elmt
  PetscInt      *n_P; //number of Press Pts on the elmt
  PetscReal     *Tow_ws, *Tow_wt; //wall shear stress of the elmt
  PetscReal     *nf_x, *nf_y, *nf_z; //normal dir
  PetscReal     *nt_x, *nt_y, *nt_z; //tangent dir
  PetscReal     *ns_x, *ns_y, *ns_z; //azimuthal dir
  PetscReal      F_x,F_y,F_z, A_tot; //Forces and Area
  PetscReal      Cp_x,Cp_y,Cp_z; //Pressure Forces
  PetscReal      Cs_x,Cs_y,Cs_z; //Surface Forces
  PetscReal      F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force

  n_elmt=ibm->n_elmt;  
  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force  %le %d %d\n",Re, n_elmt, ibm->n_elmt);

  // Allocate memory
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
    
  PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area
  PetscMalloc(n_elmt*sizeof(PetscReal), &(P));  //Press
  PetscMalloc(n_elmt*sizeof(PetscReal), &(n_P)); //no. of Press pts
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_ws)); //wall shear stress
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_wt)); //wall shear stress

  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

  // Init values
  nf_x=ibm->nf_x; nf_y=ibm->nf_y; nf_z=ibm->nf_z;
  nt_x=ibm->nt_x; nt_y=ibm->nt_y; nt_z=ibm->nt_z;
  ns_x=ibm->ns_x; ns_y=ibm->ns_y; ns_z=ibm->ns_z;
  dA  =ibm->dA;
  
  for (elmt=0; elmt<n_elmt; elmt++) {
  n_P[elmt]   =elmtinfo[elmt].n_P; 
  P[elmt]     =elmtinfo[elmt].P;
  Tow_ws[elmt]=elmtinfo[elmt].Tow_ws; 
  Tow_wt[elmt]=elmtinfo[elmt].Tow_wt;
  }

  // Calc forces
  F_x=0.;F_y=0.;F_z=0.;A_tot=0.;
  Cp_x=0.;Cp_y=0.;Cp_z=0.;
  Cs_x=0.;Cs_y=0.;Cs_z=0.;
  for (i=0; i<n_elmt; i++) {
    if (n_P[i]>0 ) {
      Cp_x+=(-P[i]*nf_x[i])*dA[i]; 
      Cp_y+=(-P[i]*nf_y[i])*dA[i]; 
      Cp_z+=(-P[i]*nf_z[i])*dA[i]; 

      Cs_x+=(Tow_ws[i]*ns_x[i]/Re + Tow_wt[i]*nt_x[i]/Re)*dA[i]; 
      Cs_y+=(Tow_ws[i]*ns_y[i]/Re + Tow_wt[i]*nt_y[i]/Re)*dA[i]; 
      Cs_z+=(Tow_ws[i]*ns_z[i]/Re + Tow_wt[i]*nt_z[i]/Re)*dA[i]; 

      //F_x+=(-P[i]*nf_x[i] + Tow_ws[i]*ns_x[i]/Re + Tow_wt[i]*nt_x[i]/Re)*dA[i]; 
      //F_y+=(-P[i]*nf_y[i] + Tow_ws[i]*ns_y[i]/Re + Tow_wt[i]*nt_y[i]/Re)*dA[i]; 
      //F_z+=(-P[i]*nf_z[i] + Tow_ws[i]*ns_z[i]/Re + Tow_wt[i]*nt_z[i]/Re)*dA[i];      
      A_tot +=dA[i];
    }
  }

  // Total Forces on each processor
  F_x=Cp_x + Cs_x; 
  F_y=Cp_y + Cs_y;
  F_z=Cp_z + Cs_z;
  
  // Global Sum
/*   PetscGlobalSum(&F_x, &F_xSum, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&F_y, &F_ySum, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&F_z, &F_zSum, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&A_tot, &A_totSum, PETSC_COMM_WORLD); */
  
/*   // Scale Check later !!!!! */
/*   F_x=F_x/A_tot*pi*2.; */
/*   F_y=F_y/A_tot*pi*2.; */
/*   F_z=F_z/A_tot*pi*2.; */
  

  // store results in FSinfo
  //FSinfo->F_x = F_xSum; FSinfo->F_y = F_ySum; FSinfo->F_z = F_zSum;
  FSinfo->A_tot = A_totSum;

  // output values
  PetscPrintf(PETSC_COMM_SELF, "F_x,F_y,F_z 2, %le %le %le %le\n",F_x,F_y,F_z,A_tot);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff2");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti, F_x, F_y, F_z,Cp_x,Cp_y,Cp_z,Cs_x,Cs_y,Cs_z, A_tot);
    fclose(f);
  }
  return(0);
}

PetscErrorCode Calc_Moments(FSInfo *FSinfo, IBMNodes *ibm, SurfElmtInfo *elmtinfo, PetscReal Re, PetscInt ti)   
{
  PetscReal      pi=3.141592654;
  PetscInt       elmt,n_elmt,n_v;
  PetscInt	*nv1, *nv2, *nv3;	// Node index of each triangle
  PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
  PetscReal     *dA ;         // area of an element
  PetscReal     *P;    //Press on the surface elmt
  PetscInt      *n_P; //number of Press Pts on the elmt
  PetscReal     *Tow_ws, *Tow_wt; //wall shear stress of the elmt
  PetscReal     *nt_x, *nt_y, *nt_z; //tangent dir
  PetscReal     *ns_x, *ns_y, *ns_z; //azimuthal dir
  PetscReal      F_x,F_y,F_z, A_tot; //Forces and Area
  PetscReal      M_x,M_y,M_z;   //Moments
  PetscReal      r_x,r_y,r_z;   //Anchor dist
  PetscReal      x, y, z;       //Mid pt of the elmt
  PetscReal      X_c,Y_c,Z_c;   //center of rotation coord
  PetscReal      M_xSum,M_ySum,M_zSum,A_totSum; //Surface Mom on all processors

  n_elmt=ibm->n_elmt; n_v=ibm->n_v;
  // Allocate memory
  PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
  PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
  PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

  PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
    
  PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area
  PetscMalloc(n_elmt*sizeof(PetscReal), &(P));  //Press
  PetscMalloc(n_elmt*sizeof(PetscReal), &(n_P)); //no. of Press pts
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_ws)); //wall shear stress
  PetscMalloc(n_elmt*sizeof(PetscReal), &(Tow_wt)); //wall shear stress

  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);
    
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

  // Init values
 /*  x_bp=ibm->x_bp; y_bp=ibm->y_bp; z_bp=ibm->z_bp; */
  nv1 =ibm->nv1 ; nv2 =ibm->nv2 ; nv3 =ibm->nv3 ;
  nf_x=ibm->nf_x; nf_y=ibm->nf_y; nf_z=ibm->nf_z;
  nt_x=ibm->nt_x; nt_y=ibm->nt_y; nt_z=ibm->nt_z;
  ns_x=ibm->ns_x; ns_y=ibm->ns_y; ns_z=ibm->ns_z;
  dA  =ibm->dA;
  
  for (elmt=0; elmt<n_elmt; elmt++) {
    n_P[elmt]   =elmtinfo[elmt].n_P; 
    P[elmt]     =elmtinfo[elmt].P;
    Tow_ws[elmt]=elmtinfo[elmt].Tow_ws; 
    Tow_wt[elmt]=elmtinfo[elmt].Tow_wt;
  }

  X_c=FSinfo->x_c; Y_c=FSinfo->y_c; Z_c=FSinfo->z_c;

  // Calc Moments Check later for /Re scaling !!!!
  M_x=0.;M_y=0.;M_z=0.;A_tot=0.;
  for (elmt=0; elmt<n_elmt; elmt++) {
    if (n_P[elmt]>0 && elmtinfo[elmt].FoundAroundcell>0) {
      x=ibm->cent_x[elmt]; 
      y=ibm->cent_y[elmt]; 
      z=ibm->cent_z[elmt]; 
      
      r_x = x-X_c;
      r_y = y-Y_c;
      r_z = z-Z_c;
      
      F_x=(-P[elmt]*nf_x[elmt] + Tow_ws[elmt]*ns_x[elmt]/Re + Tow_wt[elmt]*nt_x[elmt]/Re)*dA[elmt]; 
      F_y=(-P[elmt]*nf_y[elmt] + Tow_ws[elmt]*ns_y[elmt]/Re + Tow_wt[elmt]*nt_y[elmt]/Re)*dA[elmt]; 
      F_z=(-P[elmt]*nf_z[elmt] + Tow_ws[elmt]*ns_z[elmt]/Re + Tow_wt[elmt]*nt_z[elmt]/Re)*dA[elmt]; 
      
      M_x +=   r_y*F_z - r_z*F_y;
      M_y += -(r_x*F_z - r_z*F_x);
      M_z +=   r_x*F_y - r_y*F_x;
      
      A_tot +=dA[elmt];
    }
  }

  // scale check later for consistancy!!!!
/*   M_x=M_x/A_tot*2; */
/*   M_y=M_y/A_tot*2; */
/*   M_z=M_z/A_tot*2; */

  // Global Sum
  PetscGlobalSum(&M_x, &M_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_y, &M_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_z, &M_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&A_tot, &A_totSum, PETSC_COMM_WORLD);

  // store results in FSinfo
  FSinfo->M_x = M_xSum; FSinfo->M_y = M_ySum; FSinfo->M_z = M_zSum;
  FSinfo->A_tot = A_totSum;

  // output values
  PetscPrintf(PETSC_COMM_SELF, "M_x,M_y,M_z, %le %le %le %le\n",M_x,M_y,M_z,A_tot);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Moment_Coeff");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum, A_totSum);
    fclose(f);
  }
  return(0);

}

PetscErrorCode Calc_forces_CVM2D(UserCtx *user, FSInfo *fsi, PetscInt ti) 
{ 

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  // cell location 
  PetscInt      zin,zout,yin,yout,xin,xout;
  xin=0; xout=2;
  yin=65; yout=130;
  zin=22; zout=98;

  // center cell location
  PetscInt      zstart,zend,ystart,yend,xstart,xend;
  zstart=zin+1;
  zend  =zout+1;
  ystart=yin+1;
  yend  =yout+1;
  xstart=xin+1;
  xend  =xout+1;

  if (zstart<lzs) zstart=lzs;
  if (ystart<lys) ystart=lys;
  if (xstart<lxs) xstart=lxs;
  if (zend>lze) zend=lze;
  if (yend>lye) yend=lye;
  if (xend>lxe) xend=lxe;

  PetscPrintf(PETSC_COMM_SELF, "CV %d %d %d %d\n",xstart,xend, xin, xout);
  
  Vec           Coor;
  PetscInt	i, j, k;
  Cmpnts        ***coor, ***ucat, ***ucont;
  PetscReal     ***p;
  PetscReal     dx,dy,dz;
  PetscReal     F_x,F_y,F_z; //Forces and Area
  PetscReal     F_xSum,F_ySum,F_zSum; //Surface Force
  PetscReal     flux_z,Fp_z, Ft_z;
  PetscReal     flux_y,Fp_y, Ft_y;
  PetscReal     Re=user->ren;
  PetscReal     nj_in,nj_out,nk_in,nk_out;

  // Get Working arrays
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);

  //DAVecGetArray(fda, user->lCent,&coor);
  DAVecGetArray(fda, user->lUcat,  &ucat);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(da, user->lP, &p);

  // calculate normals of the CV
  i=xstart;
  k=zstart;
  j=yin;
  nj_in=coor[k][j][i].y-coor[k][j-1][i].y;
  nj_in=-nj_in/fabs(nj_in);
  j=yout;
  nj_out=coor[k][j+1][i].y-coor[k][j][i].y;
  nj_out=nj_out/fabs(nj_out);

  i=xstart;
  j=ystart;
  k=zin;
  nk_in=coor[k][j][i].z-coor[k-1][j][i].z;
  nk_in=-nk_in/fabs(nk_in);
  k=zout;
  nk_out=coor[k+1][j][i].z-coor[k][j][i].z;
  nk_out=nk_out/fabs(nk_out);
 
  PetscPrintf(PETSC_COMM_SELF, "CV normal %le %le %le %le\n",nk_in,nk_out,nj_in,nj_out);
 
  
  //                    z momentum 
  flux_z=0.;
  Fp_z=0.;
  Ft_z=0.;
  // 1) inflow forces
  if (zin>zs) {
    k=zin; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	dy= fabs(coor[k][j][i].y- coor[k][j-1][i].y);
	dx= fabs(coor[k][j][i].x- coor[k][j][i-1].x);
	flux_z -= ucont[k][j][i].z*
	  (ucat[k][j][i].z+ucat[k+1][j][i].z)*0.5;
	//dy      = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y -
	//	coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	Fp_z   -= 0.5*(p[k][j][i]+p[k+1][j][i])*dy*dx;
      }
    }
  }

  // 2) outflow forces
  if (zout<lze) {
    k=zout; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	dy= coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

	flux_z += ucont[k][j][i].z*
	  (ucat[k][j][i].z+ucat[k+1][j][i].z)*0.5;
	//dy      = 0.25*(coor[k][j+1][i].y+coor[k+1][j+1][i].y -
	//	coor[k][j-1][i].y-coor[k+1][j-1][i].y);
	Fp_z   += 0.5*(p[k][j][i]+p[k+1][j][i])*dy*dx;
      }
    }
  }

  // 3) upperside forces
  if (yout<lye) {
    j=yout; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	//dz     = 0.25*(coor[k+1][j][i].z+coor[k+1][j+1][i].z -
	//       coor[k-1][j][i].z-coor[k-1][j+1][i].z);
	//dy     =       coor[k][j+1][i].y-coor[k][j][i].y;
	flux_z += ucont[k][j][i].y *
	   (ucat[k][j][i].z+ucat[k][j+1][i].z)*0.5;
	dz     = coor[k][j][i].z-coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	// du_z/dy
	Ft_z  -= (ucat[k][j+1][i].z-ucat[k][j][i].z)/Re/dy*dz*dx;
      }
    }
  }

  // 4) lowerside forces
  if (yin>ys) {
    j=yin; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	//dz     = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
	//       coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	//dy     =       coor[k][j][i].y-coor[k][j-1][i].y;
	// du_z/dy
	flux_z -= ucont[k][j][i].y *
	   (ucat[k][j][i].z+ucat[k][j+1][i].z)*0.5;
	dz     = coor[k][j][i].z-coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	// du_z/dy
	Ft_z  += (ucat[k][j+1][i].z-ucat[k][j][i].z)/Re/dy*dz*dx;
      }
    }
  }

  //                    y momentum 
  flux_y=0.;
  Fp_y=0.;
  Ft_y=0.;
  // 1) inflow forces
  if (yin>ys) {
    j=yin; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	dz     = coor[k][j][i].z - coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	flux_y -= ucont[k][j][i].y*
	  (ucat[k][j][i].y+ucat[k][j+1][i].y)*0.5;
	//dz      = 0.25*(coor[k+1][j][i].z+coor[k-1][j-1][i].z -
	//	coor[k+1][j][i].z-coor[k-1][j-1][i].z);
	Fp_y   -= 0.5*(p[k][j][i]+p[k][j+1][i])*dz*dx;
      }
    }
  }

  // 2) outflow forces
  if (yout<lye) {
    j=yout; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	dz      = coor[k][j][i].z - coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	flux_y += ucont[k][j][i].y*
	  (ucat[k][j][i].y+ucat[k][j+1][i].y)*0.5;
	//dz      = 0.25*(coor[k+1][j][i].z+coor[k-1][j+1][i].z -
	//	coor[k+1][j][i].z-coor[k-1][j+1][i].z);

	Fp_y   += 0.5*(p[k][j][i]+p[k][j+1][i])*dz*dx;
      }
    }
  }

  // 3) upperside forces
  if (zout<lze) {
    k=zout; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	//dy     = 0.25*(coor[k][j+1][i].y+coor[k+1][j+1][i].y -
	//       coor[k][j-1][i].y-coor[k+1][j-1][i].y);
	//dz     =       coor[k+1][j][i].z-coor[k][j][i].z;
	flux_y += ucont[k][j][i].z*
	  (ucat[k][j][i].y+ucat[k+1][j][i].y)*0.5;
	dy   = coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	// du_y/dz
	Ft_y  -= (ucat[k+1][j][i].y-ucat[k][j][i].y)/Re/dz*dy*dx;
      }
    }
  }

  // 4) lowerside forces
  if (zin>zs) {
    k=zin; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	//dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y -
	//       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	//dz     =       coor[k][j][i].z-coor[k-1][j][i].z;
	flux_y -= ucont[k][j][i].z*
	  (ucat[k][j][i].y+ucat[k+1][j][i].y)*0.5;
	dy   = coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	// du_y/dz
	Ft_y  += (ucat[k+1][j][i].y-ucat[k][j][i].y)/Re/dz*dy*dx;
      }
    }
  }
  
  // Total Forces on each processor
  F_z = flux_z + Fp_z + Ft_z;
  F_y = flux_y + Fp_y + Ft_y;

  // Global Sum
  //PetscGlobalSum(&F_x, &F_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_y, &F_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_z, &F_zSum, PETSC_COMM_WORLD);

  // Scale Check later !!!!!
  //F_xSum=F_xSum/A_totSum*pi*2.;
  F_ySum=-2*F_ySum/0.1;///3.;
  F_zSum=-2*F_zSum/0.1;///3.;

  // store results in fsi
  fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum;

  // output values
  PetscPrintf(PETSC_COMM_SELF, "CV  F_z, %le %le %le %le  \n",F_z, flux_z, Fp_z,Ft_z);
  PetscPrintf(PETSC_COMM_SELF, "CV  F_y, %le %le %le %le  \n",F_y, flux_y, Fp_y,Ft_y);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_CVM2D");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le\n",ti, F_ySum, F_zSum);
    fclose(f);
  }

  // Restore Working arrays
  //DAVecRestoreArray(fda, user->lCent,&coor);
  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lUcat,  &ucat);
  DAVecRestoreArray(fda, user->lUcont, &ucont);
  DAVecRestoreArray(da, user->lP, &p);
  //VecDestroy(Coor);

  return(0);
}

PetscErrorCode Calc_forces_CVM2D_2(UserCtx *user, FSInfo *fsi, PetscInt ti) 
{ 

  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  // cell location 
  PetscInt      zin,zout,yin,yout,xin,xout;
  xin=0; xout=2;
  yin=60; yout=135;
  zin=22; zout=98;
/*   CV_Boundary(user,fsi); */
  yin=fsi->CV_ys;yout=fsi->CV_ye;
  zin=fsi->CV_zs;zout=fsi->CV_ze;

  if (zin<=0) zin=1;
  if (zout>=mz-2 && ze==mz) zout=mz-3;
  if (yin<=0) yin=1;
  if (yout>=my-2 && ye==my) yout=my-3;
  
  // center cell location
  PetscInt      zstart,zend,ystart,yend,xstart,xend;
  zstart=zin+1;
  zend  =zout+1;
  ystart=yin+1;
  yend  =yout+1;
  xstart=xin+1;
  xend  =xout+1;

  if (zstart<lzs) zstart=lzs;
  if (ystart<lys) ystart=lys;
  if (xstart<lxs) xstart=lxs;
  if (zend>lze) zend=lze;
  if (yend>lye) yend=lye;
  if (xend>lxe) xend=lxe;

  PetscPrintf(PETSC_COMM_SELF, "CV start end z %d %d y %d %d x %d %d %d %d\n",zstart,zend,ystart,yend,xstart,xend, lys, lye);
  PetscPrintf(PETSC_COMM_SELF, "CV in out    z %d %d y %d %d x %d %d %d %d\n",zin,zout, yin,yout, xin, xout, ys,ye);

  Vec           Coor;
  PetscInt	i, j, k;
  Cmpnts        ***coor, ***ucat, ***ucont, ***ucat_o;
  PetscReal     ***p, ***nvert;
  PetscReal     dx,dy,dz, dV;
  PetscReal     F_x,F_y,F_z; //Forces and Area
  PetscReal     F_xSum,F_ySum,F_zSum; //Surface Force
  PetscReal     Fvis_z, Fvis_y;
  PetscReal     Fvis_zSum, Fvis_ySum;
  PetscReal     Finv_z, Finv_y;
  PetscReal     Finv_zSum, Finv_ySum;
  PetscReal     flux_z,Fp_z, Ft_z, Fzz;
  PetscReal     flux_y,Fp_y, Ft_y, Fyy;
  PetscReal     massflux, massflux_z, MassFluxSum;
  PetscReal     Re=user->ren;
  PetscReal     dt=user->dt;
  PetscReal     nj_in,nj_out,nk_in,nk_out;
  PetscReal     Tzz,Tyy,Tyz;
  //intertial forces
  PetscReal     MCV_z, MCV_y, MCV_z_o, MCV_y_o;
  PetscReal     Fint_z, Fint_y;
  PetscReal     Fint_zSum, Fint_ySum;
  PetscReal     Finert_z, Finert_y;
  PetscReal     Finert_zSum, Finert_ySum;
  PetscReal     dwdt, dvdt;
  //PetscReal     Fg, FgSum;

  // Moving CV
  PetscReal     u_cv=fsi->S_new[1];
  PetscReal     v_cv=fsi->S_new[3];
  PetscReal     w_cv=fsi->S_new[5];
  Cmpnts        ***icsi,***jeta,***kzet;
  PetscReal     ***iaj,***jaj,***kaj;
  PetscReal     csi0,csi1,csi2,Jaj;
  PetscReal     Ucont_cv;

  PetscReal     dir=1., dir_y=1.;// depends on zet.z since U=zet.z*ucat.z

  // Get Working arrays
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);


  //DAVecGetArray(fda, user->lCent,&coor);
  DAVecGetArray(fda, user->lUcat,  &ucat);
  DAVecGetArray(fda, user->Ucat_o, &ucat_o);
  DAVecGetArray(fda, user->lUcont, &ucont);
  DAVecGetArray(da, user->lP, &p);
  DAVecGetArray(da, user->lNvert, &nvert);

  DAVecGetArray(fda, user->lICsi, &icsi);
  DAVecGetArray(fda, user->lJEta, &jeta);
  DAVecGetArray(fda, user->lKZet, &kzet);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

/*   DAVecGetArray(fda, user->GridSpace, &gs); */

  // calculate normals of the CV
  i=xstart;
  k=zstart;
  j=ystart;
  nj_in=coor[k][j][i].y-coor[k][j-1][i].y;
  nj_in=-nj_in/fabs(nj_in);
  j=yend;
  nj_out=coor[k][j+1][i].y-coor[k][j][i].y;
  nj_out=nj_out/fabs(nj_out);

  i=xstart;
  j=ystart;
  k=zstart;
  nk_in=coor[k][j][i].z-coor[k-1][j][i].z;
  nk_in=-nk_in/fabs(nk_in);
  k=zend;
  nk_out=coor[k+1][j][i].z-coor[k][j][i].z;
  nk_out=nk_out/fabs(nk_out);
 
  PetscPrintf(PETSC_COMM_SELF, "CV normal %le %le %le %le\n",nk_in,nk_out,nj_in,nj_out);
 
  Ucont_cv=0.;
  //                    z momentum 
  flux_z=0.;
  Fp_z=0.;
  Ft_z=0.;
  Fzz =0.;
  massflux_z=0.;
  // 1) inflow forces
  if (zin>zs) {
    k=zin; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	dy= fabs(coor[k][j][i].y- coor[k][j-1][i].y);
	dx= fabs(coor[k][j][i].x- coor[k][j][i-1].x);
	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	dz= fabs(dz);dy= fabs(dy);dx= fabs(dx);

/* 	csi0 = kzet[k][j][i].x; */
/* 	csi1 = kzet[k][j][i].y; */
/* 	csi2 = kzet[k][j][i].z; */
/* 	//Jaj = 1/ kaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_z -= dir*nk_in*(ucont[k][j][i].z-Ucont_cv)*
	  (ucat[k][j][i].z+ucat[k+1][j][i].z)*0.5;
	//dy      = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y -
	//	coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	Fp_z   -= nk_in*0.5*(p[k][j][i]+p[k+1][j][i])*dy*dx;
	// du_z/dz
	Tzz  = 2*dir*(ucat[k+1][j][i].z-ucat[k][j][i].z)/Re/dz;
	Fzz += nk_in*Tzz*dy*dx;
	massflux_z -= dir*nk_in*(ucont[k][j][i].z-Ucont_cv);
      }
    }
  }

  // 2) outflow forces
  if (zout<lze) {
    k=zout; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	dy= coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	dz= fabs(dz);
	dy= fabs(dy);
	dx= fabs(dx);

/* 	csi0 = kzet[k][j][i].x; */
/* 	csi1 = kzet[k][j][i].y; */
/* 	csi2 = kzet[k][j][i].z; */
/* 	//Jaj = 1/ kaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_z -= dir*nk_out*(ucont[k][j][i].z-Ucont_cv)*
	  (ucat[k][j][i].z+ucat[k+1][j][i].z)*0.5;
	//dy      = 0.25*(coor[k][j+1][i].y+coor[k+1][j+1][i].y -
	//	coor[k][j-1][i].y-coor[k+1][j-1][i].y);
	Fp_z   -= nk_out*0.5*(p[k][j][i]+p[k+1][j][i])*dy*dx;
	Tzz = 2*dir*(ucat[k+1][j][i].z-ucat[k][j][i].z)/Re/dz;
	Fzz  += nk_out*Tzz*dy*dx;
	massflux_z -= dir*nk_out*(ucont[k][j][i].z-Ucont_cv);
      }
    }
  }

  // 3) upperside forces
  if (yout<lye) {
    j=yout; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	//dz     = 0.25*(coor[k+1][j][i].z+coor[k+1][j+1][i].z -
	//       coor[k-1][j][i].z-coor[k-1][j+1][i].z);
	//dy     =       coor[k][j+1][i].y-coor[k][j][i].y;
/* 	csi0 = jeta[k][j][i].x; */
/* 	csi1 = jeta[k][j][i].y; */
/* 	csi2 = jeta[k][j][i].z; */
/* 	//Jaj = 1/ jaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_z -= dir_y*nj_out*(ucont[k][j][i].y-Ucont_cv) *
	   (ucat[k][j][i].z+ucat[k][j+1][i].z)*0.5;
	dz     = coor[k][j][i].z-coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	dz= fabs(dz);
	dy= fabs(dy);
	dx= fabs(dx);	
	// du_z/dy
	Tyz = (dir_y*(ucat[k][j+1][i].z-ucat[k][j][i].z)/dy +
	       dir*  (ucat[k+1][j][i].y-ucat[k-1][j][i].y +
		      ucat[k+1][j+1][i].y-ucat[k-1][j+1][i].y)*0.25/dz)/Re;
	Ft_z  += nj_out*Tyz*dz*dx;
	massflux_z -= dir_y*nj_out*(ucont[k][j][i].y-Ucont_cv);
      }
    }
  }

  // 4) lowerside forces
  if (yin>ys) {
    j=yin; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	//dz     = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
	//       coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	//dy     =       coor[k][j][i].y-coor[k][j-1][i].y;
	// du_z/dy
/* 	csi0 = jeta[k][j][i].x; */
/* 	csi1 = jeta[k][j][i].y; */
/* 	csi2 = jeta[k][j][i].z; */
/* 	//Jaj = 1/ jaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_z -= dir_y*nj_in*(ucont[k][j][i].y-Ucont_cv) *
	   (ucat[k][j][i].z+ucat[k][j+1][i].z)*0.5;
	dz     = coor[k][j][i].z-coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	dz= fabs(dz);
	dy= fabs(dy);
	dx= fabs(dx);	
	// du_z/dy
	Tyz =(dir_y*(ucat[k][j+1][i].z-ucat[k][j][i].z)/dy +
	      dir*  (ucat[k+1][j][i].y-ucat[k-1][j][i].y +
		     ucat[k+1][j+1][i].y-ucat[k-1][j+1][i].y)*0.25/dz)/Re;
	Ft_z  += nj_in*Tyz*dz*dx;
	/* 	Ft_z  -= (ucat[k][j][i].z-ucat[k][j+1][i].z)/Re/dy*dz*dx; */
	massflux_z -= dir_y*nj_in*(ucont[k][j][i].y-Ucont_cv);
      }
    }
  }

  //                    y momentum 
  flux_y=0.;
  Fp_y=0.;
  Ft_y=0.;
  Fyy=0.;
  massflux=0.;
  // 1) inflow forces
  if (zin>zs) {
    k=zin; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	//dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y -
	//       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	//dz     =       coor[k][j][i].z-coor[k-1][j][i].z;
/* 	csi0 = kzet[k][j][i].x; */
/* 	csi1 = kzet[k][j][i].y; */
/* 	csi2 = kzet[k][j][i].z; */
/* 	//Jaj = 1/ kaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_y -= dir*nk_in*(ucont[k][j][i].z-Ucont_cv)*
	  (ucat[k][j][i].y+ucat[k+1][j][i].y)*0.5;
	dy   = coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	dz= fabs(dz);
	dy= fabs(dy);
	dx= fabs(dx);	
	// du_y/dz+du_z/dy
/* 	dvdz = dir*  (ucat[k+1][j][i].y-ucat[k][j][i].y)/dz; */
/* 	dwdy = dir_y*(ucat[k][j+1][i].z-ucat[k][j-1][i].z + */
/* 		      ucat[k+1][j+1][i].z-ucat[k+1][j-1][i].z)*0.25/dy; */
	Tyz= (dir*  (ucat[k+1][j][i].y-ucat[k][j][i].y)/dz+
	      dir_y*(ucat[k][j+1][i].z-ucat[k][j-1][i].z +
		     ucat[k+1][j+1][i].z-ucat[k+1][j-1][i].z)*0.25/dy)/Re;
	Ft_y  += nk_in*Tyz*dy*dx;
	massflux += dir*nk_in*(ucont[k][j][i].z-Ucont_cv);
      }
    }
  }

/*   PetscPrintf(PETSC_COMM_SELF, "CV fluxes, %le %le %le mass %le\n", flux_y, Fp_y,Ft_y, massflux); */

  // 2) outflow forces
  if (zout<lze) {
    k=zout; 
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	//dy     = 0.25*(coor[k][j+1][i].y+coor[k+1][j+1][i].y -
	//       coor[k][j-1][i].y-coor[k+1][j-1][i].y);
	//dz     =       coor[k+1][j][i].z-coor[k][j][i].z;
/* 	csi0 = kzet[k][j][i].x; */
/* 	csi1 = kzet[k][j][i].y; */
/* 	csi2 = kzet[k][j][i].z; */
/* 	//Jaj = 1/ kaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_y -= dir*nk_out*(ucont[k][j][i].z-Ucont_cv)*
	  (ucat[k][j][i].y+ucat[k+1][j][i].y)*0.5;
	dy   = coor[k][j][i].y- coor[k][j-1][i].y;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

	dz   = 0.25*(coor[k+1][j][i].z+coor[k+1][j-1][i].z -
		     coor[k-1][j][i].z-coor[k-1][j-1][i].z);
	dz= fabs(dz);
	dy= fabs(dy);
	dx= fabs(dx);	
	// du_y/dz
	Tyz= (dir*  (ucat[k+1][j][i].y-ucat[k][j][i].y)/dz+
	      dir_y*(ucat[k][j+1][i].z-ucat[k][j-1][i].z +
		     ucat[k+1][j+1][i].z-ucat[k+1][j-1][i].z)*0.25/dy)/Re;
	Ft_y  += nk_out*Tyz*dy*dx;
/* 	Ft_y  -= (ucat[k+1][j][i].y-ucat[k][j][i].y)/Re/dz*dy*dx; */
	massflux += dir*nk_out*(ucont[k][j][i].z-Ucont_cv);
      }
    }
  }

/*   PetscPrintf(PETSC_COMM_SELF, "CV fluxes, %le %le %le mass %le\n", flux_y, Fp_y,Ft_y, massflux); */

  // 3) upperside forces
  if (yout<lye) {
    j=yout; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	dz      = coor[k][j][i].z - coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

/* 	csi0 = jeta[k][j][i].x; */
/* 	csi1 = jeta[k][j][i].y; */
/* 	csi2 = jeta[k][j][i].z; */
/* 	//Jaj = 1/ jaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_y -= dir_y*nj_out*(ucont[k][j][i].y-Ucont_cv)*
	  (ucat[k][j][i].y+ucat[k][j+1][i].y)*0.5;
	//dz      = 0.25*(coor[k+1][j][i].z+coor[k-1][j+1][i].z -
	//	coor[k+1][j][i].z-coor[k-1][j+1][i].z);
	dz= fabs(dz);	
	dx= fabs(dx);
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	dy= fabs(dy);
	Fp_y   -= nj_out*0.5*(p[k][j][i]+p[k][j+1][i])*dz*dx;
	//du_y/dy
	Tyy = 2*dir_y*(ucat[k][j+1][i].y-ucat[k][j][i].y)/Re/dy;
	Fyy    += nj_out*Tyy*dz*dx;
	massflux += dir_y*nj_out*(ucont[k][j][i].y-Ucont_cv);
      }
    }
  }

/*   PetscPrintf(PETSC_COMM_SELF, "CV fluxes, %le %le %le mass %le\n", flux_y, Fp_y,Ft_y, massflux); */

  // 4) lowerside forces
  if (yin>ys) {
    j=yin; 
    for (k=zstart; k<zend; k++) {
      for (i=xstart; i<xend; i++) {
	dz     = coor[k][j][i].z - coor[k-1][j][i].z;
	dx= coor[k][j][i].x- coor[k][j][i-1].x;

/* 	csi0 = jeta[k][j][i].x; */
/* 	csi1 = jeta[k][j][i].y; */
/* 	csi2 = jeta[k][j][i].z; */
/* 	//Jaj = 1/ jaj[k][j][i]; */
/* 	Ucont_cv = (u_cv*csi0+v_cv*csi1+w_cv*csi2); */

	flux_y -= dir_y*nj_in*(ucont[k][j][i].y-Ucont_cv)*
	  (ucat[k][j][i].y+ucat[k][j+1][i].y)*0.5;
	//dz      = 0.25*(coor[k+1][j][i].z+coor[k-1][j-1][i].z -
	//	coor[k+1][j][i].z-coor[k-1][j-1][i].z);
	dy     = 0.25*(coor[k][j+1][i].y+coor[k-1][j+1][i].y-
		       coor[k][j-1][i].y-coor[k-1][j-1][i].y);
	dy= fabs(dy);
	dz= fabs(dz);	
	dx= fabs(dx);	
	Fp_y   -= nj_in*0.5*(p[k][j][i]+p[k][j+1][i])*dz*dx;
	//du_y/dy
	Tyy = 2*dir_y*(ucat[k][j+1][i].y-ucat[k][j][i].y)/Re/dy;
	Fyy    += nj_in*Tyy*dz*dx;
/* 	Fyy    += (ucat[k][j+1][i].y-ucat[k][j][i].y)/dy*dz*dx; */
	massflux += dir_y*nj_in*(ucont[k][j][i].y-Ucont_cv);
      }
    }
  }

/*   PetscPrintf(PETSC_COMM_SELF, "CV fluxes, %le %le %le mass %le\n", flux_y, Fp_y,Ft_y, massflux); */

  // Inertial Forces
  MCV_z =0.; MCV_y=0.;
  MCV_z_o = 0.; MCV_y_o=0.;
  Finert_z = 0.; Finert_y =0.;
  Fint_z = 0.; Fint_y=0.;
  for (k=zstart; k<zend; k++) {
    for (j=ystart; j<yend; j++) {
      for (i=xstart; i<xend; i++) {
	if (nvert[k][j][i]<.2) {
	  dz= coor[k][j][i].z - coor[k-1][j][i].z;
	  dy= coor[k][j][i].y - coor[k][j-1][i].y;
	  dx= coor[k][j][i].x - coor[k][j][i-1].x;
  
	  dx= fabs(dx);dy= fabs(dy);dz= fabs(dz);	
	  
	  dV = dx*dy*dz;
	  // dV = 1/aj 
	  
	  MCV_z += ucat[k][j][i].z * dV;
	  MCV_y += ucat[k][j][i].y * dV;

	  MCV_z_o += ucat_o[k][j][i].z * dV;
	  MCV_y_o += ucat_o[k][j][i].y * dV;

	  dwdt = (ucat[k][j][i].z - ucat_o[k][j][i].z)/dt;
	  dvdt = (ucat[k][j][i].y - ucat_o[k][j][i].y)/dt;

	  Finert_z += dwdt*dV;
	  Finert_y += dvdt*dV;
	  
	  if (i==1 && j== ystart+2 && ( k==zstart || k==zstart+2))
	  PetscPrintf(PETSC_COMM_WORLD, "MCV , z %le %le y %le %le \n", ucat[k][j][i].z, ucat_o[k][j][i].z, ucat[k][j][i].y, ucat_o[k][j][i].y);	
	}
      }
    }
  }
  Fint_z = (MCV_z - MCV_z_o)/ dt;
  Fint_y = (MCV_y - MCV_y_o)/ dt;
  
  // Total Forces on each processor
  F_z = flux_z + Fp_z + Ft_z + Fzz -Fint_z;
  F_y = flux_y + Fp_y + Ft_y + Fyy -Fint_y;

  Finv_z = flux_z + Fp_z ;
  Finv_y = flux_y + Fp_y ;
  
  Fvis_z = Ft_z ;
  Fvis_y = Ft_y ;

  // Global Sum
  //PetscGlobalSum(&F_x, &F_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_y, &F_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_z, &F_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Finert_z, &Finert_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Finert_y, &Finert_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Fint_z, &Fint_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Fint_y, &Fint_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&massflux, &MassFluxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Finv_z, &Finv_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Finv_y, &Finv_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Fvis_z, &Fvis_zSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Fvis_y, &Fvis_ySum, PETSC_COMM_WORLD);

  // Scale Check later !!!!!
  F_xSum=0.;//F_xSum/A_totSum*pi*2.;
  F_ySum=2*F_ySum/0.1;///3.;
  F_zSum=2*F_zSum/0.1;///3.;

  Finv_ySum=2*Finv_ySum/0.1;///3.;
  Finv_zSum=2*Finv_zSum/0.1;///3.;

  Fvis_ySum=2*Fvis_ySum/0.1;///3.;
  Fvis_zSum=2*Fvis_zSum/0.1;///3.;

  Fint_ySum=2*Fint_ySum/0.1;///3.;
  Fint_zSum=2*Fint_zSum/0.1;///3.;

  // store results in fsi
/*   fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum; */

  // output values
/*   PetscPrintf(PETSC_COMM_SELF, "CV  F_z, %le %le %le %le %le %le %le mass %le\n",F_z, flux_z, Fp_z,Ft_z,Fzz,Fint_z,Finert_z, massflux_z); */
  PetscPrintf(PETSC_COMM_SELF, "CV  F_y, %le %le %le %le %le %le %le mass %le\n",F_y, flux_y, Fp_y,Ft_y,Fyy,Fint_y,Finert_y, massflux);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_CVM2D");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_ySum, F_zSum, Finv_ySum, Finv_zSum,Fvis_ySum, Fvis_zSum, Fint_ySum, Fint_zSum,MassFluxSum);
    fclose(f);
  }


  // Restore Working arrays
  //DAVecRestoreArray(fda, user->lCent,&coor);
  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lUcat,  &ucat);
  DAVecRestoreArray(fda, user->Ucat_o,  &ucat_o);
  DAVecRestoreArray(fda, user->lUcont, &ucont);
  DAVecRestoreArray(da, user->lP, &p);
  DAVecRestoreArray(da, user->lNvert, &nvert);

  DAVecRestoreArray(fda, user->lICsi, &icsi);
  DAVecRestoreArray(fda, user->lJEta, &jeta);
  DAVecRestoreArray(fda, user->lKZet, &kzet);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);

/*   DAVecRestoreArray(fda, user->GridSpace, &gs); */
  //VecDestroy(Coor);

  return(0);
}

PetscErrorCode CV_Boundary(UserCtx *user, FSInfo *fsi, PetscInt ibi)
{
  PetscInt	i, j, k;  
  IBMInfo       *ibminfo;
  IBMListNode   *current;

  fsi->CV_ys=9999;
  fsi->CV_zs=9999;
  fsi->CV_ye=0;
  fsi->CV_ze=0;

  current = user->ibmlist[ibi].head;
  if (!(current->next)) {
    fsi->CV_ys=0;
    fsi->CV_zs=0;
  } else {
  while (current) {
    ibminfo = &current->ibm_intp;
    current = current->next;	
 
    i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
    
    fsi->CV_ys=PetscMin(fsi->CV_ys,j);
    fsi->CV_ye=PetscMax(fsi->CV_ye,j);

    fsi->CV_zs=PetscMin(fsi->CV_zs,k);
    fsi->CV_ze=PetscMax(fsi->CV_ze,k);
  }

  fsi->CV_ys=fsi->CV_ys-radi;
  fsi->CV_ye=fsi->CV_ye+radi;
  fsi->CV_zs=fsi->CV_zs-radi;
  fsi->CV_ze=fsi->CV_ze+radi;
  }
  return(0);
}

    
PetscErrorCode Calc_forces_SI(FSInfo *fsi,UserCtx *user,
			      IBMNodes *ibm,PetscInt ti, 
			      PetscInt ibi, PetscInt bi)
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal     rei= 1./user->ren;
  Vec           Coor;
  PetscInt	i, j, k, elmt;
  PetscInt	nv1,nv2,nv3;
  Cmpnts        uelmt;
  PetscReal     sb;
  PetscReal     Tow_ws,Tow_wt,Tow_wn;
  PetscReal     Tow_x, Tow_y, Tow_z;
  Cmpnts        ***coor, ***ucat, ***cent;
  PetscReal     ***p, ***nvert;
  PetscReal     dx,dy,dz, dx1,dy1,dz1;
  PetscReal     dwdz,dwdy,dwdx;
  PetscReal     dvdz,dvdy,dvdx;
  PetscReal     dudz,dudy,dudx;
  PetscReal     Txx,Tyy,Tzz;
  PetscReal     Tzy,Tzx,Tyx;
  PetscReal     nfx,nfy,nfz;
  PetscReal     nsx,nsy,nsz;
  PetscReal     ntx,nty,ntz;
  PetscReal     F_px,F_py,F_pz,Ap_x,Ap_y,Ap_z,Ap_t; //Forces and Area
  PetscReal     F_nx,F_ny,F_nz,An_x,An_y,An_z,An_t; 
  //PetscReal     Cp_x,Cp_y,Cp_z; //Pressure Forces
  //PetscReal     Cs_x,Cs_y,Cs_z; //Surface Forces
  PetscReal     Cp_nx,Cp_ny,Cp_nz; //Pressure Forces - side
  PetscReal     Cs_nx,Cs_ny,Cs_nz; //Surface Forces - side
  PetscReal     Cp_px,Cp_py,Cp_pz; //Pressure Forces + side
  PetscReal     Cs_px,Cs_py,Cs_pz; //Surface Forces + side
  PetscReal     F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force
  PetscReal     F_pxSum,F_pySum,F_pzSum; //Surface Force
  PetscReal     F_nxSum,F_nySum,F_nzSum; //Surface Force
  PetscReal     Ap_xSum,Ap_ySum,Ap_zSum,Ap_tSum; // + side
  PetscReal     An_xSum,An_ySum,An_zSum,An_tSum; // - side
  PetscReal     A_xSum,A_ySum,A_zSum,A_tSum;
  PetscReal     Cp_xSum,Cp_ySum,Cp_zSum; //Pressure Force
  PetscReal     Cp_pxSum,Cp_pySum,Cp_pzSum; //Pressure Force
  PetscReal     Cp_nxSum,Cp_nySum,Cp_nzSum; //Pressure Force

  // Moments
  PetscReal      MF_px,MF_py,MF_pz; //Forces for Moment Calc
  PetscReal      MF_nx,MF_ny,MF_nz; //Forces for Moment Calc
  PetscReal      M_nx,M_ny,M_nz;   //Moments
  PetscReal      M_px,M_py,M_pz;   //Moments
  PetscReal      r_x,r_y,r_z;   //Anchor dist
  PetscReal      x, y, z;       //cell coord
  PetscReal      X_c,Y_c,Z_c;   //center of rotation coord
  PetscReal      M_xSum,M_ySum,M_zSum; //Surface Mom on all processors
  PetscReal      M_pxSum,M_pySum,M_pzSum; //Surface Mom on all processors
  PetscReal      M_nxSum,M_nySum,M_nzSum; //Surface Mom on all processors
  PetscReal      Iap_x,Iap_y,Iap_z;
  PetscReal      Ian_x,Ian_y,Ian_z;
  PetscReal      Iap_xSum,Iap_ySum,Iap_zSum; // + side
  PetscReal      Ian_xSum,Ian_ySum,Ian_zSum; // - side
  PetscReal      Ia_xSum,Ia_ySum,Ia_zSum;
  PetscReal      Pw_nx,Pw_ny,Pw_nz;   //Power
  PetscReal      Pw_px,Pw_py,Pw_pz;   //Power
  PetscReal      Pw_nxSum,Pw_nySum,Pw_nzSum;   //Power
  PetscReal      Pw_pxSum,Pw_pySum,Pw_pzSum;   //Power
  PetscReal      Pw_xSum,Pw_ySum,Pw_zSum;

  PetscReal      A_x,A_y,A_z,Atot;
  PetscReal      u_x,u_y,u_z;
  Cmpnts         ***csi,***eta,***zet;
  PetscReal      csi1,csi2,csi3;
  PetscReal      eta1,eta2,eta3;
  PetscReal      zet1,zet2,zet3;
  PetscReal      ***iaj,***jaj,***kaj;

  PetscReal      dr;
  PetscReal      MFdpdn_px,MFdpdn_py,MFdpdn_pz;
  PetscReal      Mdpdn_px,Mdpdn_py,Mdpdn_pz;
  PetscReal      Mdpdn_pxSum,Mdpdn_pySum,Mdpdn_pzSum;
  PetscReal      MFdpdn_nx,MFdpdn_ny,MFdpdn_nz;
  PetscReal      Mdpdn_nx,Mdpdn_ny,Mdpdn_nz;
  PetscReal      Mdpdn_nxSum,Mdpdn_nySum,Mdpdn_nzSum;
  PetscReal      Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum; 

  IBMInfo        *ibminfo;
  IBMListNode    *current;

  // Fish
  PetscReal      pi = 3.141592653589793;
  PetscReal      v_side, dzz, az0, Pw_side, Thrust, Drag;
  PetscReal      Pw_sideSum, ThrustSum, DragSum, efficiency;
  PetscReal      time=user->dt*ti;  
  PetscReal  	 ampl=0.1;
  PetscReal 	 omega=2*pi*St_exp/(2.*ampl);//2.*pi; //non-dimensionalized w
  PetscReal	 kwave= 0;//2*pi/wavelength;//8.307767;


/* ==================================================================================             */
/*   Init var */
  F_px=0.;F_py=0.;F_pz=0.;
  F_nx=0.;F_ny=0.;F_nz=0.;
  Ap_x=0.;Ap_y=0.;Ap_z=0.;Ap_t=0.;
  An_x=0.;An_y=0.;An_z=0.;An_t=0.;
  Cp_px=0.;Cp_py=0.;Cp_pz=0.;
  Cs_px=0.;Cs_py=0.;Cs_pz=0.;
  Cp_nx=0.;Cp_ny=0.;Cp_nz=0.;
  Cs_nx=0.;Cs_ny=0.;Cs_nz=0.;

  M_px=0.;M_py=0.;M_pz=0.;
  M_nx=0.;M_ny=0.;M_nz=0.;
  Iap_x=0.;Iap_y=0.;Iap_z=0.;
  Ian_x=0.;Ian_y=0.;Ian_z=0.;
  Pw_px=0.;Pw_py=0.;Pw_pz=0.;
  Pw_nx=0.;Pw_ny=0.;Pw_nz=0.;

  // Fish
  Pw_side=0.; Thrust=0.; Drag=0.;
  

  Mdpdn_px=0.;Mdpdn_py=0.;Mdpdn_pz=0.;
  Mdpdn_nx=0.;Mdpdn_ny=0.;Mdpdn_nz=0.;
  MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
  MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;

  /*   Check Later Y_c !!!!!!!!!!!!!!!!!!!!!!!! */
  X_c=fsi->x_c; Y_c=fsi->y_c; Z_c=fsi->z_c;
  X_c=x_r, Y_c=y_r, Z_c=z_r; //seokkoo

  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force  %le X_c %le %le %le fish %le %le %le\n",rei, X_c,Y_c,Z_c,omega,kwave,time);

/* ==================================================================================             */
/*   Get Working arrays */
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
  DAVecGetArray(fda, user->lCent, &cent);
  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(da, user->lP, &p);
  DAVecGetArray(da, user->lNvert, &nvert);

  DAVecGetArray(fda, user->lCsi, &csi);
  DAVecGetArray(fda, user->lEta, &eta);
  DAVecGetArray(fda, user->lZet, &zet);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

/* ==================================================================================             */
/* Loop around all ibm nodes */
  current = user->ibmlist[ibi].head;
  while (current) {
    ibminfo = &current->ibm_intp;
    current = current->next;	
    i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk;
    elmt = ibminfo->cell; // closest ibm element
    sb = ibminfo->d_s;
    //sb = ibminfo->d_i;// |dn| from node to bndry

    // normal 
    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    // 1st bi-normal of nf
    nsx=ibm->ns_x[elmt];
    nsy=ibm->ns_y[elmt];
    nsz=ibm->ns_z[elmt];

    // 2nd bi-normal of nf
    ntx=ibm->nt_x[elmt];
    nty=ibm->nt_y[elmt];
    ntz=ibm->nt_z[elmt];

    // nodes of closest ibm elmnt
    nv1=ibm->nv1[elmt];
    nv2=ibm->nv2[elmt];
    nv3=ibm->nv3[elmt];

    //velocity of the closest elmnt
    uelmt.x = (ibm->u[nv1].x+ibm->u[nv2].x+ibm->u[nv3].x)/3.;
    uelmt.y = (ibm->u[nv1].y+ibm->u[nv2].y+ibm->u[nv3].y)/3.;
    uelmt.z = (ibm->u[nv1].z+ibm->u[nv2].z+ibm->u[nv3].z)/3.;
    
    if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze) {	        

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;
      MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
      MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;


/* ==================================================================================             */
/*       Shear Stresses (2nd & 1st order) and Shear Force */

      if (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k+1<mz-1 && k-1>0) {
	zet1 = 0.25*(zet[k][j][i].x+zet[k-1][j][i].x)*
	              (kaj[k][j][i]+kaj[k-1][j][i]);
	zet2 = 0.25*(zet[k][j][i].y+zet[k-1][j][i].y)*
	              (kaj[k][j][i]+kaj[k-1][j][i]);
	zet3 = 0.25*(zet[k][j][i].z+zet[k-1][j][i].z)*
	              (kaj[k][j][i]+kaj[k-1][j][i]);
	dwdz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet1;

      } else if (nvert[k+1][j][i]<2.5 && k+1<mz-1) {
	zet1 = (zet[k][j][i].x)*kaj[k][j][i];
	zet2 = (zet[k][j][i].y)*kaj[k][j][i];
	zet3 = (zet[k][j][i].z)*kaj[k][j][i];

	dwdz = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet1;

      } else if (nvert[k-1][j][i]<2.5 && k-1>0){
	zet1 = (zet[k-1][j][i].x)*kaj[k-1][j][i];
	zet2 = (zet[k-1][j][i].y)*kaj[k-1][j][i];
	zet3 = (zet[k-1][j][i].z)*kaj[k-1][j][i];

	dwdz = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet3;
	dvdz = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet3;
	dudz = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet3;

	dwdy = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet2;
	dvdy = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet2;
	dudy = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet2;

	dwdx = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet1;
	dvdx = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet1;
	dudx = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet1;

      } else {
	dwdz = 0.;
	dvdz = 0.;
	dudz = 0.;

	dwdy = 0.;
	dvdy = 0.;
	dudy = 0.;

	dwdx = 0.;
	dvdx = 0.;
	dudx = 0.;
      }


      if (nvert[k][j+1][i]<2.5 && j+1<my-1 && nvert[k][j-1][i]<2.5 && j-1>0) {
	eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*
	              (jaj[k][j][i]+jaj[k][j-1][i]);
	eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*
	              (jaj[k][j][i]+jaj[k][j-1][i]);
	eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*
	              (jaj[k][j][i]+jaj[k][j-1][i]);
	dwdz += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta1;

      } else if (nvert[k][j+1][i]<2.5 && j+1<my-1) {
	eta1 = eta[k][j][i].x*jaj[k][j][i];
	eta2 = eta[k][j][i].y*jaj[k][j][i];
	eta3 = eta[k][j][i].z*jaj[k][j][i];

	dwdz += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta1;

      } else if (nvert[k][j-1][i]<2.5 && j-1>0){
	eta1 = eta[k][j-1][i].x*jaj[k][j-1][i];
	eta2 = eta[k][j-1][i].y*jaj[k][j-1][i];
	eta3 = eta[k][j-1][i].z*jaj[k][j-1][i];

	dwdz += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta3;
	dvdz += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta3;
	dudz += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta3;

	dwdy += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta2;
	dvdy += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta2;
	dudy += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta2;

	dwdx += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta1;
	dvdx += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta1;
	dudx += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta1;
      } 

      if (nvert[k][j][i+1]<2.5 && i+1<mx-1 && nvert[k][j][i-1]<2.5 && i-1>0) {
	csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)*
	              (iaj[k][j][i]+iaj[k][j][i-1]);
	csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)*
	              (iaj[k][j][i]+iaj[k][j][i-1]);
	csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)*
	              (iaj[k][j][i]+iaj[k][j][i-1]);

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi1;

      } else if (nvert[k][j][i+1]<2.5 && i+1<mx-1) {
	csi1 = csi[k][j][i].x*iaj[k][j][i];
	csi2 = csi[k][j][i].y*iaj[k][j][i];
	csi3 = csi[k][j][i].z*iaj[k][j][i];

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi1;

      } else if (nvert[k][j][i-1]<2.5 && i-1>0){
	csi1 = csi[k][j][i-1].x*iaj[k][j][i-1];
	csi2 = csi[k][j][i-1].y*iaj[k][j][i-1];
	csi3 = csi[k][j][i-1].z*iaj[k][j][i-1];

	dwdz += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi3;
	dvdz += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi3;
	dudz += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi3;

	dwdy += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi2;
	dvdy += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi2;
	dudy += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi2;

	dwdx += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi1;
	dvdx += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi1;
	dudx += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi1;
      }
      
      Tzz = rei * (dwdz + dwdz);
      Tyy = rei * (dvdy + dvdy);
      Txx = rei * (dudx + dudx);
      Tzy = rei * (dwdy + dvdz);
      Tzx = rei * (dwdx + dudz);
      Tyx = rei * (dvdx + dudy);


      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      r_x = 0.;
      r_y = 0.;
      r_z = 0.;
/* ==================================================================================             */
/*       K +  */

      if (nvert[k+1][j][i]<0.9 ){
	x = 0.5*(cent[k][j][i].x+cent[k+1][j][i].x);
	y = 0.5*(cent[k][j][i].y+cent[k+1][j][i].y);
	z = 0.5*(cent[k][j][i].z+cent[k+1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	A_z=zet[k  ][j][i].z;
	A_y=zet[k  ][j][i].y;
	A_x=zet[k  ][j][i].x;

	Cp_pz += -p[k+1][j][i]*A_z;
	MF_pz += -p[k+1][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

	if (MF_pz<0) 
	  Thrust += MF_pz;
	else 
	  Drag += MF_pz;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

	Cp_py += -p[k+1][j][i]*A_y;
	MF_py += -p[k+1][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */
	Cp_px += -p[k+1][j][i]*A_x;
	MF_px += -p[k+1][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	if (Tzz*A_z<0) 
	  Thrust += Tzz*A_z;
	else 
	  Drag += Tzz*A_z;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	if (Tzy*A_y<0) 
	  Thrust += Tzy*A_y;
	else 
	  Drag += Tzy*A_y;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	if (Tzx*A_x<0) 
	  Thrust += Tzx*A_x;
	else 
	  Drag += Tzx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	Ap_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k+1][j][i].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k+1][j][i].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k+1][j][i].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_px += MF_px * u_x ;//tot;
	Pw_py += MF_py * u_y ;//* Atot;
	Pw_pz += MF_pz * u_z ;//* Atot;

      }   

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       K -  */
   
      if (nvert[k-1][j][i]<0.9) {
	A_z=zet[k-1][j][i].z;
	A_y=zet[k-1][j][i].y;
	A_x=zet[k-1][j][i].x;

	x = 0.5*(cent[k][j][i].x+cent[k-1][j][i].x);
	y = 0.5*(cent[k][j][i].y+cent[k-1][j][i].y);
	z = 0.5*(cent[k][j][i].z+cent[k-1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	Cp_nz +=  p[k-1][j][i]*A_z;
	MF_nz +=  p[k-1][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

	if (MF_nz<0) 
	  Thrust += MF_nz;
	else 
	  Drag += MF_nz;

	Cp_ny +=  p[k-1][j][i]*A_y;
	MF_ny +=  p[k-1][j][i]*A_y;
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

	Cp_nx +=  p[k-1][j][i]*A_x;
	MF_nx +=  p[k-1][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	if (-Tzz*A_z<0) 
	  Thrust -= Tzz*A_z;
	else 
	  Drag -= Tzz*A_z;

	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;

	if (-Tzy*A_y<0) 
	  Thrust -= Tzy*A_y;
	else 
	  Drag -= Tzy*A_y;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	if (-Tzx*A_x<0) 
	  Thrust -= Tzx*A_x;
	else 
	  Drag -= Tzx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	An_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k-1][j][i].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k-1][j][i].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k-1][j][i].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_nx += MF_nx * u_x ;//* Atot;
	Pw_ny += MF_ny * u_y ;//* Atot;
	Pw_nz += MF_nz * u_z ;//* Atot;

      }
     
      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       j +  */

      if (nvert[k][j+1][i]<0.9 ){
	A_z=eta[k][j  ][i].z;
	A_y=eta[k][j  ][i].y;
	A_x=eta[k][j  ][i].x;

	x = 0.5*(cent[k][j][i].x+cent[k][j+1][i].x);
	y = 0.5*(cent[k][j][i].y+cent[k][j+1][i].y);
	z = 0.5*(cent[k][j][i].z+cent[k][j+1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	Cp_pz += -p[k][j+1][i]*A_z;
	MF_pz += -p[k][j+1][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

	Cp_py += -p[k][j+1][i]*A_y;
	MF_py += -p[k][j+1][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

	Cp_px += -p[k][j+1][i]*A_x;
	MF_px += -p[k][j+1][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	if (Tzz*A_z<0) 
	  Thrust += Tzz*A_z;
	else 
	  Drag += Tzz*A_z;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	if (Tzy*A_y<0) 
	  Thrust += Tzy*A_y;
	else 
	  Drag += Tzy*A_y;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	if (Tzx*A_x<0) 
	  Thrust += Tzx*A_x;
	else 
	  Drag += Tzx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	Ap_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k][j+1][i].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k][j+1][i].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k][j+1][i].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_px += MF_px * u_x ;//* Atot;
	Pw_py += MF_py * u_y ;//* Atot;
	Pw_pz += MF_pz * u_z ;//* Atot;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       j -  */

      if (nvert[k][j-1][i]<0.9) {  
	A_z=eta[k][j-1][i].z;
	A_y=eta[k][j-1][i].y;
	A_x=eta[k][j-1][i].x;

	x = 0.5*(cent[k][j][i].x+cent[k][j-1][i].x);
	y = 0.5*(cent[k][j][i].y+cent[k][j-1][i].y);
	z = 0.5*(cent[k][j][i].z+cent[k][j-1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	Cp_nz +=  p[k][j-1][i]*A_z;
	MF_nz +=  p[k][j-1][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

	Cp_ny +=  p[k][j-1][i]*A_y;
	MF_ny +=  p[k][j-1][i]*A_y;
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

	Cp_nx +=  p[k][j-1][i]*A_x;
	MF_nx +=  p[k][j-1][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	if (-Tzz*A_z<0) 
	  Thrust -= Tzz*A_z;
	else 
	  Drag -= Tzz*A_z;

	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;

	if (-Tzy*A_y<0) 
	  Thrust -= Tzy*A_y;
	else 
	  Drag -= Tzy*A_y;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	if (-Tzx*A_x<0) 
	  Thrust -= Tzx*A_x;
	else 
	  Drag -= Tzx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	An_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k][j-1][i].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k][j-1][i].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k][j-1][i].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_nx += MF_nx * u_x ;//* Atot;
	Pw_ny += MF_ny * u_y ;//* Atot;
	Pw_nz += MF_nz * u_z ;//* Atot;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       i +  */

      if (nvert[k][j][i+1]<0.9){
	A_z=csi[k][j][i].z;
	A_y=csi[k][j][i].y;
	A_x=csi[k][j][i].x;

	x = 0.5*(cent[k][j][i].x+cent[k][j][i+1].x);
	y = 0.5*(cent[k][j][i].y+cent[k][j][i+1].y);
	z = 0.5*(cent[k][j][i].z+cent[k][j][i+1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	Cp_pz += -p[k][j][i+1]*A_z;
	MF_pz += -p[k][j][i+1]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

	Cp_py += -p[k][j][i+1]*A_y;
	MF_py += -p[k][j][i+1]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

	Cp_px += -p[k][j][i+1]*A_x;
	MF_px += -p[k][j][i+1]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	if (Tzz*A_z<0) 
	  Thrust += Tzz*A_z;
	else 
	  Drag += Tzz*A_z;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	if (Tzy*A_y<0) 
	  Thrust += Tzy*A_y;
	else 
	  Drag += Tzy*A_y;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	if (Tzx*A_x<0) 
	  Thrust += Tzx*A_x;
	else 
	  Drag += Tzx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	Ap_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k][j][i+1].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k][j][i+1].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k][j][i+1].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_px += MF_px * u_x ;//* Atot;
	Pw_py += MF_py * u_y ;//* Atot;
	Pw_pz += MF_pz * u_z ;//* Atot;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       i -  */

      if  (nvert[k][j][i-1]<0.9) { 
	A_z=csi[k][j][i-1].z;
	A_y=csi[k][j][i-1].y;
	A_x=csi[k][j][i-1].x;

	x = 0.5*(cent[k][j][i].x+cent[k][j][i-1].x);
	y = 0.5*(cent[k][j][i].y+cent[k][j][i-1].y);
	z = 0.5*(cent[k][j][i].z+cent[k][j][i-1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

	Cp_nz +=  p[k][j][i-1]*A_z;
	MF_nz +=  p[k][j][i-1]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

	Cp_ny +=  p[k][j][i-1]*A_y;
	MF_ny +=  p[k][j][i-1]*A_y;
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

	Cp_nx +=  p[k][j][i-1]*A_x;
	MF_nx +=  p[k][j][i-1]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	if (-Tzz*A_z<0) 
	  Thrust -= Tzz*A_z;
	else 
	  Drag -= Tzz*A_z;

	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;

	if (-Tzy*A_y<0) 
	  Thrust -= Tzy*A_y;
	else 
	  Drag -= Tzy*A_y;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	if (-Tzx*A_x<0) 
	  Thrust -= Tzx*A_x;
	else 
	  Drag -= Tzx*A_x;
	
	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

	Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
	An_t += Atot;

/* 	u_x = 0.5*(ucat[k][j][i].x + ucat[k][j][i-1].x); */
/* 	u_y = 0.5*(ucat[k][j][i].y + ucat[k][j][i-1].y); */
/* 	u_z = 0.5*(ucat[k][j][i].z + ucat[k][j][i-1].z); */

	u_x = uelmt.x;
	u_y = uelmt.y;
	u_z = uelmt.z;

	Pw_nx += MF_nx * u_x ;//* Atot;
	Pw_ny += MF_ny * u_y ;//* Atot;
	Pw_nz += MF_nz * u_z ;//* Atot;

      }

      if (fish) {
	z = cent[k][j][i].z;
	dzz = z - 1.5;       
	az0=(ampl-0.02)*dzz+0.02;	
	v_side=az0*omega*cos(kwave*dzz-omega*time);
	
	Pw_side +=(MF_px+MF_nx)*v_side;

	//	PetscPrintf(PETSC_COMM_SELF, "MF %d %d %le %le %le %le %le %le\n",j,k,az0,dzz,M_px,M_nx, v_side,Pw_side); 
      }

      Mdpdn_px -=   r_y*MFdpdn_pz - r_z*MFdpdn_py;//
      Mdpdn_py -= -(r_x*MFdpdn_pz - r_z*MFdpdn_px);
      Mdpdn_pz -=   r_x*MFdpdn_py - r_y*MFdpdn_px;

      Mdpdn_nx -=   r_y*MFdpdn_nz - r_z*MFdpdn_ny;//
      Mdpdn_ny -= -(r_x*MFdpdn_nz - r_z*MFdpdn_nx);
      Mdpdn_nz -=   r_x*MFdpdn_ny - r_y*MFdpdn_nx;
      
/*       if (ibi==1 & i==2) */
/* 	PetscPrintf(PETSC_COMM_SELF, "MF %d %d %le %le %le %le %le %le %le %le\n",j,k,MF_py,MF_ny,MF_pz,MF_nz,M_px,M_nx,r_y,r_z); */

    } // if ibm node in CPU
/* ==================================================================================             */
/*     End of Loop ibm nodes */
  }

/* ==================================================================================             */
/*   Total Force on each processor */
  F_px = Cp_px + Cs_px; 
  F_py = Cp_py + Cs_py;
  F_pz = Cp_pz + Cs_pz;

  F_nx = Cp_nx + Cs_nx; 
  F_ny = Cp_ny + Cs_ny;
  F_nz = Cp_nz + Cs_nz;

/* ==================================================================================             */  
/*   Global Sum */
  PetscGlobalSum(&F_px, &F_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_py, &F_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_pz, &F_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&F_nx, &F_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_ny, &F_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_nz, &F_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ap_x, &Ap_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ap_y, &Ap_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ap_z, &Ap_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&An_x, &An_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_y, &An_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_z, &An_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Cp_nx, &Cp_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_ny, &Cp_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_nz, &Cp_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Cp_px, &Cp_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_py, &Cp_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_pz, &Cp_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&M_px, &M_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_py, &M_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_pz, &M_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&M_nx, &M_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_ny, &M_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_nz, &M_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Iap_x, &Iap_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Iap_y, &Iap_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Iap_z, &Iap_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ian_x, &Ian_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ian_y, &Ian_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ian_z, &Ian_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Mdpdn_px, &Mdpdn_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_py, &Mdpdn_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_pz, &Mdpdn_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Mdpdn_nx, &Mdpdn_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_ny, &Mdpdn_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_nz, &Mdpdn_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Pw_px, &Pw_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_py, &Pw_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_pz, &Pw_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Pw_nx, &Pw_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_ny, &Pw_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_nz, &Pw_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ap_t, &Ap_tSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_t, &An_tSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Thrust, &ThrustSum , PETSC_COMM_WORLD);
  PetscGlobalSum(&Drag  , &DragSum   , PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_side,&Pw_sideSum, PETSC_COMM_WORLD);

/* ==================================================================================             */
/*   Scale Check later !!!!! */

  A_xSum = 0.5 * (Ap_xSum + An_xSum);
  A_ySum = 0.5 * (Ap_ySum + An_ySum);
  A_zSum = 0.5 * (Ap_zSum + An_zSum);
  A_tSum = (Ap_tSum + An_tSum);

  F_xSum = F_pxSum + F_nxSum;
  F_ySum = F_pySum + F_nySum;
  F_zSum = F_pzSum + F_nzSum;

  if (!fish /*&& !cop && !wing*/) {
  if (fabs(A_xSum)>1e-6)
    F_xSum=F_xSum/A_xSum*2.;
  if (fabs(A_ySum)>1e-6)
    F_ySum=F_ySum/A_ySum*2.;
  if (fabs(A_zSum)>1e-6)
    F_zSum=F_zSum/A_zSum*2.;
  }


  Cp_xSum = Cp_pxSum + Cp_nxSum;
  Cp_ySum = Cp_pySum + Cp_nySum;
  Cp_zSum = Cp_pzSum + Cp_nzSum;

  if (!fish /*&& !cop && !wing*/) {
  if (fabs(A_xSum)>1e-6)
    Cp_xSum=Cp_xSum/A_xSum*2.;
  if (fabs(A_ySum)>1e-6)
    Cp_ySum=Cp_ySum/A_ySum*2.;
  if (fabs(A_zSum)>1e-6)
    Cp_zSum=Cp_zSum/A_zSum*2.;    
  }

  Ia_xSum = 0.5 * (Iap_xSum - Ian_xSum);
  Ia_ySum = 0.5 * (Iap_ySum - Ian_ySum);
  Ia_zSum = 0.5 * (Iap_zSum - Ian_zSum);

  M_xSum = M_pxSum + M_nxSum;
  M_ySum = M_pySum + M_nySum;
  M_zSum = M_pzSum + M_nzSum;

  Pw_xSum = (Pw_pxSum + Pw_nxSum);/// A_tSum;
  Pw_ySum = (Pw_pySum + Pw_nySum);/// A_tSum;
  Pw_zSum = (Pw_pzSum + Pw_nzSum);/// A_tSum;

  efficiency= Cp_zSum/(Cp_zSum+Pw_xSum);


  Mdpdn_xSum = Mdpdn_pxSum + Mdpdn_nxSum;
  Mdpdn_ySum = Mdpdn_pySum + Mdpdn_nySum;
  Mdpdn_zSum = Mdpdn_pzSum + Mdpdn_nzSum;

/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_xSum=M_xSum/A_zSum*2.; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     M_ySum=M_ySum/A_ySum*2.; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_zSum=M_zSum/A_zSum*2.; */

  A_totSum = Ap_xSum + Ap_ySum + Ap_zSum;

/* ==================================================================================             */
/*   store results in fsi */
  fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum;
  fsi->A_tot = A_totSum;
  fsi->M_x = M_xSum; fsi->M_y = M_ySum; fsi->M_z = M_zSum;
  fsi->Mdpdn_x = Mdpdn_xSum; fsi->Mdpdn_y = Mdpdn_ySum; fsi->Mdpdn_z = Mdpdn_zSum;

/* ==================================================================================             */
/*   output values */
  PetscPrintf(PETSC_COMM_WORLD, "F_x,F_y,F_z SI, %le %le %le Az %le %le Ay %le %le\n",F_xSum,F_ySum,F_zSum,Ap_zSum,An_zSum,Ap_ySum,An_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "M_x,M_y,M_z SI, %le %le %le Ia_x %le %le Ip_y %le %le\n",M_xSum,M_ySum,M_zSum,Iap_xSum,Ian_xSum,Iap_ySum,Ian_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "Mdpdn_x,Mdpdn_y,Mdpdn_z SI, %le %le %le\n",Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum);
  if (fish)
  PetscPrintf(PETSC_COMM_WORLD, "eff,P_sid,PThrust,PDrag, %le %le %le %le %le\n",efficiency,Pw_sideSum,Pw_xSum,ThrustSum,DragSum);

  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti,F_ySum, F_zSum,Cp_ySum,Cp_zSum,F_ySum-Cp_ySum,F_zSum-Cp_zSum,A_ySum,F_xSum, M_xSum); */
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum); */
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum, fsi->M_x_old, fsi->M_x_real, fsi->Mdpdn_x);
    fclose(f);

    //if (fish) {
    sprintf(filen, "Power_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f," %d %le %le %le %le %le %le %le\n",ti, efficiency,Pw_sideSum,Pw_xSum,Pw_ySum,Pw_zSum,F_zSum,A_tSum);
    fclose(f);
    //}

  }

/* ==================================================================================             */
/*   Restore Working arrays */
  DAVecRestoreArray(fda, user->lCent, &cent);
  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lP, &p);
  DAVecRestoreArray(da, user->lNvert, &nvert);

  DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);

  VecDestroy(Coor);

  return(0);
}

//add (toni) for FSI and floating structures 
PetscErrorCode Calc_forces_SI_levelset(FSInfo *fsi,UserCtx *user,
			      IBMNodes *ibm,PetscInt ti, 
			      PetscInt ibi, PetscInt bi)
{
  DA	da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal     rei= 1./user->ren;
  Vec           Coor;
  PetscInt	i, j, k, elmt;
  PetscInt	nv1,nv2,nv3;
  Cmpnts        uelmt;
  PetscReal     sb;
  PetscReal     Tow_ws,Tow_wt,Tow_wn;
  PetscReal     Tow_x, Tow_y, Tow_z;
  Cmpnts        ***coor, ***ucat, ***cent,***ucat_o;
  PetscReal     ***p, ***nvert;
	PetscReal      ***aj; 
  PetscReal     dx,dy,dz, dx1,dy1,dz1;
  PetscReal     dwdz,dwdy,dwdx;
  PetscReal     dvdz,dvdy,dvdx;
  PetscReal     dudz,dudy,dudx;
  PetscReal     Txx,Tyy,Tzz;
  PetscReal     Tzy,Tzx,Tyx;
  PetscReal     nfx,nfy,nfz;
  PetscReal     nsx,nsy,nsz;
  PetscReal     ntx,nty,ntz;
  PetscReal     F_px,F_py,F_pz,Ap_x,Ap_y,Ap_z,Ap_t; //Forces and Area
  PetscReal     F_nx,F_ny,F_nz,An_x,An_y,An_z,An_t; 
  PetscReal     Cp_nx,Cp_ny,Cp_nz; //Pressure Forces - side
  PetscReal     Cs_nx,Cs_ny,Cs_nz; //Surface Forces - side
  PetscReal     Cp_px,Cp_py,Cp_pz; //Pressure Forces + side
  PetscReal     Cs_px,Cs_py,Cs_pz; //Surface Forces + side
  PetscReal     F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force
  PetscReal     F_pxSum,F_pySum,F_pzSum; //Surface Force
  PetscReal     F_nxSum,F_nySum,F_nzSum; //Surface Force
  PetscReal     Ap_xSum,Ap_ySum,Ap_zSum,Ap_tSum; // + side
  PetscReal     An_xSum,An_ySum,An_zSum,An_tSum; // - side
  PetscReal     A_xSum,A_ySum,A_zSum,A_tSum;
  PetscReal     Cp_xSum,Cp_ySum,Cp_zSum; //Pressure Force
  PetscReal     Cp_pxSum,Cp_pySum,Cp_pzSum; //Pressure Force
  PetscReal     Cp_nxSum,Cp_nySum,Cp_nzSum; //Pressure Force
	
  // Moments
  PetscReal      MF_px,MF_py,MF_pz; //Forces for Moment Calc
  PetscReal      MF_nx,MF_ny,MF_nz; //Forces for Moment Calc
  PetscReal      M_nx,M_ny,M_nz;   //Moments
  PetscReal      M_px,M_py,M_pz;   //Moments
  PetscReal      r_x,r_y,r_z;   //Anchor dist
  PetscReal      x, y, z;       //cell coord
  PetscReal      X_c,Y_c,Z_c;   //center of rotation coord
  PetscReal      M_xSum,M_ySum,M_zSum; //Surface Mom on all processors
  PetscReal      M_pxSum,M_pySum,M_pzSum; //Surface Mom on all processors
  PetscReal      M_nxSum,M_nySum,M_nzSum; //Surface Mom on all processors
  PetscReal      Iap_x,Iap_y,Iap_z;
  PetscReal      Ian_x,Ian_y,Ian_z;
  PetscReal      Iap_xSum,Iap_ySum,Iap_zSum; // + side
  PetscReal      Ian_xSum,Ian_ySum,Ian_zSum; // - side
  PetscReal      Ia_xSum,Ia_ySum,Ia_zSum;
  PetscReal      Pw_nx,Pw_ny,Pw_nz;   //Power
  PetscReal      Pw_px,Pw_py,Pw_pz;   //Power
  PetscReal      Pw_nxSum,Pw_nySum,Pw_nzSum;   //Power
  PetscReal      Pw_pxSum,Pw_pySum,Pw_pzSum;   //Power
  PetscReal      Pw_xSum,Pw_ySum,Pw_zSum;

  PetscReal      A_x,A_y,A_z,Atot;
  PetscReal      u_x,u_y,u_z;
  Cmpnts         ***csi,***eta,***zet;
  PetscReal      csi1,csi2,csi3;
  PetscReal      eta1,eta2,eta3;
  PetscReal      zet1,zet2,zet3;
  PetscReal      ***iaj,***jaj,***kaj;

  PetscReal      dr;
  PetscReal      MFdpdn_px,MFdpdn_py,MFdpdn_pz;
  PetscReal      Mdpdn_px,Mdpdn_py,Mdpdn_pz;
  PetscReal      Mdpdn_pxSum,Mdpdn_pySum,Mdpdn_pzSum;
  PetscReal      MFdpdn_nx,MFdpdn_ny,MFdpdn_nz;
  PetscReal      Mdpdn_nx,Mdpdn_ny,Mdpdn_nz;
  PetscReal      Mdpdn_nxSum,Mdpdn_nySum,Mdpdn_nzSum;
  PetscReal      Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum; 

  IBMInfo        *ibminfo;
  IBMListNode    *current;

  // Fish
  PetscReal      pi = 3.141592653589793;
  PetscReal      v_side, dzz, az0, Pw_side, Thrust, Drag;
  PetscReal      Pw_sideSum, ThrustSum, DragSum, efficiency;
  PetscReal      time=user->dt*ti;  
  PetscReal      ampl=0.1;
  PetscReal 	   omega=2*pi*St_exp/(2.*ampl);//2.*pi; //non-dimensionalized w
  PetscReal	     kwave= 0;//2*pi/wavelength;//8.307767;


/* ==================================================================================             */
/*   Init var */
  F_px=0.;F_py=0.;F_pz=0.;
  F_nx=0.;F_ny=0.;F_nz=0.;
  Ap_x=0.;Ap_y=0.;Ap_z=0.;Ap_t=0.;
  An_x=0.;An_y=0.;An_z=0.;An_t=0.;
  Cp_px=0.;Cp_py=0.;Cp_pz=0.;
  Cs_px=0.;Cs_py=0.;Cs_pz=0.;
  Cp_nx=0.;Cp_ny=0.;Cp_nz=0.;
  Cs_nx=0.;Cs_ny=0.;Cs_nz=0.;
	
	double ro;
	double ro_air=rho_air;
	double ro_water=rho_water;

	double Cp_pz2=0.,Cp_nz2=0.;
	double dpdz_=0.,dpdz_Sum=0.,dpdz2_=0.,dpdz2_Sum=0.;
  M_px=0.;M_py=0.;M_pz=0.;
  M_nx=0.;M_ny=0.;M_nz=0.;
  Iap_x=0.;Iap_y=0.;Iap_z=0.;
  Ian_x=0.;Ian_y=0.;Ian_z=0.;
  Pw_px=0.;Pw_py=0.;Pw_pz=0.;
  Pw_nx=0.;Pw_ny=0.;Pw_nz=0.;

  // Fish
  Pw_side=0.; Thrust=0.; Drag=0.;
  

  Mdpdn_px=0.;Mdpdn_py=0.;Mdpdn_pz=0.;
  Mdpdn_nx=0.;Mdpdn_ny=0.;Mdpdn_nz=0.;
  MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
  MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;

	X_c=x_r+fsi->S_new[0], Y_c=y_r+fsi->S_new[2], Z_c=z_r+fsi->S_new[4]; //toni

  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force: %le X_c: %le Y_c: %le Z_c: %le time: %le\n",rei, X_c,Y_c,Z_c,time);

/* ==================================================================================             */
/*   Get Working arrays */
  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
  DAVecGetArray(fda, user->lCent, &cent);
  DAVecGetArray(fda, user->lUcat, &ucat);
	DAVecGetArray(fda, user->lUcat_old,&ucat_o);
  DAVecGetArray(da, user->lP, &p);
  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(fda, user->lCsi, &csi);
  DAVecGetArray(fda, user->lEta, &eta);
  DAVecGetArray(fda, user->lZet, &zet);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);
	DAVecGetArray(da, user->lAj, &aj);
	PetscReal ***mu, ***rho, ***level;
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
		DAVecGetArray(da, user->lLevelset, &level);
	}

/* ==================================================================================             */
/* Loop around all ibm nodes */
  current = user->ibmlist[ibi].head;
  while (current) {
    ibminfo = &current->ibm_intp;
    current = current->next;	
    i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk;
    elmt = ibminfo->cell; // closest ibm element
    sb = ibminfo->d_s;
    //sb = ibminfo->d_i;// |dn| from node to bndry
		int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
		int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
		int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
		double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;

    // normal 
    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    // 1st bi-normal of nf
    nsx=ibm->ns_x[elmt];
    nsy=ibm->ns_y[elmt];
    nsz=ibm->ns_z[elmt];

    // 2nd bi-normal of nf
    ntx=ibm->nt_x[elmt];
    nty=ibm->nt_y[elmt];
    ntz=ibm->nt_z[elmt];

    // nodes of closest ibm elmnt
    nv1=ibm->nv1[elmt];
    nv2=ibm->nv2[elmt];
    nv3=ibm->nv3[elmt];

    //velocity of the closest elmnt
    uelmt.x = (ibm->u[nv1].x+ibm->u[nv2].x+ibm->u[nv3].x)/3.;
    uelmt.y = (ibm->u[nv1].y+ibm->u[nv2].y+ibm->u[nv3].y)/3.;
    uelmt.z = (ibm->u[nv1].z+ibm->u[nv2].z+ibm->u[nv3].z)/3.;
		
		double Ux=uelmt.x, Uy=uelmt.y, Uz=uelmt.z;
		double Ux_old=(ibm->uold[nv1].x+ibm->uold[nv2].x+ibm->uold[nv3].x)/3.;
		double Uy_old=(ibm->uold[nv1].y+ibm->uold[nv2].y+ibm->uold[nv3].y)/3.;		
		double Uz_old=(ibm->uold[nv1].z+ibm->uold[nv2].z+ibm->uold[nv3].z)/3.;

    if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze){
    	rei= 1./user->ren;       
			if(levelset && level[k][j][i]>=0) rei = mu_water;
			if(levelset && level[k][j][i]<0) rei =mu_air;
	
      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;
      MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
      MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;

/* ==================================================================================             */
/*       Shear Stresses (2nd & 1st order) and Shear Force */

      if (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k+1<mz-1 && k-1>0) {
				zet1 = 0.25*(zet[k][j][i].x+zet[k-1][j][i].x)*(kaj[k][j][i]+kaj[k-1][j][i]);
				zet2 = 0.25*(zet[k][j][i].y+zet[k-1][j][i].y)*(kaj[k][j][i]+kaj[k-1][j][i]);
				zet3 = 0.25*(zet[k][j][i].z+zet[k-1][j][i].z)*(kaj[k][j][i]+kaj[k-1][j][i]);
				dwdz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet3;
				dvdz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet3;
				dudz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet3;
				dwdy = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet2;
				dvdy = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet2;
				dudy = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet2;
				dwdx = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet1;
				dvdx = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet1;
				dudx = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet1;
      } 
			else if (nvert[k+1][j][i]<2.5 && k+1<mz-1) {
				zet1 = (zet[k][j][i].x)*kaj[k][j][i];
				zet2 = (zet[k][j][i].y)*kaj[k][j][i];
				zet3 = (zet[k][j][i].z)*kaj[k][j][i];
				dwdz = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet3;
				dvdz = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet3;
				dudz = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet3;
				dwdy = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet2;
				dvdy = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet2;
				dudy = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet2;
				dwdx = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet1;
				dvdx = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet1;
				dudx = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet1;
      } 
			else if (nvert[k-1][j][i]<2.5 && k-1>0){
				zet1 = (zet[k-1][j][i].x)*kaj[k-1][j][i];
				zet2 = (zet[k-1][j][i].y)*kaj[k-1][j][i];
				zet3 = (zet[k-1][j][i].z)*kaj[k-1][j][i];
				dwdz = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet3;
				dvdz = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet3;
				dudz = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet3;
				dwdy = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet2;
				dvdy = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet2;
				dudy = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet2;
				dwdx = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet1;
				dvdx = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet1;
				dudx = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet1;
      } 
			else {
				dwdz = 0.;
				dvdz = 0.;
				dudz = 0.;
				dwdy = 0.;
				dvdy = 0.;
				dudy = 0.;
				dwdx = 0.;
				dvdx = 0.;
				dudx = 0.;
      }
      if (nvert[k][j+1][i]<2.5 && j+1<my-1 && nvert[k][j-1][i]<2.5 && j-1>0) {
				eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*(jaj[k][j][i]+jaj[k][j-1][i]);
				eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*(jaj[k][j][i]+jaj[k][j-1][i]);
				eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*(jaj[k][j][i]+jaj[k][j-1][i]);
				dwdz += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta3;
				dvdz += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta3;
				dudz += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta3;
				dwdy += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta2;
				dvdy += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta2;
				dudy += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta2;
				dwdx += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta1;
				dvdx += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta1;
				dudx += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta1;
      } 
			else if (nvert[k][j+1][i]<2.5 && j+1<my-1) {
				eta1 = eta[k][j][i].x*jaj[k][j][i];
				eta2 = eta[k][j][i].y*jaj[k][j][i];
				eta3 = eta[k][j][i].z*jaj[k][j][i];
				dwdz += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta3;
				dvdz += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta3;
				dudz += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta3;
				dwdy += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta2;
				dvdy += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta2;
				dudy += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta2;
				dwdx += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta1;
				dvdx += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta1;
				dudx += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta1;
      } 
			else if (nvert[k][j-1][i]<2.5 && j-1>0){
				eta1 = eta[k][j-1][i].x*jaj[k][j-1][i];
				eta2 = eta[k][j-1][i].y*jaj[k][j-1][i];
				eta3 = eta[k][j-1][i].z*jaj[k][j-1][i];
				dwdz += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta3;
				dvdz += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta3;
				dudz += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta3;
				dwdy += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta2;
				dvdy += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta2;
				dudy += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta2;
				dwdx += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta1;
				dvdx += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta1;
				dudx += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta1;
      } 
      if (nvert[k][j][i+1]<2.5 && i+1<mx-1 && nvert[k][j][i-1]<2.5 && i-1>0) {
				csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)*(iaj[k][j][i]+iaj[k][j][i-1]);
				csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)*(iaj[k][j][i]+iaj[k][j][i-1]);
				csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)*(iaj[k][j][i]+iaj[k][j][i-1]);
				dwdz += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi3;
				dvdz += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi3;
				dudz += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi3;
				dwdy += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi2;
				dvdy += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi2;
				dudy += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi2;
				dwdx += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi1;
				dvdx += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi1;
				dudx += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi1;
      } 
			else if (nvert[k][j][i+1]<2.5 && i+1<mx-1) {
				csi1 = csi[k][j][i].x*iaj[k][j][i];
				csi2 = csi[k][j][i].y*iaj[k][j][i];
				csi3 = csi[k][j][i].z*iaj[k][j][i];
				dwdz += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi3;
				dvdz += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi3;
				dudz += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi3;
				dwdy += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi2;
				dvdy += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi2;
				dudy += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi2;
				dwdx += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi1;
				dvdx += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi1;
				dudx += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi1;
      } 
			else if (nvert[k][j][i-1]<2.5 && i-1>0){
				csi1 = csi[k][j][i-1].x*iaj[k][j][i-1];
				csi2 = csi[k][j][i-1].y*iaj[k][j][i-1];
				csi3 = csi[k][j][i-1].z*iaj[k][j][i-1];
				dwdz += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi3;
				dvdz += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi3;
				dudz += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi3;
				dwdy += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi2;
				dvdy += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi2;
				dudy += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi2;
				dwdx += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi1;
				dvdx += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi1;
				dudx += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi1;
      }
      
      Tzz = rei * (dwdz + dwdz);
      Tyy = rei * (dvdy + dvdy);
      Txx = rei * (dudx + dudx);
      Tzy = rei * (dwdy + dvdz);
      Tzx = rei * (dwdx + dudz);
      Tyx = rei * (dvdx + dudy);

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      r_x = 0.;
      r_y = 0.;
      r_z = 0.;
/* ==================================================================================             */
			
/* ==================================================================================  
/*       K +  */

      if (nvert[k-1][j][i]>2.5 ){
				x = 0.5*(cent[k][j][i].x + cent[k-1][j][i].x);
				y = 0.5*(cent[k][j][i].y + cent[k-1][j][i].y);
				z = 0.5*(cent[k][j][i].z + cent[k-1][j][i].z);		
				r_x = x-X_c; r_y = y-Y_c; r_z = z-Z_c;
				A_z=zet[k][j][i].z; A_y=zet[k][j][i].y; A_x=zet[k][j][i].x;

				double rho_ = rho[kp1][jp1][ip1]*sk1 + rho[kp2][jp2][ip2]*sk2 + rho[kp3][jp3][ip3]*sk3;
				double dist=cent[k][j][i].z-cent[k-1][j][i].z;
				dist=dist*.5;
				//pressure at k=k-1/2 i.e. cell face between structure cell and ib cell
				double pp=p[k][j][i] + rho_*(Uz-Uz_old)/user->dt*dist+rho_*dist*(-gravity_z);

				Cp_pz += -(pp)*A_z;
				MF_pz += -pp*A_z;
				Ap_z  += A_z;
				Iap_x -=  r_y*A_z;
				Iap_y -= - r_x*A_z;
				MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

				if (MF_pz<0) Thrust += MF_pz;
				else Drag += MF_pz;

				Cp_py += -pp*A_y;
				MF_py += -pp*A_y;
				Ap_y  += A_y;
				Iap_x -=  - r_z*A_y;
				Iap_z -=    r_x*A_y;
				MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_px += -pp*A_x;
				MF_px += -pp*A_x;
				Ap_x  += A_x;
				Iap_z -=   - r_y*A_x;
				Iap_y -= -(- r_z*A_x);
				MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_pz += Tzz*A_z; Cs_py += Tzy*A_z; Cs_px += Tzx*A_z ;
				MF_pz += Tzz*A_z; MF_py += Tzy*A_z; MF_px += Tzx*A_z ;

				if (Tzz*A_z<0) Thrust += Tzz*A_z;
				else Drag += Tzz*A_z;

				Cs_pz += Tzy*A_y ;
				Cs_py += Tyy*A_y ;
				Cs_px += Tyx*A_y ;
					
				MF_pz += Tzy*A_y ;
				MF_py += Tyy*A_y ;
				MF_px += Tyx*A_y ;

				if (Tzy*A_y<0) Thrust += Tzy*A_y;
				else Drag += Tzy*A_y;

				Cs_pz += Tzx*A_x; Cs_py += Tyx*A_x; Cs_px += Txx*A_x;
					
				MF_pz += Tzx*A_x; MF_py += Tyx*A_x; MF_px += Txx*A_x;

				if (Tzx*A_x<0) Thrust += Tzx*A_x;
				else Drag += Tzx*A_x;

				M_px -=   r_y*MF_pz - r_z*MF_py; M_py -= -(r_x*MF_pz - r_z*MF_px); M_pz -=   r_x*MF_py - r_y*MF_px;

				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				Ap_t += Atot;

				u_x = uelmt.x; u_y = uelmt.y; u_z = uelmt.z;

				Pw_px += MF_px * u_x ; Pw_py += MF_py * u_y ; Pw_pz += MF_pz * u_z ;//* Atot;

      }   

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       K -  */
   
      if (nvert[k+1][j][i]>2.5) {
				A_z=zet[k][j][i].z; A_y=zet[k][j][i].y; A_x=zet[k][j][i].x;
				x = 0.5 * (cent[k][j][i].x + cent[k+1][j][i].x); 
				y = 0.5 * (cent[k][j][i].y + cent[k+1][j][i].y); 
				z = 0.5 * (cent[k][j][i].z + cent[k+1][j][i].z);	

				r_x = x-X_c; r_y = y-Y_c; r_z = z-Z_c;

				double dist=cent[k][j][i].z-cent[k+1][j][i].z;
				dist=dist*.5;
				double rho_   =   rho[kp1][jp1][ip1]*sk1   +   rho[kp2][jp2][ip2]*sk2   +   rho[kp3][jp3][ip3]*sk3;
				double pp=p[k][j][i] + rho_*(Uz-Uz_old)/user->dt*dist+rho_*dist*(-gravity_z);
				Cp_pz += (    pp   )*A_z;	

				MF_nz +=  pp*A_z;
				An_z  += A_z;
				Ian_x -=   r_y*A_z;
				Ian_y -= - r_x*A_z;
				MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

				if (MF_nz<0) Thrust += MF_nz;
				else Drag += MF_nz;

				Cp_ny +=  pp*A_y;
				MF_ny +=  pp*A_y;
				An_y  += A_y;
				Ian_x -=  - r_z*A_y;
				Ian_z -=    r_x*A_y;
				MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_nx +=  pp*A_x;
				MF_nx +=  pp*A_x;
				An_x  += A_x;
				Ian_z -=   - r_y*A_x;
				Ian_y -= -(- r_z*A_x);
				MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_nz -= Tzz*A_z; Cs_ny -= Tzy*A_z; Cs_nx -= Tzx*A_z;
				MF_nz -= Tzz*A_z; MF_ny -= Tzy*A_z; MF_nx -= Tzx*A_z;
				
				if (-Tzz*A_z<0) Thrust -= Tzz*A_z;
				else Drag -= Tzz*A_z;

				Cs_nz -= Tzy*A_y; Cs_ny -= Tyy*A_y; Cs_nx -= Tyx*A_y ;
				MF_nz -= Tzy*A_y; MF_ny -= Tyy*A_y; MF_nx -= Tyx*A_y;

				if (-Tzy*A_y<0) Thrust -= Tzy*A_y;
				else 	Drag -= Tzy*A_y;
				
				Cs_nz -= Tzx*A_x; Cs_ny -= Tyx*A_x; Cs_nx -= Txx*A_x;
				MF_nz -= Tzx*A_x; MF_ny -= Tyx*A_x; MF_nx -= Txx*A_x;

				if (-Tzx*A_x<0) Thrust -= Tzx*A_x;
				else 	Drag -= Tzx*A_x;

				M_nx -=   r_y*MF_nz - r_z*MF_ny; M_ny -= -(r_x*MF_nz - r_z*MF_nx); M_nz -=   r_x*MF_ny - r_y*MF_nx;
				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				An_t += Atot;
				u_x = uelmt.x; u_y = uelmt.y; u_z = uelmt.z;
				Pw_nx += MF_nx * u_x ; Pw_ny += MF_ny * u_y ; Pw_nz += MF_nz * u_z ;//* Atot;
      }
     
      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       j +  */

      if (nvert[k][j-1][i]>2.5 ){
				x = 0.5*(cent[k][j][i].x + cent[k][j-1][i].x);
				y = 0.5*(cent[k][j][i].y + cent[k][j-1][i].y);
				z = 0.5*(cent[k][j][i].z + cent[k][j-1][i].z);		
				r_x = x-X_c; r_y = y-Y_c; r_z = z-Z_c;
				A_z=eta[k][j][i].z; A_y=eta[k][j][i].y; A_x=eta[k][j][i].x;

				double rho_ = rho[kp1][jp1][ip1]*sk1 + rho[kp2][jp2][ip2]*sk2 + rho[kp3][jp3][ip3]*sk3;
				double dist=cent[k][j][i].y-cent[k][j-1][i].y;
				dist=dist*.5;
				//pressure at k=k-1/2 i.e. cell face between structure cell and ib cell
				double pp=p[k][j][i] + rho_*(Uy-Uy_old)/user->dt*dist+rho_*dist*(-gravity_y);
	
				Cp_pz += -pp*A_z;
				MF_pz += -pp*A_z;
				Ap_z  += A_z;
				Iap_x -=   r_y*A_z;
				Iap_y -= - r_x*A_z;
				MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

				Cp_py += -pp*A_y;
				MF_py += -pp*A_y;
				Ap_y  += A_y;
				Iap_x -=  - r_z*A_y;
				Iap_z -=    r_x*A_y;
				MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_px += -pp*A_x;
				MF_px += -pp*A_x;
				Ap_x  += A_x;
				Iap_z -=   - r_y*A_x;
				Iap_y -= -(- r_z*A_x);
				MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_pz += Tzz*A_z; Cs_py += Tzy*A_z; Cs_px += Tzx*A_z;
				MF_pz += Tzz*A_z; MF_py += Tzy*A_z; MF_px += Tzx*A_z;

				if (Tzz*A_z<0) Thrust += Tzz*A_z;
				else Drag += Tzz*A_z;

				Cs_pz += Tzy*A_y;	Cs_py += Tyy*A_y;	Cs_px += Tyx*A_y;
				MF_pz += Tzy*A_y;	MF_py += Tyy*A_y;	MF_px += Tyx*A_y;

				if (Tzy*A_y<0) Thrust += Tzy*A_y;
				else Drag += Tzy*A_y;

				Cs_pz += Tzx*A_x; Cs_py += Tyx*A_x; Cs_px += Txx*A_x;
				MF_pz += Tzx*A_x; MF_py += Tyx*A_x; MF_px += Txx*A_x;

				if (Tzx*A_x<0) Thrust += Tzx*A_x;
				else Drag += Tzx*A_x;

				M_px -=   r_y*MF_pz - r_z*MF_py; M_py -= -(r_x*MF_pz - r_z*MF_px); M_pz -=   r_x*MF_py - r_y*MF_px;

				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				Ap_t += Atot;

				u_x = uelmt.x; u_y = uelmt.y; u_z = uelmt.z;

				Pw_px += MF_px * u_x;	Pw_py += MF_py * u_y;	Pw_pz += MF_pz * u_z;//* Atot;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       j -  */

      if (nvert[k][j+1][i]>2.5) {
				x = 0.5*(cent[k][j][i].x + cent[k][j+1][i].x);
				y = 0.5*(cent[k][j][i].y + cent[k][j+1][i].y);
				z = 0.5*(cent[k][j][i].z + cent[k][j+1][i].z);		
				r_x = x-X_c; r_y = y-Y_c; r_z = z-Z_c;
				A_z=eta[k][j][i].z; A_y=eta[k][j][i].y; A_x=eta[k][j][i].x;

				double rho_ = rho[kp1][jp1][ip1]*sk1 + rho[kp2][jp2][ip2]*sk2 + rho[kp3][jp3][ip3]*sk3;
				double dist=cent[k][j][i].y-cent[k][j+1][i].y;
				dist=dist*.5;
				//pressure at k=k-1/2 i.e. cell face between structure cell and ib cell
				double pp=p[k][j][i] + rho_*(Uy-Uy_old)/user->dt*dist+rho_*dist*(-gravity_y);

				Cp_nz +=  pp*A_z;
				MF_nz +=  pp*A_z;
				An_z  += A_z;
				Ian_x -=   r_y*A_z;
				Ian_y -= - r_x*A_z;
				MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

				Cp_ny +=  pp*A_y;
				MF_ny +=  pp*A_y;
				An_y  += A_y;
				Ian_x -=  - r_z*A_y;
				Ian_z -=    r_x*A_y;
				MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_nx +=  pp*A_x;
				MF_nx +=  pp*A_x;
				An_x  += A_x;
				Ian_z -=   - r_y*A_x;
				Ian_y -= -(- r_z*A_x);
				MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_nz -= Tzz*A_z ;	Cs_ny -= Tzy*A_z ;	Cs_nx -= Tzx*A_z ;
				MF_nz -= Tzz*A_z ;	MF_ny -= Tzy*A_z ;	MF_nx -= Tzx*A_z ;
	
				if (-Tzz*A_z<0) Thrust -= Tzz*A_z;
				else Drag -= Tzz*A_z;

				Cs_nz -= Tzy*A_y ;	Cs_ny -= Tyy*A_y ;	Cs_nx -= Tyx*A_y ;
				MF_nz -= Tzy*A_y ;	MF_ny -= Tyy*A_y ;	MF_nx -= Tyx*A_y ;

				if (-Tzy*A_y<0) Thrust -= Tzy*A_y;
				else Drag -= Tzy*A_y;
				
				Cs_nz -= Tzx*A_x;	Cs_ny -= Tyx*A_x;	Cs_nx -= Txx*A_x;
				MF_nz -= Tzx*A_x;	MF_ny -= Tyx*A_x;	MF_nx -= Txx*A_x;

				if (-Tzx*A_x<0) Thrust -= Tzx*A_x;
				else Drag -= Tzx*A_x;

				M_nx -=   r_y*MF_nz - r_z*MF_ny;	M_ny -= -(r_x*MF_nz - r_z*MF_nx);	M_nz -=   r_x*MF_ny - r_y*MF_nx;

				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				An_t += Atot;

				u_x = uelmt.x;	u_y = uelmt.y;	u_z = uelmt.z;
				Pw_nx += MF_nx * u_x ;	Pw_ny += MF_ny * u_y ;	Pw_nz += MF_nz * u_z ;//* Atot;
      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       i +  */

      if (nvert[k][j][i-1]>2.5){
			
				A_z=csi[k][j][i].z;	A_y=csi[k][j][i].y;	A_x=csi[k][j][i].x;
				x = 0.5*(cent[k][j][i].x + cent[k][j][i-1].x);
				y = 0.5*(cent[k][j][i].y + cent[k][j][i-1].y);
				z = 0.5*(cent[k][j][i].z + cent[k][j][i-1].z);					
				
				r_x = x-X_c;	r_y = y-Y_c;	r_z = z-Z_c;

				
				double rho_ = rho[kp1][jp1][ip1]*sk1 + rho[kp2][jp2][ip2]*sk2 + rho[kp3][jp3][ip3]*sk3;
				double dist=cent[k][j][i].x-cent[k][j][i-1].x;
				dist=dist*.5;
				//pressure at k=k-1/2 i.e. cell face between structure cell and ib cell
				double pp=p[k][j][i] + rho_*(Ux-Ux_old)/user->dt*dist;				
				
				Cp_pz += -pp*A_z;
				MF_pz += -pp*A_z;
				Ap_z  += A_z;
				Iap_x -=   r_y*A_z;
				Iap_y -= - r_x*A_z;
				MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

				Cp_py += -pp*A_y;
				MF_py += -pp*A_y;
				Ap_y  += A_y;
				Iap_x -=  - r_z*A_y;
				Iap_z -=    r_x*A_y;
				MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_px += -pp*A_x;
				MF_px += -pp*A_x;
				Ap_x  += A_x;
				Iap_z -=   - r_y*A_x;
				Iap_y -= -(- r_z*A_x);
				MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_pz += Tzz*A_z ;	Cs_py += Tzy*A_z ;	Cs_px += Tzx*A_z ;
				MF_pz += Tzz*A_z ;	MF_py += Tzy*A_z ;	MF_px += Tzx*A_z ;

				if (Tzz*A_z<0) Thrust += Tzz*A_z;
				else Drag += Tzz*A_z;

				Cs_pz += Tzy*A_y ;	Cs_py += Tyy*A_y ;	Cs_px += Tyx*A_y ;
				MF_pz += Tzy*A_y ;	MF_py += Tyy*A_y ;	MF_px += Tyx*A_y ;

				if (Tzy*A_y<0) Thrust += Tzy*A_y;
				else Drag += Tzy*A_y;

				Cs_pz += Tzx*A_x;	Cs_py += Tyx*A_x;	Cs_px += Txx*A_x;
				MF_pz += Tzx*A_x;	MF_py += Tyx*A_x;	MF_px += Txx*A_x;

				if (Tzx*A_x<0) Thrust += Tzx*A_x;
				else Drag += Tzx*A_x;

				M_px -=   r_y*MF_pz - r_z*MF_py;	M_py -= -(r_x*MF_pz - r_z*MF_px);	M_pz -=   r_x*MF_py - r_y*MF_px;

				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				Ap_t += Atot;

				u_x = uelmt.x;	u_y = uelmt.y;	u_z = uelmt.z;

				Pw_px += MF_px * u_x ; Pw_py += MF_py * u_y ; Pw_pz += MF_pz * u_z ;//* Atot;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

/* ==================================================================================             */
/*       i -  */

      if  (nvert[k][j][i+1]>2.5) {
			
			
				A_z=csi[k][j][i].z;	A_y=csi[k][j][i].y;	A_x=csi[k][j][i].x;
				x = 0.5 * (cent[k][j][i].x + cent[k][j][i+1].x); 
				y = 0.5 * (cent[k][j][i].y + cent[k][j][i+1].y); 
				z = 0.5 * (cent[k][j][i].z + cent[k][j][i+1].z);	

				r_x = x-X_c; r_y = y-Y_c; r_z = z-Z_c;

				double dist=cent[k][j][i].x-cent[k][j][i+1].x;
				dist=dist*.5;
				double rho_   =   rho[kp1][jp1][ip1]*sk1   +   rho[kp2][jp2][ip2]*sk2   +   rho[kp3][jp3][ip3]*sk3;
				double pp=p[k][j][i] + rho_*(Ux-Ux_old)/user->dt*dist;			
			

				Cp_nz +=  pp*A_z;
				MF_nz +=  pp*A_z;
				An_z  += A_z;
				Ian_x -=   r_y*A_z;
				Ian_y -= - r_x*A_z;
				MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

				Cp_ny +=  pp*A_y;
				MF_ny +=  pp*A_y;
				An_y  += A_y;
				Ian_x -=  - r_z*A_y;
				Ian_z -=    r_x*A_y;
				MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

				Cp_nx +=  pp*A_x;
				MF_nx +=  pp*A_x;
				An_x  += A_x;
				Ian_z -=   - r_y*A_x;
				Ian_y -= -(- r_z*A_x);
				MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

				Cs_nz -= Tzz*A_z ;	Cs_ny -= Tzy*A_z ;	Cs_nx -= Tzx*A_z ;
				MF_nz -= Tzz*A_z ;	MF_ny -= Tzy*A_z ;	MF_nx -= Tzx*A_z ;	
				
				if (-Tzz*A_z<0) Thrust -= Tzz*A_z;
				else Drag -= Tzz*A_z;

				Cs_nz -= Tzy*A_y ;	Cs_ny -= Tyy*A_y ;	Cs_nx -= Tyx*A_y ;
				MF_nz -= Tzy*A_y ;	MF_ny -= Tyy*A_y ;	MF_nx -= Tyx*A_y ;

				if (-Tzy*A_y<0) Thrust -= Tzy*A_y;
				else Drag -= Tzy*A_y;
				
				Cs_nz -= Tzx*A_x;	Cs_ny -= Tyx*A_x;	Cs_nx -= Txx*A_x;
				MF_nz -= Tzx*A_x;	MF_ny -= Tyx*A_x;	MF_nx -= Txx*A_x;

				if (-Tzx*A_x<0) Thrust -= Tzx*A_x;
				else Drag -= Tzx*A_x;
				
				M_nx -=   r_y*MF_nz - r_z*MF_ny; M_ny -= -(r_x*MF_nz - r_z*MF_nx); M_nz -=   r_x*MF_ny - r_y*MF_nx;

				Atot = sqrt(A_x*A_x+A_y*A_y+A_z*A_z);
				An_t += Atot;

				u_x = uelmt.x;	u_y = uelmt.y;	u_z = uelmt.z;

				Pw_nx += MF_nx * u_x ;	Pw_ny += MF_ny * u_y ;	Pw_nz += MF_nz * u_z ;//* Atot;
      }

      if (fish) {
				z = cent[k][j][i].z;
				dzz = z - 1.5;       
				az0=(ampl-0.02)*dzz+0.02;	
				v_side=az0*omega*cos(kwave*dzz-omega*time);
				
				Pw_side +=(MF_px+MF_nx)*v_side;
      }

      Mdpdn_px -=   r_y*MFdpdn_pz - r_z*MFdpdn_py;//
      Mdpdn_py -= -(r_x*MFdpdn_pz - r_z*MFdpdn_px);
      Mdpdn_pz -=   r_x*MFdpdn_py - r_y*MFdpdn_px;

      Mdpdn_nx -=   r_y*MFdpdn_nz - r_z*MFdpdn_ny;//
      Mdpdn_ny -= -(r_x*MFdpdn_nz - r_z*MFdpdn_nx);
      Mdpdn_nz -=   r_x*MFdpdn_ny - r_y*MFdpdn_nx;
    
			if(level[kp1][jp1][ip1]>0 || level[kp2][jp2][ip2]>0 || level[kp3][jp3][ip3]>0) ro=ro_water;
			else ro=ro_air;
			


    } // if ibm node in CPU
/* ==================================================================================             */
/*     End of Loop ibm nodes */
  }
 


  	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
		DAVecRestoreArray(da, user->lLevelset, &level);
	}
/* ==================================================================================             */
/*   Total Force on each processor */
  F_px = Cp_px + Cs_px; 
  F_py = Cp_py + Cs_py;
  F_pz = Cp_pz + Cs_pz;

  F_nx = Cp_nx + Cs_nx; 
  F_ny = Cp_ny + Cs_ny;
  F_nz = Cp_nz + Cs_nz;

/* ==================================================================================             */  
/*   Global Sum */
  PetscGlobalSum(&F_px, &F_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_py, &F_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_pz, &F_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&F_nx, &F_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_ny, &F_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&F_nz, &F_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ap_x, &Ap_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ap_y, &Ap_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ap_z, &Ap_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&An_x, &An_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_y, &An_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_z, &An_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Cp_nx, &Cp_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_ny, &Cp_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_nz, &Cp_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Cp_px, &Cp_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_py, &Cp_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Cp_pz, &Cp_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&M_px, &M_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_py, &M_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_pz, &M_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&M_nx, &M_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_ny, &M_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&M_nz, &M_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Iap_x, &Iap_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Iap_y, &Iap_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Iap_z, &Iap_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ian_x, &Ian_xSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ian_y, &Ian_ySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Ian_z, &Ian_zSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Mdpdn_px, &Mdpdn_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_py, &Mdpdn_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_pz, &Mdpdn_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Mdpdn_nx, &Mdpdn_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_ny, &Mdpdn_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Mdpdn_nz, &Mdpdn_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Pw_px, &Pw_pxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_py, &Pw_pySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_pz, &Pw_pzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Pw_nx, &Pw_nxSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_ny, &Pw_nySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_nz, &Pw_nzSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Ap_t, &Ap_tSum, PETSC_COMM_WORLD);
  PetscGlobalSum(&An_t, &An_tSum, PETSC_COMM_WORLD);

  PetscGlobalSum(&Thrust, &ThrustSum , PETSC_COMM_WORLD);
  PetscGlobalSum(&Drag  , &DragSum   , PETSC_COMM_WORLD);
  PetscGlobalSum(&Pw_side,&Pw_sideSum, PETSC_COMM_WORLD);


	
/* ==================================================================================             */
/*   Scale Check later !!!!! */

  A_xSum = 0.5 * (Ap_xSum + An_xSum);
  A_ySum = 0.5 * (Ap_ySum + An_ySum);
  A_zSum = 0.5 * (Ap_zSum + An_zSum);
  A_tSum = (Ap_tSum + An_tSum);

  F_xSum = F_pxSum + F_nxSum;
  F_ySum = F_pySum + F_nySum;
  F_zSum = F_pzSum + F_nzSum;

  if (!levelset /*&& !cop && !wing*/) {
  if (fabs(A_xSum)>1e-6)
    F_xSum=F_xSum/A_xSum*2.;
  if (fabs(A_ySum)>1e-6)
    F_ySum=F_ySum/A_ySum*2.;
  if (fabs(A_zSum)>1e-6)
    F_zSum=F_zSum/A_zSum*2.;
  }


  Cp_xSum = Cp_pxSum + Cp_nxSum;
  Cp_ySum = Cp_pySum + Cp_nySum;
  Cp_zSum = Cp_pzSum + Cp_nzSum;

  if (!levelset /*&& !cop && !wing*/) {
  //if (fabs(A_xSum)>1e-6)
    //Cp_xSum=Cp_xSum/A_xSum*2.;
  //if (fabs(A_ySum)>1e-6)
    //Cp_ySum=Cp_ySum/A_ySum*2.;
  //if (fabs(A_zSum)>1e-6)
    //Cp_zSum=Cp_zSum/A_zSum*2.;    
  }

  Ia_xSum = 0.5 * (Iap_xSum - Ian_xSum);
  Ia_ySum = 0.5 * (Iap_ySum - Ian_ySum);
  Ia_zSum = 0.5 * (Iap_zSum - Ian_zSum);

  M_xSum = M_pxSum + M_nxSum;
  M_ySum = M_pySum + M_nySum;
  M_zSum = M_pzSum + M_nzSum;

  Pw_xSum = (Pw_pxSum + Pw_nxSum);/// A_tSum;
  Pw_ySum = (Pw_pySum + Pw_nySum);/// A_tSum;
  Pw_zSum = (Pw_pzSum + Pw_nzSum);/// A_tSum;

  efficiency= Cp_zSum/(Cp_zSum+Pw_xSum);


  Mdpdn_xSum = Mdpdn_pxSum + Mdpdn_nxSum;
  Mdpdn_ySum = Mdpdn_pySum + Mdpdn_nySum;
  Mdpdn_zSum = Mdpdn_pzSum + Mdpdn_nzSum;

/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_xSum=M_xSum/A_zSum*2.; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     M_ySum=M_ySum/A_ySum*2.; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_zSum=M_zSum/A_zSum*2.; */

  A_totSum = Ap_xSum + Ap_ySum + Ap_zSum;

/* ==================================================================================             */
/*   store results in fsi */
  fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum; 
  fsi->A_tot = A_totSum;
  fsi->M_x = M_xSum; fsi->M_y = M_ySum; fsi->M_z = M_zSum;
  fsi->Mdpdn_x = Mdpdn_xSum; fsi->Mdpdn_y = Mdpdn_ySum; fsi->Mdpdn_z = Mdpdn_zSum;
	
	
/* ==================================================================================             */
/*   output values */
  PetscPrintf(PETSC_COMM_WORLD, "F_x,F_y,F_z SI, %le %le %le Az %le %le Ay %le %le\n",F_xSum,F_ySum,F_zSum,Ap_zSum,An_zSum,Ap_ySum,An_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "M_x,M_y,M_z SI, %le %le %le Ia_x %le %le Ip_y %le %le\n",M_xSum,M_ySum,M_zSum,Iap_xSum,Ian_xSum,Iap_ySum,Ian_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "Mdpdn_x,Mdpdn_y,Mdpdn_z SI, %le %le %le\n",Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum);
  if (fish)
  PetscPrintf(PETSC_COMM_WORLD, "eff,P_sid,PThrust,PDrag, %le %le %le %le %le\n",efficiency,Pw_sideSum,Pw_xSum,ThrustSum,DragSum);

  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum,F_zSum ,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti,F_ySum, F_zSum,Cp_ySum,Cp_zSum,F_ySum-Cp_ySum,F_zSum-Cp_zSum,A_ySum,F_xSum, M_xSum); */
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum); */
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum, fsi->M_x_old, fsi->M_x_real, fsi->Mdpdn_x);
    fclose(f);

    if (fish) {
    sprintf(filen, "Power_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f," %d %le %le %le %le %le %le %le\n",ti, efficiency,Pw_sideSum,Pw_xSum,Pw_ySum,Pw_zSum,F_zSum,A_tSum);
    fclose(f);
    }
  }

/* ==================================================================================             */
/*   Restore Working arrays */
  DAVecRestoreArray(fda, user->lCent, &cent);
  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
	DAVecRestoreArray(fda, user->lUcat_old,&ucat_o);
  DAVecRestoreArray(da, user->lP, &p);
  DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);
  DAVecRestoreArray(da, user->lAj, &aj);	
  VecDestroy(Coor);
  return(0);
}

//end (toni)
