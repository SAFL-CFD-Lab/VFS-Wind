/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "petsc.h"
#include "variables.h"

#define EPSILON 0.00000001
#define CROSS(dest, v1, v2) \
	dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
	dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
	dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define SUB(dest, v1, v2) \
	dest[0] = v1[0] - v2[0]; \
	dest[1] = v1[1] - v2[1]; \
	dest[2] = v1[2] - v2[2];

PetscInt intsect_triangle(PetscReal orig[3], PetscReal dir[3],
			  PetscReal vert0[3], PetscReal vert1[3], PetscReal vert2[3],
			  PetscReal *t, PetscReal *u, PetscReal *v)
{
  PetscReal edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  PetscReal det, inv_det;

  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);
  
  CROSS(pvec, dir, edge2);
  
  det = DOT(edge1, pvec);

  if (det > -EPSILON && det < EPSILON)
    return 0;
  inv_det = 1.0 / det;

  SUB(tvec, orig, vert0);
  
  *u = DOT(tvec, pvec) * inv_det;

  if (*u < 0.0 || *u > 1.0)
    return 0;

  CROSS(qvec, tvec, edge1);

  *v = DOT(dir, qvec) * inv_det;

  if (*v < 0.0 || *u + *v > 1.0)
    return 0;

  *t = DOT(edge2, qvec) * inv_det;

  return 1;
}


PetscInt ISSameSide2D(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3)
     /* Check whether 2D point p is located on the same side of line p1p2
	with point p3. Returns:
	-1	different side
	1	same side (including the case when p is located
	        right on the line)
	If p and p3 is located on the same side to line p1p2, then
	the (p-p1) X (p2-p1) and (p3-p1) X (p2-p1) should have the same sign
     */
{
  PetscReal t1, t2, t3;
  PetscReal	epsilon = 1.e-10;

  PetscReal A, B, C;
  A = p2.y - p1.y;
  B = -(p2.x - p1.x);
  C = (p2.x - p1.x) * p1.y - (p2.y - p1.y) * p1.x;

  t3 = fabs(A * p.x + B * p.y + C) / sqrt(A*A + B*B);

/*   if (t3<1.e-3) return(1); */
  if (t3 < 1.e-3) {
    t1 = A * p.x + B * p.y + C;
    t2 = A * p3.x + B * p3.y + C;
    //    if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", t1, t2, t3, A, B, C);
  }
  else {
    t1 = (p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x);
    t2 = (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x);
  }

  //!!!!!!!!!!!!1 Change t1, t2 & lt !!!!!!!
  t1 = (p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x);
  t2 = (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x);
  PetscReal lt;
  lt = sqrt((p1.x - p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
  //  if(flagprint) PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", p1.x, p2.x, p3.x, p1.y, p2.y, p3.y);
  //if (fabs(t1) < epsilon) { // Point is located along the line of p1p2
  if (fabs(t1/lt) < epsilon) { // Point is located along the line of p1p2
    return(1);
  }
  // End of change !!!!!!!!!!!!!1

  if (t1 > 0) {
    if (t2 > 0) return (1); // same side
    else return(-1);  // not
  }
  else {
    if (t2 < 0) return(1); // same side
    else return(-1);
  }
}

PetscInt ISInsideTriangle2D(Cpt2D p, Cpt2D pa, Cpt2D pb, Cpt2D pc)
{
  // Check if point p and p3 is located on the same side of line p1p2
  PetscInt 	ls;

  ls = ISSameSide2D(p, pa, pb, pc);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "aaa, %d\n", ls);
  if (ls < 0) {
    return (ls);
  }
  ls = ISSameSide2D(p, pb, pc, pa);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "bbb, %d\n", ls);
  if (ls < 0) {
    return (ls);
  }
  ls = ISSameSide2D(p, pc, pa, pb);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "ccc, %d\n", ls);
  if (ls <0) {
    return(ls);
  }
  return (ls);
}

PetscInt ISPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3,
			   PetscReal nfx, PetscReal nfy, PetscReal nfz) 
{
  PetscInt flag;
  Cpt2D	pj, pj1, pj2, pj3;
  if (fabs(nfz) >= fabs(nfx) && fabs(nfz) >= fabs(nfy) ) {
    pj.x = p.x; pj.y = p.y;
    pj1.x = p1.x; pj1.y = p1.y;
    pj2.x = p2.x; pj2.y = p2.y;
    pj3.x = p3.x; pj3.y = p3.y;
  }
  else if (fabs(nfx) >= fabs(nfy) && fabs(nfx) >= fabs(nfz)) {
    pj.x = p.z; pj.y = p.y;
    pj1.x = p1.z; pj1.y = p1.y;
    pj2.x = p2.z; pj2.y = p2.y;
    pj3.x = p3.z; pj3.y = p3.y;
  }
  else {
    pj.x = p.x; pj.y = p.z;
    pj1.x = p1.x; pj1.y = p1.z;
    pj2.x = p2.x; pj2.y = p2.z;
    pj3.x = p3.x; pj3.y = p3.z;
  }
  flag = ISInsideTriangle2D(pj, pj1, pj2, pj3);
  //  if (flag > 0)
  //  PetscPrintf(PETSC_COMM_WORLD, "%d %e %e %e dddd", flag, nfx, nfy, nfz);
  return(flag);
}

PetscErrorCode Dis_P_Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d)
{
  PetscReal	dmin;
  PetscReal	dx21, dy21, dz21, dx31, dy31, dz31, t;

  dx21 = p2.x - p1.x; dy21 = p2.y - p1.y; dz21 = p2.z - p1.z;
  dx31 = p.x  - p1.x; dy31 = p.y  - p1.y; dz31 = p.z  - p1.z;

  t = (dx31 * dx21 + dy31 * dy21 + dz31 * dz21) / (dx21*dx21 + dy21*dy21 +
						   dz21 * dz21);
  if (t<0) { // The closet point is p1
    po->x = p1.x; po->y = p1.y; po->z = p1.z;
    *d = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
  }
  else if (t>1) { // The closet point is p2
    po->x = p2.x; po->y = p2.y; po->z = p2.z;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
	     (p.z - po->z) * (p.z - po->z));
  }
  else { // The closet point lies between p1 & p2
    po->x = p1.x + t * dx21; po->y = p1.y + t * dy21; po->z = p1.z + t*dz21;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
	     (p.z - po->z) * (p.z - po->z));
  }
  return(0);
}

/* Temp solution to get the value on surfaces */
PetscErrorCode triangle_intp2(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo->cs1 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo->cs2 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo->cs3 = 1. - ibminfo->cs1 - ibminfo->cs2;

  if (ibminfo->cs1 != ibminfo->cs1) {
    PetscPrintf(PETSC_COMM_SELF, "Aa %e %e %e %e %e\n", a, x13, y23, x23, y13);
    PetscPrintf(PETSC_COMM_SELF, "XS %e %e %e %e %e %e\n", p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
  }
  return 0;
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);
  
}

PetscErrorCode triangle_intp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo)
{
	
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo->cr1 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo->cr2 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo->cr3 = 1. - ibminfo->cr1 - ibminfo->cr2;
  return 0;
}

// seokkoo

double tri_area(double x1, double y1, double z1, 
			double x2, double y2, double z2,
			double x3, double y3, double z3)
{
	double xA, yA, zA;
	double xB, yB, zB;
	
	xA = x2 - x1;
	yA = y2 - y1;
	zA = z2 - z1;
	
	xB = x3 - x1;
	yB = y3 - y1;
	zB = z3 - z1;
	
	/*
		i	j	k
		xA	yA	zA
		xB	yB	zB
	*/
	double p = (yA*zB - zA*yB);
	double q = - (xA*zB - zA*xB);
	double r = (xA*yB - yA*xB);
	return sqrt ( p*p + q*q + r*r ) * 0.5;
};

// seokkoo
PetscErrorCode triangle_intp3D(double x, double y, double z, 
						double x1, double y1, double z1, 
						double x2, double y2, double z2,
						double x3, double y3, double z3, IBMInfo *ibminfo)
{
	double Aa = tri_area(x, y, z, x1, y1, z1, x2, y2, z2);
	double Ab = tri_area(x, y, z, x2, y2, z2, x3, y3, z3);
	double Ac = tri_area(x, y, z, x3, y3, z3, x1, y1, z1);
	
	double w1 = Ab;
	double w2 = Ac;
	double w3 = Aa;
	
	double wsum = w1 + w2 + w3;
	w1 /= wsum;
	w2 /= wsum;
	w3 /= wsum;
	
	ibminfo->cr1 = w1;
	ibminfo->cr2 = w2;
	ibminfo->cr3 = w3;
}



PetscErrorCode triangle_intpp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			 IBMInfo *ibminfo)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;
  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a =x13 * y23 - x23 * y13;
  ibminfo->cr11 = (y23 * xp3 - x23 * yp3) / a;
  ibminfo->cr22 = (-y13 * xp3 + x13 * yp3) / a;
  ibminfo->cr33 = 1. - ibminfo->cr11 - ibminfo->cr22;

  if (ibminfo->cr1 < -1.e-6 || ibminfo->cr1 > 1.+1.e-6) {
    x13 = 0.;
    y13 = 0.;
  }
  //  if (fabs(a)<1.-5)
  //  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e\n", y23*xp3, x23*yp3, ibminfo[number].cr1, ibminfo[number].cr2);

  return 0;
}

PetscTruth ISLineTriangleIntp(Cmpnts p1, Cmpnts p2, IBMNodes *ibm, PetscInt ln_v)
{
  PetscInt cutthrough;

 /*  PetscInt ic1, jc1, kc1, ic2, jc2, kc2, ics, ice, jcs, jce, kcs, kce; */

/*   PetscInt i, j, k; */

/*   PetscInt celln; */
  PetscInt n1e, n2e, n3e;
  
  PetscReal *x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
  PetscReal *nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;

  PetscReal xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3;

  PetscReal nfx, nfy, nfz, nftx, nfty, nftz;

  PetscReal tf1, tf2, dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3;

  PetscReal d;
  Cmpnts    pint, pa, pb, pc;

  n1e = ibm->nv1[ln_v]; n2e = ibm->nv2[ln_v]; n3e = ibm->nv3[ln_v];
  
  xf1 = x_bp[n1e]; yf1 = y_bp[n1e]; zf1 = z_bp[n1e];
  xf2 = x_bp[n2e]; yf2 = y_bp[n2e]; zf2 = z_bp[n2e];
  xf3 = x_bp[n3e]; yf3 = y_bp[n3e]; zf3 = z_bp[n3e];

  nfx = nf_x[ln_v]; nfy = nf_y[ln_v]; nfz = nf_z[ln_v];
  
  tf1 = (p1.x - xf1) * nfx + (p1.y - yf1) * nfy + (p1.z - zf1) * nfz;
  tf2 = (p2.x - xf1) * nfx + (p2.y - yf1) * nfy + (p2.z - zf1) * nfz;
  
  if (tf1 * tf2 < 0) { // if p1 & p2 on different sides of triangle
    
    dx1 = p1.x - xf1; dy1 = p1.y - yf1; dz1 = p1.z - zf1;
    
    nftx = p2.x - p1.x; nfty = p2.y - p1.y; nftz = p2.z - p1.z;
    
    dx2 = xf2 - xf1; dy2 = yf2 - yf1; dz2 = zf2 - zf1;
    dx3 = xf3 - xf1; dy3 = yf3 - yf1; dz3 = zf3 - zf1;
    
    d = - (dx1 * (dy2 * dz3 - dz2 * dy3) - 
	   dy1 * (dx2 * dz3 - dz2 * dx3) + 
	   dz1 * (dx2 * dy3 - dy2 * dx3)) / 
         (nftx * (dy2 * dz3 - dz2 * dy3) - 
	  nfty * (dx2 * dz3 - dz2 * dx3) + 
	  nftz * (dx2 * dy3 - dy2 * dx3));

    pint.x = p1.x + d * nftx;
    pint.y = p1.y + d * nfty;
    pint.z = p1.z + d * nftz;

    if (d > 0 && d < 1.) {
      pa.x = xf1; pa.y = yf1; pa.z = zf1;
      pb.x = xf2; pb.y = yf2; pb.z = zf2;
      pc.x = xf3; pc.y = yf3; pc.z = zf3;
      
      cutthrough = ISPointInTriangle(pint, pa, pb, pc,nfx, nfy, nfz);
      
      if (cutthrough == 1) {
	return (PETSC_TRUE);
      }
    }
    
  }
  
  return (PETSC_FALSE);
}
