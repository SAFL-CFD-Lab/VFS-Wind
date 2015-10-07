/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

/*  Start and End nodes on each part of cop */
//PetscInt L_dim;

// legs
PetscInt n1_min,n1_max;
PetscInt n2_min,n2_max;
PetscInt n3_min,n3_max;
PetscInt n4_min,n4_max;
PetscInt n5_min,n5_max;
PetscInt n6_min,n6_max;
PetscInt n7_min,n7_max;
PetscInt n8_min,n8_max;
PetscInt n9_min,n9_max;
PetscInt n10_min,n10_max;
PetscInt n11_min,n11_max;
PetscInt n12_min,n12_max;

// antene
PetscInt n_antene_mn,n_antene_mx;
PetscInt n_antene_0,n_antene_1,n_antene_2,n_antene_3;
PetscReal x_antene_2,z_antene_0,z_antene_1;

// tail 
PetscInt n_tail_0;
PetscReal x_tail_0;

// bot
PetscInt n_bot_0,n_bot_1,n_bot_2;
PetscReal z_bot_0,z_bot_1,z_bot_2;

// points of rotation
PetscInt n_leg_1_rot,n_leg_2_rot,n_leg_3_rot;
PetscInt n_leg_4_rot,n_leg_5_rot,n_leg_6_rot;
PetscInt n_leg_7_rot,n_leg_8_rot,n_leg_9_rot;
PetscInt n_leg_10_rot,n_leg_11_rot,n_leg_12_rot;

// points of ratation coor
PetscReal sx_tail,sy_tail,sz_tail;
PetscReal sx_antene,sy_antene,sz_antene;
PetscReal sx_l1,sy_l1,sz_l1;
PetscReal sx_l2,sy_l2,sz_l2;
PetscReal sx_l3,sy_l3,sz_l3;
PetscReal sx_l4,sy_l4,sz_l4;
PetscReal sx_l5,sy_l5,sz_l5;
PetscReal sx_l6,sy_l6,sz_l6;
PetscReal sx_l7,sy_l7,sz_l7;
PetscReal sx_l8,sy_l8,sz_l8;
PetscReal sx_l9,sy_l9,sz_l9;
PetscReal sx_l10,sy_l10,sz_l10;
PetscReal sx_l11,sy_l11,sz_l11;
PetscReal sx_l12,sy_l12,sz_l12;

/*  total number of nodes on each part of cop */
PetscInt n_tail;
PetscInt n_antene;
PetscInt n_leg1,n_leg2,n_leg3;
PetscInt n_leg4,n_leg5,n_leg6;
PetscInt n_leg7,n_leg8,n_leg9;
PetscInt n_leg10,n_leg11,n_leg12;

/*  Node aray of each part of cop */
PetscInt *nv_tail, *nv_antene;
PetscInt *nv_leg1, *nv_leg2, *nv_leg3;
PetscInt *nv_leg4, *nv_leg5, *nv_leg6;
PetscInt *nv_leg7, *nv_leg8, *nv_leg9;
PetscInt *nv_leg10, *nv_leg11, *nv_leg12;


/* ==================================================================================             */
/*  time variable */

PetscReal nfrq,T_period,N_period;

// tail
PetscReal t1_t,t2_t,t3_t,t4_t;
// antene
PetscReal t1_a,t2_a,t3_a,t4_a;
// legs
PetscReal t1_l12,t2_l12,t3_l12,t4_l12;
PetscReal t1_l34,t2_l34,t3_l34,t4_l34;
PetscReal t1_l56,t2_l56,t3_l56,t4_l56;
PetscReal t1_l78,t2_l78,t3_l78,t4_l78;
PetscReal t1_l910,t2_l910,t3_l910,t4_l910;
PetscReal t1_l1112,t2_l1112,t3_l1112,t4_l1112;

PetscReal tet_t, tet_a;
PetscReal tet_l12,tet_l34,tet_l56,tet_l78,tet_l910,tet_l1112;
PetscReal tet_tmin,tet_tmax;
PetscReal tet_amin,tet_amax;
PetscReal tet_l12min,tet_l34min,tet_l56min,tet_l78min,tet_l910min,tet_l1112min;
PetscReal tet_l12max,tet_l34max,tet_l56max,tet_l78max,tet_l910max,tet_l1112max;

