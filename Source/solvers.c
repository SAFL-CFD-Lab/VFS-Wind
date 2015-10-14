/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
extern PetscInt   block_number, ti, tiout, NumberOfBodies;
extern PetscInt   movefsi, rotatefsi, immersed, STRONG_COUPLING;
extern PetscInt   cop, regime, fish, MHV;
extern PetscReal  max_angle, Flux_in;
extern int averaging, les, wallfunction, freesurface, rans;
extern int tistart, time_marching, poisson;
extern void free_surafe_BC(UserCtx *user);
extern void Update_Velocity_by_Gravity(UserCtx *user);
extern void pseudo_periodic_BC(UserCtx *user);
extern void K_Omega_Set_Constant(UserCtx *user);
extern void Solve_K_Omega(UserCtx *user);
extern void K_Omega_IC(UserCtx *user);
void Init_LevelSet_Vectors(UserCtx *user);
void Destroy_LevelSet_Vectors(UserCtx *user);
void Distance_Function_IC(UserCtx *user);
void Solve_Distance(UserCtx *user);
extern void Compute_Distance_Function(UserCtx *user);

PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, 
			    FSInfo *fsi, PetscInt itr_sc,
			    PetscInt tistart, 
			    PetscTruth *DoSCLoop)
{
  PetscReal     dS_sc, dS_MIN=1e-5, dSmax;
  UserCtx	*user;
  PetscInt	i,bi,ibi, level, Add_dUndt=1,MHV_stuck=0 ;

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;

/* ==================================================================================             */
/*     Store old values to determine SC convergence */
  if (movefsi || rotatefsi || MHV || fish || cop) {
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
			for (i=0;i<6;i++){
				fsi[ibi].S_old[i] = fsi[ibi].S_new[i];
				fsi[ibi].S_ang_o[i]=fsi[ibi].S_ang_n[i];      
				if (itr_sc==1) {
					fsi[ibi].dS[i]=0.;
					fsi[ibi].atk=0.3;	
				}
				fsi[ibi].dS_o[i]=fsi[ibi].dS[i];
				fsi[ibi].atk_o=fsi[ibi].atk;
			}
			if (itr_sc==2) fsi[ibi].atk_o=0.298;
			fsi[ibi].F_x_old=fsi[ibi].F_x;
			fsi[ibi].F_y_old=fsi[ibi].F_y;
			fsi[ibi].F_z_old=fsi[ibi].F_z;
			fsi[ibi].M_x_old=fsi[ibi].M_x;
			fsi[ibi].M_y_old=fsi[ibi].M_y;
			fsi[ibi].M_z_old=fsi[ibi].M_z;
    }
  }
/* ==================================================================================             */
/*     Calculating Forces! */
  if (MHV) Add_dUndt=0;
  if (immersed){
    for (bi=0; bi<block_number; bi++){
      for (ibi=0;ibi<NumberOfBodies;ibi++){
				if (levelset) Calc_forces_SI_levelset(&fsi[ibi],&(user[bi]),&ibm[ibi], ti, ibi, bi);
				else Calc_forces_SI(&fsi[ibi],&(user[bi]),&ibm[ibi], ti, ibi, bi);
      }
      if (itr_sc==1) VecCopy(user[bi].Ucat, user[bi].Ucat_o);
      if ((MHV || movefsi || rotatefsi || fish) && itr_sc>1){
				PetscPrintf(PETSC_COMM_WORLD, "Corrector Step itr # %d\n", itr_sc);
				VecCopy(user[bi].Ucont_o, user[bi].Ucont);
				VecCopy(user[bi].P_o, user[bi].P);
				DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
				DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
				DAGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
				DAGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
				Contra2Cart(&(user[bi]));
      }
      PetscBarrier(PETSC_NULL);	  
    }
  }
/* ==================================================================================             */
/*     Find The new Position & Move the BODY */
  if (movefsi && !fsi_6dof){
		for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--){
			user = usermg->mgctx[level].user;
			for (bi=0; bi<block_number; bi++){
				if (immersed){
					for (ibi=0;ibi<NumberOfBodies;ibi++){
						if(forced_motion){
							Forced_Motion(&fsi[ibi], 0.5,user->dt);	  
						}
						else{
							if(levelset){
								Calc_FSI_pos_SC_levelset(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt, user->ren);
							}
							else {
								Calc_FSI_pos_SC(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt, user->ren);			
							}
						}
					}
					CollisionDetectionOfCylinders(fsi,NumberOfBodies);
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi]);
					}
					PetscBarrier(PETSC_NULL);
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
						ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
						PetscBarrier(PETSC_NULL);	  
					}
				}
			}
		}
  }

  if (movefsi && fsi_6dof){
		for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--){
			user = usermg->mgctx[level].user;
			for (bi=0; bi<block_number; bi++){
				if (immersed){
					for (ibi=0;ibi<NumberOfBodies;ibi++){
						if(forced_motion){
							Forced_Motion(&fsi[ibi], 0.5,user->dt);	  
						}
						else{
							if(levelset){
								Calc_FSI_pos_6dof_levelset(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt, user->ren);
							}
							else {
								Calc_FSI_pos_SC(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt, user->ren);			
							}
						}
					}
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						Elmt_Move_FSI_ROT_TRANS(&fsi[ibi], &ibm[ibi],user->dt,0);					
					}
					PetscBarrier(PETSC_NULL);
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
						ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
						PetscBarrier(PETSC_NULL);	  
					}
				}
			}
		}
  }
  
  if (rotatefsi) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
			user = usermg->mgctx[level].user;
			for (bi=0; bi<block_number; bi++) {
				if (immersed) {
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						//Calc_FSI_Ang(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt,ibi,&user[bi]) ;
						//Forced_Rotation(fsi, 20*pi0180,2.*pi*0.185*0.97,user->dt);
						fsi[ibi].x_c = x_r;
						fsi[ibi].y_c = y_r;
						fsi[ibi].z_c = z_r;
						if(ibi==0) {
							Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
						}
						PetscBarrier(PETSC_NULL);
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA solver.c begin\n");
						ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA solver.c end\n");
						PetscBarrier(PETSC_NULL);	  
					}
				}
			}
		}
  }

  if (MHV && ti>-1) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      PetscReal dir=1.;
      for (bi=0; bi<block_number; bi++) {
				if (immersed) {
					// for leaflets ibi = 1 & 2
					for (ibi=1;ibi<NumberOfBodies;ibi++) {
						dir = -1*dir;
						Calc_FSI_Ang_intg(&fsi[ibi], &ibm[ibi], user->dt,itr_sc,ibi,&user[bi]) ;
						if ((dir*fsi[ibi].S_ang_n[0])> -max_angle) {
							fsi[ibi].S_ang_n[0]= -dir*max_angle;
							fsi[ibi].S_ang_n[1]= 0.;
							MHV_stuck=1;
						}
						if ((dir*fsi[ibi].S_ang_n[0])< 0.0) {
							fsi[ibi].S_ang_n[0]= dir*0.0;
							fsi[ibi].S_ang_n[1]= 0.;
							MHV_stuck=1;
						}
					}
					for (ibi=1;ibi<NumberOfBodies;ibi++) {
						Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
					}
					PetscBarrier(PETSC_NULL);
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
						ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
						PetscBarrier(PETSC_NULL);	  
					}
				}
      }
    }
  }





/* ==================================================================================             */
/*   Convergence of the SC Loop */
/* ==================================================================================             */

  *DoSCLoop = PETSC_FALSE;
  dSmax=1e-10;
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
		if ((movefsi || fish) && STRONG_COUPLING) {
			for (i=0;i<6;i++) {
				dS_sc = fabs(fsi[ibi].S_new[i]-fsi[ibi].S_old[i]);
				if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
				if (dS_sc > dSmax) dSmax=dS_sc;
			}
		}

		if ((rotatefsi||MHV) && STRONG_COUPLING) {
			dS_sc = fabs(fsi[ibi].S_ang_n[0]-fsi[ibi].S_ang_o[0]);
			if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
			dS_sc = fabs(fsi[ibi].S_ang_n[1]-fsi[ibi].S_ang_o[1]);
			if(fabs(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_o[1])>2.) 
			dS_sc /= 0.5*fabs(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_o[1]);
			if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
			if (dS_sc > dSmax) dSmax=dS_sc;
		}
  
  if ((movefsi || rotatefsi || MHV || fish) && STRONG_COUPLING && itr_sc<2) 
    *DoSCLoop = PETSC_TRUE;  
  } // ibi  

  if (itr_sc>10) *DoSCLoop = PETSC_FALSE;
  if (MHV_stuck && itr_sc==1) *DoSCLoop = PETSC_FALSE;

  PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence %d %le %le %le\n", itr_sc, dSmax,fsi[1].S_ang_n[1],fsi[1].S_ang_o[1]);
  
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    if ((movefsi || rotatefsi || cop || MHV) && !(*DoSCLoop)) FSI_DATA_Output(&fsi[ibi], ibi);
  }
/* ==================================================================================             */
  return(0);
}
/* ==================================================================================             */

PetscErrorCode Struc_predictor(UserMG *usermg,IBMNodes *ibm, 
			       FSInfo *fsi, PetscInt itr_sc,
			       PetscInt tistart, 
			       PetscTruth *DoSCLoop)
{
  UserCtx	*user;
  PetscInt	bi,ibi, level, MHV_stuck=0 ;
  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;
  if (MHV && ti>-1) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      PetscReal dir=1.;

      for (bi=0; bi<block_number; bi++) {
				if (immersed) {
					// for leaflets ibi = 1 & 2
					for (ibi=1;ibi<NumberOfBodies;ibi++) {
						dir = -1*dir;
						fsi[ibi].S_ang_n[1]=2.*fsi[ibi].S_ang_r[1]-fsi[ibi].S_ang_rm1[1];
						fsi[ibi].S_ang_n[0]=fsi[ibi].S_ang_r[0]+0.5*(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_r[1])*user->dt;

						if ((dir*fsi[ibi].S_ang_n[0])> -max_angle) {
							fsi[ibi].S_ang_n[0]= -dir*max_angle;
							fsi[ibi].S_ang_n[1]= 0.;
							MHV_stuck=1;
						}
						if ((dir*fsi[ibi].S_ang_n[0])< 0.0) {
							fsi[ibi].S_ang_n[0]= dir*0.0;
							fsi[ibi].S_ang_n[1]= 0.;
							MHV_stuck=1;
						}
					}
					for (ibi=1;ibi<NumberOfBodies;ibi++) {
						Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
					}
					PetscBarrier(PETSC_NULL);
					for (ibi=0;ibi<NumberOfBodies;ibi++) {
						PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
						ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
						PetscBarrier(PETSC_NULL);	  
					}
				}
      }
    }
  }    
  return(0);
}
/* ==================================================================================             */


/* ==================================================================================             */
/*     Flow Solver! */
/* ==================================================================================             */
PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi, PetscInt itr_sc, IBMNodes *wtm, ACL *acl, FSInfo *fsi_wt, IBMNodes *ibm_ACD, FSInfo *fsi_IBDelta,IBMNodes *ibm_IBDelta, IBMNodes *ibm_acl2ref, FSInfo *fsi_acl2ref, IBMNodes *ibm_nacelle, FSInfo *fsi_nacelle)
{
  UserCtx	*user;
  PetscInt	bi, level;
	PetscReal ts, te, cputime;     // xiaolei
	int ibi;
	/* ==================================================================================             */
  if (immersed) {
    for (level=usermg->mglevels-1; level>0; level--) {
      for (bi=0; bi<block_number; bi++) {
				MyNvertRestriction(&(usermg->mgctx[level].user[bi]), &(usermg->mgctx[level-1].user[bi]));
      }
    }
  }
	/* ==================================================================================             */

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;
	
	Calc_Minimum_dt(user);	// momentum.c
  
	#ifdef DIRICHLET
	if(freesurface && ti!=tistart) {
		void update_free_surface_position(UserCtx *user);
		update_free_surface_position(&user[0]);
	}
	#endif

	/* ==================================================================================             */
	/*   Momentum Solver! */
	/* ==================================================================================             */
	if(ti==tistart || movefsi || rotatefsi) {
		for (bi=0; bi<block_number; bi++) {
			if (immersed) {
				ibm_interpolation_advanced(&user[bi]);
			}
			IB_BC(&user[bi]);
			DALocalToGlobal(user[bi].fda, user[bi].lUcont, INSERT_VALUES, user[bi].Ucont);
		}
	}

	if(pseudo_periodic || k_periodic || kk_periodic) {
		pseudo_periodic_BC(user);
		if(save_inflow) save_inflow_section(user);
	}

	if (inletprofile==100) {	// read inflow data
		read_inflow_section(user);
	}

	Calc_k_Flux(&user[0]);

	if(les){
		DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
		DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
		Contra2Cart(user);
		if(ti%dynamic_freq==0 || ti==tistart) Compute_Smagorinsky_Constant_1(user, user->lUcont, user->lUcat);
		Compute_eddy_viscosity_LES(user);
	}

	if (levelset && itr_sc==1) {
		if(ti==tistart && ti==0) {
			Levelset_Function_IC(&user[0]);
		}
		if(!fix_level && ti!=0) {
			double dt = 1.*user[0].dt;
			VecCopy(user[0].Levelset, user[0].Levelset_o);
			DAGlobalToLocalBegin(user[0].da, user[0].Levelset, INSERT_VALUES, user[0].lLevelset);
			DAGlobalToLocalEnd(user[0].da, user[0].Levelset, INSERT_VALUES, user[0].lLevelset);
			Levelset_BC(&user[0]);
			Advect_Levelset(&user[0],dt);
			Levelset_BC(&user[0]);
			Reinit_Levelset(&user[0]);
			Levelset_BC(&user[0]);
		}
		Compute_Density(&user[0]);
		if(surface_tension) Compute_Surface_Tension(&user[0]);
	}

	if(rans) {
	  extern char path[256];
	  char filen[256];
	  PetscViewer     viewer;
	  if( ti==tistart && rans==3 ) {
	    bi=0;
	    sprintf(filen, "%s/Distance_%1d.dat", path, user->_this);
	    if(!ti || !file_exist(filen)) {
	      Compute_Distance_Function(user);
	      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	      VecView(user->Distance, viewer);
	      PetscViewerDestroy(viewer);
	    }
	  }
		if(ti==0) {
			PetscPrintf(PETSC_COMM_WORLD, "\nInitializing K-omega ... \n\n");
			K_Omega_IC(user);
			VecSet(user->lNu_t, user->ren);
			bi=0;
			VecCopy(user[bi].K_Omega, user[bi].K_Omega_o);
			DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
			DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
			DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
			DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
		}
		else {
		  bi=0;
		  if(ti==tistart) {
		    VecCopy(user[bi].K_Omega, user[bi].K_Omega_o);
		    DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
		    DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
		    DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
		    DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
		  }
		  K_Omega_Set_Constant(user);
		}
	}

	VecDuplicate(user[0].lUcont, &user[0].Fp);
	VecDuplicate(user[0].lUcont, &user[0].Div1);
	VecDuplicate(user[0].lUcont, &user[0].Div2);
	VecDuplicate(user[0].lUcont, &user[0].Div3);
	VecDuplicate(user[0].lUcont, &user[0].Visc1);
	VecDuplicate(user[0].lUcont, &user[0].Visc2);
	VecDuplicate(user[0].lUcont, &user[0].Visc3);

	Pressure_Gradient(&user[0], user[0].dP);

	//add (Toni)	
	//for wave_momentum_source
	if(wave_momentum_source && itr_sc==1){
		if(wave_momentum_source==1){
			extern void WAVE_DATA_input(UserCtx *user);
			if (ti == (ti/wave_skip)*wave_skip || ti==tistart) {
				for (bi=0; bi<block_number; bi++){
					WAVE_DATA_input(&user[bi]);	
				}
			}		
		}
		else if(wave_momentum_source==2) {
			extern void WAVE_DATA_input(UserCtx *user);
			for (bi=0; bi<block_number; bi++){
				WAVE_DATA_input(&user[bi]);	
			}	
		}
		else if(wave_momentum_source==3 || wave_momentum_source==4) {
			extern void WAVE_DATA_input(UserCtx *user);
			for (bi=0; bi<block_number; bi++){
				WAVE_DATA_input(&user[bi]);	
			}	
		}	
		extern void WAVE_Formfuction2(UserCtx *user);
		for (bi=0; bi<block_number; bi++){
			WAVE_Formfuction2(&user[bi]);
		}	
	}	
	if(wave_sponge_layer && itr_sc==1){
		extern void WAVE_SL_Formfuction2(UserCtx *user);
		for (bi=0; bi<block_number; bi++){
			WAVE_SL_Formfuction2(&user[bi]);
		}
	}			
	//for air_flow_levelset
	if(air_flow_levelset==2 && itr_sc==1){
		extern void WIND_DATA_input(UserCtx *user);
		if (ti == (ti/wind_skip)*wind_skip){
			for (bi=0; bi<block_number; bi++){
				WIND_DATA_input(&user[bi]);	
			}
		}
	}
	//end (Toni)
	if (nacelle_model || IB_delta || rotor_model) VecSet(user[0].lF_eul,0.0);
	if (nacelle_model || IB_delta || rotor_model) VecSet(user[0].F_eul,0.0);
	if (rotor_model && temperature_rotormodel && temperature) VecSet(user[0].lFtmprt_eul,0.0);
	if (rotor_model && temperature_rotormodel && temperature) VecSet(user[0].Ftmprt_eul,0.0);
	if (rotor_model) {
		int ibi;
		if (rotor_model == 1) {
			Calc_F_lagr(user, wtm, fsi_wt, NumberOfTurbines);
			char fname[80];
			sprintf(fname,"Turbine_AD01");
			Calc_forces_rotor(user, wtm, fsi_wt, 0, fname, NumberOfTurbines);
			int df = deltafunc;
			Calc_F_eul(user, wtm, fsi_wt, NumberOfTurbines, 1.0, df);
		}
		if (rotor_model == 3) {
			PetscPrintf(PETSC_COMM_WORLD, "Calculate U_ref  \n");
			Uref_ACL(user, wtm, ibm_ACD, fsi_wt, NumberOfTurbines);
			Calc_turbineangvel(user->dt, wtm, fsi_wt);
			char fname[80];
			sprintf(fname,"line");
			for(ibi=0;ibi<NumberOfTurbines;ibi++) {
				int rotor_6dof=1;
				if(rotor_6dof) {
					//need to set turbine motion
					//rotor_Rot_6dof_fsi(&fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
					rotor_Rot_6dof_fsi(&fsi[ibi], &fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
				}
				else {
					rotor_Rot(&fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
				}
			}
			Pre_process(user, wtm, NumberOfTurbines);
			Calc_U_lagr(user, wtm, fsi_wt, NumberOfTurbines);
			Calc_F_lagr_ACL(user, wtm, fsi_wt, acl, NumberOfTurbines);
			Calc_forces_ACL(user, wtm, fsi_wt, 0);
			if(floating_turbine_case){
				fsi[0].Force_rotor_z=fsi_wt[0].Force_rotor_z;
				fsi[0].Mom_rotor_x=fsi_wt[0].Mom_rotor_x;
			}
			int df = deltafunc;
			Calc_F_eul(user, wtm, fsi_wt, NumberOfTurbines, 1.0, df);
		}
		if (rotor_model == 4) {
			PetscPrintf(PETSC_COMM_WORLD, "Calculate U_ref  \n");
			Uref_ACL(user, wtm, ibm_ACD, fsi_wt, NumberOfTurbines);
			Calc_F_lagr(user, wtm, fsi_wt, NumberOfTurbines);
			char fname[80];
			sprintf(fname,"Turbine_AD04");
			Calc_forces_rotor(user, wtm, fsi_wt, 0, fname, NumberOfTurbines);
			int df = deltafunc;
			Calc_F_eul(user, wtm, fsi_wt, NumberOfTurbines, 1.0, df);
		}
		if (rotor_model == 6) {
			PetscPrintf(PETSC_COMM_WORLD, "Calculate U_ref  \n");
			Uref_ACL(user, wtm, ibm_ACD, fsi_wt, NumberOfTurbines);
			Calc_turbineangvel(user->dt, wtm, fsi_wt);
			char fname[80];
			sprintf(fname,"line");
			for(ibi=0;ibi<NumberOfTurbines;ibi++) rotor_Rot(&fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
			Pre_process(user, wtm, NumberOfTurbines);
			for(ibi=0;ibi<NumberOfTurbines;ibi++) refAL_Rot(&fsi_wt[ibi], &wtm[ibi], &ibm_acl2ref[ibi], ibi);
			Pre_process(user, ibm_acl2ref, NumberOfTurbines);
			PetscPrintf(PETSC_COMM_WORLD, "AL: U_lagr\n");
			Calc_U_lagr(user, ibm_acl2ref, fsi_acl2ref, NumberOfTurbines);
			int i;
			for(ibi=0;ibi<NumberOfTurbines;ibi++) {
				for (i=0; i<wtm[ibi].n_elmt; i++) {
					wtm[ibi].U_lagr_x[i]= ibm_acl2ref[ibi].U_lagr_x[i];
					wtm[ibi].U_lagr_y[i]= ibm_acl2ref[ibi].U_lagr_y[i];
					wtm[ibi].U_lagr_z[i]= ibm_acl2ref[ibi].U_lagr_z[i];
				}
			}
			PetscPrintf(PETSC_COMM_WORLD, "lagrange force  \n");
			Calc_F_lagr_ACL(user, wtm, fsi_wt, acl, NumberOfTurbines);
			PetscPrintf(PETSC_COMM_WORLD, "write turbin_AL data \n");
			Calc_forces_ACL(user, wtm, fsi_wt, 0);
			int df = deltafunc;
			PetscPrintf(PETSC_COMM_WORLD, "distribute to Eulerian grid \n");
			Calc_F_eul(user, wtm, fsi_wt, NumberOfTurbines, 1.0, df);
		}
		PetscBarrier(PETSC_NULL);
	}
	if (IB_delta) {
		PetscReal ts1, te1;
		int my_rank;
		MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
		PetscGetTime(&ts);  // xiaolei
		// Be care, the geometry is assumed no changing during rotation
		if (rotate_IBdelta) {
			Cmpnts nr, na, nt; 
			int i, ibi;
			for(ibi=0;ibi<NumberOfIBDelta;ibi++) {
				for (i=0; i<ibm_IBDelta[ibi].n_v; i++) {
					na.x=fsi_IBDelta[ibi].nx_tb;	
					na.y=fsi_IBDelta[ibi].ny_tb;	
					na.z=fsi_IBDelta[ibi].nz_tb;	
					double rx = ibm_IBDelta[ibi].x_bp[i]-fsi_IBDelta[ibi].x_c;
					double ry = ibm_IBDelta[ibi].y_bp[i]-fsi_IBDelta[ibi].y_c;
					double rz = ibm_IBDelta[ibi].z_bp[i]-fsi_IBDelta[ibi].z_c;
					double rr = sqrt(rx*rx+ry*ry+rz*rz)+1.e-19;
					nr.x = rx/rr; 
					nr.y = ry/rr; 
					nr.z = rz/rr;
					nt.x=na.y*nr.z-na.z*nr.y;
					nt.y=na.z*nr.x-na.x*nr.z;
					nt.z=na.x*nr.y-na.y*nr.x;
					double Ut=fsi_IBDelta[ibi].angvel_axis*rr;
					ibm_IBDelta[ibi].u[i].x = Ut*nt.x;
					ibm_IBDelta[ibi].u[i].y = Ut*nt.y;
					ibm_IBDelta[ibi].u[i].z = Ut*nt.z;
				}
			}
		}
		Calc_F_lagr_noslip(user, ibm_IBDelta, fsi_IBDelta, NumberOfIBDelta);
		char fname[80];
		sprintf(fname,"IBDeltaForce");
		Calc_forces_rotor(user, ibm_IBDelta, fsi_IBDelta, 0, fname, NumberOfIBDelta);
		double dh=1.0;
		int df = deltafunc;
		Calc_F_eul(user, ibm_IBDelta, fsi_IBDelta, NumberOfIBDelta, dh, df); 
  	PetscBarrier(PETSC_NULL);
		PetscGetTime(&te);  // xiaolei
		PetscPrintf(PETSC_COMM_WORLD, "Time for IB_delta  %le\n", te-ts);
 }
	if(rotor_model){
		if (rotor_model == 3) {
			char fname[80];
			sprintf(fname,"line");
			for(ibi=0;ibi<NumberOfTurbines;ibi++) Export_lagrdata(&fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
		}
		else if (rotor_model == 6) {
			char fname[80];
			sprintf(fname,"line");
			for(ibi=0;ibi<NumberOfTurbines;ibi++) Export_lagrdata(&fsi_wt[ibi], &wtm[ibi], user->dt, ibi, fname, 1);
			sprintf(fname,"refline");
			for(ibi=0;ibi<NumberOfTurbines;ibi++) Export_lagrdata(&fsi_acl2ref[ibi], &ibm_acl2ref[ibi], user->dt, ibi, fname, 1);
		}
		PetscBarrier(PETSC_NULL);
	}
	if(inletprofile==20){}
	else if (implicit==1) ImplicitMomentumSolver(user, ibm, fsi);
	else if (implicit==2) ImplicitMomentumSolver1(user, ibm, fsi);
	else if (implicit==3) ImpRK(user, ibm, fsi);
	else if (implicit==4) {
	  VecSet (user[0].RHS_o, 0.);
	  Formfunction_2 (&user[0], user[0].RHS_o, 1.0);
		extern PetscErrorCode Implicit_MatrixFree(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);	  
		Implicit_MatrixFree(user, ibm, fsi);
	}
  else {
		COEF_TIME_ACCURACY=1.0;
		RungeKutta(user, ibm, fsi);
	}
	VecDestroy(user[0].Fp);
	VecDestroy(user[0].Div1);
	VecDestroy(user[0].Div2);
	VecDestroy(user[0].Div3);
	VecDestroy(user[0].Visc1);
	VecDestroy(user[0].Visc2);
	VecDestroy(user[0].Visc3);

	/* ==================================================================================             */
	/*    Poisson Solver! */
	/* ==================================================================================             */
	PetscBarrier(PETSC_NULL);

	for(bi=0; bi<block_number; bi++) {
		if(inletprofile==20){}
		else if(poisson==-1) PoissonSolver_MG_original(usermg, ibm, user[bi].ibm_intp);
		else if(poisson==0) PoissonSolver_MG(usermg, ibm, user[bi].ibm_intp);
		else if(poisson==1) PoissonSolver_Hypre(&user[bi], ibm, user[bi].ibm_intp);
	}

	/* ==================================================================================             */
	/*    Velocity Correction! */
	/* ==================================================================================             */
	if(inletprofile!=20)
	for (bi=0; bi<block_number; bi++) {
		UpdatePressure(&user[bi]);
		Projection(&(user[bi]));
		DAGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
		DAGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
	}

	/* ==================================================================================             */
	/*   BC!!    */
	/* ==================================================================================             */
	if (block_number>1) {
		Block_Interface_U(user);
	}
  for (bi=0; bi<block_number; bi++) {
    if (immersed) {
			ibm_interpolation_advanced(&user[bi]);
    }
	}
	
	/* ==================================================================================             */
  /*    Checking Convergence! */
	/* ==================================================================================             */
	bi = 0;
	Divergence(&(user[bi]));
	for (bi=0; bi<block_number; bi++) {
		IB_BC(&user[bi]);
		DALocalToGlobal(user[bi].fda, user[bi].lUcont, INSERT_VALUES, user[bi].Ucont);
		Contra2Cart(&(user[bi]));
	}
	bi=0;
	Calc_ShearStress(&user[0]);
 
 	/* ==================================================================================             */
  /*    Second step of the leveset method if necessary */
	/* ==================================================================================             */
	if (levelset && itr_sc==1 && !fix_level && levelset_solve2) {
		VecAXPBY(user[0].Levelset,.5,.5,user[0].Levelset_o);
		VecCopy(user[0].Levelset, user[0].Levelset_o);
		DAGlobalToLocalBegin(user[0].da, user[0].Levelset, INSERT_VALUES, user[0].lLevelset);
		DAGlobalToLocalEnd(user[0].da, user[0].Levelset, INSERT_VALUES, user[0].lLevelset);
		if(!fix_level && ti!=0) {
			double dt = .5*user[0].dt;
			Levelset_BC(&user[0]);
			Advect_Levelset(&user[0],dt);
			Levelset_BC(&user[0]);
			Reinit_Levelset(&user[0]);
			Levelset_BC(&user[0]);
		}		
	}
	
	if(levelset && export_FS_elev){
		if (ti == (ti/tiout)*tiout) {
			Calc_free_surface_location(&user[0]);
		}
	}
	
	if(levelset)write_data(&user[0]);

	/* ==================================================================================             */
	/*     OUTPUT Values! */
	/* ==================================================================================             */
	if(averaging) {
		extern PetscErrorCode Do_averaging(UserCtx *user);
		Do_averaging(&user[0]);
	}

	//add (Toni)
	//for wave_momentum_source in the case of reading a wave frequencies information file from far-field domain. It is the time elapsed since the external file was imported. When a new file is imported it is set back to zero in the function WAVE_DATA_input.
	if(wave_momentum_source==1){
		user[0].wave_inf[0].WAVE_time_since_read+=user->dt;
	}	
	//end (Toni)

	extern PetscErrorCode KE_Output(UserCtx *user);
	KE_Output(user);
	int i;
	for (i = 0; i < (rand() % 3000); i++) (rand() % 3000);	//seokkoo
	if(rans) {
		K_Omega_Set_Constant(user);
		Solve_K_Omega(user);
		VecCopy(user[bi].K_Omega, user[bi].K_Omega_o);
		DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
		DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
		DAGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
		DAGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
	}
	bi=0;
	if (ti == (ti/tiout) * tiout) {
		Ucont_P_Binary_Output(&(user[bi]));
	}
	else if (tiout_ufield>0 && ti == (ti/tiout_ufield) * tiout_ufield && ti<=tiend_ufield) {
		Ucat_Binary_Output(&(user[bi]));
	}
	/* ==================================================================================             */
	PetscPrintf(PETSC_COMM_WORLD, "Finish Flow Solver  \n");
	/*     End of Flow Solver! */
	/* ==================================================================================             */
  return(0);
}
