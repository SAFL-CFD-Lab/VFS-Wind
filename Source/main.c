/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"
static char help[] = "Testing programming!";
PetscReal COEF_TIME_ACCURACY=1.5;
PetscInt ti,tistart=0;
PetscReal	Flux_in = 4.104388e-04, angle = 0;
PetscInt tiout = 10;
PetscInt block_number;
PetscReal FluxInSum, FluxOutSum;
PetscReal FluxInSum_gas, FluxOutSum_gas;	// seokkoo
PetscInt immersed = 0;
PetscInt inviscid = 0;
PetscInt movefsi = 0, rotatefsi=0;
PetscInt implicit = 0;
PetscInt imp_MAX_IT = 50; 
PetscInt radi=10;
PetscInt inletprofile=1;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 0;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0;
PetscReal max_angle = -54.*3.1415926/180.;// for MHV; min=0
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=0,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=0 ;
PetscInt NumberOfBodies=1;
PetscReal L_dim;
PetscInt averaging=0;
PetscInt binary_input=0;
PetscInt xyz_input=0;
PetscInt les=0;
PetscInt inlet_buffer_k=1;
PetscInt wallfunction=0;
PetscInt slipbody=-1000;
PetscInt central=0, second_order=0;
//PetscInt initialzero=0;
PetscInt freesurface=0;
PetscInt rans=0, lowRe=0;
PetscInt cross_diffusion=1;
PetscInt surface_tension=0;
PetscInt sloshing=0;
PetscInt export_FS_elev_center=0;
PetscInt export_FS_elev=0;
PetscInt forced_motion=0;
PetscInt fall_cyll_case=0;
double sloshing_a=0.05, sloshing_b=2, sloshing_d=1;


// add (Toni)
	//Variables for levelset
double level_in_height=0.;	
int level_in=0;
int levelset_it=10;
double levelset_tau=0.01;
int levelset_weno=0;
	//Variables for wave_momentum_source
int wave_momentum_source=0;
int wave_sponge_layer=0;
double wave_sponge_zs=4.;//length of the sponge layer
double wave_sponge_z01=-10.;//start of the sponge layer at the left x boundary
double wave_sponge_z02=16.;//start of the sponge layer at the right x boundary	
double wave_sponge_xs=4.;//length of the sponge layer
double wave_sponge_x01=-10.;//start of the sponge layer at the left x boundary
double wave_sponge_x02=16.;//start of the sponge layer at the right x boundary	
double wave_angle_single=0.;
double wave_K_single=1.0;
double wave_depth=1.0;
double wave_a_single=0.0;	
int wave_ti_start=0;
	//End variables for wave_momentum_source
	//Variables for air_flow_levelset
int air_flow_levelset=0;
int air_flow_levelset_periodic=0;
int wave_average_k=0;
int wave_skip=0;
int wind_skip=0;
int wind_start_read=0;
int wave_start_read=0;
int wind_recicle=10000;
int wave_recicle=10000;
double wave_wind_reflength=1.0;
double wave_wind_refvel=1.0;
double wave_wind_yshift=0.0;
int floating_turbine_case=0;
int wave_k_ave=0, wave_i_ave=0;
int wave_ti_startave=1000;
int freesurface_wallmodel=0;
int viscosity_wallmodel=0;
double channel_height=1.0;
	//End variables for air_flow_levelset
	//variables for FSI
int fsi_6dof=0;
double body_mass=0.;
double body_inertia_x=0., body_inertia_y=0., body_inertia_z=0.;
double body_alpha_rot_x=0., body_alpha_rot_y=0., body_alpha_rot_z=0.;
double body_alpha_lin_x=0., body_alpha_lin_y=0., body_alpha_lin_z=0.;
double body_beta_rot_x=0., body_beta_rot_y=0., body_beta_rot_z=0.;
double body_beta_lin_x=0., body_beta_lin_y=0., body_beta_lin_z=0.;
double angle_x0=0.,angle_y0=0.,angle_z0=0.;
// End (Toni)


// add (xiaolei)
PetscInt forcewidthfixed = 0; // fix the width for force distribution
PetscReal dhi_fixed, dhj_fixed, dhk_fixed;

PetscInt temperature_rotormodel = 0; // scalar source/sink on the blades
PetscReal tmprt_initval=0.; // the initial value for temperature field 
PetscReal u_settling=0., v_settling=0., w_settling=0.; // settling velocities in x-, y- and z-directions  

PetscReal poisson_threshold = 0.1; // xiaolei test
PetscInt fractional_Max_IT = 1; 
PetscInt imin_wm=0, imax_wm=0, jmin_wm=0, jmax_wm=0, kmin_wm=0, kmax_wm=0;
PetscInt imin_wmtmprt=0, imax_wmtmprt=0, jmin_wmtmprt=0, jmax_wmtmprt=0, kmin_wmtmprt=0, kmax_wmtmprt=0;
PetscInt powerlawwallmodel=0;

PetscInt NumberOfTurbines=1; //xyang
PetscInt NumberOfIBDelta=1; //xyang
PetscInt NumberOfNacelle=1; //xyang
PetscInt IB_delta = 0; // xyang
PetscInt NumIBPerLoc = 1;
PetscInt NumNacellePerLoc = 1;

PetscInt IB_wm = 0;
PetscInt IB_wmtmprt = 0;
PetscInt i_periodicIB = 0, j_periodicIB = 0, k_periodicIB = 0;
PetscInt ApproxBC_wm = 0;
PetscInt Shear_wm = 1;
PetscInt Vel_wm = 0;
PetscInt Force_wm = 0;
PetscInt infRe = 0;
//int num_innergrid = 41; // the number of grid for wall model xyang:
int dymmatch_wm = 1; // xyang
PetscInt wallmodel_test = 0;
PetscInt dp_wm = 0;
PetscReal alfa_wm = 0.0;

double xmin,xmax,ymin,ymax,zmin,zmax; // the range of domain used for moving frame with wind turbines

PetscInt nacelle_model = 0;  
PetscInt rotate_IBdelta = 0;  
PetscInt rotate_nacelle = 0;  
PetscInt rotor_model = 0;  // xyang 12-7-2010 1: actuator disk model 2, 3: actuator line model
PetscInt forceavgperiod_AL = 1;  // force averaing period 
PetscReal indf_ax = 0.25; // xyang 12-16-2010
PetscReal percent_weno = 0.0; // xyang 12-16-2010
PetscInt surface_p_out = 0; //xyang
PetscInt num_blade = 3; // xyang 1-21-2011
PetscInt num_foiltype = 2; // The number of airfoil types used in the blade xyang 03-18-2011
PetscReal dh1_wm = 0.001; // the first off-wall grid spacing for wall model 4-11-2011
PetscReal dhratio_wm = 1.05; // the ratio of grid spacing in wall model 
PetscReal reflength_wt = 1.0;  
PetscReal reflength_nacelle = 1.0;  
PetscReal refvel_wt = 1.0;  
PetscReal refvel_cfd = 1.0;  
PetscReal reflength_IBDelta = 1.0; 
PetscInt temperature = 0; 
PetscInt deltafunc = 10;
PetscInt add_fluctuations = 0;
PetscInt add_fluctuations_tmprt = 0;
PetscReal halfwidth_dfunc = 4.0;
PetscReal tipspeedratio = 4.1;
PetscReal r_nacelle = 0.0;
PetscReal L_nacelle = 1.0;
PetscReal dh_nacelle = 1.0;
PetscReal loc_refvel = 1;
PetscReal X_control = 3;

PetscInt les_prt = 0;

PetscInt AL_Noslip=0;
PetscReal u_frame, v_frame, w_frame;
PetscInt MoveFrame = 0;

PetscInt ii_periodicWT=0, jj_periodicWT=0, kk_periodicWT=0; // periodic WT, a row/column of ghost wind turbines needs to be added
PetscReal Sx_WT=1.0, Sy_WT=1.0, Sz_WT=1.0;
PetscInt Nx_WT, Ny_WT, Nz_WT;

// Idealized water wave
PetscReal a_iww, lamda_iww, C_iww;

char path_inflow[256];

PetscReal prt_eps=1.e-20;

PetscInt New_wallmodel;

PetscInt tisteps;

// Sponge layer at outlet for levelset
//
PetscInt SpongeLayer=0;
PetscInt SpongeDistance=50;

PetscInt MoveCylinderTest=0;

PetscInt inflow_levelset=0;
PetscInt sublevel=0;
PetscReal smoothlevel=0;

// TUrbine model

PetscReal dfunc_wd=2.0;

PetscInt FixTipSpeedRatio=1;
PetscInt fixturbineangvel=0;
PetscInt fixturbinegeneratortorque=0;

PetscInt rstart_turbinerotation=0;
PetscInt turbinetorquecontrol=0;
PetscInt Shen_AL=0;
PetscReal c0_CL=1;
PetscReal c1_CH=2.2;
PetscReal c2_CH=1;
PetscReal c3_CH=4;
PetscInt correction3D_CH=0;
PetscInt correction3D_CL=0;
PetscInt Shen1_AL=0;
PetscReal correction_ALShen1=1;
PetscReal relax_AL=1;
PetscReal a_shen=0.125;
PetscReal b_shen=21;
PetscReal c_shen=0.1;
PetscInt Prandtl_AL=0;
PetscInt smoothforce_AL=0;

PetscReal refangle_AL=20;

PetscInt radialforce_AL=0;
PetscReal cf_radialforce_AL=0.1;
PetscReal count_AL;

PetscInt turbine_prescrived_motion=0, turbine_prescrived_motion_heave=0, turbine_prescrived_motion_pitch=0;
PetscInt turbine_6dof_fsi_motion=0;
// Mobile bed sediment

PetscReal dt_bed=1;  // for every dt_bed flow time steps, bed is computed once. 
PetscReal particle_dens=1992.0;  // for every dt_bed flow time steps, bed is computed once. 
PetscReal particle_D50=0.0018;  // for every dt_bed flow time steps, bed is computed once. 
PetscInt smooth_bed=200;  // smooth coefficient for mobile bed. 
PetscInt bed_avg=200;  // smooth coefficient for mobile bed. 
PetscInt particlevel_model=0;  // model for parivle velocity. 
PetscInt LS_bedConstr=1;  // using least squares for constructing xp, yp, zp 
PetscReal Angle_repose=50; 
PetscReal bed_porosity=0.5; 
PetscReal sb_bed=0.01; 
PetscInt number_smooth=2;
PetscReal deltab=0.01;
PetscReal rs_bed=1;
 
// FSI

PetscInt prescribed_rotation=1;


PetscReal scale_velocity=1;
// end add (xiaolei)

double dt_inflow;

int levelset=0;
int levelset_solve2=0;
int fix_level=0;
int laplacian=0;
int qcr=0;
int poisson=1;
int amg_agg=1;
double amg_thresh=0.75;
int periodic=0;
int i_periodic=0;
int ii_periodic=0;
int j_periodic=0;
int jj_periodic=0;
int k_periodic=0;
int kk_periodic=0;
int pseudo_periodic=0;
double inlet_flux=-1;
int delete_previous_file=0;
int mixed=0;
int clark=0;
int vorticity=0;
int initial_perturbation=0;
int skew=0;
int dynamic_freq=1;
int my_rank;
int ib_bctype[128];
char path[256], gridfile[256];
int i_proc=PETSC_DECIDE, j_proc=PETSC_DECIDE, k_proc=PETSC_DECIDE;
double imp_free_tol=1.e-4;
double poisson_tol=5.e-9;	// relative tolerance
double les_eps=1.e-7;
//double Fr=1.0;

double mean_pressure_gradient=0;	// relative tolerance
PetscReal max_cs=0.5;
PetscTruth dpdz_set=PETSC_FALSE;
int save_inflow=0;
int save_inflow_period=100;
int save_inflow_minus=0;
int ti_lastsave=0;
int localstep=1;
int inflow_recycle_perioid=20000;
int save_memory=1;
int ibm_search=0;
PetscTruth rough_set=PETSC_FALSE;
double roughness_size=0.0;
int save_point[3000][10];
int nsave_points=0;
int testfilter_ik=0;
int testfilter_1d=0;
int i_homo_filter=0;
int j_homo_filter=0;
int k_homo_filter=0;
int poisson_it=10;
int tiout_ufield = -1, tiend_ufield = 10000000;
//int display_implicit_count=0;
double dx_min, di_min, dj_min, dk_min;
double di_max, dj_max, dk_max;

//double rho_water=1000., rho_air=1.204;	// bubble
//double mu_water=1, mu_air=0.1;
double dthick=1.5;
PetscTruth dthick_set=PETSC_FALSE;
double rho_water=1., rho_air=0.001;	// sloshing
double mu_water=1.e-3, mu_air=1.78e-5;
double angvel=3.141592;
//double rho_water=10., rho_air=0.01;

double gravity_x=0, gravity_y=0, gravity_z=0;
double inlet_y=0, outlet_y=0;
double inlet_z=0, outlet_z=0;
int fix_outlet=0, fix_inlet=0;
int rotdir=2; // 0: rotate around the x-axis, 1:y-axis, 2:z-axis
double x_r=0, y_r=0, z_r=0; // center of rotation of rfsi

PetscTruth	rstart_flg=PETSC_FALSE;
PetscTruth	inlet_y_flag=PETSC_FALSE;
PetscTruth	inlet_z_flag=PETSC_FALSE;

IBMNodes	*ibm_ptr;
FSInfo        *fsi_ptr;
UserCtx	*user_ptr;

int file_exist(char *str)
{
  int r=0;

  if(!my_rank) {
	FILE *fp=fopen(str, "r");
	if(!fp) {
	  r=0;
	  printf("\nFILE !!! %s does not exist !!!\n", str); 
	}
	else {
		fclose(fp);
		r=1;
	}
  }
  MPI_Bcast(&r, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  return r;
};

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
	PetscViewer	viewer;
	char filen[90];
	
	PetscInt N;
	VecGetSize(user->Ucont, &N);
	PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  
	sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucont));
	PetscViewerDestroy(viewer);

	PetscBarrier(PETSC_NULL);

	PetscOptionsClearValue("-vecload_block_size");

	sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->P));
	PetscViewerDestroy(viewer);

	sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Nvert_o));
	PetscViewerDestroy(viewer);
  
	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);	// Seokkoo Kang
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoadIntoVector(viewer, (user->Ucat));
	PetscViewerDestroy(viewer);
  
	if(!immersed) {
          VecSet(user->Nvert, 0.);
          VecSet(user->Nvert_o, 0.);
	}

	DAGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
	DAGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
	
	DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

	VecCopy(user->Ucont, user->Ucont_o);
	DAGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
	DAGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
  
	DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
	DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
  
	DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);
        DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);

	DAGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);
	DAGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);
  
	if(averaging) {	// Seokkoo Kang
		sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
		FILE *fp=fopen(filen, "r");
	  
		VecSet(user->Ucat_sum, 0);
		VecSet(user->Ucat_cross_sum, 0);
		VecSet(user->Ucat_square_sum, 0);
		VecSet(user->P_sum, 0);
	  
		if(les) {
		  VecSet(user->Nut_sum, 0.);
		}

		if(rans) {
		  VecSet(user->K_sum, 0.);
		}

		if(averaging>=2) {
			VecSet(user->P_square_sum, 0);
			//VecSet(user->P_cross_sum, 0);
		}
		if(averaging>=3) {
			if(les) {
				VecSet(user->tauS_sum, 0);
			}
			VecSet(user->Udp_sum, 0);
			VecSet(user->dU2_sum, 0);
			VecSet(user->UUU_sum, 0);
			VecSet(user->Vort_sum, 0);
			VecSet(user->Vort_square_sum, 0);
		}
		
		if(fp==NULL) {
			PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and contiues the computation ... ***\n\n", filen);
		}
		else {
			fclose(fp);
			PetscBarrier(PETSC_NULL);
			sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, user->Ucat_sum);
			PetscViewerDestroy(viewer);
			
			sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, user->Ucat_cross_sum);
			PetscViewerDestroy(viewer);
		  
			sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, (user->Ucat_square_sum));
			PetscViewerDestroy(viewer);
			
			sprintf(filen, "%s/sp_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, user->P_sum);
			PetscViewerDestroy(viewer);
			
			if(les) {
			  sprintf(filen, "%s/snut_%06d_%1d.dat", path, ti, user->_this);
			  if( file_exist(filen) ) {
			    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			    VecLoadIntoVector(viewer, user->Nut_sum);
			    PetscViewerDestroy(viewer);
			  }
			}
			
			if(rans) {
			  sprintf(filen, "%s/sk_%06d_%1d.dat", path, ti, user->_this);
			  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			  VecLoadIntoVector(viewer, user->K_sum);
			  PetscViewerDestroy(viewer);
			}

			if(averaging>=2) {
				sprintf(filen, "%s/sp2_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->P_square_sum);
				PetscViewerDestroy(viewer);
				/*
				sprintf(filen, "%s/sp1_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->P_cross_sum);
				PetscViewerDestroy(viewer);
				*/
			}
			
			if(averaging>=3) {
				if(les) {
					sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
					if( file_exist(filen) ) {
					  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					  VecLoadIntoVector(viewer, user->tauS_sum);
					  PetscViewerDestroy(viewer);
					}
					else PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !\n", filen);
				}
				
				sprintf(filen, "%s/su3_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->Udp_sum);
				PetscViewerDestroy(viewer);

				sprintf(filen, "%s/su4_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->dU2_sum);
				PetscViewerDestroy(viewer);

				sprintf(filen, "%s/su5_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->UUU_sum);
				PetscViewerDestroy(viewer);

				sprintf(filen, "%s/svo_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->Vort_sum);
				PetscViewerDestroy(viewer);
				
				sprintf(filen, "%s/svo2_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoadIntoVector(viewer, user->Vort_square_sum);
				PetscViewerDestroy(viewer);
			}
			
			PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
		}
	}
  
	if(levelset) {
		
		sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti, user->_this);
		FILE *fp=fopen(filen, "r");

		if(fp==NULL) {
			PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, terminates ... ***\n\n", filen);
			PetscFinalize();
			exit(0);
		}
		else {
			fclose(fp);
		
			PetscBarrier(PETSC_NULL);
			
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, user->Levelset);
			PetscViewerDestroy(viewer);
			
			VecCopy (user->Levelset, user->Levelset_o);
			DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
			DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
		}
		
	}
	
	if(les) {
		Vec Cs;
		VecDuplicate(user->P, &Cs);
		
		sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti, user->_this);
		FILE *fp=fopen(filen, "r");

		if(fp==NULL) {
			PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
			VecSet(Cs, 0);
		}
		else {
			fclose(fp);
		
			PetscBarrier(PETSC_NULL);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer, Cs);
			PetscViewerDestroy(viewer);
		}
		
		DAGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
		DAGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
		
		VecDestroy(Cs);
	}
  
	if(rans) {
		// K-Omega
		sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
		if( file_exist(filen) ) {
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer,user->K_Omega);
			PetscViewerDestroy(viewer);
		}
		else {
		  //K_Omega_IC(user);
			PetscPrintf(PETSC_COMM_WORLD, "\nInitializing K-omega ... \n\n");
                        K_Omega_IC(user);
                        VecSet(user->lNu_t, user->ren);
		}
		
		VecCopy(user->K_Omega, user->K_Omega_o);
			
		DAGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
		DAGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
			
		DAGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
		DAGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
		
		if(rans==3) {
		  // distance
		  sprintf(filen, "%s/Distance_%1d.dat", path, user->_this);
		  if( file_exist(filen) ) {
		    PetscPrintf(PETSC_COMM_WORLD, "Reading %s !\n", filen);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoadIntoVector(viewer,user->Distance);
			PetscViewerDestroy(viewer);
		  }
		  
		  else {
		    PetscPrintf(PETSC_COMM_WORLD, "File %s does not exist. Recalculating distance function !\n", filen);
		    /*Compute_Distance_Function(user);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Distance, viewer);
			PetscViewerDestroy(viewer);
		    */
		  }
		  
		}
	}
	
	return 0;
}

PetscErrorCode Ucat_Binary_Output(UserCtx *user)
{
	PetscViewer	viewer;
	char filen[80];
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	PetscBarrier(PETSC_NULL);
	
	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucat, viewer);
	PetscViewerDestroy(viewer);
	sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(PETSC_NULL);
  
	return 0;
}

int delete_count=0;

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user)
{
	PetscViewer	viewer;
	char filen[80];
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	
	PetscBarrier(PETSC_NULL);
	
	sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucont, viewer);
	PetscViewerDestroy(viewer);
	sprintf(filen, "%s/vfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(PETSC_NULL);

	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucat, viewer);
	PetscViewerDestroy(viewer);
	sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(PETSC_NULL);

	sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->P, viewer);
	PetscViewerDestroy(viewer);
	sprintf(filen, "%s/pfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(PETSC_NULL);
	
	if(qcr) {
		Vec Q;
		VecDuplicate(user->P, &Q);
		Compute_Q(user,  Q);
		
		sprintf(filen, "%s/qfield%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(Q, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/qfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		VecDestroy(Q);
		PetscBarrier(PETSC_NULL);
	}
	
	sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Nvert, viewer);
	PetscViewerDestroy(viewer);
	sprintf(filen, "%s/nvfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(PETSC_NULL);
  
	if(averaging && ti!=0) {	// Seokkoo Kang
		sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_sum, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/su0_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(PETSC_NULL);
		  
		sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_cross_sum, viewer);
		PetscViewerDestroy(viewer);  
		sprintf(filen, "%s/su1_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(PETSC_NULL);
		
		sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_square_sum, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/su2_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(PETSC_NULL);
		  
		sprintf(filen, "%s/sp_%06d_%1d.dat",path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->P_sum, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/sp_%06d_%1d.dat.info",path, ti, user->_this);	if(!rank) unlink(filen);
		
		if(les) {
		  sprintf(filen, "%s/snut_%06d_%1d.dat",path, ti, user->_this);
		  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		  VecView(user->Nut_sum, viewer);
		  PetscViewerDestroy(viewer);
		  sprintf(filen, "%s/snut_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}

		if(rans) {
		  sprintf(filen, "%s/sk_%06d_%1d.dat",path, ti, user->_this);
		  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		  VecView(user->K_sum, viewer);
		  PetscViewerDestroy(viewer);
		  sprintf(filen, "%s/sk_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}









		if(averaging>=2) {
			sprintf(filen, "%s/sp2_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->P_square_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/sp2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			/*
			sprintf(filen, "%s/sp1_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->P_cross_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/sp1_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			*/
		}
		
		if(averaging>=3) {
			if(les) {
				sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
				VecView(user->tauS_sum, viewer);
				PetscViewerDestroy(viewer);
				sprintf(filen, "%s/stauS_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);
			}
				
			sprintf(filen, "%s/su3_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Udp_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/su3_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/su4_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->dU2_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/su4_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/su5_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->UUU_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/su5_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/svo_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Vort_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/svo_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			
			sprintf(filen, "%s/svo2_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Vort_square_sum, viewer);
			PetscViewerDestroy(viewer);
			sprintf(filen, "%s/svo2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}
		
		PetscBarrier(PETSC_NULL);
	}
  
	if(levelset) {
		sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti, user->_this);
		
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Levelset, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/lfield%06d_%1d.dat.info",path, ti, user->_this);	if(!rank) unlink(filen);
	}
	
	if(les) {
		Vec Cs;
		
		VecDuplicate(user->P, &Cs);
		DALocalToGlobal(user->da, user->lCs, INSERT_VALUES, Cs);
		
		sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(Cs, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/cs_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(PETSC_NULL);
		VecDestroy(Cs);
	}
	
	if(rans) {
		sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->K_Omega, viewer);
		PetscViewerDestroy(viewer);
		sprintf(filen, "%s/kfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(PETSC_NULL);
	}
  
	if(!rank && delete_previous_file && delete_count++>=2 && ti-tiout*2!=0) {
		sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		
		if(!(tiout_ufield>0 && ti == (ti/tiout_ufield) * tiout_ufield && ti<=tiend_ufield)) {
			sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		}
		sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		if(averaging) {
			sprintf(filen, "%s/sp_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			if(averaging>=2) {
			  //sprintf(filen, "%s/sp1_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
				sprintf(filen, "%s/sp2_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			}
			if(averaging>=3) {
			  sprintf(filen, "%s/su3_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);
			  sprintf(filen, "%s/su4_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);
			  sprintf(filen, "%s/su5_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);

				sprintf(filen, "%s/svo_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
				sprintf(filen, "%s/svo2_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
			}
			if(les) {
			  sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
			  sprintf(filen, "%s/snut_%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
			}
		}
		if(les>=2) {
			sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		}
		if(rans) {
			sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti-tiout*2, user->_this);	if(!rank) unlink(filen);
		}
		if(levelset) {
			sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
		}

	}
	PetscBarrier(PETSC_NULL);
  
	return 0;
}

PetscErrorCode Divergence(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DAVecGetArray(fda,user->lUcont, &ucont);
  DAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DAVecGetArray(da, Div, &div);
  DAVecGetArray(da, user->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x + ucont[k][j][i].y - ucont[k][j-1][i].y + ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	//if(i==mx-2) printf("%f %f %f %f %f %f\n", ucont[k][j][i].x, ucont[k][j][i-1].x, ucont[k][j][i].y, ucont[k][j-1][i].y, ucont[k][j][i].z, ucont[k-1][j][i].z);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] + nvert[k][j+1][i] + nvert[k][j-1][i] + nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	if (air_flow_levelset && k==lze-1) maxdiv = 0.;
	
	
	div[k][j][i] = maxdiv;
	
      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv %d %d %e\n", ti, i, maxdiv);
  PetscInt mi;
  

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (lidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k);
	}
      }
    }
  }
  
	
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
	FILE *f;
	char filen[80];
	sprintf(filen, "%s/Converge_dU", path);
	f = fopen(filen, "a");
	PetscFPrintf(PETSC_COMM_WORLD, f, " Maxdiv=%.2e\n", maxdiv);
	fclose(f);
  }
  
  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lUcont, &ucont);
  DAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(Div);
  return(0);
}

void write_data(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;

	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt	lxs, lys, lzs, lxe, lye, lze;
	PetscInt	i, j, k;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	double lvol=0, vol=0;

	PetscReal ***p, ***aj, ***nvert, ***rho, ***level;
	Cmpnts	***ucat, ***ucont, ***csi, ***eta, ***zet, ***cent;
  
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lLevelset, &level);
	}
	DAVecGetArray(da,user->lAj, &aj);
	DAVecGetArray(da, user->lP, &p);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda,user->lCsi, &csi);
	DAVecGetArray(fda,user->lEta, &eta);
	DAVecGetArray(fda,user->lZet, &zet);
	DAVecGetArray(fda,user->lUcat, &ucat);
	DAVecGetArray(fda,user->lUcont, &ucont);
	DAVecGetArray(fda,user->lCent, &cent);
	
	
	// for sloshing recording
	int ci=mx/2, ck=mz/2;
	double lz_sloshing=-10, z_sloshing;
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>0.1) continue;
		if(levelset && (sloshing || export_FS_elev_center)) {
			lvol += rho[k][j][i] / aj[k][j][i];
			if ( i==ci && k==ck) {
				if( level[k][j][i]>=0 && level[k][j+1][i]<0 ) {	// water surface is above my cell center
					if(sloshing || level_in==1){
						lz_sloshing = cent[k][j][i].z + level[k][j][i];					
					}
					else if(level_in==2) {
						lz_sloshing = cent[k][j][i].y + level[k][j][i];
					}
				}
				/*else if( level[k][j][i]<0 && level[k][j-1][i]>=0 ) {	// water surface is below my cell center
					lz_sloshing = cent[k][j][i].z + level[k][j][i];
				}*/
			}
		}
		
		for(int m=0; m<nsave_points; m++) {
			if(i==save_point[m][0] && j==save_point[m][1] && k==save_point[m][2]) {
				double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
				double dpdc, dpde, dpdz;
				double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
				double dp_dx, dp_dy, dp_dz;
			
				double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
				double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
				double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
				double ajc = aj[k][j][i];
				
				double Ai = sqrt ( csi0*csi0 + csi1*csi1 + csi2*csi2 );
				double Aj = sqrt ( eta0*eta0 + eta1*eta1 + eta2*eta2 );
				double Ak = sqrt ( zet0*zet0 + zet1*zet1 + zet2*zet2 );
				
				double U = 0.5*(ucont[k][j][i].x+ucont[k][j][i-1].x) / Ai;
				double V = 0.5*(ucont[k][j][i].y+ucont[k][j-1][i].y) / Aj;
				double W = 0.5*(ucont[k][j][i].z+ucont[k-1][j][i].z) / Ak;
			
				Compute_du_center ( i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
				
				Compute_dscalar_center ( i, j, k, mx, my, mz, p, nvert, &dpdc, &dpde, &dpdz );
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz );
				
				double vort_x = dw_dy - dv_dz,  vort_y = du_dz - dw_dx, vort_z = dv_dx - du_dy;
				
				FILE *f;
				char filen[80];
				
				sprintf(filen, "%s/Flow0_%04d_%04d_%04d_dt_%g.dat", path, i, j, k, user->dt);
				f = fopen(filen, "a");
				if(ti==tistart) fprintf(f, "\n");//fprintf(f, "\n\n**** Beginning to log time history of u v w p U V W dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz at point (%d,%d,%d) ****\n", i, j, k);
				fprintf(f, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n", ti, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z, p[k][j][i], U, V, W, vort_x, vort_y, vort_z);
				fclose(f);
				
				sprintf(filen, "%s/Flow1_%04d_%04d_%04d_dt_%g.dat", path, i, j, k, user->dt);
				f = fopen(filen, "a");
				if(ti==tistart) fprintf(f, "\n");
				fprintf(f, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n", ti, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dp_dx, dp_dy, dp_dz);
				fclose(f);
				
				break;
			}
		}
	}
	
	if(levelset && export_FS_elev_center) {
		PetscBarrier(PETSC_NULL);
		PetscGlobalSum(&lvol, &vol, PETSC_COMM_WORLD);
		PetscGlobalMax(&lz_sloshing, &z_sloshing, PETSC_COMM_WORLD);
		if(!my_rank) {
			char filen[256];
			sprintf(filen, "%s/mass.dat", path);
			FILE *fp = fopen(filen, "a");
			fprintf(fp, "%d %.10e\n", ti, vol );
			fclose(fp);
			double _time = (ti+1)*user->dt;
			sprintf(filen, "%s/fs_elev_center.dat", path);
			fp = fopen(filen, "a");
			fprintf(fp, "%e %.7e\n", _time,z_sloshing);
			fclose(fp);
		}
	} 
	if(levelset && sloshing) {
		PetscBarrier(PETSC_NULL);
		PetscGlobalSum(&lvol, &vol, PETSC_COMM_WORLD);
		PetscGlobalMax(&lz_sloshing, &z_sloshing, PETSC_COMM_WORLD);
		if(!my_rank) {
			char filen[256];
			sprintf(filen, "%s/mass.dat", path);
			FILE *fp = fopen(filen, "a");
			fprintf(fp, "%d %.10e\n", ti, vol );
			fclose(fp);
			double z_sloshing_exact;
			double _time = (ti+1)*user->dt;
			double a = 0.05, b = 2.0, d = 1.0, g = 1., x = 1.0, y;
			double k2 = 2*M_PI/b, k4 = 4*M_PI/b;
			double w2 = sqrt ( k2 * g * tanh (k2*d) ), w4 = sqrt ( k4 * g * tanh (k4*d) );
			
			if(inviscid && sloshing==1) {
				z_sloshing_exact = a * cos(w2*_time) * cos(k2*x);
				z_sloshing_exact += 0.125/g * ( 2*pow(w2*a,2.) * cos(2.*w2*_time) + pow(a/w2,2.) * ( pow(k2*g,2.) + pow(w2,4.) ) - pow(a/w2,2.) * ( pow(k2*g,2.) + 3.*pow(w2,4.) ) * cos(w4*_time) );
			}			
			else if(sloshing==2) {
			  i=ci;k=ck;j=5;
				x=10.;
				y=10.;
				double L = 20., d=1., Beta=0.25, g=fabs(gravity_z), a=0.1; // L = width of the 3D tank
				printf("x=%f, y=%f\n", x, y );
				double eta0 = a * exp ( -Beta * ( pow(x-L/2, 2) + pow(y-L/2, 2) ) );
				std::complex<double> I(0,1);
				std::complex<double> one(1,0);
				std::complex<double> val(0,0);
				for(int m=0; m<40; m++)
				for(int n=0; n<40; n++) {
					double coeff_m=2., coeff_n=2.;
					double k_mn = sqrt( pow(M_PI/L,2.) * (m*m + n*n) );
					double omega_mn = sqrt ( g * k_mn * tanh (k_mn * d) );
					if(m==0) coeff_m=1.;
					if(n==0) coeff_n=1.;
					std::complex<double> Integral;
					Integral   = M_PI * a / (16*Beta) * exp ( -0.25*M_PI * ( 2.*(m+n)*I + (m*m+n*n)*M_PI/(Beta*L*L)*one ) );
					Integral *= ( one + exp(m*M_PI*I) ) * ( one + exp(n*M_PI*I) );
					Integral *= ERF ( (Beta*L*L - m*M_PI*I) * pow(2*sqrt(Beta)*L, -1) ) + ERF ( (Beta*L*L + m*M_PI*I) * pow(2*sqrt(Beta)*L, -1) );
					Integral *= ERF ( (Beta*L*L - n*M_PI*I) * pow(2*sqrt(Beta)*L, -1) ) + ERF ( (Beta*L*L + n*M_PI*I) * pow(2*sqrt(Beta)*L, -1) );
		
					double eta_mn = (1./(L*L)) * coeff_m * coeff_n * Integral.real();
					val += eta_mn * exp ( -I * omega_mn * _time ) * cos (n * M_PI / L * x) * cos (m * M_PI / L * y);
				}
				z_sloshing_exact = val.real();
			}			
			
			
			sprintf(filen, "%s/sloshing.dat", path);
			fp = fopen(filen, "a");
			fprintf(fp, "%e %.7e %.7e %.7e\n", _time, -1.0+z_sloshing,  z_sloshing_exact, fabs(z_sloshing_exact-(-1.0+z_sloshing)));
			fclose(fp);
		}
	} 	

	PetscBarrier(PETSC_NULL);
	if(user->bctype[0]==11 && ti) {
		PetscGlobalSum(&user->lA_cyl, &user->A_cyl, PETSC_COMM_WORLD);
		PetscGlobalSum(&user->lA_cyl_x, &user->A_cyl_x, PETSC_COMM_WORLD);
		PetscGlobalSum(&user->lA_cyl_z, &user->A_cyl_z, PETSC_COMM_WORLD);
		
		PetscGlobalSum(&user->lFpx_cyl, &user->Fpx_cyl, PETSC_COMM_WORLD);
		PetscGlobalSum(&user->lFpz_cyl, &user->Fpz_cyl, PETSC_COMM_WORLD);
		PetscGlobalSum(&user->lFvx_cyl, &user->Fvx_cyl, PETSC_COMM_WORLD);
		PetscGlobalSum(&user->lFvz_cyl, &user->Fvz_cyl, PETSC_COMM_WORLD);
		
		//double Cpx=user->Fpx_cyl/user->A_cyl*2., Cpz=user->Fpz_cyl/user->A_cyl*2.;
		//double Cvx=user->Fvx_cyl/user->A_cyl*2., Cvz=user->Fvz_cyl/user->A_cyl*2.;
		
		double Cpx=user->Fpx_cyl/user->A_cyl_x*4., Cpz=user->Fpz_cyl/user->A_cyl_z*4.;
		double Cvx=user->Fvx_cyl/user->A_cyl_x*4., Cvz=user->Fvz_cyl/user->A_cyl_z*4.;
		
		double Cdx=Cpx+Cvx, Cdz=Cpz+Cvz;
		
		if(!my_rank) {
			char filen[256];
			sprintf(filen, "%s/Force_Cylinder.dat", path);
			FILE *fp = fopen(filen, "a");
			if(ti==tistart) fprintf(fp, "\n\n***\n\n");
			fprintf(fp, "%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", ti, Cdx, 0., Cdz, Cpx, 0., Cpz, user->A_cyl_x, 0., user->A_cyl_z, user->A_cyl );
			fclose(fp);
		}
	}
	PetscBarrier(PETSC_NULL);
	
	if(/*inletprofile==13 &&*/ ti) {
		
		if(!my_rank) {
			char filen[256];
			sprintf(filen, "%s/shear_velocity.dat", path);
			FILE *fp = fopen(filen, "a");
			//if(ti==tistart) fprintf(fp, "\n\n***\n\n");
			fprintf(fp, "%d ", ti);
			fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e\n", 
				user->ustar_now[0],user->ustar_now[1],user->ustar_now[2],user->ustar_now[3],user->ustar_now[4],user->ustar_now[5]);
			fclose(fp);
		}
	}
  
	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lLevelset, &level);
	}
	DAVecRestoreArray(da,user->lAj, &aj);
	DAVecRestoreArray(da, user->lP, &p);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(fda,user->lCsi, &csi);
	DAVecRestoreArray(fda,user->lEta, &eta);
	DAVecRestoreArray(fda,user->lZet, &zet);
	DAVecRestoreArray(fda,user->lUcat, &ucat);
	DAVecRestoreArray(fda,user->lUcont, &ucont);
	DAVecRestoreArray(fda,user->lCent, &cent);
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
	Vec	ResidualT;
	UserCtx	*user;
	PetscErrorCode	ierr;
	PetscInt	i, bi, ibi;
	PetscReal	norm;
	IBMNodes	*ibm, *ibm0, *ibm1;
	IBMInfo	*ibm_intp;

	// Added for fsi
	FSInfo        *fsi;
	PetscTruth    DoSCLoop;
	PetscInt      itr_sc;
	
	PetscInt level;
	UserMG usermg;
	PetscTruth	flg;
	//PetscInt tistart = 0;

	// begin add (xiaolei)
	ACL             *acl;   //xyang 3-18-2011
	IBMNodes        *wtm;   // xyang 04212011
	FSInfo        *fsi_wt;
	IBMNodes        *ibm_ACD;
	IBMNodes        *ibm_IBDelta;
	FSInfo        *fsi_IBDelta;
	IBMNodes        *ibm_acl2ref; // 20140807
	FSInfo        *fsi_acl2ref; // 20140807
	IBMNodes        *ibm_nacelle; 
	FSInfo        *fsi_nacelle; 
	// end add (xiaolei)
	
	PetscInitialize(&argc, &argv, (char *)0, help);
	PetscBarrier(PETSC_NULL);
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
	srand( time(NULL)) ;	// Seokkoo Kang

	PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
	
	PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-tiou", &tiout_ufield, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-tieu", &tiend_ufield, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-inv", &inviscid, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
	PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-reg", &regime, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-thin", &thin, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_z", &dgf_z, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_y", &dgf_y, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_x", &dgf_x, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_az", &dgf_az, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_ay", &dgf_ay, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-dgf_ax", &dgf_ax, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0
	PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, PETSC_NULL);	// Seokkoo Kang: if 1 binary PLOT3D file, if 0 ascii.
	PetscOptionsGetInt(PETSC_NULL, "-xyz", &xyz_input, PETSC_NULL);			// Seokkoo Kang: if 1 text xyz format, useful for very big (>1GB) Cartesian grid
	PetscOptionsGetInt(PETSC_NULL, "-les", &les, PETSC_NULL);				// Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
	PetscOptionsGetInt(PETSC_NULL, "-inlet_buffer_k", &inlet_buffer_k, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wallfunction", &wallfunction, PETSC_NULL);	// Seokkoo Kang: 1 or 2
	PetscOptionsGetInt(PETSC_NULL, "-slipbody", &slipbody, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-central", &central, PETSC_NULL);//central differencing
	PetscOptionsGetInt(PETSC_NULL, "-second_order", &second_order, PETSC_NULL);
	//PetscOptionsGetInt(PETSC_NULL, "-initialzero", &initialzero, PETSC_NULL);	// Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-freesurface", &freesurface, PETSC_NULL);	// Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);			// Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-cross_diffusion", &cross_diffusion, PETSC_NULL);                     // Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-lowRe", &lowRe, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-stension", &surface_tension, PETSC_NULL);			// Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-delete", &delete_previous_file, PETSC_NULL);	// Seokkoo Kang: delete previous time step's saved file for saving disk space
	PetscOptionsGetInt(PETSC_NULL, "-mixed", &mixed, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
	PetscOptionsGetInt(PETSC_NULL, "-clark", &clark, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
	PetscOptionsGetInt(PETSC_NULL, "-vorticity", &vorticity, PETSC_NULL);			// Seokkoo Kang: vorticity form for viscous terms
	PetscOptionsGetInt(PETSC_NULL, "-pseudo", &pseudo_periodic, PETSC_NULL);	// Seokkoo Kang: pseudo periodic BC in k-direction for genenration of inflow condition

	PetscOptionsGetInt(PETSC_NULL, "-levelset", &levelset, PETSC_NULL);     // Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-levelset_solve2", &levelset_solve2, PETSC_NULL);     // Seokkoo Kang	
	PetscOptionsGetInt(PETSC_NULL, "-fix_level", &fix_level, PETSC_NULL);     // Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-levelset_it", &levelset_it, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL, "-rotdir", &rotdir, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-laplacian", &laplacian, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-qcrout", &qcr, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ii_periodic", &ii_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-jj_periodic", &jj_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-kk_periodic", &kk_periodic, PETSC_NULL);	
	periodic = i_periodic+j_periodic+k_periodic+ii_periodic+jj_periodic+kk_periodic;
	PetscOptionsGetInt(PETSC_NULL, "-perturb", &initial_perturbation, PETSC_NULL);	// Seokkoo Kang: give a random perturbation for initial condition
	PetscOptionsGetInt(PETSC_NULL, "-skew", &skew, PETSC_NULL);				// Seokkoo Kang: skew symmetric form of advection term
	PetscOptionsGetInt(PETSC_NULL, "-dynamic_freq", &dynamic_freq, PETSC_NULL);		// Seokkoo Kang: LES dynamic compute frequency 
	if(dynamic_freq<1) dynamic_freq=1;

	PetscOptionsGetInt(PETSC_NULL, "-save_inflow", &save_inflow, PETSC_NULL);		// Seokkoo Kang: save infow BC to files; should be used in conjunction wiht -pseudo 1
	PetscOptionsGetInt(PETSC_NULL, "-save_inflow_period", &save_inflow_period, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-save_inflow_minus", &save_inflow_minus, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ti_lastsave", &ti_lastsave, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL, "-localstep", &localstep, PETSC_NULL);		// Seokkoo Kang: localstep ( explict + implicit momentum solver )
	PetscOptionsGetInt(PETSC_NULL, "-recycle", &inflow_recycle_perioid, PETSC_NULL);	// Seokkoo Kang, set recycling period of the inflow data
	PetscOptionsGetInt(PETSC_NULL, "-save_memory", &save_memory, PETSC_NULL);	// Seokkoo Kang, save_memory
	PetscOptionsGetInt(PETSC_NULL, "-ibm_search", &ibm_search, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ip", &i_proc, PETSC_NULL);			// Seokkoo Kang: number of processors in i direction
	PetscOptionsGetInt(PETSC_NULL, "-jp", &j_proc, PETSC_NULL);			// Seokkoo Kang: number of processors in j direction
	PetscOptionsGetInt(PETSC_NULL, "-kp", &k_proc, PETSC_NULL);			// Seokkoo Kang: number of processors in k direction
	PetscOptionsGetInt(PETSC_NULL, "-poisson", &poisson, PETSC_NULL); 	// Seokkoo Kang
	PetscOptionsGetInt(PETSC_NULL, "-amg_agg", &amg_agg, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-testfilter_ik", &testfilter_ik, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-testfilter_1d", &testfilter_1d, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-poisson_it", &poisson_it, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-i_homo_filter", &i_homo_filter, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-j_homo_filter", &j_homo_filter, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-k_homo_filter", &k_homo_filter, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-fix_outlet", &fix_outlet, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-fix_inlet", &fix_inlet, PETSC_NULL);
	
	// add (Toni)
	//for levelset
	PetscOptionsGetInt(PETSC_NULL, "-levelset_it", &levelset_it, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-levelset_tau", &levelset_tau, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-level_in_height", &level_in_height, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-level_in", &level_in, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-sloshing", &sloshing, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-levelset_weno", &levelset_weno, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-export_FS_elev_center", &export_FS_elev_center, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-export_FS_elev", &export_FS_elev, PETSC_NULL);	
	
	//for IB FSI
	PetscOptionsGetInt(PETSC_NULL, "-forced_motion", &forced_motion, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-fall_cyll_case", &fall_cyll_case, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-fsi_6dof", &fsi_6dof, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-body_mass", &body_mass, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_inertia_x", &body_inertia_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_inertia_y", &body_inertia_y, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_inertia_z", &body_inertia_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_rot_x", &body_alpha_rot_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_rot_y", &body_alpha_rot_y, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_rot_z", &body_alpha_rot_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_lin_x", &body_alpha_lin_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_lin_y", &body_alpha_lin_y, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_alpha_lin_z", &body_alpha_lin_z, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_rot_x", &body_beta_rot_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_rot_y", &body_beta_rot_y, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_rot_z", &body_beta_rot_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_lin_x", &body_beta_lin_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_lin_y", &body_beta_lin_y, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-body_beta_lin_z", &body_beta_lin_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-angle_x0", &(angle_x0), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-angle_y0", &(angle_y0), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-angle_z0", &(angle_z0), PETSC_NULL);

	//Variables for wave_momentum_source
	PetscOptionsGetInt(PETSC_NULL, "-wave_momentum_source", &wave_momentum_source, PETSC_NULL); //For momentum source 1 in x direction and 2 both directions
	PetscOptionsGetInt(PETSC_NULL, "-wave_sponge_layer", &wave_sponge_layer, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_zs", &wave_sponge_zs, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_z01", &wave_sponge_z01, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_z02", &wave_sponge_z02, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_xs", &wave_sponge_xs, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_x01", &wave_sponge_x01, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_x02", &wave_sponge_x02, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_angle_single", &wave_angle_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_K_single", &wave_K_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_depth", &wave_depth, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_a_single", &wave_a_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_reflength", &wave_wind_reflength, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_refvel", &wave_wind_refvel, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_yshift", &wave_wind_yshift, PETSC_NULL);	
	//for air_flow_levelset
	PetscOptionsGetInt(PETSC_NULL, "-air_flow_levelset", &air_flow_levelset, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-air_flow_levelset_periodic", &air_flow_levelset_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_average_k", &wave_average_k, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_skip", &wave_skip, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_ti_start", &wave_ti_start, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wind_skip", &wind_skip, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wind_start_read", &wind_start_read, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_start_read", &wave_start_read, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wind_recicle", &wind_recicle, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_recicle", &wave_recicle, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_k_ave", &wave_k_ave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_i_ave", &wave_i_ave, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_ti_startave", &wave_ti_startave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-freesurface_wallmodel", &freesurface_wallmodel, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-viscosity_wallmodel", &viscosity_wallmodel, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-channel_height", &channel_height, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-floating_turbine_case", &floating_turbine_case, PETSC_NULL);
	//End variables for wave_momentum_source
	// End (Toni)	
	
	// add begin (xiaolei)
	PetscOptionsGetInt(PETSC_NULL, "-forcewidthfixed", &forcewidthfixed, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetReal(PETSC_NULL, "-dhi_fixed", &dhi_fixed, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetReal(PETSC_NULL, "-dhj_fixed", &dhj_fixed, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetReal(PETSC_NULL, "-dhk_fixed", &dhk_fixed, PETSC_NULL); // xyang 6-3-2011
	if(forcewidthfixed) {
		PetscPrintf(PETSC_COMM_WORLD, "\n Actuator model: the width for force distribution is fixed at dhi=%le dhj=%le dhk=%le \n\n", dhi_fixed, dhj_fixed, dhj_fixed);
	}
	PetscOptionsGetInt(PETSC_NULL, "-temperature_rotormodel", &temperature_rotormodel, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetReal(PETSC_NULL, "-tmprt_initval", &tmprt_initval, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-u_settling", &u_settling, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-v_settling", &v_settling, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-w_settling", &w_settling, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-fractional_Max_IT", &fractional_Max_IT, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-wallmodel_test", &wallmodel_test, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-dp_wm", &dp_wm, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-alfa_wm", &alfa_wm, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-forceavgperiod_AL", &forceavgperiod_AL, PETSC_NULL); // xyang 12-7-2010
	PetscOptionsGetInt(PETSC_NULL, "-rotor_modeled", &rotor_model, PETSC_NULL); // xyang 12-7-2010
	PetscOptionsGetInt(PETSC_NULL, "-rotate_IBdelta", &rotate_IBdelta, PETSC_NULL); // xyang 12-7-2010
	PetscOptionsGetInt(PETSC_NULL, "-rotate_nacelle", &rotate_nacelle, PETSC_NULL); // xyang 12-7-2010
	PetscOptionsGetInt(PETSC_NULL, "-nacelle_model", &nacelle_model, PETSC_NULL); // xyang 12-7-2010
	PetscOptionsGetReal(PETSC_NULL, "-indf_ax", &indf_ax, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-percent_weno", &percent_weno, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetInt(PETSC_NULL, "-imin_wm", &imin_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-imax_wm", &imax_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-jmin_wm", &jmin_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-jmax_wm", &jmax_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-kmin_wm", &kmin_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-kmax_wm", &kmax_wm, PETSC_NULL); // xyang 1-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-imin_wmtmprt", &imin_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-imax_wmtmprt", &imax_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-jmin_wmtmprt", &jmin_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-jmax_wmtmprt", &jmax_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-kmin_wmtmprt", &kmin_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-kmax_wmtmprt", &kmax_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-surface_p_out", &surface_p_out, PETSC_NULL); // xyang 3-2-2011
	PetscOptionsGetInt(PETSC_NULL, "-num_blade", &num_blade, PETSC_NULL); // xyang 3-2-2011
	PetscOptionsGetInt(PETSC_NULL, "-num_foiltype", &num_foiltype, PETSC_NULL); // xyang 3-18-2011
	PetscOptionsGetReal(PETSC_NULL, "-dh1_wm", &dh1_wm, PETSC_NULL); //xyang 4-11-2011
	PetscOptionsGetReal(PETSC_NULL, "-dhratio_wm", &dhratio_wm, PETSC_NULL); //xyang 4-11-2011
	PetscOptionsGetInt(PETSC_NULL, "-IB_delta", &IB_delta, PETSC_NULL); // xyang 3-18-2011
	PetscOptionsGetInt(PETSC_NULL, "-NumIBPerLoc", &NumIBPerLoc, PETSC_NULL); // xyang 3-18-2011
	PetscOptionsGetInt(PETSC_NULL, "-NumNacellePerLoc", &NumNacellePerLoc, PETSC_NULL); // xyang 3-18-2011
	PetscOptionsGetReal(PETSC_NULL, "-reflength_wt", &reflength_wt, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-reflength_nacelle", &reflength_nacelle, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-refvel_wt", &refvel_wt, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-refvel_cfd", &refvel_cfd, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-reflength_IBDelta", &reflength_IBDelta, PETSC_NULL); // xyang 12-16-2010
	PetscOptionsGetReal(PETSC_NULL, "-tipspeedratio", &tipspeedratio, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetReal(PETSC_NULL, "-r_nacelle", &r_nacelle, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetReal(PETSC_NULL, "-L_nacelle", &L_nacelle, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetReal(PETSC_NULL, "-dh_nacelle", &dh_nacelle, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetReal(PETSC_NULL, "-loc_refvel", &loc_refvel, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetReal(PETSC_NULL, "-X_control", &X_control, PETSC_NULL); // xyang 2-12-2012
	PetscOptionsGetInt(PETSC_NULL, "-IB_wm", &IB_wm, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetInt(PETSC_NULL, "-IB_wmtmprt", &IB_wmtmprt, PETSC_NULL); // xyang 10-22-2012
	PetscOptionsGetInt(PETSC_NULL, "-temperature", &temperature, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetInt(PETSC_NULL, "-deltafunc", &deltafunc, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetReal(PETSC_NULL, "-prt_eps", &prt_eps, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-halfwidth_dfunc", &halfwidth_dfunc, PETSC_NULL); 
	PetscOptionsGetInt(PETSC_NULL, "-i_periodicIB", &i_periodicIB, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetInt(PETSC_NULL, "-j_periodicIB", &j_periodicIB, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetInt(PETSC_NULL, "-k_periodicIB", &k_periodicIB, PETSC_NULL); // xyang 6-3-2011
	PetscOptionsGetInt(PETSC_NULL, "-Force_wm", &Force_wm, PETSC_NULL); // xyang 10-24-2012
	PetscOptionsGetInt(PETSC_NULL, "-Shear_wm", &Shear_wm, PETSC_NULL); // xyang 10-24-2012
	PetscOptionsGetInt(PETSC_NULL, "-Vel_wm", &Vel_wm, PETSC_NULL); // xyang 10-24-2012
	PetscOptionsGetInt(PETSC_NULL, "-add_fluc", &add_fluctuations, PETSC_NULL); // xyang 11-03-2011
	PetscOptionsGetInt(PETSC_NULL, "-add_fluc_tmprt", &add_fluctuations_tmprt, PETSC_NULL); // xyang 11-03-2011
	PetscOptionsGetInt(PETSC_NULL, "-les_prt", &les_prt, PETSC_NULL); // xyang 11-03-2011
	PetscOptionsGetInt(PETSC_NULL, "-infRe", &infRe, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-MoveFrame", &MoveFrame, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-u_frame", &u_frame, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-v_frame", &v_frame, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-w_frame", &w_frame, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-C_iww", &C_iww, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-a_iww", &a_iww, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-lamda_iww", &lamda_iww, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-xmin", &xmin, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-xmax", &xmax, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-ymin", &ymin, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-ymax", &ymax, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-zmin", &zmin, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-zmax", &zmax, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ii_periodicWT", &ii_periodicWT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-jj_periodicWT", &jj_periodicWT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-kk_periodicWT", &kk_periodicWT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-Nx_WT", &Nx_WT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-Ny_WT", &Ny_WT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-Nz_WT", &Nz_WT, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-Sx_WT", &Sx_WT, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-Sy_WT", &Sy_WT, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-Sz_WT", &Sz_WT, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-New_wallmodel", &New_wallmodel, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-turbine", &NumberOfTurbines, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-NumberOfIBDelta", &NumberOfIBDelta, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-NumberOfNacelle", &NumberOfNacelle, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-SpongeLayer", &SpongeLayer, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-SpongeDistance", &SpongeDistance, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-MoveCylinderTest", &MoveCylinderTest, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-inflow_levelset", &inflow_levelset, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-sublevel", &sublevel, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-smoothlevel", &smoothlevel, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-dfunc_wd", &dfunc_wd, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-FixTipSpeedRatio", &FixTipSpeedRatio, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-rstart_turbinerotation", &rstart_turbinerotation, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-turbinetorquecontrol", &turbinetorquecontrol, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-turbine_prescrived_motion", &turbine_prescrived_motion, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-turbine_6dof_fsi_motion", &turbine_6dof_fsi_motion, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-turbine_prescrived_motion_heave", &turbine_prescrived_motion_heave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-turbine_prescrived_motion_pitch", &turbine_prescrived_motion_pitch, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-fixturbineangvel", &fixturbineangvel, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-fixturbinegeneratortorque", &fixturbinegeneratortorque, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-Shen_AL", &Shen_AL, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-correction3D_CH", &correction3D_CH, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-correction3D_CL", &correction3D_CL, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-smoothforce_AL", &smoothforce_AL, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-radialforce_AL", &radialforce_AL, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-AL_Noslip", &AL_Noslip, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-Shen1_AL", &Shen1_AL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-correction_ALShen1", &correction_ALShen1, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-relax_AL", &relax_AL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-cf_radialforce_AL", &cf_radialforce_AL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-c0_CL", &c0_CL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-c1_CH", &c1_CH, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-c2_CH", &c2_CH, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-c3_CH", &c3_CH, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-a_shen", &a_shen, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-b_shen", &b_shen, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-c_shen", &c_shen, PETSC_NULL);  // xyang
	PetscOptionsGetInt(PETSC_NULL, "-Prandtl_AL", &Prandtl_AL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-refangle_AL", &refangle_AL, PETSC_NULL);  // xyang
	PetscOptionsGetReal(PETSC_NULL, "-dt_bed", &dt_bed, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-particle_dens", &particle_dens, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-particle_D50", &particle_D50, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-smooth_bed", &smooth_bed, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-bed_avg", &bed_avg, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-particlevel_model", &particlevel_model, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-LS_bedConstr", &LS_bedConstr, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-Angle_repose", &Angle_repose, PETSC_NULL);
	Angle_repose	  *= M_PI;
	Angle_repose	  /= 180.;
	PetscOptionsGetReal(PETSC_NULL, "-bed_porosity", &bed_porosity, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-sb_bed", &sb_bed, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-number_smooth", &number_smooth, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-deltab", &deltab, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-rs_bed", &rs_bed, PETSC_NULL);
	// FSI
	PetscOptionsGetInt(PETSC_NULL, "-powerlawwallmodel", &powerlawwallmodel, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-prescribed_rotation", &prescribed_rotation, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-scale_velocity", &scale_velocity, PETSC_NULL);
	// add end (xiaolei)
	
	if(movefsi || rotatefsi) save_memory=0;
	sprintf(path, ".");
	PetscOptionsGetString(PETSC_NULL,"-path", path, 256, PETSC_NULL);		//  Seokkoo Kang: path for saving output; grid.dat, bcs,dat should be put there, but control.dat should exist in the current directory where job is submitted
	
	sprintf(gridfile, "grid.dat");
	PetscOptionsGetString(PETSC_NULL,"-grid", gridfile, 256, PETSC_NULL);	//  Seokkoo Kang: the name of the grid file other than grid.dat if you want

	sprintf(path_inflow, "./inflow");
	PetscOptionsGetString(PETSC_NULL,"-path_inflow", path_inflow, 256, PETSC_NULL);         //  Xiaolei Yang: path for inflow field

	int len=strlen(path);
	if(path[len-1]=='/') path[len-1]=0;
  
	extern void read_grid();
  //  read_grid();
  
	char saveoption_file[400];
	sprintf(saveoption_file, "%s/savepoints", path);
	FILE *fp=fopen(saveoption_file, "r");
	
	if(fp!=NULL) {
		i=0;
		do {
			fscanf(fp, "%d %d %d\n", &save_point[i][0], &save_point[i][1], &save_point[i][2]);
			i++;
		} while(!feof(fp));
		nsave_points=i;
		fclose(fp);
	}
	
	if(TwoD) PetscPrintf(PETSC_COMM_WORLD, "\n\n!!! 2D computation !!! \n\n");
	if(i_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nI-Periodic\n");
	if(ii_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nII-Periodic\n");
	if(j_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nJ-Periodic \n");
	if(jj_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nJJ-Periodic \n");
	if(k_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nK-Periodic \n");
	if(kk_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nKK-Periodic \n");
 
	//PetscOptionsGetReal(PETSC_NULL, "-Fr", &Fr, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-flux", &inlet_flux, PETSC_NULL);			// Seokkoo Kang: the amount of inlet flux, if not set mean bulk velocity is set to 1
	PetscOptionsGetReal(PETSC_NULL, "-imp_tol", &imp_free_tol, PETSC_NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(PETSC_NULL, "-poisson_tol", &poisson_tol, PETSC_NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(PETSC_NULL, "-les_eps", &les_eps, PETSC_NULL);		// Seokkoo Kang: small value for preventing very large Cs values in les>1
	PetscOptionsGetReal(PETSC_NULL, "-dpdz", &mean_pressure_gradient, &dpdz_set);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(PETSC_NULL, "-roughness", &roughness_size, &rough_set);	// Seokkoo Kang: roughness_size
	PetscOptionsGetReal(PETSC_NULL, "-amg_thresh", &amg_thresh, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-rho0", &rho_water, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-rho1", &rho_air, PETSC_NULL);		
	PetscOptionsGetReal(PETSC_NULL, "-mu0", &mu_water, PETSC_NULL);	
	PetscOptionsGetReal(PETSC_NULL, "-mu1", &mu_air, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-dthick", &dthick, &dthick_set);
	PetscPrintf(PETSC_COMM_WORLD, "\nrho0=%f, rho1=%f, mu0=%f, mu1=%f\n", rho_water, rho_air, mu_water, mu_air);
	PetscOptionsGetReal(PETSC_NULL, "-gx", &gravity_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-gy", &gravity_y, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-gz", &gravity_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-x_r", &(x_r), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-y_r", &(y_r), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-z_r", &(z_r), PETSC_NULL);

	PetscOptionsGetReal(PETSC_NULL, "-inlet_y", &inlet_y, &inlet_y_flag);
	PetscOptionsGetReal(PETSC_NULL, "-outlet_y", &outlet_y, PETSC_NULL);
	
	PetscOptionsGetReal(PETSC_NULL, "-inlet_z", &inlet_z, &inlet_z_flag);
	PetscOptionsGetReal(PETSC_NULL, "-outlet_z", &outlet_z, PETSC_NULL);
	
	PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);

	PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);
	
	PetscOptionsGetReal(PETSC_NULL, "-angvel", &angvel, PETSC_NULL);

  if (fish) {
    PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);
  }
  PetscPrintf(PETSC_COMM_WORLD, "tiout %d %le %le thin %d!\n",tiout, imp_atol,imp_rtol,thin);

  if (MHV) 
    L_dim=1./25.4;//.005;
  else
    L_dim=1.;
  
  if (MHV) NumberOfBodies=3;//2;
  
  if (immersed) {
    PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
    PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
	  ibm_ptr = ibm;
	  fsi_ptr = fsi;
  }

	// add begin (xiaolei)
	if (IB_delta) {
		PetscMalloc(NumberOfIBDelta*sizeof(IBMNodes), &ibm_IBDelta);
		PetscMalloc(NumberOfIBDelta*sizeof(FSInfo), &fsi_IBDelta);
	}

	if (rotor_model)  {
		PetscPrintf(PETSC_COMM_WORLD, "Plot discrete Delta functions \n");
		deltafunc_test( );
		if (rotor_model == 1)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
		}
		if (rotor_model == 2)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
			PetscMalloc(num_foiltype*sizeof(ACL), &acl);
		}
		if (rotor_model == 3)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
			PetscMalloc(num_foiltype*sizeof(ACL), &acl);
			
		}
		if (rotor_model == 4)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
		}
		if (rotor_model == 5)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_acl2ref);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_acl2ref);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_acl2ref[i].angvel_x=0;
				fsi_acl2ref[i].angvel_y=0;
				fsi_acl2ref[i].angvel_z=0;
				fsi_acl2ref[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
			PetscMalloc(num_foiltype*sizeof(ACL), &acl);
		}
		if (rotor_model == 6)  {
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_wt[i].angvel_x=0;
				fsi_wt[i].angvel_y=0;
				fsi_wt[i].angvel_z=0;
				fsi_wt[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_acl2ref);
			PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_acl2ref);
			for (i=0;i<NumberOfTurbines;i++) {
				fsi_acl2ref[i].angvel_x=0;
				fsi_acl2ref[i].angvel_y=0;
				fsi_acl2ref[i].angvel_z=0;
				fsi_acl2ref[i].angvel_axis=0.0;
			}
			PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
			PetscMalloc(num_foiltype*sizeof(ACL), &acl);
		}
	}
	if (nacelle_model) {
		PetscMalloc(NumberOfNacelle*sizeof(IBMNodes), &ibm_nacelle);
		PetscMalloc(NumberOfNacelle*sizeof(FSInfo), &fsi_nacelle);
		for (i=0;i<NumberOfNacelle;i++) {
			fsi_nacelle[i].angvel_x=0;
			fsi_nacelle[i].angvel_y=0;
			fsi_nacelle[i].angvel_z=0;
			fsi_nacelle[i].angvel_axis=0.0;
		}
	}
	// add end (xiaolei)	
	
  MG_Initial(&usermg, ibm);
  
   // Seokkoo Kang
  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
 
  VecDuplicate(user->lP, &user->lUstar);
 
// add (Toni)
	//Initialitzation of wave_momentum_source
	if(wave_momentum_source){
		extern void Initialize_wave(UserCtx *user);
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)Initialize_wave(&user[bi]);	
	}
	//End Initialitzation of wave_momentum_source	
	if(air_flow_levelset==2){
		extern void Initialize_wind(UserCtx *user);
		extern void WIND_DATA_input(UserCtx *user);				
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)Initialize_wind(&user[bi]);	
		for (bi=0; bi<block_number; bi++)WIND_DATA_input(&user[bi]);	
	}		
// End (Toni)
	if(averaging) {	// Seokkoo Kang
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		VecDuplicate(user->Ucat, &user->Ucat_sum);
		VecDuplicate(user->Ucat, &user->Ucat_cross_sum);
		VecDuplicate(user->Ucat, &user->Ucat_square_sum);
		VecSet(user->Ucat_sum,0);
		VecSet(user->Ucat_cross_sum,0);
		VecSet(user->Ucat_square_sum,0);
		  
		VecDuplicate(user->P, &user->P_sum);
		VecSet(user->P_sum,0);
		
		if(rans) {
			VecDuplicate(user->P, &user->K_sum);
			VecSet(user->K_sum, 0.);
		}
		
		if(les) {
			VecDuplicate(user->P, &user->Nut_sum);
			VecSet(user->Nut_sum, 0.);
		}

		if(averaging>=2) {
			VecDuplicate(user->P, &user->P_square_sum);
			VecSet(user->P_square_sum,0);
		}
		
		if(averaging>=3) {
			if(les) {
				VecDuplicate(user->P, &user->tauS_sum);
				VecSet(user->tauS_sum, 0);
			}
			
			VecDuplicate(user->P, &user->Udp_sum);
			VecSet(user->Udp_sum, 0);
			
			VecDuplicate(user->Ucont, &user->dU2_sum);
			VecSet(user->dU2_sum, 0);
			
			VecDuplicate(user->Ucont, &user->UUU_sum);
			VecSet(user->UUU_sum, 0);
			
			VecDuplicate(user->Ucont, &user->Vort_sum);
			VecSet (user->Vort_sum, 0);
			
			VecDuplicate(user->Ucont, &user->Vort_square_sum);
			VecSet (user->Vort_square_sum, 0);
		}
	}

  // Seokkoo Kang
  #ifdef DIRICHLET
	if(freesurface) {
		extern void Initialize_free_surface(UserCtx *user);
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)  Initialize_free_surface(&user[bi]);
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD, "\n**************** Warning !! Freesurface option not set *********************************!\n");
	}
  #endif
  
  if (immersed) {
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist));
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
				InitIBMList(&(user[bi].ibmlist[ibi]));
      }
    }
    if (MHV) {
      i=0;
      // read casing
      CMz_c=0.;//1.105;
      ibm_read_ucd(&ibm[i], i);
      // read valves
      CMz_c=4.49+0.31;
      CMy_c=.0;
      L_dim=1./28.;
      for (ibi=1; ibi<NumberOfBodies; ibi++) {
				if (ibi==2) CMy_c=-CMy_c;
				ibm_read_ucd(&ibm[ibi], ibi);
				PetscPrintf(PETSC_COMM_WORLD, "Ibm read MHV!\n");

				FsiInitialize(0, &fsi[ibi], ibi);
      }
      
      fsi[1].y_c = -0.0878; fsi[1].z_c = 4.143;//4.21;
      fsi[2].y_c =  0.0878; fsi[2].z_c = 4.143;//4.21;

      fsi[1].S_ang_n[0]= max_angle; fsi[1].S_ang_r[0]= max_angle; fsi[1].S_ang_r[0]= max_angle;  fsi[1].S_ang_rm1[0]= max_angle;
      fsi[2].S_ang_n[0]= -max_angle; fsi[2].S_ang_r[0]= -max_angle; fsi[2].S_ang_r[0]= -max_angle;  fsi[2].S_ang_rm1[0]= -max_angle;

      for (ibi=1; ibi<NumberOfBodies; ibi++) {	
				Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],0.,ibi);
      }
      PetscBarrier(PETSC_NULL);
      for (ibi=0; ibi<NumberOfBodies; ibi++) {
				ibm_surface_out(&ibm[ibi], 0, ibi);
      }
    } 
		else {
      for (i=0;i<NumberOfBodies;i++) {
				PetscPrintf(PETSC_COMM_WORLD, "Ibm read!\n");
				/*     ibm_read(ibm0); */
				ibm_read_ucd(&ibm[i], i);
				PetscBarrier(PETSC_NULL);	
				// init for fsi
				/* FsiInitialize(ibm[i].n_elmt, &fsi[i], i); */
				FsiInitialize(0, &fsi[i], i);
      }
    }
    ti = 0;
    if (rstart_flg) ti = tistart;		
  }

	if (rotor_model) {
		PetscReal cl = 1.;
		PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
		if (!my_rank) {
			FILE *fd;
			char str[256];
			sprintf(str, "%s/Turbine.inp", path);
			fd = fopen(str, "r");
			if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);
			char string[256];
			fgets(string, 256, fd);
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				fscanf(fd, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le ", &(fsi_wt[ibi].nx_tb), &(fsi_wt[ibi].ny_tb), &(fsi_wt[ibi].nz_tb), &(fsi_wt[ibi].x_c), &(fsi_wt[ibi].y_c), &(fsi_wt[ibi].z_c), &(wtm[ibi].indf_axis), &(wtm[ibi].Tipspeedratio), &(fsi_wt[ibi].J_rotation) , &(fsi_wt[ibi].r_rotor), &(fsi_wt[ibi].CP_max), &(fsi_wt[ibi].TSR_max), &(fsi_wt[ibi].angvel_fixed), &(fsi_wt[ibi].Torque_generator), &(wtm[ibi].pitch[0]), &(wtm[ibi].CT));
				double rr=sqrt(pow(fsi_wt[ibi].nx_tb,2)+pow(fsi_wt[ibi].ny_tb,2)+pow(fsi_wt[ibi].nz_tb,2));
				fsi_wt[ibi].nx_tb=fsi_wt[ibi].nx_tb/rr; 
				fsi_wt[ibi].ny_tb=fsi_wt[ibi].ny_tb/rr; 
				fsi_wt[ibi].nz_tb=fsi_wt[ibi].nz_tb/rr;
				MPI_Bcast(&(fsi_wt[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The rotating center %f %f %f \n", (fsi_wt[ibi].nx_tb), (fsi_wt[ibi].ny_tb), (fsi_wt[ibi].nz_tb));
				fsi_wt[ibi].x_c=fsi_wt[ibi].x_c/cl;
				fsi_wt[ibi].y_c=fsi_wt[ibi].y_c/cl;
				fsi_wt[ibi].z_c=fsi_wt[ibi].z_c/cl;
				MPI_Bcast(&(fsi_wt[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].Tipspeedratio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].J_rotation), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].r_rotor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].CP_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].TSR_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].angvel_fixed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].Torque_generator), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].pitch[0]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].CT), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The rotating center for %d th turbne %f %f %f \n", ibi, (fsi_wt[ibi].x_c), (fsi_wt[ibi].y_c), (fsi_wt[ibi].z_c));
				PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (wtm[ibi].indf_axis));
				PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (wtm[ibi].Tipspeedratio));
				PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_wt[ibi].J_rotation));
				PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_wt[ibi].r_rotor));
				PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_wt[ibi].CP_max));
				PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_wt[ibi].TSR_max));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_wt[ibi].angvel_fixed));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed Torque for %d th turbine  %f \n", ibi, (fsi_wt[ibi].Torque_generator));
				PetscPrintf(PETSC_COMM_WORLD, "pitch for %d th turbine  %f \n", ibi, (wtm[ibi].pitch[0]));
				PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (wtm[ibi].CT));
			}
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				fsi_wt[ibi].x_c0=fsi_wt[ibi].x_c;
				fsi_wt[ibi].y_c0=fsi_wt[ibi].y_c;
				fsi_wt[ibi].z_c0=fsi_wt[ibi].z_c;
			}
			fclose(fd);
		}
		else {
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				MPI_Bcast(&(fsi_wt[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].Tipspeedratio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].J_rotation), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].r_rotor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].CP_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].TSR_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].angvel_fixed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_wt[ibi].Torque_generator), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].pitch[0]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(wtm[ibi].CT), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			}
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				fsi_wt[ibi].x_c0=fsi_wt[ibi].x_c;
				fsi_wt[ibi].y_c0=fsi_wt[ibi].y_c;
				fsi_wt[ibi].z_c0=fsi_wt[ibi].z_c;
			}
		}
		if (rotor_model == 5) {
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				fsi_acl2ref[ibi].x_c = fsi_wt[ibi].x_c;
				fsi_acl2ref[ibi].y_c = fsi_wt[ibi].y_c;
				fsi_acl2ref[ibi].z_c = fsi_wt[ibi].z_c;
				ibm_acl2ref[ibi].indf_axis = wtm[ibi].indf_axis;
				ibm_acl2ref[ibi].Tipspeedratio = wtm[ibi].Tipspeedratio;
				fsi_acl2ref[ibi].J_rotation = fsi_wt[ibi].J_rotation;
				fsi_acl2ref[ibi].r_rotor = fsi_wt[ibi].r_rotor;
				fsi_acl2ref[ibi].CP_max = fsi_wt[ibi].CP_max;
				fsi_acl2ref[ibi].TSR_max = fsi_wt[ibi].TSR_max;
				fsi_acl2ref[ibi].angvel_fixed = fsi_wt[ibi].angvel_fixed;
				fsi_acl2ref[ibi].Torque_generator = fsi_wt[ibi].Torque_generator;
				ibm_acl2ref[ibi].pitch[0] = wtm[ibi].pitch[0];
				ibm_acl2ref[ibi].CT = wtm[ibi].CT;
				fsi_acl2ref[ibi].x_c0=fsi_wt[ibi].x_c0;
				fsi_acl2ref[ibi].y_c0=fsi_wt[ibi].y_c0;
				fsi_acl2ref[ibi].z_c0=fsi_wt[ibi].z_c0;
				fsi_acl2ref[ibi].nx_tb=fsi_wt[ibi].nx_tb;
				fsi_acl2ref[ibi].ny_tb=fsi_wt[ibi].ny_tb;
				fsi_acl2ref[ibi].nz_tb=fsi_wt[ibi].nz_tb;
				PetscPrintf(PETSC_COMM_WORLD, "Reference line file\n");
				PetscPrintf(PETSC_COMM_WORLD, "The rotating center for %d th turbne %f %f %f \n", ibi, (fsi_acl2ref[ibi].x_c), (fsi_acl2ref[ibi].y_c), (fsi_acl2ref[ibi].z_c));
				PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].indf_axis));
				PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].Tipspeedratio));
				PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].J_rotation));
				PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].r_rotor));
				PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].CP_max));
				PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].TSR_max));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_fixed));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed Torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Torque_generator));
				PetscPrintf(PETSC_COMM_WORLD, "pitch for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[0]));
				PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].CT));
			}
		}
		if (rotor_model == 6) {
			for (ibi=0;ibi<NumberOfTurbines;ibi++) {
				fsi_acl2ref[ibi].x_c = fsi_wt[ibi].x_c;
				fsi_acl2ref[ibi].y_c = fsi_wt[ibi].y_c;
				fsi_acl2ref[ibi].z_c = fsi_wt[ibi].z_c;
				ibm_acl2ref[ibi].indf_axis = wtm[ibi].indf_axis;
				ibm_acl2ref[ibi].Tipspeedratio = wtm[ibi].Tipspeedratio;
				fsi_acl2ref[ibi].J_rotation = fsi_wt[ibi].J_rotation;
				fsi_acl2ref[ibi].r_rotor = fsi_wt[ibi].r_rotor;
				fsi_acl2ref[ibi].CP_max = fsi_wt[ibi].CP_max;
				fsi_acl2ref[ibi].TSR_max = fsi_wt[ibi].TSR_max;
				fsi_acl2ref[ibi].angvel_fixed = fsi_wt[ibi].angvel_fixed;
				fsi_acl2ref[ibi].Torque_generator = fsi_wt[ibi].Torque_generator;
				ibm_acl2ref[ibi].pitch[0] = wtm[ibi].pitch[0];
				ibm_acl2ref[ibi].CT = wtm[ibi].CT;
				fsi_acl2ref[ibi].x_c0=fsi_wt[ibi].x_c0;
				fsi_acl2ref[ibi].y_c0=fsi_wt[ibi].y_c0;
				fsi_acl2ref[ibi].z_c0=fsi_wt[ibi].z_c0;
				fsi_acl2ref[ibi].nx_tb=fsi_wt[ibi].nx_tb;
				fsi_acl2ref[ibi].ny_tb=fsi_wt[ibi].ny_tb;
				fsi_acl2ref[ibi].nz_tb=fsi_wt[ibi].nz_tb;
				PetscPrintf(PETSC_COMM_WORLD, "Reference line file\n");
				PetscPrintf(PETSC_COMM_WORLD, "The rotating center for %d th turbne %f %f %f \n", ibi, (fsi_acl2ref[ibi].x_c), (fsi_acl2ref[ibi].y_c), (fsi_acl2ref[ibi].z_c));
				PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].indf_axis));
				PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].Tipspeedratio));
				PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].J_rotation));
				PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].r_rotor));
				PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].CP_max));
				PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].TSR_max));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_fixed));
				PetscPrintf(PETSC_COMM_WORLD, "Fixed Torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Torque_generator));
				PetscPrintf(PETSC_COMM_WORLD, "pitch for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[0]));
				PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].CT));
			}
		}
		PetscBarrier(PETSC_NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Turbines read!\n");
		if (rotor_model == 1) {
			double reflength = reflength_wt;
			char fname[80];
			sprintf(fname,"acddata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
		}
		if (rotor_model == 2) {
			double reflength = reflength_wt;
			char fname[80];
			sprintf(fname,"acl2data000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
			sprintf(fname,"Urefdata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
			Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
			airfoil_ACL(acl, wtm,  fsi_wt);
		}
		if (rotor_model == 3) {
			double reflength = reflength_wt;
			for (i=0;i<NumberOfTurbines;i++) ACL_read_ucd(&wtm[i], i, &fsi_wt[i], reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
			char fname[80];
			sprintf(fname,"Urefdata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
			Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
			airfoil_ACL(acl, wtm,  fsi_wt);
		}
		if (rotor_model == 4) {
			double reflength = reflength_wt;
			char fname[80];
			sprintf(fname,"acddata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
			sprintf(fname,"Urefdata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
			Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
		}
		if (rotor_model == 5) {
			double reflength = reflength_wt;
			//ACD_read(&wtm[i], i, &fsi_wt[i], 0);
			char fname[80];
			sprintf(fname,"acsdata000");
			//for (i=0;i<NumberOfTurbines;i++) ACL_read_ucd(&wtm[i], i, &fsi_wt[i]);
			for (i=0;i<NumberOfTurbines;i++) surface_read_xpatch(&wtm[i], i, &fsi_wt[i], fname, reflength);
			//for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname);
			PetscPrintf(PETSC_COMM_WORLD, "acl2 pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines);  
			//disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname);
			sprintf(fname,"Urefdata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
			Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
			PetscPrintf(PETSC_COMM_WORLD, "Read uref line!\n");
			for (i=0;i<NumberOfTurbines;i++) ACL_read_ucd(&ibm_acl2ref[i], i, &fsi_acl2ref[i], reflength);  // 20140807
			PetscPrintf(PETSC_COMM_WORLD, "Uref line pre-processing!\n");
			Pre_process(&(user[0]), ibm_acl2ref, NumberOfTurbines); // 20140807
			calc_s2l(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref, NumberOfTurbines);
			airfoil_ACL(acl, wtm,  fsi_wt);
			airfoil_ACL(acl, ibm_acl2ref,  fsi_acl2ref);
		}
		if (rotor_model == 6) {
			double reflength = reflength_wt;
			for (i=0;i<NumberOfTurbines;i++) ACL_read_ucd(&wtm[i], i, &fsi_wt[i], reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
			char fname[80];
			sprintf(fname,"Urefdata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
			Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
			airfoil_ACL(acl, wtm,  fsi_wt);
			PetscPrintf(PETSC_COMM_WORLD, "Read uref line!\n");
			for (i=0;i<NumberOfTurbines;i++) ACL_read_ucd(&ibm_acl2ref[i], i, &fsi_acl2ref[i], reflength);  // 20140807
			PetscPrintf(PETSC_COMM_WORLD, "Uref line pre-processing!\n");
			Pre_process(&(user[0]), ibm_acl2ref, NumberOfTurbines); // 20140807
		}
		PetscBarrier(PETSC_NULL);
		ti = 0;
		if (rstart_flg) ti = tistart;
	}	
	if (IB_delta) {
		PetscReal cl = 1.;
		PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
		if (!my_rank) {
			FILE *fd;
			char str[256];
			sprintf(str, "%s/IBDelta.inp", path);
			fd = fopen(str, "r");
			if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);
			char string[256];
			fgets(string, 256, fd);
			for (ibi=0;ibi<NumberOfIBDelta;ibi++) {
				fscanf(fd, "%le %le %le %le %le %le %le %le %le %le", &(fsi_IBDelta[ibi].nx_tb), &(fsi_IBDelta[ibi].ny_tb), &(fsi_IBDelta[ibi].nz_tb), &(fsi_IBDelta[ibi].x_c), &(fsi_IBDelta[ibi].y_c), &(fsi_IBDelta[ibi].z_c), &(ibm_IBDelta[ibi].CD_bluff), &(ibm_IBDelta[ibi].indf_axis), &(ibm_IBDelta[ibi].indf_tangent), &(fsi_IBDelta[ibi].angvel_axis));
				double rr=sqrt(pow(fsi_IBDelta[ibi].nx_tb,2)+pow(fsi_IBDelta[ibi].ny_tb,2)+pow(fsi_IBDelta[ibi].nz_tb,2));
				fsi_IBDelta[ibi].nx_tb=fsi_IBDelta[ibi].nx_tb/rr; 
				fsi_IBDelta[ibi].ny_tb=fsi_IBDelta[ibi].ny_tb/rr; 
				fsi_IBDelta[ibi].nz_tb=fsi_IBDelta[ibi].nz_tb/rr;
				MPI_Bcast(&(fsi_IBDelta[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The directions of IBDelta %f %f %f \n", (fsi_IBDelta[ibi].nx_tb), (fsi_IBDelta[ibi].ny_tb), (fsi_IBDelta[ibi].nz_tb));
				fsi_IBDelta[ibi].x_c=fsi_IBDelta[ibi].x_c/cl;
				fsi_IBDelta[ibi].y_c=fsi_IBDelta[ibi].y_c/cl;
				fsi_IBDelta[ibi].z_c=fsi_IBDelta[ibi].z_c/cl;
				MPI_Bcast(&(fsi_IBDelta[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].CD_bluff), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].indf_tangent), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The Locations for %d th  IBDelta %f %f %f \n", ibi, (fsi_IBDelta[ibi].x_c), (fsi_IBDelta[ibi].y_c), (fsi_IBDelta[ibi].z_c));
				PetscPrintf(PETSC_COMM_WORLD, "The drag coefficient for %d th IBDelta body %f \n", ibi, (ibm_IBDelta[ibi].CD_bluff));
				PetscPrintf(PETSC_COMM_WORLD, "The axial induction factor for %d th IBDelta body %f \n", ibi, (ibm_IBDelta[ibi].indf_axis));
				PetscPrintf(PETSC_COMM_WORLD, "The tangential induction factor for %d th IBDelta body %f \n", ibi, (ibm_IBDelta[ibi].indf_tangent));
				PetscPrintf(PETSC_COMM_WORLD, "The angular velocity for %d th IBDelta body %f \n", ibi, (fsi_IBDelta[ibi].angvel_axis));
			}
			for (ibi=0;ibi<NumberOfIBDelta;ibi++) {
				fsi_IBDelta[ibi].x_c0=fsi_IBDelta[ibi].x_c;
				fsi_IBDelta[ibi].y_c0=fsi_IBDelta[ibi].y_c;
				fsi_IBDelta[ibi].z_c0=fsi_IBDelta[ibi].z_c;
			}
			fclose(fd);
		}
		else {
			for (ibi=0;ibi<NumberOfIBDelta;ibi++) {
				MPI_Bcast(&(fsi_IBDelta[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

				MPI_Bcast(&(fsi_IBDelta[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].CD_bluff), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_IBDelta[ibi].indf_tangent), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_IBDelta[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			}

			for (ibi=0;ibi<NumberOfIBDelta;ibi++) {
				fsi_IBDelta[ibi].x_c0=fsi_IBDelta[ibi].x_c;
				fsi_IBDelta[ibi].y_c0=fsi_IBDelta[ibi].y_c;
				fsi_IBDelta[ibi].z_c0=fsi_IBDelta[ibi].z_c;
			}
		}
		PetscPrintf(PETSC_COMM_WORLD, "IBDeltas read  11!\n");
		double reflength = reflength_IBDelta;
		int ipt;
		int NumLoc=NumberOfIBDelta/NumIBPerLoc;
		for (ibi=0;ibi<NumLoc;ibi++) 
		for (ipt=0;ipt<NumIBPerLoc;ipt++) {
			int iname=ibi*NumIBPerLoc+ipt; 
			char fname[80];
			sprintf(fname,"ibmDelta%3.3d", ipt);
			disk_read_ucd(&ibm_IBDelta[iname], iname, &fsi_IBDelta[iname], 0, fname, reflength);	
		}
		Pre_process(&(user[0]), ibm_IBDelta, NumberOfIBDelta); // xyang 12-13-2010
		ti = 0;
		if (rstart_flg) ti = tistart;
		for (bi=0; bi<block_number; bi++) {
			for (ibi=0;ibi<NumberOfIBDelta;ibi++) {
				//Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
			}
		}
  }
	if (nacelle_model) {
		PetscReal cl = 1.;
		PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
		if (!my_rank) {
			FILE *fd;
			char str[256];
			sprintf(str, "%s/Nacelle.inp", path);
			fd = fopen(str, "r");
			if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);
			char string[256];
			fgets(string, 256, fd);
			for (ibi=0;ibi<NumberOfNacelle;ibi++) {
				fscanf(fd, "%le %le %le %le %le %le %le %d %le %le %le %le %le %le %le", &(fsi_nacelle[ibi].nx_tb), &(fsi_nacelle[ibi].ny_tb), &(fsi_nacelle[ibi].nz_tb), &(fsi_nacelle[ibi].x_c), &(fsi_nacelle[ibi].y_c), &(fsi_nacelle[ibi].z_c), &(fsi_nacelle[ibi].angvel_axis), &(fsi_nacelle[ibi].rotate_alongaxis), &(ibm_nacelle[ibi].axialforcecoefficient), &(ibm_nacelle[ibi].tangentialforcecoefficient), &(ibm_nacelle[ibi].axialprojectedarea), &(ibm_nacelle[ibi].tangentialprojectedarea), &(ibm_nacelle[ibi].pressure_factor), &(ibm_nacelle[ibi].friction_factor), &(ibm_nacelle[ibi].dh));
				double rr=sqrt(pow(fsi_nacelle[ibi].nx_tb,2)+pow(fsi_nacelle[ibi].ny_tb,2)+pow(fsi_nacelle[ibi].nz_tb,2));
				fsi_nacelle[ibi].nx_tb=fsi_nacelle[ibi].nx_tb/rr; 
				fsi_nacelle[ibi].ny_tb=fsi_nacelle[ibi].ny_tb/rr; 
				fsi_nacelle[ibi].nz_tb=fsi_nacelle[ibi].nz_tb/rr;
				MPI_Bcast(&(fsi_nacelle[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The directions of nacelle %f %f %f \n", (fsi_nacelle[ibi].nx_tb), (fsi_nacelle[ibi].ny_tb), (fsi_nacelle[ibi].nz_tb));
				fsi_nacelle[ibi].x_c=fsi_nacelle[ibi].x_c/cl;
				fsi_nacelle[ibi].y_c=fsi_nacelle[ibi].y_c/cl;
				fsi_nacelle[ibi].z_c=fsi_nacelle[ibi].z_c/cl;
				MPI_Bcast(&(fsi_nacelle[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].rotate_alongaxis), 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].axialforcecoefficient), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].tangentialforcecoefficient), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].axialprojectedarea), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].tangentialprojectedarea), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].pressure_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].friction_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].dh), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				PetscPrintf(PETSC_COMM_WORLD, "The Locations for %d th  nacelle %f %f %f \n", ibi, (fsi_nacelle[ibi].x_c), (fsi_nacelle[ibi].y_c),
				(fsi_nacelle[ibi].z_c));
				PetscPrintf(PETSC_COMM_WORLD, "The angular velocity for %d th nacelle body %f \n", ibi, (fsi_nacelle[ibi].angvel_axis));
				PetscPrintf(PETSC_COMM_WORLD, "Rotate the nacelle along the axis? %d \n",(fsi_nacelle[ibi].rotate_alongaxis) );
				PetscPrintf(PETSC_COMM_WORLD, "The axial force coefficient of %d th nacelle body %f \n", ibi, (ibm_nacelle[ibi].axialforcecoefficient));
				PetscPrintf(PETSC_COMM_WORLD, "The tangential force coefficientt of %d th nacelle body %f \n", ibi, (ibm_nacelle[ibi].tangentialforcecoefficient));
				PetscPrintf(PETSC_COMM_WORLD, "The axial projected area of %d th nacelle body %f \n", ibi, (ibm_nacelle[ibi].axialprojectedarea));
				PetscPrintf(PETSC_COMM_WORLD, "The tangential projected area of %d th nacelle body %f \n", ibi, (ibm_nacelle[ibi].tangentialprojectedarea));
				PetscPrintf(PETSC_COMM_WORLD, "The pressure factor of %d th nacelle %f \n", ibi, (ibm_nacelle[ibi].pressure_factor));
				PetscPrintf(PETSC_COMM_WORLD, "The friction factor of %d th nacelle %f \n", ibi, (ibm_nacelle[ibi].friction_factor));
				PetscPrintf(PETSC_COMM_WORLD, "The wall-normal thickness of %d th nacelle mesh %f \n", ibi, (ibm_nacelle[ibi].dh));
			}
			for (ibi=0;ibi<NumberOfNacelle;ibi++) {
				fsi_nacelle[ibi].x_c0=fsi_nacelle[ibi].x_c;
				fsi_nacelle[ibi].y_c0=fsi_nacelle[ibi].y_c;
				fsi_nacelle[ibi].z_c0=fsi_nacelle[ibi].z_c;
			}
			fclose(fd);
		}
		else {
			for (ibi=0;ibi<NumberOfNacelle;ibi++) {
				MPI_Bcast(&(fsi_nacelle[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(fsi_nacelle[ibi].rotate_alongaxis), 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].axialforcecoefficient), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].tangentialforcecoefficient), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].axialprojectedarea), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].tangentialprojectedarea), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].pressure_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].friction_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&(ibm_nacelle[ibi].dh), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			}
			for (ibi=0;ibi<NumberOfNacelle;ibi++) {
				fsi_nacelle[ibi].x_c0=fsi_nacelle[ibi].x_c;
				fsi_nacelle[ibi].y_c0=fsi_nacelle[ibi].y_c;
				fsi_nacelle[ibi].z_c0=fsi_nacelle[ibi].z_c;
			}
		}
   	PetscPrintf(PETSC_COMM_WORLD, "nacells read  11!\n");
		double reflength = reflength_nacelle;
		int ipt;
		int NumLoc=(int)NumberOfNacelle/(int)NumNacellePerLoc;
		for (ibi=0;ibi<NumLoc;ibi++) 
		for (ipt=0;ipt<(int)NumNacellePerLoc;ipt++) {
			int iname=ibi*(int)NumNacellePerLoc+ipt; 
			char fname[80];
			sprintf(fname,"nacelle%3.3d", ipt);
			disk_read_ucd(&ibm_nacelle[iname], iname, &fsi_nacelle[iname], 0, fname, reflength);	
		}
		Pre_process(&(user[0]), ibm_nacelle, NumberOfNacelle); // xyang 12-13-2010
		ti = 0;
		if (rstart_flg) ti = tistart;
  }

  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (rstart_flg) {
    ti = tistart; tistart++;
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi]));
      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);
      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DAGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
      DAGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
      DAGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      DAGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      Contra2Cart(&(user[bi]));
			if (rstart_fsi) {
				for (ibi=0;ibi<NumberOfBodies;ibi++) {
					if(!rotatefsi) FSI_DATA_Input(&fsi[ibi],ibi);
					if (movefsi && !fsi_6dof) {
						Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi]);	
						for (i=0;i<6;i++){
							fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
							fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
						}
						for (i=0; i<ibm[ibi].n_v; i++) {
							ibm[ibi].uold[i].x = fsi[ibi].S_real[1];
							ibm[ibi].uold[i].y = fsi[ibi].S_real[3];
							ibm[ibi].uold[i].z = fsi[ibi].S_real[5];
						}
						for (i=0; i<ibm[ibi].n_v; i++) {
							ibm[ibi].urm1[i].x = fsi[ibi].S_realm1[1];
							ibm[ibi].urm1[i].y = fsi[ibi].S_realm1[3];
							ibm[ibi].urm1[i].z = fsi[ibi].S_realm1[5];
						}
					}
					else if (movefsi && fsi_6dof) {
						Elmt_Move_FSI_ROT_TRANS(&fsi[ibi], &ibm[ibi],user->dt,0);			
						for (i=0; i<ibm[ibi].n_v; i++) {
							double rot_angle;
							//rotate
							rotate_xyz6dof (ti, user->dt, fsi[ibi].S_ang_r[0],fsi[ibi].S_ang_r[2],fsi[ibi].S_ang_r[4], x_r, y_r, z_r, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &ibm->x_bp_o[i], &ibm->y_bp_o[i], &ibm->z_bp_o[i], &rot_angle);
							//translate
							ibm->x_bp_o[i] = ibm->x_bp_o[i]+(fsi[ibi].S_real[0]);//-FSinfo->S_real[0]);
							ibm->y_bp_o[i] = ibm->y_bp_o[i]+(fsi[ibi].S_real[2]);//-FSinfo->S_real[2]);
							ibm->z_bp_o[i] = ibm->z_bp_o[i]+(fsi[ibi].S_real[4]);//-FSinfo->S_real[4]); 							
							ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt ;
							ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt ;
							ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt ;
							ibm[ibi].uold[i] = ibm[ibi].u[i];
							ibm[ibi].urm1[i] = ibm[ibi].u[i];					
						}	
						PetscInt rank;
						MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
						if (!rank) {
							FILE *f;
							char filen[80];
							sprintf(filen, "RESTARTsurface%3.3d.dat",ti);
							f = fopen(filen, "w");
							PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
							PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", ibm[ibi].n_v, ibm[ibi].n_elmt);
							for (i=0; i<ibm[ibi].n_v; i++) {
								PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
							}
							for (i=0; i<ibm[ibi].n_elmt; i++) {
								PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
							}
							fclose(f);					
						}
						for (i=0;i<6;i++){
							fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
							fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
							fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
							fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];							
						}
						fsi[ibi].F_x_real=fsi[ibi].F_x;
						fsi[ibi].F_y_real=fsi[ibi].F_y;
						fsi[ibi].F_z_real=fsi[ibi].F_z;	   
						fsi[ibi].M_x_rm3=fsi[ibi].M_x;
						fsi[ibi].M_y_rm3=fsi[ibi].M_y;
						fsi[ibi].M_z_rm3=fsi[ibi].M_z;
						fsi[ibi].M_x_rm2=fsi[ibi].M_x;
						fsi[ibi].M_y_rm2=fsi[ibi].M_y;
						fsi[ibi].M_z_rm2=fsi[ibi].M_z;
						fsi[ibi].M_x_real=fsi[ibi].M_x;
						fsi[ibi].M_y_real=fsi[ibi].M_y;
						fsi[ibi].M_z_real=fsi[ibi].M_z;						
					}
					if (rotatefsi|| MHV) {
						fsi[ibi].x_c = x_r;
						fsi[ibi].y_c = y_r;
						fsi[ibi].z_c = z_r;
						if(ibi==0) {
							Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
						}
						else {
							for (i=0; i<ibm[ibi].n_v; i++) {
								ibm[ibi].u[i].x = 0;
								ibm[ibi].u[i].y = 0;
								ibm[ibi].u[i].z = 0;
								ibm[ibi].uold[i] = ibm[ibi].u[i];
								ibm[ibi].urm1[i] = ibm[ibi].u[i];
							}
						}
						// if read ti, then will start for ti+1
						for (i=0;i<6;i++){
							fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
							fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
						}
						fsi[ibi].F_x_real=fsi[ibi].F_x;
						fsi[ibi].F_y_real=fsi[ibi].F_y;
						fsi[ibi].F_z_real=fsi[ibi].F_z;	   
						fsi[ibi].M_x_rm3=fsi[ibi].M_x;
						fsi[ibi].M_y_rm3=fsi[ibi].M_y;
						fsi[ibi].M_z_rm3=fsi[ibi].M_z;
						fsi[ibi].M_x_rm2=fsi[ibi].M_x;
						fsi[ibi].M_y_rm2=fsi[ibi].M_y;
						fsi[ibi].M_z_rm2=fsi[ibi].M_z;
						fsi[ibi].M_x_real=fsi[ibi].M_x;
						fsi[ibi].M_y_real=fsi[ibi].M_y;
						fsi[ibi].M_z_real=fsi[ibi].M_z;
						/*
						PetscReal rx,ry,rz;
						for (i=0; i<ibm[ibi].n_v; i++) {
							rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
							ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
							rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;      
							ibm[ibi].u[i].x =   ry*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[3]*rz  ;
							ibm[ibi].u[i].y =-( rx*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[1]*rz );
							ibm[ibi].u[i].z =   rx*fsi[ibi].S_ang_n[3]-fsi[ibi].S_ang_n[1]*ry  ;     
							ibm[ibi].uold[i].x =   ry*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[3]*rz  ;
							ibm[ibi].uold[i].y =-( rx*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[1]*rz );
							ibm[ibi].uold[i].z =   rx*fsi[ibi].S_ang_r[3]-fsi[ibi].S_ang_r[1]*ry  ;      
							ibm[ibi].urm1[i].x =   ry*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[3]*rz  ;
							ibm[ibi].urm1[i].y =-( rx*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[1]*rz );
							ibm[ibi].urm1[i].z =   rx*fsi[ibi].S_ang_rm1[3]-fsi[ibi].S_ang_rm1[1]*ry  ;
						}*/
					}
				}//ibi
			} // if rstart fsi
		}// bi
	} // if rstart

// do the search once if elmt is not moving!
  if (immersed) {
    for (level = usermg.mglevels-1; level>=usermg.mglevels-1; level--) {
      user = usermg.mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {
				for (ibi=0;ibi<NumberOfBodies;ibi++) {
					PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA %d \n", ibi);
					ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
				}
				PetscBarrier(PETSC_NULL);
				PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP\n");
				ibm_interpolation_advanced(&user[bi]);
      }
    }
  }

  // Copy Ucont to Ucont_o for the finest level
  for (bi=0; bi<block_number; bi++) {
    //VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o));
    ti = 0;
    if (rstart_flg) ti = tistart;
    if(ti==tistart && ti==0 && levelset) {
			Levelset_Function_IC(&user[bi]);
			DAGlobalToLocalBegin(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
			DAGlobalToLocalEnd(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
			VecCopy(user[bi].Levelset, user[bi].Levelset_o);
    }
    if(ti==tistart) Calc_Inlet_Area(&user[bi]);
    if (ti==0) {
			VecSet(user[bi].Ucont,0.);
			VecSet(user[bi].lUcont,0.);
			VecSet(user[bi].Ucont_o,0.);
			VecSet(user[bi].lUcont_o,0.);
			VecSet(user[bi].Ucat,0.);
			VecSet(user[bi].lUcat,0.);
			VecSet(user[bi].P,0.);
			VecSet(user[bi].lP,0.);
			//if(initialzero) PetscPrintf(PETSC_COMM_WORLD, "\nInitial Guess is Zero !\n");
			//else 
			SetInitialGuessToOne(&(user[bi]));
			Contra2Cart(&(user[bi]));
			DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
			DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
    }

    VecCopy(user[bi].Ucont, user[bi].Ucont_o);
    //VecCopy(user[bi].Ucont, user[bi].Ucont_rm2);	// allocate at init.c
    VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);
    VecCopy(user[bi].Ucat, user[bi].Ucat_o);
    VecCopy(user[bi].P, user[bi].P_o);
    DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
  }
  PetscBarrier(PETSC_NULL);
  PetscInt tisteps = 1000000;
  PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg);
  if (tistart==0) tisteps ++;
	
/* ==================================================================================             */
/*   pysical time Step Loop */
  for (ti = tistart; ti<tistart + tisteps; ti++) {
    PetscPrintf(PETSC_COMM_WORLD, "Time %d\n", ti);
    if (inletprofile==3) {
			if (MHV && (fsi[1].S_ang_n[0]<0.8*max_angle || fsi[2].S_ang_n[0]>-0.8*max_angle)) 
			angle=angle+1;
			else
			angle=0.;
			fluxin(&(usermg.mgctx[usermg.mglevels-1].user[0]));
    }
    /* ==================================================================================             */
    /*     Strong-Coupling (SC) Loop */
    DoSCLoop= PETSC_TRUE ; itr_sc = 0;
    while (DoSCLoop) {
      itr_sc++;
      PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);
      /*     Structral Solver! */
      if (immersed)
      Struc_Solver(&usermg, ibm, fsi, itr_sc,tistart, &DoSCLoop);
      else
      DoSCLoop = PETSC_FALSE;
			/*       /\*     Structral Solver! *\/ */
			/*       if (immersed) */
			/*       Struc_predictor(&usermg, ibm, fsi, itr_sc,tistart, &DoSCLoop); */
			/*       else */
			/*       DoSCLoop = PETSC_FALSE; */
      /*     Flow Solver! */
      if(levelset) Calc_Inlet_Area(&(usermg.mgctx[usermg.mglevels-1].user[0]));
      //Flow_Solver(&usermg, ibm, fsi,itr_sc);
			//Flow_Solver(&usermg, ibm, fsi,itr_sc, wtm, acl, fsi_wt, ibm_ACD, fsi_IBDelta, ibm_IBDelta);
      Flow_Solver(&usermg, ibm, fsi, itr_sc, wtm, acl, fsi_wt, ibm_ACD, fsi_IBDelta, ibm_IBDelta, ibm_acl2ref, fsi_acl2ref, ibm_nacelle, fsi_nacelle);
		
      if(rotatefsi || movefsi || NumberOfBodies==2) for (ibi=0;ibi<NumberOfBodies;ibi++) ibm_surface_out_with_pressure(&ibm[ibi], ibi);
      
    }// End of while SC loop
    /* ==================================================================================             */
		/*  put the time accuracy coefficient back to 1.5 
				after the 1st real-time step */
		/*     COEF_TIME_ACCURACY=1.5; */
		/* ==================================================================================             */
		/*     Save the old values (at ti) for later */
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      if (immersed) {
				VecCopy(user[bi].Nvert, user[bi].Nvert_o);
				DAGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
				DAGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      }
      //VecCopy(user[bi].Ucont_rm1, user[bi].Ucont_rm2);
			if(levelset) VecCopy(user[bi].Levelset, user[bi].Levelset_o);
      VecCopy(user[bi].Ucont_o, user[bi].Ucont_rm1);
      VecCopy(user[bi].Ucont, user[bi].Ucont_o);
      VecCopy(user[bi].P, user[bi].P_o);
      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      //seokkoo
      DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
      DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
    }

    if (immersed && (movefsi || rotatefsi || cop || fish || MHV)){
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
				for (i=0; i<ibm[ibi].n_v; i++) {
					ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
					ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
					ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];
					ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;
					ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
					ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;
					ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
					ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
					ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
				}
				for (i=0;i<6;i++){
					fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
					fsi[ibi].S_real[i]=fsi[ibi].S_new[i];

					fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
					fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
				}
				fsi[ibi].F_x_real=fsi[ibi].F_x;
				fsi[ibi].F_y_real=fsi[ibi].F_y;
				fsi[ibi].F_z_real=fsi[ibi].F_z;
				fsi[ibi].M_x_rm3=fsi[ibi].M_x_rm2;
				fsi[ibi].M_y_rm3=fsi[ibi].M_y_rm2;
				fsi[ibi].M_z_rm3=fsi[ibi].M_z_rm2;
				fsi[ibi].M_x_rm2=fsi[ibi].M_x_real;
				fsi[ibi].M_y_rm2=fsi[ibi].M_y_real;
				fsi[ibi].M_z_rm2=fsi[ibi].M_z_real;
				fsi[ibi].M_x_real=fsi[ibi].M_x;
				fsi[ibi].M_y_real=fsi[ibi].M_y;
				fsi[ibi].M_z_real=fsi[ibi].M_z;
      } //ibi
    }

		/* ==================================================================================             */
				
		////////////////////////////////////---------------------------
		/*     if (ti == (ti/100)*100) */
		/*       Ucont_P_Binary_Output(&user); */

  } // ti (physical time) loop
	/* ==================================================================================             */
  PetscPrintf(PETSC_COMM_WORLD, "\n\n ******* Finished computation ti=%d ******* \n\n", ti);
  MG_Finalize(&usermg);
  PetscFinalize();
	/* ==================================================================================             */
  return(0);
}