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
PetscBool dpdz_set=PETSC_FALSE;
int save_inflow=0;
int save_inflow_period=100;
int save_inflow_minus=0;
int ti_lastsave=0;
int localstep=1;
int inflow_recycle_perioid=20000;
int save_memory=1;
int ibm_search=0;
PetscBool rough_set=PETSC_FALSE;
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
PetscBool dthick_set=PETSC_FALSE;
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

PetscBool	rstart_flg=PETSC_FALSE;
PetscBool	inlet_y_flag=PETSC_FALSE;
PetscBool	inlet_z_flag=PETSC_FALSE;

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
	VecLoad(user->Ucont, viewer);
	PetscViewerDestroy(&viewer);

	PetscBarrier(NULL);

	PetscOptionsClearValue(NULL, "-vecload_block_size");

	sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->P, viewer);
	PetscViewerDestroy(&viewer);

	sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->Nvert_o, viewer);
	PetscViewerDestroy(&viewer);
  
	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);	// Seokkoo Kang
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->Ucat, viewer);
	PetscViewerDestroy(&viewer);
  
	if(!immersed) {
          VecSet(user->Nvert, 0.);
          VecSet(user->Nvert_o, 0.);
	}

	DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
	DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
	
	DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

	VecCopy(user->Ucont, user->Ucont_o);
	DMGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
	DMGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
  
	DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
	DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
  
	DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);
        DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);

	DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);
	DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);
  
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
			PetscBarrier(NULL);
			sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->Ucat_sum, viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->Ucat_cross_sum, viewer);
			PetscViewerDestroy(&viewer);
		  
			sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->Ucat_square_sum, viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%s/sp_%06d_%1d.dat", path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->P_sum, viewer);
			PetscViewerDestroy(&viewer);
			
			if(les) {
			  sprintf(filen, "%s/snut_%06d_%1d.dat", path, ti, user->_this);
			  if( file_exist(filen) ) {
			    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			    VecLoad(user->Nut_sum, viewer);
			    PetscViewerDestroy(&viewer);
			  }
			}
			
			if(rans) {
			  sprintf(filen, "%s/sk_%06d_%1d.dat", path, ti, user->_this);
			  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			  VecLoad(user->K_sum, viewer);
			  PetscViewerDestroy(&viewer);
			}

			if(averaging>=2) {
				sprintf(filen, "%s/sp2_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->P_square_sum, viewer);
				PetscViewerDestroy(&viewer);
				/*
				sprintf(filen, "%s/sp1_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->P_cross_sum, viewer);
				PetscViewerDestroy(&viewer);
				*/
			}
			
			if(averaging>=3) {
				if(les) {
					sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
					if( file_exist(filen) ) {
					  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					  VecLoad(user->tauS_sum, viewer);
					  PetscViewerDestroy(&viewer);
					}
					else PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !\n", filen);
				}
				
				sprintf(filen, "%s/su3_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->Udp_sum, viewer);
				PetscViewerDestroy(&viewer);

				sprintf(filen, "%s/su4_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->dU2_sum, viewer);
				PetscViewerDestroy(&viewer);

				sprintf(filen, "%s/su5_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->UUU_sum, viewer);
				PetscViewerDestroy(&viewer);

				sprintf(filen, "%s/svo_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->Vort_sum, viewer);
				PetscViewerDestroy(&viewer);
				
				sprintf(filen, "%s/svo2_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user->Vort_square_sum, viewer);
				PetscViewerDestroy(&viewer);
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
		
			PetscBarrier(NULL);
			
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->Levelset, viewer);
			PetscViewerDestroy(&viewer);
			
			VecCopy (user->Levelset, user->Levelset_o);
			DMGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
			DMGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
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
		
			PetscBarrier(NULL);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(Cs, viewer);
			PetscViewerDestroy(&viewer);
		}
		
		DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
		DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
		
		VecDestroy(&Cs);
	}
  
	if(rans) {
		// K-Omega
		sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
		if( file_exist(filen) ) {
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->K_Omega, viewer);
			PetscViewerDestroy(&viewer);
		}
		else {
		  //K_Omega_IC(user);
			PetscPrintf(PETSC_COMM_WORLD, "\nInitializing K-omega ... \n\n");
                        K_Omega_IC(user);
                        VecSet(user->lNu_t, user->ren);
		}
		
		VecCopy(user->K_Omega, user->K_Omega_o);
			
		DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
		DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
			
		DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
		DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
		
		if(rans==3) {
		  // distance
		  sprintf(filen, "%s/Distance_%1d.dat", path, user->_this);
		  if( file_exist(filen) ) {
		    PetscPrintf(PETSC_COMM_WORLD, "Reading %s !\n", filen);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user->Distance, viewer);
			PetscViewerDestroy(&viewer);
		  }
		  
		  else {
		    PetscPrintf(PETSC_COMM_WORLD, "File %s does not exist. Recalculating distance function !\n", filen);
		    /*Compute_Distance_Function(user);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Distance, viewer);
			PetscViewerDestroy(&viewer);
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
	
	PetscBarrier(NULL);
	
	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucat, viewer);
	PetscViewerDestroy(&viewer);
	sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(NULL);
  
	return 0;
}

int delete_count=0;

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user)
{
	PetscViewer	viewer;
	char filen[80];
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	
	PetscBarrier(NULL);
	
	sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucont, viewer);
	PetscViewerDestroy(&viewer);
	sprintf(filen, "%s/vfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(NULL);

	sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Ucat, viewer);
	PetscViewerDestroy(&viewer);
	sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(NULL);

	sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->P, viewer);
	PetscViewerDestroy(&viewer);
	sprintf(filen, "%s/pfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(NULL);
	
	if(qcr) {
		Vec Q;
		VecDuplicate(user->P, &Q);
		Compute_Q(user,  Q);
		
		sprintf(filen, "%s/qfield%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(Q, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/qfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		VecDestroy(&Q);
		PetscBarrier(NULL);
	}
	
	sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
	VecView(user->Nvert, viewer);
	PetscViewerDestroy(&viewer);
	sprintf(filen, "%s/nvfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
	
	PetscBarrier(NULL);
  
	if(averaging && ti!=0) {	// Seokkoo Kang
		sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_sum, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/su0_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(NULL);
		  
		sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_cross_sum, viewer);
		PetscViewerDestroy(&viewer);  
		sprintf(filen, "%s/su1_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(NULL);
		
		sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Ucat_square_sum, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/su2_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(NULL);
		  
		sprintf(filen, "%s/sp_%06d_%1d.dat",path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->P_sum, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/sp_%06d_%1d.dat.info",path, ti, user->_this);	if(!rank) unlink(filen);
		
		if(les) {
		  sprintf(filen, "%s/snut_%06d_%1d.dat",path, ti, user->_this);
		  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		  VecView(user->Nut_sum, viewer);
		  PetscViewerDestroy(&viewer);
		  sprintf(filen, "%s/snut_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}

		if(rans) {
		  sprintf(filen, "%s/sk_%06d_%1d.dat",path, ti, user->_this);
		  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		  VecView(user->K_sum, viewer);
		  PetscViewerDestroy(&viewer);
		  sprintf(filen, "%s/sk_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}









		if(averaging>=2) {
			sprintf(filen, "%s/sp2_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->P_square_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/sp2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			/*
			sprintf(filen, "%s/sp1_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->P_cross_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/sp1_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			*/
		}
		
		if(averaging>=3) {
			if(les) {
				sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
				VecView(user->tauS_sum, viewer);
				PetscViewerDestroy(&viewer);
				sprintf(filen, "%s/stauS_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);
			}
				
			sprintf(filen, "%s/su3_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Udp_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/su3_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/su4_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->dU2_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/su4_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/su5_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->UUU_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/su5_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

			sprintf(filen, "%s/svo_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Vort_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/svo_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
			
			sprintf(filen, "%s/svo2_%06d_%1d.dat",path, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
			VecView(user->Vort_square_sum, viewer);
			PetscViewerDestroy(&viewer);
			sprintf(filen, "%s/svo2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
		}
		
		PetscBarrier(NULL);
	}
  
	if(levelset) {
		sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti, user->_this);
		
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->Levelset, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/lfield%06d_%1d.dat.info",path, ti, user->_this);	if(!rank) unlink(filen);
	}
	
	if(les) {
		Vec Cs;
		
		VecDuplicate(user->P, &Cs);
		DMLocalToGlobal(user->da, user->lCs, INSERT_VALUES, Cs);
		
		sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(Cs, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/cs_%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(NULL);
		VecDestroy(&Cs);
	}
	
	if(rans) {
		sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
		VecView(user->K_Omega, viewer);
		PetscViewerDestroy(&viewer);
		sprintf(filen, "%s/kfield%06d_%1d.dat.info", path, ti, user->_this);	if(!rank) unlink(filen);
		
		PetscBarrier(NULL);
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
	PetscBarrier(NULL);
  
	return 0;
}

PetscErrorCode Divergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

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

  DMDAVecGetArray(fda,user->lUcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
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
  DMDAVecRestoreArray(da, Div, &div);
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
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}

void write_data(UserCtx *user)
{
	DM		da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;

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
		DMDAVecGetArray(da, user->lDensity, &rho);
		DMDAVecGetArray(da, user->lLevelset, &level);
	}
	DMDAVecGetArray(da,user->lAj, &aj);
	DMDAVecGetArray(da, user->lP, &p);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(fda,user->lCsi, &csi);
	DMDAVecGetArray(fda,user->lEta, &eta);
	DMDAVecGetArray(fda,user->lZet, &zet);
	DMDAVecGetArray(fda,user->lUcat, &ucat);
	DMDAVecGetArray(fda,user->lUcont, &ucont);
	DMDAVecGetArray(fda,user->lCent, &cent);
	
	
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
					if(level_in==1){
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
		PetscBarrier(NULL);
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
		PetscBarrier(NULL);
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
				double zz=10.;
				double L = 20., d=1., Beta=0.25, g=fabs(gravity_y), a=0.1; // L = width of the 3D tank
				printf("x=%f, z=%f\n", x, zz );
				double eta0 = a * exp ( -Beta * ( pow(x-L/2, 2) + pow(zz-L/2, 2) ) );
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
					val += eta_mn * exp ( -I * omega_mn * _time ) * cos (n * M_PI / L * x) * cos (m * M_PI / L * zz);
				}
				z_sloshing_exact = val.real();
			}			
			
			
			sprintf(filen, "%s/sloshing.dat", path);
			fp = fopen(filen, "a");
			fprintf(fp, "%e %.7e %.7e %.7e\n", _time, -1.0+z_sloshing,  z_sloshing_exact, fabs(z_sloshing_exact-(-1.0+z_sloshing)));
			fclose(fp);
		}
	} 	

	PetscBarrier(NULL);
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
	PetscBarrier(NULL);
	
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
		DMDAVecRestoreArray(da, user->lDensity, &rho);
		DMDAVecRestoreArray(da, user->lLevelset, &level);
	}
	DMDAVecRestoreArray(da,user->lAj, &aj);
	DMDAVecRestoreArray(da, user->lP, &p);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(fda,user->lCsi, &csi);
	DMDAVecRestoreArray(fda,user->lEta, &eta);
	DMDAVecRestoreArray(fda,user->lZet, &zet);
	DMDAVecRestoreArray(fda,user->lUcat, &ucat);
	DMDAVecRestoreArray(fda,user->lUcont, &ucont);
	DMDAVecRestoreArray(fda,user->lCent, &cent);
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
	PetscBool    DoSCLoop;
	PetscInt      itr_sc;
	
	PetscInt level;
	UserMG usermg;
	PetscBool	flg;
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
	PetscBarrier(NULL);
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
	srand( time(NULL)) ;	// Seokkoo Kang

	PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);
	
	PetscOptionsGetInt(NULL, NULL, "-tio", &tiout, NULL);
	PetscOptionsGetInt(NULL, NULL, "-tiou", &tiout_ufield, NULL);
	PetscOptionsGetInt(NULL, NULL, "-tieu", &tiend_ufield, NULL);

	PetscOptionsGetInt(NULL, NULL, "-imm", &immersed, NULL);
	PetscOptionsGetInt(NULL, NULL, "-inv", &inviscid, NULL);
	PetscOptionsGetInt(NULL, NULL, "-rstart", &tistart, &rstart_flg);
	PetscOptionsGetInt(NULL, NULL, "-imp", &implicit, NULL);
	PetscOptionsGetInt(NULL, NULL, "-imp_MAX_IT", &imp_MAX_IT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-fsi", &movefsi, NULL);
	PetscOptionsGetInt(NULL, NULL, "-rfsi", &rotatefsi, NULL);
	PetscOptionsGetInt(NULL, NULL, "-radi", &radi, NULL);
	PetscOptionsGetInt(NULL, NULL, "-inlet", &inletprofile, NULL);
	PetscOptionsGetInt(NULL, NULL, "-str", &STRONG_COUPLING, NULL);
	PetscOptionsGetInt(NULL, NULL, "-rs_fsi", &rstart_fsi, NULL);
	PetscOptionsGetInt(NULL, NULL, "-cop", &cop, NULL);
	PetscOptionsGetInt(NULL, NULL, "-fish", &fish, NULL);
	PetscOptionsGetInt(NULL, NULL, "-mhv", &MHV, NULL);
	PetscOptionsGetInt(NULL, NULL, "-reg", &regime, NULL);
	PetscOptionsGetInt(NULL, NULL, "-twoD", &TwoD, NULL);
	PetscOptionsGetInt(NULL, NULL, "-thin", &thin, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_z", &dgf_z, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_y", &dgf_y, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_x", &dgf_x, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_az", &dgf_az, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_ay", &dgf_ay, NULL);
	PetscOptionsGetInt(NULL, NULL, "-dgf_ax", &dgf_ax, NULL);
	PetscOptionsGetInt(NULL, NULL, "-body", &NumberOfBodies, NULL);

	PetscOptionsGetInt(NULL, NULL, "-averaging", &averaging, NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0
	PetscOptionsGetInt(NULL, NULL, "-binary", &binary_input, NULL);	// Seokkoo Kang: if 1 binary PLOT3D file, if 0 ascii.
	PetscOptionsGetInt(NULL, NULL, "-xyz", &xyz_input, NULL);			// Seokkoo Kang: if 1 text xyz format, useful for very big (>1GB) Cartesian grid
	PetscOptionsGetInt(NULL, NULL, "-les", &les, NULL);				// Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
	PetscOptionsGetInt(NULL, NULL, "-inlet_buffer_k", &inlet_buffer_k, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wallfunction", &wallfunction, NULL);	// Seokkoo Kang: 1 or 2
	PetscOptionsGetInt(NULL, NULL, "-slipbody", &slipbody, NULL);
	PetscOptionsGetInt(NULL, NULL, "-central", &central, NULL);//central differencing
	PetscOptionsGetInt(NULL, NULL, "-second_order", &second_order, NULL);
	//PetscOptionsGetInt(NULL, NULL, "-initialzero", &initialzero, NULL);	// Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-freesurface", &freesurface, NULL);	// Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-rans", &rans, NULL);			// Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-cross_diffusion", &cross_diffusion, NULL);                     // Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-lowRe", &lowRe, NULL);
	PetscOptionsGetInt(NULL, NULL, "-stension", &surface_tension, NULL);			// Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-delete", &delete_previous_file, NULL);	// Seokkoo Kang: delete previous time step's saved file for saving disk space
	PetscOptionsGetInt(NULL, NULL, "-mixed", &mixed, NULL);			// Seokkoo Kang: mixed model option for LES
	PetscOptionsGetInt(NULL, NULL, "-clark", &clark, NULL);			// Seokkoo Kang: mixed model option for LES
	PetscOptionsGetInt(NULL, NULL, "-vorticity", &vorticity, NULL);			// Seokkoo Kang: vorticity form for viscous terms
	PetscOptionsGetInt(NULL, NULL, "-pseudo", &pseudo_periodic, NULL);	// Seokkoo Kang: pseudo periodic BC in k-direction for genenration of inflow condition

	PetscOptionsGetInt(NULL, NULL, "-levelset", &levelset, NULL);     // Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-levelset_solve2", &levelset_solve2, NULL);     // Seokkoo Kang	
	PetscOptionsGetInt(NULL, NULL, "-fix_level", &fix_level, NULL);     // Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-levelset_it", &levelset_it, NULL);

	PetscOptionsGetInt(NULL, NULL, "-rotdir", &rotdir, NULL);
	PetscOptionsGetInt(NULL, NULL, "-i_periodic", &i_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-j_periodic", &j_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-k_periodic", &k_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-laplacian", &laplacian, NULL);
	PetscOptionsGetInt(NULL, NULL, "-qcrout", &qcr, NULL);
	PetscOptionsGetInt(NULL, NULL, "-ii_periodic", &ii_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-jj_periodic", &jj_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-kk_periodic", &kk_periodic, NULL);	
	periodic = i_periodic+j_periodic+k_periodic+ii_periodic+jj_periodic+kk_periodic;
	PetscOptionsGetInt(NULL, NULL, "-perturb", &initial_perturbation, NULL);	// Seokkoo Kang: give a random perturbation for initial condition
	PetscOptionsGetInt(NULL, NULL, "-skew", &skew, NULL);				// Seokkoo Kang: skew symmetric form of advection term
	PetscOptionsGetInt(NULL, NULL, "-dynamic_freq", &dynamic_freq, NULL);		// Seokkoo Kang: LES dynamic compute frequency 
	if(dynamic_freq<1) dynamic_freq=1;

	PetscOptionsGetInt(NULL, NULL, "-save_inflow", &save_inflow, NULL);		// Seokkoo Kang: save infow BC to files; should be used in conjunction wiht -pseudo 1
	PetscOptionsGetInt(NULL, NULL, "-save_inflow_period", &save_inflow_period, NULL);
	PetscOptionsGetInt(NULL, NULL, "-save_inflow_minus", &save_inflow_minus, NULL);
	PetscOptionsGetInt(NULL, NULL, "-ti_lastsave", &ti_lastsave, NULL);

	PetscOptionsGetInt(NULL, NULL, "-localstep", &localstep, NULL);		// Seokkoo Kang: localstep ( explict + implicit momentum solver )
	PetscOptionsGetInt(NULL, NULL, "-recycle", &inflow_recycle_perioid, NULL);	// Seokkoo Kang, set recycling period of the inflow data
	PetscOptionsGetInt(NULL, NULL, "-save_memory", &save_memory, NULL);	// Seokkoo Kang, save_memory
	PetscOptionsGetInt(NULL, NULL, "-ibm_search", &ibm_search, NULL);
	PetscOptionsGetInt(NULL, NULL, "-ip", &i_proc, NULL);			// Seokkoo Kang: number of processors in i direction
	PetscOptionsGetInt(NULL, NULL, "-jp", &j_proc, NULL);			// Seokkoo Kang: number of processors in j direction
	PetscOptionsGetInt(NULL, NULL, "-kp", &k_proc, NULL);			// Seokkoo Kang: number of processors in k direction
	PetscOptionsGetInt(NULL, NULL, "-poisson", &poisson, NULL); 	// Seokkoo Kang
	PetscOptionsGetInt(NULL, NULL, "-amg_agg", &amg_agg, NULL);
	PetscOptionsGetInt(NULL, NULL, "-testfilter_ik", &testfilter_ik, NULL);
	PetscOptionsGetInt(NULL, NULL, "-testfilter_1d", &testfilter_1d, NULL);
	PetscOptionsGetInt(NULL, NULL, "-poisson_it", &poisson_it, NULL);
	PetscOptionsGetInt(NULL, NULL, "-i_homo_filter", &i_homo_filter, NULL);
	PetscOptionsGetInt(NULL, NULL, "-j_homo_filter", &j_homo_filter, NULL);
	PetscOptionsGetInt(NULL, NULL, "-k_homo_filter", &k_homo_filter, NULL);
	PetscOptionsGetInt(NULL, NULL, "-fix_outlet", &fix_outlet, NULL);
	PetscOptionsGetInt(NULL, NULL, "-fix_inlet", &fix_inlet, NULL);
	
	// add (Toni)
	//for levelset
	PetscOptionsGetInt(NULL, NULL, "-levelset_it", &levelset_it, NULL);
	PetscOptionsGetReal(NULL, NULL, "-levelset_tau", &levelset_tau, NULL);
	PetscOptionsGetReal(NULL, NULL, "-level_in_height", &level_in_height, NULL);
	PetscOptionsGetInt(NULL, NULL, "-level_in", &level_in, NULL);
	PetscOptionsGetInt(NULL, NULL, "-sloshing", &sloshing, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-levelset_weno", &levelset_weno, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-export_FS_elev_center", &export_FS_elev_center, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-export_FS_elev", &export_FS_elev, NULL);	
	
	//for IB FSI
	PetscOptionsGetInt(NULL, NULL, "-forced_motion", &forced_motion, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-fall_cyll_case", &fall_cyll_case, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-fsi_6dof", &fsi_6dof, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-body_mass", &body_mass, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_inertia_x", &body_inertia_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_inertia_y", &body_inertia_y, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_inertia_z", &body_inertia_z, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_rot_x", &body_alpha_rot_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_rot_y", &body_alpha_rot_y, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_rot_z", &body_alpha_rot_z, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_lin_x", &body_alpha_lin_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_lin_y", &body_alpha_lin_y, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_alpha_lin_z", &body_alpha_lin_z, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-body_beta_rot_x", &body_beta_rot_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_beta_rot_y", &body_beta_rot_y, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-body_beta_rot_z", &body_beta_rot_z, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_beta_lin_x", &body_beta_lin_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-body_beta_lin_y", &body_beta_lin_y, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-body_beta_lin_z", &body_beta_lin_z, NULL);
	PetscOptionsGetReal(NULL, NULL, "-angle_x0", &(angle_x0), NULL);
	PetscOptionsGetReal(NULL, NULL, "-angle_y0", &(angle_y0), NULL);
	PetscOptionsGetReal(NULL, NULL, "-angle_z0", &(angle_z0), NULL);

	//Variables for wave_momentum_source
	PetscOptionsGetInt(NULL, NULL, "-wave_momentum_source", &wave_momentum_source, NULL); //For momentum source 1 in x direction and 2 both directions
	PetscOptionsGetInt(NULL, NULL, "-wave_sponge_layer", &wave_sponge_layer, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_zs", &wave_sponge_zs, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_z01", &wave_sponge_z01, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_z02", &wave_sponge_z02, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_xs", &wave_sponge_xs, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_x01", &wave_sponge_x01, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_sponge_x02", &wave_sponge_x02, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_angle_single", &wave_angle_single, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_K_single", &wave_K_single, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_depth", &wave_depth, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_a_single", &wave_a_single, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_wind_reflength", &wave_wind_reflength, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_wind_refvel", &wave_wind_refvel, NULL);
	PetscOptionsGetReal(NULL, NULL, "-wave_wind_yshift", &wave_wind_yshift, NULL);	
	//for air_flow_levelset
	PetscOptionsGetInt(NULL, NULL, "-air_flow_levelset", &air_flow_levelset, NULL);
	PetscOptionsGetInt(NULL, NULL, "-air_flow_levelset_periodic", &air_flow_levelset_periodic, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-wave_average_k", &wave_average_k, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wave_skip", &wave_skip, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wave_ti_start", &wave_ti_start, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-wind_skip", &wind_skip, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wind_start_read", &wind_start_read, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wave_start_read", &wave_start_read, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-wind_recicle", &wind_recicle, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wave_recicle", &wave_recicle, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-wave_k_ave", &wave_k_ave, NULL);
	PetscOptionsGetInt(NULL, NULL, "-wave_i_ave", &wave_i_ave, NULL);	
	PetscOptionsGetInt(NULL, NULL, "-wave_ti_startave", &wave_ti_startave, NULL);
	PetscOptionsGetInt(NULL, NULL, "-freesurface_wallmodel", &freesurface_wallmodel, NULL);
	PetscOptionsGetInt(NULL, NULL, "-viscosity_wallmodel", &viscosity_wallmodel, NULL);
	PetscOptionsGetReal(NULL, NULL, "-channel_height", &channel_height, NULL);
	PetscOptionsGetInt(NULL, NULL, "-floating_turbine_case", &floating_turbine_case, NULL);
	//End variables for wave_momentum_source
	// End (Toni)	
	
	// add begin (xiaolei)
	PetscOptionsGetInt(NULL, NULL, "-forcewidthfixed", &forcewidthfixed, NULL); // xyang 6-3-2011
	PetscOptionsGetReal(NULL, NULL, "-dhi_fixed", &dhi_fixed, NULL); // xyang 6-3-2011
	PetscOptionsGetReal(NULL, NULL, "-dhj_fixed", &dhj_fixed, NULL); // xyang 6-3-2011
	PetscOptionsGetReal(NULL, NULL, "-dhk_fixed", &dhk_fixed, NULL); // xyang 6-3-2011
	if(forcewidthfixed) {
		PetscPrintf(PETSC_COMM_WORLD, "\n Actuator model: the width for force distribution is fixed at dhi=%le dhj=%le dhk=%le \n\n", dhi_fixed, dhj_fixed, dhj_fixed);
	}
	PetscOptionsGetInt(NULL, NULL, "-temperature_rotormodel", &temperature_rotormodel, NULL); // xyang 6-3-2011
	PetscOptionsGetReal(NULL, NULL, "-tmprt_initval", &tmprt_initval, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-u_settling", &u_settling, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-v_settling", &v_settling, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-w_settling", &w_settling, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-fractional_Max_IT", &fractional_Max_IT, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-wallmodel_test", &wallmodel_test, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-dp_wm", &dp_wm, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-alfa_wm", &alfa_wm, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-forceavgperiod_AL", &forceavgperiod_AL, NULL); // xyang 12-7-2010
	PetscOptionsGetInt(NULL, NULL, "-rotor_modeled", &rotor_model, NULL); // xyang 12-7-2010
	PetscOptionsGetInt(NULL, NULL, "-rotate_IBdelta", &rotate_IBdelta, NULL); // xyang 12-7-2010
	PetscOptionsGetInt(NULL, NULL, "-rotate_nacelle", &rotate_nacelle, NULL); // xyang 12-7-2010
	PetscOptionsGetInt(NULL, NULL, "-nacelle_model", &nacelle_model, NULL); // xyang 12-7-2010
	PetscOptionsGetReal(NULL, NULL, "-indf_ax", &indf_ax, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-percent_weno", &percent_weno, NULL); // xyang 12-16-2010
	PetscOptionsGetInt(NULL, NULL, "-imin_wm", &imin_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-imax_wm", &imax_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-jmin_wm", &jmin_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-jmax_wm", &jmax_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-kmin_wm", &kmin_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-kmax_wm", &kmax_wm, NULL); // xyang 1-11-2011
	PetscOptionsGetInt(NULL, NULL, "-imin_wmtmprt", &imin_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-imax_wmtmprt", &imax_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-jmin_wmtmprt", &jmin_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-jmax_wmtmprt", &jmax_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-kmin_wmtmprt", &kmin_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-kmax_wmtmprt", &kmax_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-surface_p_out", &surface_p_out, NULL); // xyang 3-2-2011
	PetscOptionsGetInt(NULL, NULL, "-num_blade", &num_blade, NULL); // xyang 3-2-2011
	PetscOptionsGetInt(NULL, NULL, "-num_foiltype", &num_foiltype, NULL); // xyang 3-18-2011
	PetscOptionsGetReal(NULL, NULL, "-dh1_wm", &dh1_wm, NULL); //xyang 4-11-2011
	PetscOptionsGetReal(NULL, NULL, "-dhratio_wm", &dhratio_wm, NULL); //xyang 4-11-2011
	PetscOptionsGetInt(NULL, NULL, "-IB_delta", &IB_delta, NULL); // xyang 3-18-2011
	PetscOptionsGetInt(NULL, NULL, "-NumIBPerLoc", &NumIBPerLoc, NULL); // xyang 3-18-2011
	PetscOptionsGetInt(NULL, NULL, "-NumNacellePerLoc", &NumNacellePerLoc, NULL); // xyang 3-18-2011
	PetscOptionsGetReal(NULL, NULL, "-reflength_wt", &reflength_wt, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-reflength_nacelle", &reflength_nacelle, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-refvel_wt", &refvel_wt, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-refvel_cfd", &refvel_cfd, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-reflength_IBDelta", &reflength_IBDelta, NULL); // xyang 12-16-2010
	PetscOptionsGetReal(NULL, NULL, "-tipspeedratio", &tipspeedratio, NULL); // xyang 2-12-2012
	PetscOptionsGetReal(NULL, NULL, "-r_nacelle", &r_nacelle, NULL); // xyang 2-12-2012
	PetscOptionsGetReal(NULL, NULL, "-L_nacelle", &L_nacelle, NULL); // xyang 2-12-2012
	PetscOptionsGetReal(NULL, NULL, "-dh_nacelle", &dh_nacelle, NULL); // xyang 2-12-2012
	PetscOptionsGetReal(NULL, NULL, "-loc_refvel", &loc_refvel, NULL); // xyang 2-12-2012
	PetscOptionsGetReal(NULL, NULL, "-X_control", &X_control, NULL); // xyang 2-12-2012
	PetscOptionsGetInt(NULL, NULL, "-IB_wm", &IB_wm, NULL); // xyang 6-3-2011
	PetscOptionsGetInt(NULL, NULL, "-IB_wmtmprt", &IB_wmtmprt, NULL); // xyang 10-22-2012
	PetscOptionsGetInt(NULL, NULL, "-temperature", &temperature, NULL); // xyang 6-3-2011
	PetscOptionsGetInt(NULL, NULL, "-deltafunc", &deltafunc, NULL); // xyang 6-3-2011
	PetscOptionsGetReal(NULL, NULL, "-prt_eps", &prt_eps, NULL);
	PetscOptionsGetReal(NULL, NULL, "-halfwidth_dfunc", &halfwidth_dfunc, NULL); 
	PetscOptionsGetInt(NULL, NULL, "-i_periodicIB", &i_periodicIB, NULL); // xyang 6-3-2011
	PetscOptionsGetInt(NULL, NULL, "-j_periodicIB", &j_periodicIB, NULL); // xyang 6-3-2011
	PetscOptionsGetInt(NULL, NULL, "-k_periodicIB", &k_periodicIB, NULL); // xyang 6-3-2011
	PetscOptionsGetInt(NULL, NULL, "-Force_wm", &Force_wm, NULL); // xyang 10-24-2012
	PetscOptionsGetInt(NULL, NULL, "-Shear_wm", &Shear_wm, NULL); // xyang 10-24-2012
	PetscOptionsGetInt(NULL, NULL, "-Vel_wm", &Vel_wm, NULL); // xyang 10-24-2012
	PetscOptionsGetInt(NULL, NULL, "-add_fluc", &add_fluctuations, NULL); // xyang 11-03-2011
	PetscOptionsGetInt(NULL, NULL, "-add_fluc_tmprt", &add_fluctuations_tmprt, NULL); // xyang 11-03-2011
	PetscOptionsGetInt(NULL, NULL, "-les_prt", &les_prt, NULL); // xyang 11-03-2011
	PetscOptionsGetInt(NULL, NULL, "-infRe", &infRe, NULL);
	PetscOptionsGetInt(NULL, NULL, "-MoveFrame", &MoveFrame, NULL);
	PetscOptionsGetReal(NULL, NULL, "-u_frame", &u_frame, NULL);
	PetscOptionsGetReal(NULL, NULL, "-v_frame", &v_frame, NULL);
	PetscOptionsGetReal(NULL, NULL, "-w_frame", &w_frame, NULL);
	PetscOptionsGetReal(NULL, NULL, "-C_iww", &C_iww, NULL);
	PetscOptionsGetReal(NULL, NULL, "-a_iww", &a_iww, NULL);
	PetscOptionsGetReal(NULL, NULL, "-lamda_iww", &lamda_iww, NULL);
	PetscOptionsGetReal(NULL, NULL, "-xmin", &xmin, NULL);
	PetscOptionsGetReal(NULL, NULL, "-xmax", &xmax, NULL);
	PetscOptionsGetReal(NULL, NULL, "-ymin", &ymin, NULL);
	PetscOptionsGetReal(NULL, NULL, "-ymax", &ymax, NULL);
	PetscOptionsGetReal(NULL, NULL, "-zmin", &zmin, NULL);
	PetscOptionsGetReal(NULL, NULL, "-zmax", &zmax, NULL);
	PetscOptionsGetInt(NULL, NULL, "-ii_periodicWT", &ii_periodicWT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-jj_periodicWT", &jj_periodicWT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-kk_periodicWT", &kk_periodicWT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-Nx_WT", &Nx_WT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-Ny_WT", &Ny_WT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-Nz_WT", &Nz_WT, NULL);
	PetscOptionsGetReal(NULL, NULL, "-Sx_WT", &Sx_WT, NULL);
	PetscOptionsGetReal(NULL, NULL, "-Sy_WT", &Sy_WT, NULL);
	PetscOptionsGetReal(NULL, NULL, "-Sz_WT", &Sz_WT, NULL);
	PetscOptionsGetInt(NULL, NULL, "-New_wallmodel", &New_wallmodel, NULL);
	PetscOptionsGetInt(NULL, NULL, "-turbine", &NumberOfTurbines, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-NumberOfIBDelta", &NumberOfIBDelta, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-NumberOfNacelle", &NumberOfNacelle, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-SpongeLayer", &SpongeLayer, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-SpongeDistance", &SpongeDistance, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-MoveCylinderTest", &MoveCylinderTest, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-inflow_levelset", &inflow_levelset, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-sublevel", &sublevel, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-smoothlevel", &smoothlevel, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-dfunc_wd", &dfunc_wd, NULL);
	PetscOptionsGetInt(NULL, NULL, "-FixTipSpeedRatio", &FixTipSpeedRatio, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-rstart_turbinerotation", &rstart_turbinerotation, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-turbinetorquecontrol", &turbinetorquecontrol, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-turbine_prescrived_motion", &turbine_prescrived_motion, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-turbine_6dof_fsi_motion", &turbine_6dof_fsi_motion, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-turbine_prescrived_motion_heave", &turbine_prescrived_motion_heave, NULL);
	PetscOptionsGetInt(NULL, NULL, "-turbine_prescrived_motion_pitch", &turbine_prescrived_motion_pitch, NULL);
	PetscOptionsGetInt(NULL, NULL, "-fixturbineangvel", &fixturbineangvel, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-fixturbinegeneratortorque", &fixturbinegeneratortorque, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-Shen_AL", &Shen_AL, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-correction3D_CH", &correction3D_CH, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-correction3D_CL", &correction3D_CL, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-smoothforce_AL", &smoothforce_AL, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-radialforce_AL", &radialforce_AL, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-AL_Noslip", &AL_Noslip, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-Shen1_AL", &Shen1_AL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-correction_ALShen1", &correction_ALShen1, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-relax_AL", &relax_AL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-cf_radialforce_AL", &cf_radialforce_AL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-c0_CL", &c0_CL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-c1_CH", &c1_CH, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-c2_CH", &c2_CH, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-c3_CH", &c3_CH, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-a_shen", &a_shen, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-b_shen", &b_shen, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-c_shen", &c_shen, NULL);  // xyang
	PetscOptionsGetInt(NULL, NULL, "-Prandtl_AL", &Prandtl_AL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-refangle_AL", &refangle_AL, NULL);  // xyang
	PetscOptionsGetReal(NULL, NULL, "-dt_bed", &dt_bed, NULL);
	PetscOptionsGetReal(NULL, NULL, "-particle_dens", &particle_dens, NULL);
	PetscOptionsGetReal(NULL, NULL, "-particle_D50", &particle_D50, NULL);
	PetscOptionsGetInt(NULL, NULL, "-smooth_bed", &smooth_bed, NULL);
	PetscOptionsGetInt(NULL, NULL, "-bed_avg", &bed_avg, NULL);
	PetscOptionsGetInt(NULL, NULL, "-particlevel_model", &particlevel_model, NULL);
	PetscOptionsGetInt(NULL, NULL, "-LS_bedConstr", &LS_bedConstr, NULL);
	PetscOptionsGetReal(NULL, NULL, "-Angle_repose", &Angle_repose, NULL);
	Angle_repose	  *= M_PI;
	Angle_repose	  /= 180.;
	PetscOptionsGetReal(NULL, NULL, "-bed_porosity", &bed_porosity, NULL);
	PetscOptionsGetReal(NULL, NULL, "-sb_bed", &sb_bed, NULL);
	PetscOptionsGetInt(NULL, NULL, "-number_smooth", &number_smooth, NULL);
	PetscOptionsGetReal(NULL, NULL, "-deltab", &deltab, NULL);
	PetscOptionsGetReal(NULL, NULL, "-rs_bed", &rs_bed, NULL);
	// FSI
	PetscOptionsGetInt(NULL, NULL, "-powerlawwallmodel", &powerlawwallmodel, NULL);
	PetscOptionsGetInt(NULL, NULL, "-prescribed_rotation", &prescribed_rotation, NULL);
	PetscOptionsGetReal(NULL, NULL, "-scale_velocity", &scale_velocity, NULL);
	// add end (xiaolei)
	
	if(movefsi || rotatefsi) save_memory=0;
	sprintf(path, ".");
	PetscOptionsGetString(NULL, NULL, "-path", path, 256, NULL);		//  Seokkoo Kang: path for saving output; grid.dat, bcs,dat should be put there, but control.dat should exist in the current directory where job is submitted
	
	sprintf(gridfile, "grid.dat");
	PetscOptionsGetString(NULL, NULL, "-grid", gridfile, 256, NULL);	//  Seokkoo Kang: the name of the grid file other than grid.dat if you want

	sprintf(path_inflow, "./inflow");
	PetscOptionsGetString(NULL, NULL, "-path_inflow", path_inflow, 256, NULL);         //  Xiaolei Yang: path for inflow field

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
 
	//PetscOptionsGetReal(NULL, NULL, "-Fr", &Fr, NULL);
	PetscOptionsGetReal(NULL, NULL, "-max_cs", &max_cs, NULL);
	PetscOptionsGetReal(NULL, NULL, "-flux", &inlet_flux, NULL);			// Seokkoo Kang: the amount of inlet flux, if not set mean bulk velocity is set to 1
	PetscOptionsGetReal(NULL, NULL, "-imp_tol", &imp_free_tol, NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(NULL, NULL, "-poisson_tol", &poisson_tol, NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(NULL, NULL, "-les_eps", &les_eps, NULL);		// Seokkoo Kang: small value for preventing very large Cs values in les>1
	PetscOptionsGetReal(NULL, NULL, "-dpdz", &mean_pressure_gradient, &dpdz_set);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
	PetscOptionsGetReal(NULL, NULL, "-roughness", &roughness_size, &rough_set);	// Seokkoo Kang: roughness_size
	PetscOptionsGetReal(NULL, NULL, "-amg_thresh", &amg_thresh, NULL);
	PetscOptionsGetReal(NULL, NULL, "-rho0", &rho_water, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-rho1", &rho_air, NULL);		
	PetscOptionsGetReal(NULL, NULL, "-mu0", &mu_water, NULL);	
	PetscOptionsGetReal(NULL, NULL, "-mu1", &mu_air, NULL);
	PetscOptionsGetReal(NULL, NULL, "-dthick", &dthick, &dthick_set);
	PetscPrintf(PETSC_COMM_WORLD, "\nrho0=%f, rho1=%f, mu0=%f, mu1=%f\n", rho_water, rho_air, mu_water, mu_air);
	PetscOptionsGetReal(NULL, NULL, "-gx", &gravity_x, NULL);
	PetscOptionsGetReal(NULL, NULL, "-gy", &gravity_y, NULL);
	PetscOptionsGetReal(NULL, NULL, "-gz", &gravity_z, NULL);
	PetscOptionsGetReal(NULL, NULL, "-x_r", &(x_r), NULL);
	PetscOptionsGetReal(NULL, NULL, "-y_r", &(y_r), NULL);
	PetscOptionsGetReal(NULL, NULL, "-z_r", &(z_r), NULL);

	PetscOptionsGetReal(NULL, NULL, "-inlet_y", &inlet_y, &inlet_y_flag);
	PetscOptionsGetReal(NULL, NULL, "-outlet_y", &outlet_y, NULL);
	
	PetscOptionsGetReal(NULL, NULL, "-inlet_z", &inlet_z, &inlet_z_flag);
	PetscOptionsGetReal(NULL, NULL, "-outlet_z", &outlet_z, NULL);
	
	PetscOptionsGetReal(NULL, NULL, "-x_c", &(CMx_c), NULL);
	PetscOptionsGetReal(NULL, NULL, "-y_c", &(CMy_c), NULL);
	PetscOptionsGetReal(NULL, NULL, "-z_c", &(CMz_c), NULL);

	PetscOptionsGetReal(NULL, NULL, "-imp_atol", &(imp_atol), NULL);
	PetscOptionsGetReal(NULL, NULL, "-imp_rtol", &(imp_rtol), NULL);
	PetscOptionsGetReal(NULL, NULL, "-imp_stol", &(imp_stol), NULL);
	
	PetscOptionsGetReal(NULL, NULL, "-angvel", &angvel, NULL);

  if (fish) {
    PetscOptionsGetReal(NULL, NULL, "-St_exp", &(St_exp), NULL);
    PetscOptionsGetReal(NULL, NULL, "-wlngth", &(wavelength), NULL);
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
      PetscBarrier(NULL);
      for (ibi=0; ibi<NumberOfBodies; ibi++) {
				ibm_surface_out(&ibm[ibi], 0, ibi);
      }
    } 
		else {
      for (i=0;i<NumberOfBodies;i++) {
				PetscPrintf(PETSC_COMM_WORLD, "Ibm read!\n");
				/*     ibm_read(ibm0); */
				ibm_read_ucd(&ibm[i], i);
				PetscBarrier(NULL);	
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
		PetscOptionsGetReal(NULL, NULL, "-chact_leng", &cl, NULL);
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
		PetscBarrier(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Turbines read!\n");
		if (rotor_model == 1) {
			double reflength = reflength_wt;
			char fname[80];
			sprintf(fname,"acddata000");
			for (i=0;i<NumberOfTurbines;i++) disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
			PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
			Pre_process(&(user[0]), wtm, NumberOfTurbines); 
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
		PetscBarrier(NULL);
		ti = 0;
		if (rstart_flg) ti = tistart;
	}	
	if (IB_delta) {
		PetscReal cl = 1.;
		PetscOptionsGetReal(NULL, NULL, "-chact_leng", &cl, NULL);
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

  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (rstart_flg) {
    ti = tistart; tistart++;
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi]));
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
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
				PetscBarrier(NULL);
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
			DMGlobalToLocalBegin(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
			DMGlobalToLocalEnd(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
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
			DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
			DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
    }

    VecCopy(user[bi].Ucont, user[bi].Ucont_o);
    //VecCopy(user[bi].Ucont, user[bi].Ucont_rm2);	// allocate at init.c
    VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);
    VecCopy(user[bi].Ucat, user[bi].Ucat_o);
    VecCopy(user[bi].P, user[bi].P_o);
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
  }
  PetscBarrier(NULL);
  PetscInt tisteps = 1000000;
  PetscOptionsGetInt(NULL, NULL, "-totalsteps", &tisteps, &flg);
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
				DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
				DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      }
      //VecCopy(user[bi].Ucont_rm1, user[bi].Ucont_rm2);
			if(levelset) VecCopy(user[bi].Levelset, user[bi].Levelset_o);
      VecCopy(user[bi].Ucont_o, user[bi].Ucont_rm1);
      VecCopy(user[bi].Ucont, user[bi].Ucont_o);
      VecCopy(user[bi].P, user[bi].P_o);
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      //seokkoo
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
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
