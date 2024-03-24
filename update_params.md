# Modified parameters

Below only shows the modified parameters. The updated full `MOM_input` can be found in `reorder_configs/MOM6-CICE6_0.25deg_ryf_ml/MOM_input`

```fortran
! === module MOM ===
DIABATIC_FIRST = False 					
                                ! default=False hence is deleted in the updated MOM_input
dttherm = 3600 									
                                ! test with 4x of dt (baroclinic dynamic ts)
HFREEZE = 10.0 									
                                ! for calculation of melt potential

! below BAD_VAL are set back to default values, hence are removed in the updated MOM_input
BAD_VAL_SSH_MAX = 50.0          
                                !   [m] default = 20.0
                                ! The value of SSH above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SSS_MAX = 75.0          
                                !   [PPT] default = 45.0
                                ! The value of SSS above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SST_MAX = 55.0          
                                !   [deg C] default = 45.0
                                ! The value of SST above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SST_MIN = -3.0          
                                !   [deg C] default = -2.1
                                ! The value of SST below which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
! above BAD_VAL are set back to default values, hence are deleted in the updated MOM_input

MIN_SALINITY = 0.0
                                ! "[PPT] default = 0.01
                                !  The minimum value of salinity when BOUND_SALINITY=True. The default is 0.01
                                !  for backward compatibility but ideally should be 0."

! === module MOM_grid_init ===                       
MINIMUM_DEPTH = 0.0 						
                                ! the minimum depth is controled via `topog.nc`, default=0.0, 
																! hence is deleted in MOM_input.
CHANNEL_CONFIG = "None" 				
                                ! default=None, hence is deleted in MOM_input
! === module MOM_verticalGrid ===
NK = 50 												
                                !nk set to 75 at both 1° and 0.25° to match ACCESS-OM2-01 - use KDS75 z* 
																! - test hybrid adaptive coords later;
! === module MOM_EOS ===
DTFREEZE_DP = -7.75E-08         
                                !   [deg C Pa-1] default = 0.0
                                ! When TFREEZE_FORM=LINEAR, this is the derivative of the freezing potential
                                ! temperature with pressure.
                                ! differs from default 0 and would be redundant if 
                                ! we used TFREEZE_FORM = “TEOS10” - should we move to TEOS10 for EOS as well?         
! === module MOM_dynamics_split_RK2 ===
TIDES = False                   
                                !   [Boolean] default = False
                                ! If true, apply tidal momentum forcing.
                                ! remove all TIDERS associated parameters. 
                                ! Hence is removed in the updated MOM_input             
! === module MOM_tidal_forcing ===
TIDE_M2 = False                 
                                !   [Boolean] default = False
                                ! If true, apply tidal momentum forcing at the M2 frequency. This is only used
                                ! if TIDES is true.
                                ! Hence is removed in the updated MOM_input
TIDE_SAL_SCALAR_VALUE = 0.094   
                                !   [m m-1]
                                ! The constant of proportionality between sea surface height (really it should
                                ! be bottom pressure) anomalies and bottom geopotential anomalies. 
                                ! This is only used if TIDES and TIDE_USE_SAL_SCALAR are true.
                                ! Hence is removed in the updated MOM_input
! === module MOM_diag_mediator ===
NUM_DIAG_COORDS 	= 2
DIAG_COORDS 			= "z 01 ZSTAR", "rho2 02 RHO"
DIAG_COORD_DEF_01 = "FILE:ocean_vgrid.nc,interfaces=zeta"
DIAG_COORD_DEF_02 = "FILE:diag_rho2.nc,interfaces=rho2"


! === module MOM_MEKE ===
! not use MEKE, test parameters.
USE_MEKE = False
                                ![Boolean] default = False
                                ! If true, turns on the MEKE scheme which calculates a sub-grid mesoscale eddy
                                ! kinetic energy budget. 
																! Hence is removed in the updated MOM_input

! === module MOM_lateral_mixing_coeffs ===
DEPTH_SCALED_KHTH = False
                                ! "[Boolean] default = False
                                !  If true,  KHTH is scaled away when the depth is shallower
                                !  than a reference depth: KHTH = MIN(1,H/H0)**N * KHTH,
                                !  where H0 is a reference depth, controlled via
                                !  DEPTH_SCALED_KHTH_H0, and theexponent (N) is
                                !  controlled via DEPTH_SCALED_KHTH_EXP."
                                ! default is False, Hence is removed in the updated MOM_input.
KH_RES_SCALE_COEF = 0.7
                                ![nondim] default = 1.0
                                ! A coefficient that determines how KhTh is scaled away if RESOLN_SCALED_... is
                                ! true, as F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).
VISC_RES_SCALE_COEF = 4
                                ![nondim] default = 0.4
                                ! A coefficient that determines how Kh is scaled away if RESOLN_SCALED_... is
                                ! true, as F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER). This
                                ! function affects lateral viscosity, Kh, and not KhTh.

! === module MOM_thickness_diffuse ===
USE_KH_IN_MEKE = False
                                ![Boolean] default = False
                                ! If true, uses the thickness diffusivity calculated here to diffuse MEKE.
																! Hence is removed in the updated MOM_input
KHTH_MAX_CFL = 0.1
                                ![nondimensional] default = 0.8
                                ! The maximum value of the local diffusive CFL ratio that is permitted for the
                                ! thickness diffusivity. 1.0 is the marginally unstable value in a pure layered
                                ! model, but much smaller numbers (e.g. 0.1) seem to work better for ALE-based
                                ! models.

! === module MOM_energetic_PBL === 
! (parameters of EPBL are copied from GFDL-OM4-025 and mom6-panan)
ML_OMEGA_FRAC = 0.001
                                ![nondim] default = 0.0
                                ! When setting the decay scale for turbulence, 
                                ! use this fraction of the absolute
                                ! rotation rate blended with the local value of f, as sqrt((1-of)*f^2 +
                                ! of*4*omega^2).
TKE_DECAY = 0.01
                                ![nondim] default = 2.5
                                ! TKE_DECAY relates the vertical rate of decay of the TKE available for
                                ! mechanical entrainment to the natural Ekman depth.
EPBL_MSTAR_SCHEME = "OM4"
                                !default = "CONSTANT"
                                ! EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are:
                                !    CONSTANT   - Use a fixed mstar given by MSTAR
                                !    OM4        - Use L_Ekman/L_Obukhov in the sabilizing limit, as in OM4
                                !    REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.
MSTAR_CAP = 10.0
                                ![nondim] default = -1.0
                                ! If this value is positive, it sets the maximum value of mstar allowed in ePBL.
                                ! (This is not used if EPBL_MSTAR_SCHEME = CONSTANT).
MSTAR2_COEF1 = 0.29
                                ![nondim] default = 0.3
                                ! Coefficient in computing mstar when rotation and stabilizing effects are both
                                ! important (used if EPBL_MSTAR_SCHEME = OM4).
MSTAR2_COEF2 = 0.152
                                ![nondim] default = 0.085
                                ! Coefficient in computing mstar when only rotation limits the total mixing
                                ! (used if EPBL_MSTAR_SCHEME = OM4)
NSTAR = 0.06
                                ![nondim] default = 0.2
                                ! The portion of the buoyant potential energy imparted by surface fluxes that is
                                ! available to drive entrainment at the base of mixed layer when that energy is
                                ! positive.
MSTAR_CONV_ADJ = 0.667
                                ![nondim] default = 0.0
                                ! Coefficient used for reducing mstar during convection due to reduction of
                                ! stable density gradient.
USE_MLD_ITERATION = True
                                ![Boolean] default = False
                                ! A logical that specifies whether or not to use the distance to the bottom of
                                ! the actively turbulent boundary layer to help set the EPBL length scale.
EPBL_TRANSITION_SCALE = 0.01
                                ![nondim] default = 0.1
                                ! A scale for the mixing length in the transition layer at the edge of the
                                ! boundary layer as a fraction of the boundary layer thickness.
MIX_LEN_EXPONENT = 1.0
                                ![nondim] default = 2.0
                                ! The exponent applied to the ratio of the distance to the MLD and the MLD depth
                                ! which determines the shape of the mixing length. This is only used if
                                ! USE_MLD_ITERATION is True.
USE_LA_LI2016 = True
                                ![nondim] default = False
                                ! A logical to use the Li et al. 2016 (submitted) formula to determine the
                                ! Langmuir number.
EPBL_LANGMUIR_SCHEME = "ADDITIVE"
                                !default = "NONE"
                                ! EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence.
                                ! Valid values are:
                                !    NONE     - Do not do any extra mixing due to Langmuir turbulence
                                !    RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir
                                !      turbulence
                                !    ADDITIVE - Add a Langmuir turblence contribution to mstar to other
                                !      contributions
LT_ENHANCE_COEF = 0.044
                                ![nondim] default = 0.447
                                ! Coefficient for Langmuir enhancement of mstar
LT_ENHANCE_EXP = -1.5
                                ![nondim] default = -1.33
                                ! Exponent for Langmuir enhancementt of mstar
LT_MOD_LAC1 = 0.0
                                ![nondim] default = -0.87
                                ! Coefficient for modification of Langmuir number due to MLD approaching Ekman
                                ! depth.
LT_MOD_LAC4 = 0.0
                                ![nondim] default = 0.95
                                ! Coefficient for modification of Langmuir number due to ratio of Ekman to
                                ! stable Obukhov depth.
LT_MOD_LAC5 = 0.22
                                ![nondim] default = 0.95
                                ! Coefficient for modification of Langmuir number due to ratio of Ekman to
                                ! unstable Obukhov depth.

! === module MOM_hor_visc ===
SMAGORINSKY_KH = False
                                ! Smagorinsky nonlinear eddy viscosity. default is false, so it is removed in the MOM_input
SMAGORINSKY_AH = True
                                ![Boolean] default = False
                                ! If true, use a biharmonic Smagorinsky nonlinear eddy viscosity.

! === module MOM_neutral_diffusion ===
USE_NEUTRAL_DIFFUSION = False
                                              ![Boolean] default = False
                                              ! If true, enables the neutral diffusion module.
                                              ! default is false, so it is removed in the MOM_input
! === module MOM_tidal_mixing ===
INT_TIDE_DECAY_SCALE = 300.3003003003003
                                              ! "[m] default = 500.0
                                              ! The decay scale away from the bottom for tidal TKE with
                                              ! the new coding when INT_TIDE_DISSIPATION is used."
                                              ! LEE_WAVE_DISSIPATION=False and USE_CVMix_TIDAL=False, this parameter has no meaning, removed


! === module MOM_bkgnd_mixing ===
HORIZ_VARYING_BACKGROUND = False
                                              ![Boolean] default = False
                                              ! If true, apply vertically uniform, latitude-dependent background diffusivity,
                                              ! as described in Danabasoglu et al., 2012
                                              ! default is false, so it is removed in the MOM_input

PRANDTL_BKGND = 1.0
                                              ![nondim] default = 1.0
                                              ! Turbulent Prandtl number used to convert vertical background diffusivities
                                              ! into viscosities.
                                              ! default is 1.0, so it is removed in the MOM_input
KD_MAX = -1.0
                                              ![m2 s-1] default = -1.0
                                              ! The maximum permitted increment for the diapycnal diffusivity from TKE-based
                                              ! parameterizations, or a negative value for no limit.
                                              ! default is -1.0, so it is removed in the MOM_input
! === module ocean_model_init ===
RESTART_CONTROL = 1
                                              !default = 1
                                              ! An integer whose bits encode which restart files are written. Add 2 (bit 1)
                                              ! for a time-stamped file, and odd (bit 0) for a non-time-stamped file.  A
                                              ! restart file will be saved at the end of the run segment for any non-negative
                                              ! value.
                                              ! default is 1, so it is removed in the MOM_input
RESTORE_SALINITY = False
                                              ![Boolean] default = False
                                              ! If true, the coupled driver will add a globally-balanced fresh-water flux that
                                              ! drives sea-surface salinity toward specified values.
                                              ! CESM default False, removed.                                              
! === module MOM_surface_forcing_nuopc ===
LATENT_HEAT_FUSION = 3.337E+05
                                              ![J/kg] default = 3.34E+05
                                              ! The latent heat of fusion.
                                              ! leave it to default, removed 
LATENT_HEAT_VAPORIZATION = 2.501E+06
                                              ![J/kg] default = 2.5E+06
                                              ! The latent heat of fusion.
                                              ! leave it to default, removed
ENTHALPY_FROM_COUPLER = False
                                              ![Boolean] default = False
                                              ! If True, the heat (enthalpy) associated with mass entering/leaving the ocean
                                              ! is provided via coupler.
                                              ! default is false, so it is removed in the MOM_input
SRESTORE_AS_SFLUX = False
                                              ![Boolean] default = False
                                              ! If true, the restoring of salinity is applied as a salt flux instead of as a
                                              ! freshwater flux.
                                              ! default is false, so it is removed in the MOM_input
MAX_DELTA_SRESTORE = 999.0
                                              ![PSU or g kg-1] default = 999.0
                                              ! The maximum salinity difference used in restoring terms.     
                                              ! default is 999, so it is removed in the MOM_input                                     
! === module MOM_kappa_shear ===
! Parameterization of shear-driven turbulence following Jackson, Hallberg and Legg, JPO 2008
USE_JACKSON_PARAM = True      

​												                !   [Boolean] default = False
​                               			 	! If true, use the Jackson-Hallberg-Legg (JPO 2008) shear mixing
​                                				! parameterization.
! === module MOM_CVMix_shear ===
USE_LMD94 = False
                                              ![Boolean] default = False
                                              ! If true, use the Large-McWilliams-Doney (JGR 1994) shear mixing
                                              ! parameterization.
                                              ! default is false, so it is removed in the MOM_input
N_SMOOTH_RI = 0
                                              !default = 0
                                              ! If > 0, vertically smooth the Richardson number by applying a 1-2-1 filter
                                              ! N_SMOOTH_RI times.
                                              ! Parameterization of mixing due to double diffusion processes via CVMix
                                              ! default is 0, so it is removed in the MOM_input
! === module MOM_CVMix_ddiff ===
USE_CVMIX_DDIFF = False
                                              ![Boolean] default = False
                                              ! If true, turns on double diffusive processes via CVMix. Note that double
                                              ! diffusive processes on viscosity are ignored in CVMix, see
                                              ! http://cvmix.github.io/ for justification.
                                              ! default is false, so it is removed in the MOM_input

! === module MOM_tracer_hor_diff ===
MAX_TR_DIFFUSION_CFL = -1.0
                                              ![nondim] default = -1.0
                                              ! If positive, locally limit the along-isopycnal tracer diffusivity to keep the
                                              ! diffusive CFL locally at or below this value.  The number of diffusive
                                              ! iterations is often this value or the next greater integer.
                                              ! default is -1.0, so it is removed in the MOM_input
KHTR_MIN = 0.0
                                              ![m2 s-1] default = 0.0
                                              ! The minimum along-isopycnal tracer diffusivity.
                                              ! never used in the source code (checked), removed.
! === module MOM_hor_bnd_diffusion ===
USE_HORIZONTAL_BOUNDARY_DIFFUSION = False
                                              ![Boolean] default = False
                                              ! If true, enables the horizonal boundary tracer's diffusion module.
                                              ! default is false, so it is removed in the MOM_input
HBD_LINEAR_TRANSITION = False
                                              ![Boolean] default = False
                                              ! If True, apply a linear transition at the base/top of the boundary.
                                              ! The flux will be fully applied at k=k_min and zero at k=k_max.   
                                              ! default is false, so it is removed in the MOM_input                                           

! === module MOM_entrain_diffusive ===
MAX_ENT_IT = 5
                                ! "default = 5
                                ! The maximum number of iterations that may be used to
                                ! calculate the interior diapycnal entrainment."
                                ! default value hence is removed in the updated MOM_input
TOLERANCE_ENT = 2.683281572999748E-05
                                ! "[m] default = 2.683281572999748E-05
                                ! The tolerance with which to solve for entrainment values."
                                ! default value hence is removed in the updated MOM_input

```

