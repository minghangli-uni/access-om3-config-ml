# Modified parameters

Below only shows the modified parameters. The updated full `MOM_input` can be found in `reorder_configs/MOM6-CICE6_0.25deg_ryf_ml/MOM_input`

```fortran
! === module MOM ===
DIABATIC_FIRST = False 					! default=False hence is deleted in the updated MOM_input
dttherm = 3600 									! test with 4x of dt (baroclinic dynamic ts)
HFREEZE = 10.0 									! for calculation of melt potential

! below BAD_VAL are set back to default values, hence are removed in the updated MOM_input
BAD_VAL_SSH_MAX = 50.0          !   [m] default = 20.0
                                ! The value of SSH above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SSS_MAX = 75.0          !   [PPT] default = 45.0
                                ! The value of SSS above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SST_MAX = 55.0          !   [deg C] default = 45.0
                                ! The value of SST above which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
BAD_VAL_SST_MIN = -3.0          !   [deg C] default = -2.1
                                ! The value of SST below which a bad value message is triggered, if
                                ! CHECK_BAD_SURFACE_VALS is true.
! above BAD_VAL are set back to default values, hence are deleted in the updated MOM_input

! === module MOM_grid_init ===                       
MINIMUM_DEPTH = 0.0 						! the minimum depth is controled via `topog.nc`, default=0.0, 
																! hence is deleted in MOM_input.
CHANNEL_CONFIG = "None" 				! default=None, hence is deleted in MOM_input

! === module MOM_verticalGrid ===
NK = 75 												!nk set to 75 at both 1° and 0.25° to match ACCESS-OM2-01 - use KDS75 z* 
																! - test hybrid adaptive coords later;

! === module MOM_EOS ===
DTFREEZE_DP = -7.75E-08         !   [deg C Pa-1] default = 0.0
                                ! When TFREEZE_FORM=LINEAR, this is the derivative of the freezing potential
                                ! temperature with pressure.
                                ! differs from default 0 and would be redundant if 
                                ! we used TFREEZE_FORM = “TEOS10” - should we move to TEOS10 for EOS as well?
                                
! === module MOM_dynamics_split_RK2 ===
TIDES = False                   !   [Boolean] default = False
                                ! If true, apply tidal momentum forcing.
                                ! remove all TIDERS associated parameters. 
                                ! Hence is removed in the updated MOM_input
                                
! === module MOM_tidal_forcing ===
TIDE_M2 = False                 !   [Boolean] default = False
                                ! If true, apply tidal momentum forcing at the M2 frequency. This is only used
                                ! if TIDES is true.
                                ! Hence is removed in the updated MOM_input
TIDE_SAL_SCALAR_VALUE = 0.094   !   [m m-1]
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
                                ! Smagorinsky nonlinear eddy viscosity.
SMAGORINSKY_AH = True
                                ![Boolean] default = False
                                ! If true, use a biharmonic Smagorinsky nonlinear eddy viscosity.

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

