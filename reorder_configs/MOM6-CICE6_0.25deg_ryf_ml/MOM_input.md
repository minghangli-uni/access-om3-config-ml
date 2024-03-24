! === module MOM ===
USE_REGRIDDING = True
                                              ![Boolean] default = False
                                              ! If True, use the ALE algorithm (regridding/remapping). If False, use the
                                              ! layered isopycnal algorithm.
THICKNESSDIFFUSE = True
                                              ![Boolean] default = False
                                              ! If true, isopycnal surfaces are diffused with a Laplacian coefficient of KHTH.
THICKNESSDIFFUSE_FIRST = True
                                              ![Boolean] default = False
                                              ! If true, do thickness diffusion or interface height smoothing before dynamics.
                                              ! This is only used if THICKNESSDIFFUSE or APPLY_INTERFACE_FILTER is true.
DT = 900.0
                                              ![s]
                                              ! The (baroclinic) dynamics time step.  The time-step that is actually used will
                                              ! be an integer fraction of the forcing time-step (DT_FORCING in ocean-only mode
                                              ! or the coupling timestep in coupled mode.)
DT_THERM = 3600.0
                                              ![s] default = 1800.0
                                              ! The thermodynamic and tracer advection time step. Ideally DT_THERM should be
                                              ! an integer multiple of DT and less than the forcing or coupling time-step,
                                              ! unless THERMO_SPANS_COUPLING is true, in which case DT_THERM can be an integer
                                              ! multiple of the coupling timestep.  By default DT_THERM is set to DT.
HFREEZE = 10.0
                                              ![m] default = -1.0
                                              ! If HFREEZE > 0, melt potential will be computed. The actual depth over which
                                              ! melt potential is computed will be min(HFREEZE, OBLD), where OBLD is the
                                              ! boundary layer depth. If HFREEZE <= 0 (default), melt potential will not be
                                              ! computed.
DTBT_RESET_PERIOD = 0.0
                                              ![s] default = 3600.0
                                              ! The period between recalculations of DTBT (if DTBT <= 0). If DTBT_RESET_PERIOD
                                              ! is negative, DTBT is set based only on information available at
                                              ! initialization.  If 0, DTBT will be set every dynamics time step. The default
                                              ! is set by DT_THERM.  This is only used if SPLIT is true.
FRAZIL = True
                                              ![Boolean] default = False
                                              ! If true, water freezes if it gets too cold, and the accumulated heat deficit
                                              ! is returned in the surface state.  FRAZIL is only used if
                                              ! ENABLE_THERMODYNAMICS is true.
BOUND_SALINITY = True
                                              ![Boolean] default = False
                                              ! If true, limit salinity to being positive. (The sea-ice model may ask for more
                                              ! salt than is available and drive the salinity negative otherwise.)

MIN_SALINITY = 0.0
                                			  ! "[PPT] default = 0.01
                               		 	      !  The minimum value of salinity when BOUND_SALINITY=True. The default is 0.01
                                			  !  for backward compatibility but ideally should be 0."

C_P = 3992.0
                                              ![J kg-1 K-1] default = 3991.86795711963
                                              ! The heat capacity of sea water, approximated as a constant. This is only used
                                              ! if ENABLE_THERMODYNAMICS is true. The default value is from the TEOS-10
                                              ! definition of conservative temperature.

CHECK_BAD_SURFACE_VALS = True
                                              ![Boolean] default = False
                                              ! If true, check the surface state for ridiculous values.
SAVE_INITIAL_CONDS = True
                                              ![Boolean] default = False
                                              ! If true, write the initial conditions to a file given by IC_OUTPUT_FILE.
! === module MOM_domains ===
TRIPOLAR_N = True
                                              ![Boolean] default = False
                                              ! Use tripolar connectivity at the northern edge of the domain.  With
                                              ! TRIPOLAR_N, NIGLOBAL must be even.
NIGLOBAL = 1440
                                              !
                                              ! The total number of thickness grid points in the x-direction in the physical
                                              ! domain. With STATIC_MEMORY_ this is set in MOM_memory.h at compile time.
NJGLOBAL = 1080
                                              !
                                              ! The total number of thickness grid points in the y-direction in the physical
                                              ! domain. With STATIC_MEMORY_ this is set in MOM_memory.h at compile time.
! === module MOM_hor_index ===

! === module MOM_grid ===

! === module MOM_fixed_initialization ===
INPUTDIR = "./input/"
                                              !default = "."
                                              ! The directory in which input files are found.
! === module MOM_grid_init ===
GRID_CONFIG = "mosaic"
                                              !
                                              ! A character string that determines the method for defining the horizontal
                                              ! grid.  Current options are:
                                              !     mosaic - read the grid from a mosaic (supergrid)
                                              !              file set by GRID_FILE.
                                              !     cartesian - use a (flat) Cartesian grid.
                                              !     spherical - use a simple spherical grid.
                                              !     mercator - use a Mercator spherical grid.
GRID_FILE = "ocean_hgrid.nc"
                                              !
                                              ! Name of the file from which to read horizontal grid data.
TOPO_CONFIG = "file"
                                              !
                                              ! This specifies how bathymetry is specified:
                                              !     file - read bathymetric information from the file
                                              !       specified by (TOPO_FILE).
                                              !     flat - flat bottom set to MAXIMUM_DEPTH.
                                              !     bowl - an analytically specified bowl-shaped basin
                                              !       ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH.
                                              !     spoon - a similar shape to 'bowl', but with an vertical
                                              !       wall at the southern face.
                                              !     halfpipe - a zonally uniform channel with a half-sine
                                              !       profile in the meridional direction.
                                              !     bbuilder - build topography from list of functions.
                                              !     benchmark - use the benchmark test case topography.
                                              !     Neverworld - use the Neverworld test case topography.
                                              !     DOME - use a slope and channel configuration for the
                                              !       DOME sill-overflow test case.
                                              !     ISOMIP - use a slope and channel configuration for the
                                              !       ISOMIP test case.
                                              !     DOME2D - use a shelf and slope configuration for the
                                              !       DOME2D gravity current/overflow test case.
                                              !     Kelvin - flat but with rotated land mask.
                                              !     seamount - Gaussian bump for spontaneous motion test case.
                                              !     dumbbell - Sloshing channel with reservoirs on both ends.
                                              !     shelfwave - exponential slope for shelfwave test case.
                                              !     Phillips - ACC-like idealized topography used in the Phillips config.
                                              !     dense - Denmark Strait-like dense water formation and overflow.
                                              !     USER - call a user modified routine.
TOPO_FILE = "topog.nc"
                                              !default = "topog.nc"
                                              ! The file from which the bathymetry is read.
MAXIMUM_DEPTH = 6000.0
                                              ![m]
                                              ! The maximum depth of the ocean.
                                              ! Controls where open boundaries are located, what kind of boundary condition to impose, and what data to apply,
                                              ! if any.
! === module MOM_open_boundary ===

!all OBC params irrelevant

! === module MOM_verticalGrid ===
NK = 50
                                              ![nondim]
                                              ! The number of model layers.
! === module MOM_tracer_registry ===

! === module MOM_EOS ===

DTFREEZE_DP = -7.75E-08
                                              ![deg C Pa-1] default = 0.0
                                              ! When TFREEZE_FORM=LINEAR, this is the derivative of the freezing potential
                                              ! temperature with pressure.
! === module MOM_restart ===

! === module MOM_tracer_flow_control ===
USE_IDEAL_AGE_TRACER = True
                                              ![Boolean] default = False
                                              ! If true, use the ideal_age_example tracer package.

! === module ideal_age_example ===

! === module MOM_coord_initialization ===

REGRIDDING_COORDINATE_MODE = "ZSTAR"
                                              !default = "LAYER"
                                              ! Coordinate mode for vertical regridding. Choose among the following
                                              ! possibilities:  LAYER - Isopycnal or stacked shallow water layers
                                              !  ZSTAR, Z* - stretched geopotential z*
                                              !  SIGMA_SHELF_ZSTAR - stretched geopotential z* ignoring shelf
                                              !  SIGMA - terrain following coordinates
                                              !  RHO   - continuous isopycnal
                                              !  HYCOM1 - HyCOM-like hybrid coordinate
                                              !  HYBGEN - Hybrid coordinate from the Hycom hybgen code
                                              !  SLIGHT - stretched coordinates above continuous isopycnal
                                              !  ADAPTIVE - optimize for smooth neutral density surfaces
ALE_COORDINATE_CONFIG = "FILE:ocean_vgrid.nc,interfaces=zeta"
                                              !default = "UNIFORM"
                                              ! Determines how to specify the coordinate resolution. Valid options are:
                                              !  PARAM       - use the vector-parameter ALE_RESOLUTION
                                              !  UNIFORM[:N] - uniformly distributed
                                              !  FILE:string - read from a file. The string specifies
                                              !                the filename and variable name, separated
                                              !                by a comma or space, e.g. FILE:lev.nc,dz
                                              !                or FILE:lev.nc,interfaces=zw
                                              !  WOA09[:N]   - the WOA09 vertical grid (approximately)
                                              !  FNC1:string - FNC1:dz_min,H_total,power,precision
                                              !  HYBRID:string - read from a file. The string specifies
                                              !                the filename and two variable names, separated
                                              !                by a comma or space, for sigma-2 and dz. e.g.
                                              !                HYBRID:vgrid.nc,sigma2,dz
                                              !ALE_RESOLUTION = 2.303499698638916, 2.6903486251831055, 3.1421399116516113, 3.6697616577148438, 4.285917282104492, 5.005424499511719, 5.845563888549805, 6.826459884643555, 7.971549987792969, 9.308074951171875, 10.867660522460938, 12.686931610107422, 14.808158874511719, 17.279945373535156, 20.157821655273438, 23.504684448242188, 27.390975952148438, 31.894271850585938, 37.097900390625, 43.088226318359375, 49.94970703125, 57.757049560546875, 66.56375122070312, 76.386962890625, 87.18865966796875, 98.85760498046875, 111.1953125, 123.914794921875, 136.6578369140625, 149.03271484375, 160.6646728515625, 171.2481689453125, 180.5816650390625, 188.5797119140625, 195.2608642578125, 200.720703125, 205.10205078125, 208.565185546875, 211.2705078125, 213.363525390625, 214.97119140625, 216.198974609375, 217.13232421875, 217.83984375, 218.37451171875, 218.7783203125, 219.08203125, 219.310546875, 219.482421875, 219.6123046875 !   [m]
REMAPPING_SCHEME = "PPM_H4"
                                              !default = "PLM"
                                              ! This sets the reconstruction scheme used for vertical remapping for all
                                              ! variables. It can be one of the following schemes: PCM         (1st-order
                                              ! accurate)
                                              ! PLM         (2nd-order accurate)
                                              ! PLM_HYBGEN  (2nd-order accurate)
                                              ! PPM_H4      (3rd-order accurate)
                                              ! PPM_IH4     (3rd-order accurate)
                                              ! PPM_HYBGEN  (3rd-order accurate)
                                              ! WENO_HYBGEN (3rd-order accurate)
                                              ! PQM_IH4IH3  (4th-order accurate)
                                              ! PQM_IH6IH5  (5th-order accurate)
! === module MOM_state_initialization ===
INIT_LAYERS_FROM_Z_FILE = True
                                              ! "[Boolean] default = False
                                              ! If true, intialize the layer thicknesses, temperatures,
                                              ! and salnities from a Z-space file on a latitude-
                                              ! longitude grid."
TEMP_SALT_Z_INIT_FILE = ocean_temp_salt.res.nc
                                              ! "default = 'temp_salt_z.nc'
                                              ! The name of the z-space input file used to initialize
                                              ! the layer thicknesses, temperatures and salinities."
Z_INIT_FILE_PTEMP_VAR = temp
                                              ! "default = 'ptemp'
                                              ! The name of the potential temperature variable in
                                              ! TEMP_SALT_Z_INIT_FILE."

TEMP_SALT_INIT_VERTICAL_REMAP_ONLY = True
                                              ! "[Boolean] default = False"
Z_INIT_REMAP_OLD_ALG = False
                                              ! "[Boolean] default = True
                                              ! If false, uses the preferred remapping algorithm for initialization. If true,
                                              ! use an older, less robust algorithm for remapping."
Z_INIT_ALE_REMAPPING = True
                                              ! "[Boolean] default = False
                                              ! If True, then remap straight to model coordinate from file."
! === module MOM_diag_mediator ===

NUM_DIAG_COORDS   = 2

DIAG_COORDS       = "z 01 ZSTAR", "rho2 02 RHO"

DIAG_COORD_DEF_01 = "FILE:ocean_vgrid.nc,interfaces=zeta"

DIAG_COORD_DEF_02 = "FILE:diag_rho2.nc,interfaces=rho2"

! === module MOM_MEKE ===

! === module MOM_lateral_mixing_coeffs ===
USE_VARIABLE_MIXING = True
                                              ![Boolean] default = False
                                              ! If true, the variable mixing code will be called.  This allows diagnostics to
                                              ! be created even if the scheme is not used.  If KHTR_SLOPE_CFF>0 or
                                              ! KhTh_Slope_Cff>0, this is set to true regardless of what is in the parameter
                                              ! file.
RESOLN_SCALED_KH = True
                                              ![Boolean] default = False
                                              ! If true, the Laplacian lateral viscosity is scaled away when the first
                                              ! baroclinic deformation radius is well resolved.
RESOLN_SCALED_KHTH = True
                                              ![Boolean] default = False
                                              ! If true, the interface depth diffusivity is scaled away when the first
                                              ! baroclinic deformation radius is well resolved.
KHTH_USE_EBT_STRUCT = True
                                              ![Boolean] default = False
                                              ! If true, uses the equivalent barotropic structure as the vertical structure of
                                              ! thickness diffusivity.
KHTH_SLOPE_CFF = 0.01
                                              ![nondim] default = 0.0
                                              ! The nondimensional coefficient in the Visbeck formula for the interface depth
                                              ! diffusivity
USE_STORED_SLOPES = True
                                              ![Boolean] default = False
                                              ! If true, the isopycnal slopes are calculated once and stored for re-use. This
                                              ! uses more memory but avoids calling the equation of state more times than
                                              ! should be necessary.
KH_RES_SCALE_COEF = 0.7
                                              ![nondim] default = 1.0
                                              ! A coefficient that determines how KhTh is scaled away if RESOLN_SCALED_... is
                                              ! true, as F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).
VISC_RES_SCALE_COEF = 4
                                              ![nondim] default = 0.4
                                              ! A coefficient that determines how Kh is scaled away if RESOLN_SCALED_... is
                                              ! true, as F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER). This
                                              ! function affects lateral viscosity, Kh, and not KhTh.

! === module MOM_set_visc ===
CHANNEL_DRAG = True
                                              ![Boolean] default = False
                                              ! If true, the bottom drag is exerted directly on each layer proportional to the
                                              ! fraction of the bottom it overlies.
HBBL = 10.0
                                              ![m]
                                              ! The thickness of a bottom boundary layer with a viscosity increased by
                                              ! KV_EXTRA_BBL if BOTTOMDRAGLAW is not defined, or the thickness over which
                                              ! near-bottom velocities are averaged for the drag law if BOTTOMDRAGLAW is
                                              ! defined but LINEAR_DRAG is not.
DRAG_BG_VEL = 0.1
                                              ![m s-1] default = 0.0
                                              ! DRAG_BG_VEL is either the assumed bottom velocity (with LINEAR_DRAG) or an
                                              ! unresolved  velocity that is combined with the resolved velocity to estimate
                                              ! the velocity magnitude.  DRAG_BG_VEL is only used when BOTTOMDRAGLAW is
                                              ! defined.
BBL_THICK_MIN = 0.1
                                              ![m] default = 0.0
                                              ! The minimum bottom boundary layer thickness that can be used with
                                              ! BOTTOMDRAGLAW. This might be Kv/(cdrag*drag_bg_vel) to give Kv as the minimum
                                              ! near-bottom viscosity.
KV = 1.0E-04
                                              ![m2 s-1]
                                              ! The background kinematic viscosity in the interior. The molecular value, ~1e-6
                                              ! m2 s-1, may be used.
! === module MOM_thickness_diffuse ===
KHTH_USE_FGNV_STREAMFUNCTION = True
                                              ![Boolean] default = False
                                              ! If true, use the streamfunction formulation of Ferrari et al., 2010, which
                                              ! effectively emphasizes graver vertical modes by smoothing in the vertical.
FGNV_C_MIN = 0.01
                                              ![m s-1] default = 0.0
                                              ! A minium wave speed used in the Ferrari et al., 2010, streamfunction
                                              ! formulation.

KHTH_MAX_CFL = 0.1
                                              ![nondimensional] default = 0.8
                                              ! The maximum value of the local diffusive CFL ratio that is permitted for the
                                              ! thickness diffusivity. 1.0 is the marginally unstable value in a pure layered
                                              ! model, but much smaller numbers (e.g. 0.1) seem to work better for ALE-based
                                              ! models.

! === module MOM_porous_barriers ===

! === module MOM_dynamics_split_RK2 ===

! === module MOM_continuity ===

! === module MOM_continuity_PPM ===
ETA_TOLERANCE = 1.0E-06
                                              ![m] default = 2.5E-09
                                              ! The tolerance for the differences between the barotropic and baroclinic
                                              ! estimates of the sea surface height due to the fluxes through each face.  The
                                              ! total tolerance for SSH is 4 times this value.  The default is
                                              ! 0.5*NK*ANGSTROM, and this should not be set less than about
                                              ! 10^-15*MAXIMUM_DEPTH.
VELOCITY_TOLERANCE = 1.0E-04
                                              ![m s-1] default = 3.0E+08
                                              ! The tolerance for barotropic velocity discrepancies between the barotropic
                                              ! solution and  the sum of the layer thicknesses.
! === module MOM_CoriolisAdv ===
BOUND_CORIOLIS = True
                                              ![Boolean] default = False
                                              ! If true, the Coriolis terms at u-points are bounded by the four estimates of
                                              ! (f+rv)v from the four neighboring v-points, and similarly at v-points.  This
                                              ! option would have no effect on the SADOURNY Coriolis scheme if it were
                                              ! possible to use centered difference thickness fluxes.
! === module MOM_tidal_forcing ===

! === module MOM_PressureForce ===

! === module MOM_PressureForce_FV ===
MASS_WEIGHT_IN_PRESSURE_GRADIENT = True
                                              ![Boolean] default = False
                                              ! If true, use mass weighting when interpolating T/S for integrals near the
                                              ! bathymetry in FV pressure gradient calculations.
! === module MOM_hor_visc ===

LAPLACIAN = True                              
                                              !   [Boolean] default = False
                               			      ! If true, use a Laplacian horizontal viscosity.

SMAGORINSKY_AH = True

​                                              !   [Boolean] default = False
​                                              ! If true, use a biharmonic Smagorinsky nonlinear eddy viscosity.
AH_VEL_SCALE = 0.01
​					                           ! The velocity scale which is multiplied by the cube of the grid spacing to calculate the biharmonic viscosity. 
​					                           ! The final viscosity is the largest of this scaled viscosity, the Smagorinsky and Leith viscosities, and AH.
SMAG_BI_CONST = 0.06
​                                              ! "[nondim] default = 0.0
​                                              ! The nondimensional biharmonic Smagorinsky constant,
​                                              ! typically 0.015 - 0.06."
! === module MOM_vert_friction ===
HMIX_FIXED = 0.5
​                                              ![m]
​                                              ! The prescribed depth over which the near-surface viscosity and diffusivity are
​                                              ! elevated when the bulk mixed layer is not used.
MAXVEL = 6.0
​                                              ![m s-1] default = 3.0E+08
​                                              ! The maximum velocity allowed before the velocity components are truncated.
U_TRUNC_FILE = U_velocity_truncations
​                                              ! "default = ''
​                                              ! The absolute path to a file into which the accelerations
​                                              ! leading to zonal velocity truncations are written.
​                                              ! Undefine this for efficiency if this diagnostic is not
​                                              ! needed."
V_TRUNC_FILE = V_velocity_truncations
​                                              ! "default = ''
​                                              ! The absolute path to a file into which the accelerations
​                                              ! leading to meridional velocity truncations are written.
​                                              ! Undefine this for efficiency if this diagnostic is not
​                                              ! needed."
! === module MOM_barotropic ===
BOUND_BT_CORRECTION = True
​                                              ![Boolean] default = False
​                                              ! If true, the corrective pseudo mass-fluxes into the barotropic solver are
​                                              ! limited to values that require less than maxCFL_BT_cont to be accommodated.
BT_PROJECT_VELOCITY = True
​                                              ![Boolean] default = False
​                                              ! If true, step the barotropic velocity first and project out the velocity
​                                              ! tendency by 1+BEBT when calculating the transport.  The default (false) is to
​                                              ! use a predictor continuity step to find the pressure field, and then to do a
​                                              ! corrector continuity step using a weighted average of the old and new
​                                              ! velocities, with weights of (1-BEBT) and BEBT.
BEBT = 0.2
​                                              ![nondim] default = 0.1
​                                              ! BEBT determines whether the barotropic time stepping uses the forward-backward
​                                              ! time-stepping scheme or a backward Euler scheme. BEBT is valid in the range
​                                              ! from 0 (for a forward-backward treatment of nonrotating gravity waves) to 1
​                                              ! (for a backward Euler treatment). In practice, BEBT must be greater than about
​                                              ! 0.05.
DTBT = -0.95
​                                              ![s or nondim] default = -0.98
​                                              ! The barotropic time step, in s. DTBT is only used with the split explicit time
​                                              ! stepping. To set the time step automatically based the maximum stable value
​                                              ! use 0, or a negative value gives the fraction of the stable value. Setting
​                                              ! DTBT to 0 is the same as setting it to -0.98. The value of DTBT that will
​                                              ! actually be used is an integer fraction of DT, rounding down.
! === module MOM_mixed_layer_restrat ===
MIXEDLAYER_RESTRAT = True
​                                              ![Boolean] default = False
​                                              ! If true, a density-gradient dependent re-stratifying flow is imposed in the
​                                              ! mixed layer. Can be used in ALE mode without restriction but in layer mode can
​                                              ! only be used if BULKMIXEDLAYER is true.
FOX_KEMPER_ML_RESTRAT_COEF = 1.0
​                                              ![nondim] default = 0.0
​                                              ! A nondimensional coefficient that is proportional to the ratio of the
​                                              ! deformation radius to the dominant lengthscale of the submesoscale mixed layer
​                                              ! instabilities, times the minimum of the ratio of the mesoscale eddy kinetic
​                                              ! energy to the large-scale geostrophic kinetic energy or 1 plus the square of
​                                              ! the grid spacing over the deformation radius, as detailed by Fox-Kemper et al.
​                                              ! (2010)
MLE_FRONT_LENGTH = 5000.0
​                                              ![m] default = 0.0
​                                              ! If non-zero, is the frontal-length scale used to calculate the upscaling of
​                                              ! buoyancy gradients that is otherwise represented by the parameter
​                                              ! FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is non-zero, it is recommended
​                                              ! to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.
MLE_MLD_DECAY_TIME = 2592000.0
​                                              ![s] default = 0.0
​                                              ! The time-scale for a running-mean filter applied to the mixed-layer depth used
​                                              ! in the MLE restratification parameterization. When the MLD deepens below the
​                                              ! current running-mean the running-mean is instantaneously set to the current
​                                              ! MLD.
! === module MOM_diagnostics ===

! === module MOM_diabatic_driver ===
USE_LEGACY_DIABATIC_DRIVER = False
                                              ![Boolean] default = True
                                              ! If true, use a legacy version of the diabatic subroutine. This is temporary
                                              ! and is needed to avoid change in answers.

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

! === module MOM_CVMix_KPP ===

! === module MOM_CVMix_conv ===

! === module MOM_set_diffusivity ===
SIMPLE_TKE_TO_KD = True
                                              ! "[Boolean] default = False
                                              ! If true, uses a simple estimate of Kd/TKE that will
                                              ! work for arbitrary vertical coordinates. If false,
                                              ! calculates Kd/TKE and bounds based on exact
                                              ! energetics/nfor an isopycnal layer-formulation."
! === module MOM_tidal_mixing ===
BBL_MIXING_AS_MAX = False
                                              ![Boolean] default = True
                                              ! If true, take the maximum of the diffusivity from the BBL mixing and the other
                                              ! diffusivities. Otherwise, diffusivity from the BBL_mixing is simply added.
! === module MOM_bkgnd_mixing ===
KD = 2.0E-05
                                              ![m2 s-1] default = 0.0
                                              ! The background diapycnal diffusivity of density in the interior. Zero or the
                                              ! molecular value, ~1e-7 m2 s-1, may be used.
KD_MIN = 2.0E-06
                                              ![m2 s-1] default = 2.0E-07
                                              ! The minimum diapycnal diffusivity.
! === module MOM_kappa_shear ===
USE_JACKSON_PARAM = True      

​												!   [Boolean] default = False
​                               			 	! If true, use the Jackson-Hallberg-Legg (JPO 2008) shear mixing
​                                				! parameterization.

MAX_RINO_IT = 25
                                              ! "[nondim] default = 50
                                              ! The maximum number of iterations that may be used to
                                              ! estimate the Richardson number driven mixing."
! === module MOM_CVMix_shear ===

! === module MOM_CVMix_ddiff ===

! === module MOM_diabatic_aux ===
PRESSURE_DEPENDENT_FRAZIL = True
                                              ![Boolean] default = False
                                              ! If true, use a pressure dependent freezing temperature when making frazil. The
                                              ! default is false, which will be faster but is inappropriate with ice-shelf
                                              ! cavities.
! === module MOM_regularize_layers ===

! === module MOM_opacity ===
PEN_SW_SCALE = 15.0
                                              ![m] default = 0.0
                                              ! The vertical absorption e-folding depth of the penetrating shortwave
                                              ! radiation.
PEN_SW_FRAC = 0.42
                                              ![nondim] default = 0.0
                                              ! The fraction of the shortwave radiation that penetrates below the surface.
! === module MOM_tracer_advect ===
TRACER_ADVECTION_SCHEME = "PPM:H3"
                                              !default = "PLM"
                                              ! The horizontal transport scheme for tracers:
                                              !   PLM    - Piecewise Linear Method
                                              !   PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)
                                              !   PPM    - Piecewise Parabolic Method (Colella-Woodward)
! === module MOM_tracer_hor_diff ===
CHECK_DIFFUSIVE_CFL = True
                                              ![Boolean] default = False
                                              ! If true, use enough iterations the diffusion to ensure that the diffusive
                                              ! equivalent of the CFL limit is not violated.  If false, always use the greater
                                              ! of 1 or MAX_TR_DIFFUSION_CFL iteration.
! === module MOM_neutral_diffusion ===

! === module MOM_hor_bnd_diffusion ===

! === module MOM_sum_output ===

! === module ocean_stochastics_init ===

! === module ocean_model_init ===
OCEAN_SURFACE_STAGGER = "A"
                                              !default = "C"
                                              ! A case-insensitive character string to indicate the staggering of the surface
                                              ! velocity field that is returned to the coupler.  Valid values include 'A',
                                              ! 'B', or 'C'.
! === module MOM_surface_forcing_nuopc ===
ADJUST_NET_FRESH_WATER_TO_ZERO = True
                                              ![Boolean] default = False
                                              ! If true, adjusts the net fresh-water forcing seen by the ocean (including
                                              ! restoring) to zero.
FLUXCONST = 0.11
                                              ![m day-1] default = 0.0
                                              ! The constant that relates the restoring surface fluxes to the relative surface
                                              ! anomalies (akin to a piston velocity).  Note the non-MKS units.
SALT_RESTORE_FILE = "salt_sfc_restore.nc"
                                              !default = "salt_restore.nc"
                                              ! A file in which to find the surface salinity to use for restoring.
GUST_CONST = 0.02
                                              ![Pa] default = 0.0
                                              ! The background gustiness in the winds.
! === module MOM_file_parser ===

! === module MOM_shared_initialization ===

! === module MOM_surface_forcing ===

! === module MOM_debugging ===

! === module MOM_ALE ===

! === module MOM_regridding ===

! === module MOM_ocean_model_nuopc ===

! === module MOM_set_viscosity ===

! === module MOM_entrain_diffusive ===
! === module MOM_mixed_layer ===
HMIX_MIN = 2.0
                                              ! "[m] default = 0.0
                                              ! The minimum mixed layer depth if the mixed layer depth
                                              ! is determined dynamically."