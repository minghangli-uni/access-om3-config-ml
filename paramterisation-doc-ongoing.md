## MOM parameters - ongoing

> "!" listed in `Parameters` indicates the usage in 1deg configuration, while not used in 025deg configuration. 

| Parameters             | Value  | Reasons for selection                                        |
| ---------------------- | ------ | ------------------------------------------------------------ |
| DIABATIC_FIRST         | True   | ! If true, apply diabatic and thermodynamic processes, including buoyancy forcing and mass gain or loss, before stepping the dynamics forward. |
| USE_REGRIDDING         | True   | use the ALE algorithm (regridding/remapping)                 |
| THICKNESSDIFFUSE       | True   | If true, isopycnal surfaces are diffused with a Laplacian coefficient of KHTH |
| THICKNESSDIFFUSE_FIRST | True   | do thickness diffusion or interface height smoothing before dynamics. |
| DT                     | 900    | The (baroclinic) dynamics time step.                         |
| DT_THERM               | 1800   | The thermodynamic and tracer advection time step             |
| HFREEZE                | 10.0   | !melt potential is computed will be min(HFREEZE, OBLD), where OBLD is the boundary layer depth<br />.src/core/MOM.F90<br />!comment<br />3633         ! Here it is assumed that p=0 is OK, since HFrz ~ 10 to 20m, but under ice-shelves this<br/>3634         ! can be a very bad assumption.  ###To fix this, uncomment the following...<br/>3635         !   pres(i) = p_surface(i) + 0.5*(GV%g_Earth*GV%H_to_RZ)*h(i,j,1) |
| DTBT_RESET_PERIOD      | 0.0    | ! The period between recalculations of DTBT (if DTBT <= 0)<br />If 0, DTBT will be set every dynamics time step. |
| FRAZIL                 | TRUE   | ENABLE_THERMODYNAMICS is TRUE in default.<br />src/core/MOM.F90<br /> 160 character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR" !< A string for specifying the freezing point expression<br/> 161 character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78" !< A string for specifying<br/> 162 !! freezing point expression<br/> 163 character*(10), parameter :: TFREEZE_TEOS10_STRING = "TEOS10" !< A string for specifying the freezing point expression<br/> 164 character*(10), parameter :: TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING !< The default freezing point expression<br />default TFREEZE_DEFAULT is in linear form, hence another parameter in DTFREEZE_DP is set to -7.75E-08. |
| BOUND_SALINITY         | TRUE   | ! comment: limit salinity to positive                        |
| C_P                    | 3992.0 |                                                              |
| CHECK_BAD_SURFACE_VALS | True   |                                                              |
| !BAD_VAL_SSH_MAX       | 50.0   | ! The value of SSH above which a bad value message is triggered, |
| !BAD_VAL_SSS_MAX       | 75.0   | ! The value of SSS above which a bad value message is triggered, |
| !BAD_VAL_SST_MAX       | 55.0   | ! The value of SST above which a bad value message is triggered, |
| !BAD_VAL_SST_MIN       | -3.0   | ! The value of SST below which a bad value message is triggered, |
| !SAVE_INITIAL_CONDS    | True   | ! If true, write the initial conditions to a file given by IC_OUTPUT_FILE. |



## module MOM_domains

| Parameters | Value | Reasons for selection                                        |
| ---------- | ----- | ------------------------------------------------------------ |
| TRIPOLAR_N | True  | ! Use tripolar connectivity at the northern edge of the domain.  With TRIPOLAR_N, NIGLOBAL must be even. |
| NIGLOBAL   | 1440  | ! The total number of thickness grid points in the x-direction in the physical domain. |
| NJGLOBAL   | 1080  | ! The total number of thickness grid points in the y-direction in the physical domain. |



## module MOM_hor_index

## module MOM_grid

## module MOM_fixed_initialization

| Parameters | Value      | Reasons for selection                           |
| ---------- | ---------- | ----------------------------------------------- |
| INPUTDIR   | "./input/" | ! The directory in which input files are found. |



## module MOM_grid_init 

| Parameters    | Value          | Reasons for selection                                        |
| ------------- | -------------- | ------------------------------------------------------------ |
| GRID_CONFIG   | mosaic         | ! mosaic - read the grid from a mosaic (supergrid) file set by GRID_FILE. |
| GRID_FILE     | ocean_hgrid.nc | read horizontal grid data                                    |
| TOPO_CONFIG   | file           | ! file - read bathymetric information from the file specified by (TOPO_FILE). |
| TOPO_FILE     | topog.nc       | ! The file from which the bathymetry is read.                |
| MINIMUM_DEPTH | 0.5            | ! anything shallower than MINIMUM_DEPTH is assumed to be land and all fluxes are masked out. |
| MAXIMUM_DEPTH | 6000.0         | ! The maximum depth of the ocean. Controls where open boundaries are located, what kind of boundary condition to impose, and what data to apply, |

## module MOM_open_boundary

| Parameters     | Value       | Reasons for selection                                        |
| -------------- | ----------- | ------------------------------------------------------------ |
| CHANNEL_CONFIG | global_1deg | !     global_1deg - Sets 16 specific channels appropriate for a 1-degree model, as used in CM2G. |

## module MOM_verticalGrid

| Parameters | Value | Reasons for selection         |
| ---------- | ----- | ----------------------------- |
| NK         | 100   | ! The number of model layers. |

```bash
$ pwd
/g/data/ik11/inputs/access-om2/input_20230515_025deg_topog/mom_025deg
[ml0072@gadi-login-03 mom_025deg]$ ncdump -c ocean_vgrid.nc
netcdf ocean_vgrid {
dimensions:
	nzv = 101 ;
variables:
	double zeta(nzv) ;
		zeta:units = "meters" ;
		zeta:standard_name = "vertical_grid_vertex" ;
		zeta:long_name = "vgrid" ;
		zeta:author = "Kial Stewart" ;

// global attributes:
		:history = " | Updated on Mon Nov  2 16:01:19 AEDT 2020 using https://github.com/COSIMA/make_025deg_topo/tree/f262fbf" ;
data:
}
```



## module MOM_tracer_registry



## module MOM_EOS

| Parameters  | Value     | Reasons for selection                                        |
| ----------- | --------- | ------------------------------------------------------------ |
| DTFREEZE_DP | -7.75E-08 | !   [deg C Pa-1] default = 0.0<br/> ! When TFREEZE_FORM=LINEAR, this is the derivative of the freezing potential temperature with pressure. |

## module MOM_restart



## module MOM_tracer_flow_control

| Parameters           | Value | Reasons for selection                     |
| -------------------- | ----- | ----------------------------------------- |
| USE_IDEAL_AGE_TRACER | True  | use the ideal_age_example tracer package. |

## module ideal_age_example



## module MOM_coord_initialization

| Parameters                 | Value                               | Reasons for selection                                        |
| -------------------------- | ----------------------------------- | ------------------------------------------------------------ |
| REGRIDDING_COORDINATE_MODE | ZSTAR                               | !  ZSTAR, Z* - stretched geopotential z*                     |
| ALE_COORDINATE_CONFIG      | FILE:ocean_vgrid.nc,interfaces=zeta | Determines how to specify the coordinate resolution.         |
| REMAPPING_SCHEME           | PPM_H4                              | vertical remapping for all variable! <br />PPM_H4      (3rd-order accurate) |



## module MOM_state_initialization

| Parameters                         | Value                  | Reasons for selection                                        |
| ---------------------------------- | ---------------------- | ------------------------------------------------------------ |
| INIT_LAYERS_FROM_Z_FILE            | True                   | ! If true, initialize the layer thicknesses, temperatures, and salinities from a<br/>! Z-space file on a latitude-longitude grid. |
| TEMP_SALT_Z_INIT_FILE              | ocean_temp_salt.res.nc | ! The name of the z-space input file used to initialize the layer thicknesses, temperatures and salinities." |
| Z_INIT_FILE_PTEMP_VAR              | temp                   | ! The name of the potential temperature variable in                                     TEMP_SALT_Z_INIT_FILE." |
| TEMP_SALT_INIT_VERTICAL_REMAP_ONLY | True                   | ! If true, initial conditions are on the model horizontal grid. Extrapolation over missing ocean values is done using an ICE-9 procedure with vertical ALE remapping . |
| Z_INIT_REMAP_OLD_ALG               | False                  | uses the preferred remapping algorithm for initialization.   |
| Z_INIT_ALE_REMAPPING               | True                   | remap straight to model coordinate from file.                |



## module MOM_diag_mediator

| Parameters       | Value                               | Reasons for selection |
| ---------------- | ----------------------------------- | --------------------- |
| DIAG_COORD_DEF_Z | FILE:ocean_vgrid.nc,interfaces=zeta |                       |



## module MOM_MEKE

$$
\partial_{\tilde{t}}E = \dot{E}_b + \gamma_\eta \dot{E}_{\eta} + \gamma_v \dot{E}_v - (\lambda + C_d |\mathbf{U}_d| + \gamma^2_b)E + \nabla \cdot \left( (\kappa_E + \gamma_M \kappa_M \tilde{E}_M) \nabla E - \kappa_4 \nabla^3 E \right)\\
$$

where, 

- $\tilde{t}=at$, $a$ is scaled time ($a\ge1$  is used to accelerate towards equilibrium) $a$ (MEKE_DTSCALE)
- $\dot{E}_b + \gamma_\eta \dot{E}_{\eta} + \gamma_v \dot{E}_v$ : sources
  - where, $\dot{E_b}$ is a constant background source of energy to avoid $E->0$
  - $\dot{E}_{\eta}$: "GM" source term, equals the mean PE removed by the Gent-Williams closure, and is included/excluded in the MEKE budget by the efficient paramater $\gamma_\eta\in[0,1]$ (MEKE_GMCOEFF default -1)
  - $\dot{E}_v$: friction source, equals the mean kinetic energy removed by lateral viscous fluxes, and is included/excluded in the MEKE budget by the efficiency parameter $\gamma_v\in[0,1]$ (MEKE_FrCOEFF default -1)
- $(\lambda + C_d |\mathbf{U}_d| + \gamma^2_b)E$: local dissipation
- $\nabla \cdot \left( (\kappa_E + \gamma_M \kappa_M \tilde{E}_M) \nabla E - \kappa_4 \nabla^3 E \right)$: smoothing
  - $E$ is laterally diffused by a diffusivity $\kappa_E + \gamma_M \kappa_M$, where
    - $\kappa_E$ is a constant diffusivity (MEKE_KH)
    - $\gamma_M\kappa_M$ is a self diffusion using the diffusivity calculated 

The MEKE equation is 2D and obtained by depth averaging the 3D eddy energy equation. 

| Parameters                 | Value   | Reasons for selection                                        |
| -------------------------- | ------- | ------------------------------------------------------------ |
| USE_MEKE                   | True    | turns on the MEKE scheme which calculates a sub-grid mesoscale eddy kinetic energy budget. |
| MEKE_GMCOEFF               | 1.0     | Efficiency of conversion of Potential Energy into MEKE       |
| MEKE_GEOMETRIC             | True    | uses the GM coefficient formulation from the GEOMETRIC framework (Marshall et al., 2012). |
| MEKE_EQUILIBRIUM_ALT       | True    | use an alternative formula for computing the (equilibrium)initial value of MEKE.<br /> 791     if (CS%MEKE_equilibrium_alt) then<br/> 792       MEKE%MEKE(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * US%Z_to_L*depth_tot(i,j))**2 / cd2<br />!MEKE_GEOMETRIC_alpha: The nondimensional coefficient governing the efficiency of the GEOMETRIC thickness diffusion, default=0.05<br />!SN: Eady growth rate (1/T)<br />!US%Z_to_L<br />!depth_tot, &    ! The depth of the water column [Z ~> m].<br />!cdrag         !< The bottom drag coefficient for MEKE [nondim]. ($C_d$ )<br />!cd2: drag^2 |
| MEKE_EQUILIBRIUM_RESTORING | True    | ! restore MEKE back to its equilibrium value, which is calculated at each time step. |
| MEKE_RESTORING_TIMESCALE   | 1.0E+07 | !   [s] default = 1.0E+06                                    |
|                            |         |                                                              |
| MEKE_KHTH_FAC              | 0.0     | !parameterizations/lateral/MOM_thickness_diffuse.F90<br /> ! A factor that maps MEKE%Kh to KhTh<br /> 232   if (allocated(MEKE%Kh)) then<br/> 233     if (CS%MEKE_GEOMETRIC) then<br/> 234       !$OMP do<br/> 235       do j=js,je ; do I=is-1,ie<br/> 236         Khth_loc_u(I,j) = Khth_loc_u(I,j) + G%OBCmaskCu(I,j) * CS%MEKE_GEOMETRIC_alpha * &<br/> 237                           0.5*(MEKE%MEKE(i,j)+MEKE%MEKE(i+1,j)) / &<br/> 238                           (VarMix%SN_u(I,j) + CS%MEKE_GEOMETRIC_epsilon)<br/> 239       enddo ; enddo<br/> 240     else<br/> 241       do j=js,je ; do I=is-1,ie<br/> 242         Khth_loc_u(I,j) = Khth_loc_u(I,j) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))<br/> 243       enddo ; enddo<br/> 244     endif<br/> 245   endif<br />! |
| MEKE_KHTR_FAC              | 0.0     | Similar to the above. <br />! A factor that maps MEKE%Kh to KhTr. |
| MEKE_KHMEKE_FAC            | 0.5     | ! A factor that maps MEKE%Kh to Kh for MEKE itself.<br />527       do j=js,je ; do I=is-1,ie<br/> 528         ! Limit Kh to avoid CFL violations.<br/> 529         if (allocated(MEKE%Kh)) &<br/> 530           Kh_here = max(0., CS%MEKE_Kh) + &<br/> 531               CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i+1,j)) |
| MEKE_MIN_LSCALE            | True    | ! If true, use a strict minimum of provided length scales rather than harmonic mean. |
| MEKE_VISCOSITY_COEFF_KU    | 0.0     | ! Can be negative torepresent backscatter from the unresolved eddies. cause instabilities issues. |
| MEKE_ALPHA_DEFORM          | 1.0     |                                                              |
| MEKE_ALPHA_RHINES          | 1.0     |                                                              |
| MEKE_ALPHA_EADY            | 1.0     |                                                              |
| MEKE_ALPHA_FRICT           | 1.0     |                                                              |
| MEKE_ALPHA_GRID            | 1.0     |                                                              |
| MEKE_ADVECTION_FACTOR      | 1.0     |                                                              |

Harmonic mean: 
$$
I_m = \frac{1}{\frac{\alpha_d}{L_d}+\frac{\alpha_f}{L_f}+\frac{\alpha_R}{L_R}+\frac{\alpha_\Delta}{L_\Delta}+\frac{\delta[L_c]}{L_c}}
$$

- $\alpha_d$: MEKE_ALPHA_DEFORM

- $\alpha_f$: MEKE_ALPHA_FRICT

- $\alpha_R$: MEKE_ALPHA_RHINES

- $\alpha_\Delta$: MEKE_ALPHA_GRID

- $L_c$: MEKE_FIXED_MIXING_LENGTH

  

### Diffusivity derived from MEKE

The predicted eddy velocity scale $U_e$, can be combined with a mixing length scale to form a diffusivity. 

The primary use of a MEKE derived diffusivity is for use in thickness diffusion (`mom_thickness_diffuse()`)

optionally in along isopycnal mixing of tracers (`mom_tracer_hor_diff()`)

### Viscosity derived from MEKE

As for $\kappa_M$, the predicted eddy viscosity velocity scale can be used to from a harmonic eddy viscosity, and a bi-harmonic eddy viscosity.

```fortran
K_u ! harmonic eddy viscosity
K_4 ! biharmonic eddy viscosity
```

- $\kappa_u=\gamma_u\sqrt{U^2_eA_\Delta}$: harmonic eddy viscosity
- $\kappa_4=\gamma_4\sqrt{U^2_eA^3_\Delta}$: bi-harmonic eddy viscosity

## MOM_lateral_mixing_coeffs

| Parameters          | Value | Reasons for selection                                        |
| ------------------- | ----- | ------------------------------------------------------------ |
| USE_VARIABLE_MIXING | True  | turns on the MEKE scheme which calculates a sub-grid mesoscale eddy kinetic energy budget. |
| RESOLN_SCALED_KH    | True  | Laplacian lateral viscosity is scaled away when the first baroclinic deformation radius is well resolved. |
| RESOLN_SCALED_KHTH  | True  | KHTH is scaled away when the depth is shallowerthan a reference depth: KHTH = MIN(1,H/H0)**N * KHTH, where H0 is a reference depth, |
| KHTH_SLOPE_CFF      | 0.0   |                                                              |
| USE_STORED_SLOPES   | True  | the isopycnal slopes are calculated once and stored for re-use. This uses more memory but avoids calling the equation of state more times than should be necessary. |
| KH_RES_SCALE_COEF   | 0.7   | ! A coefficient that determines how KhTh is scaled away if RESOLN_SCALED_... is true, as F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER). |
| KH_RES_FN_POWER     | 4     |                                                              |
|                     |       |                                                              |

### Resolution function 

the resolution function is expressed in terms of the ratio of grid-spacing to deformation radius. The square of the reolution parameter is, (in MOM6, $L_d$ is $R_d$)
$$
R^2=\frac{L_d^2}{\Delta^2}=\frac{c_g^2}{f^2\Delta^2+c_g\beta\Delta^2}
$$
where, grid spacing $\Delta$,
$$
\Delta^2=\Delta x^2+\Delta y^2
$$

- $L_d$: the first baroclinic deformation radius,
- $c_g$: first mode internal gravity wave speed,
- $f$: Coriolis parameter.
- $\beta=\frac{\partial f}{\partial y}$: meridional gradient. 

1. When $R^2$ is small, the baroclinic eddy dynamics are not resolved by the grid spacing, so an eddy parameterisation is required.  

2. When $R^2$ is larger, parameterisation is unnecessary and may be counterproductive. 

3. The above two considerations suggest that the eddy diffusivity should be multiplied by a resolution function $r(\Delta,L_d)$, which 

   - decreases from 1 for small $R$,

   - to 0 for large $R$.

   - One candidate resolution function that meets these criteria is,
     $$
     F_4(R^2)=\frac{1}{1+(0.5R^2)^2}
     $$
     

The **resolution function** used in scaling diffusivities is,
$$
r(\Delta,L_d) = \frac{1}{1+(\alpha R)^p}
$$
It can be applied independently to,

- thickness diffusion (module `mom_thickness_diffuse()`)
- tracer diffusion (`mom_tracer_hordiff`)
- lateral viscosity (`mom_hor_visc`)

Parametrers:

- $\alpha$: KH_RES_SCALE_COEF (for thickness and tracer diffusivity) Rd: Deformation radius
- $p$: KH_RES_FN_POWER (for thickness and tracer diffusivity)
- $\alpha$: VISC_RES_SCALE_COEF (for lateral viscosity)
- $p$: VISC_RES_FN_POWER (for lateral viscosity)

**Robert Hallberg 2013 using a resolution function to regulate parameterisations of oceanic mesoscale eddy effects. **

### Visbeck diffusivity

thickness diffusivity: describes how quickly layers of different density within the ocean mix and spread horizontally. 

This module calculates factors used in setting the **thickness diffusivity** similar to a Visbeck et al 1997. The factors are combined in `mom_thickness_diffuse()`,
$$
\kappa_h=\alpha_sL_s^2SN
$$
where, 

- $S$ is the magnitude of the isoneutral slope,
- $N$ is the Brunt-Vaisala frequency.

Parametrers:

- $\alpha_s$: KHTH_SLOPE_CFF
- $L_s$: VISBECK_L_SCALE
- $𝑆_{𝑚𝑎𝑥}$: VISBECK_MAX_SLOPE

### Vertical structure function of KhTh

Thickness diffusivity can be prescribed a vertical distribution with the shape of the equivalent **barotropic** velocity mode. 

The structure function is stored in the control structure for this module (`varmix_cs`) but is calculated using subroutines in `mom_wave_speed`.

thickness diffusivity uses the equivalent barotropic structure as the vertical structure of thickness diffusivity: KHTH_USE_EBT_STRUCT



## MOM_set_visc (related to bbl, such as viscosity and thickness of bbl)

> Calculates various values related to the bottom boundary layer, such as the viscosity and thickness of the BBL (set_viscous_BBL).

| Parameters    | Value  | Reasons for selection                                        |
| ------------- | ------ | ------------------------------------------------------------ |
| CHANNEL_DRAG  | True   | the bottom drag is exerted directly on each layer            |
| PRANDTL_TURB  | 1.0    | applied to shear instability                                 |
| HBBL          | 10.0   | The thickness of a bottom boundary layer with a viscosity increased by KV_EXTRA_BBL, if BOTTOMDRAGLAW is not defined, or the thickness over which near-bottom velocities are averaged for the drag law if BOTTOMDRAGLAW is defined but LINEAR_DRAG is not. |
| BBL_THICK_MIN | 0.1    | ! The minimum bottom boundary layer thickness that can be used with BOTTOMDRAGLAW. This might be Kv/(cdrag*drag_bg_vel) to give Kv as the minimum near-bottom viscosity. |
| KV            | 1.0E-4 | The background kinematic viscosity in the interior.          |

## module MOM_thickness_diffuse

| Parameters                   | Value | Reasons for selection                                        |
| ---------------------------- | ----- | ------------------------------------------------------------ |
| KHTH_USE_FGNV_STREAMFUNCTION | True  | use the streamfunction formulation of Ferrari et al., 2010, which effectively emphasizes graver vertical modes by smoothing in the vertical. |
| FGNV_C_MIN                   | 0.01  | ! A minium wave speed used in the Ferrari et al., 2010, streamfunction  formulation. |
| USE_KH_IN_MEKE               | True  | uses the thickness diffusivity calculated here to diffuse MEKE. |

## module MOM_porous_barriers



## module MOM_continuity

## module MOM_continuity_PPM

| Parameters         | Value | Reasons for selection                                        |
| ------------------ | ----- | ------------------------------------------------------------ |
| ETA_TOLERANCE      | 1e-6  | The tolerance for the differences between the barotropic and baroclinic estimates of the sea surface height due to the fluxes through each face.  The total tolerance for SSH is 4 times this value.  The default is 0.5*NK*ANGSTROM, and this should not be set less than about 10^-15*MAXIMUM_DEPTH. |
| VELOCITY_TOLERANCE | 1e-4  | ! The tolerance for barotropic velocity discrepancies between the barotropic solution and  the sum of the layer thicknesses. |



## module MOM_CoriolisAdv

| Parameters     | Value | Reasons for selection                                        |
| -------------- | ----- | ------------------------------------------------------------ |
| BOUND_CORIOLIS | True  | the Coriolis terms at u-points are bounded by the four estimates of (f+rv)v from the four neighboring v-points, and similarly at v-points.  This option would have no effect on the SADOURNY Coriolis scheme if it were possible to use centered difference thickness fluxes. |

## module MOM_tidal_forcing

| Parameters            | Value | Reasons for selection                                        |
| --------------------- | ----- | ------------------------------------------------------------ |
| TIDE_M2               | True  | apply tidal momentum forcing at the M2 frequency. This is only used  if TIDES is true. |
| TIDE_SAL_SCALAR_VALUE | 0.094 | The constant of proportionality between sea surface height (really it should be bottom pressure) anomalies and bottom geopotential anomalies. This is only used if TIDES and TIDE_USE_SAL_SCALAR are true. |

## module MOM_PressureForce



## module MOM_PressureForce_FV

| Parameters                       | Value | Reasons for selection                                        |
| -------------------------------- | ----- | ------------------------------------------------------------ |
| MASS_WEIGHT_IN_PRESSURE_GRADIENT | True  | ! If true, use mass weighting when interpolating T/S for integrals near the bathymetry in FV pressure gradient calculations. |

## module MOM_hor_visc

| Parameters     | Value                       | Reasons for selection                                        |
| -------------- | --------------------------- | ------------------------------------------------------------ |
| LAPLACIAN      | True                        | Laplacian horizontal viscosity.                              |
| SMAG_BI_CONST  | 0.06                        | The nondimensional biharmonic Smagorinsky constant           |
| SMAGORINSKY_KH | True                        | Smagorinsky nonlinear eddy viscosity.                        |
| SMAG_LAP_CONST | 0.15                        | nondimensional Laplacian Smagorinsky constant, often 0.15    |
| AH_VEL_SCALE   | 0.01(both 1deg and 0.25deg) | The velocity scale which is multiplied by the cube of the grid spacing to calculate the biharmonic viscosity. The final viscosity is the largest of this scaled viscosity, the Smagorinsky and Leith viscosities, and AH. |
| KH_VEL_SCALE   | 0.01 - 1deg / 0.0 - 0.25deg | ! The velocity scale which is multiplied by the grid spacing to calculate the Laplacian viscosity. The final viscosity is the largest of this scaled viscosity, the Smagorinsky and Leith viscosities, and KH. |

## module  MOM_vert_friction 

| Parameters             | Value                  | Reasons for selection                                        |
| ---------------------- | ---------------------- | ------------------------------------------------------------ |
| HMIX_FIXED             | 0.5                    | The prescribed depth over which the near-surface viscosity and diffusivity.are elevated when the bulk mixed layer is not used. |
| MAXVEL                 | 6.0                    | The maximum velocity allowed before the velocity components are truncated. |
| CFL_TRUNCATE_RAMP_TIME | 7200                   | The time over which the CFL truncation value is ramped up at the beginning of the run. |
| U_TRUNC_FILE           | U_velocity_truncations | ! The absolute path to a file into which the accelerations leading to zonal velocity truncations are written.                     ! Undefine this for efficiency if this diagnostic is not                         ! needed." |
| V_TRUNC_FILE           | V_velocity_truncations | ! The absolute path to a file into which the accelerations   leading to meridional velocity truncations are written. Undefine this for efficiency if this diagnostic is not needed." |

## module MOM_barotropic

| Parameters          | Value | Reasons for selection                                        |
| ------------------- | ----- | ------------------------------------------------------------ |
| BOUND_BT_CORRECTION | True  | the corrective pseudo mass-fluxes into the barotropic solver are limited to values that require less than maxCFL_BT_cont to be accommodated. |
| BT_PROJECT_VELOCITY | True  | step the barotropic velocity first and project out the velocity tendency by 1+BEBT when calculating the transport.  The default (false) is to use a predictor continuity step to find the pressure field, and then to do a corrector continuity step using a weighted average of the old and new velocities, with weights of (1-BEBT) and BEBT. |
| BEBT                | 0.2   | BEBT determines whether the barotropic time stepping uses the forward-backward time-stepping scheme or a backward Euler scheme. BEBT is valid in the range from 0 (for a forward-backward treatment of nonrotating gravity waves) to 1 (for a backward Euler treatment). In practice, BEBT must be greater than about  0.05. |
| DTBT                | -0.95 | The barotropic time step, in s. DTBT is only used with the split explicit time stepping. To set the time step automatically based the maximum stable value use 0, or a negative value gives the fraction of the stable value. Setting DTBT to 0 is the same as setting it to -0.98. The value of DTBT that will actually be used is an integer fraction of DT, rounding down. |

## module MOM_mixed_layer_restrat

| Parameters                 | Value                  | Reasons for selection                                        |
| -------------------------- | ---------------------- | ------------------------------------------------------------ |
| MIXEDLAYER_RESTRAT         | True                   | a density-gradient dependent re-stratifying flow is imposed in the mixed layer. Can be used in ALE mode without restriction but in layer mode can only be used if BULKMIXEDLAYER is true. |
| FOX_KEMPER_ML_RESTRAT_COEF | 1.0                    | A nondimensional coefficient that is proportional to the ratio of the deformation radius to the dominant lengthscale of the submesoscale mixed layer instabilities, times the minimum of the ratio of the mesoscale eddy kinetic energy to the large-scale geostrophic kinetic energy or 1 plus the square of the grid spacing over the deformation radius, as detailed by Fox-Kemper et al. (2010)<br />re-stratification in the surface mixed layer due to submesoscale eddies. This parameterization applies an overturning circulation dependent on the horizontal buoyancy gradients within the mixed layer. |
| MLE_FRONT_LENGTH           | 500 -1deg/ 250-0.25deg | is the frontal-length scale used to calculate the upscaling of buoyancy gradients that is otherwise represented by the parameter  FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.<br />**a strong sensitivity of SST biases to the MLE parameterization parameters, with the frontal length the most efficient parameter found to optimize for reducing biases.** |
| MLE_MLD_DECAY_TIME         | 2592000                | The time-scale for a running-mean filter applied to the mixed-layer depth used in the MLE restratification parameterization. When the MLD deepens below the current running-mean the running-mean is instantaneously set to the current MLD. |
|                            |                        |                                                              |

## module MOM_diagnostics



## module MOM_diabatic_driver

| Parameters                 | Value | Reasons for selection                                        |
| -------------------------- | ----- | ------------------------------------------------------------ |
| USE_LEGACY_DIABATIC_DRIVER | False | ! If true, use a legacy version of the diabatic subroutine. This is temporary and is needed to avoid change in answers. |

## module MOM_CVMix_KPP

| Parameters      | Value         | Reasons for selection                                        |
| --------------- | ------------- | ------------------------------------------------------------ |
| USE_KPP         | True          | to calculate diffusivities and non-local transport in the OBL. |
| KPP%            |               |                                                              |
| N_SMOOTH        | 3             | The number of times the 1-1-4-1-1 Laplacian filter is applied on OBL depth. |
| MATCH_TECHNIQUE | MatchGradient | MatchGradient  = sigma*(1-sigma)^2 for NLT; diffusivity profile |
| KPP_IS_ADDITIVE | False         | ! If true, adds KPP diffusivity to diffusivity from other schemes.<br />If false, KPP is the only diffusivity wherever KPP is non-zero. |

## module MOM_CVMix_conv

> Vertical viscosity evaluated by 
>
> parameterizations/vertical/MOM_CVMix_conv.F90

| Parameters           | Value | variable_incode | Reasons for selection                                        |
| -------------------- | ----- | --------------- | ------------------------------------------------------------ |
| USE_CVMix_CONVECTION | True  | CVMix_conv_init | turns on the enhanced mixing due to convection via CVMix. This scheme increases diapycnal diffs./viscs. at statically unstable interfaces. Relevant   parameters are contained in the CVMix_CONVECTION% parameter block. |

```fortran
 60   real    :: prandtl_conv !< Turbulent Prandtl number used in convective instabilities.
 67   call get_param(param_file, mdl, "USE_CVMix_CONVECTION", CVMix_conv_init, default=.false., do_not_log=.true.)
 
 97   call get_param(param_file, mdl, "PRANDTL_CONV", prandtl_conv, &
 98                  "The turbulent Prandtl number applied to convective "//&
 99                  "instabilities (i.e., used to convert KD_CONV into KV_CONV)", &
100                  units="nondim", default=1.0)

102   call get_param(param_file, mdl, 'KD_CONV', CS%kd_conv_const, &
103                  "Diffusivity used in convective regime. Corresponding viscosity "//&
104                  "(KV_CONV) will be set to KD_CONV * PRANDTL_CONV.", &
105                  units='m2/s', default=1.00)

114   ! set kv_conv_const based on kd_conv_const and prandtl_conv
115   CS%kv_conv_const = CS%kd_conv_const * prandtl_conv
```

## module MOM_set_diffusivity

| Parameters       | Value | Reasons for selection                                        |
| ---------------- | ----- | ------------------------------------------------------------ |
| SIMPLE_TKE_TO_KD | True  | uses a simple estimate of Kd/TKE that will work for arbitrary vertical coordinates. If false, calculates Kd/TKE and bounds based on exact energetics/nfor an isopycnal layer-formulation." |

## module MOM_tidal_mixing

| Parameters           | Value             | Reasons for selection                                        |
| -------------------- | ----------------- | ------------------------------------------------------------ |
| BBL_MIXING_AS_MAX    | False             | ! If true, take the maximum of the diffusivity from the BBL mixing and the other diffusivities. Otherwise, diffusivity from the BBL_mixing is simply added. |
| INT_TIDE_DECAY_SCALE | 300.3003003003003 | The decay scale away from the bottom for tidal TKE with the new coding when INT_TIDE_DISSIPATION is used." |

>  source: CVMix/Chapter 6
>
> Mixing arises when mechanical energy dissipates at the small scales in the presence of a nonzero gradient of tracer and/or momentum.
>
> 1. breaking internal waves in ocean interior
> 2. tidal wave interacting with continental shelves (friction bottom drag)
>
> general practice to assume unit Prandtl number. 







## MOM_bkgnd_mixing

> ! Adding static vertical background mixing coefficients
>
> ! parameterizations/vertical/MOM_bkgnd_mixing.F90

| Parameters               | Value   | variable_incode          | Reasons for selection                                        |
| ------------------------ | ------- | ------------------------ | ------------------------------------------------------------ |
| PRANDTL_BKGND            | 5.0     | prandtl_bkgnd            | Turbulent Prandtl number used to convert vertical background diffusivities into viscosities. |
| KD                       | 2.0E-5  | Kd                       | The background diapycnal diffusivity of density in the interior. molecular value, ~1e-7 m2 s-1, may be used. |
| HORIZ_VARYING_BACKGROUND | False   | horiz_varying_background | If true, apply vertically uniform, latitude-dependent background diffusivity, as described in Danabasoglu et al., 2012<br />https://mom6.readthedocs.io/en/main/api/generated/pages/Vertical_Viscosity.html?highlight=rayleigh#channel-drag |
| KD_MIN                   | 2.0E-06 | Kd_min                   | The minimum diapycnal diffusivity.                           |
| KD_MAX                   | 0.1     | Kd_max                   | The maximum permitted increment for the diapycnal diffusivity from TKE-based parameterizations, or a negative value for no limit. |

```fortran
 55   real    :: Kd_min                 !< minimum diapycnal diffusivity [Z2 T-1 ~> m2 s-1]
 56   real    :: Kd                     !< interior diapycnal diffusivity [Z2 T-1 ~> m2 s-1]
```

code snippet related to static vertical background mixing coefficients

```fortran
117   real :: Kv                    ! The interior vertical viscosity [Z2 T-1 ~> m2 s-1] - read to set Prandtl
118                                 ! number unless it is provided as a parameter
119   real :: prandtl_bkgnd_comp    ! Kv/CS%Kd. Gets compared with user-specified prandtl_bkgnd.
250     prandtl_bkgnd_comp = CS%prandtl_bkgnd
251     if (CS%Kd /= 0.0) prandtl_bkgnd_comp = Kv / CS%Kd
253     if ( abs(CS%prandtl_bkgnd - prandtl_bkgnd_comp)>1.e-14) then
254       call MOM_error(FATAL, "bkgnd_mixing_init: The provided KD, KV and PRANDTL_BKGND values "//&
255                             "are incompatible. The following must hold: KD*PRANDTL_BKGND==KV")
256     endif
```

## module MOM_CVMix_shear

| Parameters  | Value | Reasons for selection                                        |
| ----------- | ----- | ------------------------------------------------------------ |
| USE_LMD94   | True  | Large-McWilliams-Doney (JGR 1994) shear mixing parameterization. |
| N_SMOOTH_RI | 1     | If > 0, vertically smooth the Richardson number by applying a 1-2-1 filter N_SMOOTH_RI times. |

$$
\kappa_{kpp\ shear} = \left\{

\begin{aligned}

&\kappa_0\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  Ri<0\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{gravitational instability regime}\\
&\kappa_0[1-(\frac{Ri}{Ri_0})^2]^3\ 0<Ri<Ri_0\ \ \ \ \text{shear instability regime}\\
&0\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  Ri\ge Ri_0\ \ \ \ \ \ \ \ \ \ \ \ \text{stable regime}


\end{aligned}

\right.
$$

The form in the shear instability regime falls most rapidly near $Ri = 0.4 Ri_0$, which aims to parameterize the onset of shear instability. In this neighborhood, rapid changes in Ri can cause gravitational instabilities to develop in the vertical, but these are largely controlled by vertically smoothing Ri profiles with a 1 − 2 − 1 smoother. Unlike Pacanowski and Philander (1981), Large et al. (1994) chose a unit Prandtl number for shear induced mixing; i.e., the shear induced viscosity is the same as the shear induced diffusivity,
$$
Pr_{kpp} = \frac{\nu_{kpp\ shear}}{\kappa_{kpp\ shear}}
$$


## module MOM_CVMix_ddiff

| Parameters      | Value | Reasons for selection                                        |
| --------------- | ----- | ------------------------------------------------------------ |
| USE_CVMIX_DDIFF | True  | If true, turns on double diffusive processes via CVMix. Note that double diffusive processes on viscosity are ignored in CVMix, see http://cvmix.github.io/ for justification. |

>  [Heat diffuses about 100 times faster than salt](https://www.sciencedirect.com/topics/earth-and-planetary-sciences/double-diffusion).
>
> 



## module MOM_diabatic_aux

| Parameters                | Value | Reasons for selection                                        |
| ------------------------- | ----- | ------------------------------------------------------------ |
| PRESSURE_DEPENDENT_FRAZIL | True  | If true, use a pressure dependent freezing temperature when making frazil. The default is false, which will be faster but is inappropriate with ice-shelf cavities. |

## module MOM_opacity

| Parameters   | Value | Reasons for selection                                        |
| ------------ | ----- | ------------------------------------------------------------ |
| PEN_SW_SCALE | 15.0  | The vertical absorption e-folding depth of the penetrating shortwave radiation. |
| PEN_SW_FRAC  | 0.42  | The fraction of the shortwave radiation that penetrates below the surface. |

## module MOM_tracer_advect

| Parameters              | Value  | Reasons for selection                        |
| ----------------------- | ------ | -------------------------------------------- |
| TRACER_ADVECTION_SCHEME | PPM:H3 | Piecewise Parabolic Method (Huyhn 3rd order) |

## module MOM_tracer_hor_diff

| Parameters           | Value | Reasons for selection                                        |
| -------------------- | ----- | ------------------------------------------------------------ |
| KHTR_MIN             | 50.0  | The minimum along-isopycnal tracer diffusivity.              |
| CHECK_DIFFUSIVE_CFL  | True  | ! If true, use enough iterations the diffusion to ensure that the diffusive equivalent of the CFL limit is not violated.  If false, always use the greater of 1 or MAX_TR_DIFFUSION_CFL iteration. |
| MAX_TR_DIFFUSION_CFL | 2.0   | ! If positive, locally limit the along-isopycnal tracer diffusivity to keep the diffusive CFL locally at or below this value.  The number of diffusive iterations is often this value or the next greater integer. |

## module MOM_neutral_diffusion

| Parameters            | Value | Reasons for selection                                        |
| --------------------- | ----- | ------------------------------------------------------------ |
| USE_NEUTRAL_DIFFUSION | True  | enables the neutral diffusion module.                        |
| NDIFF_INTERIOR_ONLY   | True  | only applies neutral diffusion in the ocean interior.That is, the algorithm will exclude the surface and bottomboundary layers. |

## module MOM_hor_bnd_diffusion

| Parameters                        | Value | Reasons for selection                                        |
| --------------------------------- | ----- | ------------------------------------------------------------ |
| USE_HORIZONTAL_BOUNDARY_DIFFUSION | True  | enables the horizonal boundary tracer's diffusion module.    |
| HBD_LINEAR_TRANSITION             | True  | apply a linear transition at the base/top of the boundary. The flux will be fully applied at k=k_min and zero at k=k_max. |

## module MOM_sum_output 

## module ocean_stochastics_init

## module ocean_model_init

| Parameters            | Value | Reasons for selection                                        |
| --------------------- | ----- | ------------------------------------------------------------ |
| RESTART_CONTROL       | 3     | An integer whose bits encode which restart files are written. Add 2 (bit 1) for a time-stamped file, and odd (bit 0) for a non-time-stamped file.  A restart file will be saved at the end of the run segment for any non-negative value. |
| OCEAN_SURFACE_STAGGER | A     | A case-insensitive character string to indicate the staggering of the surface velocity field that is returned to the coupler.  Valid values include 'A','B', or 'C'. |
| RESTORE_SALINITY      | True  | the coupled driver will add a globally-balanced fresh-water flux that drives sea-surface salinity toward specified values. |



## module MOM_surface_forcing_nuopc

| Parameters                     | Value               | Reasons for selection                                        |
| ------------------------------ | ------------------- | ------------------------------------------------------------ |
| LATENT_HEAT_FUSION             | 3.337E+05           | The latent heat of fusion.                                   |
| LATENT_HEAT_VAPORIZATION       | 2.501E+06           | The latent heat of fusion.                                   |
| ADJUST_NET_SRESTORE_TO_ZERO    | True                | If true, adjusts the salinity restoring seen to zero whether restoring is via a salt flux or virtual precip. |
| ADJUST_NET_FRESH_WATER_TO_ZERO | True                | If true, adjusts the net fresh-water forcing seen by the ocean (including restoring) to zero. |
| ENTHALPY_FROM_COUPLER          | True                | ! If True, the heat (enthalpy) associated with mass entering/leaving the ocean is provided via coupler. |
| FLUXCONST                      | 0.11                | ! The constant that relates the restoring surface fluxes to the relative surface anomalies (akin to a piston velocity).  Note the non-MKS units. |
| SALT_RESTORE_FILE              | salt_sfc_restore.nc | to find the surface salinity to use for restoring            |
| SRESTORE_AS_SFLUX              | True                | If true, the restoring of salinity is applied as a salt flux instead of as a freshwater flux. |
| MAX_DELTA_SRESTORE             | 0.5                 | The maximum salinity difference used in restoring terms.     |
| GUST_CONST                     | 0.02                | The background gustiness in the winds.                       |
| FIX_USTAR_GUSTLESS_BUG         | True                | If true correct a bug in the time-averaging of the gustless wind friction velocity |

## module MOM_entrain_diffusive

| Parameters    | Value | Reasons for selection                                        |
| ------------- | ----- | ------------------------------------------------------------ |
| MAX_ENT_IT    | 20    | ! "default = 5 The maximum number of iterations that may be used to calculate the interior diapycnal entrainment." |
| TOLERANCE_ENT | 1e-5  | ! The tolerance with which to solve for entrainment values." |

## module MOM_mixed_layer

| Parameters | Value | Reasons for selection                                        |
| ---------- | ----- | ------------------------------------------------------------ |
| HMIX_MIN   | 2.0   | ! The minimum mixed layer depth if the mixed layer depths determined dynamically." |
