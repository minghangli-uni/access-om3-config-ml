# access-om3-config-ml

This branch (`EPBL_MEKE_FALSE`) implements one test of MOM6 parameterisation. 

Most of the updates are sourced from discussions in [namelist-discussion](https://forum.access-hive.org.au/t/namelist-configuration-discussion-meeting/1917/9?u=minghangli). Some significant updates are highlighted below,

- surface boundary layer parameterisation: 
  - implements an energetic constrained parameterisation of the surface boundary layer (EPBL), providing vertically diffusivity and viscosity, and the depth of active mixing (BL thickness).
  - replaces the KPP parameterisation used in CVmix project. 

- Removal of mesoscale eddy mixing parametersiations (e.g., `USE_MEKE=False`)
- Removal of `BAD_VAL_XXX_MAX` nor `BAD_VAL_XXX_MIN`
- Set `NK=75` to match `ACCESS-OM2-01` - use `KDS75 z*`
- No parameters associated with `TIDES`
- `NUM_DIAG_COORDS=2` includes `z_star` and `rho_2`

---

For a more detailed breakdown of only the changed parameters, please refer to the `only_changed_params.md` file.

---

A complete input file with updated parameters can be found in `reorder_configs/MOM6-CICE6_0.25deg_ryf_ml/MOM_input`

