[comment]: <> (this is a markdown document and the special characters are for formatting)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="man/logo/snl-pstess_logo_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="man/logo/snl-pstess_logo.png">
  <img src="man/logo/snl-pstess_logo.png" alt="PSTess logo" width=300px margin="auto" />
</picture>

# Master branch changelog
Welcome to the Power and Energy Storage Systems Toolbox. This is the
**PSTess** CHANGELOG. Please refer to the table of contents below to
help you navigate.

## Table of contents
<a id="toc"></a>
- [Initial release, Version 1.0](#initial_release)
    - [Bug fixes](#bug_fixes)
    - [Miscellaneous](#miscellaneous)
    - [File modifications](#file_mods)

## Initial release, Version 1.0
<a id="initial_release"></a>

This section describes the changes, on a per-file basis, made to
existing **PST** functions and scripts in the creation of **PSTess**.
Please see the README for context and additional information about new
features added to the toolbox.

Note: The **PST v3.0** source code on which **PSTess** is based is
available in the main `pstess` directory under the folder titled
`archive_pstV3p0`.

### Bug fixes
<a id="bug_fixes"></a>

* **exc_st3**: the angle variable `theta` was mistakenly indexed by
  the machine number `n` rather than the bus number `n_bus`; this has
  been corrected
* **inv_lf**: the angle variable `gamma` was mistakenly reported in
  radians rather than degrees; this has been corrected for consistency
  with the rest of the HVDC code, including **rec_lf**
* **tcsc**: this model is evidently based on **svc** because in the
  dynamics calculation the state derivative was accidentally called
  `dB_cv` instead of `dB_tcsc`; this `dB_cv` was subsequently never
  used, so this has been corrected

### Miscellaneous
<a id="miscellaneous"></a>

* Imaginary numbers: the original `jay = sqrt(-1)` convention has been
  eliminated in favor of `1j`
* Dummy variables: functions no longer return unnecessary dummy
  variables (the most common was `f`)
* Unused function inputs: it appears that many functions were based on
  a common template and had unused inputs as a result; these have been
  removed to avoid confusion

### File modifications
<a id="file_mods"></a>

* **calc:** formatting, removal of `jay`, edited the header to remove reference to unused inputs
* **cdps:** rewritten to accommodate batch processing and **get_path**
* **chq_lim:** formatting, removal of dummy variable `f`, reorganized globals, changed inputs and outputs to avoid unnecessary variables
* **contents:** updated the table of contents explaining the purpose of each file
* **contents_global:** new file listing all global variables and field names
* **dbcage:** formatting
* **dc_cont:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements
* **dc_cur:** formatting, removal of `jay`, reorganized globals, removal of unused input `k`
* **dc_indx:** formatting, removal of dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **dc_lf:** formatting, removal of `jay`, reorganized globals, revised error statements
* **dc_line:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused inputs `k` and `bus`
* **dc_load:** formatting, removal of `jay`, typo correction `Yiii` to `Yii`, reorganized globals, removed unused input `k`
* **dc_sim:** formatting, function call changes due to removal of dummy variable `f`, reorganized globals, added `dcrd_sig` and `dcid_sig` as inputs to avoid unnecessary globals
* **dc_vidc:** formatting, removal of dummy variable `f`, reorganized globals
* **dci_sud:** formatting, corrections to deprecated `get()` syntax, revised error and warning statements
* **dcr_sud:** formatting, revised error and warning statements
* **deepbar:** formatting
* **dessat:** formatting, changed the name of the output variable to `gsat` to avoid conflict with the global struct `g`
* **dpwf:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **dpwf_indx:** formatting, removal of dummy variable `f`, reorganized globals, revised error statements
* **ess:** new energy storage device model with reactive modulation capability
* **ess_indx:** new index function for energy storage models
* **ess_sat:** new auxiliary function for implementing saturation characteristics for constant-current inverter-based resources (includes low-voltage power logic)
* **ess_sud:** new user-defined energy storage control model
* **exc_dc12:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **exc_indx:** formatting, removal of dummy variable `f`, reorganized globals
* **exc_st3:** formatting, removal of `jay` and dummy variable `f`, patched `n` vs. `n_bus` indexing bug,
reorganized globals, revised error statements, removed unused input `bus`
* **form_jac:** formatting, removal of `jay`, reorganized globals, removed unused inputs `ang_red` and `volt_red`
* **freqcalc:** frequency measurement model developed by Felipe Wilches-Bernal
* **get_dypar:** new function that stores information about **PSTess** application settings
* **get_nstate:** new function that stores data about the number of states associated with each model
* **get_path:** new function that stores information about the **PSTess** data path and working file
* **i_simu:** formatting, removal of `jay`, changed **nc_load** solution tolerance to `tol = 1e-6`, reorganized globals
* **import_var:** a new script to import dynamic model parameters into the global struct `g` and clean up the workspace
* **ind_ldto:** formatting, removal of dummy variable `f`, reorganized globals
* **insimit:** formatting, removal of `dummy` in function call
* **inv_lf:** formatting, reorganized globals, revised error statements, fixed radian vs. degree bug in `gamma` specification
* **ivmmod_dyn:** new function that specifies the dynamics for voltage behind reactance converter models **mac_ivm**
* **lfdc:** formatting, batch processing implementation, removal of `jay`, reorganized globals, revised error statements, added `clear all` statement at the beginning of the script
* **lfdcs:** formatting, removal of `jay`, reorganized globals, revised error statements
* **lfdemo:** formatting, batch processing implementation, changed loadflow solution tolerance to `tol = 1e-8`, reorganized globals, revised error statements
* **lftap:** formatting, reorganized globals, revised warning statements
* **line_cur:** formatting, removal of `jay`, changed `nline` to `n_line` for consistency
* **line_pq:** formatting, removal of `jay`, changed `nline` to `n_line` for consistency
* **lmod:** formatting, removal of dummy variable `f`, corrections to anti-windup reset behavior, reorganized globals, removed unused input `bus`
* **lmod_indx:** formatting, renamed from **lm_indx** for clarity, removal of dummy variable `f`, reorganized globals, revised error statements
* **loadflow:** formatting, removal of `jay`, used `nargout` to control Jacobian output behavior, reorganized globals, revised error and warning statements
* **lsc:** new LTV synchronizing torque controller for use with inverter-based resources modeled using **ess**
* **lsc_indx:** new index function for LTV synchronizing torque controllers
* **mac_em:** formatting, removal of `jay` and dummy variable `f`, reorganized globals
* **mac_ib:** formatting, removal of `jay` and dummy variable `f`, reorganized globals
* **mac_igen:** formatting, removal of `jay`, reorganized globals, revised error statements
* **mac_ind:** formatting, removal of `jay`, reorganized globals, revised error statements
* **mac_ind_single_cage:** single-cage induction motor model (alternate version of **mac_ind**)
* **mac_indx:** formatting, removal of dummy variable `f`, reorganized globals
* **mac_ivm:** new voltage behind reactance converter model
* **mac_sub:** replaced the original file with a model developed by Dan Trudnowski and John Undrill based on **genrou** in PSLF
* **mac_sub_pst:** original **mac_sub** model from **PST v3.0**
* **mac_tra:** formatting, removal of `jay` and dummy variable `f`, correction of `vex` vs `vec` typo, reorganized globals, revised warning statements
* **mac_trip_logic:** logic for tripping generators (e.g., overspeed or out-of-step protection)
* **mdc_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **mess_sig:** specifies the auxiliary input for the energy storage-based control model
* **mexc_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **mlmod_sig:** formatting, renamed from **ml_sig** for clarity, removal of dummy variable `f`, reorganized globals
* **mpm_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **mrlmod_sig:** formatting, renamed from **rml_sig** for clarity, removal of dummy variable `f`, reorganized globals
* **msvc_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **mtcsc_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **mtg_sig:** formatting, removal of dummy variable `f`, reorganized globals
* **nc_load:** formatting, removal of `jay`, additions related to **ess**, rearranged mismatch calculation for clarity, reorganized globals, revised error statements
* **nm_if:** formatting, removal of dummy variable `f`
* **ns_file:** formatting, integration of `nss` struct stored by **get_nstate**, reorganized globals, revised error statements, added support for **ess** and **lsc** models
* **p_cont:** formatting, set `p_ratio = 1e-4`, switched to the central difference method, reorganized globals, added support for **ess** and **lsc** models
* **p_dpw:** formatting, removed hardcoded peturbation magnitude, switched to the central difference method, reorganized globals
* **p_exc:** formatting, corrected mislabeled state `Efd` to `st_name(k,j) = 10`, removed hardcoded peturbation magnitude, switched to the central difference method, reorganized globals
* **p_file:** formatting, switched to the central difference method, reorganized globals, revised error statements, added support for **ess** and **lsc** models
* **p_m_file:** formatting, rearranged repetitive code into loops, reorganized globals, added support for **ess** and **lsc** models to linearization routine
* **p_pss:** formatting, removed hardcoded peturbation magnitude, switched to the central difference method,
reorganized globals
* **p_tg:** formatting, removed hardcoded peturbation magnitude, switched to the central difference method, reorganized globals
* **pss:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **pss_des:** formatting, batch processing implementation
* **pss_indx:** formatting, removal of dummy variable `f`, reorganized globals, revised error statements
* **pss_phse:** formatting, removal of `jay` (in this case `i`)
* **pwrmod_dyn:** new function that specifies the dynamics for constant power/current injection models
* **pwrmod_dyn_alt:** new alternate example of **pwrmod_dyn**
* **pwrmod_indx:** new index function for constant power/current injection models
* **pwrmod_p:** new function for specifying real power injections
* **pwrmod_q:** new function for specifying reactive power injections
* **rbus_ang:** formatting, removal of `jay`, changed `nline` to `n_line` for consistency
* **rec_lf:** formatting, reorganized globals, revised error and warning statements
* **red_ybus:** formatting, optional additions related to **ess**, reorganized globals
* **rlmod:** formatting, removal of dummy variable `f`, corrections to anti-windup reset behavior, reorganized globals, removed unused input `bus`
* **rlmod_indx:** formatting, renamed from **rlm_indx** for clarity, removal of dummy variable `f`, reorganized globals, revised error statement
* **rltf:** formatting
* **rootloc:** formatting
* **s_simu:** formatting, batch processing implementation, removal of `jay`, multi-terminal hvdc additions by Ruisheng Diao, reorganized globals, revised error and warning statements, added support for **ess** and **lsc** models
* **sd_torque:** formatting
* **set_path:** new companion function that alters the file and/or path information stored in **get_path**
* **smpexc:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **smppi:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **stab_d:** formatting, batch processing implementation
* **stab_f:** formatting, removal of `jay` (in this case `i`)
* **state_f:** formatting, removal of `jay` (in this case `i`)
* **step_res:** formatting, revised error statements
* **svc:** formatting, removal of `jay`, reorganized globals, revised error statements
* **svc_indx:** formatting, removal of dummy variable `f`, reorganized globals, revised error statements
* **svc_sud:** formatting, corrections to deprecated `get()` syntax, revised error and warning statements
* **svm_mgen:** formatting, batch processing implementation, removal of `jay`, switched to the central difference method, made the code follow the outline of **s_simu** more closely, reorganized globals, revised error statements, added support for **ess** and **lsc** models
* **swcap:** formatting
* **tcsc:** formatting, removal of `jay` and dummy variable `f`, correction of `dB_cv` typo, reorganized globals, removed unused input `bus`
* **tcsc_indx:** formatting, removal of dummy variable `f`, reorganized globals, revised error_statements
* **tcsc_sud:** formatting, corrections to deprecated `get()` syntax, revised error and warning statements
* **tg:** formatting, removal of dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **tg_hydro:** formatting, removal of `jay` and dummy variable `f`, reorganized globals, revised error statements, removed unused input `bus`
* **tg_indx:** formatting, removal of dummy variable `f`, defined new variable `n_tg_tot` that stores the total number of turbine governors, reorganized globals
* **time_stamp:** formatting
* **vsdemo:** formatting, batch processing implementation, reorganized globals, revised error statements
* **y_sparse:** formatting, removal of `jay`, reorganized globals
* **y_switch:** formatting, implementation of switching sequence for resistive brake insertions, reorganized globals, revised error statements
* **ybus:** formatting, removal of `jay`, reorganized globals

[comment]: <> (eof)
