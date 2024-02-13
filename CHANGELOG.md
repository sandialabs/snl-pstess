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
- [Version 1.1](#version_1p1)
    - [Bug fixes](#bug_fixes)
    - [Miscellaneous](#miscellaneous)
    - [File modifications](#file_mods)

## Version 1.1
<a id="version_1p1"></a>

This section describes the changes, on a per-file basis, made to
existing functions and scripts in the creation of **PSTess** Version 1.1.
Please see the README for context and additional information about new
features added to the toolbox.

Note: The **PST v3.0** source code on which **PSTess** is based is
available in the main `pstess` directory under the folder titled
`archive_pstV3p0`.

### Bug fixes
<a id="bug_fixes"></a>

* **mac_ivm**: the original function mistakenly scaled the generator
  terminal current onto the system base rather than the generator base;
  this has been corrected. See the file modifications section for
  additional changes.

### Miscellaneous
<a id="miscellaneous"></a>

* The main purpose of this update is to add support for the
  grid-following inverter model REEC_D and the grid-forming
  model REGFM_A1. For details of their implementation, see
  `reec`, `ivm`, and `gfma`.

### File modifications
<a id="file_mods"></a>

* **gfma:** new model for representing REGFM_A1 droop-controlled grid-forming inverters
* **gfma_indx:** new function for indexing `gfma` instances
* **ivm_indx:** new function for indexing `ivm` instances
* **ivm_sud:** new function for implementing user-defined grid-forming inverter control
* **ess:** implemented the voltage-dependent current limit (VDL) curves specified in REEC_D; can now turn off SOC limits by setting `ess_con(,7) <= 0`; can now turn off LVPL function when `ess_con(,17) <= 0`
* **ess_indx:** implemented the voltage-dependent current limit (VDL) curves specified in REEC_D
* **ess_sat:** implemented the voltage-dependent current limit (VDL) curves and protection functions specified in REEC_D
* **get_nstate:** added support for `ivm`, `gfma`, and `reec` models
* **import_var:** added support for `ivm`, `gfma`, and `reec` models
* **lsc:** corrected typo in header
* **mac_em:** allowed `pm_sig` modulation when `i > 0`
* **mac_ib:** fixed bug where `psikd`, `psikq`, `edprime`, and `eqprime` were not being updated for `mac_sub` models in vectorized computation; allowed `mac_ivm` generators to be modeled as infinite buses (Beta feature, experimental)
* **mac_indx:** implemented error checking for infinite bus dimensions; changed how `ivm` generators are detected to account for VSMs with nonzero inertia
* **mac_ivm:** updated interface model for representing grid-forming inverters; the voltage magnitude state is now `eqprime` for consistency with other generator models; the commanded voltage magnitude is set in `vex` and the commanded angle in `fldcur`
* **mgfma_sig.m:** new function for perturbing `gfma` models
* **mreec_sig.m:** new function for perturbing `reec` models
* **nc_load:** added a new hook to allow tripping of non-conforming loads (Beta feature, experimental)
* **ns_file:** added support for `ivm`, `gfma`, and `reec` models
* **p_cont:** added support for `ivm`, `gfma`, and `reec` models
* **p_file:** added support for `ivm`, `gfma`, and `reec` models
* **p_m_file:** added support for `ivm`, `gfma`, and `reec` models
* **reec:** new model for representing REEC_D grid-following inverters
* **reec_indx:** new function for indexing `reec` instances
* **s_simu:** added support for `ivm`, `gfma`, and `reec` models
* **svm_mgen:** added support for `ivm`, `gfma`, and `reec` models
* **trip_handler.m** new helper function for implementing protection and RAS logic
* **trip_indx.m** new function for indexing trippable loads, e.g., for UFLS
* **trip_logic.m** new user-defined function for implementing protection and RAS logic (Beta feature, experimental)
