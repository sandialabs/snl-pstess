[comment]: <> (this is a markdown document and the special characters are for formatting)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="man/logo/snl-pstess_logo_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="man/logo/snl-pstess_logo.png">
  <img src="man/logo/snl-pstess_logo.png" alt="PSTess logo" width=300px margin="auto" />
</picture>

# PSTess: The Power and Energy Storage Systems Toolbox

Welcome to the Power and Energy Storage Systems Toolbox. This is the
**PSTess** README. Please refer to the table of contents below to help
you navigate.

Current release version: 1.1

Release date: February, 2024

## Contact
For issues and feedback, we would appreciate it if you could use the
"Issues" feature of this repository. This helps others join the
discussion and helps us keep track of and document issues.

### Email
Entity account `@sandia.gov: snl-pstess`<br />
Project maintainer (Ryan Elliott) `@sandia.gov: rtellio`<br />
Project co-maintainer (Hyungjin Choi) `@sandia.gov: hchoi`

### Partners
The development of **PSTess** was a collaborative effort between
Sandia National Laboratories and Montana Tech. We thank Dr. Daniel
Trudnowski and Tam Nguyen for their contributions.

## Table of contents
- [Introduction](#intro)
    - [Models](#models)
    - [Features](#features)
    - [Acknowledgment](#acknowledgment)
- [Getting started](#getting_started)
    - [Setting the path](#set_path)
    - [Running an example](#first_example)
- [Frequently Asked Questions](#faq)
    - [General](#general)
- [References](#references)

## What is it?
<a id="intro"></a>

**PSTess** is an open-source, MATLAB-based toolbox for dynamic
simulation and analysis of power systems with utility-scale,
inverter-based energy storage systems (ESS). Of course, it can also be
used to study conventional power systems. **PSTess** is a fork of the
Power System Toolbox, called **PST** for short. It is based on **PST
v3.0**, released by Rensselaer Polytechnic Institute (RPI) in August
of 2020. The major differences between the two packages are described
below. For a detailed description of the modifications made to **PST**
source files, see the CHANGELOG.

### Models
<a id="models"></a>

* **ess:** energy storage device model with reactive modulation
  capability
* **ess_sud:** user-defined energy storage control model; the pipeline
  has been customized by Montana Tech
* **lsc:** linear time-varying (LTV) synchronizing torque controller
  (for use with **ess**)
* **pwrmod:** constant power (or current) injection model
* **ivmmod:** experimental voltage behind reactance converter model
  (primarily for studying grid-forming inverters)

### Features
<a id="features"></a>

* Global variables: all global variables have been reorganized into a
  common struct, called `g`
* Batch processing: a new mode of operation has been added to modify
  application-wide behavior when performing multiple simulation or
  analysis runs (see **get_dypar**)
* Error statements: all errors and warnings issued by the application
  have been edited to be more specific, consistent, and useful
* Dynamic brake insertions: **y_switch** has been updated to support
  three-phase resistive shunt insertions (e.g., for modeling the Chief
  Joseph Brake); this option is now available as `f_type = 8`
* Linearization routine: **svm_mgen** has been redesigned around the
  central difference method for improved accuracy; the surrounding
  code has also been modified to better track the time-domain
  simulation procedure in **s_simu**

### Acknowledgment
<a id="acknowledgment"></a>

This work was supported by the U.S. DOE Energy Storage Program within
the Office of Electricity. The authors would like to thank Dr. Imre
Gyuk, Director of Energy Storage Research.

This toolbox owes a great deal to the original developers of **PST**:
Drs. Joe Chow, Kwok Cheung, and Graham Rogers. We thank them for
their contributions and for their generosity in releasing **PST**
under the MIT license, enabling others to build upon the foundation
they created. See the NOTICE file in this repository for details about
the **PST** license.

## Getting started
<a id="getting_started"></a>

Before working with **PSTess**, please see the documentation listed
under the `man` directory of this repository. In addition, the
CHANGELOG includes detailed information about changes made to **PST**
functions and scripts in the creation of **PSTess**.

### Setting the path
<a id="set_path"></a>

To set the **PSTess** path and or working data file utilized in batch
mode, they may be entered directly in `get_path.m`, or set
programmatically via `set_path.m`. Note that calls to `set_path.m`
re-write the `get_path.m` file, so they are durable in
response to `clear` calls. The path and data file specified in this
way are used by **PSTess** when `dypar.batch_mode = true`, as
specified in `get_dypar.m`. To turn batch mode off, simply set
`dypar.batch_mode = false`. When batch mode is turned off, the
application will prompt the user for information, as necessary.

### Running an example
<a id="first_example"></a>

Please see the `analysis` folder in this repository for two examples
based on a version of the Kundur 2-area system augmented with energy storage.

* `main_ess_example1.m` -- Perturbs the real power command of the
  energy storage systems, as specified in `mess_sig_example1.m`.
* `main_ess_example2.m` -- Perturbs the reactive power command of the
  energy storage systems, as specified in `mess_sig_example2.m`.

Before running either of the examples, read the information in the
header of the main script. The instructions for `main_ess_example1.m`
are reprinted below for reference.

    % Note: You must run this example from the 'analysis' folder. Before
    %       running, rename the original /pstess/mess_sig.m file to
    %       mess_sig_original.m. Then move mess_sig_example1.m into /pstess
    %       and rename it mess_sig.m. Make sure the working directory is
    %       configured in get_path.m. Before running double check that sw_con
    %       in d2asbegp_ess.m corresponds to ftype=6 (do nothing).

After running `main_ess_example1.m`, a small collection of figures
will be automatically generated and saved in the `analysis/fig`
directory. The figure on the left below shows the bus frequency
deviation vs. time in response to the input perturbation. The solid
blue trace shows the output of the full nonlinear simulation, and the
dashed magenta trace the output of the linearized simulation. The plot
on the right shows the same comparison but for the real power
modulation by the ESSs. The subscripts on the y-axis labels correspond
to the ESS index numbers. Both examples include the ESS-based
transient stability control strategy described in a recent
[IET journal article](https://digital-library.theiet.org/content/journals/10.1049/iet-gtd.2020.1319).

Bus frequency deviation    |  Real power modulation
:-------------------------:|:-------------------------:
![Bus frequency](/analysis/fig/ess_example1_bus_freq_web.png) | ![Real power command](/analysis/fig/ess_example1_pwr_cmd_web.png)

## Frequently asked questions
<a id="faq"></a>

### General
<a id="general"></a>

> Where can I find documentation for PST and PSTess?

Please see the `man` folder under the top level of the repository. It
contains documentation for **PST v3.0** and a supplemental document
describing the features of **PSTess**. For users without prior **PST**
experience, we strongly recommend reading both of these documents
before working with **PSTess**.

> I think I found a bug, what should I do?

Please use the "Issues" feature of this repository. This helps others
join the discussion and helps us keep track of and document issues. We
will do what we can to address these issues; however, the **PST** base
application has both known and unknown issues that users will
encounter. Users of **PSTess** should be comfortable enough with
MATLAB programming to address these issues on their own as they arise.

## References
<a id="references"></a>

R. Elliott, D. Trudnowski, H. Choi, T. Nguyen, "The Power and Energy Storage
Systems Toolbox – PSTess Version 1.0" September 2021, SAND2021-11259. Available:
[https://www.sandia.gov/ess-ssl/](https://www.sandia.gov/ess-ssl/wp-content/uploads/2021/09/sand2021_11259_pstess_elliott_choi.pdf)

J. Chow and K. Cheung, "A toolbox for power system dynamics and
control engineering education and research," IEEE Trans. Power Syst.,
vol. 7, no. 4, pp. 1559–1564, 1992.

J. Chow, "Power System Toolbox Version 3.0 manual." Available:
https://www.ecse.rpi.edu/~chowj/PSTMan.pdf

R. Elliott, A. Ellis, P. Pourbeik, J. Sanchez-Gasca, J. Senthil,
and J. Weber, "Generic photovoltaic system models for WECC–A status
report," in IEEE Power Energy Soc. Gen. Meeting, pp. 1–5, 2015.

R. T. Elliott, P. Arabshahi, and D. S. Kirschen, "Stabilizing
transient disturbances with utility-scale inverter-based resources,"
IET Gen. Tran. Dist., vol. 14, pp. 6534–6544, Dec. 2020. Available:
[https://ietresearch.onlinelibrary.wiley.com/doi/](https://ietresearch.onlinelibrary.wiley.com/doi/epdf/10.1049/iet-gtd.2020.1319)

CAISO, "Inverter-based interconnection requirements (ER19-1153)."
Available: [https://www.caiso.com/Documents/](http://www.caiso.com/Documents/Feb28-2019-InterconnectionProcessEnhancements-Inverter-BasedGeneratorInterconnectionRequirements-ER19-1153.pdf)

M. Klein, G. Rogers, and P. Kundur, "A fundamental study of
inter-area oscillations in power systems," IEEE Trans. Power Syst.,
vol. 6, no. 3, pp. 914–921, 1991.

D. Trudnowski, D. Kosterev, and J. Undrill, "PDCI damping control
analysis for the western North American power system," in IEEE Power
Energy Soc. Gen. Meeting, pp. 1–5, 2013.

[comment]: <> (eof)
