# User guide 
## Compiling the model and creating a project
The only prerequisite to install LAKE model is Fortran compiler.
The model can be built from source code by GNU and Intel Fortran compilers. GNU compiler is set by default in the model install configuration. Install configuration is given in files [Makefile](/../blob/master/Makefile) and [source/Makefile](/../blob/master/source/Makefile).  
If you want to use Intel compiler, you have to change variable `FC` value to `FC = ifort` in these files.  
Once you unpack the LAKE model archive or cloned its git project just type in the model root folder  
```
user@computer:~/LAKE/lake$ make all
```
to build the model from the source code. The next step is to configure the
model for you needs. This means that you create a project, identified
by the project name hereafter we will refer to it as ```project_name```,
it may be any text, e.g. combination of the lake name and the time span you need to simulate this lake, for instance, "Valkea-Kotinen2006"). Once you've chosen the project name, type
```
user@computer:~/LAKE/lake$ ./crproj <project_name>
```

This creates the folder structure in `./results`, to which the output
files will be written. After that, create `./setup/project_name_driver.dat` and `./setup/project_name_setup.dat` files. These files must have
a certain structure, so that we recommend that you create them by copying from
sample files already existing in `./setup`. Once copied, you can modify them appropriately. This is what the next two sections are devoted to.

To run the model you have to provide a file with meteorological data. It must be an ASCII table, where each column contains a time series of a certain meteorological variable. A sequence of columns is set in `/setup/project_name_driver.dat`. Shortly after start, the model will seek the file with meteorological data as `./data/project_name.dat`.

## Configuring the driver file
Driver file is used to set the outer environment of the LAKE model
that is normally the atmospheric model. This means, that atmospheric
model calculates the atmospheric forcing for the LAKE model, and some
external parameters are also being passed to the LAKE model through
the driver interface. When the LAKE model is coupled to atmospheric 
model the driver file is not needed.  
The format of the driver file is following. Each string of the file
starts with a comment (indicated by `#`) or with a keyword. The last
string of the file contains the keyword `end`. The keywords may
be listed in any order, but we recommend to keep them in thematic groups,
as it is done in sample files. Each keyword, excepting `end`, must 
have a value that is set after the keyword at the same string separated
by at least one space. Some keywords are followed by a group of numbers.  
### Parameters in driver file
| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*dataname* | The project name, character, enclosed by `' '` |
|*forc_format* | The input file format: 0 - ASCII(text), 1 - netcdf |
|*npoints* | Number of lakes to be simulated |
|*select_call*, opt | The list of up to 20 numbers, each at new line, denoting particular lakes, that will be simulated. Relevant only if *npoints* \> 1. Other lakes are omitted in the simulation. |
|*lakinterac* | The number of interacting lakes: 1 - no interaction, 2 - two lakes are linked by a channel, \> 2 - not operational. |
|*form* | obsolete, put 0 (means adjustable ASCII format, see below) |
|*height_T_q* | The height above the lake surface at which temperature and humidity are given (calculated by the atmospheric model or measured), m|
|*height_u* | The height above the lake surface at which the wind speed is given (calculated by the atmospheric model or measured), m|
|*interval* | The time step of meteorological data, s |
|*rad* | Defines, atmospheric radiation (1) or net radiation (2) is in the relevant column of datafile, relevant if input file is ASCII (forc_format = 0) |

**The set of parameters for adjustable format of input text file, relevant if forc_format = 0, form = 0**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*N_header_lines* | The number of lines in a file header  |
|*N_coloumns* | The total number of coloumns in a file |
|*N_Year* | The number of coloumn with the number of a year (not used in the model) |
|*N_Month* | The number of coloumn with the number of a month (not used in the model) |
|*N_Day* | The number of coloumn with the number of a day (not used in the model) |
|*N_Hour* | The number of coloumn with the number of a hour (not used in the model) |
|*N_Uspeed* | The number of coloumn with x-component speed, m/s |
|*N_Vspeed* | The number of coloumn with y-component speed, m/s |
|*N_Temp* | The number of coloumn with air temperature, K |
|*N_Hum* | The number of coloumn with air humidity, kg/kg |
|*N_Pres* | The number of coloumn with atmospheric pressure, Pa |
|*N_SWdown* | The number of coloumn with total solar radiation, W/m^2 |
|*N_LWdown* | The number of coloumn with atmospheric radiation, W/m^2) |
|*N_Precip* | The number of coloumn with precipitation, m/s |
|*N_SensFlux* | The number of coloumn with sensible heat flux, W/m^2 |
|*N_LatentFlux* | The number of coloumn with latent heat flux, W/m^2 |
|*N_Ustar* | The number of coloumn with friction velocity, m/s |
|*N_surfrad* | The number of coloumn with radiation emitted by the surface, W/m^2 |
|*N_cloud* | The number of coloumn with total cloudiness, fraction |

**The parameters of time integration**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*year0* | Julian year of initial moment of integration |
|*month0* | Julian month of initial moment of integration |
|*day0* | Julian day of initial moment of integration |
|*hour0* | Julian hour (may be fraction) of initial moment of integration |
|*tinteg* | The duration of integration, days |
|*spinup_times* | Number of spinup periods |
|*spinup_period* | The duration of spinup period, days |
|*dt* | Time step of integration, s |
|*call_Flake* | The switch for integrating FLake model: 1 - on, 0 - off |

**External physical parameters of a lake**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*extwat* | Solar radiation extinction in a water body, m^-1 |
|*extice* | Solar radiation extinction in an ice, m^-1 |
|*alphax* | Slope of water surface in the x-direction, deg |
|*alphay* | Slope of water surface in the y-direction, deg |
|*a_veg* | Effective cross-section of vegetation, m^2/m^2 |
|*c_veg* | Friction coefficient of vegetation in water, n/d |
|*h_veg* | The height of vegetation in the lake, m |
|*kor* | Coriolis parameter, s^-1; if -999., then it is computed from latitude |
|*phi* | Latitude, deg |
|*lam* | Longitude, deg |
|*fetch* | Average wind fetch, m |
|*area_lake* | Lake area, m^2 |
|*cellipt* | The ratio of x dimension to y dimension of a lake, n/d |
|*lakeform* | The form of the lake's horizontal cross-section: 1 -- ellipse, 2 -- rectangle |
|*trib_inflow* | Total tributaries' inflow, m^3/s; if -9999., then
it is dynamically adjusted to keep the lake depth nearly constant |
|*morphometry* | Depth - lake cross-section area table, m^2; if not given, the area is constant with depth |
|*effl_outflow* | effluent discharge parameters group, specifying polynomial dependence of discharge on water level; the last value is relative altitude of effluent bottom over lake bottom (at deepest points, respectively); the first N values are polynomial coefficients, where N stands after *effl_outflow* keyword |

**Initial conditions keywords**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*l10* | Initial ice thickness, m |
|*h10* | Initial water thickness (lake depth), m |
|*select_h10*, opt | The list of up to 20 depths, each at new line, specifying lakes' depths, for lakes selected by *select_call*. |
|*ls10* | Initial bottom ice thickness, m |
|*hs10* | Initial snow thickness, m |
|*Ts0* | Initial temperature of the top mixed layer, °C |
|*Tb0* | Initial bottom temperature, °C |
|*Tbb0* | Initial temperature at the bottom of soil layer, °C |
|*Tbb0* | Initial mean temperature of the water column, °C |
|*h_ML0* | Initial mixed-layer depth, m |
|*Sals0* | Initial salinity in the mixed layer, kg/kg |
|*Salb0* | Initial bottom salinity, kg/kg |
|*us0* | Initial x-component of velocity in the mixed layer, m/s |
|*vs0* | Initial y-component of velocity in the mixed layer, m/s |
|*init_T* | The type of profiles initialization: 1 - using *h_ML0*, *Ts0* and *Tb0*; 2 - using Tm,    Ts0 and Tb0; 3 - using the temperature profile, specified in `\*_setup.dat` or   file specified after keyword *T_profile*; The temperature, salinity, CH_4,CO_2,O_2 will be initialized via profiles. *init_T*=2 works only for temperature |

**Output controlling parameters**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*nstep_ncout* | the interval of netcdf output from driver, timesteps (if -1 no netcdf output from driver) |
|*nstep_out_Flake* | the interval of output of Flake variables from driver, timesteps (if -1 the output of FLake variables from driver is not implemented), relevant if *call_Flake* = 1 |
|*moving_average_window* | the interval of moving average, netcdf output steps (intervals) |
|*mean_cycle_period* | the length of mean cycle, netcdf output steps (intervals) |

## Configuring the setup file
Setup file has the same format, as the driver file. The list of parameters is given in a table below. 

### The parameters of setup file

**General parameters**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*path* | The full path to the model root, needed on some systems; normally, specify it as `''` |
|*runmode* | The mode of LAKE model: 1 - standalone, 2 - driven by atmospheric model|
|*omp* | The OpenMP: 0 - not used, 1 - used; normally OpenMP does not give significant speedup |

**Spatio-temporal resolution group**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*M* | Number of layers in water |
|*ns* | Number of layers in soil |
|*nsoilcols* | Number of soil columns (1 - single soil column contacting with the lowest numerical layer (of water or deep ice) only, \> 1 - additional soil columns at the bottom slope, contacting with all numerical layers of water, ice and deep ice; morphometry information is necessary in driver file, if *nsoilcols* \> 1) |
|*Mice* | Number of layers in bottom ice |
|*d_surf* | Grid zooming parameter at the surface, n/d |
|*d_bot* | Grid zooming parameter at the bottom, n/d |

**Physical parameterizations**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*varalb* | Water albedo: 0 - constant, 1 - zenith-angle dependent |
|*PBLpar* | Surface fluxes scheme: -1 - sensible, latent heat and momentum fluxes are given as input for the model; 0 - the latent heat flux is set to zero, while sensible heat and momentum fluxes are constant in time, specified by *sensflux0* and *momflux0*; 1 - Businger-Dayer formulas (Monin-Obukhov theory) for exchange coefficients, the code from INM RAS climate model; 2 - formulation from NH3d (Louis, 1979); 3 - formulation from FLake;
4 -  formulation implemented by D.Chechin |
|*waveenh* | The shallow water correction of surface fluxes (Panin et al., 1996) (0 - off, 1 - on) |
|*momflxpart* | Momentum flux treatment: 0 - all momentum flux from the atmosphere is consumed by currents acceleration; 2 - momentum flux from the atmosphere is partitioned between wave development (controlled by fetch) and currents acceleration (Stepanenko et al., 2014) |
|*c_d* | The momentum exchange coefficient, n/d; if -999, momentum flux is calculated by surface flux scheme |
|*kwe* | The factor of turbulence enhancement by wave breaking (wave energy factor), n/d (used in boundary conditions of $`K-\epsilon`$ model) |
|*relwind* | In momentum flux calculation the currents speed is substracted from wind speed: 1 - on, 0 - off |
|*eos* | Equation of state: 1 - from Hostetler model (freshwater), 2 - TEOS-2010 (currently not operational), 3 - for Kivu lake (Schmid et al., 2003), 4 - UNESCO formula (freshwater), 5 - McCatcheon et al., 1993  |
|*lindens* | Switch for using linearized equation of state respective to temperature and salinity: 1 - on, 0 - off |
|*lindens* | Linear dependency of density on temperature and salinity: 0 - off, 1 - on  |
|*nmeltpoint* | Melting point: 0 - constant, 0 °C, 1 - linearly dependent on salinity, 2 - TEOS-2010 formula |
|*Turbpar* | Turbulence model: 2 - $` K-\epsilon `$ closure; 8 - Hendersson-Sellers diffusivity (1985); 9 - Henderson-Sellers diffusivity adopted for shallow lakes |
|*stabfunc* | Stability function in $` K-\epsilon `$ : 1 - constant coefficients (standard K-\epsilon model); 2 - stability functions according to (Canuto et al., 2001); 3 - stability functions according to (Galperin et al., 1988) |
|*kepsbc* | Boundary conditions in K-\epsilon model: 1 - Neuman boundary conditions for unstratified sheared flow (Burchard, 2002); 2 -  Neuman boundary conditions for unstratified non-sheared flow with wave breaking; 3 - Neuman boundary conditions for unstratified sheared flow with wave breaking; 4 - Neuman boundary conditions for free convection |
|*soiltype* | The soil type: 1 - sand, 2 - loamy sand, 3 - sandy loam, 4 - loam, 5 - silt loam, 6 - sandy clay loam, 7 - clay loam, 8 - silty clay loam, 9 - sandy clay, 10 - silty clay, 11 - clay |
|*soil_depth* | The depth of the soil layer, m |
|*soilswitch* | The switch for soil model: 0 - off, 1 - on |
|*tricemethhydr* | The ice treatment in the soil pores: 0. - treat as pure ice; 1. - treat as methane hydrate |
|*skin* | The skin temperature parameterization: 0 - off, 1 - on |
|*massflux* | The massflux parameterization of convection (Siebesma et al., 2007): 0 - off, 1 - on. |
|*ifrad* | Radiation fluxes at the water surface: 0 - set to zero; 1 - full treatment |
|*ifbubble* | Switch for the bubble model: 0 - off (bubbles do not change with depth); 1 - on |
|*carbon_model* | Switch for the carbon model: 1 - Stefan and Fang model; 2 - Hanson et al. model |
|*sedim* | Gravitational sedimentation of tracer in water: 0 - off, 1 - on |
|*salsoil* | Switch for salinity transport in soil: 0 - off (zero flux at the bottom), 1 - on |
|*dyn_pgrad* | Dynamic pressure gradient: 0 - off, 1 - barotropic pressure gradient, 2 - two-layer approximation for pressure gradient, 3 - multilayer model for pressure gradient |
|*zero_model* | Switch for zero-dimensional model: 0 -off, 1 - on |
|
*thermokarst_meth_prod* | Switch for old organics methane production under thermokarst lakes: 1. - on, 0. - off |
|*soil_meth_prod* | Switch for new organics methane production under lakes: 1. - on, 0. - off |
|*outflpar* | The treatment of effluent temperature: 0 - value at the outflow = cross-section mean; 1 - the cross-section mean = 0.5*(inflow value + outflow value); 2 - variables at the outflow are calculated using Lagrangian approach |
|*sensflux0* | Sensible heat flux upwards, constant in time (relevant if PBLpar = 0), W/m^2 |
|*momflux0* | Momentum flux downwards (positive), constant in time (relevant if PBLpar = 0), N/m^2 |
|*soilbotflx* | The constant heat flux at the soil or lake bottom (depending on if soil is switched on or off), W/m^2 |
|*cuette* | The switch for Cuette flow boundary conditions: 0 - no, 1 - Dirichlet conditions for temperature and velocities at both top and bottom boundaries, 11 - the same as 1, but momentum flux is imposed at the top boundary; **Note!** Other switches for physical parameterizations and initial conditions must be adjusted to the Cuette flow configuration |
|*deadvol* | Variables at the outflow are calculated using Lagrangian approach |

**Initial conditions**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*T_profile* | Initial temperature (°C, col.2), salinity (kg/kg, col.3), CH_4 (mol/m^3, col.4), CO_2 (mol/m^3, col.5), O_2 (mol/m^3, col.6) profiles, at depth levels (m, col.1) |

**Tributaries and effluents**

| Keyword | Description and acceptable values set |
|:---------------|:-----------------------------------------------------------|
|*tribheat* | The switch for thermal effect of tributaries and effluents: 0 - off, 1 - on, constant inflow profiles, 2 - on, time-dependent inflow profiles |
|*N_tribin* | The number of tributaries, the next line should be a sequence of *N_tribin* numbers, from 1 to 4, describing the location of each tributary |
|*N_triblev* | The number of levels in tributaries data, relevant if tribheat = 1 or 2 |
|*iefflloc* | from 1 to 4, the location of an effluent |

## Running the model
The model is launched by typing in the terminal
```
user@computer:~/LAKE/lake$ ./lake.out
```

However, it usually happen that you need to run the model with different configurations (input parameters). In this case we recommend to use the script `launch` in a root folder. To run the model you invoke this script with an argument specifying the experiment label
```
user@computer:~/LAKE/lake$ ./launch <project_name>_<spec>
```
for instance
```
user@computer:~/LAKE/lake$ ./launch Valkea-Kotinen2006_base
```
This script does the following. It runs the model, the model writes output data to `results/Valkea-Kotinen2006`, then the script creates `results/Valkea-Kotinen2006_base` and copies there the content of `results/Valkea-Kotinen2006`, `setup/Valkea-Kotinen2006_setup.dat` and `setup/Valkea-Kotinen2006_driver.dat`. As a result, both model output and configuration files of an experiment labeled as `base` are written in a specific directory, so you will not mix them with that from other experiments. 

## Model output
The model output is written in `./results/project_name` folder.
```
user@computer:~/LAKE/lake/results/project\_name$ ls
daily  everystep  hourly  monthly  netcdf  time_series
user@computer:~/LAKE/lake/results/project\_name$
```
The folders `hourly`, `daily`, `monthly` contain vertical profiles
of different variables in a lake, averaged over respective time intervals. Output files in netCDF format are being placed in `netcdf` folder, and files 
with time series - in `time_series` folder. Sometimes, for debugging purposes,
it is advantageous to use the output at every time step of the model's numerical scheme, the respective files reside in `time_series`

 