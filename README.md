Author: Dominic Brass (dombra@ceh.ac.uk)
Date: 05/03/24

A repository of files used to produce the results shown in the paper "Phenotypic plasticity in vector traits drives trends in global disease incidence".

Contents:

Dens_Temp_WL.xlsx   - An excel sheet containing all the life-history traits extracted from published laboratory experiments describing the 
		      life-history of Aedes albopictus. (Used in Parametrisation.Rmd)

Parametrisation.Rmd - An R markdown file that parametrises reaction norms describing the relationship between temperature, food availability
		      and the life-history traits of Aedes albopictus (Requires Dens_Temp_WL.xlsx).

LDgam.csv           - CSV of predictions made by the GAMM describing the relationship between larval rearing temperature and food availability on
		      the duration of the larval stage. (Output from Parametrisation.Rmd and used as input in AO2_public.jl & SO2_public.jl)

LSgam.csv           - CSV of predictions made by the GAMM describing the relationship between larval rearing temperature and food availability on
		      the survival through the larval stage. (Output from Parametrisation.Rmd and used as input in AO2_public.jl & SO2_public.jl)

WL_re.csv           - CSV of predictions made by the GAMM describing the relationship between larval rearing temperature and food availability on
		      the wing length of adults. (Output from Parametrisation.Rmd and used as input in AO2_public.jl & SO2_public.jl)

AMgam.csv           - CSV of predictions made by the GAMM describing the relationship between the current temperature and adult wing length on
		      the survival of adult mosquitoes (Output from Parametrisation.Rmd and used as input in AO2_public.jl & SO2_public.jl)

RO2_example.jl      - Example Julia script uses the files here to produce the predictions for Rimini, Italy. Note that this requires the user to supply
		      climate data from ERA5 which cannot be shared here due to liscence agreements. The form that this data should take is indicated within
		      this script with links to the repositories used. (Requires LDgam.csv, LSgam.csv, WL_re.csv, AMgam.csv, RMSE_fit.jl, Rim_Clim.csv*, Car_dat_Rim.csv*
		      and AO2_public.jl) 

AO2_public.jl       - Julia script called by RO2_public.jl, describes the stage/phenotypically structured delay-differential equation model 
		      uesed to predict the population dynamics of Aedes albopictus. These dynamics are the same as those used in SO2_public.jl but do 
	              not include the dengue dynamics for speed of computation. (Requires LDgam.csv, LSgam.csv, WL_re.csv, AMgam.csv, and functions 
		      such as are described in RO2_public.jl) 

RMSE_fit.jl         - Julia script defining a function the compares the fit between model predictions and field observations.

R_t_comp.jl         - Julia script used to produce the results shown in Figure 2. Note that this requires the user to supply climate data from ERA5 which 
                      cannot be shared here due to liscence agreements. The form that this data should take is indicated within this script with links to 
                      the repositories used. (Requires LDgam.csv, LSgam.csv, WL_re.csv, AMgam.csv, R_t_calc.jl, Gua_Clim.csv*, Cag_Clim.csv*,Tok_Clim.csv*, 
                      SB_Clim.csv*, Kon_Clim.csv*, and SO2_public.jl) 

R_t_calc.jl         - Julia script that computes R_t from output of SO2_public.jl

SO2_public.jl       - Julia script called by R_t_comp.jl, used to run the stage/phenotypically structured delay-differential equation model 
		      describing the population dynamics of Aedes albopictus. (Requires LDgam.csv, LSgam.csv, WL_re.csv, AMgam.csv, and functions 
		      described in RO2_public.jl)                                 
 
* - Not provided here
