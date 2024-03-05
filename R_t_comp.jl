## R_t_comp.jl
#Author: Dominic Brass (dombra@ceh.ac.uk)
#Date: 04/03/2024

#Code accompanying the paper, "Phenotypic plasticity in vector traits drives trends in global disease incidence".
#Code used to produce the predictions of disease dynamics in Figure 2. Due to the ERA5 license agreement
#the climate files used to produce the paper outputs cannot be shared here but this data is freely available. To perform the comparison between 
#model prediction observed cases additional files are required as specified in the Supplementary Materials. 

#Requires:
#SO2_public.jl
#Gua_Clim.csv
#Cag_Clim.csv
#Tok_Clim.csv
#SB_Clim.csv
#Kon_Clim.csv

#Loads packages
using(DifferentialEquations)
using(Dierckx)
using(Plots)
using(CSV)
using(Interpolations)
using(QuadGK)
using(Statistics)
using(Dates)
using(DataFrames)
using(LaTeXStrings)

pal = get_color_palette(:auto, plot_color(:white)) #Defines colour pallette

##Guangzhou 2013 & 2014

lat   = 23.1291        #Latitude of the simulation
tspan = (0.0, 1461.0); #Time span to simulate
pop   = 32113.22 * 4;  # Population density of Guangzhou per 4km^2

#Introduction of infected humans to the population
function infect(time)
  if mod(time,365) >= 151  && mod(time,365) <= 190 && time > 1095
     impulse = 0.25
 elseif mod(time,365) >= 150  && mod(time,365) <= 151  && time < 1095 && time > 730
     impulse =0.25
   else
     impulse = 0
   end
   return(impulse)
end;
 
Guang_Clim = CSV.read(raw"Gua_Clim.csv", header=true,DataFrame) #Loads climate of Guangzhou
##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.

Guang_Temp = Guang_Clim[:,3]
Guang_Prec = Guang_Clim[:,4]
Guang_EV   = Guang_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)), Guang_Prec[4019:5479])
Temps_spl  = Spline1D(collect(range(1,1461,length = 1461)) ,Guang_Temp[4019:5479])
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)) ,  Guang_EV[4019:5479])

#Gua_D_13   = CSV.read(raw"C:\Users\dombra\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Gua_13_R.csv", header=false,DataFrame)
#Gua_D_14   = CSV.read(raw"C:\Users\dombra\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\G_D_14.csv", header=false,DataFrame)

include("SO2_public.jl") #Runs SIR model
 
#Processes model outputs
solgua14       = sol.t
inf_sum_2014   = inf_sum
adult_sum_2014 = adult_sum
H_I_2014       = out[numb_H_I,:]
H_R_2014       = out[numb_H_R,:]
H_S_2014       = out[numb_H_S,:]
out_gua        = out

#Number of new infections produced each day
Daily_inf_2014 = (vector_human.(Temp.(sol.t .- 4)) .* bite.(Temp.(sol.t .- 4)))   .*   40 .* 4 .* inf_sum .* out[numb_H_S,:] ./ (out[numb_H_S,:] .+ out[numb_H_I,:] .+ out[numb_H_R,:] .+ 15000   )

#Plots outputs
dt = DateTime(2011,1,1)
      min_date = DateTime(2014,1,1)
      max_date = DateTime(2015,1,1)
      time_axis = add_float2datetime.(solgua14,dt)
      time_axis2 = add_float2datetime.(1:1461,dt)

gH_I_2014 = Plots.plot(time_axis, H_I_2014 ,label = "Simulated dengue cases", xlabel = "Time (Days)", ylabel = "Dengue Cases")
                      plot!(xlims = Dates.value.([min_date,max_date]))
#                      plot!(add_float2datetime.(Gua_D_14[:,1]  .+ 1095 ,dt),(Gua_D_14[:,2]), label = "Observed dengue cases")
                      plot!(time_axis, H_S_2014, legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")
                      plot!(time_axis, H_R_2014 , legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")

daily_2014 = plot(time_axis, (Daily_inf_2014 .* 13.5),label = "Predicted cases" , ylabel = "New cases per day", linewidth = 1)
                  plot!(xlims = Dates.value.([min_date,max_date]))
 #                 plot!(add_float2datetime.(Gua_D_14[:,1] .+ 1095    ,dt),(Gua_D_14[:,2]) ,label = "Observed cases" , linewidth = 1, color = pal[2])
 #                 plot!(add_float2datetime.(Gua_D_13[:,1] .+ 730 .+ 14,dt),exp.(Gua_D_13[:,2]) ./ (31) , label = nothing, linewidth =1 , color = pal[2])

 #Calculates R_t for Guangzhou
include("R_t_calc.jl")                 
                  
R0_spl_gua = Spline1D(collect(range(1,1461,length = 1461)),sum(R_0',dims = 2)[:,1])

plot(time_axis, R0_spl_gua.(sol.t))

#Alpes-Martimes 2020

lat = 43.67            # Latitude of simulation
tspan = (0.0, 1461.0); # Number of days to be simulated
pop = 2088.012 .* 4;   # Population density per 4km^2

#Introduction of infected humans to the population
function infect(time)
 if mod(time,365) >= 140  && mod(time,365) <= 141  && time > 1095
    impulse = 1
  else
    impulse = 0
  end
  return(impulse)
end;

Cag_DP_Clim = CSV.read(raw"Cag_Clim.csv", header=true, DataFrame) #Loads climate data
##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.

Cag_Temp = Cag_DP_Clim[:,3]
Cag_Prec = Cag_DP_Clim[:,4]
Cag_EV   = Cag_DP_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Cag_Prec[6210:7670])
Temps_spl  = Spline1D(collect(range(1,1461,length = 1461)),Cag_Temp[6210:7670 ])
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)),Cag_EV[6210:7670 ])

include("SO2_public.jl") #Runs SIR model

#Processes outputs
sol_cag20 = sol.t

inf_sum_cag20   = inf_sum
adult_sum_cag20 = adult_sum
H_I_cag20       = out[numb_H_I,:]
H_R_cag20       = out[numb_H_R,:]
H_S_cag20       = out[numb_H_S,:]

#New infections produced per day
Daily_inf_cag20 = (vector_human.(Temp.(sol.t .- 4)) .* bite.(Temp.(sol.t .- 4)))   .*  40 .* 4 .* inf_sum .* out[numb_H_S,:] ./ (out[numb_H_S,:] .+ out[numb_H_I,:] .+ out[numb_H_R,:]  .+ 15000   )

#Makes plots
dt = DateTime(2016,1,1)
     min_date = DateTime(2019,4,1)
     max_date = DateTime(2019,12,31)
    time_axis = add_float2datetime.(sol_cag20,dt)
   time_axis2 = add_float2datetime.(1:1461,dt)

gH_I_cag20 = Plots.plot(time_axis, H_I_cag20  ,label = "Simulated dengue cases", xlabel = "Time (Days)", ylabel = "Dengue Cases")
                plot!(xlims = Dates.value.([min_date,max_date]))
                Plots.plot!(time_axis,H_S_cag20, legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")
                Plots.plot(time_axis, H_R_cag20, legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")

daily_cag20 = plot(time_axis, (Daily_inf_cag20  ) , label ="Simulated dengue cases" , ylabel = "Number of new cases", xlabel = "Time (days)")
                  plot!(xlims = Dates.value.([min_date,max_date]))

##Tokyo 2014

lat   = 35.6762
tspan = (0.0, 5844.0);
pop   = 17254.96 .* 4;

function infect(time)
  if mod(time,365) >= 170  && mod(time,365) <= 185  && time > 5110
     impulse = 2
   else
     impulse = 0
   end
   return(impulse)
 end;
 
Tok_Clim = CSV.read(raw"Tok_Clim.csv", header=true,DataFrame)
##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.

#Tok_Den_14 = CSV.read(raw"Tok_Den_14.csv", header=true,DataFrame)

Tok_Temp = Tok_Clim[:,3]
Tok_Prec = Tok_Clim[:,4]
Tok_EV =   Tok_Clim[:,5]

Precip_spl    = Spline1D(collect(range(1,5844,length = 5844)),Tok_Prec[1:5844] )
   Temps_spl  = Spline1D(collect(range(1,5844,length = 5844)) ,Tok_Temp[1:5844])
   EV_spl     = Spline1D(collect(range(1,5844,length = 5844)) ,Tok_EV[1:5844])
   
include("SO2_public.jl")

solTok14       = sol.t
inf_sum_2014   = inf_sum
adult_sum_2014 = adult_sum
H_I_T_2014     = out[numb_H_I,:]
H_R_T_2014     = out[numb_H_R,:]
H_S_T_2014     = out[numb_H_S,:]

Daily_inf_T_2014= (vector_human.(Temp.(sol.t .- 4)) .* bite.(Temp.(sol.t .- 4)))   .*   40 .* 4 .* inf_sum .* out[numb_H_S,:] ./ (out[numb_H_S,:] .+ out[numb_H_I,:] .+ out[numb_H_R,:]  .+ 15000   )

tick_years = DateTime.(["2014-08-01", "2014-09-01", "2014-10-01"])
DateTick = ["August", "September", "October"]

dt = DateTime(2000,1,1)
        min_date = DateTime(2013,5,1)
        max_date = DateTime(2015,12,31)
       time_axis = add_float2datetime.(solTok14,dt)
      time_axis2 = add_float2datetime.(1:5844,dt)

gH_I_T_2014 = Plots.plot(time_axis, H_I_T_2014,label = "Simulated dengue cases", xlabel = "Time (Days)", ylabel = "New cases per day")
                        plot!(xlims = Dates.value.([min_date,max_date]))
                  Plots.plot!(time_axis, H_S_T_2014, legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")
                  Plots.plot!(time_axis, H_R_T_2014, legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")

daily_T_2014 = plot(time_axis, (Daily_inf_T_2014 ), ylabel = "New cases per day", linewidth = 1, label = "Predicted cases")
                  plot!(xlims = Dates.value.([min_date,max_date]))
                  plot!(add_float2datetime.(Tok_Den_14[:,1] .+ 5110   ,dt),(Tok_Den_14[:,2]),  ylabel = "New cases per day", linewidth = 1, label = "",xticks = false)
                  plot!(ylims = (-0.45,15))
                  plot!(xticks=(tick_years,DateTick), xtickfontsize=7)

#Saint Paul, la Reunion

lat = -21.012870920008957
tspan = (0.0, 2192.0);
pop = 4 .*  4000;

function infect(time)
  if  time > 440 && time < 800
     impulse = 1/12
   else
     impulse = 0
   end
   return(impulse)
 end;
 
SP_Clim   = CSV.read(raw"SB_Clim.csv", header=true,DataFrame)
##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.


#RE_Den_18 = CSV.read(raw"RE_Den_18.csv", header=false,DataFrame)
#RE_Den_19 = CSV.read(raw"RE_Den_1820.csv", header=false,DataFrame)

SP_Temp = SP_Clim[:,3]
SP_Prec = SP_Clim[:,4]
SP_EV   = SP_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,2192 ,length = 2192 )),SP_Prec[2923:5114])
Temps_spl  = Spline1D(collect(range(1,2192 ,length =  2192 )),SP_Temp[2923:5114])
EV_spl     = Spline1D(collect(range(1,2192 ,length = 2192 )),SP_EV[2923:5114])

include("SO2_public.jl")

dt = DateTime(2016,1,1)
     min_date = DateTime(2017,3,1)
     max_date = DateTime(2021,9,1)
    time_axis = add_float2datetime.(1:2192,dt)
   time_axis2 = add_float2datetime.(sol.t,dt)

solsp19        = sol.t
inf_sum_sp19   = inf_sum
adult_sum_sp19 = adult_sum
H_I_sp19       = out[numb_H_I,:]
H_R_sp19       = out[numb_H_R,:]
H_S_sp19       = out[numb_H_S,:]
wing_sp_13     = avg_wing

Daily_inf_sp19= (vector_human.(Temp.(solsp19 .- 4)) .* bite.(Temp.(solsp19 .- 4)))   .*   40 .* 4  .* inf_sum_sp19 .* H_S_sp19 ./ (H_S_sp19 .+ H_I_sp19 .+H_R_sp19 .+ 15000  )

gH_I_RE = Plots.plot(time_axis2, H_I_sp19 ,label = "Simulated dengue cases", xlabel = "Time (Days)", ylabel = "Dengue Cases")
                  plot!(time_axis2, out[numb_H_R,:])

Daily_inf_spl = Spline1D(sol.t, Daily_inf_sp19)

Inf_vec1 = Daily_inf_spl.(1:2192)
Inf_vec = repeat([0.0],2192)

for i = 7:2192
    Inf_vec[i] = sum(Inf_vec1[i-6:i])
end

Inf_spl = movingaverage(Daily_inf_spl.(1:2192),1)

daily_RE = plot(time_axis, (Inf_vec .* 55.3 ),label ="Plastic model" , ylabel = "New cases per week", linewidth = 2, xtickfontsize=12, ytickfontsize=12, yguidefontsize = 16, color = pal[1])
               # plot!(add_float2datetime.(RE_Den_19[:,1]  .+ 730  .- 20  ,dt),(RE_Den_19[:,2]) ,label ="Observed cases", linewidth = 2, color = pal[2])
               # plot!(add_float2datetime.(RE_Den_18[1:24,1]  .+ 365  .-20  ,dt),(RE_Den_18[1:24,2]) , label = nothing, linewidth = 2, color = pal[2])
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims = (-99,3300))

#Hawai'i America

lat = 19.64
tspan = (0.0, 1461.0);
pop = 4 .*  1323;

function infect(time)
  if  time > 950 && time < 1030
       impulse = 0.3
   else
       impulse = 0
   end
   return(impulse)
end;

Kona_Clim   = CSV.read(raw"Kon_Clim.csv", header=true,DataFrame)
##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.

Kona_Temp = Kona_Clim[:,3]
Kona_Prec = Kona_Clim[:,4]
Kona_EV   = Kona_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Kona_Prec)
Temps_spl  = Spline1D(collect(range(1,1461,length =  1461)),Kona_Temp)
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)),Kona_EV)

#Kona_Deng_W = CSV.read(raw"Kon_Den_W.csv", header=false,DataFrame)

include("SO2_public.jl")

dt = DateTime(2013,1,1)
     min_date = DateTime(2015,1,1)
     max_date = DateTime(2017,1,1)
    time_axis = add_float2datetime.(1:1461,dt)
   time_axis2 = add_float2datetime.(sol.t,dt)

solkon        = sol.t
inf_sum_kon   = inf_sum
adult_sum_kon = adult_sum
H_I_kon       = out[numb_H_I,:]
H_R_kon       = out[numb_H_R,:]
H_S_kon       = out[numb_H_S,:]

Daily_inf_kon= (vector_human.(Temp.(sol.t .- 4)) .* bite.(Temp.(sol.t .- 4)))   .*   40 * 4  .* inf_sum .* out[numb_H_S,:] ./ (out[numb_H_S,:] .+ out[numb_H_I,:] .+ out[numb_H_R,:]  .+ 15000  )

H_S_spl = Spline1D(sol.t, out[numb_H_S,:])

new_infections = zeros(1461)

new_infections[1] = 0
for i in 2:1461
    new_infections[i] =   H_S_spl(i-1) - H_S_spl(i)
end

weekly_infection = zeros(1461)

for i in 7:1461
    weekly_infection[i] =   sum(new_infections[i-6:i])
end

plot(time_axis, weekly_infection .* 56,  label ="Predicted cases" , ylabel = "New cases per week", linewidth = 1)
   # plot!(add_float2datetime.(Kona_Deng_W[:,1] .+ 730 .- 20, dt) , Kona_Deng_W[:,2],label = "Observed cases", linewidth = 1)
    plot!(xlims = Dates.value.([min_date,max_date]))