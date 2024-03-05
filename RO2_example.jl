## RO2_example.jl
#Author: Dominic Brass (dombra@ceh.ac.uk)
#Date: 04/03/2024

#Code accompanying the paper, "Phenotypic plasticity in vector traits drives trends in global disease incidence".
#Provides an example of how predictions of mosquito population dynamics are produced for Rimini Italy. Due to the ERA5 license agreement
#the climate files used to produce the paper outputs cannot be shared here but this data is freely available. To perform the comparison between 
#model fit and field data additional files are required commented with where the field data used to produce our figures can be obtained. 

#Requires:
#AO2_public.jl
#Rim_clim.csv

#Loads packages
using(Dates)
using(DifferentialEquations)
using(Dierckx)
using(Plots)
using(CSV)
using(Interpolations)
using(QuadGK)
using(Statistics)
using(DataFrames)

cd(@__DIR__) #Sets working directory to location of file 

#Defines function used to define moving average
function movingaverage(X::Vector,numofele::Int)
    BackDelta = div(numofele,2)
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(X)
    Y = similar(X)
    for n = 1:len
        lo = max(1,n - BackDelta)
        hi = min(len,n + ForwardDelta)
        Y[n] = mean(X[lo:hi])
    end
    return Y
end

#Defines function that coputes R^2 between model predictions and field observations
include("RMSE_fit.jl") 

pal = get_color_palette(:auto, plot_color(:white)) #Colour palette

####Europe#####

#Rimini - Carrieri 2008

lat = 44.0678          # Latitude of location 
tspan = (0.0, 4000.0); # Number of days to be simulated

Rim_Clim = CSV.read(raw"Rim_Clim.csv", header=true,DataFrame)  #Loads climae data for Rimini, Italy

##NOTE this file is not supplied in the repository, to run model download the following environmental variables 
#from the ERA5-Land hourly data from 1950 to present (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form):
#2m temperature
#Evaporation from open water surfaces excluding oceans
#Total precipitation
#The following code expects a CSV of daily average climate variables from 2000-2021 with date in the second column
#temperature in the thrid column, precipitation in the fourth column, and evaporation in the fifth colummn.

Rim_Temp = Rim_Clim[:,3]
Rim_Prec = Rim_Clim[:,4]
Rim_EV   = Rim_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_Prec[1828:6000])
Temps_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_Temp[1828:6000])
EV_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_EV[1828:6000])

#Car_dat_Rim = CSV.read(raw"Car_dat_Rim.csv", header=false,DataFrame) #Loads observed oviposition activity in Rimini Italy taken from Carrieri et al. (2011)
#"Aedes albopictus (Diptera: Culicidae) Population Size Survey in the 2007 Chikungunya Outbreak Area in Italy. I. Characterization of
#Breeding Sites and Evaluation of Sampling Methodologies"

include("AO2_public.jl") #Runs the model

#Processes the model outputs
Ovi_Obs    = repeat([0.0],Int(tspan[2]))

sampling_period = 7

for i = 21:Int(tspan[2])
  if tau_E_daily(i) < sampling_period
    Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:i - tau_E_daily(i)))
  else
    Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:Int(i - sampling_period)))
  end
end

Egg_sum_rim1 = Ovi_Obs
adult_sum_rim = adult_sum
temp_sim_rim = Temp_sim
larv_sum_rim  = larv_sum
solrim = 1:4000
soltrim = sol.t

wing_rim = avg_wing

#Defines times for plotting
dt = DateTime(2005,1,1)
              min_date = DateTime(2008,1,1)
              max_date = DateTime(2010,1,1)
             time_axis = add_float2datetime.(solrim,dt)
            time_axis2 = add_float2datetime.(soltrim,dt)

##########Plots##############

######Climate
Clim  = plot(time_axis2,Temp_sim  , label = nothing, xlabel = "Date", ylabel = "Temperature (°C)")
            plot!(xlims = Dates.value.([min_date,max_date]))

Precip  = plot(time_axis2,Precip_check.(soltrim) ./ surface_area  , label = nothing, xlabel = "Date", ylabel = "Precipiation (mm)")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-0.4, 40))

Photo  = plot(time_axis2,PP.(soltrim)   , label = nothing, xlabel = "Date", ylabel = "Photoperiod (hours)")
            plot!(xlims = Dates.value.([min_date,max_date]))

Food_plot  = plot(time_axis2, Water.(sol.t) , label = nothing, xlabel = "Date", ylabel = "Food availability (mg/day)")
            plot!(xlims = Dates.value.([min_date,max_date]))

#####Eggs
E =  plot(time_axis2,Egg_sum  , label = "Total eggs", xlabel = "Date", ylabel = "Number of Eggs")
    plot!(time_axis2, out[1,:], label = "Active eggs")
    plot!(time_axis2, out[2,:], label = "Diapausing eggs")
    plot!(time_axis2,out[3,:], label = "Quiescent eggs")
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(ylims = (-2.5 * 10^2, 2.5*10^4))

E_tau = Plots.plot(time_axis2, out[numb_t_E,:], legend = false, xlabel = "Date", ylabel = "Duration of active egg stage (Days)")
        plot!(xlims = Dates.value.([min_date,max_date]))

E_surv = Plots.plot(time_axis2, out[numb_P_E,:], legend = false, xlabel = "Date", ylabel = "Through stage survival of active eggs")
        plot!(xlims = Dates.value.([min_date,max_date]))

#####Larvae
L = Plots.plot(time_axis2, out[numb_L_1,:], legend = false, xlabel = "Date", ylabel = "Number of Larvae")
    plot!(xlims = Dates.value.([min_date,max_date]))
L_tau = Plots.plot(time_axis2, out[numb_t_L,:], legend = false, xlabel = "Date", ylabel = "Duration of larval stage (Days)")
    plot!(xlims = Dates.value.([min_date,max_date]))
L_surv = Plots.plot(time_axis2, out[numb_P_L,:], legend = false, xlabel = "Date", ylabel = "Through stage survival of larvae")
    plot!(xlims = Dates.value.([min_date,max_date]))

#####Pupae
P_tau = Plots.plot(time_axis2, out[numb_t_P,:], legend = false, xlabel = "Date", ylabel = "Duration of pupal stage (Days)")
    plot!(xlims = Dates.value.([min_date,max_date]))
P_surv = Plots.plot(time_axis2, out[numb_P_P,:], legend = false, xlabel = "Date", ylabel = "Through stage survival of pupae")
    plot!(xlims = Dates.value.([min_date,max_date]))

#####Adults
A_plot = plot(time_axis2, out[numb_A_1:numb_A_2,:]', label = nothing, color = :inferno, line_z = Wing_Vals', xlabel = "Date", ylabel = "Number of Adults", colorbar_title = "Wing length (mm)")
      Plots.plot!(time_axis2,movingaverage(adult_sum ,1), label = "Total Adults", color = :black, xlabel = "Date", ylabel = "Number of Adults")
      plot!(xlims = Dates.value.([min_date,max_date]))

######Compares model prediction to field observations, requires Car_dat_Rim.csv.      
#rim_fit1 = RMSE_fun(movingaverage(Egg_sum_rim1,7) ,solrim,Car_Tot_1,7, 1095, 0.01,0.0001)

#grim1 =plot(time_axis,movingaverage(Egg_sum_rim1,7) .* rim_fit1[3]  , label = "Simulated oviposition activity", ylabel = "Oviposition activity")
#    plot!(add_float2datetime.(Car_Tot_1[:,1] .+ 1095 .+ rim_fit1[2],dt) , Car_Tot_1[:,2]       ,   label = "Observed oviposition activity" )
#    plot!(add_float2datetime.(Car_Tot_1[:,1] .+ 1095 .+ rim_fit1[2],dt)  , Car_Tot_1[:,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
#        plot!(xlims = Dates.value.([min_date,max_date]))
