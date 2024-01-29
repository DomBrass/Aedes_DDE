
using(Dates)
using(DifferentialEquations)
using(Dierckx)
using(Plots)
using(CSV)
using(Interpolations)
using(QuadGK)
using(Statistics)
using(DataFrames)

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

include("RMSE_fit.jl")
pal = get_color_palette(:auto, plot_color(:white))

####Europe#####

##Germany


# ##Freiburg DONE

lat = 47.9987
alt  = 3
tspan = (0.0, 1461.0);                    # Number of days to be simulated

Metz_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Frei_Clim.csv", header=true,DataFrame)

Metz_Temp = Metz_DP_Clim[:,3]
Metz_Prec = Metz_DP_Clim[:,4]
Metz_EV   = Metz_DP_Clim[:,5]
Metz_DP   = Metz_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461 ,length = 1461 )),Metz_Prec[6211:7671])
Temps_spl = Spline1D(collect(range(1,1461 ,length = 1461 )),Metz_Temp[6211:7671])
DP_spl = Spline1D(collect(range(1,1461 ,length = 1461 )),Metz_DP[6211:7671])
EV_spl = Spline1D(collect(range(1,1461 ,length = 1461 )),Metz_EV[6211:7671])
Metz_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Metz_Eggs.csv", header=false,DataFrame)

RMSE_fun( Temps_spl.(1:1461),Metz_DP_Clim[6211:7671,1],[Metz_DP_Clim[6211:7671,1] Metz_Temp[6211:7671]],0, 0, 1,1)
RMSE_fun( Precip_spl.(1:1461),Metz_DP_Clim[6211:7671,1],[Metz_DP_Clim[6211:7671,1] Metz_Prec[6211:7671]],0, 0, 1,1)

plot(Metz_DP_Clim[6211:7671,1],Metz_Temp[6211:7671]) 
  plot!(Metz_DP_Clim[6211:7671,1],Temps_spl.(1:1461))

include("AO2.jl")

Egg_sum_Metz = Ovi_Obs
temp_Metz = Temp_sim
solMetz = 1:1461

dt = DateTime(2017,1,1)
              min_date = DateTime(2020,6,1)
              max_date = DateTime(2020,11,1)
              time_axis = add_float2datetime.(solMetz,dt)

Met_E = [Metz_Egg[3:6,1]  Metz_Egg[3:6,7]]

met_fit = RMSE_fun(movingaverage(Egg_sum_Metz,1) ,solMetz,Met_E,7, 1095, 0.2,0.001)

gMetz =plot(time_axis,movingaverage(Egg_sum_Metz,7) .* met_fit[3]  ,  label = "Simulated oviposition activity", ylabel = "Oviposition activity")
        plot!(add_float2datetime.(Metz_Egg[:,1] .+ 1095 .+met_fit[2] ,dt) , Metz_Egg[:,7]       ,   label = "Oviposition in non-intervention area" )
        plot!(add_float2datetime.(Metz_Egg[:,1] .+ 1095  .+met_fit[2] ,dt)  , Metz_Egg[:,7]       ,   label = nothing, seriestype =:scatter, color=:orange)
        plot!(add_float2datetime.(Metz_Egg[:,1] .+ 1095 .+ met_fit[2]  ,dt) , Metz_Egg[:,4]       ,   label = "Oviposition in intervention area" )
        plot!(add_float2datetime.(Metz_Egg[:,1] .+ 1095 .+met_fit[2]  ,dt)  , Metz_Egg[:,4]       ,   label = nothing, seriestype =:scatter, color=:purple)
                  plot!(xlims = Dates.value.([min_date,max_date]))

##Ludwigschafen DONE

lat = 49.48776449468759
alt  = 3
tspan = (0.0, 1461.0);                    # Number of days to be simulated

Lud_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Lud_Clim.csv", header=true,DataFrame)

Lud_Temp = Lud_DP_Clim[:,3]
Lud_Prec = Lud_DP_Clim[:,4]
Lud_EV   = Lud_DP_Clim[:,5]
Lud_DP   = Lud_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Lud_Prec[6211:7671])
Temps_spl = Spline1D(collect(range(1,1461,length = 1461)),Lud_Temp[6211:7671])
DP_spl = Spline1D(collect(range(1,1461,length = 1461)),Lud_DP[6211:7671])
EV_spl = Spline1D(collect(range(1,1461,length = 1461)),Lud_EV[6211:7671])
Lud_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Lud_Egg.csv", header=false,DataFrame)

include("AO2.jl")

Egg_sum_Lud = Ovi_Obs
temp_Lud = Temp_sim
solLud = 1:1461

dt = DateTime(2017,1,1)
            min_date = DateTime(2020,6,1)
            max_date = DateTime(2020,11,1)
            time_axis = add_float2datetime.(solLud,dt)

lud_fit = RMSE_fun(movingaverage(Egg_sum_Lud,7) ,solLud,Lud_Egg,7, 1095, 0.1,0.001)


gLud =plot(time_axis,movingaverage(Egg_sum_Lud,1) .* lud_fit[3]  ,  label = "Simulated oviposition activity", ylabel = "Oviposition activity")
      plot!(add_float2datetime.(Lud_Egg[:,1] .+ 1088 .+ lud_fit[2] ,dt) , Lud_Egg[:,2]       ,   label = "Observed oviposition activity" )
      plot!(add_float2datetime.(Lud_Egg[:,1] .+ 1088 .+ lud_fit[2]  ,dt)  , Lud_Egg[:,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
                plot!(xlims = Dates.value.([min_date,max_date]))



# ##Budva DONE

lat = 42.288432232332674
alt  = 3
tspan = (0.0, 2000.0);                    # Number of days to be simulated

Bud_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Bud_Clim.csv", header=true,DataFrame)

Bud_Temp = Bud_DP_Clim[:,3]
Bud_Prec = Bud_DP_Clim[:,4]
Bud_EV   = Bud_DP_Clim[:,5]
Bud_DP   = Bud_DP_Clim[:,6]


Precip_spl = Spline1D(collect(range(1,4200,length = 4200)),Bud_Prec[3654:7853])
Temps_spl = Spline1D(collect(range(1,4200,length = 4200)),Bud_Temp[3654:7853])
DP_spl = Spline1D(collect(range(1,4200,length = 4200)),Bud_DP[3654:7853])
EV_spl = Spline1D(collect(range(1,4200,length = 4200)),Bud_EV[3654:7853])

Bud_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Bud_O1.csv", header=false,DataFrame;  dateformat="dd/mm/yyyy")

include("AO2.jl")

Egg_sum_bud = Ovi_Obs
temp_bud = Temp_sim
solbud = 1:2000

dt = DateTime(2010,1,1)
              min_date = DateTime(2012,4,1)
              max_date = DateTime(2012,12,1)
              time_axis = add_float2datetime.(solbud,dt)

Temp_bud = zeros(18)

for i in 1:18
Temp_bud[i] = (Bud_Egg[i,1] .-  Date(2010,1,1)).value
end

Temp_bud1 = [Temp_bud[1:18,1] Bud_Egg[1:18,2]]

bud_fit = RMSE_fun(movingaverage(Egg_sum_bud,7) ,solbud,Temp_bud1,7, 0, 0.2,0.001)

gbud =plot(time_axis,movingaverage(Egg_sum_bud,7) .* bud_fit[3] ,  label = "Simulated oviposition activity", ylabel = "Oviposition activity")
                plot!(add_float2datetime.(Temp_bud1[:,1] .+ bud_fit[3],dt) , Bud_Egg[:,2]       ,   label = "Observed oviposition activity" )
                plot!(add_float2datetime.(Temp_bud1[:,1]  .+ bud_fit[3] ,dt)  , Bud_Egg[:,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
                  plot!(xlims = Dates.value.([min_date,max_date]))
                  plot!(ylims = (-3.6,360))


# ##Podgarica DONE
#
lat = 42.46024640300247
alt  = 40
tspan = (0.0, 2000.0);                    # Number of days to be simulated

Pod_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Pod_Clim.csv", header=true,DataFrame)

Pod_Temp = Pod_DP_Clim[:,3]
Pod_Prec = Pod_DP_Clim[:,4]
Pod_EV   = Pod_DP_Clim[:,5]
Pod_DP   = Pod_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,4200,length = 4200)),Pod_Prec[3654:7853])
Temps_spl = Spline1D(collect(range(1,4200,length = 4200)),Pod_Temp[3654:7853])
DP_spl = Spline1D(collect(range(1,4200,length = 4200)),Pod_DP[3654:7853])
EV_spl = Spline1D(collect(range(1,4200,length = 4200)),Pod_EV[3654:7853])

Pod_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Brad_Egg.csv", header=true,DataFrame)

include("AO2.jl")

Egg_sum_pod= Ovi_Obs

temp_pod= Temp_sim
solpod = sol.t
solpod2 = 1:2000

dt = DateTime(2010,1,1)
            min_date = DateTime(2013,4,1)
            max_date = DateTime(2014,1,1)
            time_axis = add_float2datetime.(solpod2,dt)

Pod_Egg = [Pod_Egg[19:37,1]      Pod_Egg[19:37,4]]

pod_fit = RMSE_fun(movingaverage(Egg_sum_pod ,7),1:2000,Pod_Egg,14, 0,0.1,0.0001)


gpod =plot(time_axis,movingaverage(Egg_sum_pod,7)  .* pod_fit[3] , ylabel = "Oviposition activity", label = "Predicted oviposition activity")
              plot!(add_float2datetime.(Pod_Egg[:,1] .+ 730 .+ pod_fit[2] ,dt) , Pod_Egg[:,2]       ,   label = "Observed oviposition activity" )
              plot!(add_float2datetime.(Pod_Egg[:,1]  .+ 730 .+ pod_fit[2]   ,dt)  , Pod_Egg[:,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
                plot!(xlims = Dates.value.([min_date,max_date]))

# ##Zambelici DONE
#
lat = 42.42917290324098
alt  = 40
tspan = (0.0, 2009.0);                    # Number of days to be simulated

Zam_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Zam_Clim.csv", header=true,DataFrame)

Zam_Temp = Zam_DP_Clim[:,3]
Zam_Prec = Zam_DP_Clim[:,4]
Zam_EV   = Zam_DP_Clim[:,5]
Zam_DP   = Zam_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2009,length = 2009)),Zam_Prec[5845:7853])
Temps_spl = Spline1D(collect(range(1,2009,length = 2009)),Zam_Temp[5845:7853])
DP_spl = Spline1D(collect(range(1,2009,length = 2009)),Zam_DP[5845:7853])
EV_spl = Spline1D(collect(range(1,2009,length = 2009)),Zam_EV[5845:7853])


Zam_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Zam_O1.csv", header=false,DataFrame; dateformat ="dd/mm/yyyy")

include("AO2.jl")

Egg_sum_zam= Ovi_Obs
temp_zam= Temp_sim
solzam = 1:2009

dt = DateTime(2016,1,1)
            min_date = DateTime(2018,1,1)
            max_date = DateTime(2020,1,1)
            time_axis = add_float2datetime.(solzam,dt)

Temp_zam = zeros(54)

for i in 1:54
    Temp_zam[i] = (Zam_Egg[i,1] .-  Date(2016,1,1)).value
end

Temp_zam1 = [Temp_zam[1:54,1] Zam_Egg[1:54,2] ./ 15]

zam_fit = RMSE_fun(movingaverage(Egg_sum_zam ,16),1:2009,Temp_zam1,14, 0,0.1,0.001)

gzam =plot(time_axis,movingaverage(Egg_sum_zam,7)  .* zam_fit[3]   , ylabel = "Oviposition activity", label = "Predicted oviposition activity")
              plot!(add_float2datetime.(Temp_zam1[:,1] .+ zam_fit[2] ,dt) , Zam_Egg[1:54,2]  ./ 15    ,   label = "Observed oviposition activity" )
              plot!(add_float2datetime.(Temp_zam1[:,1] .+ zam_fit[2] ,dt)  , Zam_Egg[1:54,2]  ./ 15    ,   label = nothing, seriestype =:scatter, color=:orange)
                plot!(xlims = Dates.value.([min_date,max_date]))

###Serbia

# #Novi Sad DONE

lat = 45.250315
alt  = 80
tspan = (0.0, 2071.0);                    # Number of days to be simulated

Ser_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Nov_Clim.csv", header=true,DataFrame)

Ser_Temp = Ser_DP_Clim[:,3]
Ser_Prec = Ser_DP_Clim[:,4]
Ser_EV   = Ser_DP_Clim[:,5]
Ser_DP   = Ser_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2071,length = 2071)),Ser_Prec[5845:7915])
Temps_spl = Spline1D(collect(range(1,2071,length = 2071)),Ser_Temp[5845:7915])
DP_spl = Spline1D(collect(range(1,2071,length = 2071)),Ser_DP[5845:7915])
EV_spl = Spline1D(collect(range(1,2071,length = 2071)),Ser_EV[5845:7915])

Ser_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Nov_S_O.csv", header=false,DataFrame;  dateformat="dd/mm/yyyy")

include("AO2.jl")

Egg_sum_ser = Ovi_Obs
temp_ser = Temp_sim
solser = 1:2071

dt = DateTime(2016,1,1)
              min_date = DateTime(2020,5,1)
              max_date = DateTime(2021,9,1)
              time_axis = add_float2datetime.(solser,dt)

Temp_ser = zeros(30)

for i in 1:30
    Temp_ser[i] = (Ser_Egg[i,1] .-  Date(2016,1,1)).value
end

Temp_ser1 = [Temp_ser[1:30,1] Ser_Egg[1:30,2]]

ser_fit = RMSE_fun(movingaverage(Egg_sum_ser ,7),1:2071,Temp_ser1[1:27,:],14, 0,0.1,0.0001)

gser =plot(time_axis,movingaverage(Egg_sum_ser,7) .* ser_fit[3] , ylabel = "Oviposition activity", label = "Predicted oviposition activity")
  plot!(add_float2datetime.(Temp_ser1[1:27,1]  .+ ser_fit[2] ,dt) , Ser_Egg[1:27,2]       ,   label = "Observed oviposition activity" )
  plot!(add_float2datetime.(Temp_ser1[1:27,1]  .+ ser_fit[2],dt)  , Ser_Egg[1:27,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(ylims = (-4,400))

##Loznica DONE

lat = 44.5338
alt  = 121
tspan = (0.0, 2071.0);                    # Number of days to be simulated

Loz_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Loz_Clim.csv", header=true,DataFrame)

Loz_Temp = Loz_DP_Clim[:,3]
Loz_Prec = Loz_DP_Clim[:,4]
Loz_EV   = Loz_DP_Clim[:,5]
Loz_DP   = Loz_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2009,length = 2009)),Loz_Prec[5845:7853])
Temps_spl = Spline1D(collect(range(1,2009,length = 2009)),Loz_Temp[5845:7853])
DP_spl = Spline1D(collect(range(1,2009,length = 2009)),Loz_DP[5845:7853])
EV_spl = Spline1D(collect(range(1,2009,length = 2009)),Loz_EV[5845:7853])

Loz_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Loz_O3.csv", header=false,DataFrame;  dateformat="dd/mm/yyyy")
Loz_A = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Loz_A3.csv", header=false,DataFrame;  dateformat="dd/mm/yyyy")

include("AO2.jl")

Egg_sum_loz = Ovi_Obs
temp_loz = Temp_sim
solloz = sol.t
solloz2 = 1:2071
adult_sum_L = adult_sum

dt = DateTime(2016,1,1)
              min_date = DateTime(2019,5,1)
              max_date = DateTime(2021,12,1)
              time_axis = add_float2datetime.(solloz,dt)
              time_axis2 = add_float2datetime.(solloz2,dt)

Temp_loz = zeros(20)

for i in 1:20
    Temp_loz[i] = (Loz_Egg[i,1] .-  Date(2016,1,1)).value
end

Temp_loz1 = [Temp_loz[1:20,1] Loz_Egg[1:20,2]]

Temp_lozA = zeros(21)

for i in 1:21
    Temp_lozA[i] = (Loz_A[i,1] .-  Date(2016,1,1)).value
end

Temp_lozA1 = [Temp_lozA[1:21,1] Loz_A[1:21,4]]


loz_fit = RMSE_fun(movingaverage(Egg_sum_loz ,14),1:2071,Temp_loz1[1:10,:],14, 0,2,0.01)
loz_fitA = RMSE_fun(movingaverage(adult_sum_L ,14),solloz,Temp_lozA1,2, 0,0.5,0.01)

gloz =plot(time_axis2,movingaverage(Egg_sum_loz,7) .* loz_fit[3]  , ylabel = "Oviposition activity", label = "Predicted oviposition activity")
    plot!(add_float2datetime.(Temp_loz1[1:10,1] .+ loz_fit[2]   ,dt) , Loz_Egg[1:10,2]       ,   label = "Observed oviposition activity" )
    plot!(add_float2datetime.(Temp_loz1[1:10,1] .+ loz_fit[2] ,dt)  , Loz_Egg[1:10,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
    plot!(add_float2datetime.(Temp_loz1[11:20,1] .+ loz_fit[2]   ,dt) , Loz_Egg[11:20,2]       ,   label = nothing , color=palette(:default)[2])
    plot!(add_float2datetime.(Temp_loz1[11:20,1] .+ loz_fit[2] ,dt)  , Loz_Egg[11:20,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(ylims = (-5,1200))
    vline!(add_float2datetime.([182,189,292,222,249,267,536,546,572,590,600].+1460,dt), color=palette(:default)[3], line=:dash , label = "Control", linewidth = 0.5)

gloz_A =plot(time_axis,movingaverage(adult_sum_L,10) .* loz_fitA[3]  , label = "Simulated oviposition activity", ylabel = "Predicted oviposition activity")
          plot!(add_float2datetime.(Temp_lozA1[1:21,1]     ,dt) , Loz_A[1:21,4]       ,   label = "Observed oviposition activity" )
        plot!(add_float2datetime.(Temp_lozA1[1:21,1]  ,dt)  , Loz_A[1:21,4]       ,   label = nothing, seriestype =:scatter, color=:orange)
        vline!(add_float2datetime.([182,189,292,222,249,267,536,546,572,590,600].+1460,dt), color=palette(:default)[3], line=:dash , label = "Control", linewidth = 0.5)
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylims = (-0.5,60))

#Italy###

##Catiania - Bella 2018 DONE

lat = 37.5079
alt  = 5
tspan = (0.0, 4000.0);                    # Number of days to be simulated

Cat_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Cat_Clim.csv", header=true,DataFrame)
Cat_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Cat_Egg.csv", header=true,DataFrame)

Cat_Temp = Cat_Clim[:,3]
Cat_Prec = Cat_Clim[:,4]
Cat_EV   = Cat_Clim[:,5]
Cat_DP   = Cat_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,4173,length = 4173)),Cat_Prec[1828:6000])
Temps_spl = Spline1D(collect(range(1,4173,length = 4173)),Cat_Temp[1828:6000])
DP_spl = Spline1D(collect(range(1,4173,length = 4173)),Cat_DP[1828:6000])
EV_spl = Spline1D(collect(range(1,4173,length = 4173)),Cat_EV[1828:6000])

include("AO2.jl")

Egg_sum_cat = Ovi_Obs
temp_cat = Temp_sim
solcat = 1:4000

dt = DateTime(2005,1,1)
              min_date = DateTime(2008,1,1)
              max_date = DateTime(2013,12,1)
              time_axis = add_float2datetime.(solcat,dt)

cat_fit1 = RMSE_fun(movingaverage(Egg_sum_cat,14) ,solcat,Cat_Egg[9:17,:],30, 1095, 0.4,0.001)


gcat =plot(time_axis,movingaverage(Egg_sum_cat,14) .* cat_fit1[3] , label = "Simulated oviposition activity", ylabel = "Oviposition activity")
  plot!(add_float2datetime.(Cat_Egg[1:8,1] .+ 1095  .+ cat_fit1[2]   ,dt) , Cat_Egg[1:8,2]       ,   label = "Observed oviposition activity" )
  plot!(add_float2datetime.(Cat_Egg[1:8,1] .+ 1095  .+ cat_fit1[2],dt)  , Cat_Egg[1:8,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
  plot!(add_float2datetime.(Cat_Egg[9:17,1] .+ 1095  .+ cat_fit1[2]  ,dt) , Cat_Egg[9:17,2]       ,   label = nothing, color=pal[2])
  plot!(add_float2datetime.(Cat_Egg[9:17,1] .+ 1095 .+ cat_fit1[2],dt)  , Cat_Egg[9:17,2]       ,   label = nothing, color=:orange, seriestype =:scatter)
  plot!(xlims = Dates.value.([min_date,max_date]))
  plot!(ylims = (-0.5,4800))


#Rimini - Carrieri 2008, Carrieri 2017## DONE

lat = 44.0678
alt  = 5
tspan = (0.0, 4000.0);                    # Number of days to be simulated

Rim_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Rim_Clim.csv", header=true,DataFrame)

Rim_Temp = Rim_DP_Clim[:,3]
Rim_Prec = Rim_DP_Clim[:,4]
Rim_EV   = Rim_DP_Clim[:,5]
Rim_DP   = Rim_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_Prec[1828:6000])
Temps_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_Temp[1828:6000])
DP_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_DP[1828:6000])
EV_spl = Spline1D(collect(range(1,4173,length = 4173)),Rim_EV[1828:6000])

Rim_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Rim_Egg_1.csv", header=false,DataFrame)
Car_Tot_1 = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Car_Tot_1.csv", header=false,DataFrame)
Car_dat_Rim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Car_dat_Rim.csv", header=false,DataFrame)

include("AO2.jl")

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

dt = DateTime(2005,1,1)
              min_date = DateTime(2008,1,1)
              max_date = DateTime(2010,1,1)
              time_axis = add_float2datetime.(solrim,dt)
              time_axis2 = add_float2datetime.(soltrim,dt)

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


L = Plots.plot(time_axis2, out[numb_L_1,:], legend = false, xlabel = "Date", ylabel = "Number of Larvae")
    plot!(xlims = Dates.value.([min_date,max_date]))
L_tau = Plots.plot(time_axis2, out[numb_t_L,:], legend = false, xlabel = "Date", ylabel = "Duration of larval stage (Days)")
    plot!(xlims = Dates.value.([min_date,max_date]))
L_surv = Plots.plot(time_axis2, out[numb_P_L,:], legend = false, xlabel = "Date", ylabel = "Through stage survival of larvae")
    plot!(xlims = Dates.value.([min_date,max_date]))
P_tau = Plots.plot(time_axis2, out[numb_t_P,:], legend = false, xlabel = "Date", ylabel = "Duration of pupal stage (Days)")
    plot!(xlims = Dates.value.([min_date,max_date]))
P_surv = Plots.plot(time_axis2, out[numb_P_P,:], legend = false, xlabel = "Date", ylabel = "Through stage survival of pupae")
    plot!(xlims = Dates.value.([min_date,max_date]))


Clim  = plot(time_axis2,Temp_sim  , label = nothing, xlabel = "Date", ylabel = "Temperature (°C)")
            plot!(xlims = Dates.value.([min_date,max_date]))

Precip  = plot(time_axis2,Precip_check.(soltrim) ./ surface_area  , label = nothing, xlabel = "Date", ylabel = "Precipiation (mm)")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-0.4, 40))

Photo  = plot(time_axis2,PP.(soltrim)   , label = nothing, xlabel = "Date", ylabel = "Photoperiod (hours)")
            plot!(xlims = Dates.value.([min_date,max_date]))

Food_plot  = plot(time_axis2, Water.(sol.t) , label = nothing, xlabel = "Date", ylabel = "Food availability (mg/day)")
            plot!(xlims = Dates.value.([min_date,max_date]))

A_plot = plot(time_axis2, out[numb_A_1:numb_A_2,:]', label = nothing, color = :inferno, line_z = Wing_Vals', xlabel = "Date", ylabel = "Number of Adults", colorbar_title = "Wing length (mm)")
      Plots.plot!(time_axis2,movingaverage(adult_sum ,1), label = "Total Adults", color = :black, xlabel = "Date", ylabel = "Number of Adults")
      plot!(xlims = Dates.value.([min_date,max_date]))

I_plot = plot(time_axis2, out[numb_I_1:numb_I_2,:]', label = nothing, color = :inferno, line_z = Wing_Vals', xlabel = "Date", ylabel = "Number of infected adults", colorbar_title = "Wing length (mm)")
      Plots.plot!(time_axis2,movingaverage(inf_sum ,1), label = "Total Adults", color = :black, xlabel = "Date", ylabel = "Number of infected adults")
      plot!(xlims = Dates.value.([min_date,max_date]))

Mort_Comp  = plot(time_axis2, 1 ./ movingaverage(avg_mort,60), label = "Model prediction", ylabel = "Adult longeivity (Days)")
                plot!(time_axis2, 1 ./ movingaverage(AL(Temp_sim),60), label = "Laboratory model")
                plot!(time_axis2, 1 ./  movingaverage(AM(Temp_sim),60), label = "Field model")
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims =(-4.5,150))
          

#plot(time_axis2, out[numb_H_I,:], xlabel = "Date", ylabel = "Number of people", label = "Infected")
#                plot!(time_axis2, out[numb_H_S,:], label = "Susceptible")
#                plot!(time_axis2, out[numb_H_R,:], label = "Recovered")
#                plot!(xlims = Dates.value.([min_date,max_date]))
#                plot!(ylims = (-100,12000))


#rim_fit1 = RMSE_fun(Egg_sum_rim ,solrim,Car_Tot_1,2, 1095, 0.2,0.01)
rim_fit1 = RMSE_fun(movingaverage(Egg_sum_rim1,7) ,solrim,Car_Tot_1,7, 1095, 0.01,0.0001)

grim1 =plot(time_axis,movingaverage(Egg_sum_rim1,7) .* rim_fit1[3]  , label = "Simulated oviposition activity", ylabel = "Oviposition activity")
    plot!(add_float2datetime.(Car_Tot_1[:,1] .+ 1095 .+ rim_fit1[2],dt) , Car_Tot_1[:,2]       ,   label = "Observed oviposition activity" )
    plot!(add_float2datetime.(Car_Tot_1[:,1] .+ 1095 .+ rim_fit1[2],dt)  , Car_Tot_1[:,2]       ,   label = nothing, seriestype =:scatter, color=:orange)
        plot!(xlims = Dates.value.([min_date,max_date]))

min_date=DateTime(2014,1,1)
max_date =  DateTime(2016,1,1)
 
Ovi_Obs    = repeat([0.0],Int(tspan[2]))

sampling_period = 14

for i = 21:Int(tspan[2])
  if tau_E_daily(i) < sampling_period
    Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:i - tau_E_daily(i)))
  else
    Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:Int(i - sampling_period)))
  end
end


Egg_sum_rim2 = Ovi_Obs

Rim_1415 = [Rim_Egg[:,1]    Rim_Egg[:,2]
            Car_dat_Rim[:,1] .+ 365  Car_dat_Rim[:,2]  ]

rim_fit2 = RMSE_fun(movingaverage(Egg_sum_rim2 ,7),solrim,Rim_1415,14, 3650 - 365 , 0.1,0.001)
#rim_fit3 = RMSE_fun(movingaverage(Egg_sum_rim2 ,14),solrim,Car_dat_Rim,7, 4015 .- 365 ,4,0.01)

grim2 =plot(time_axis,movingaverage(Egg_sum_rim2 ,7) .* rim_fit2[3] , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
  plot!(add_float2datetime.(Rim_Egg[:,1] .+ 3650 .- 365 .+ rim_fit2[2],dt),Rim_Egg[:,2], label = "Observed oviposition activity" , color=:orange  )
  plot!(add_float2datetime.(Rim_Egg[:,1] .+ 3650 .- 365 .+   rim_fit2[2],dt),Rim_Egg[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
  plot!(add_float2datetime.(Car_dat_Rim[:,1] .+ 4015 .- 365  .+ rim_fit2[2],dt ) , Car_dat_Rim[:,2]     ,   label = nothing, color=:orange )
  plot!(add_float2datetime.(Car_dat_Rim[:,1] .+ 4015 .- 365  .+ rim_fit2[2],dt ) , Car_dat_Rim[:,2] ,   label = nothing, seriestype =:scatter, color=:orange )
  plot!(xlims = Dates.value.([min_date,max_date]))
  plot!(ylims = (-1,1500))


min_date = DateTime(2008,4,1)
max_date = DateTime(2008,11,1)

rim_nuis = [-60   0
            0   0.54
            30  1.22
            60   1.4]

grimbite = plot(time_axis2, movingaverage(adult_sum_rim,100) .* bite.(temp_sim_rim),right_margin = 10Plots.mm, color = pal[1], label =nothing, legend=:topleft, xlabel = "Time", ylabel = "Nuisance factor")
            plot!(1, repeat([NaN],1), color = pal[1], label = "Simluated nuisance factor")
            plot!(1, repeat([NaN],1), color = pal[2], label = "Observed bites per day")
            plot!(xlims = Dates.value.([min_date,max_date]))
            subplot = twinx()
            plot!(subplot,add_float2datetime.(rim_nuis[2:4,1].+5*31 .+ 1095 ,dt),rim_nuis[2:4,2], label = nothing, color=pal[2], ylabel = "Bites per day")
            plot!(subplot,add_float2datetime.(rim_nuis[2:4,1].+5*31 .+ 1095 ,dt),rim_nuis[2:4,2] , label = nothing, seriestype =:scatter, color=:orange)
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(subplot,ylims = (-0.05,1.8))


# Rome Manica/Toma## DONE

lat = 41.9028
tspan = (0.0, 5000.0);                    # Number of days to be simulated

Rom_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Rom_Clim.csv", header=true, DataFrame)

Rom_Temp = Rom_DP_Clim[:,3]
Rom_Prec = Rom_DP_Clim[:,4]
Rom_EV   = Rom_DP_Clim[:,5]
Rom_DP   = Rom_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,5000,length = 5000)),Rom_Prec[1:5000])
Temps_spl = Spline1D(collect(range(1,5000,length = 5000)),Rom_Temp[1:5000])
DP_spl = Spline1D(collect(range(1,5000,length = 5000)),Rom_DP[1:5000])
EV_spl = Spline1D(collect(range(1,5000,length = 5000)),Rom_EV[1:5000])

Rome_Adults= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Rome_Ad.csv", header=false, DataFrame)
Rome_Eggs_2000= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Rome_Adults_2000.csv", header=true, DataFrame)

include("AO2.jl")

Egg_sum_rom = Ovi_Obs
adult_sum_rom = adult_sum
solrom = sol.t

rom_fit1 = RMSE_fun(movingaverage(Egg_sum_rom ,7),1:5000,Rome_Eggs_2000,8, 0, 0.1,0.0001)

dt = DateTime(2000,1,1)
      min_date = DateTime(2000,4,1)
      max_date = DateTime(2001,1,1)
      time_axis = add_float2datetime.(solrom,dt)
      time_axis2 = add_float2datetime.(1:5000,dt)

grom1 = plot(time_axis2,movingaverage(Egg_sum_rom ,7) .* rom_fit1[3] , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
      plot!(xlims = Dates.value.([min_date,max_date]))
      plot!(add_float2datetime.(Rome_Eggs_2000[:,1] .+ rom_fit1[2],dt ) , Rome_Eggs_2000[:,2]     ,   label = "Predicted oviposition activity" )
      plot!(add_float2datetime.(Rome_Eggs_2000[:,1] .+ rom_fit1[2],dt ) , Rome_Eggs_2000[:,2] ,   label = nothing, seriestype =:scatter, color=:orange )
      plot!(ylims= (0,180))

rom_fit2 = RMSE_fun(movingaverage(adult_sum_rom ,100),solrom,Rome_Adults,7, 4380, 0.1,0.001)

min_date = DateTime(2012,4,1)
max_date = DateTime(2012,12,31)

grom2 = plot(time_axis,movingaverage(adult_sum_rom ,1) .* rom_fit2[3], label = "Simulated adults", xlabel = "Time (Days)")
               plot!(xlims = Dates.value.([min_date,max_date]))
               plot!(add_float2datetime.(Rome_Adults[:,1] .+ 4380 .+ rom_fit2[2],dt ) , Rome_Adults[:,2]     ,   label = "Observed adults" )
               plot!(add_float2datetime.(Rome_Adults[:,1] .+ 4380 .+ rom_fit2[2],dt ) , Rome_Adults[:,2] ,   label = nothing, seriestype =:scatter, color=:orange )


# #Cosenza - Bonnacci et al. 2015## DONE
#
lat = 39.2983
alt = 238

tspan = (0.0, 1501.0);                    # Number of days to be simulated

Cose_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Cos_Clim.csv", header=true, DataFrame)

Cos_Temp = Cose_DP_Clim[:,3]
Cos_Prec = Cose_DP_Clim[:,4]
Cos_EV   = Cose_DP_Clim[:,5]
Cos_DP   = Cose_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1500,length = 1500)),Cos_Prec[3654:5153])
Temps_spl = Spline1D(collect(range(1,1500,length = 1500)),Cos_Temp[3654:5153])
DP_spl = Spline1D(collect(range(1,1500,length = 1500)),Cos_DP[3654:5153])
EV_spl = Spline1D(collect(range(1,1500,length = 1500)),Cos_EV[3654:5153])

Bon_Eggs =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Bon_O.csv", header=false, DataFrame)

Bon_Oggs = [Bon_Eggs[:,1]   Bon_Eggs[:,2] ./ 20]

include("AO2.jl")

Egg_sum_cos = Ovi_Obs
solcos = 1:1501

dt = DateTime(2010,1,1)
    min_date = DateTime(2013,1,1)
    max_date = DateTime(2014,3,5)
    time_axis = add_float2datetime.(solcos,dt)

cos_fit1 = RMSE_fun(movingaverage(Egg_sum_cos ,7),solcos,Bon_Oggs,40, 1200, 0.1,0.001)

gcos = plot(time_axis,movingaverage(Egg_sum_cos   ,7) .* cos_fit1[3]    , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
  plot!(add_float2datetime.(Bon_Oggs[:,1] .+ 1200 .+ cos_fit1[2]   ,dt), Bon_Oggs[:,2] ,   label = "Observed oviposition activity" )
  plot!(add_float2datetime.(Bon_Oggs[:,1] .+ 1200  .+ cos_fit1[2],dt), Bon_Oggs[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
  plot!(xlims = Dates.value.([min_date,max_date]))


##Como -- Suter## Done

lat = 43.8081
alt = 201

tspan = (0.0, 1826.0);                    # Number of days to be simulated

Como_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Com_Clim.csv", header=true, DataFrame)

Com_Temp = Como_DP_Clim[:,3]
Com_Prec = Como_DP_Clim[:,4]
Com_EV   = Como_DP_Clim[:,5]
Com_DP   = Como_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1826,length = 1826)),Com_Prec[3289:5114])
Temps_spl = Spline1D(collect(range(1,1826,length = 1826)),Com_Temp[3289:5114])
DP_spl = Spline1D(collect(range(1,1826,length = 1826)),Com_DP[3289:5114])
EV_spl = Spline1D(collect(range(1,1826,length = 1826)),Com_EV[3289:5114])

Como_Ovi= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Como_Ovi.csv", header=true, DataFrame;  dateformat="dd/mm/yyyy")
Como_Test= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Com_test.csv", header=false, DataFrame)

include("AO2.jl")

solcom = sol.t
solcom2 = 1:1826
Ovi_com = Ovi_Obs

dt = DateTime(2009,1,1)
    min_date = DateTime(2012,4,1)
    max_date = DateTime(2013,12,1)
    time_axis = add_float2datetime.(solcom,dt)
    time_axis2 = add_float2datetime.(1:1826,dt)



Como_Test1 = [Como_Test[1:23,1] Como_Test[1:23,2]./ 70]

com_fit2 = RMSE_fun(movingaverage(Ovi_com ,7),1:1826,Como_Test1 ,7, 1095, 0.1,0.001)

gcom = plot(time_axis2,movingaverage(Ovi_com   ,7) .* com_fit2[3], label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
    plot!(add_float2datetime.(Como_Test[:,1] .+ 1095 .+ com_fit2[2]  ,dt), Como_Test[:,2] ./ 70,   label = "Observed oviposition activity" )
    plot!(add_float2datetime.(Como_Test[:,1]  .+ 1095 .+ com_fit2[2] ,dt), Como_Test[:,2] ./ 70,   label = nothing, seriestype =:scatter, color=:orange  )
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(ylims = (-1,300))

##Trento - Roiz Marini Done

lat = 46.0748
alt = 194

tspan = (0.0, 5479.0);                    # Number of days to be simulated

Tre_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Tre_Clim.csv", header=true, DataFrame)

Tre_Temp = Tre_DP_Clim[:,3]
Tre_Prec = Tre_DP_Clim[:,4]
Tre_EV   = Tre_DP_Clim[:,5]
Tre_DP   = Tre_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,5479,length = 5479)),Tre_Prec[2193:7671])
Temps_spl = Spline1D(collect(range(1,5479,length = 5479)),Tre_Temp[2193:7671])
DP_spl = Spline1D(collect(range(1,5479,length = 5479)),Tre_DP[2193:7671])
EV_spl = Spline1D(collect(range(1,5479,length = 5479)),Tre_EV[2193:7671])

Trento_Adults =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Trento_Adults.csv", header=true,DataFrame)
Trento_Eggs =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Trento_Eggs.csv", header=true,DataFrame)
Mar_ad =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Marini_ad.csv", header=true,DataFrame)

#https://www.muse.it/it/La-Ricerca/Zoologia-invertebrati-idrobiologia/Azioni-sul-territorio/Documents/Andamento%20stagionale%20del%20numero%20di%20uova%20%282010-2020%29.pdf

include("AO2.jl")

soltre2 =sol.t
soltre =1:5479
adult_sum_tre = adult_sum
Egg_sum_tre = Ovi_Obs

dt = DateTime(2006,1,1)
    min_date = DateTime(2008,5,1)
    max_date = DateTime(2008,12,5)
    time_axis = add_float2datetime.(1:5479,dt)
    time_axis2 = add_float2datetime.(soltre2,dt)

tre_fit1 = RMSE_fun(movingaverage(adult_sum_tre ,100),soltre2,Trento_Adults,7,  .+ 730 , 0.3,0.001)

gtre1 = plot(time_axis2,movingaverage(adult_sum_tre   ,100) .* tre_fit1[3]  , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
            plot!(add_float2datetime.(Trento_Adults[:,1]  .+ 730 .+ tre_fit1[2],dt), Trento_Adults[:,2],   label = "Observed adults" )
            plot!(add_float2datetime.(Trento_Adults[:,1]  .+ 730 .+ tre_fit1[2],dt), Trento_Adults[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))

min_date = DateTime(2010,1,1)
max_date = DateTime(2020,3,5)

tre_fit2 = RMSE_fun(movingaverage(Egg_sum_tre ,7),soltre,Trento_Eggs,7, 1460 ,0.03,0.0001)

gtre2 = plot(time_axis,movingaverage(Egg_sum_tre   ,7)  .* tre_fit2[3]  , label = "Simulated oviposition activity", ylabel = "Oviposition activity")
            plot!(add_float2datetime.(Trento_Eggs[:,1] .+ 1460  .+ tre_fit2[2] ,dt), Trento_Eggs[:,2],   label = "Observed oviposition activity" )
            plot!(add_float2datetime.(Trento_Eggs[:,1]  .+ 1460 .+ tre_fit2[2],dt), Trento_Eggs[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(size=(1200,600))

min_date = DateTime(2014,1,1)
max_date = DateTime(2016,3,5)

tre_fit3 = RMSE_fun(movingaverage(adult_sum_tre ,100),soltre2,Mar_ad,7, 2555 .+ 365  , 0.2,0.001)

gtre3 = plot(time_axis2,movingaverage(adult_sum_tre   ,100) .* tre_fit3[3] ,   label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
            plot!(add_float2datetime.(Mar_ad[:,1] .+ 2555 .+ 365 .+ tre_fit3[2] ,dt), Mar_ad[:,2],   label = "Observed adults" )
            plot!(add_float2datetime.(Mar_ad[:,1]  .+2555  .+ 365  .+ tre_fit3[2],dt), Mar_ad[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-1.5,50))

###France##~

##Cagnes-sur-mer - Lacour et al. 2015 Done

lat = 43.6637
alt = 0
tspan = (0.0, 2508.0);
                    # Number of days to be simulated

Cag_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Cag_Clim.csv", header=true, DataFrame)
Cagnes_Eggs = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Cagnes_Eggs.csv", header=false, DataFrame)
Cagnes_O2 = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Cag_O2.csv", header=false, DataFrame)

Cag_Temp = Cag_DP_Clim[:,3]
Cag_Prec = Cag_DP_Clim[:,4]
Cag_EV   = Cag_DP_Clim[:,5]
Cag_DP   = Cag_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2508,length = 2508)),Cag_Prec[2193:4700])
Temps_spl = Spline1D(collect(range(1,2508,length = 2508)),Cag_Temp[2193:4700])
DP_spl = Spline1D(collect(range(1,2508,length = 2508)),Cag_DP[2193:4700])
EV_spl = Spline1D(collect(range(1,2508,length = 2508)),Cag_EV[2193:4700])


include("AO2.jl")

Egg_sum_cag = Ovi_Obs
solcag = 1:2508

cag_fit = RMSE_fun(movingaverage(Egg_sum_cag ,7),solcag,Cagnes_O2,8, 1095, 0.045,0.001)

dt = DateTime(2006,1,1)
              min_date = DateTime(2008,1,1)
              max_date = DateTime(2012,12,1)
                  time_axis = add_float2datetime.(1:2508,dt)

gcag = plot(time_axis,movingaverage(Egg_sum_cag,7) .* cag_fit[3]  , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
        plot!(add_float2datetime.(Cagnes_O2[:,1] .+ 1095 .+ cag_fit[2],dt)  , Cagnes_O2[:,2]     ,   label = "Observed oviposition activity" )
        plot!(add_float2datetime.(Cagnes_O2[:,1] .+ 1095 .+ cag_fit[2],dt)  , Cagnes_O2[:,2]     ,   label = nothing, seriestype =:scatter, color=:orange)
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylim = (-1,550))


cag_fit = RMSE_fun(movingaverage(Egg_sum_cag ,7),solcag,Cagnes_Eggs,8, 1460, 0.1,0.001)

dt = DateTime(2006,1,1)
              min_date = DateTime(2010,4,1)
              max_date = DateTime(2012,8,1)
            time_axis = add_float2datetime.(1:2508,dt)

gcag = plot(time_axis,movingaverage(Egg_sum_cag,7) .* cag_fit[3]  , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
        plot!(add_float2datetime.(Cagnes_Eggs[:,1] .+ 1460 .+ cag_fit[2],dt)  , Cagnes_Eggs[:,2]     ,   label = "Observed oviposition activity" )
        plot!(add_float2datetime.(Cagnes_Eggs[:,1] .+ 1460 .+ cag_fit[2],dt)  , Cagnes_Eggs[:,2]     ,   label = nothing, seriestype =:scatter, color=:orange)
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylim = (-10.5,350))

#
###Spain###

# ##Baix Llobregat DONE

lat = 37.613289786633416

alt = 20
tspan = (0.0, 5479.0);


Cart_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Cat_Clim.csv", header=true,DataFrame)
Cart_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Cart_Egg.csv", header=false,DataFrame)
Cata_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Cata_Egg.csv", header=true,DataFrame)

Cart_Temp = Cart_DP_Clim[:,3]
Cart_Prec = Cart_DP_Clim[:,4]
Cart_EV   = Cart_DP_Clim[:,5]
Cart_DP   = Cart_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,5479,length = 5479)),Cart_Prec[1:5479])
Temps_spl = Spline1D(collect(range(1,5479,length = 5479)),Cart_Temp[1:5479])
DP_spl = Spline1D(collect(range(1,5479,length = 5479)),Cart_DP[1:5479])
EV_spl = Spline1D(collect(range(1,5479,length = 5479)),Cart_EV[1:5479])

include("AO2.jl")

Egg_sum_Cart = Ovi_Obs



solCart = 1:365

days = 1:365
cart_Avg = zeros(365)
cart_Out = zeros(365)

for i in 1462:5479
    if days[mod(i,365) + 1] == solCart[mod(i, 365) + 1]

        global cart_Avg[mod(i,365)+1] = cart_Avg[mod(i,365)+1] + Egg_sum_Cart[i]
    end

    global cart_Out = cart_Avg./9
end

Cart_fit = RMSE_fun(movingaverage(cart_Out ,7),solCart,Cata_Egg,7, 0, 0.1,0.0001)

dt = DateTime(2000,1,1)
    min_date = DateTime(2014,1,1)
    max_date = DateTime(2014,12,31)
    time_axis = add_float2datetime.(solCart,dt)

gCart = plot(1:365,movingaverage(cart_Out ,7) .* Cart_fit[3], label = "Simulated oviposition activity", xlabel = "Day of year", ylabel = "Average oviposition activity")
    plot!(Cata_Egg[:,1]  .+ Cart_fit[2],Cata_Egg[:,2], label = "Observed oviposition activity"   )
    plot!(Cata_Egg[:,1]  .+ Cart_fit[2] ,Cata_Egg[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
    plot!(ylims = (0,160))
    plot!(xlims = (120,365))


#Majorca

lat = 39.6953

alt = 20
tspan = (0.0, 1461.0);


Maj_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Maj_Clim.csv", header=true,DataFrame)
Maj_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Maj_Egg.csv", header=true,DataFrame)

Maj_Temp = Maj_DP_Clim[:,3]
Maj_Prec = Maj_DP_Clim[:,4]
Maj_EV   = Maj_DP_Clim[:,5]
Maj_DP   = Maj_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Maj_Prec[4019:5479])
Temps_spl = Spline1D(collect(range(1,1461,length = 1461)),Maj_Temp[4019:5479])
DP_spl = Spline1D(collect(range(1,1461,length = 1461)),Maj_DP[4019:5479])
EV_spl = Spline1D(collect(range(1,1461,length = 1461)),Maj_EV[4019:5479])

include("AO2.jl")

Egg_sum_Maj = Ovi_Obs
solMaj = 1:1461

maj_fit = RMSE_fun(movingaverage(Egg_sum_Maj ,7),solMaj,Maj_Egg,7, 1095, 0.1,0.001)

dt = DateTime(2011,1,1)
    min_date = DateTime(2014,1,1)
    max_date = DateTime(2014,12,31)
    time_axis = add_float2datetime.(solMaj,dt)

gMaj = plot(time_axis,movingaverage(Egg_sum_Maj ,1) .* maj_fit[3], label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(add_float2datetime.(Maj_Egg[:,1] .+ 1095 .+ maj_fit[2],dt),Maj_Egg[:,2], label = "Observed oviposition activity"   )
    plot!(add_float2datetime.(Maj_Egg[:,1] .+ 1095 .+ maj_fit[2],dt),Maj_Egg[:,2],   label = nothing, seriestype =:scatter, color=:orange  )

#
##Goiri et al. 2020## DONE
#
lat = 43.3381
alt = 20
tspan = (0.0, 2922.0);


Iru_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Iru_Clim.csv", header=true,DataFrame)
Irun_Eggs1= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Iru_OM.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")
Irun_Eggs2= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Iru_OW.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")

Iru_Temp = Iru_DP_Clim[:,3]
Iru_Prec = Iru_DP_Clim[:,4]
Iru_EV   = Iru_DP_Clim[:,5]
Iru_DP   = Iru_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2922,length = 2922)),Iru_Prec[4750:7671])
Temps_spl = Spline1D(collect(range(1,2922,length = 2922)),Iru_Temp[4750:7671])
DP_spl = Spline1D(collect(range(1,2922,length = 2922)),Iru_DP[4750:7671])
EV_spl = Spline1D(collect(range(1,2922,length = 2922)),Iru_EV[4750:7671])

include("AO2.jl")

Egg_sum_iru = Ovi_Obs
soliru = 1:2922

Temp_Iru1 = zeros(20)

for i in 1:20
    Temp_Iru1[i] = (Irun_Eggs1[i,1] .-  Date(2013,1,1)).value
end

Temp_Iru11 = [Temp_Iru1[1:20,1] movingaverage(Irun_Eggs1[1:20,2],1)]



iru_fit1 = RMSE_fun(movingaverage(Egg_sum_iru ,7),soliru,Temp_Iru11[16:20,:],14, 0, 0.005,0.00001)
iru_fit2 = RMSE_fun(movingaverage(Egg_sum_iru ,7),soliru,Temp_Iru11[1:20,:],14, 0, 0.005,0.00001)


dt = DateTime(2013,1,1)
    min_date = DateTime(2015,5,1)
    max_date = DateTime(2018,11,5)
    time_axis = add_float2datetime.(soliru,dt)

giru = plot(time_axis,movingaverage(Egg_sum_iru ,7) .* iru_fit1[3] , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
    plot!(add_float2datetime.(Temp_Iru11[1:20,1] .+ iru_fit1[2],dt),Temp_Iru11[1:20,2], label = "Observed oviposition activity"   )
    plot!(add_float2datetime.(Temp_Iru11[1:20,1] .+ iru_fit1[2],dt),Temp_Iru11[1:20,2],   label = nothing, seriestype =:scatter, color=:orange  )
    plot!(xlims = Dates.value.([min_date,max_date]))
   plot!(ylims = (-1,50))
    Plots.vline!(add_float2datetime.([12725,13065,13072,13079,13397,13437,13707,13747,13777,13797,13817].-11680,dt),  line=:dash , label = "Control", linewidth = 0.5)
#
###Greece### DONE

##GIATROPOULOS et al. 2012
#
lat = 37.9838
alt = 60
tspan = (0.0, 4018.0);

Ath_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Ath_Clim.csv", header=true,DataFrame)
Gia_Dat =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Gia_O.csv", header=false,DataFrame);


Ath_Temp = Ath_DP_Clim[:,3]
Ath_Prec = Ath_DP_Clim[:,4]
Ath_EV   = Ath_DP_Clim[:,5]
Ath_DP   = Ath_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,4018,length = 4018)),Ath_Prec[1:4018])
Temps_spl = Spline1D(collect(range(1,4018,length = 4018)),Ath_Temp[1:4018])
DP_spl = Spline1D(collect(range(1,4018,length = 4018)),Ath_DP[1:4018])
EV_spl = Spline1D(collect(range(1,4018,length = 4018)),Ath_EV[1:4018])

include("AO2.jl")

solath = 1:4018
Egg_sum_ath = Ovi_Obs

ath_fit = RMSE_fun(movingaverage(Egg_sum_ath ,7),solath,Gia_Dat,7, 3289, 0.15,0.001)

dt = DateTime(2000,1,1)
    min_date = DateTime(2009,1,1)
    max_date = DateTime(2011,3,5)
    time_axis = add_float2datetime.(solath,dt)

gath = plot(time_axis,movingaverage(Egg_sum_ath,7) .*  ath_fit[3]  , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Predicted oviposition activity")
    plot!(add_float2datetime.(Gia_Dat[:,1]   .+ ath_fit[2] .+ 3289,dt),Gia_Dat[:,2], label = "Observed oviposition activity"   )
    plot!(add_float2datetime.(Gia_Dat[:,1]  .+ ath_fit[2] .+ 3289,dt),Gia_Dat[:,2],   label = nothing, seriestype =:scatter, color=:orange  )
    plot!(xlims = Dates.value.([min_date,max_date]))
    plot!(ylims = (-1,150))


#Croatia###

##Zitko and Merdic 2014## DONE

lat = 43.5081
alt = 0
tspan = (0.0, 1943.0);

Spl_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Spl_Clim.csv", header=true,DataFrame)
Zitka = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Zitka_Dat.csv", header=false,DataFrame)

Spl_Temp = Spl_DP_Clim[:,3]
Spl_Prec = Spl_DP_Clim[:,4]
Spl_EV   = Spl_DP_Clim[:,5]
Spl_DP   = Spl_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1943,length = 1943)),Spl_Prec[2558:4500])
Temps_spl = Spline1D(collect(range(1,1943,length = 1943)),Spl_Temp[2558:4500])
DP_spl = Spline1D(collect(range(1,1943,length = 1943)),Spl_DP[2558:4500])
EV_spl = Spline1D(collect(range(1,1943,length = 1943)),Spl_EV[2558:4500])

include("AO2.jl")

solspl2 = 1:1943
Ovi_spl = Ovi_Obs

zit_fit2 = RMSE_fun(movingaverage(Ovi_spl ,1),solspl2,Zitka,7, 730,3,0.001)

dt = DateTime(2007,1,1)
    min_date = DateTime(2009,5,1)
    max_date = DateTime(2010,11,1)
    time_axis2 = add_float2datetime.(solspl2,dt)


gspli =  plot(time_axis2,movingaverage(Ovi_spl ,1) .* zit_fit2[3] , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
        plot!(add_float2datetime.(Zitka[:,1] .+ 730 .+ zit_fit2[2],dt), Zitka[:,2]  ,   label = "Observed oviposition activity" )
        plot!(add_float2datetime.(Zitka[:,1] .+ 730 .+ zit_fit2[2],dt), Zitka[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylims = (-3,1000))

##Portugal### DONE
#Loule -  Osorio et al. (2020)

lat = 37.1379
alt = 49
tspan = (0.0, 1145.0);

Lou_DP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Lou_Clim.csv", header=true,DataFrame)

Lou_Temp = Lou_DP_Clim[:,3]
Lou_Prec = Lou_DP_Clim[:,4]
Lou_EV   = Lou_DP_Clim[:,5]
Lou_DP   = Lou_DP_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1144,length = 1144)),Lou_Prec[6211:7354])
Temps_spl = Spline1D(collect(range(1,1144,length = 1144)),Lou_Temp[6211:7354])
DP_spl = Spline1D(collect(range(1,1144,length = 1144)),Lou_DP[6211:7354])
EV_spl = Spline1D(collect(range(1,1144,length = 1144)),Lou_EV[6211:7354])

Os_Eggs= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Lou_Ov_M.csv", header=false,DataFrame)
Os_Ads= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Lou_Av_o.csv", header=false,DataFrame)

include("AO2.jl")

Egg_sum_lou = Ovi_Obs
adult_sum_lou = adult_sum

sollou = sol.t
sollou2 = 1:1145

lou_fit = RMSE_fun(movingaverage(Egg_sum_lou ,1),sollou2,Os_Eggs,15, 730, 1,0.001)
lou_fit2 = RMSE_fun(movingaverage(adult_sum_lou ,1),sollou,Os_Ads,15, 730, 1,0.001)

dt = DateTime(2017,1,1)
    min_date = DateTime(2019,2,1)
    max_date = DateTime(2020,1,1)
    time_axis = add_float2datetime.(sollou,dt)
    time_axis2 = add_float2datetime.(sollou2,dt)


glou=  plot(time_axis2,movingaverage(Egg_sum_lou ,1) .* lou_fit[3]  , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
   plot!(xlims = Dates.value.([min_date,max_date]))
   plot!(add_float2datetime.(Os_Eggs[:,1] .+ 730 .+ lou_fit[2],dt), Os_Eggs[:,2]  ,   label = "Observed Oviposition activity" )
   plot!(add_float2datetime.(Os_Eggs[:,1] .+ 730 .+ lou_fit[2],dt), Os_Eggs[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )

glouad=  plot(time_axis,movingaverage(adult_sum_lou ,100) .* lou_fit2[3]   , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
      plot!(xlims = Dates.value.([min_date,max_date]))
      plot!(add_float2datetime.(Os_Ads[:,1] .+ 730 .+ lou_fit2[2],dt), Os_Ads[:,2]  ,   label = "Observed adults" )
      plot!(add_float2datetime.(Os_Ads[:,1] .+ 730 .+ lou_fit2[2],dt), Os_Ads[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
      plot!(ylims = (-1,120))

###America####
# #
## Lake Charles - Willis and Nasci DONE
#
lat = 30.3322
alt = 4.5

tspan = (0.0, 1461.0);

Wil_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\LC_Clim.csv", header=true,DataFrame)
Wil_Dat = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wil_Dat.csv", header=false,DataFrame)
willis_time_2= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\wil_time_2.csv", header=false,DataFrame)
Wil_Ad = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wil_P.csv", header=false,DataFrame)

Wil_Temp = Wil_Clim[:,3]
Wil_Prec = Wil_Clim[:,4]
Wil_EV =   Wil_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Wil_Prec[1:1461] )
Temps_spl = Spline1D(collect(range(1,1461,length =  1461)),Wil_Temp[1:1461])
EV_spl = Spline1D(collect(range(1,1461,length = 1461)),Wil_EV[1:1461])

include("AO2.jl")


solwil = sol.t
avg_wing_wil = avg_wing
out_wing_wil = out_wing
adult_sum_wil = adult_sum
big_a_wil = big_adult_sum

dt = DateTime(1988,1,1)
    min_date = DateTime(1990,3,1)
    max_date = DateTime(1991,2,1)
    time_axis = add_float2datetime.(solwil,dt)

wil_fit = RMSE_fun(movingaverage(adult_sum_wil ,100),solwil,Wil_Ad,30, 730, 0.05,0.001)

gwil =  heatmap(time_axis,
        Wing_Vals, out_wing_wil,
        c=cgrad(:oslo, categorical = true),
        xlabel="Time (Days)", ylabel="Wing length (mm)", colorbar_title = "Proportion of adults", background_color_inside = "black")
        plot!(clims = (0,1))
        plot!(ylims = (1.5,4))
        plot!(time_axis,movingaverage(avg_wing_wil,100), label = "Simulated average wing length", color = :green ,linewidth = 2)
        plot!(add_float2datetime.(Wil_Dat[:,1]  .+ 730  ,dt) ,Wil_Dat[:,2], label = "Average wing length of captured adults", color = :orange ,linewidth = 2)
        plot!(add_float2datetime.(Wil_Dat[:,1]  .+ 730  ,dt) ,Wil_Dat[:,2], label = nothing, color = :orange, seriestype =:scatter)
        plot!(add_float2datetime.(willis_time_2[:,1]  .+ 365 ,dt) ,willis_time_2[:,2], label = nothing, color = :purple, seriestype =:scatter)
        plot!(add_float2datetime.(willis_time_2[:,1]  .+ 365 ,dt) ,willis_time_2[:,2], label = "Average wing length of emerging adults", color = :purple ,linewidth = 2)
        plot!(xlims = Dates.value.([min_date,max_date]))


plot(time_axis,movingaverage(avg_wing_wil,100), label = "Simulated average wing length" ,linewidth = 2, ylabel = "Average wing length (mm)")
  plot!(add_float2datetime.(Wil_Dat[:,1]  .+ 730  ,dt) ,Wil_Dat[:,2], label = "Field observed wing length" ,linewidth = 2)
  plot!(add_float2datetime.(Wil_Dat[:,1]  .+ 730  ,dt) ,Wil_Dat[:,2], label = nothing, color = :orange, seriestype =:scatter)
  plot!(xlims = Dates.value.([min_date,max_date]))
  plot!(ylims = (2.0,4))


gwilad=  plot(time_axis,movingaverage(adult_sum_wil ,100) .* wil_fit[3], label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
                plot!(add_float2datetime.(Wil_Ad[:,1] .+ 730 .+ wil_fit[2] ,dt), Wil_Ad[:,2]  ,   label = "Observed adults" )
                plot!(add_float2datetime.(Wil_Ad[:,1] .+ 730 .+ wil_fit[2],dt), Wil_Ad[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims = (-0.2,9))

                gwilad=  plot(time_axis,movingaverage(adult_sum_wil ,100) .* wil_fit[3], label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults",linewidth = 2)
                plot!(add_float2datetime.(Wil_Ad[:,1] .+ 730 .+ wil_fit[2] ,dt), Wil_Ad[:,2]  ,   label = "Field observed adults" ,linewidth = 2)
                plot!(add_float2datetime.(Wil_Ad[:,1] .+ 730 .+ wil_fit[2],dt), Wil_Ad[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims = (-0.2,9))


#VectorBase - Fort Worth## DONE

lat = 30.3322
alt = 199
tspan = (0.0, 3347.0);

FW_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\FW_Clim.csv", header=true,DataFrame)

FW_Temp = FW_Clim[:,3]
FW_Prec = FW_Clim[:,4]
FW_EV =   FW_Clim[:,5]
FW_DP =   FW_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,3347,length = 3347)),FW_Prec[3654:7000] )
Temps_spl = Spline1D(collect(range(1,3347,length =  3347)),FW_Temp[3654:7000])
DP_spl = Spline1D(collect(range(1,3347,length = 3347)),FW_DP[3654:7000])
EV_spl = Spline1D(collect(range(1,3347,length = 3347)),FW_EV[3654:7000])

FW_Adults_CDC = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FW_CDC.csv", header=true,DataFrame;  dateformat="dd/mm/yyyy")
FW_Adults_BG = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FW_BG.csv", header=true,DataFrame;  dateformat="dd/mm/yyyy")


include("AO2.jl")

adult_sum_fw = adult_sum
big_adult_sum_fw = big_adult_sum
solfw = sol.t

dt = DateTime(2010,1,1)
  min_date = DateTime(2015,1,1)
  max_date = DateTime(2018,1,1)
  time_axis = add_float2datetime.(solfw,dt)

Temp_West = zeros(122)

for i in 1:122
    Temp_West[i] = (FW_Adults_CDC[i,1] .-  Date(2010,1,1)).value
end

Temp_West1 = [Temp_West[1:122,1] movingaverage(FW_Adults_CDC[1:122,2],1)]

Temp_BG = zeros(71)

for i in 1:71
    Temp_BG[i] = (FW_Adults_CDC[i,1] .-  Date(2010,1,1)).value
end

Temp_BG1 = [Temp_BG[1:71,1] movingaverage(FW_Adults_CDC[1:71,2],1)]

fw_fit1 = RMSE_fun(movingaverage(adult_sum_fw ,1),solfw,Temp_BG1,14, 0,2,0.01)
fw_fit2 = RMSE_fun(movingaverage(adult_sum_fw ,1),solfw,Temp_West1,14, 0,0.05,0.0001)
fw_fit3 = RMSE_fun(movingaverage(big_adult_sum_fw ,1),solfw,Temp_West1,14, 0,0.1,0.001)

gfwcdc=  plot(time_axis,movingaverage(adult_sum_fw ,1) .* fw_fit1[3]  , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
        plot!(DateTime.(FW_Adults_BG[:,1] ) , FW_Adults_BG[:,2]  ,   label = "Observed adults" )
        plot!((DateTime.(FW_Adults_BG[:,1])),FW_Adults_BG[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))

gfwcdc2=  plot(time_axis,movingaverage(adult_sum_fw ,1) .* fw_fit2[3]  , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
        plot!(DateTime.(FW_Adults_CDC[:,1] ) , FW_Adults_CDC[:,2]  ,   label = "Observed adults" )
        plot!((DateTime.(FW_Adults_CDC[:,1])),FW_Adults_CDC[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))

gfwcdc3=  plot(time_axis,movingaverage(big_adult_sum_fw ,1) .* fw_fit3[3]  , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
        plot!(DateTime.(FW_Adults_CDC[:,1] ) , FW_Adults_CDC[:,2]  ,   label = "Observed adults" )
        plot!((DateTime.(FW_Adults_CDC[:,1])),FW_Adults_CDC[:,2] ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))


##New Orleans - Comiskey## DONE
#
Orl_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Orl_Clim.csv", header=true,DataFrame)

Orl_Temp = Orl_Clim[:,3]
Orl_Prec = Orl_Clim[:,4]
Orl_EV =   Orl_Clim[:,5]

Precip_spl = Spline1D(collect(range(1,1827,length = 1827)),Orl_Prec[1:1827] )
Temps_spl  = Spline1D(collect(range(1,1827,length =  1827)),Orl_Temp[1:1827])
DP_spl     = Spline1D(collect(range(1,1827,length = 1827)),Orl_EV[1:1827])
EV_spl     = Spline1D(collect(range(1,1827,length = 1827)),Orl_EV[1:1827])

lat = 30
alt = 0

tspan = (0.0,1827)

include("AO2.jl")

solorl = sol.t
avg_wing_orl = avg_wing
out_wing_orl = out_wing
larv_orl = out[4,:]
ads_orl = adult_sum

dt = DateTime(1992,1,1)
    min_date = DateTime(1995,1,1)
    max_date = DateTime(1996,3,1)
    time_axis = add_float2datetime.(solorl,dt)

com_time = [91,4*31,5*31,6*31,7*31,8*31,9*31,10*31,11*31] .+ 15
com_dat  = [2.68,2.56,2.46,2.43,2.49,2.51,2.49,2.51,2.67]
com_min = [2.25,2.07,1.8,2.04,2.1,1.98,1.92,2.34,2.37]
com_max = [3.2,2.88,2.88,2.76,2.91,3,2.97,2.61,2.91]
com_larv = [16,15,24,28,24,16,15,15,7]
com_ads = [27,19,15,24,39,21,26,3,13]

com_Dat = [com_time     com_dat]

gorl =  heatmap(time_axis,
    Wing_Vals, out_wing_orl,
    c=cgrad(:oslo, categorical = true),
    xlabel="Time (Days)", ylabel="Wing length (mm)", colorbar_title = "Proportion of adults", background_color_inside = "black")
    plot!(clims = (0,1))
    plot!(ylims = (1.7,4))
    plot!(time_axis,movingaverage(avg_wing_orl,1000), label = "Simulated average wing length", color = :green ,linewidth = 2)
    plot!(add_float2datetime.(com_time[:,1] .+ 1095  ,dt)  ,com_dat, label ="Observed average wing length")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 ,dt)  ,com_max, label ="Maximum observed wing length")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 ,dt) ,com_min, label ="Minimum observed wing length")
    plot!(xlims = Dates.value.([min_date,max_date]))

orl_L = [com_time   com_larv]
orl_A = [com_time   com_ads]

orl_fitl = RMSE_fun(movingaverage(larv_orl ,100),solorl,orl_L,14, 1095,0.02,0.0001)
orl_fitA = RMSE_fun(movingaverage(ads_orl ,100),solorl,orl_A,14, 1095,0.1,0.001)


gorl_L = Plots.plot(time_axis,movingaverage(larv_orl ,100) .* orl_fitl[3]     , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 .+ orl_fitl[2] ,dt)  ,com_larv, label ="Observed larvae")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 .+ orl_fitl[2] ,dt)  ,com_larv  ,   label = nothing, seriestype =:scatter , color = :orange)
    plot!(xlims = Dates.value.([min_date,max_date]))

gorl_A = Plots.plot(time_axis,movingaverage(ads_orl ,100) .* orl_fitA[3]    , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 .+ orl_fitA[2] ,dt)  ,com_ads, label ="Observed adults")
    plot!(add_float2datetime.(com_time[:,1] .+ 1095 .+ orl_fitA[2] ,dt)  ,com_ads,   label = nothing, seriestype =:scatter , color = :orange)
    plot!(xlims = Dates.value.([min_date,max_date]))

# ##Stratford, Conneticut - Armstrong DONE
#
CT_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Str_Clim.csv", header=true,DataFrame)
CT_Larv= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Arm_LBW.csv", header=true,DataFrame)
CT_Ad = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Arm_ABW.csv", header=true,DataFrame)

CT_Temp = CT_Clim[:,3]
CT_Prec = CT_Clim[:,4]
CT_EV =   CT_Clim[:,5]
CT_DP =   CT_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2592,length = 2592)),CT_Prec[3654:6245] )
Temps_spl  = Spline1D(collect(range(1,2592,length =  2592)),CT_Temp[3654:6245])
DP_spl     = Spline1D(collect(range(1,2592,length = 2592)),CT_DP[3654:6245])
EV_spl     = Spline1D(collect(range(1,2592,length = 2592)),CT_EV[3654:6245])
alt = 50
lat = 41.1845
tspan = (0.0, 2592.0);

include("AO2.jl")

solct = sol.t
adult_sum_ct = adult_sum
larv_sum_ct = larv_sum

dt = DateTime(2010,1,1)
    min_date = DateTime(2013,1,1)
    max_date = DateTime(2017,1,1)
    time_axis = add_float2datetime.(solct,dt)

CT_Ad2 = CT_Ad[1:33,:]

CT_L2 =
        CT_Larv[1:43,:]


ct_fitA = RMSE_fun(movingaverage(adult_sum_ct ,100),solct,CT_Ad2,7, 1095,1,0.001)
ct_fitL = RMSE_fun(movingaverage(larv_sum_ct ,100),solct,CT_L2,7, 1095,0.01,0.0001)
ct_fitA = RMSE_fun(movingaverage(adult_sum_ct ,100),solct,CT_Ad2[26:33,:],7, 1095,1,0.001)
ct_fitL = RMSE_fun(movingaverage(larv_sum_ct ,100),solct,CT_L2[35:43,:],7, 1095,0.1,0.001)


gct_A = Plots.plot(time_axis,movingaverage(adult_sum_ct ,100) .* ct_fitA[3]      , label = "Predicted adults", xlabel = "Time (Days)", ylabel = "Number of Adults")
            plot!(add_float2datetime.(CT_Ad[1:10,1] .+ 1095  .+ ct_fitA[2]   ,dt), CT_Ad[1:10,2]    ,   label = "Observed adults", color=:orange  )
            plot!(add_float2datetime.(CT_Ad[1:10,1] .+ 1095  .+ ct_fitA[2] ,dt), CT_Ad[1:10,2]     ,   label = nothing, seriestype =:scatter, color=:orange)
            plot!(add_float2datetime.(CT_Ad[11:18,1] .+ 1095 .+ ct_fitA[2]   ,dt), CT_Ad[11:18,2]    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.(CT_Ad[11:18,1] .+ 1095  .+ ct_fitA[2] ,dt), CT_Ad[11:18,2]     ,   label = nothing, color=:orange)
            plot!(add_float2datetime.(CT_Ad[19:25,1] .+ 1095 .+ ct_fitA[2]   ,dt), CT_Ad[19:25,2]    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.(CT_Ad[19:25,1] .+ 1095  .+ ct_fitA[2] ,dt), CT_Ad[19:25,2]     ,   label = nothing, color=:orange)
            plot!(add_float2datetime.(CT_Ad[26:33,1] .+ 1095 .+ ct_fitA[2]   ,dt), CT_Ad[26:33,2]    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.(CT_Ad[26:33,1] .+ 1095  .+ ct_fitA[2] ,dt), CT_Ad[26:33,2]     ,   label = nothing, color=:orange)
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-1.5,30))

gct_L = Plots.plot(time_axis,movingaverage(larv_sum_ct ,100) .* ct_fitL[3]  , label = "Predicted larvae",  xlabel = "Time (Days)", ylabel = "Number of larvae")
            plot!(add_float2datetime.(CT_Larv[1:13,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[1:13,2]   ,  color=:orange,    label = "Observed larvae"  )
            plot!(add_float2datetime.(CT_Larv[1:13,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[1:13,2]   , color=:orange,   label = nothing , seriestype =:scatter)
            plot!(add_float2datetime.(CT_Larv[14:26,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[14:26,2]   , color=:orange,   label = nothing, seriestype =:scatter  )
            plot!(add_float2datetime.(CT_Larv[14:26,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[14:26,2]   , color=:orange,   label = nothing )
            plot!(add_float2datetime.(CT_Larv[27:34,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[27:34,2]   , color=:orange,   label = nothing, seriestype =:scatter  )
            plot!(add_float2datetime.(CT_Larv[27:34,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[27:34,2]   , color=:orange,   label = nothing )
            plot!(add_float2datetime.(CT_Larv[35:43,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[35:43,2]   , color=:orange,   label = nothing, seriestype =:scatter  )
            plot!(add_float2datetime.(CT_Larv[35:43,1] .+ 1095 .+ ct_fitL[2],dt), CT_Larv[35:43,2]   , color=:orange,   label = nothing )
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims = (-3.6,120))

#                 plot!(xlims = Dates.value.([min_date,max_date]))

#
# ##Monmouth, New Jersey - Fonesca DONE

NJ_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Mon_Clim.csv", header=true,DataFrame)

NJ_Temp = NJ_Clim[:,3]
NJ_Prec = NJ_Clim[:,4]
NJ_EV =   NJ_Clim[:,5]
NJ_DP =   NJ_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7000,length = 7000)),NJ_Prec[1:7000] )
Temps_spl  = Spline1D(collect(range(1,7000,length =  7000)),NJ_Temp[1:7000])
DP_spl     = Spline1D(collect(range(1,7000,length = 7000)),NJ_DP[1:7000])
EV_spl     = Spline1D(collect(range(1,7000,length = 7000)),NJ_EV[1:7000])

NJ_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Fon_O3.csv", header=false,DataFrame)
NJ_Ad= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Fon_Ad.csv", header=false,DataFrame)

lat = 40.2589
alt = 20

tspan = (0.0, 4000.0);

include("AO2.jl")

solnj = sol.t
adult_sum_nj = adult_sum
egg_sum_nj = Ovi_Obs
big_adult_sum_nj = big_adult_sum
dt = DateTime(2000,1,1)
  min_date = DateTime(2009,4,1)
  max_date = DateTime(2010,1,1)
  time_axis = add_float2datetime.(solnj,dt)
  time_axis2 = add_float2datetime.(1:4000,dt)

nj_fit = RMSE_fun(movingaverage(egg_sum_nj ,7),1:4000,NJ_Egg,7, 3285,0.01,0.0001)

gnj1  = Plots.plot(time_axis2,movingaverage(egg_sum_nj ,7) .* nj_fit[3]  , label = "Simulated ovipostition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
        plot!(add_float2datetime.(NJ_Egg[:,1] .+ 3285 .+   nj_fit[2]  ,dt), NJ_Egg[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(add_float2datetime.(NJ_Egg[:,1] .+ 3285 .+   nj_fit[2],dt), NJ_Egg[:,2]   ,   label ="Observed ovipostition activity", color=pal[2]  )
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylims = (-1,40))

nj_fit2 = RMSE_fun(movingaverage(adult_sum_nj ,1),solnj,NJ_Ad,7, 3285,0.2,0.001)


gnj2 = Plots.plot(time_axis,movingaverage(adult_sum_nj ,1).*  nj_fit2[3]    , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
                plot!(add_float2datetime.(NJ_Ad[:,1] .+ 3285 .+  nj_fit2[2]   ,dt), NJ_Ad[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
                plot!(add_float2datetime.(NJ_Ad[:,1] .+ 3285 .+  nj_fit2[2],dt), NJ_Ad[:,2]   ,   label ="Observed adults", color=pal[2]  )
                plot!(xlims = Dates.value.([min_date,max_date]))

nj_fit3 = RMSE_fun(movingaverage(big_adult_sum_nj ,100),solnj,NJ_Ad,7, 3285,6,0.01)


gnj3 = Plots.plot(time_axis,movingaverage(big_adult_sum_nj ,100).*  nj_fit3[3]    , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
            plot!(add_float2datetime.(NJ_Ad[:,1] .+ 3285 .+  nj_fit3[2]   ,dt), NJ_Ad[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.(NJ_Ad[:,1] .+ 3285 .+  nj_fit3[2],dt), NJ_Ad[:,2]   ,   label ="Observed adults", color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))


#TTU, Lubbock, Texas - Erickson

TTU_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Lub_Clim.csv", header=true,DataFrame)

TTU_Temp = TTU_Clim[:,3]
TTU_Prec = TTU_Clim[:,4]
TTU_EV =   TTU_Clim[:,5]
TTU_DP =   TTU_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7000,length = 7000)),TTU_Prec[1:7000] )
Temps_spl  = Spline1D(collect(range(1,7000,length =  7000)),TTU_Temp[1:7000])
DP_spl     = Spline1D(collect(range(1,7000,length = 7000)),TTU_DP[1:7000])
EV_spl     = Spline1D(collect(range(1,7000,length = 7000)),TTU_EV[1:7000])

TTU_Adults =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\TTU_Adults.csv", header=true,DataFrame)
TTU_Good =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\TTU_2007.csv", header=false,DataFrame)


lat = 33.5843
alt = 976

tspan = (0.0, 7000.0);

include("AO2.jl")

solttu = sol.t
adult_sum_ttu = adult_sum

dt = DateTime(2000,1,1)
  min_date = DateTime(2004,1,1)
  max_date = DateTime(2009,1,1)
  time_axis = add_float2datetime.(solttu,dt)

ttu_fit = RMSE_fun(movingaverage(adult_sum_ttu ,7),solttu,TTU_Adults[11:74,:],14, 730,0.01,0.0001)


gttu = Plots.plot(time_axis,movingaverage(adult_sum_ttu ,1) .* ttu_fit[3]  , label = "Predicted Adults", color = pal[1], xlabel = "Time (Days)", ylabel = "Number of Adults")
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(add_float2datetime.(TTU_Adults[11:74,1] .+730   ,dt), TTU_Adults[11:74,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(add_float2datetime.(TTU_Adults[11:74,1] .+ 730, dt), TTU_Adults[11:74,2]   , label = "Observed Adults" , color = pal[2])

#VectorBase- Asheville## DONE

Ashe_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Ash_Clim.csv", header=true,DataFrame)

Ashe_Temp = Ashe_Clim[:,3]
Ashe_Prec = Ashe_Clim[:,4]
Ashe_EV =   Ashe_Clim[:,5]
Ashe_DP =   Ashe_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1886,length = 1886)),Ashe_Prec[5115:7000] )
Temps_spl  = Spline1D(collect(range(1,1886,length =  1886)),Ashe_Temp[5115:7000])
DP_spl     = Spline1D(collect(range(1,1886,length = 1886)),Ashe_DP[5115:7000])
EV_spl     = Spline1D(collect(range(1,1886,length = 1886)),Ashe_EV[5115:7000])

Ashe_Ovi =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Ash_OviW.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")

lat = 35.5951
alt = 650

tspan = (0.0, 1886.0);

include("AO2.jl")

solash = sol.t
Ovi_ash = Ovi_Obs

dt = DateTime(2014,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,1,1)
  time_axis = add_float2datetime.(1:1886,dt)


Temp_ash = zeros(28)

for i in 1:28
  Temp_ash[i] = (Ashe_Ovi[i,1] .-  Date(2014,1,1)).value
end

Temp_ash1 = [Temp_ash[1:28,1] Ashe_Ovi[1:28,2]]

ash_fit_O = RMSE_fun(movingaverage(Ovi_Obs ,1),1:1886,Temp_ash1,14, 0,0.002,0.00001)


gash =  plot(time_axis,movingaverage(Ovi_ash ,1) .* ash_fit_O[3]   , ylabel = "Oviposition activity", xlabel = "Time (Days)", label = "Predicted oviposition activity")
   plot!(add_float2datetime.(Temp_ash[:,1]  .+ ash_fit_O[2],dt), Ashe_Ovi[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
   plot!(add_float2datetime.(Temp_ash[:,1] .+ ash_fit_O[2] ,dt), Ashe_Ovi[:,2]   ,    color=:orange, label = "Observed oviposition activity")
   plot!(xlims = Dates.value.([min_date,max_date]))
   plot!(ylims = (-1,30))

##VecotrBase - Charlotte Mundis, VB, Whiteman DONE
#
Char_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Cha_Clim.csv", header=true,DataFrame)

Char_Temp = Char_Clim[:,3]
Char_Prec = Char_Clim[:,4]
Char_EV =   Char_Clim[:,5]
Char_DP =   Char_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1886,length = 1886)),Char_Prec[5115:7000] )
Temps_spl  = Spline1D(collect(range(1,1886,length =  1886)),Char_Temp[5115:7000])
DP_spl     = Spline1D(collect(range(1,1886,length = 1886)),Char_DP[5115:7000])
EV_spl     = Spline1D(collect(range(1,1886,length = 1886)),Char_EV[5115:7000])

Char_Adults =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Char_Ad1.csv", header=true,DataFrame; dateformat="dd/mm/yyyy")
Char_Ovi =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Char_Ovi.csv", header=true,DataFrame; dateformat="dd/mm/yyyy")
Char_T=  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Char_T.csv", header=false,DataFrame)


lat = 35.2271
alt = 232

tspan = (0.0, 1501.0);

include("AO2.jl")

solcha = sol.t
adult_sum_cha = adult_sum
big_adult_sum_cha = big_adult_sum
Ovi_cha = Ovi_Obs

dt = DateTime(2014,1,1)
  min_date = DateTime(2017,1,1)
  max_date = DateTime(2018,1,1)
  time_axis = add_float2datetime.(solcha,dt)
  time_axis2 = add_float2datetime.(1:1501,dt)

Temp_CH = zeros(12)

for i in 1:12
  Temp_CH[i] = (Char_Adults[i,1] .-  Date(2014,1,1)).value
end

Temp_CH1 = [Temp_CH[1:12,1] movingaverage(Char_Adults[1:12,2],1)]

Temp_OV = zeros(24)

for i in 1:24
  Temp_OV[i] = (Char_Ovi[i,1] .-  Date(2014,1,1)).value
end

Temp_Ov1 = [Temp_OV[1:24,1] movingaverage(Char_Ovi[1:24,2],1)]

cha_fit_A1 = RMSE_fun(movingaverage(adult_sum_cha ,100),solcha,Temp_CH1,10, 0,3,0.001)
cha_fit_O = RMSE_fun(movingaverage(Ovi_Obs ,7),1:1501,Temp_Ov1,14, 0,0.02,0.0001)

gcha_A =  plot(time_axis,movingaverage(adult_sum_cha ,1) .* cha_fit_A1[3], label = "Predicted adults", xlabel = "Time (Days)", ylabel = "Number of adults")
     plot!(xlims = Dates.value.([min_date,max_date]))
     plot!(add_float2datetime.(Temp_CH1[:,1]  .+ cha_fit_A1[2] ,dt), Char_Adults[:,2]   ,   label = "Observed adults")
     plot!(add_float2datetime.(Temp_CH1[:,1] .+ cha_fit_A1[2],dt), Char_Adults[:,2]   ,   label = nothing, color=:orange , seriestype =:scatter )
     plot!(ylims = (-0.1,8))


min_date = DateTime(2016,1,1)
max_date = DateTime(2017,1,1)

gcha_O =  plot(time_axis2,movingaverage(Ovi_cha ,7) .* cha_fit_O[3], label = "Predicted oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
      plot!(xlims = Dates.value.([min_date,max_date]))
      plot!(add_float2datetime.(Temp_Ov1[:,1]  .+ cha_fit_O[2] ,dt), Char_Ovi[:,2]   ,   label = nothing, seriestype =:scatter , color=:orange )
      plot!(add_float2datetime.(Temp_Ov1[:,1] .+ cha_fit_O[2] ,dt), Char_Ovi[:,2]   ,   label = "Observed oviposition activity", color = pal[2])
      plot!(ylims = (-3,100))

min_date = DateTime(2017,1,1)
max_date = DateTime(2018,1,1)

gCT =  heatmap(time_axis,
         Wing_Vals, out_wing,
         c=cgrad(:oslo, categorical = true),
         xlabel="Time (Days)", ylabel="Wing length (mm)", colorbar_title = "Proportion of adults", background_color_inside = "black")
         plot!(clims = (0,1))
         plot!(ylims = (1.7,4))
         plot!(time_axis,movingaverage(avg_wing,100), label = "Simulated average wing length", color = :green ,linewidth = 2)
         plot!(add_float2datetime.(Char_T[:,1] .+ 1095 .+ 10 ,dt), Char_T[:,2]  , label = "Average wing length of emerging adults", color = :orange ,linewidth = 2 )
         plot!(add_float2datetime.(Char_T[:,1] .+ 1095 .+ 10,dt), Char_T[:,2]   , label = nothing,color = :orange, seriestype =:scatter)
         plot!(xlims = Dates.value.([min_date,max_date]))


##VecotrBase - Raliegh DONE

Ragh_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Rag_Clim.csv", header=true,DataFrame)

Ragh_Temp = Ragh_Clim[:,3]
Ragh_Prec = Ragh_Clim[:,4]
Ragh_EV =   Ragh_Clim[:,5]
Ragh_DP =   Ragh_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1201,length = 1201)),Ragh_Prec[5115:6315] )
Temps_spl  = Spline1D(collect(range(1,1201,length =  1201)),Ragh_Temp[5115:6315])
DP_spl     = Spline1D(collect(range(1,1201,length = 1201)),Ragh_DP[5115:6315])
EV_spl     = Spline1D(collect(range(1,1201,length = 1201)),Ragh_EV[5115:6315])

Rag_Ovi =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Rag_Ovi.csv", header=true,DataFrame; dateformat="dd/mm/yyyy")

lat = 35.7796
alt = 96

tspan = (0.0, 1201.0);

include("AO2.jl")

solrag = sol.t
Ovi_rag = Ovi_Obs

dt = DateTime(2014,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,3,1)
  time_axis = add_float2datetime.(1:1201,dt)

Temp_rag = zeros(28)

for i in 1:28
  Temp_rag[i] = (Rag_Ovi[i,1] .-  Date(2014,1,1)).value
end

Temp_rag1 = [Temp_rag[1:28,1] Rag_Ovi[1:28,2]]

rag_fit = RMSE_fun(movingaverage(Ovi_Obs ,1),1:1201,Temp_rag1,14, 0,0.02,0.0001)


grag =  plot(time_axis,movingaverage(Ovi_rag ,1) .*  rag_fit[3]   , ylabel = "Oviposition activity", xlabel = "Time (Days)", label = "Predicted oviposition activity")
     plot!(xlims = Dates.value.([min_date,max_date]))
     plot!(add_float2datetime.(Temp_rag1[:,1] .+ 0 .+ rag_fit[2],dt),Temp_rag1[:,2], label = "Observed oviposition activity")
     plot!(add_float2datetime.(Temp_rag1[:,1] .+ 0 .+ rag_fit[2] ,dt),Temp_rag1[:,2],  label = nothing, seriestype =:scatter, color=:orange  )
     plot!(ylim = (-5,120))


##Vectorbase Washington

Wash_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Was_Clim.csv", header=true,DataFrame)

Wash_Temp = Wash_Clim[:,3]
Wash_Prec = Wash_Clim[:,4]
Wash_EV =   Wash_Clim[:,5]
Wash_DP =   Wash_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7853,length = 7853)),Wash_Prec)
Temps_spl  = Spline1D(collect(range(1,7853,length =  7853)),Wash_Temp)
DP_spl     = Spline1D(collect(range(1,7853,length = 7853)),Wash_DP)
EV_spl     = Spline1D(collect(range(1,7853,length = 7853)),Wash_EV)

#Wash_Adults_BG =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wash_Ad_BG_lure_carbon.csv", DataFrame,header=true; dateformat="yyyy-mm-dd")
#Wash_Adults_CDC =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wash_Ad_CDC_light_cd.csv", DataFrame,header=true; dateformat="yyyy-mm-dd")
Wash_Adults_CDC =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wash_ADM.csv", DataFrame,header=true; dateformat="dd/mm/yyyy")

lat = 38.834758449896114
alt = 520

tspan = (0.0, 7853.0);

include("AO2.jl")

solwas = sol.t
adult_sum_was = adult_sum

dt = DateTime(2000,1,1)
    min_date = DateTime(2010,3,1)
    max_date = DateTime(2020,1,1)
    time_axis = add_float2datetime.(solwas,dt)
    time_axis2 = add_float2datetime.(1:7853,dt)
Wash_Ad = Wash_Adults_CDC

Temp_wash = zeros(61)

for i in 1:61
  Temp_wash[i] = (Wash_Adults_CDC[i,1] .-  Date(2000,1,1)).value
end

Temp_wash1 = [Temp_wash[1:61,1] movingaverage(Wash_Adults_CDC[1:61,2],1)]


was_fit2 = RMSE_fun(movingaverage(adult_sum_was ,1),solwas,Temp_wash1[11:61,:],14, 0,0.1,0.001)

gwas_CDC = plot(time_axis,movingaverage(adult_sum_was ,1) .* was_fit2[3], label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
         plot!(add_float2datetime.(Temp_wash1[11:15,1].+ was_fit2[2],dt )   ,Float64.(Temp_wash1[11:15,2])  ,   label = nothing, seriestype =:scatter,  color=:orange  )
         plot!(add_float2datetime.(Temp_wash1[11:15,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[11:15,2])   ,   label = "Observed adults", color=pal[2]  )
         plot!(add_float2datetime.(Temp_wash1[16:19,1].+ was_fit2[2],dt )   ,Float64.(Temp_wash1[16:19,2])  ,   label = nothing, seriestype =:scatter, color=:orange  )
         plot!(add_float2datetime.(Temp_wash1[16:19,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[16:19,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[20:24,1].+ was_fit2[2],dt )   ,Float64.(Temp_wash1[20:24,2])  ,   label = nothing, seriestype =:scatter, color=:orange  )
         plot!(add_float2datetime.(Temp_wash1[20:24,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[20:24,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[25:30,1].+ was_fit2[2],dt )   ,Float64.(Temp_wash1[25:30,2])  ,   label = nothing, seriestype =:scatter, color=:orange  )
         plot!(add_float2datetime.(Temp_wash1[25:30,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[25:30,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[31:35,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[31:35,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[36:40,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[36:40,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[41:46,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[41:46,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[47:51,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[47:51,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[52:55,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[52:55,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[56:61,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[56:61,2])   ,   label = nothing, color=pal[2] )
         plot!(add_float2datetime.(Temp_wash1[31:61,1].+ was_fit2[2] ,dt)   ,Float64.(Temp_wash1[31:61,2])   ,   label = nothing, seriestype =:scatter, color=:orange  )
         plot!(xlims = Dates.value.([min_date,max_date]))
         plot!(ylims = (-1,90))

#VectorBase - Fort Myers

FM_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\FM_Clim.csv", header=true,DataFrame)

FM_Temp = FM_Clim[:,3]
FM_Prec = FM_Clim[:,4]
FM_EV =   FM_Clim[:,5]
FM_DP =   FM_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7000,length = 7000)),FM_Prec[1:7000] )
Temps_spl  = Spline1D(collect(range(1,7000,length =  7000)),FM_Temp[1:7000])
DP_spl     = Spline1D(collect(range(1,7000,length = 7000)),FM_DP[1:7000])
EV_spl     = Spline1D(collect(range(1,7000,length = 7000)),FM_EV[1:7000])

# FM_CDCg =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FM_CDCg.csv", header=true,DataFrame)
# FM_CDCl =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FM_CDCl.csv", header=true,DataFrame)
# FM_BG =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FM_BG.csv", header=true,DataFrame)
FM_BG =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\FM_BG.csv", header=true,DataFrame)

lat = 26.614607620070455
alt = 3

tspan = (0.0, 7000.0);

include("AO2.jl")

solfm = sol.t
adult_sum_fm = adult_sum

dt = DateTime(2000,1,1)
  min_date = DateTime(2008,1,1)
  max_date = DateTime(2019,1,1)
  time_axis = add_float2datetime.(solfm,dt)



gfmbg = plot(time_axis,movingaverage(adult_sum_fm ,100) ./60, label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Predicted oviposition activity")
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(DateTime.(FM_BG[:,2]), movingaverage(FM_BG[:,3],3)    ,   label = nothing, seriestype =:scatter, color=:orange  )
                plot!(DateTime.(FM_BG[:,2]), movingaverage(FM_BG[:,3],3)   ,   label = nothing, color=:orange )
                plot!(ylims = (0,25))

#Indianapolous - VB

lat = 39.7684
tspan = (0.0, 5661.0);
alt = 219

Ind_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Ind_Clim.csv", header=true,DataFrame)

Ind_Temp = Ind_Clim[:,3]
Ind_Prec = Ind_Clim[:,4]
Ind_EV =   Ind_Clim[:,5]
Ind_DP =   Ind_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,5661,length = 5661)),Ind_Prec[2193:7853] )
Temps_spl  = Spline1D(collect(range(1,5661,length =  5661)),Ind_Temp[2193:7853])
DP_spl     = Spline1D(collect(range(1,5661,length = 5661)),Ind_DP[2193:7853])
EV_spl     = Spline1D(collect(range(1,5661,length = 5661)),Ind_EV[2193:7853])

Ind_ad =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\ind_ad_rev1.csv", header=true,DataFrame;  dateformat = "yyyy-mm-dd")

include("AO2.jl")

solind= sol.t
adult_sum_ind = adult_sum

dt = DateTime(2006,1,1)
  min_date = DateTime(2018,1,1)
  max_date = DateTime(2021,1,1)
  time_axis = add_float2datetime.(solind,dt)

Ind2 = Ind_ad[72:91,2:3]

Temp_ind = zeros(20)

for i in 1:20
  Temp_ind[i] = (Ind2[i,1] .-  Date(2006,1,1)).value
end

Ind3 = [Temp_ind    Ind2[:,2]]

ind_fit = RMSE_fun(movingaverage(adult_sum_ind ,100),solind,Ind3,14, 0,0.01,0.0001)

gind = plot(time_axis,movingaverage(adult_sum_ind ,100)  .* ind_fit[3]  , ylabel = "Number of adults", xlabel = "Time (Days)", label = "Predicted adults")
        plot!(DateTime.(Ind_ad[72:77,2] ) ,movingaverage(Ind_ad[72:77,3],1),  label = "Observed adults")
        plot!(DateTime.(Ind_ad[72:77,2] ),movingaverage(Ind_ad[72:77,3],1),   label = nothing, seriestype =:scatter, color=:orange )
        plot!(DateTime.(Ind_ad[78:84,2] ) ,movingaverage(Ind_ad[78:84,3],1),  label = nothing , color = pal[2])
        plot!(DateTime.(Ind_ad[78:84,2] ),movingaverage(Ind_ad[78:84,3],1),   label = nothing, seriestype =:scatter, color=:orange )
        plot!(DateTime.(Ind_ad[85:91,2] ) ,movingaverage(Ind_ad[85:91,3],1),  label = nothing , color = pal[2])
        plot!(DateTime.(Ind_ad[85:91,2] ),movingaverage(Ind_ad[85:91,3],1),   label = nothing, seriestype =:scatter, color=:orange )
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylims = (-0.018,0.6))


# Colombus - VB

lat = 39.9612
tspan = (0.0, 2739.0);
alt = 275

Col_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Col_Clim.csv", header=true,DataFrame)

Col_Temp = Col_Clim[:,3]
Col_Prec = Col_Clim[:,4]
Col_EV =   Col_Clim[:,5]
Col_DP =   Col_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,2739,length = 2739)),Col_Prec[5115:7853] )
Temps_spl  = Spline1D(collect(range(1,2739,length =  2739)),Col_Temp[5115:7853])
DP_spl     = Spline1D(collect(range(1,2739,length = 2739)),Col_DP[5115:7853])
EV_spl     = Spline1D(collect(range(1,2739,length = 2739)),Col_EV[5115:7853])

Col_ad_bg =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Col_Ad_2.csv", header=true,DataFrame)

include("AO2.jl")

solcol= sol.t
adult_sum_col= adult_sum

dt = DateTime(2014,1,1)
  min_date = DateTime(2018,4,1)
  max_date = DateTime(2019,1,1)
  time_axis = add_float2datetime.(solcol,dt)
  time_axis2 = add_float2datetime.(1:2739,dt)

Temp_col= zeros(31)

for i in 1:31
  Temp_col[i] = (Col_ad_bg[i,2] .-  Date(2014,1,1)).value
end

Temp_col1 = [Temp_col[1:31] Col_ad_bg[1:31,3]]


col_fit = RMSE_fun(movingaverage(adult_sum_col ,1),solcol,Temp_col1[13:31,:],8, 0, 0.1,0.001)


bcbg = plot(time_axis,movingaverage(adult_sum_col ,1) .* col_fit[3], ylabel = "Number of adults", xlabel = "Time (Days)", label = "Predicted adults")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(DateTime.(Col_ad_bg[13:31,2]) ,movingaverage(Col_ad_bg[13:31,3],1),  label = "Observed adults")
            plot!(DateTime.(Col_ad_bg[13:31,2])  , Col_ad_bg[13:31,3]  ,   label = nothing, seriestype =:scatter, color=:orange  )

#VB - Norfolk DONE

lat = 36.8508
tspan = (0.0, 5661.0);
alt = 2.13

Nor_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Suf_Clim.csv", header=true,DataFrame)

Nor_Temp = Nor_Clim[:,3]
Nor_Prec = Nor_Clim[:,4]
Nor_EV =   Nor_Clim[:,5]
Nor_DP =   Nor_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,5661,length = 5661)),Nor_Prec[2193:7853] )
Temps_spl  = Spline1D(collect(range(1,5661,length =  5661)),Nor_Temp[2193:7853])
DP_spl     = Spline1D(collect(range(1,5661,length = 5661)),Nor_DP[2193:7853])
EV_spl     = Spline1D(collect(range(1,5661,length = 5661)),Nor_EV[2193:7853])

#Nor_ad_BG =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Nor_BGW.csv", header=true,DataFrame;  dateformat="dd/mm/yyyy")
Nor_ad_BG =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Nor_BGM.csv", header=true,DataFrame;  dateformat="dd/mm/yyyy")

include("AO2.jl")

solnor= sol.t
adult_sum_nor = adult_sum

dt = DateTime(2006,1,1)
    min_date = DateTime(2009,3,1)
    max_date = DateTime(2019,1,1)
    time_axis = add_float2datetime.(solnor,dt)

Temp_nor = zeros(71)

for i in 1:71
    Temp_nor[i] = (Nor_ad_BG[i,1] .-  Date(2006,1,1)).value
end

Temp_nor1 = [Temp_nor[1:71,1] Nor_ad_BG[1:71,2]]


nor_fit = RMSE_fun(movingaverage(adult_sum_nor ,100),solnor,Temp_nor1,8, 0, 0.2,0.001)


gnor1 = plot(time_axis,movingaverage(adult_sum_nor ,10) .* nor_fit[3]      , ylabel = "Number of adults", xlabel = "Time (Days)", label = "Predicted adults")
              plot!(xlims = Dates.value.([min_date,max_date]))
              plot!(add_float2datetime.(Temp_nor1[1:7,1] .+ nor_fit[2],dt), Nor_ad_BG[1:7,2], label = "Observed adults", color=pal[2] )
              plot!(add_float2datetime.(Temp_nor1[1:7,1] .+ nor_fit[2],dt )  , Nor_ad_BG[1:7,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[8:14,1] .+ nor_fit[2],dt), Nor_ad_BG[8:14,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[8:14,1] .+ nor_fit[2],dt )  , Nor_ad_BG[8:14,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[15:21,1] .+ nor_fit[2],dt), Nor_ad_BG[15:21,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[15:21,1] .+ nor_fit[2],dt )  , Nor_ad_BG[15:21,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[22:28,1] .+ nor_fit[2],dt), Nor_ad_BG[22:28,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[22:28,1] .+ nor_fit[2],dt )  , Nor_ad_BG[22:28,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[29:35,1] .+ nor_fit[2],dt), Nor_ad_BG[29:35,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[29:35,1] .+ nor_fit[2],dt )  , Nor_ad_BG[29:35,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[36:42,1] .+ nor_fit[2],dt), Nor_ad_BG[36:42,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[36:42,1] .+ nor_fit[2],dt )  , Nor_ad_BG[36:42,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[43:49,1] .+ nor_fit[2],dt), Nor_ad_BG[43:49,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[43:49,1] .+ nor_fit[2],dt )  , Nor_ad_BG[43:49,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[50:57,1] .+ nor_fit[2],dt), Nor_ad_BG[50:57,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[50:57,1] .+ nor_fit[2],dt )  , Nor_ad_BG[50:57,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[58:64,1] .+ nor_fit[2],dt), Nor_ad_BG[58:64,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[58:64,1] .+ nor_fit[2],dt )  , Nor_ad_BG[58:64,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(add_float2datetime.(Temp_nor1[65:71,1] .+ nor_fit[2],dt), Nor_ad_BG[65:71,2] ,   label = nothing, color=pal[2])
              plot!(add_float2datetime.(Temp_nor1[65:71,1] .+ nor_fit[2],dt )  , Nor_ad_BG[65:71,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(ylims = (-1,120))

## Greenville - VB DONE

lat = 34.8526
tspan = (0.0, 1165.0);
alt = 294

Gre_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Gre_Clim.csv", header=true,DataFrame)

Gre_Temp = Gre_Clim[:,3]
Gre_Prec = Gre_Clim[:,4]
Gre_EV =   Gre_Clim[:,5]
Gre_DP =   Gre_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1165,length = 1165)),Gre_Prec[5115:6279] )
Temps_spl  = Spline1D(collect(range(1,1165,length =  1165)),Gre_Temp[5115:6279])
DP_spl     = Spline1D(collect(range(1,1165,length = 1165)),Gre_DP[5115:6279])
EV_spl     = Spline1D(collect(range(1,1165,length = 1165)),Gre_EV[5115:6279])

Gre_ov =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Gre_O.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")

include("AO2.jl")

solgre= sol.t

gre_OV = Ovi_Obs

dt = DateTime(2014,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,1,1)
  time_axis = add_float2datetime.(1:1165,dt)

Temp_gre = zeros(29)

for i in 1:29
  Temp_gre[i] = (Gre_ov[i,1] .-  Date(2014,1,1)).value
end

Temp_gre1 = [Temp_gre[1:29,1] Gre_ov[1:29,2]]



gre_fit = RMSE_fun(movingaverage(gre_OV ,1),1:1165,Temp_gre1,7, 0, 0.01,0.00001)


ggre = plot(time_axis,movingaverage(gre_OV ,1) .* gre_fit[3]   , label = "Predicted oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
              plot!(xlims = Dates.value.([min_date,max_date]))
              plot!(add_float2datetime.(Temp_gre1[:,1].+ gre_fit[2] ,dt),Temp_gre1[:,2]  , label = "Observed oviposition activity")
              plot!(add_float2datetime.(Temp_gre1[:,1].+ gre_fit[2] ,dt )  , Temp_gre1[:,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(xlims = Dates.value.([min_date,max_date]))
              plot!(ylims = (-1,80))


#Fayetteville - VB

lat = 36.0662
tspan = (0.0, 800.0);
alt = 31

Fay_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Fay_Clim.csv", header=true,DataFrame)

Fay_Temp = Fay_Clim[:,3]
Fay_Prec = Fay_Clim[:,4]
Fay_EV =   Fay_Clim[:,5]
Fay_DP =   Fay_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,800,length = 800)),Fay_Prec[5480:6279] )
Temps_spl  = Spline1D(collect(range(1,800,length =  800)),Fay_Temp[5480:6279])
DP_spl     = Spline1D(collect(range(1,800,length = 800)),Fay_DP[5480:6279])
EV_spl     = Spline1D(collect(range(1,800,length = 800)),Fay_EV[5480:6279])

Fay_lar_ov =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Fay_Ovi.csv", header=true,DataFrame; dateformat="dd/mm/yyyy")

include("AO2.jl")

solfay= 1:800

Ovi_Fay = Ovi_Obs

dt = DateTime(2015,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,1,1)
  time_axis = add_float2datetime.(solfay,dt)

Temp_fay = zeros(22)

for i in 1:22
Temp_fay[i] = (Fay_lar_ov[i,1] .-  Date(2015,1,1)).value
end

Temp_fay1 = [Temp_fay[1:22,1] movingaverage(Fay_lar_ov[1:22,2],1)]


fay_fit = RMSE_fun(movingaverage(Ovi_Fay ,7),1:800,Temp_fay1,7, 0, 0.5,0.0005)


gfay = plot(time_axis,movingaverage(Ovi_Fay ,7) .* 0.001 , ylabel = "Oviposition activity", xlabel = "Time (Days)", label = "Predicted oviposition activity")
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(DateTime.(Fay_lar_ov[:,1]),Fay_lar_ov[:,2], label = "Observed oviposition activity")
                plot!(DateTime.(Fay_lar_ov[:,1] )  , Fay_lar_ov[:,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )


# Wilmington - VB

lat = 34.2104
tspan = (0.0, 800.0);
alt = 28

Wilm_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Wil_Clim.csv", header=true,DataFrame)

Wilm_Temp = Wilm_Clim[:,3]
Wilm_Prec = Wilm_Clim[:,4]
Wilm_EV =   Wilm_Clim[:,5]
Wilm_DP =   Wilm_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,800,length = 800)),Wilm_Prec[5480:6279])
Temps_spl  = Spline1D(collect(range(1,800,length = 800)),Wilm_Temp[5480:6279])
DP_spl     = Spline1D(collect(range(1,800,length = 800)),Wilm_DP[5480:6279])
EV_spl     = Spline1D(collect(range(1,800,length = 800)),Wilm_EV[5480:6279])

Wilm_lar_ov =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Wilm_L.csv", header=true,DataFrame)

include("AO2.jl")

solwilm= 1:800
egg_sum_wilm = Ovi_Obs
lar_sum_wilm = larv_sum

dt = DateTime(2015,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,1,1)
  time_axis = add_float2datetime.(solwilm,dt)

  Temp_wil = zeros(12)

  for i in 1:12
    Temp_wil[i] = (Wilm_lar_ov[i,2] .-  Date(2015,1,1)).value
  end

  Temp_wil1 = [Temp_wil[1:12] Wilm_lar_ov[1:12,3]]

wilm_fit = RMSE_fun(movingaverage(egg_sum_wilm ,7),solwilm,Temp_wil1,7, 0, 0.5,0.01)

gwilm = plot(time_axis,movingaverage(egg_sum_wilm ,7)  .* wilm_fit[3]  , ylabel = "Number of larvae", xlabel = "Time (Days)", label = "Predicted larvae")
                  plot!(xlims = Dates.value.([min_date,max_date]))
                  plot!(add_float2datetime.(Temp_wil1[:,1] .+ wilm_fit[3],dt ) ,Temp_wil1[:,2],label = "Observed larvae")
                  plot!(add_float2datetime.(Temp_wil1[:,1]  .+ wilm_fit[3],dt )   , Temp_wil1[:,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )

#Santa Rosa beach - VB DONE

lat = 38.4404
tspan = (0.0, 3470.0);
alt = 0.91

SR_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\San_Clim.csv", header=true,DataFrame)

SR_Temp = SR_Clim[:,3]
SR_Prec = SR_Clim[:,4]
SR_EV =   SR_Clim[:,5]
SR_DP =   SR_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,3470,length = 3470)),SR_Prec[4384:7853] )
Temps_spl  = Spline1D(collect(range(1,3470,length =  3470)),SR_Temp[4384:7853])
DP_spl     = Spline1D(collect(range(1,3470,length = 3470)),SR_DP[4384:7853])
EV_spl     = Spline1D(collect(range(1,3470,length = 3470)),SR_EV[4384:7853])

SR_ad =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\SR_AdW.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")
SR_ad2 =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\SR_AdM.csv", header=true,DataFrame; dateformat = "dd/mm/yyyy")

include("AO2.jl")


solsr= sol.t
adult_sum_sr = adult_sum

dt = DateTime(2012,1,1)
  min_date = DateTime(2015,1,1)
  max_date = DateTime(2018,1,1)
  time_axis = add_float2datetime.(solsr,dt)

Temp_sr = zeros(156)

for i in 1:156
    Temp_sr[i] = (SR_ad[i,1] .-  Date(2012,1,1)).value
end

Temp_sr1 = [Temp_sr[1:156,1] SR_ad[1:156,2]]


SR_fit = RMSE_fun(movingaverage(adult_sum_sr ,100),solsr,Temp_sr1[1:156,:],7, 0, 0.005,0.0001)

gsr = plot(time_axis,movingaverage(adult_sum_sr ,1) .* SR_fit[3]   , ylabel = "Number of adults", xlabel = "Time (Days)", label = "Predicted adults")
                plot!(add_float2datetime.(Temp_sr1[1:156,1] .+ SR_fit[2], dt),SR_ad[1:156,2] ,   label = "Observed adults")
                plot!(add_float2datetime.(Temp_sr1[1:156,1].+ SR_fit[2], dt )  , movingaverage(SR_ad[1:156,2],1)  ,   label = nothing, seriestype =:scatter, color=:orange  )
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims = (-0.1,4))




# St. Augustine

lat = 29.9012
tspan = (0.0, 7479.0);
alt = 0

Sta_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Sta_Clim.csv", header=true,DataFrame)

Sta_Temp = Sta_Clim[:,3]
Sta_Prec = Sta_Clim[:,4]
Sta_EV =   Sta_Clim[:,5]
Sta_DP =   Sta_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7479,length = 7479)),Sta_Prec[1:7479] )
Temps_spl  = Spline1D(collect(range(1,7479,length =  7479)),Sta_Temp[1:7479])
DP_spl     = Spline1D(collect(range(1,7479,length = 7479)),Sta_DP[1:7479])
EV_spl     = Spline1D(collect(range(1,7479,length = 7479)),Sta_EV[1:7479])

sta_ad_bg =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\sta_ad_bg.csv", header=true,DataFrame)
sta_ad_bl =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\sta_ad_bl.csv", header=true,DataFrame)

include("AO2.jl")

solsta= sol.t
adult_sum_sta = adult_sum

dt = DateTime(2000,1,1)
  min_date = DateTime(2016,1,1)
  max_date = DateTime(2017,1,1)
  time_axis = add_float2datetime.(solsta,dt)


 gsta = Plots.plot(time_axis,movingaverage(adult_sum_sta ,50) .* 0.1    , label = "Simulated eggs", xlabel = "Time (Days)", ylabel = "Number of eggs")
     plot!(DateTime.(sta_ad_bl[:,2] )   ,sta_ad_bl[:,3] , label ="Observed adults")
     plot!(DateTime.(sta_ad_bl[:,2] )  , sta_ad_bl[:,3]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(xlims = Dates.value.([min_date,max_date]))

#Windsor, Ontario


lat = 42.26
tspan = (0.0, 1461.0);
alt = 0

Ont_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Win_Clim.csv", header=true,DataFrame)
Ont_BGS = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Ont_BGS.csv", header=false,DataFrame)
Ont_O = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Ont_O.csv", header=false,DataFrame)

Ont_Temp = Ont_Clim[:,3]
Ont_Prec = Ont_Clim[:,4]
Ont_EV =   Ont_Clim[:,5]
Ont_DP =   Ont_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Ont_Prec[5480:6940] )
Temps_spl  = Spline1D(collect(range(1,1461,length =  1461)),Ont_Temp[5480:6940])
DP_spl     = Spline1D(collect(range(1,1461,length = 1461)),Ont_DP[5480:6940])
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)),Ont_EV[5480:6940])

include("AO2.jl")

solont= sol.t
adult_sum_ont = adult_sum
Ovi_Obs_ont = Ovi_Obs
dt = DateTime(2015,1,1)
  min_date = DateTime(2018,1,1)
  max_date = DateTime(2019,1,1)
  time_axis = add_float2datetime.(solont,dt)
  time_axis2 = add_float2datetime.(1:1461,dt)

ont_fit = RMSE_fun(movingaverage(adult_sum_ont ,100),solont,Ont_BGS,14, 1095, 20,0.1)
ont_fit2 = RMSE_fun(movingaverage(Ovi_Obs_ont ,7),1:1461,Ont_O,14, 1095, 1,0.01)


gont = Plots.plot(time_axis,movingaverage(adult_sum_ont ,50) .* 0.2     , label = "Simulated eggs", xlabel = "Time (Days)", ylabel = "Observed adults of eggs")
                plot!(add_float2datetime.(Ont_BGS[:,1] .+ 1095, dt),Ont_BGS[:,2] ,   label = "Observed adults")
                plot!(add_float2datetime.(Ont_BGS[:,1] .+ 1095 ,dt)  , Ont_BGS[:,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
              plot!(xlims = Dates.value.([min_date,max_date]))
              plot!(ylims = (-1,100))


gont_O = Plots.plot(time_axis2,movingaverage(Ovi_Obs_ont ,7) .* 0.01   , label = "Simulated eggs", xlabel = "Time (Days)", ylabel = "Observed adults of eggs")
          plot!(add_float2datetime.(Ont_O[:,1] .+ 1095, dt),Ont_O[:,2] ,   label = "Observed adults")
          plot!(add_float2datetime.(Ont_O[:,1] .+ 1095 ,dt)  , Ont_O[:,2]  ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-1,500))


##Asia####

#China##

# Khiri Rat Nikham - VB

lat = 28.9005
tspan = (0.0, 1461.0);
alt = 10

Khi_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Khi_Clim.csv", header=true,DataFrame)

Khi_Temp = Khi_Clim[:,3]
Khi_Prec = Khi_Clim[:,4]
Khi_EV =   Khi_Clim[:,5]
Khi_DP =   Khi_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Khi_Prec[6211:7671])
Temps_spl  = Spline1D(collect(range(1,1461,length =  1461)),Khi_Temp[6211:7671])
DP_spl     = Spline1D(collect(range(1,1461,length = 1461)),Khi_DP[6211:7671])
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)),Khi_EV[6211:7671])

khi_ad_man =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Khi_Ad.csv", header=true,DataFrame)

include("AO2.jl")

solkhir= sol.t
adult_sum_khir= adult_sum

dt = DateTime(2017,1,1)
  min_date = DateTime(2019,1,1)
  max_date = DateTime(2020,1,1)
  time_axis = add_float2datetime.(solkhir,dt)

Temp_khi = zeros(11)

for i in 1:11
  Temp_khi[i] = (khi_ad_man[i,2] .-  Date(2017,1,1)).value
end

Temp_khi1 = [Temp_khi[1:11,1] khi_ad_man[1:11,3]]

khi_fit = RMSE_fun(movingaverage(adult_sum_khir ,100),solkhir,Temp_khi1 ,21, 0, 1,0.001)


gkhir=Plots.plot(time_axis,movingaverage(adult_sum_khir ,1) .* khi_fit[3], label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
          plot!(DateTime.(khi_ad_man[:,2]  .+Day(khi_fit[2]))   ,khi_ad_man[:,3] ,label ="Observed adults")
          plot!(DateTime.(khi_ad_man[:,2] .+ Day(khi_fit[2]))  , khi_ad_man[:,3] ,   label = nothing, seriestype =:scatter, color=:orange  )
          plot!(xlims = Dates.value.([min_date,max_date]))
          plot!(ylims = (0,20))


#Gunagzahi - ZHang Jia Xia Xu DONE

lat = 23.1291
alt = 21
tspan = (0.0, 7671.0);

Guang_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Gua_Clim.csv", header=true,DataFrame)

Guang_Temp = Guang_Clim[:,3]
Guang_Prec = Guang_Clim[:,4]
Guang_EV =   Guang_Clim[:,5]
Guang_DP =   Guang_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7671,length = 7671)),Guang_Prec[1:7671] )
Temps_spl  = Spline1D(collect(range(1,7671,length =  7671)),Guang_Temp[1:7671])
DP_spl     = Spline1D(collect(range(1,7671,length = 7671)),Guang_DP[1:7671])
EV_spl     = Spline1D(collect(range(1,7671,length = 7671)),Guang_EV[1:7671])

Guang_Eggs =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Guang_Eggs_2016_2017.csv", header=false,DataFrame)
Guang_M_Adults =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Guang_More_Adults.csv", header=false,DataFrame)

include("AO2.jl")

solgua= sol.t
Egg_sum_gua = Ovi_Obs
adult_sum_gua = adult_sum
larv_sum_gua = larv_sum

gua_fit1 = RMSE_fun(movingaverage(Egg_sum_gua ,3),1:7671,Guang_Eggs,14, 6205 -365, 0.02,0.0001)

dt = DateTime(2000,1,1)
  min_date = DateTime(2016,6,1)
  max_date = DateTime(2018,1,1)
  time_axis = add_float2datetime.(solgua,dt)
  time_axis2 = add_float2datetime.(1:7671,dt)

ggua1 = Plots.plot(time_axis2,movingaverage(Egg_sum_gua ,7) .* gua_fit1[3]    , label = "Oviposition activity", xlabel = "Time (Days)", ylabel = "Simulated oviposition activity")
            plot!(add_float2datetime.(Guang_Eggs[:,1]  .+ gua_fit1[2] .+ 6205 .-365,dt)  ,Guang_Eggs[:,2], label ="Observed oviposition activity")
            plot!(add_float2datetime.(Guang_Eggs[:,1]  .+ gua_fit1[2].+ 6205 .-365,dt), Guang_Eggs[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-10,730))

min_date = DateTime(2004,1,1)
max_date = DateTime(2014,12,1)

gua_fit2 = RMSE_fun(movingaverage(adult_sum_gua ,100),solgua,Guang_M_Adults[40:110,:],16, 2193 .- 365, 0.007,0.0001)


ggua2 =Plots.plot(time_axis,movingaverage(adult_sum_gua ,1) .* gua_fit2[3], label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
        plot!(add_float2datetime.(Guang_M_Adults[1:110,1] .+ 2193 .- 365 ,dt)  ,Guang_M_Adults[1:110,2], label ="Observed adults")
        plot!(add_float2datetime.(Guang_M_Adults[1:110,1] .+ 2193 .- 365 ,dt), Guang_M_Adults[1:110,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))

Mort_Comp  = plot(time_axis, 1 ./ movingaverage(avg_mort,60), label = "Model prediction", ylabel = "Adult longeivity (Days)")
                plot!(time_axis, 1 ./ movingaverage(AL(Temp_sim),60), label = "Laboratory model")
                plot!(time_axis, 1 ./  movingaverage(AM(Temp_sim),60), label = "Field model")
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims =(-4.5,150))

#Japan ###

#Naha - - Toma

Naha_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Nah_Clim.csv", header=true,DataFrame)

Naha_Temp = Naha_Clim[:,3]
Naha_Prec = Naha_Clim[:,4]
Naha_EV =   Naha_Clim[:,5]
Naha_DP =   Naha_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,730,length = 730)),Naha_Prec[1:730])
Temps_spl  = Spline1D(collect(range(1,730,length =  730)),Naha_Temp[1:730])
DP_spl     = Spline1D(collect(range(1,730,length = 730)),Naha_DP[1:730])
EV_spl     = Spline1D(collect(range(1,730,length = 730)),Naha_EV[1:730])

Tom_Ad = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tom_Ad.csv", header=false,DataFrame)
Tomo_Larv =  CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tomo_Larv.csv", header=false,DataFrame);
Tomo_Egg = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tomo_Egg.csv", header=false,DataFrame)

lat = 26.2126
at = 3
tspan = (0.0, 730.0);

include("AO2.jl")

solnah= sol.t
Egg_sum_nah = Ovi_Obs
adult_sum_nah = adult_sum
larv_sum_nah = larv_sum

dt = DateTime(1978,1,1)
    min_date = DateTime(1978,3,1)
    max_date = DateTime(1979,6,1)
    time_axis = add_float2datetime.(solnah,dt)
    time_axis2 = add_float2datetime.(1:730,dt)

Tom_Ad2 = [0,57,112,34,93,97,73,72,21,30,14,5,22]

Tom_Ad3 = [Tom_Ad[1:13,1]  Tom_Ad2]

nah_fit1 = RMSE_fun(movingaverage(adult_sum_nah ,100),solnah,Tom_Ad  ,14, 0, 0.05,0.001)

gnahad= Plots.plot(time_axis,movingaverage(adult_sum_nah ,100) .* nah_fit1[3]  , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
                        plot!(xlims = Dates.value.([min_date,max_date]))
                        plot!(add_float2datetime.(Tom_Ad[2:13,1] .+ nah_fit1[2] ,dt), Tom_Ad[2:13,2] , label = "Observed adults")
                        plot!(add_float2datetime.(Tom_Ad[2:13,1] .+ nah_fit1[2] ,dt), Tom_Ad[2:13,2] ,   label = nothing, seriestype =:scatter, color=:orange  )

nah_fit2 = RMSE_fun(movingaverage(larv_sum_nah ,100),solnah,Tomo_Larv,14, 0, 0.1,0.001)

gnahlar= Plots.plot(time_axis,movingaverage(larv_sum_nah ,1)  .* nah_fit2[3]     , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
                        plot!(xlims = Dates.value.([min_date,max_date]))
                        plot!(add_float2datetime.(Tomo_Larv[:,1]  .+ nah_fit2[2],dt), Tomo_Larv[:,2]  , label = "Observed larvae")
                        plot!(add_float2datetime.(Tomo_Larv[:,1]  .+ nah_fit2[2],dt), Tomo_Larv[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )

nah_fit3 = RMSE_fun(movingaverage(Egg_sum_nah ,7),1:730,Tomo_Egg,14, 0, 0.04,0.0001)

gnahegg= Plots.plot(time_axis2,movingaverage(Egg_sum_nah ,7)  .* nah_fit3[3]    , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Number of eggs")
                        plot!(xlims = Dates.value.([min_date,max_date]))
                        plot!(add_float2datetime.(Tomo_Egg[:,1] .+ nah_fit3[2],dt), Tomo_Egg[:,2]  , label = "Observed  oviposition activity")
                        plot!(add_float2datetime.(Tomo_Egg[:,1] .+ nah_fit3[2],dt), Tomo_Egg[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )

##Nagasaki - - Suzuki DONE

Nag_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Nag_Clim.csv", header=true,DataFrame)
Nag_WL = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Nag_wl_NP.csv", header=false,DataFrame)
Nag_Ad = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Nag_Ad.csv", header=false,DataFrame)

Nag_Temp = Nag_Clim[:,3]
Nag_Prec = Nag_Clim[:,4]
Nag_EV =   Nag_Clim[:,5]
Nag_DP =   Nag_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1826,length = 1826)),Nag_Prec)
Temps_spl  = Spline1D(collect(range(1,1826,length =  1826)),Nag_Temp)
DP_spl     = Spline1D(collect(range(1,1826,length = 1826)),Nag_DP)
EV_spl     = Spline1D(collect(range(1,1826,length = 1826)),Nag_EV)

lat = 26.2126
alt  = 15
tspan = (0.0, 1826.0);


include("AO2.jl")

solnag= sol.t
Egg_sum_nag = Egg_sum
adult_sum_nag = adult_sum
avg_wing_nag = avg_wing
out_wing_nag = out_wing

dt = DateTime(1987,1,1)
    min_date = DateTime(1990,4,1)
    max_date = DateTime(1991,1,30)
    time_axis = add_float2datetime.(solnag,dt)

nag_fit = RMSE_fun(movingaverage(adult_sum_nag ,7),solnag,Nag_Ad,14, 1095, 0.2,0.001)


gnag= Plots.plot(time_axis,movingaverage(adult_sum_nag ,7) .* nag_fit[3]    , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of Adults")
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(add_float2datetime.(Nag_Ad[:,1]  .+ 1095 .+  nag_fit[2] ,dt) ,Nag_Ad[:,2], label = "Observed adults", color = pal[2] ,linewidth = 2)
        plot!(add_float2datetime.(Nag_Ad[:,1]  .+ 1095 .+  nag_fit[2] ,dt) ,Nag_Ad[:,2], label = nothing, color = :orange, seriestype =:scatter)
        plot!(ylims = (-1,90))


gnag_wl =  heatmap(time_axis,
        Wing_Vals , out_wing_nag ,
        c=cgrad(:oslo, categorical = true),
        xlabel="Time (Days)", ylabel="Wing length (unspecified units)", colorbar_title = "Proportion of adults", background_color_inside = "black")
        plot!(clims = (0,1))
        plot!(ylims = (1.5,4.5))
        plot!(time_axis,movingaverage(avg_wing_nag,1)  , label = "Simulated average wing length", color = :green ,linewidth = 2)
        plot!(add_float2datetime.(Nag_WL[:,1]  .-15 .+ 1095    ,dt) ,(Nag_WL[:,2] ) , label = "Average wing length of emerging adults", color = :orange ,linewidth = 2)
        plot!(add_float2datetime.(Nag_WL[:,1]  .-15 .+ 1095  ,dt) ,(Nag_WL[:,2] ), label = nothing, color = :orange, seriestype =:scatter)
        plot!(xlims = Dates.value.([min_date,max_date]))

#Tokyo - - Kori

Tok_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Tok_Clim.csv", header=true,DataFrame)
Tok_L= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tok_L.csv", header=false,DataFrame)
Tok_A= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tok_A.csv", header=false,DataFrame)
Tok_A2= CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Tok_Ad2.csv", header=false,DataFrame)

Tok_Temp = Tok_Clim[:,3]
Tok_Prec = Tok_Clim[:,4]
Tok_EV =   Tok_Clim[:,5]
Tok_DP =   Tok_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,7671,length = 7671)),Tok_Prec)
Temps_spl  = Spline1D(collect(range(1,7671,length =  7671)),Tok_Temp)
DP_spl     = Spline1D(collect(range(1,7671,length = 7671)),Tok_DP)
EV_spl     = Spline1D(collect(range(1,7671,length = 7671)),Tok_EV)

lat = 35.6762
at = 40
tspan = (0.0, 7671.0);

include("AO2.jl")

solTok= sol.t
Egg_sum_Tok = Egg_sum
adult_sum_Tok = adult_sum
larv_sum_Tok = larv_sum

dt = DateTime(2000,1,1)
    min_date = DateTime(2010,3,1)
    max_date = DateTime(2019,1,1)
    time_axis = add_float2datetime.(solTok,dt)

tok_A_fit = RMSE_fun(movingaverage(adult_sum_Tok ,100),solTok,Tok_A,14, 5480, 0.015,0.0001)
tok_A2_fit = RMSE_fun(movingaverage(adult_sum_Tok ,100),solTok,Tok_A2,14, 3654, 0.25,0.001)
tok_L_fit = RMSE_fun(movingaverage(larv_sum_Tok ,1),solTok,Tok_L,14, 5480, 0.004,0.0001)

gTok_A= Plots.plot(time_axis,movingaverage(adult_sum_Tok ,1)  .* tok_A_fit[3]   , label = "Predicted adults", xlabel = "Time (Days)", ylabel = "Number of Adults")
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(add_float2datetime.( (Tok_A[1:8,1] ).+ 5480 .+ tok_A_fit[2] ,dt), Tok_A[1:8,2]  , label = "Observed adults", color = pal[2])
        plot!(add_float2datetime.( (Tok_A[1:8,1] ).+ 5480 .+ tok_A_fit[2],dt), Tok_A[1:8,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(add_float2datetime.( (Tok_A[9:16,1] ).+ 5480 .+ tok_A_fit[2] ,dt), Tok_A[9:16,2]  , label = nothing, color = pal[2])
        plot!(add_float2datetime.( (Tok_A[9:16,1] ).+ 5480 .+ tok_A_fit[2],dt), Tok_A[9:16,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(add_float2datetime.( (Tok_A[17:24,1] ).+ 5480 .+ tok_A_fit[2] ,dt), Tok_A[17:24,2]  , label = nothing, color = pal[2])
        plot!(add_float2datetime.( (Tok_A[17:24,1] ).+ 5480 .+ tok_A_fit[2],dt), Tok_A[17:24,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(ylims = (-0.1,12))

gTok_A2= Plots.plot(time_axis,movingaverage(adult_sum_Tok ,1) .* tok_A2_fit[3]  , label = "Predicted Adults", xlabel = "Time (Days)", ylabel = "Number of Adults")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(add_float2datetime.( (Tok_A2[1:10,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[1:10,2]  , label = "Observed adults")
            plot!(add_float2datetime.( (Tok_A2[1:10,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[1:10,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.( (Tok_A2[11:18,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[11:18,2]  , label = nothing, color = pal[2])
            plot!(add_float2datetime.( (Tok_A2[11:18,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[11:18,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.( (Tok_A2[19:28,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[19:28,2]  , label = nothing, color = pal[2])
            plot!(add_float2datetime.( (Tok_A2[19:28,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[19:28,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.( (Tok_A2[29:37,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[29:37,2]  , label = nothing, color = pal[2])
            plot!(add_float2datetime.( (Tok_A2[29:37,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[29:37,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(add_float2datetime.( (Tok_A2[38:46,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[38:46,2]  , label = nothing, color = pal[2])
            plot!(add_float2datetime.( (Tok_A2[38:46,1] ).+ 3654 .+ tok_A2_fit[2],dt), Tok_A2[38:46,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )

gTok_L= Plots.plot(time_axis,movingaverage(larv_sum_Tok ,1) .* tok_L_fit[3]   , label = "Predicted larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(add_float2datetime.( (Tok_L[:,1] ).+ 5480 .+ tok_L_fit[2],dt), Tok_L[:,2]  , label = "Observed larvae")
        plot!(add_float2datetime.( (Tok_L[:,1] ).+ 5480 .+ tok_L_fit[2],dt), Tok_L[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )

#Penang, Malaysia - - Rozilwalta

Pen_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Pen_Clim.csv", header=true,DataFrame)


Pen_Temp = Pen_Clim[:,3]
Pen_Prec = Pen_Clim[:,4]
Pen_EV =   Pen_Clim[:,5]
Pen_DP =   Pen_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1827,length = 1827)),Pen_Prec[366:2192])
Temps_spl  = Spline1D(collect(range(1,1827,length =  1827)),Pen_Temp[366:2192])
DP_spl     = Spline1D(collect(range(1,1827,length = 1827)),Pen_DP[366:2192])
EV_spl     = Spline1D(collect(range(1,1827,length = 1827)),Pen_EV[366:2192])

Rozi_Eggs = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Rozi_Eggs.csv", header=false,DataFrame)

lat = 5.44
alt = 2
tspan = (0.0, 1827.0);

include("AO2.jl")

solpen= 1:1827
Eggs_pen = Ovi_Obs

dt = DateTime(2000,1,1)
        min_date = DateTime(2003,3,1)
        max_date = DateTime(2004,6,1)
        time_axis = add_float2datetime.(solpen,dt)

pen_fit = RMSE_fun(movingaverage(Eggs_pen ,7),solpen,Rozi_Eggs,7, 1095, 0.002,0.00001)


gpen= Plots.plot(time_axis,movingaverage(Eggs_pen ,7) .* pen_fit[3] , label = "Simulated Eggs", xlabel = "Time (Days)", ylabel = "Number of Adults")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(add_float2datetime.( (Rozi_Eggs[:,1] ).+ 1095 ,dt), Rozi_Eggs[:,2]   , label = "Observed Eggs")
            plot!(add_float2datetime.( (Rozi_Eggs[:,1] ).+ 1095 ,dt), Rozi_Eggs[:,2]   ,   label = nothing, seriestype =:scatter, color=:orange  )


#Suwon, Korea - - Hwang DONE

Suw_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\Suw_Clim.csv", header=true,DataFrame)
Suw_Eggs = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Suw_Ad_M.csv", header=true,DataFrame)


Suw_Temp = Suw_Clim[:,3]
Suw_Prec = Suw_Clim[:,4]
Suw_EV =   Suw_Clim[:,5]
Suw_DP =   Suw_Clim[:,6]

Precip_spl = Spline1D(collect(range(1,1461,length = 1461)),Suw_Prec[4750:6210])
Temps_spl  = Spline1D(collect(range(1,1461,length = 1461)),Suw_Temp[4750:6210])
DP_spl     = Spline1D(collect(range(1,1461,length = 1461)),Suw_DP[4750:6210])
EV_spl     = Spline1D(collect(range(1,1461,length = 1461)),Suw_EV[4750:6210])


lat = 37.278621373072795
alt = 2
tspan = (0.0, 1461.0);

include("AO2.jl")

solSuw= sol.t
Ad_Suw = adult_sum

dt = DateTime(2013,1,1)
        min_date = DateTime(2016,3,1)
        max_date = DateTime(2017,1,1)
        time_axis = add_float2datetime.(solSuw,dt)

hwa_fit = RMSE_fun(movingaverage(Ad_Suw ,100),solSuw,Suw_Eggs,14, 1095 , 0.1,0.001)

gSuw= Plots.plot(time_axis,movingaverage(Ad_Suw ,100) .* hwa_fit[3]   , label = "Simulated adults", xlabel = "Time (Days)", ylabel = "Number of adults")
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(add_float2datetime.(Suw_Eggs[:,1] .+ 1095 .+ hwa_fit[2],dt),movingaverage(Suw_Eggs[:,2],1)  , label = "Observed adults")
            plot!(add_float2datetime.(Suw_Eggs[:,1] .+ 1095 .+ hwa_fit[2],dt), movingaverage(Suw_Eggs[:,2],1)   ,   label = nothing, seriestype =:scatter, color=:orange  )


#Africa###

##La reunion - St. Paul

SP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\SP_Clim.csv", header=true,DataFrame)

SP_Temp = SP_Clim[:,3]
SP_Prec = SP_Clim[:,4]
SP_EV =   SP_Clim[:,5]
SP_DP =   SP_Clim[:,6]


Precip_spl = Spline1D(collect(range(1,4749,length = 4749)),SP_Prec)
Temps_spl  = Spline1D(collect(range(1,4749,length =  4749)),SP_Temp)
DP_spl     = Spline1D(collect(range(1,4749,length = 4749)),SP_DP)
EV_spl     = Spline1D(collect(range(1,4749,length = 4749)),SP_EV)

SP_Larv = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\RE_SP_Larv.csv", header=true,DataFrame)

lat = -21.012870920008957
alt = 4
tspan = (0.0, 4749.0);

include("AO2.jl")

solsp= sol.t
larv_sum_SP = larv_sum

dt = DateTime(2008,1,1)
        min_date = DateTime(2012,1,1)
        max_date = DateTime(2014,1,1)
        time_axis = add_float2datetime.(solsp,dt)

sp_fit = RMSE_fun(movingaverage(larv_sum_SP ,100),solsp,SP_Larv,7, 1460, 0.1,0.001)

gre_sp= Plots.plot(time_axis,movingaverage(larv_sum_SP ,100) .* sp_fit[3]  , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
          plot!(add_float2datetime.(SP_Larv[:,1] .+ 1460 .+  sp_fit[2] ,dt),movingaverage(SP_Larv[:,2],1)  , label = "Observed larvae")
          plot!(add_float2datetime.(SP_Larv[:,1] .+ 1460 .+  sp_fit[2],dt), movingaverage(SP_Larv[:,2],1)    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (0,80))

Mort_Comp  = plot(time_axis, 1 ./ movingaverage(avg_mort,60), label = "Model prediction", ylabel = "Adult longeivity (Days)")
                plot!(time_axis, 1 ./ movingaverage(AL(Temp_sim),60), label = "Laboratory model")
                plot!(time_axis, 1 ./  movingaverage(AM(Temp_sim),60), label = "Field model")
                plot!(xlims = Dates.value.([min_date,max_date]))
                plot!(ylims =(-4.5,150))


#la Possession

LP_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\LP_Clim.csv", header=true,DataFrame)

LP_Temp = LP_Clim[:,3]
LP_Prec = LP_Clim[:,4]
LP_EV =   LP_Clim[:,5]
LP_DP =   LP_Clim[:,6]


Precip_spl = Spline1D(collect(range(1,4749,length = 4749)),LP_Prec)
Temps_spl  = Spline1D(collect(range(1,4749,length =  4749)),LP_Temp)
DP_spl     = Spline1D(collect(range(1,4749,length = 4749)),LP_DP)
EV_spl     = Spline1D(collect(range(1,4749,length = 4749)),LP_EV)

LP_Larv = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\RE_LP_Larv.csv", header=true,DataFrame)

lat = -21.012870920008957
alt = 4
tspan = (0.0, 4749.0);

include("AO2.jl")

sollp= sol.t
larv_sum_lp = larv_sum

dt = DateTime(2008,1,1)
        min_date = DateTime(2012,1,1)
        max_date = DateTime(2014,1,1)
        time_axis = add_float2datetime.(sollp,dt)

LP_fit = RMSE_fun(movingaverage(larv_sum_lp ,100),sollp,LP_Larv,7, 1460, 0.1,0.001)

gre_lp= Plots.plot(time_axis,movingaverage(larv_sum_lp ,100) .*LP_fit[3]  , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
          plot!(add_float2datetime.(LP_Larv[:,1] .+ 1460 .+  LP_fit[2] ,dt),movingaverage(LP_Larv[:,2],1)  , label = "Observed larvae")
          plot!(add_float2datetime.(LP_Larv[:,1] .+ 1460 .+  LP_fit[2],dt), movingaverage(LP_Larv[:,2],1)    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (0,80))

#Saint-Bennoit

SB_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\SB_Clim.csv", header=true,DataFrame)

SB_Temp = SB_Clim[:,3]
SB_Prec = SB_Clim[:,4]
SB_EV =   SB_Clim[:,5]
SB_DP =   SB_Clim[:,6]


Precip_spl = Spline1D(collect(range(1,5114,length = 5114)),SB_Prec)
Temps_spl  = Spline1D(collect(range(1,5114,length =  5114)),SB_Temp)
DP_spl     = Spline1D(collect(range(1,5114,length = 5114)),SB_DP)
EV_spl     = Spline1D(collect(range(1,5114,length = 5114)),SB_EV)

SB_Larv = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\RE_SB_Larv.csv", header=true,DataFrame)

lat = -21.012870920008957
alt = 4
tspan = (0.0, 4749.0);

include("AO2.jl")

solSB= sol.t
larv_sum_SB = larv_sum

dt = DateTime(2008,1,1)
        min_date = DateTime(2011,1,1)
        max_date = DateTime(2013,1,1)
        time_axis = add_float2datetime.(solSB,dt)

SB_fit = RMSE_fun(movingaverage(larv_sum_SB ,40),solSB,SB_Larv,7, 1460, 0.31,0.001)

gre_sb= Plots.plot(time_axis,movingaverage(larv_sum_SB ,40) .* 0.01 , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
          plot!(add_float2datetime.(SB_Larv[:,1] .+ 1460 .+  SB_fit[2] ,dt),movingaverage(SB_Larv[:,2],1)  , label = "Observed larvae")
          plot!(add_float2datetime.(SB_Larv[:,1] .+ 1460 .+  SB_fit[2],dt), movingaverage(SB_Larv[:,2],1)    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (0,45))

#Sainte-Marie

SM_Clim = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\New_Clims\SM_Clim.csv", header=true,DataFrame)

SM_Temp = SM_Clim[:,3]
SM_Prec = SM_Clim[:,4]
SM_EV =   SM_Clim[:,5]
SM_DP =   SM_Clim[:,6]


Precip_spl = Spline1D(collect(range(1,4749,length = 4749)),SM_Prec)
Temps_spl  = Spline1D(collect(range(1,4749,length =  4749)),SM_Temp)
DP_spl     = Spline1D(collect(range(1,4749,length = 4749)),SM_DP)
EV_spl     = Spline1D(collect(range(1,4749,length = 4749)),SM_EV)

SM_Larv = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\RE_SM_Larv.csv", header=true,DataFrame)
SM_Ov = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Re_Ov.csv", header=false,DataFrame;  dateformat="dd/mm/yyyy")

lat = -21.012870920008957
alt = 4
tspan = (0.0, 4749.0);

include("AO2.jl")

solSM= sol.t
larv_sum_SM = larv_sum
SM_Obs = Ovi_Obs

dt = DateTime(2008,1,1)
        min_date = DateTime(2013,1,1)
        max_date = DateTime(2014,9,1)
        time_axis = add_float2datetime.(solSM,dt)
        time_axis1 = add_float2datetime.(1:4749,dt)

SM_fit = RMSE_fun(movingaverage(larv_sum_SM ,100),solSM,SM_Larv,7, 1460, 0.1,0.001)

gre_sm= Plots.plot(time_axis,movingaverage(larv_sum_SM ,100) .* SM_fit[3] , label = "Simulated larvae", xlabel = "Time (Days)", ylabel = "Number of larvae")
          plot!(add_float2datetime.(SM_Larv[:,1] .+ 1460 .+  SM_fit[2] ,dt),movingaverage(SM_Larv[:,2],1)  , label = "Observed larvae")
          plot!(add_float2datetime.(SM_Larv[:,1] .+ 1460 .+  SM_fit[2],dt), movingaverage(SM_Larv[:,2],1)    ,   label = nothing, seriestype =:scatter, color=:orange  )
            plot!(xlims = Dates.value.([min_date,max_date]))
            plot!(ylims = (-5.4,110))

Temp_OV = zeros(86)

for i in 1:86
  Temp_OV[i] = (SM_Ov[i,1] .-  Date(2008,1,1)).value
end

Temp_Ov1 = [Temp_OV[1:64,1] movingaverage(SM_Ov[1:64,2],1)]


SM_fit2 = RMSE_fun(movingaverage(SM_Obs ,1),1:4383,Temp_Ov1,7, 0, 0.1,0.001)

gre_sm2 = Plots.plot(time_axis1,movingaverage(SM_Obs ,1) .* SM_fit2[3]   , label = "Predicted oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")
      plot!(add_float2datetime.(Temp_Ov1[1:64,1]  .+  SM_fit2[2] ,dt),movingaverage(SM_Ov[1:64,2],1)  , label = "Observed oviposition activity")
      plot!(add_float2datetime.(Temp_Ov1[1:64,1]  .+  SM_fit2[2],dt), movingaverage(SM_Ov[1:64,2],1)    ,   label = nothing, seriestype =:scatter, color=:orange  )
        plot!(xlims = Dates.value.([min_date,max_date]))
        plot!(ylims = (-6.6,220))
