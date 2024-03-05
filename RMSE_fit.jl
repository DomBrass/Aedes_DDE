#RMSE_fit.jl
#Author: Dominic Brass (dombra@ceh.ac.uk)

#Code accompanying the paper, "Phenotypic plasticity in vector traits drives trends in global disease incidence".
#Compares model predictions to field observations.

function RMSE_fun(pred,time,dat1,period, year, max_sca,step_sca)

    time_adjust = collect(range(-period,period,step = 1))
    scaling_factor = collect(range(0,max_sca,step = step_sca))
    dat = dat1
    time_points = dat[:,1] .+ year

    nmb_pnt = length(dat[:,1])

    RMSE = zeros(length(time_adjust),length(scaling_factor))
    R2 = zeros(length(time_adjust),length(scaling_factor))
    sqr_error = zeros(nmb_pnt)
    points = zeros(nmb_pnt)
    ss_tot = zeros(nmb_pnt)


    for i in 1:length(time_adjust)
        points = time_points .+ time_adjust[i]
        for j in 1:length(scaling_factor)
            for z in 1:nmb_pnt
                sol_loc = argmin(abs.(points[z] .- time))
                predict = pred[sol_loc]
                sqr_error[z] = ( dat[z,2] - scaling_factor[j] * predict)^2
                ss_tot[z] = (dat[z,2]- mean(dat[:,2]))^2
            end
            R2[i,j] = 1 - sum(sqr_error) / sum(ss_tot)
            RMSE[i,j] = mean(sqr_error)^0.5
        end
    end

    RMSE_min = findmin(RMSE)
    R2_max = findmax(R2)


    return(RMSE[R2_max[2]]/mean(dat[:,2]), time_adjust[R2_max[2][1]],scaling_factor[R2_max[2][2]], RMSE[R2_max[2]],R2[R2_max[2]])
end
