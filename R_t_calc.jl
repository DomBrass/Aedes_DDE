#R_t_calc.jl
#Author: Dominic Brass (dombra@ceh.ac.uk)

#Code accompanying the paper, "Phenotypic plasticity in vector traits drives trends in global disease incidence".
#Predicts R_t for model output of SO2_public.jl.

R_0 = reshape(repeat(Float64[0], Int64(tspan[2]) * nWing) , (nWing,Int64(tspan[2])))
ADU = reshape(repeat(Float64[0], Int64(tspan[2]) * nWing) , (nWing,Int64(tspan[2])))
TIME_R_0 = Spline1D(sol.t, sol.t .- out[numb_EIP,:])
TAU_EIP = Spline1D(sol.t, out[numb_EIP,:])

for j in 1:nWing
    global store_ad = Spline1D(sol.t,ad_out[j,:])
    global store_del = Spline1D(sol.t, 1 ./ Delta_sim[j,:])
    global store_surv = Spline1D(sol.t, out[numb_P_EIP_1 + j - 1,:])
    global store_H_S = Spline1D(sol.t, out[numb_H_S,:])
    global store_H_I = Spline1D(sol.t, out[numb_H_I,:])
    global store_H_R = Spline1D(sol.t, out[numb_H_R,:])
    global store_A_S = Spline1D(sol.t, sum(ad_out .* out[numb_P_EIP_2:numb_P_EIP_2,:], dims = 1)[1,:])

    times = collect(range(1,tspan[2],length = Int64(tspan[2])))
    
    print(j)
    print(" ")

    for i in 1:Int64(tspan[2])
        low = Float64(times[i] ) 
        high = Float64(times[i] + 4)
         R_0[j,i] = quadgk(t -> 40 * 4 * ((1/EIP(Temp(t))) ./ (1/EIP(Temp(t - TAU_EIP(t))))) * human_vector(Temp(t - TAU_EIP(t))) * bite(Temp(t - TAU_EIP(t))) * store_ad(t - TAU_EIP(t)) * store_surv(t) * (1 ./  (y0[numb_H_S]  + 15000 ))
         *  quadgk(u -> bite(Temp(u)) * vector_human(Temp(u)) * (y0[numb_H_S] ./  (y0[numb_H_S]  + 15000 )),  t, t +  store_del(t), rtol = 1)[1], low, high, rtol = 1)[1]
         ADU[j,i] = store_ad(i)


    end
end

R_0_spl = Spline1D(TIME_R_0.(collect(range(1,Int(tspan[2]),length = Int(tspan[2])))), sum(R_0',dims = 2)[:,1])
