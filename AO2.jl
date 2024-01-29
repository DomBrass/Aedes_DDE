#Code to run the population model,

#Required inputs:

#Temp_spl(time)   : function that takes the day of the year and outputs the average temperature on that day
#Precip_spl(time) : function that takes the day of the year and outputs the total precipitation on that day
#EV_spl(time)     : function that takes the day of the year and outputs the total evaporation from open water bodies

#tspan            : unit step range defining the time period over which the simulation takes place, e.g. tspan = (0.0, 4000.0)
#lat              : latitude of the population

#Packages

using(DifferentialEquations)
using(Dierckx)
using(Plots)
using(CSV)
using(Interpolations)
using(QuadGK)
using(Statistics)
using(Dates)
using(DataFrames)

#Function to convert integer dates to date time format

function add_float2datetime(x,datetime)
    d = Int(floor(x))
    e = x-d
    h = Int(floor(e*24))
    e = e-h/24
    m = Int(floor(e*24*60))
    e = e-m/24/60
    s = Int(floor(e*24*60*60))
    e = e-s/24/60/60
    ms = Int(floor(e*24*60*60*1e3))
    return datetime+Day(d)+Hour(h)+Minute(m)+Second(s)+Millisecond(ms)
end

#Loads reaction norms

WL     = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\WL_re.csv", header=true, DataFrame)      #Wing length as predicted by larval temperature and density
A_Long = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\AMgam.csv", header=true, DataFrame)      #Adult longeivity as predicted by wing length and adult temperature
L_Surv = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\Fin_LSurv.csv", header=true, DataFrame)  #Larval survival as predicted by larval temperature and density
L_Dur  = CSV.read(raw"C:\Users\dombr\OneDrive - UKCEH\Dom Mosquito Project\Project\Data\LDgam.csv", header=true, DataFrame)      #Larval duration as predicted by larval temperature and density

#Misc parameters

nTemp = 200; nDens = 200; nWing = 64;     # Number of environmental classes
n     = nTemp * nDens;                    # Number of combinations of temp and density
ooo = 16                                  # Time taken to evaporate

#Container dimensions

H_width  = 114                            # Width of container
H_length = 114                            # Length of container
height   = 38                             # Height of container
surface_area = H_width * H_length         # Surface area of container
m_precip     = height * surface_area      # Volume of container

#Function that calculates moving averages

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

#Defines the temperature, requires "Temps_spl" which should be a function taking the
#day of the year and outputting temperature

function Temp(t)
    time = t   ;
    out = Temps_spl(time);
    return(out)
end

#Function that calculates precipitation accumulating in habitat at time t by scaling the precipitation per unit area to match the surface area of the larval habitat

function Precip_check(t)
    out = Precip_spl.(t) .* surface_area
    if out < 0
      out = 0
    end
    return out
end

#Defines a vector that describes the volume of water in habitat for each day of the simulation, accounting for evaporation, accumlation, and overspill

Water_Vect = repeat([0.0],Int(tspan[2]))    #Initialises precipitation vector
Water_Vect[1:ooo] = repeat([m_precip],ooo)  #Starting conditions

for i in (1 + ooo):Int(tspan[2])
  if Water_Vect[i-1] > 0
    Water_Vect[i] = Water_Vect[i-1] +    Precip_check(i)  .+  EV_spl(i) .* surface_area
  else
      Water_Vect[i] = Water_Vect[i-1] +  Precip_check(i)
  end
  if Water_Vect[i] .> m_precip
    Water_Vect[i] = m_precip
  end
end

#Defines spline that calculates volume of water in habitat

Water_spl = Spline1D(collect(range(1,Int(tspan[2]),length = Int(tspan[2]))),Water_Vect)

#Function that ensures the water level is never 0

function Water(t)
    out =Water_spl.( t)
    if out <0.1
      out = 0.1
    end
    return(( out ))
end

#Function defining the amount of food produced by the larval habitat each day through respiration and decomposition

function resp(TEMP,time)
    10 .^ (0.45 + 0.095 * TEMP)   .* Water.(time) .* 1e-6 .+ 80
end

#Function that releases eggs from quiescence upon inundation

function R(t)
      if Water(t) > Water(t-1)
        out  = 1 -  (m_precip  - Water(t))./(m_precip)
      else
        out = 0
      end

      if out < 0
        out = 0
      end
      if out > 1
        out = 1
      end
      return(out)
end;

#Intiates vector for storing derivatives and adult recruitment

derivs   = zeros(nWing + 10) ; R_A      = zeros(nWing) ;

#Defines numbers for referencing derivatives

numb_E_1   = Int(1)     ; numb_E_D = 2         ; numb_E_Q = 3 ;
numb_L_1   = 4          ;
numb_t_E   = 5          ; numb_P_E = 6         ;
numb_t_L   = 7          ; numb_P_L = 8         ;
numb_t_P   = 9          ; numb_P_P = 10        ;
numb_A_1   = 11         ; numb_A_2 = nWing + 10;

#Creates discretisastion of environmental cues and defines environmental classes

Temp_Disc    = collect(range(-15,40,length = nTemp+1)) ; Temp_Vals          = (Temp_Disc[1:nTemp].+Temp_Disc[2:nTemp+1])./2 ;
Temp_Disc[1] = -100                                  ; Temp_Disc[nTemp+1] =  100

Dens_Disc = (collect(range(-8,1,length = nDens+1))); Dens_Vals = ((Dens_Disc[1:nDens]+Dens_Disc[2:nDens+1])/2) ;
Dens_Disc[nDens+1] = exp(14); Dens_Vals[1] = -8; Dens_Disc[1] = -exp(14);

Wing_Disc    = collect(range(1.5,4,length = nWing+1)) ; Wing_Vals          = (Wing_Disc[1:nWing].+Wing_Disc[2:nWing+1])./2 ;
Wing_Disc[1] = 0                               ; Wing_Disc[nWing+1] =  4

### Adult parameters ###

#Defines adult wing lengths based on larval experience of temperature and density

T_Vals = unique(WL[:,3])
D_Vals = unique(WL[:,2])

Wing_W       = transpose(reshape(WL[:,4], (400,400))); nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,Wing_W,Gridded(Linear()));
wing_spl     = extrapolate(sGird,Flat());

function wing_func(Temp,Dens)
  ifelse( Temp > 15, wing_spl(Temp,Dens),  wing_spl(15,Dens))
end

#Defines the mortality of adults based on wing length and adult temperature

T_Vals = unique(A_Long[:,2])
W_Vals = unique(A_Long[:,3])

Long    = transpose(reshape(A_Long[:,4], (400,400))); nodes = (W_Vals,T_Vals); sGird = interpolate(nodes,Long,Gridded(Linear()));
del_spl = extrapolate(sGird,Flat());

function delta_AT(Temp)

  del = del_spl.(Wing_Vals,Temp)
  del =ifelse.(del .< 0.00001,0.00001,del)
  out = -log(0.5) ./ del
  out = ifelse.(out .> 0.5,  0.5, out)

  return(out)
end

#Defines the length of the gonotrophic cycle

function G_Time(Temp)

   if Temp > 38.3
      Temp = 38.3
   end
   if Temp < 12
     Temp =  12
   end
   out = 1 ./ (1.93e-04  * Temp * (Temp - 10.25 ) * (38.32 - Temp) ^ (1/2))

   return(out)
end

#Defines function for eggs per gonotrophic cycle based on wing length

function q_func(Temp)

    out =  exp.(2.35 .+ 0.69 .* Wing_Vals)./ 2

    return(out)
end;

#Proportion of eggs that are quiescent

function Qui(t)

    if Water(t) < Water(t-1)
      out = 1 .* (m_precip   - Water(t))./(m_precip )
      if out > 1
        out = 1
      end
      if out < 0
        out = 0
      end
    else
      out = 0
    end

    return(out)
end

###Egg parameters###

#Egg development rate

function g_E(Temp)

    if Temp < 6
      out =    -0.0008256     * 6 ^2  +0.0334072    * 6-0.0557825
    elseif Temp > 38
      out =    -0.0008256     * 38 ^2  +0.0334072    * 38 -0.0557825
    else
      out =    -0.0008256     * Temp ^2  +0.0334072    *Temp -0.0557825
    end

    return(out)
end;

#Length of egg stage

function tau_E(Temp)

    out = 1/g_E(Temp)

    return(out)
end;

#Through egg stage survival

function P_E(Temp)

    out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((Temp - 24.672 )/6.115) ^ 2)

    if Temp < 6
      out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((6 - 24.672 )/6.115) ^ 2)
    end
    if Temp > 40
      out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((40 - 24.672 )/6.115) ^ 2)
    end

    return(out)
end;

#Active egg mortality rate

function delta_E(Temp)

    del = -log.(P_E(Temp)) ./ tau_E(Temp)
    if del > 0.99
      del = 0.99
    end

    return(del)
end;

#Quiescent egg mortality rate

function delta_E_Q(Temp)

    out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((Temp - 24.672 )/6.115) ^ 2)

    if out <= 0.01
        out = 0.01
    end

    del = -log.(out) ./ tau_E(Temp)
    if del > 1
      del = 1
    end

    return(del)
end;

#Diapausing egg mortality rate

function delta_E_D(Temp)
    if Temp > -12
      out = 0.01
    else
      out = 0.1
    end
    return(out)
end;

###Pupal parameters###

#Pupal development rate

function g_P(TEMP)

    if TEMP > 40
        out = 2.916e-05  * 40 * (40- 1.008e+01 ) * (4.768e+01 - 40) ^ (1/8.317e-01 )
    elseif TEMP < 12.3
        out = 2.916e-05  *12.3 * (12.3- 1.008e+01 ) * (4.768e+01- 12.3) ^ (1/8.317e-01 )
    else
        out = 2.916e-05  * TEMP * (TEMP - 1.008e+01 ) * (4.768e+01 - TEMP) ^ (1/8.317e-01 )
    end

    return(out)
end;

#Pupal development time

function tau_P(TEMP)

  tau = 1/g_P(TEMP)

  return(tau)
end;

#Through pupal-stage survival

function P_P(TEMP)

    out =     -0.0070628      * TEMP ^2  +0.3331028    *TEMP -2.9878761
    if TEMP > 34
      out =  -0.0070628      *  34 ^2  +0.3331028    * 34 -2.9878761
    elseif TEMP < 12.3
      out =  -0.0070628      *  12.3 ^2  +0.3331028    * 12.3 -2.9878761
    end

    if out <= 0.000001
      out = 0.000001
    end
    if out > 0.99
      out = 0.99
    end
    return(out)
end;

#Pupal mortality rate

function delta_P(TEMP,time)

    del =   (- log.(P_P.(TEMP))./tau_P.(TEMP))

    if Precip_check.(time) .>=  m_precip  && del < 0.2 && Water(time) == m_precip
      del = 0.2
    end

    if del > 0.99
      del = 0.99
    end
    if Water(time) <= 6000
      del = 0.5
    end

    return(del)
end;

###Larval parameters###

#Larval development time

D_Vals = unique(L_Dur[:,2])
T_Vals = unique(L_Dur[:,3])
L_D       = transpose(reshape(L_Dur[:,4], (400,400))); nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,L_D, Gridded(Linear()));
tau_L_spl = extrapolate(sGird,Flat());

#Larval development rate

function g_L(TEMP,DENS)

      out = 1 ./tau_L_spl(TEMP,DENS);
      if TEMP < 12
        out = 1 ./ tau_L_spl(12, DENS)
      end

      return(out)
end;

#Through larval-stage survival

D_Vals = unique(L_Surv[:,3])
T_Vals = unique(L_Surv[:,2])
L_S    = (reshape(L_Surv[:,4], (400,400))); nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,L_S, Gridded(Linear()));
surv_L = extrapolate(sGird,Flat());

function P_L(Temp,Dens)

    out = surv_L(Temp,Dens)

    if out < 0.01
      out = 0.01
    end
    if out > 0.99
      out = 0.99
    end

    return(out)
end

#Larval mortality rate


#The biotic component of larval mortality, from temperature and copmetition
function delta_L_spl(Temp,Dens)

    s = P_L(Temp,Dens)
    del = (- log.(s)./tau_L_spl(Temp,Dens))
    if del > 0.99
      del = 0.99
    end

    return(del)
end

#Modification of the larval mortality rate to account for abiotic factors such as overspill, evaporation, and overcrowding
function delta_L(dens,temp,larv,time)

    mort = delta_L_spl(temp,dens)
    if Water(time) <= 6000 #evaporation
      mort = 0.5
    end

    out = mort .+    exp.( -   exp.( 1.0 .-  0.33.* abs.(larv + 1) ./ (Water.(time) ./ 1000))) #overcrowding

    if Precip_check.(time) .>=   m_precip  && out < 0.2 && Water(time) == m_precip #overspill induced flushing
      out = 0.2
    end
    if out > 0.99
      out = 0.99
    end
    if out < 0.01
      out = 0.01
    end


    return(out)
end ;

#Function for initial impulse

function on_off(t,a)

  impulse = zeros(length(a));

  for i = 1:length(a)
    if t <= (a[i] .+ 10) && t .>= (a[i] .- 10)
      impulse[i] = 1
    end
  end

  return(impulse)
end;

#Function for assigning adults to environmental classes

function w_func(TEMP,DENS,VECT)

  VECT_TEMP = VECT

  inp  = wing_func(TEMP,DENS)

  for z in 1:nWing
          if Wing_Disc[z] <= inp < Wing_Disc[z + 1]
              VECT_TEMP[z] = 1
          end
  end

  return(VECT_TEMP)
end;

#Function that calculates photoperiod at time t

function PP(T)

  L = lat
  EPS = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (T-183.5)))))
  NUM = sin(0.8333 *pi/180) + (sin(L*pi/180) * sin(EPS))
  DEN = cos(L*pi/180) * cos(EPS)
  DAYLIGHT = 24 - (24/pi) * acos(NUM / DEN)

  return(DAYLIGHT)
end

#Defines critical photoperiod

a = 1.5 * abs(lat)/14 + 12.5 - 28*1.5/14

#Function for production of diapausing eggs

function Dia(t)
  if PP.(t) < PP(t-1)
    if Temp(t) < 18
      out =  1- 1/(1 + 15*exp( (PP(t) - a)))
    else
      out = 1
    end
  else
    out = 1
  end
  return(out)
end

#Function for release from diapause

function H(t)
  if Temp(t) > 14.5
    if PP(t) > a
      if PP(t) > PP(t-1)
        out = 1
      else
        out = 0
      end
    else
      out = 0
    end
  else
    out = 0
  end
end;

#Defines intial conditions

y0           = zeros(nWing + 10) ; y0[numb_E_1] = 1; y0[numb_L_1] = 150; y0[numb_E_D] = 100; y0[numb_E_Q] = 1;
y0[numb_t_P] = tau_P(Temp(0))
y0[numb_t_L] = tau_L_spl(Temp(0),log(Water(0)./y0[numb_L_1]))
y0[numb_t_E] = tau_E(Temp(0))
y0[numb_P_P] = exp(-y0[numb_t_P] * delta_P(Temp(0),0))
y0[numb_P_L] = exp(-y0[numb_t_L] * delta_L(log(Water(0)./y0[numb_L_1]),Temp(0),y0[numb_L_1],0))
y0[numb_P_E] = exp(-y0[numb_t_E] * delta_E(Temp(0)))

lags    = y0          ;
w_A      = zeros(nWing)     ;
h(p, t) = lags;

#Function that defines the model

function aedes_model(dy,y,h,p,t)

  #Temperatures at developmental milestones
  TEMP_NOW = Temp(t              )[1] ;
  TEMP_E   = Temp(t - y[numb_t_E])[1] ;
  TEMP_L   = Temp(t - y[numb_t_L])[1] ;
  TEMP_P   = Temp(t - y[numb_t_P])[1] ;

  #State of variables at developmental milestones
  ylag_E   = h(p,t -  y[numb_t_E])    ;
  ylag_L   = h(p,t -  y[numb_t_L])    ;
  ylag_P   = h(p,t -  y[numb_t_P])    ;

  #Average temperature over the larval period
  t1 = t-y[numb_t_P]
  t2 = t-y[numb_t_P]-ylag_P[numb_t_L]

  TEMP_P_AVG = mean(Temp.(range(t2,t1 , length = 20)))

  #Temperatures at developmental milestones
  TEMP_E_L   = Temp(t - y[numb_t_L] - ylag_L[numb_t_E])[1] ;
  TEMP_L_P   = Temp(t - y[numb_t_P] - ylag_P[numb_t_L])[1] ;

  #State of variables at developmental milestones
  ylag_E_L   = h(p, t - y[numb_t_L] - ylag_L[numb_t_E])    ;
  ylag_L_P   = h(p, t - y[numb_t_P] - ylag_P[numb_t_L])    ;
  ylag_L_L   = h(p, t - y[numb_t_L] - ylag_L[numb_t_L])    ;

  #Temperatures at developmental milestones
  TEMP_E_L_P = Temp(t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E])[1] ;

  #State of variables at developmental milestones
  ylag_E_L_P = h( p, t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E])    ;

  v         = log(             ( resp(TEMP_NOW,t)) / abs(       y[numb_L_1] + 1)) ;
  v_L       = log(             (  resp(TEMP_L,t - y[numb_t_L])  ) / abs(  ylag_L[numb_L_1] + 1)) ;
  v_L_P     = log(             (  resp(TEMP_L_P, t - ylag_P[numb_t_L] - y[numb_t_P])) / abs(ylag_L_P[numb_L_1] + 1)) ;
  v_P       = log(             (  resp(TEMP_P, t - y[numb_t_P])  ) / abs(  ylag_P[numb_L_1] + 1)) ;
  val_P     = log( quadgk(t -> (  resp(Temp(t),t) ) / abs(  h(p,t)[numb_L_1] + 1)  , t - y[numb_t_P] - ylag_P[numb_t_L], t - y[numb_t_P], rtol = 1)[1] ./ ylag_P[numb_t_L] ) ;

  #Proportion of adults producing diapausing eggs
  Dia_Now   = Dia( t                                                       ) ;
  Dia_E     = Dia( t - y[numb_t_E]                                         ) ;
  Dia_E_L   = Dia( t - y[numb_t_L] - ylag_L[numb_t_E]                      ) ;
  Dia_E_L_P = Dia( t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E] ) ;

  #Defines transition functions
  w_A_temp  = zeros(nWing)                              ;
  w_A       = w_func(TEMP_P_AVG, val_P, w_A_temp) ;

  #Recruitement terms
  R_E       = Dia_Now .* sum(q_func(TEMP_NOW) .* y[numb_A_1:numb_A_2] ./ (G_Time(TEMP_NOW)))  .+ on_off(t,0)[1] ;
  R_E_lag   = Dia_E   .* ( g_E(TEMP_NOW) / g_E(TEMP_E) ) .* sum( q_func(TEMP_E)  .* ylag_E[numb_A_1:numb_A_2] ./ (G_Time(TEMP_E)) ) .* y[numb_P_E] .+
                         ( g_E(TEMP_NOW) / g_E(TEMP_E) ) .* on_off(t-y[numb_t_E],0)[1] .* y[numb_P_E]      ;

  R_E_D     = (1-Dia_Now) .* sum(q_func(TEMP_NOW)  .* y[numb_A_1:numb_A_2] ./ (G_Time(TEMP_NOW)))                ;
  R_E_D_lag =  H(t) * y[numb_E_D] ;

  R_E_Q     = Qui(t) .* (R_E_lag .+ R_E_D_lag)
  R_E_Q_mat =  R(t) * y[numb_E_Q]

  R_L     = (1-Qui(t)) * Dia_E .* (g_E(TEMP_NOW) /  g_E(TEMP_E) ) * sum(q_func(TEMP_E)  .* ylag_E[numb_A_1:numb_A_2] ./ (G_Time(TEMP_E))) .* y[numb_P_E]  .+
            (1-Qui(t)) * H(t) * y[numb_E_D] .+
             R(t) * y[numb_E_Q] .+
            (g_E(TEMP_NOW) /  g_E(TEMP_E) ) * on_off(t-y[numb_t_E],0)[1] .* y[numb_P_E]

  R_L_lag = (1-Qui(t - y[numb_t_L])) * Dia_E_L .* (g_E(TEMP_L ) / g_E(TEMP_E_L)) * (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * sum(q_func(TEMP_E_L)  .* ylag_E_L[numb_A_1:numb_A_2] ./ (G_Time(TEMP_E_L))) .* ylag_L[numb_P_E] .* y[numb_P_L] .+
            (1-Qui(t - y[numb_t_L])) *                                             (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L ))  * H(t - y[numb_t_L]) * ylag_L[numb_E_D] .* y[numb_P_L] .+
                                                                                                      (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) *  R(t - y[numb_t_L]) * ylag_L[numb_E_Q] .* y[numb_P_L] .+
                                                                     (g_E(TEMP_L ) / g_E(TEMP_E_L)) * (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * on_off(t-ylag_L[numb_t_E] - y[numb_t_L],0)[1].*  ylag_L[numb_P_E]  * y[numb_P_L]

  R_A   = (1-Qui(t - ylag_P[numb_t_L] - y[numb_t_P])) * Dia_E_L_P .* w_A .* (g_E(TEMP_L_P)/ g_E(TEMP_E_L_P)) * (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*             sum(q_func(TEMP_E_L_P)  .*   ylag_E_L_P[numb_A_1:numb_A_2] ./ (G_Time(TEMP_E_L_P))) .* ylag_L_P[numb_P_E] .* ylag_P[numb_P_L] .* y[numb_P_P] .+
          (1-Qui(t - ylag_P[numb_t_L] - y[numb_t_P]))*               w_A .*                                    (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*                H(t - ylag_P[numb_t_L] - y[numb_t_P]) .* ylag_L_P[numb_E_D] .* ylag_P[numb_P_L] .* y[numb_P_P] .+
                                                                     w_A .*                                    (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*                R(t - ylag_P[numb_t_L] - y[numb_t_P]) .* ylag_L_P[numb_E_Q] .* ylag_P[numb_P_L] .* y[numb_P_P] .+
                                                                     w_A .* (g_E(TEMP_L_P)/ g_E(TEMP_E_L_P)) * (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .* on_off(t-ylag_L_P[numb_t_E] - ylag_P[numb_t_L] - y[numb_t_P],repeat([0],nWing)) .* ylag_L_P[numb_P_E] * ylag_P[numb_P_L] * y[numb_P_P]

  #Differential equations to be solved
  dy[numb_E_1]          = R_E   - R_E_lag   -  delta_E(TEMP_NOW)            * y[numb_E_1]                                                                                   ;
  dy[numb_E_D]          = R_E_D - R_E_D_lag -  delta_E_D(TEMP_NOW)            * y[numb_E_D]                                                                                   ;
  dy[numb_E_Q]          = R_E_Q - R_E_Q_mat -  delta_E_Q(TEMP_NOW)            * y[numb_E_Q]                                                                                   ;
  dy[numb_L_1]          = R_L   - R_L_lag   -  delta_L(v,TEMP_NOW,y[numb_L_1],t)  * y[numb_L_1]                                                                                   ;
  dy[numb_A_1:numb_A_2] = R_A              .-  delta_AT(TEMP_NOW )         .* y[numb_A_1:numb_A_2]                                                                          ;
  dy[numb_t_E]          = 1 - (g_E(TEMP_NOW    ) ./ g_E(TEMP_E)     )                                                                                                         ;
  dy[numb_t_L]          = 1 - (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L))                                                                                                         ;
  dy[numb_t_P]          = 1 - (g_P(TEMP_NOW    ) ./ g_P(TEMP_P)     )                                                                                                         ;
  dy[numb_P_E]          = y[numb_P_E]   * ((g_E(TEMP_NOW)    * delta_E(TEMP_E    )) / g_E(TEMP_E)     - delta_E(TEMP_NOW))            ;
  dy[numb_P_L]          = y[numb_P_L]   * ((g_L(TEMP_NOW, v) * delta_L(v_L,TEMP_L, ylag_L[numb_L_1], t - y[numb_t_L])) / g_L(TEMP_L,v_L) - delta_L(v,TEMP_NOW, y[numb_L_1],t)) ;
  dy[numb_P_P]          = y[numb_P_P]   * ((g_P(TEMP_NOW)    * delta_P(TEMP_P      , t- y[numb_t_P]                                    )) / g_P(TEMP_P)     - delta_P(TEMP_NOW,t))              ;

print(t)
print(" ")
end;

p = ( n,Water, tau_E, tau_P, delta_E, delta_L, P_E, numb_L_1, numb_A_1, numb_A_2, nDens) #Defines named parameters
prob = DDEProblem(aedes_model,y0,h,tspan,p;
                  depedent_lags = ((y,p,t) ->  [y[numb_t_E],y[numb_t_L], y[numb_t_P], y[numb_t_L] + ylag_L[numb_t_E], y[numb_t_P] + ylag_P[numb_t_L], y[numb_t_L] + ylag_L[numb_t_L], y[numb_t_P] + ylag_P[numb_t_L] + ylag_L_P[numb_t_E]], ))  #Defines problem to be solved
alg = MethodOfSteps(RK4())                 #Algorithm used to solve problem

sol = solve(prob,alg, rtol = 1e-3, etol = 1e-3)  #Runs the model

#Processing of model outputs

out = hcat(sol.u...);
thing = out[numb_A_1:numb_A_2,:]
adult_sum = sum(thing, dims = 1)[1,:]
big_adult_sum = sum(thing[30:60,:], dims = 1)[1,:]
larv_sum = out[4,:]

out_wing = ( (thing ) .* Wing_Vals) ./ (sum((thing ) .* Wing_Vals , dims = 1)[1,:])'
avg_wing = sum((thing ) .* Wing_Vals, dims = 1)[1,:] ./ (adult_sum )
avg_wing1 = ifelse.(isnan.(avg_wing) == true, 0, avg_wing)
for i in 1:length(avg_wing)
  avg_wing[i] = ifelse(isnan(avg_wing[i]) == true, 0, avg_wing[i])
  avg_wing[i] = ifelse((adult_sum[i]) <= 0.2, NaN, avg_wing[i])
  out_wing[:,i] = ifelse((adult_sum[i]) <= 0.2, repeat([NaN],nWing), out_wing[:,i])
end

Egg_sum = out[1,:] .+ out[3,:]

Temp_sim = Temp.(sol.t)
Temp_sim = ifelse.(Temp_sim .> 37 ,37 ,Temp_sim)
temp_mort = 0.8692 .- 0.159 .*Temp_sim .+ 0.01116 .*Temp_sim.^2 .- 3.408 .*10. ^(.-4).*Temp_sim.^3 + 3.809 .* 10 .^(.-6.00) .* Temp_sim.^4

function bite(Temp)
   0.0043 .* Temp .+ 0.0943
end

function human_vector(Temp)
  ifelse.(Temp .<= 26.1, 0.0729.*Temp .- 0.9037,1.0)
end

function vector_human(Temp)
  if Temp > 37.6
    Temp = 37.6
  end
  0.001044.*Temp.*(Temp .- 12.286) .* (37.6.-Temp).^0.5
end;

function EIP(Temp)
  4 .+ exp.(5.15 .- 0.123 .* Temp)
end

VC_nondynamic = (bite(Temp_sim)  .^ 2 .* human_vector(Temp_sim) .* vector_human.(Temp_sim) .* exp.(.- temp_mort  .* EIP(Temp_sim)))./ temp_mort

Temp_sim = float(Temp.(sol.t));

Delta_sim = zeros(nWing,length(sol.t))
for i = 1:length(sol.t)
    Delta_sim[:,i] .= delta_AT(Temp_sim[i])
end;
mort = Delta_sim .* thing

avg_mort = (sum(Delta_sim.*thing,dims=1)[1,:]) ./adult_sum

for i in 1:length(avg_mort)
  avg_mort[i] = ifelse(isnan(avg_mort[i]) == true, 0, avg_mort[i])
  avg_mort[i] = ifelse((adult_sum[i]) <= 0.2, NaN, avg_mort[i])
end

VC = (bite(Temp_sim)'  .^ 2.0 .* human_vector(Temp_sim)' .* vector_human.(Temp_sim)' .* exp.(.-EIP(Temp_sim)' .* Delta_sim))./Delta_sim
VCS = (VC.*thing)./adult_sum'

VA = VectorOfArray(q_func.(Temp.(sol.t)))
q_store = convert(Array,VA)

Eggs_In = sum(q_store  .* out[numb_A_1:numb_A_2,:] ./ transpose((G_Time.(Temp.(sol.t) ))), dims = 1)[1,:]
Eggs_daily = Spline1D(sol.t,Eggs_In)
tau_E_daily = Spline1D(sol.t,out[5,:])

Ovi_Obs    = repeat([0.0],Int(tspan[2]))

sampling_period = 14

for i = 21:Int(tspan[2])
   if tau_E_daily(i) < sampling_period
     Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:i - tau_E_daily(i)))
   else
    Ovi_Obs[i] =   sum(Eggs_daily(1:i)) - sum(Eggs_daily(1:Int(i - sampling_period)))
   end
end

Ads_daily = Spline1D(sol.t,adult_sum)

Ad_Obs    = repeat([0.0],Int(tspan[2]))

sampling_period = 14

for i = 21:Int(tspan[2])
   if tau_E_daily(i) < sampling_period
    Ad_Obs[i] =   sum(Ads_daily(1:i)) - sum(Ads_daily(1:Int(i - sampling_period)))
   end
end

A_Temp           = [00,5,10.0,15,20.0,25.0,30.0,35.0,40];
Adult_Mort_Field = [1,12,22,35,40,41,40,21,1]
Adult_Mort_Lab   = [1,17,38,76,94,112,96,64,1]
delta_A_max      = Spline1D(A_Temp,Adult_Mort_Lab)
delta_A_min      = Spline1D(A_Temp,Adult_Mort_Field)

AL(Temp) = 1 ./ delta_A_max(Temp)
AM(Temp) = 1 ./ delta_A_min(Temp)

E =  plot(sol.t,Egg_sum  , label = "Total eggs", xlabel = "Time (Days)", ylabel = "Number of Eggs")
     plot!(sol.t, out[1,:], label = "Active eggs")
     plot!(sol.t, out[2,:], label = "Diapausing eggs")
     plot!(sol.t,out[3,:], label = "Quiescent eggs")

L = Plots.plot(sol, vars = (numb_L_1), legend = false, xlabel = "Time (Days)", ylabel = "Number of Larvae")

A = plot(sol, vars = (numb_A_1:numb_A_2), label = nothing, color = :inferno, line_z = Wing_Vals', xlabel = "Time (Days)", ylabel = "Number of Adults")
      Plots.plot!(sol.t,movingaverage(adult_sum ,1), label = "Simulated Adults", color = :black, xlabel = "Time (Days)", ylabel = "Number of Adults")

plot(sol.t, VCS' ,palette=:inferno, legend  = false)
plot(sol.t,movingaverage(sum(VCS, dims =1)[1,:],300), colour = :black, label = "Simulated rVC", xlabel = "Time (Days)", ylabel = "rVC")
  plot!(sol.t,movingaverage( VC_nondynamic, 300) , label = "Estimated rVC")

Temp_Liu = [10,12,14,16,18,20,22,24,25,26,27,28.2,29.3,30.2,31,31.5,32,32.3,32.5]
m_Liu = [1.32,1.5,1.43,1.27,1.13,1.04,1.01,1.03,1.05,1.08,1.11,1.15,1.16,1.15,1.12,1.09,1.05,1.02,1]

  m_spl = Spline1D(Temp_Liu,m_Liu)

plot(sol.t,movingaverage(sum(VCS , dims =1)[1,:] .* (adult_sum ./ 100),30), colour = :black, label = "Simulated rVC")
    plot!(sol.t,movingaverage( VC_nondynamic .*  (m_spl.(Temp_sim)), 30), label = "Estimated rVC")

plot(sol.t,movingaverage(avg_mort,60))
  plot!(sol.t, movingaverage(AL(Temp_sim),60))
  plot!(sol.t, movingaverage(AM(Temp_sim),60))

plot(sol.t, movingaverage(Eggs_In   ,7)       , label = "Simulated oviposition activity", xlabel = "Time (Days)", ylabel = "Oviposition activity")

plot(movingaverage(Ovi_Obs ,7))

W =  heatmap(sol.t,
      Wing_Vals, out_wing,
      c=cgrad(:oslo, categorical = true, scale = :exp),
      xlabel="Time (Days)", ylabel="Wing length (mm)", colorbar_title = "Proportion of adults", background_color_inside = "black")
      plot!(clims = (0,1))
      plot!(ylims = (1,4.5))
      plot!(sol.t,movingaverage(avg_wing,100), label = "Simulated average wing length", linewidth = 2)
      plot!(xlims = (0,730))
