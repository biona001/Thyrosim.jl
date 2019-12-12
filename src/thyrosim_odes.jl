"""
Initial conditions for original thyrosim model
"""
function initialize_original_thyrosim(
    dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88];
    scale_Vp::Bool = true,
    height=1.77,
    weight=70,
    sex=true #true = male, false = female
    )

    # Initial conditions
    ic    = zeros(Float64, 19)
    ic[1] = 0.322114215761171
    ic[2] = 0.201296960359917
    ic[3] = 0.638967411907560
    ic[4] = 0.00663104034826483
    ic[5] = 0.0112595761822961
    ic[6] = 0.0652960640300348
    ic[7] = 1.78829584764370
    ic[8] = 7.05727560072869
    ic[9] = 7.05714474742141
    ic[10] = 0
    ic[11] = 0
    ic[12] = 0
    ic[13] = 0
    ic[14] = 3.34289716182018
    ic[15] = 3.69277248068433
    ic[16] = 3.87942133769244
    ic[17] = 3.90061903207543
    ic[18] = 3.77875734283571
    ic[19] = 3.55364471589659

    # Parameter values
    p = zeros(Float64, 48)
    p[1] = 0.00174155 * dial[1]     #S4
    p[2] = 8              #tau
    p[3] = 0.868          #k12
    p[4] = 0.108          #k13
    p[5] = 584            #k31free
    p[6] = 1503           #k21free
    p[7] = 0.000289       #A
    p[8] = 0.000214       #B
    p[9] = 0.000128       #C
    p[10] = -8.83*10^-6    #D
    p[11] = 0.88           #k4absorb; originally 0.881
    p[12] = 0.0189         #k02
    p[13] = 0.00998996     #VmaxD1fast
    p[14] = 2.85           #KmD1fast
    p[15] = 6.63*10^-4     #VmaxD1slow
    p[16] = 95             #KmD1slow
    p[17] = 0.00074619     #VmaxD2slow
    p[18] = 0.075          #KmD2slow
    p[19] = 3.3572*10^-4 * dial[3]   #S3
    p[20] = 5.37           #k45
    p[21] = 0.0689         #k46
    p[22] = 127            #k64free
    p[23] = 2043           #k54free
    p[24] = 0.00395        #a
    p[25] = 0.00185        #b
    p[26] = 0.00061        #c
    p[27] = -0.000505      #d
    p[28] = 0.88           #k3absorb; originally 0.882
    p[29] = 0.207          #k05
    p[30] = 1166           #Bzero
    p[31] = 581            #Azero
    p[32] = 2.37           #Amax
    p[33] = -3.71          #phi
    p[34] = 0.53           #kdegTSH-HYPO
    p[35] = 0.037          #VmaxTSH
    p[36] = 23             #K50TSH
    p[37] = 0.118          #k3
    p[38] = 0.29           #T4P-EU
    p[39] = 0.006          #T3P-EU
    p[40] = 0.037          #KdegT3B
    p[41] = 0.0034         #KLAG-HYPO
    p[42] = 5              #KLAG
    p[43] = 1.3            #k4dissolve
    p[44] = 0.12 * dial[2] #k4excrete; originally 0.119 (change with dial 2)
    p[45] = 1.78           #k3dissolve
    p[46] = 0.12 * dial[4] #k3excrete; originally 0.118 (change with dial 4)
    p[47] = 3.2            #Vp
    p[48] = 5.2            #VTSH

    if scale_Vp
        Vp, Vtsh = plasma_volume(height, weight, sex)
        p[47] = Vp
        p[48] = Vtsh
    end

    return ic, p
end

"""
initial conditions for new thyrosim model
"""
function initialize(
    dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88],
    scale_Vp::Bool = true,
    height=1.77,
    weight=70,
    sex=true #true = male, false = female
    )

    # TODO: need to calculate initial steady state
    ic    = zeros(Float64, 19)
    ic[1] = 0.322114215761171
    ic[2] = 0.201296960359917
    ic[3] = 0.638967411907560
    ic[4] = 0.00663104034826483
    ic[5] = 0.0112595761822961
    ic[6] = 0.0652960640300348
    ic[7] = 1.78829584764370
    ic[8] = 7.05727560072869
    ic[9] = 7.05714474742141
    ic[10] = 0
    ic[11] = 0
    ic[12] = 0
    ic[13] = 0
    ic[14] = 3.34289716182018
    ic[15] = 3.69277248068433
    ic[16] = 3.87942133769244
    ic[17] = 3.90061903207543
    ic[18] = 3.77875734283571
    ic[19] = 3.55364471589659

    # Parameter values
    p = zeros(Float64, 60)
    p[1] = 0.00174155 * dial[1]     #S4
    p[2] = 8               #tau
    p[3] = 0.868           #k12
    p[4] = 0.108           #k13
    p[5] = 584             #k31free
    p[6] = 1503            #k21free
    p[7] = 0.000289        #A
    p[8] = 0.000214        #B
    p[9] = 0.000128        #C
    p[10] = -8.83*10^-6    #D
    p[11] = 0.88           #k4absorb; originally 0.881
    p[12] = 0.0189         #k02
    p[13] = 0.00998996     #VmaxD1fast
    p[14] = 2.85           #KmD1fast
    p[15] = 6.63*10^-4     #VmaxD1slow
    p[16] = 95             #KmD1slow
    p[17] = 0.00074619     #VmaxD2slow
    p[18] = 0.075          #KmD2slow
    p[19] = 3.3572*10^-4 * dial[3]   #S3
    p[20] = 5.37           #k45
    p[21] = 0.0689         #k46
    p[22] = 127            #k64free
    p[23] = 2043           #k54free
    p[24] = 0.00395        #a
    p[25] = 0.00185        #b
    p[26] = 0.00061        #c
    p[27] = -0.000505      #d
    p[28] = 0.4978                #k3absorb; fitted to jonklaas
    p[29] = 0.207          #k05
    p[30] = 101                   #Bzero; fitted to blakesley
    p[31] = 47.64                 #Azero; fitted to blakesley
    p[32] = 0                     #Amax;  should be around 0 because 1976 weeke says hypothyroid patients should have no oscillations.
    p[33] = -3.71          #phi
    p[34] = 0.53           #kdegTSH-HYPO
    p[35] = 0.226                 #VmaxTSH; originally it's 0.037 but this is probably a typo because eq4 of 2010 eigenberg it not a real hill function
    p[36] = 23             #K50TSH
    p[37] = 0.118          #k3
    p[38] = 0.29           #T4P-EU
    p[39] = 0.006          #T3P-EU
    p[40] = 0.037          #KdegT3B
    p[41] = 0.0034         #KLAG-HYPO
    p[42] = 5              #KLAG
    p[43] = 1.3            #k4dissolve
    p[44] = 0.12 * dial[2] #k4excrete; originally 0.119 (change with dial 2)
    p[45] = 1.78           #k3dissolve
    p[46] = 0.12 * dial[4] #k3excrete; originally 0.118 (change with dial 4)
    p[47] = 3.2            #Vp
    p[48] = 5.2            #VTSH

    #parameters for hill functions in f_circ and SRtsh
    p[49] = 4.57           #K_circ -> fitted to Blakesley data (this will be recalculated in the ODE equations
    p[50] = 3.90           #K_SR_tsh -> fitted to Blakesley data (this will be recalculated in the ODE equations
    p[51] = 11.0 #6.91     #hill exponent in f_circ
    p[52] = 5.0  #7.66     #hill exponent in SR_tsh
    p[53] = 3.5            #Km for f4
    p[54] = 8.0            #hill exponent for f4

    if scale_Vp
        Vp, Vtsh = plasma_volume(height, weight, sex)
        p[47] = Vp
        p[48] = Vtsh
    end

    #TODO: setup p[60] = phase

    return ic, p
end

"""
    plasma_volume(height, weight, sex)

# Parameters used to get reference plasma volume (Vp) values:
## MCL: NEED TO DOUBLE-CHECK HEIGHT/WEIGHT
## Blakesley data: half male, half female all of "normal weight and height" (but no values given in paper).
## height: Average height in USA.
## weight: I think we used approximate values from back-transforming it from BMI = 22.5?
    male_height   = 1.77
    female_height = 1.63
    male_weight   = 70.0
    female_weight = 59.0

# Inputs
+ `height`: measured in meters
+ `weight`: measured in KG 
+ `sex`: 1 = male, 0 = female

# Outputs 
+ `Vp_new`: Scaled plasma volume (liters)
+ `Vtsh_new`: Scaled TSH distribution volume (liters)
"""
function plasma_volume(h, w, sex::Bool)
    Hem = 0.40 + 0.05 * sex #.45 for male and .4 for females (by default)
    BMI = w / h^2

    # some Vp reference volume. 
    male_ref_vp   = 2.932691217834299
    female_ref_vp = 2.5137479938183125

    # calculate Ideal Weight fitted to Feldschush's data
    if sex == 1
        iw = 176.3 - 220.6 * h + 93.5 * h^2
    elseif sex == 0
        iw = 145.8 - 182.7 * h + 79.55 * h^2
    end

    # power law fitted to Feldchush data
    a, n = 1.26975706e+03, 3.72981228e-01
    Δiw = (w - iw) / iw * 100  #deviation from ideal weight, in percentage
    Vb_per_kg = a * (100.0 + Δiw)^(n - 1)
    Vb = Vb_per_kg * w / 1000
    
    # calculate new Vp
    Vp_new = Vb*(1-Hem);       #Vp (orginally 3.2)
    if sex == 1
        Vp_new = Vp_new * 3.2 / male_ref_vp;    
    else
        Vp_new = Vp_new * 3.2 / female_ref_vp;   
    end

    # scale Vtsh according to new Vp
    Vtsh_new = 5.2 + Vp_new - 3.2 # Vtsh_old + (Vp_new - Vp_old) 

    return Vp_new, Vtsh_new
end

"""
Original thyrosim ODEs.

Source: https://bitbucket.org/DistefanoLab/thyrosim/src/master/resource/matlab/thyrosim_core.m
"""
function original_thyrosim(dq, q, p, t)
    kdelay = 5/8

    # Auxillary equations
    q4F = (p[24]+ p[25] * q[1] + p[26] * q[1]^2 + p[27] *q[1]^3) * q[4] #FT3p
    q1F = (p[7] + p[8] * q[1] + p[9] * q[1]^2 + p[10] * q[1]^3) * q[1]  #FT4p
    SR3 = (p[19] * q[19])                                        #Brain delay (dial 3)
    SR4 = (p[1] * q[19])                                         #Brain delay (dial 1)
    
    fCIRC = 1 + (p[32] / (p[31] * exp(-q[9])) - 1) * (1 / (1 + exp(10*q[9] - 55)))
    SRTSH = (p[30] + p[31] * fCIRC * sin(pi/12 * t - p[33])) * exp(-q[9])
    fdegTSH = p[34] + p[35] / (p[36] + q[7])
    fLAG = p[41] + 2*q[8]^11 / (p[42]^11 + q[8]^11)
    f4 = p[37] + 5 * p[37] / (1 + exp(2 * q[8] - 7))
    NL = p[13] / (p[14] + q[2])

    # ODEs
    dq[1]  = SR4 + p[3] * q[2] + p[4] * q[3] - (p[5] + p[6]) * q1F + p[11] * q[11] #T4dot (need to remove u1)
    dq[2]  = p[6] * q1F - (p[3] + p[12] + NL) * q[2]                                    #T4fast
    dq[3]  = p[5] * q1F -(p[4] + p[15] / (p[16] + q[3]) + p[17] /(p[18] + q[3])) *q[3]  #T4slow
    dq[4]  = SR3 + p[20] * q[5] + p[21] * q[6] - (p[22] + p[23]) * q4F + p[28] * q[13]  #T3pdot
    dq[5]  = p[23] * q4F + NL * q[2] - (p[20] + p[29]) * q[5]                         #T3fast
    dq[6]  = p[22] * q4F + p[15] * q[3] / (p[16] + q[3]) + p[17] * q[3] / (p[18] + q[3]) -(p[21])*q[6] #T3slow
    dq[7]  = SRTSH - fdegTSH * q[7]                                           #TSHp
    dq[8]  = f4 / p[38] * q[1] + p[37] / p[39] * q[4] - p[40] * q[8]          #T3B
    dq[9]  = fLAG * (q[8] - q[9])                                             #T3B LAG
    dq[10] = -p[43] * q[10]                                                   #T4PILLdot
    dq[11] =  p[43] * q[10] - (p[44] + p[11]) * q[11]                         #T4GUTdot
    dq[12] = -p[45] * q[12]                                                   #T3PILLdot
    dq[13] =  p[45] * q[12] - (p[46] + p[28]) * q[13]                         #T3GUTdot

    # Delay ODEs
    # dq[14] = -kdelay * q[14] + q[7]                                           #delay1 CHECK, might be wrong 
    dq[14] = kdelay * (q[7] - q[14]) 
    dq[15] = kdelay * (q[14] - q[15])                                         #delay2
    dq[16] = kdelay * (q[15] - q[16])                                         #delay3
    dq[17] = kdelay * (q[16] - q[17])                                         #delay4
    dq[18] = kdelay * (q[17] - q[18])                                         #delay5
    dq[19] = kdelay * (q[18] - q[19])                                         #delay6
end

"""
ODEs for latest thyrosim model. 
"""
function thyrosim(dq, q, p, t)
    kdelay = 5/8

    # Auxillary equations
    q4F = (p[24]+ p[25] * q[1] + p[26] * q[1]^2 + p[27] *q[1]^3) * q[4] #FT3p
    q1F = (p[7] + p[8] * q[1] + p[9] * q[1]^2 + p[10] * q[1]^3) * q[1]  #FT4p
    SR3 = (p[19] * q[19])                                        #Brain delay (dial 3)
    SR4 = (p[1] * q[19])                                         #Brain delay (dial 1)
    fCIRC = q[9]^p[51] / (q[9]^p[51] + p[49]^p[51])
    SRTSH = (p[30]+p[31]*fCIRC*sin(pi/12*(t+p[60])-p[33]))*(p[50]^p[52]/(p[50]^p[52] + q[9]^p[52]))
    fdegTSH = p[34] + p[35] / (p[36] + q[7])
    fLAG = p[41] + 2*q[8]^11 / (p[42]^11 + q[8]^11)
    f4 = p[37]*(1 + 5*(p[53]^p[54]) / (p[53]^p[54]+q[8]^p[54]))
    NL = p[13] / (p[14] + q[2])

    # ODEs
    dq[1]  = SR4 + p[3] * q[2] + p[4] * q[3] - (p[5] + p[6]) * q1F + p[11] * q[11] #T4dot (need to remove u1)
    dq[2]  = p[6] * q1F - (p[3] + p[12] + NL) * q[2]                                    #T4fast
    dq[3]  = p[5] * q1F -(p[4] + p[15] / (p[16] + q[3]) + p[17] /(p[18] + q[3])) *q[3]  #T4slow
    dq[4]  = SR3 + p[20] * q[5] + p[21] * q[6] - (p[22] + p[23]) * q4F + p[28] * q[13]  #T3pdot
    dq[5]  = p[23] * q4F + NL * q[2] - (p[20] + p[29]) * q[5]                         #T3fast
    dq[6]  = p[22] * q4F + p[15] * q[3] / (p[16] + q[3]) + p[17] * q[3] / (p[18] + q[3]) -(p[21])*q[6] #T3slow
    dq[7]  = SRTSH - fdegTSH * q[7]                                           #TSHp
    dq[8]  = f4 / p[38] * q[1] + p[37] / p[39] * q[4] - p[40] * q[8]          #T3B
    dq[9]  = fLAG * (q[8] - q[9])                                             #T3B LAG
    dq[10] = -p[43] * q[10]                                                   #T4PILLdot
    dq[11] =  p[43] * q[10] - (p[44] + p[11]) * q[11]                         #T4GUTdot
    dq[12] = -p[45] * q[12]                                                   #T3PILLdot
    dq[13] =  p[45] * q[12] - (p[46] + p[28]) * q[13]                         #T3GUTdot

    # Delay ODEs
    dq[14] = kdelay * (q[7] - q[14]) 
    dq[15] = kdelay * (q[14] - q[15])                                         #delay2
    dq[16] = kdelay * (q[15] - q[16])                                         #delay3
    dq[17] = kdelay * (q[16] - q[17])                                         #delay4
    dq[18] = kdelay * (q[17] - q[18])                                         #delay5
    dq[19] = kdelay * (q[18] - q[19])                                         #delay6
end

function output_equations(sol, p)
    return [777.0 * sol[1, :] / p[47], #T4
            651.0 * sol[4, :] / p[47], #T3
            5.6 * sol[7, :] / p[48]] #TSH
end

"""
Set initial conditions from data. Options to set other compartments to steady state,
optionally including the TSH lag compartments.
"""
function set_patient_ic!(ic, p, t4, t3, tsh;
        steady_state::Bool=false, set_tsh_lag::Bool=false)
    # Set IC for observed compartments. 
    ic[1] = (p[47] * t4) / 777.0
    ic[4] = (p[47] * t3) / 651.0
    ic[7] = (p[48] * tsh) / 5.6
    
    if steady_state
        q4F = (p[24]+ p[25] * ic[1] + p[26] * ic[1]^2 + p[27] *ic[1]^3) * ic[4] #FT3p
        q1F = (p[7] + p[8] * ic[1] + p[9] * ic[1]^2 + p[10] * ic[1]^3) * ic[1]  #FT4p
        
        B = p[6] * q1F - p[14] * (p[3] + p[12]) - p[13]
        A = -(p[3] + p[12])
        C = p[6] * p[14] * q1F
        ic[2] = (-B - sqrt(B^2 - 4.0 * A *C)) / (2.0 * A)
        
        B = p[5] * q1F - (p[4] + p[15] / p[16]) * p[18] - p[17]
        A = -(p[4] + p[15] / p[16])
        C = p[5] * p[18] * q1F
        ic[3] = (-B - sqrt(B^2 - 4.0 * A *C)) / (2.0 * A)
        
        ic[5] = (p[23] * q4F + (p[13] / (p[14] + ic[2])) * ic[2]) / (p[20] + p[29])
        ic[6] = (p[22] * q4F + p[15] * (ic[3] / (p[16] + ic[3]))
            + p[17] * (ic[3] / (p[18] + ic[3]))) / p[21]
    end
    
    if set_tsh_lag
        # Probably not 100% correct since they're supposed to be lagged, but probably better than the default.
        ic[14:19] .= ic[7]
    end
end
