# Initial conditions
function initialize(
    dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88],
    scale_Vp::Bool = true,
    height=1.77,
    weight=70,
    sex=true #true = male, false = female
    )

    # Initial conditions
    ic = [0.322114215761171;
          0.201296960359917;
          0.638967411907560;
          0.00663104034826483;
          0.0112595761822961;
          0.0652960640300348;
          1.78829584764370;
          7.05727560072869;
          7.05714474742141;
          0;
          0;
          0;
          0;
          3.34289716182018;
          3.69277248068433;
          3.87942133769244;
          3.90061903207543;
          3.77875734283571;
          3.55364471589659]

    # Parameter values
    p = [0.00174155 * dial[1];     #S4
         8;              #tau
         0.868;          #k12
         0.108;          #k13
         584;            #k31free
         1503;           #k21free
         0.000289;       #A
         0.000214;       #B
         0.000128;       #C
         -8.83*10^-6;    #D
         0.88;           #k4absorb; originally 0.881
         0.0189;         #k02
         0.00998996;     #VmaxD1fast
         2.85;           #KmD1fast
         6.63*10^-4;     #VmaxD1slow
         95;             #KmD1slow
         0.00074619;     #VmaxD2slow
         0.075;          #KmD2slow
         3.3572*10^-4 * dial[3];   #S3
         5.37;           #k45
         0.0689;         #k46
         127;            #k64free
         2043;           #k54free
         0.00395;        #a
         0.00185;        #b
         0.00061;        #c
         -0.000505;      #d
         0.88;           #k3absorb; originally 0.882
         0.207;          #k05
         1166;           #Bzero
         581;            #Azero
         2.37;           #Amax
         -3.71;          #phi
         0.53;           #kdegTSH-HYPO
         0.037;          #VmaxTSH
         23;             #K50TSH
         0.118;          #k3
         0.29;           #T4P-EU
         0.006;          #T3P-EU
         0.037;          #KdegT3B
         0.0034;         #KLAG-HYPO
         5;              #KLAG
         1.3;            #k4dissolve
         0.12 * dial[2]; #k4excrete; originally 0.119 (change with dial 2)
         1.78;           #k3dissolve
         0.12 * dial[4]; #k3excrete; originally 0.118 (change with dial 4)
         3.2;            #Vp
         5.2]            #VTSH

    if scale_Vp
        Vp, Vtsh = plasma_volume(height, weight, sex)
        p[47] = Vp
        p[48] = Vtsh
    end

    return ic, p
end

"""
    plasma_volume(height, weight, sex)

# Parameters used to get reference plasma volume (Vp) values:

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
    male_ref_vp   = 2.92;
    female_ref_vp = 2.48;

    # calculate Ideal Weight fitted to Feldschush's data
    if sex == 1
        iw = 176.3 - 220.6 * h + 93.5 * h^2
    elseif sex == 0
        iw = 145.8 - 182.7 * h + 79.55 * h^2
    end

    # hill function fitted to feldchush data
    a, n, K = 1.42634 * 10^4, 0.745964264, 100.0
    Δiw = (w - iw) / iw * 100  #deviation from ideal weight, in percentage
    Vb_per_kg = a * (100.0 + Δiw)^(n - 1) / ((100 + Δiw)^n + K^n)
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

# original thyrosim ODEs: https://bitbucket.org/DistefanoLab/thyrosim/src/master/resource/matlab/thyrosim_core.m
function original_thyrosim(dq, q, p, t)
    kdelay = 5/8

    # Auxillary equations
    q4F = (p[24]+ p[25] * q[1] + p[26] * q[1]^2 + p[27] *q[1]^3) * q[4]; #FT3p
    q1F = (p[7] + p[8] * q[1] + p[9] * q[1]^2 + p[10] * q[1]^3) * q[1];   #FT4p
    SR3 = (p[19] * q[19]);                                        #Brain delay (dial 3)
    SR4 = (p[1] * q[19]);                                         #Brain delay (dial 1)
    
    fCIRC = 1 + (p[32] / (p[31] * exp(-q[9])) - 1) * (1 / (1 + exp(10*q[9] - 55)));
    SRTSH = (p[30] + p[31] * fCIRC * sin(pi/12 * t - p[33])) * exp(-q[9]);
    fdegTSH = p[34] + p[35] / (p[36] + q[7]);
    fLAG = p[41] + 2*q[8]^11 / (p[42]^11 + q[8]^11);
    f4 = p[37] + 5 * p[37] / (1 + exp(2 * q[8] - 7));
    NL = p[13] / (p[14] + q[2]);

    # ODEs
    dq[1]  = SR4 + p[3] * q[2] + p[4] * q[3] - (p[5] + p[6]) * q1F + p[11] * q[11]; #T4dot (need to remove u1)
    dq[2]  = p[6] * q1F - (p[3] + p[12] + NL) * q[2];                                    #T4fast
    dq[3]  = p[5] * q1F -(p[4] + p[15] / (p[16] + q[3]) + p[17] /(p[18] + q[3])) *q[3];  #T4slow
    dq[4]  = SR3 + p[20] * q[5] + p[21] * q[6] - (p[22] + p[23]) * q4F + p[28] * q[13];  #T3pdot
    dq[5]  = p[23] * q4F + NL * q[2] - (p[20] + p[29]) * q[5];                         #T3fast
    dq[6]  = p[22] * q4F + p[15] * q[3] / (p[16] + q[3]) + p[17] * q[3] / (p[18] + q[3]) -(p[21])*q[6]; #T3slow
    dq[7]  = SRTSH - fdegTSH * q[7];                                           #TSHp
    dq[8]  = f4 / p[38] * q[1] + p[37] / p[39] * q[4] - p[40] * q[8];          #T3B
    dq[9]  = fLAG * (q[8] - q[9]);                                             #T3B LAG
    dq[10] = -p[43] * q[10];                                                   #T4PILLdot
    dq[11] =  p[43] * q[10] - (p[44] + p[11]) * q[11];                         #T4GUTdot
    dq[12] = -p[45] * q[12];                                                   #T3PILLdot
    dq[13] =  p[45] * q[12] - (p[46] + p[28]) * q[13];                         #T3GUTdot

    # Delay ODEs
    dq[14] = -kdelay * q[14] + q[7];                                           #delay1 CHECK, might be wrong 
    # dq[14] = kdelay * (q[7] - q[14]); 
    dq[15] = kdelay * (q[14] - q[15]);                                         #delay2
    dq[16] = kdelay * (q[15] - q[16]);                                         #delay3
    dq[17] = kdelay * (q[16] - q[17]);                                         #delay4
    dq[18] = kdelay * (q[17] - q[18]);                                         #delay5
    dq[19] = kdelay * (q[18] - q[19]);                                         #delay6
end

function output_equations(sol, p)
    return [777.0 * sol[1, :] / p[47], #T4
            651.0 * sol[4, :] / p[47], #T3
            5.6 * sol[7, :] / p[48]] #TSH
end
