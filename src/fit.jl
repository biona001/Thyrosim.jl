"""
Helper function that sets the parameter of `p` that are not in `fitting_index`
to its original value. 
"""
function reset_p!(
	p::Vector{Float64}, 
	fitting_index::Vector{Int},
	dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88]
	)

	ic, p0 = initialize(dial)

	for i in 1:length(p)
		if i in fitting_index 
			continue
		else
			p[i] = p0[i]
		end
	end

	return nothing
end

"""
In place version of initialize; avoids initializing ic and p in every loop. 
Used in fitting functions only.
"""
function initialize!(
    ic::Vector,
    p::Vector,
    dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88],
    scale_Vp::Bool = true,
    height=1.77,
    weight=70,
    sex=true #true = male, false = female
    )

    # TODO: need to calculate initial steady state
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
    p[51] = 6.91 #11.0     #hill exponent in f_circ -> fitted to Blakesley data
    p[52] = 7.66 #5.0      #hill exponent in SR_tsh -> fitted to Blakesley data
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
	compute_schneider_error(train)

Computes training error on the `train` dataset

# Error is defined as:
- When given an euthyroid T4 dose, if any TSH values not in [0.5, 4.5] in the last 
24h of simulation, then error + 1 (i.e. **patients receiving correct dose should 
have normal TSH**)
- When given the initial T4 dose, if the initial T4 dose is not equal to euthyroid 
T4 dose, and all TSH values is in [0.5, 4.5], then error + 1 (i.e. **patients not 
receiving correct dose should NOT have normal TSH**)

# Parameter definition:
- `p[55]:` Daily T4 oral dose
- `p[56]:` Daily T3 oral dose
"""
function compute_schneider_error(train_data)
    dial = [0.0; 0.88; 0.0; 0.88]
    scale_Vp = true
    tot_loss = 0.0
    
    # define function for adding dose
    function add_dose!(integrator)
        integrator.u[10] += integrator.p[55]
        integrator.u[12] += integrator.p[56]
    end
    cbk = PeriodicCallback(add_dose!, 24.0);

    #loop over all patients
    for i in 1:size(train_data, 1)
        height = train_data[i, Symbol("Ht.m")]
        weight = train_data[i, Symbol("Wt.kg")]
        sex    = Bool(train_data[i, Symbol("Sex")])
        ic, p  = initialize(dial, scale_Vp, height, weight, sex)
        ic[7]  = train_data[i, Symbol("TSH.preop")] #set initial TSH value
        tspan  = (0.0, 24.0train_data[i, Symbol("Days.to.euthyroid")]) #(0, total hours)
        
        # calculate error for euthyroid dose
        euthyroid_dose = train_data[i, Symbol("LT4.euthyroid.dose")] / 777.0
        p[55] = euthyroid_dose
        prob  = ODEProblem(thyrosim,ic,tspan,p,callback=cbk)
        sol   = solve(prob, save_idxs=7)
        tot_loss += compute_euthyroid_dose_error(sol)
        
        # when initial dose != euthyroid dose, calculate error
        initial_dose = train_data[i, Symbol("LT4.initial.dose")] / 777.0
        if initial_dose != euthyroid_dose
            p[55] = initial_dose
            prob  = ODEProblem(thyrosim,ic,tspan,p,callback=cbk)
            sol   = solve(prob, save_idxs=7)
            tot_loss += compute_initial_dose_error(sol)
        end
    end
    
    return tot_loss
end
