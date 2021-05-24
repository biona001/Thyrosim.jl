using DifferentialEquations
using Thyrosim
using DiffEqCallbacks
using DiffEqParamEstim
using Optim
using Statistics
using LinearAlgebra
BLAS.set_num_threads(1)

function objective(
    p_being_optimized::Vector, # first n elem is parameters, n+1:end is T4/T3 secrete rates
    fitting_index::Vector, 
    lowerbound::Vector,
    upperbound::Vector,
    blakesley_time::Matrix,
    blakesley_my400_data::Matrix, 
    blakesley_my450_data::Matrix, 
    blakesley_my600_data::Matrix,
    jonklaas_time::Vector, 
    jonklaas_patient_t4::Matrix, 
    jonklaas_patient_t3::Matrix, 
    jonklaas_patient_tsh::Matrix, 
    jonklaas_patient_param::Matrix, 
    jonklaas_patient_dose::Matrix,
    jonklaas_exclude_idx::Vector,
    jonklaas_secrete_rate_clusters::Vector,
    schneider_height::Vector, 
    schneider_weight::Vector, 
    schneider_sex::Vector, 
    schneider_tspan::Vector, 
    schneider_init_tsh::Vector, 
    schneider_euthy_dose::Vector, 
    schneider_init_dose::Vector,
    schneider_postTSH::Vector;      
    verbose::Bool = false, #set to true to display intermediate errors,
    blakesley_tsh_penalty::Float64 = 0.0,
    scale_plasma_ode=false,
    scale_slow_ode=false,
    scale_fast_ode=false,
    scale_allometric_exponent::Bool = false,
    scale_clearance::Bool = false
    )
    total_neg_logl = 0.0
    # quick return
    for i in eachindex(p_being_optimized)
        if !(lowerbound[i] ≤ p_being_optimized[i] ≤ upperbound[i])
            return Inf
        end
    end
    #
    # Blakesley male patient with BMI = p[65] (assuming height 1.77m)
    #
    bmi = p_being_optimized[findfirst(x -> x == 65, fitting_index)]
    w = bmi * 1.77^2 # BMI * h^2
    ic, p = initialize([1.0; 0.88; 1.0; 0.88], true, 1.77, w, true, 
        fitting_index=fitting_index, p_being_optimized=p_being_optimized, 
        scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
        scale_fast_ode=scale_fast_ode, scale_allometric_exponent=scale_allometric_exponent,
        scale_clearance=scale_clearance) 
    tspan = (0.0, 120.0)
    cbk   = ContinuousCallback(blakesley_condition, add_dose!); 
    p_400 = copy(p)
    p_450 = copy(p)
    p_600 = copy(p)
    find_patient_ic!(ic, p_400, 30)
    # solve different ODE problems for varying doses
    p_400[55] = 400.0 / 777.0
    p_450[55] = 450.0 / 777.0
    p_600[55] = 600.0 / 777.0
    prob_400 = ODEProblem(thyrosim,ic,tspan,p_400,callback=cbk)
    prob_450 = ODEProblem(thyrosim,ic,tspan,p_450,callback=cbk)
    prob_600 = ODEProblem(thyrosim,ic,tspan,p_600,callback=cbk)
    sol_400 = solve(prob_400, save_idxs=[1, 4, 7])
    sol_450 = solve(prob_450, save_idxs=[1, 4, 7])
    sol_600 = solve(prob_600, save_idxs=[1, 4, 7])
    T4_error = blakesley_t4_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[47], p[61]) + 
               blakesley_t4_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[47], p[61]) + 
               blakesley_t4_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[47], p[61])
    T3_error = blakesley_t3_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[47], p[62]) + 
               blakesley_t3_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[47], p[62]) + 
               blakesley_t3_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[47], p[62])
    TSH_error = blakesley_tsh_neg_logl(sol_400,blakesley_time,blakesley_my400_data,p[48],p[63],blakesley_tsh_penalty) + 
                blakesley_tsh_neg_logl(sol_450,blakesley_time,blakesley_my450_data,p[48],p[63],blakesley_tsh_penalty) + 
                blakesley_tsh_neg_logl(sol_600,blakesley_time,blakesley_my600_data,p[48],p[63],blakesley_tsh_penalty)
    blakesley_male_err = T4_error + T3_error + TSH_error
    verbose && println("blakesley male neg logl: T4 = $T4_error, T3 = $T3_error, TSH = $TSH_error")
    total_neg_logl += blakesley_male_err
    #
    # Blakesley female with BMI = p[66] (assuming height 1.63m)
    #
    bmi = p_being_optimized[findfirst(x -> x == 66, fitting_index)]
    w = bmi * 1.63^2 # BMI * h^2
    ic, p = initialize([1.0; 0.88; 1.0; 0.88], true, 1.63, w, false, 
        fitting_index=fitting_index, p_being_optimized=p_being_optimized,
        scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
        scale_fast_ode=scale_fast_ode, scale_allometric_exponent=scale_allometric_exponent,
        scale_clearance=scale_clearance)  
    tspan = (0.0, 120.0)
    cbk   = ContinuousCallback(blakesley_condition, add_dose!); 
    p_400 = copy(p)
    p_450 = copy(p)
    p_600 = copy(p)
    find_patient_ic!(ic, p_400, 30)
    # solve different ODE problems for varying doses
    p_400[55] = 400.0 / 777.0
    p_450[55] = 450.0 / 777.0
    p_600[55] = 600.0 / 777.0
    prob_400 = ODEProblem(thyrosim,ic,tspan,p_400,callback=cbk)
    prob_450 = ODEProblem(thyrosim,ic,tspan,p_450,callback=cbk)
    prob_600 = ODEProblem(thyrosim,ic,tspan,p_600,callback=cbk)
    sol_400 = solve(prob_400, save_idxs=[1, 4, 7])
    sol_450 = solve(prob_450, save_idxs=[1, 4, 7])
    sol_600 = solve(prob_600, save_idxs=[1, 4, 7])
    T4_error = blakesley_t4_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[47], p[61]) + 
               blakesley_t4_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[47], p[61]) + 
               blakesley_t4_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[47], p[61])
    T3_error = blakesley_t3_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[47], p[62]) + 
               blakesley_t3_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[47], p[62]) + 
               blakesley_t3_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[47], p[62])
    TSH_error = blakesley_tsh_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[48], p[63],blakesley_tsh_penalty) + 
                blakesley_tsh_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[48], p[63],blakesley_tsh_penalty) + 
                blakesley_tsh_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[48], p[63],blakesley_tsh_penalty)
    blakesley_err = T4_error + T3_error + TSH_error
    verbose && println("blakesley female neg logl: T4 = $T4_error, T3 = $T3_error, TSH = $TSH_error")
    total_neg_logl += blakesley_err
    #
    # Old Jonklaas data
    #
#     tspan = (0.0, 8.0)
#     cbk   = ContinuousCallback(jonklaas_condition, add_dose!); 
#     jonklaas_err = 0.0
#     for i in 1:size(jonklaas_patient_dose, 1)
#         # initialize ODE problem for patient 1
#         height, weight, sex = jonklaas_patient_param[i, :]
#         T4init, T3init, TSHinit = jonklaas_patient_t4[i, 1], jonklaas_patient_t3[i, 1], jonklaas_patient_tsh[i, 1]
#         ic, p = initialize([0.0; 0.88; 0.0; 0.88], true, height, weight, Bool(sex)) 
#         set_patient_ic!(ic, p, T4init, T3init, TSHinit, steady_state=true, set_tsh_lag=true)
#         # set parameters being fitted
#         p[fitting_index] .= @view(p_being_optimized[1:length(fitting_index)])
#         # solve different ODE problems for varying doses
#         p[56] = jonklaas_patient_dose[i] / 651.0
#         prob = ODEProblem(thyrosim,ic,tspan,p,callback=cbk)
#         sol  = solve(prob, save_idxs=4)
#         jonklaas_err += jonklaas_t3_neg_logl(sol, jonklaas_time, 
#             jonklaas_patient_t3[i, :], p[47], p[62])
#     end
#     verbose && println("jonklaas's negative loglikelihood = $jonklaas_err")
#     total_neg_logl += jonklaas_err
    #
    # Updated Jonklaas data
    #
    jonklaas_err = T4_error = T3_error = TSH_error = 0.0
    dial = [1.0; 0.88; 1.0; 0.88]
    nparams = length(fitting_index)
    nsamples = size(jonklaas_patient_param, 1)
    cbk = PeriodicCallback(add_dose!, 24.0) 
    weight_w1 = jonklaas_patient_param[:, 1] # week 1 weight
    height = jonklaas_patient_param[:, 3]
    sex = convert(BitVector, jonklaas_patient_param[:, 4])
    for i in 1:nsamples
        i in jonklaas_exclude_idx && continue
        # run model to steady state before actual simulation
        dial[1] = dial[3] = 1.0
        sol = simulate(height[i], weight_w1[i], sex[i], days=50, dial=dial, warmup=false, 
            fitting_index=fitting_index, parameters=p_being_optimized[1:length(fitting_index)])
        _, p = initialize(dial, true, height[i], weight_w1[i], sex[i], fitting_index=fitting_index,
            p_being_optimized=p_being_optimized, scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
            scale_fast_ode=scale_fast_ode, scale_allometric_exponent=scale_allometric_exponent,
            scale_clearance=scale_clearance)
        T4_error += jonklaas_T4_neg_logl(sol, jonklaas_patient_t4[i, 2], p[47], p[64])
        T3_error += jonklaas_T3_neg_logl(sol, jonklaas_patient_t3[i, 2], p[47], p[62])
        TSH_error += jonklaas_TSH_neg_logl(sol, jonklaas_patient_tsh[i, 2], p[47], p[63])
        # run first 8 week simulations, interpolate weight weekly
        weight_diff = (jonklaas_patient_param[i, 2] - jonklaas_patient_param[i, 1]) / 16.0
        dial[1] = dial[3] = 0.0
        for week in 1:8
            # reset parameters using new weight
            ic, p = initialize(dial, true, height[i], weight_w1[i] + week*weight_diff, sex[i],
                fitting_index=fitting_index, p_being_optimized=p_being_optimized,
                scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
                scale_fast_ode=scale_fast_ode, scale_allometric_exponent=scale_allometric_exponent,
                scale_clearance=scale_clearance)
            p[55] = jonklaas_patient_dose[i, 1] / 777.0
            p[fitting_index] .= @view(p_being_optimized[1:length(fitting_index)])
            # use last week's end value
            ic .= sol[end]
            ic[10] += p[55] # manually add dose for first day of the week
            prob = ODEProblem(thyrosim,ic,(0.0, 168.0),p,callback=cbk)
            sol  = solve(prob)
        end
        T4_error += jonklaas_T4_neg_logl(sol, jonklaas_patient_t4[i, 3], p[47], p[64])
        T3_error += jonklaas_T3_neg_logl(sol, jonklaas_patient_t3[i, 3], p[47], p[62])
        TSH_error += jonklaas_TSH_neg_logl(sol, jonklaas_patient_tsh[i, 3], p[47], p[63])
        # run next 8 week, interpolate weight weekly
        for week in 9:16
            # reset parameters using new weight
            ic, p = initialize(dial, true, height[i], weight_w1[i] + week*weight_diff, sex[i],
                fitting_index=fitting_index, p_being_optimized=p_being_optimized,
                scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
                scale_fast_ode=scale_fast_ode, scale_allometric_exponent=scale_allometric_exponent,
                scale_clearance=scale_clearance)
            p[55] = jonklaas_patient_dose[i, 2] / 777.0
            # use last week's end value
            ic .= sol[end]
            ic[10] += p[55] # manually add dose for first day of the week
            prob = ODEProblem(thyrosim,ic,(0.0, 168.0),p,callback=cbk)
            sol  = solve(prob)
        end
        T4_error += jonklaas_T4_neg_logl(sol, jonklaas_patient_t4[i, 3], p[47], p[64])
        T3_error += jonklaas_T3_neg_logl(sol, jonklaas_patient_t3[i, 3], p[47], p[62])
        TSH_error += jonklaas_TSH_neg_logl(sol, jonklaas_patient_tsh[i, 3], p[47], p[63])
    end
    verbose && println("jonklaas neg logl: T4 = $T4_error, T3 = $T3_error, TSH = $TSH_error")
    jonklaas_err = T4_error + T3_error + TSH_error
    total_neg_logl += jonklaas_err
    # 
    # Schneider
    #
#     num_params = length(fitting_index)
#     num_sample = length(schneider_height)
#     num_params == length(fitting_index) || error("check parameter length bro")
#     schneider_logl = zeros(Threads.nthreads())
#     Threads.@threads for i in 1:num_sample
#         # get current patient's parameters
#         sex = schneider_sex[i]
#         h = schneider_height[i]
#         w = schneider_weight[i]
#         euthyroid_dose = schneider_euthy_dose[i]
#         initial_dose = schneider_init_dose[i]
#         initial_tsh = schneider_init_tsh[i]
#         id = Threads.threadid()
#         # calculate observations neg logl (units in mcg/kg)
#         μ = schneider_postTSH[i] # post operative TSH values after 6 weeks treatment
#         σ = p[62]
#         y = schneider_end_tsh(p_being_optimized, fitting_index, h, w, sex, 
#             initial_tsh, initial_dose)
#         schneider_logl[id] += log(2π) / 2 + log(σ) + (y - μ)^2 / 2 / σ^2
#     end
#     verbose && println("schneider's negative loglikelihood = $(sum(schneider_logl))")
#     total_neg_logl += sum(schneider_logl)
    #
    # Return final error
    #
    return total_neg_logl
end

# gives 400/450/600 mcg of oral T4 at hour 24
function blakesley_condition(u, t, integrator)
    return t - 24.0
end
# gives T3 dose at hour 0
function jonklaas_condition(u, t, integrator)
    return t - 0.01 #cannot make this exactly 0
end
# gives T4 dose at hour 24
function new_jonklaas_condition(u, t, integrator)
    return t - 24.0
end
# define function for adding dose
function add_dose!(integrator)
    integrator.u[10] += integrator.p[55]
    integrator.u[12] += integrator.p[56]
end
function blakesley_t4_neg_logl(sol, time, data, Vp, σ) # sol includes T4/T3/TSH only
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        n = length(time)
        for i in 1:n
            T4_predicted = sol(time[i])[1] * 777.0 / Vp
            tot_loss += (T4_predicted - data[i, 1])^2
        end
        tot_loss /= 2σ^2
        tot_loss += n * log(2π) / 2 + n * log(σ)
    end
    return tot_loss
end
function blakesley_t3_neg_logl(sol, time, data, Vp, σ) # sol includes T4/T3/TSH only
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        n = length(time)
        for i in 1:n
            T3_predicted = sol(time[i])[2] * 651.0 / Vp
            tot_loss += (T3_predicted - data[i, 2])^2
        end
        tot_loss /= 2σ^2
        tot_loss += n * log(2π) / 2 + n * log(σ)
    end
    return tot_loss
end
function blakesley_tsh_neg_logl(sol, time, data, Vtsh, σ, λ = 0.0) # sol includes T4/T3/TSH only
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        n = length(time)
        for i in 1:length(time)
            predicted_tsh = sol(time[i])[3] * 5.6 / Vtsh
            tot_loss += (predicted_tsh - data[i, 3])^2
        end
        tot_loss /= 2σ^2
        tot_loss += n * log(2π) / 2 + n * log(σ)
        # add penalty for 1st and 2nd high/low peak in TSH data (no penalty if λ=0)           
        for i in [9, 13, 24, 28]
            predicted_tsh = sol(time[i])[3] * 5.6 / Vtsh
            tot_loss += λ * (predicted_tsh - data[i, 3])^2
        end
    end
    return tot_loss
end
"""
Converts TT4 (μmol) to FT4 (μmol) using Thyrosim's internal 4th order polynomial
"""
function FT4(TT4::Float64)
    return (0.000289 + (0.000214 * TT4) + (0.000128 * TT4^2) - (8.83e-6 * TT4^3)) * TT4
end

"""
Converts TT4 (μmol, Thyrosim output q1) to FT4 (ng/dL, jonklaas's units) 
"""
function TT4_to_FT4(TT4::Float64, Vp::Float64)
    return FT4(TT4) * 777 / Vp * 1000 * 0.45 / 10
end
# calculate new jonklaas error for T4 at end of sol
function jonklaas_T4_neg_logl(sol, data, Vp, σ) # sol includes all comparments
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        predicted_FT4 = TT4_to_FT4(sol[end][1], Vp)
        tot_loss += (predicted_FT4 - data)^2 / 2σ^2 + log(2π) / 2 + log(σ)
    end
    return tot_loss
end
# calculate new jonklaas error for T3 at end of sol
function jonklaas_T3_neg_logl(sol, data, Vp, σ) # sol includes all comparments
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        predicted = sol[end][4] * 651.0 / Vp
        tot_loss += (predicted - data)^2 / 2σ^2 + log(2π) / 2 + log(σ)
    end
    return tot_loss
end
# calculate new jonklaas error for TSH at end of sol
function jonklaas_TSH_neg_logl(sol, data, Vtsh, σ) # sol includes all comparments
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        predicted = sol[end][7] * 5.6 / Vtsh
        tot_loss += (predicted - data)^2 / 2σ^2 + log(2π) / 2 + log(σ)
    end
    return tot_loss
end
function schneider_end_tsh( 
    current_iter::Vector,
    fitting_index::Vector,
    height::Float64, 
    weight::Float64, 
    sex::Bool,
    initial_tsh::Float64,
    initial_dose::Float64
    )

    #initialize simulation parameters
    scale_Vp = true
    dial  = [0.0; 0.88; 0.0; 0.88]
    ic, p = initialize(dial, scale_Vp, height, weight, sex)
    ic[7] = initial_tsh
    cbk = PeriodicCallback(add_dose!, 24.0)# function to add dose

    # update parameter for fitting 
    p[fitting_index] .= @view(current_iter[1:length(fitting_index)])

    # run ODE simulation
    p[55] = initial_dose / 777.0
    prob  = ODEProblem(thyrosim,ic,(0.0, 1008),p,callback=cbk) # simulate for 6 weeks
    sol   = solve(prob, Tsit5(), save_idxs=[1, 7])

    # return observed value
    return sol.u[end][2] * 5.6 / p[48]
end
function schneider_TSH_err(initTSH, sol, Vtsh)
    tsh = sol.u[end][2] * 5.6 / Vtsh
    return true
end
function TSH_within_interval(sol, Vtsh)
    tot_loss = 0
    tsh = sol.u[end][2] * 5.6 / Vtsh
    if 0.5 ≤ tsh ≤ 4.5
        return true
    end
    return false
end
# Calculate mean/std of eythyroid dose (mcg/kg) for 
# male and female patients in different BMI categories 
function compute_patient_categories(
    sex::AbstractVector, 
    bmi::AbstractVector,
    euthyroid_dose::AbstractVector
    )
    categories = Dict{Symbol, Tuple{Float64, Float64}}()
    
    # get index for different cateories
    male_normal_idx = intersect(findall(iszero, sex), findall(x -> x < 24.9, bmi))
    male_overweight_idx = intersect(findall(iszero, sex), findall(x -> 24.9 <= x < 29.9, bmi))
    male_obese_idx = intersect(findall(iszero, sex), findall(x -> 29.9 <= x, bmi))
    female_normal_idx = intersect(findall(isone, sex), findall(x -> x < 24.9, bmi))
    female_overweight_idx = intersect(findall(isone, sex), findall(x -> 24.9 <= x < 29.9, bmi))
    female_obese_idx = intersect(findall(isone, sex), findall(x -> 29.9 <= x, bmi))
    
    # compute mean and var of euthyroid dose. If empty, set both as 0
    categories[:male_normal] = isempty(male_normal_idx) ? (0, 0) : 
        (mean(euthyroid_dose[male_normal_idx]), std(euthyroid_dose[male_normal_idx]))
    categories[:male_overweight] = isempty(male_overweight_idx) ? (0, 0) : 
        (mean(euthyroid_dose[male_overweight_idx]), std(euthyroid_dose[male_overweight_idx]))
    categories[:male_obese] = isempty(male_obese_idx) ? (0, 0) : 
        (mean(euthyroid_dose[male_obese_idx]), std(euthyroid_dose[male_obese_idx]))
    categories[:female_normal] = isempty(female_normal_idx) ? (0, 0) : 
        (mean(euthyroid_dose[female_normal_idx]), std(euthyroid_dose[female_normal_idx]))
    categories[:female_overweight] = isempty(female_overweight_idx) ? (0, 0) : 
        (mean(euthyroid_dose[female_overweight_idx]), std(euthyroid_dose[female_overweight_idx]))
    categories[:female_obese] = isempty(female_obese_idx) ? (0, 0) : 
        (mean(euthyroid_dose[female_obese_idx]), std(euthyroid_dose[female_obese_idx]))

    return categories
end
# male = 0, female = 1. Height in meters and weight in kg
function get_category(categories::Dict, sex, height, weight)
    bmi = weight / height^2
    if sex == 0
        if bmi < 24.9 
            return categories[:male_normal]
        elseif 24.9 <= bmi < 29.9
            return categories[:male_overweight]
        elseif 29.9 <= bmi
            return categories[:male_obese]
        else
            error("male patient without assigned category!")
        end
    elseif sex == 1
        if bmi < 24.9 
            return categories[:female_normal]
        elseif 24.9 <= bmi < 29.9
            return categories[:female_overweight]
        elseif 29.9 <= bmi
            return categories[:female_obese]
        else
            error("female patient without assigned category!")
        end
    else
        error("undefined sex!")
    end
end
function update_logl_by_category!(logl_by_category::Vector, logl, sex, height, weight)
    bmi = weight / height^2
    if sex == 0 # male
        if bmi < 24.9 
            logl_by_category[8] += logl
        elseif 24.9 <= bmi < 29.9
            logl_by_category[16] += logl
        elseif 29.9 <= bmi
            logl_by_category[24] += logl
        else
            error("male patient without assigned category!")
        end
    elseif sex == 1
        if bmi < 24.9 
            logl_by_category[32] += logl
        elseif 24.9 <= bmi < 29.9
            logl_by_category[40] += logl
        elseif 29.9 <= bmi
            logl_by_category[48] += logl
        else
            error("female patient without assigned category!")
        end
    else
        error("undefined sex!")
    end
end

function fit_all()
    fitting_index =
        [1; 13;                  # S4, VtshMax
        30; 31; 37               # A0, B0, k3
        49; 50; 51; 52; 53; 54;  # hill function parameters
        65; 66; 72; 73]          # reference male/female BMI, fat-free and fat constant
    initial_guess = [0.0019892210815454564, 0.012318557740933649, 78.03368752668696, 63.079747932889816,
        0.06578735870878696, 3.3739342983833187, 4.39393376334155, 7.183642942358456, 8.91034232003827,
        6.863194346722813, 18.848701766376884, 23.929032682987728, 22.5, 0.5, 0.5]
    lowerbound = zeros(length(initial_guess))
    upperbound = initial_guess .* 10.0
    lowerbound[findall(x -> x == 65, fitting_index)] .= 20.0
    lowerbound[findall(x -> x == 66, fitting_index)] .= 20.0
    upperbound[findall(x -> x == 65, fitting_index)] .= 25.0
    upperbound[findall(x -> x == 66, fitting_index)] .= 25.0
    upperbound[findall(x -> x == 54, fitting_index)] .= 20.0
    upperbound[findall(x -> x == 72, fitting_index)] .= 1.0
    upperbound[findall(x -> x == 73, fitting_index)] .= 1.0

    # whether to scale plasma compartments by the Vp ratio
    scale_plasma_ode = false
    scale_slow_ode = false
    scale_fast_ode = false
    scale_allometric_exponent = false
    scale_clearance = true
    
    # blakesley setup 
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    blakesley_tsh_penalty = 100.0 # penalize the peak TSH values 

    # jonklaas setup
#     jonklaas_exclude_idx = [4, 7, 11, 18, 22, 23, 26, 36, 44, 45, 46] # TSH > 2 at final time point
#     jonklaas_exclude_idx = [7, 11, 18, 23, 26, 36, 44, 45, 46] # TSH > 3 at final time point
    jonklaas_exclude_idx = [11, 18, 44] # TSH > 4 at final time point
    jonklaas_secrete_rate_clusters = [4,2,2,1,3,1,3,1,2,1,2,4,1,2,1,3,1,2,3,1,3,4,4,1,1,3,4,1,1,1,2,1,1]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data_new()
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    
    # schneider setup
    alldata, train_data = schneider_data();
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")]) # 1 is female
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])
    postTSH = convert(Vector{Float64}, train_data[!, Symbol("6 week TSH")])
    
    return optimize(p -> objective(p, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, patient_t4, 
        patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        jonklaas_exclude_idx, jonklaas_secrete_rate_clusters, height, weight, sex, tspan, 
        init_tsh, euthy_dose, init_dose, postTSH, verbose=false, 
        blakesley_tsh_penalty=blakesley_tsh_penalty, scale_plasma_ode=scale_plasma_ode,
        scale_slow_ode=scale_slow_ode, scale_fast_ode=scale_fast_ode, 
        scale_allometric_exponent=scale_allometric_exponent, scale_clearance=scale_clearance),
            initial_guess, NelderMead(),
            Optim.Options(time_limit = 72000.0, iterations = 10000, g_tol=1e-5,
            show_trace = true, allow_f_increases=true))
end

function prefit_error()
    fitting_index =
        [1; 13;                  # S4, VtshMax
        30; 31; 37               # A0, B0, k3
        49; 50; 51; 52; 53; 54;  # hill function parameters
        65; 66; 72; 73]          # reference male/female BMI, fat-free and fat constant
    initial_guess = [0.0019892210815454564, 0.012318557740933649, 78.03368752668696, 63.079747932889816,
        0.06578735870878696, 3.3739342983833187, 4.39393376334155, 7.183642942358456, 8.91034232003827,
        6.863194346722813, 18.848701766376884, 23.929032682987728, 22.5, 0.5, 0.5]
    lowerbound = zeros(length(initial_guess))
    upperbound = initial_guess .* 10.0

    # whether to scale plasma compartments by the Vp ratio
    scale_plasma_ode = false
    scale_slow_ode = false
    scale_fast_ode = false
    scale_allometric_exponent = false
    scale_clearance = true

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    blakesley_tsh_penalty = 100.0
    # jonklaas setup
    # jonklaas_exclude_idx = [4, 7, 11, 18, 22, 23, 26, 36, 44, 45, 46] # TSH > 2 at final time point
    # jonklaas_exclude_idx = [7, 11, 18, 23, 26, 36, 44, 45, 46] # TSH > 3 at final time point
    jonklaas_exclude_idx = [11, 18, 44] # TSH > 4 at final time point
    jonklaas_secrete_rate_clusters = [4,2,2,1,3,1,3,1,2,1,2,4,1,2,1,3,1,2,3,1,3,4,4,1,1,3,4,1,1,1,2,1,1]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data_new()
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
#     jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    alldata, train_data = schneider_data()
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])
    postTSH = convert(Vector{Float64}, train_data[!, Symbol("6 week TSH")])

    return objective(initial_guess, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, patient_t4, 
        patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        jonklaas_exclude_idx, jonklaas_secrete_rate_clusters, height, weight, sex, 
        tspan, init_tsh, euthy_dose, init_dose, postTSH, verbose=true, 
        blakesley_tsh_penalty=blakesley_tsh_penalty, scale_plasma_ode=scale_plasma_ode, 
        scale_slow_ode=scale_slow_ode, scale_fast_ode=scale_fast_ode, 
        scale_allometric_exponent = scale_allometric_exponent,scale_clearance=scale_clearance)
end

function postfit_error(minimizer)
    fitting_index =
        [1; 13;                  # S4, VtshMax
        30; 31; 37               # A0, B0, k3
        49; 50; 51; 52; 53; 54;  # hill function parameters
        65; 66; 72; 73]          # reference male/female BMI, fat-free and fat constant
    initial_guess = [0.0019892210815454564, 0.012318557740933649, 78.03368752668696, 63.079747932889816,
        0.06578735870878696, 3.3739342983833187, 4.39393376334155, 7.183642942358456, 8.91034232003827,
        6.863194346722813, 18.848701766376884, 23.929032682987728, 22.5, 0.5, 0.5]
    lowerbound = zeros(length(minimizer))
    upperbound = Inf .* ones(length(minimizer))

    # whether to scale plasma compartments by the Vp ratio
    scale_plasma_ode = false
    scale_slow_ode = false
    scale_fast_ode = false
    scale_allometric_exponent = false
    scale_clearance = true

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    blakesley_tsh_penalty = 100.0
    # jonklaas setup
    # jonklaas_exclude_idx = [4, 7, 11, 18, 22, 23, 26, 36, 44, 45, 46] # TSH > 2 at final time point
    # jonklaas_exclude_idx = [7, 11, 18, 23, 26, 36, 44, 45, 46] # TSH > 3 at final time point
    jonklaas_exclude_idx = [11, 18, 44] # TSH > 4 at final time point
    jonklaas_secrete_rate_clusters = [4,2,2,1,3,1,3,1,2,1,2,4,1,2,1,3,1,2,3,1,3,4,4,1,1,3,4,1,1,1,2,1,1]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data_new()
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
#     jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    alldata, train_data = schneider_data()
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])
    postTSH = convert(Vector{Float64}, train_data[!, Symbol("6 week TSH")])

    return objective(minimizer, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, patient_t4, 
        patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        jonklaas_exclude_idx, jonklaas_secrete_rate_clusters, height, weight, sex,
        tspan, init_tsh, euthy_dose, init_dose, postTSH, verbose=true,
        blakesley_tsh_penalty=blakesley_tsh_penalty, scale_plasma_ode=scale_plasma_ode, 
        scale_slow_ode=scale_slow_ode, scale_fast_ode=scale_fast_ode, 
        scale_allometric_exponent=scale_allometric_exponent,scale_clearance=scale_clearance)
end

println("Threads = ", Threads.nthreads())
flush(stdout)
println("prefit error: ")
prefit = prefit_error()
println("total prefit error = $prefit \n")
flush(stdout)

result = fit_all()

println("postfit error:")
post = postfit_error(result.minimizer)
println("total postfit error = $post \n")

println("result:")
println(result)
println(result.minimizer)
