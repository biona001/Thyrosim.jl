using DifferentialEquations
using Thyrosim
using DiffEqCallbacks
using DiffEqParamEstim
using SharedArrays
using Optim
using Statistics
function objective(
    p_being_optimized::Vector, 
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
    schneider_height::Vector, 
    schneider_weight::Vector, 
    schneider_sex::Vector, 
    schneider_tspan::Vector, 
    schneider_init_tsh::Vector, 
    schneider_euthy_dose::Vector, 
    schneider_init_dose::Vector;      
    verbose::Bool = false #set to true to display intermediate errors
    )
    
    total_neg_logl = 0.0
    # quick return
    for (i, val) in enumerate(p_being_optimized)
        if val < lowerbound[i] || val > upperbound[i]
            return Inf
        end
    end
    #
    # Blakesley
    #
    ic, p = initialize([1.0; 0.88; 1.0; 0.88]) 
    p[fitting_index] .= p_being_optimized
    tspan = (0.0, 120.0)
    cbk   = ContinuousCallback(blakesley_condition, add_dose!); 
    p_400 = p_450 = p_600 = copy(p)
    find_patient_ic!(ic, p_400, 30)
#     T4_init_400, T3_init_400, TSH_init_400 = blakesley_my400_data[1, :]
#     T4_init_450, T3_init_450, TSH_init_450 = blakesley_my450_data[1, :]
#     T4_init_600, T3_init_600, TSH_init_600 = blakesley_my600_data[1, :]
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
    TSH_error = blakesley_tsh_neg_logl(sol_400, blakesley_time, blakesley_my400_data, p[48], p[63]) + 
                blakesley_tsh_neg_logl(sol_450, blakesley_time, blakesley_my450_data, p[48], p[63]) + 
                blakesley_tsh_neg_logl(sol_600, blakesley_time, blakesley_my600_data, p[48], p[63])
    blakesley_err = T4_error + T3_error + TSH_error
    verbose && println("blakesley negative loglikelihood: T4 = $T4_error, T3 = $T3_error, TSH = $TSH_error")
    total_neg_logl += blakesley_err
    #
    # Jonklaas
    #
    tspan = (0.0, 8.0)
    cbk   = ContinuousCallback(jonklaas_condition, add_dose!); 
    jonklaas_err = 0.0
    for i in 1:size(jonklaas_patient_dose, 1)
        # initialize ODE problem for patient 1
        height, weight, sex = jonklaas_patient_param[i, :]
        T4init, T3init, TSHinit = jonklaas_patient_t4[i, 1], jonklaas_patient_t3[i, 1], jonklaas_patient_tsh[i, 1]
        ic, p = initialize([0.0; 0.88; 0.0; 0.88], true, height, weight, Bool(sex)) 
        set_patient_ic!(ic, p, T4init, T3init, TSHinit, steady_state=true, set_tsh_lag=true)
        # set parameters being fitted
        p[fitting_index] .= p_being_optimized
        # solve different ODE problems for varying doses
        p[56] = jonklaas_patient_dose[i] / 651.0
        prob = ODEProblem(thyrosim,ic,tspan,p,callback=cbk)
        sol  = solve(prob, save_idxs=4)
        jonklaas_err += jonklaas_t3_neg_logl(sol, jonklaas_time, 
            jonklaas_patient_t3[i, :], p[47], p[62])
    end
    verbose && println("jonklaas's negative loglikelihood = $jonklaas_err")
    total_neg_logl += jonklaas_err
    #
    # Schneider
    #
    num_params = length(p_being_optimized)
    num_sample = length(schneider_height)
    num_params == length(fitting_index) || error("check parameter length bro")
    schneider_logl = 0.0
    Threads.@threads for i in 1:num_sample
        # get current patient's parameters
        sex = schneider_sex[i]
        h = schneider_height[i]
        w = schneider_weight[i]
        euthyroid_dose = schneider_euthy_dose[i]
        initial_dose = schneider_init_dose[i]
        initial_tsh = schneider_init_tsh[i]
        
        # only calculate for patients with euthyroid dose ≥ initial dose
        if euthyroid_dose < initial_dose
            continue
        elseif euthyroid_dose == initial_dose
            μ = initial_tsh # assume patients maintain TSH.preop
        else # initial dose too low
            μ = 5.0 # assume patient had TSH value 5.0 
        end
        #calculate each observations neg likelihood (units in mcg/kg)
        σ = p[64]
        y = schneider_end_tsh(p_being_optimized, fitting_index, 
            h, w, sex, schneider_tspan[i], initial_tsh, euthyroid_dose, initial_dose)
        schneider_logl += log(2π) / 2 + log(σ) + (y - μ)^2 / 2 / σ^2
    end
    verbose && println("schneider's negative loglikelihood = $schneider_logl")
    total_neg_logl += schneider_logl
    #
    # Return final error
    #
    return total_neg_logl
end
# gives 400 mcg of oral T4 at hour 24
function blakesley_condition(u, t, integrator)
    return t - 24.0
end
# gives T3 dose at hour 0
function jonklaas_condition(u, t, integrator)
    return t - 0.01 #cannot make this exactly 0
end
# define function for adding dose
function add_dose!(integrator)
    integrator.u[10] += integrator.p[55]
    integrator.u[12] += integrator.p[56]
end
function blakesley_t4_neg_logl(sol, time, data, Vp, σ)
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
function blakesley_t3_neg_logl(sol, time, data, Vp, σ)
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
function blakesley_tsh_neg_logl(sol, time, data, Vtsh, σ)
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
    end
    return tot_loss
end
function jonklaas_t3_neg_logl(sol, time, data, Vp, σ)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        n = length(time)
        for i in 1:n
            T3_predicted = sol(time[i]) * 651.0 / Vp
            tot_loss += (T3_predicted - data[i])^2
        end
        tot_loss /= 2σ^2
        tot_loss += n * log(2π) / 2 + n * log(σ)
    end
    return tot_loss
end
# calculates mean passing dose (mcg/kg) and returns negative loglikelihood of this observation
function schneider_dose_logobs( 
    current_iter::Vector,
    fitting_index::Vector, 
    μ::Float64, # mean of euthyroid dose for a group
    σ::Float64, # std  of euthyroid dose for a group
    height::Float64, 
    weight::Float64, 
    sex::Bool,
    tspan::Float64,
    initial_tsh::Float64,
    euthyroid_dose::Float64,
    initial_dose::Float64
    )

    possible_doses = 50.0:12.5:325.0
    mindose = 0.0
    maxdose = 0.0

    #initialize simulation parameters
    scale_Vp = true
    dial  = [0.0; 0.88; 0.0; 0.88]
    ic, p = initialize(dial, scale_Vp, height, weight, sex)
    ic[7] = initial_tsh
    tot_loss = zero(Int)
    cbk = PeriodicCallback(add_dose!, 24.0)# function to add dose

    # update parameter for fitting 
    p[fitting_index] .= current_iter

    # check doses from bottom up
    for dose in possible_doses
        p[55] = dose / 777.0
        prob  = ODEProblem(thyrosim,ic,(0.0, tspan),p,callback=cbk)
        sol   = solve(prob, save_idxs=[1, 7])

        # save doses that work
        if TSH_within_interval(sol, p[48])
            mindose = dose
            break
        end
    end
    if mindose == 0 
        mindose = Inf
    end

    # check doses from top to bottom
    for dose in Iterators.reverse(possible_doses)
        p[55] = dose / 777.0
        prob  = ODEProblem(thyrosim,ic,(0.0, tspan),p,callback=cbk)
        sol   = solve(prob, save_idxs=[1, 7])

        # save doses that work
        if TSH_within_interval(sol, p[48])
            maxdose = dose
            break
        end
    end
    if maxdose == 0 
        maxdose = Inf
    end

    y = (maxdose + mindose) / 2 / weight
    return log(2π) / 2 + log(σ) + (y - μ)^2 / 2 / σ^2
end
function schneider_end_tsh( 
    current_iter::Vector,
    fitting_index::Vector,
    height::Float64, 
    weight::Float64, 
    sex::Bool,
    tspan::Float64,
    initial_tsh::Float64,
    euthyroid_dose::Float64,
    initial_dose::Float64
    )
    # quick return
    euthyroid_dose < initial_dose && return 0.0

    #initialize simulation parameters
    scale_Vp = true
    dial  = [0.0; 0.88; 0.0; 0.88]
    ic, p = initialize(dial, scale_Vp, height, weight, sex)
    ic[7] = initial_tsh
    cbk = PeriodicCallback(add_dose!, 24.0)# function to add dose

    # update parameter for fitting 
    p[fitting_index] .= current_iter

    # run ODE simulation
    p[55] = initial_dose / 777.0
    prob  = ODEProblem(thyrosim,ic,(0.0, tspan),p,callback=cbk)
    sol   = solve(prob, save_idxs=[1, 7])

    # return observed value
    return sol.u[end][2] * 5.6 / p[48]
#     return log(2π) / 2 + log(σ) + (y - initial_tsh)^2 / 2 / σ^2
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
    # fitting all parameters suggested by mauricio
    fitting_index = 
        [1; 
        13; 15; 17; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;     
        37; 40; 41; 42; 
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters
#     initial_guess = [ # original thyrosim params
#         0.00174155; 
#         0.00998996; 6.63*10^-4; 0.00074619; 3.3572*10^-4;
#         101; 47.64; 0.0; 0.53; 0.226; 
#         0.118; 0.037; 0.0034; 5; 
#         4.57; 3.90; 11.0; 5.0; 3.5; 8.0;
#         1.0; 1.0; 1.0]
    initial_guess = [ # start with good blakesley fit
        0.0023492162046331597, 
        0.015112938102610458, 0.013859639854710276, 0.00030393692427532434, 0.0003794439939511424, 
        34.17522329125075, 49.14968552336458, 0.04177085881286718, 0.08492645933069559, 0.15237583045798547, 
        0.21771547250949363, 0.04702635602348462, 0.03715391750986074, 7.285777789184459, 
        3.7112194422859854, 4.069201730462423, 0.6010144440857, 5.873363376140577, 2.7040412393908153, 9.782254820851314, 
        3.4320284272461308, 0.1328734399674661, 0.23159077159261604, 1.0]
    lowerbound = zeros(length(initial_guess))
    upperbound = initial_guess .* 10.0

    # blakesley setup 
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    
    # schneider setup
    train, test, toy = schneider_data();
    train_data = train
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")]) # 1 is female
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])
    
    return optimize(p -> objective(p, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, patient_t4, 
        patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=false), 
        initial_guess, NelderMead(), 
        Optim.Options(time_limit = 54000.0, iterations = 10000, g_tol=1e-5))
end

function prefit_error()
    # fitting all parameters suggested by mauricio
    fitting_index = 
        [1; 
        13; 15; 17; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;     
        37; 40; 41; 42; 
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters
#     initial_guess = [ # original thyrosim params
#         0.00174155; 
#         0.00998996; 6.63*10^-4; 0.00074619; 3.3572*10^-4;
#         101; 47.64; 0.0; 0.53; 0.226; 
#         0.118; 0.037; 0.0034; 5; 
#         4.57; 3.90; 11.0; 5.0; 3.5; 8.0;
#         1.0; 1.0; 1.0]
    initial_guess = [ # start with good blakesley fit
        0.0023492162046331597, 
        0.015112938102610458, 0.013859639854710276, 0.00030393692427532434, 0.0003794439939511424, 
        34.17522329125075, 49.14968552336458, 0.04177085881286718, 0.08492645933069559, 0.15237583045798547, 
        0.21771547250949363, 0.04702635602348462, 0.03715391750986074, 7.285777789184459, 
        3.7112194422859854, 4.069201730462423, 0.6010144440857, 5.873363376140577, 2.7040412393908153, 9.782254820851314, 
        3.4320284272461308, 0.1328734399674661, 0.23159077159261604, 1.0]

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data();
    train_data = train
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])

    lowerbound = zeros(length(initial_guess))
    upperbound = initial_guess .* 10.0

    return objective(initial_guess, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, 
        patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=true)
end

function postfit_error(minimizer)
    # need to know fitting index
    fitting_index = 
        [1; 
        13; 15; 17; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;     
        37; 40; 41; 42; 
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data();
    train_data = train
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])

    lowerbound = zeros(length(fitting_index))
    upperbound = minimizer .* 10.0

    return objective(minimizer, fitting_index, lowerbound, upperbound, 
        blakesley_time, my400_data, my450_data, my600_data, jonklaas_time, 
        patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
        height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=true)
end

println("Threads = ", Threads.nthreads())
println("prefit error: ")
prefit = prefit_error()
println("total prefit error = $prefit \n")

result = fit_all()

println("postfit error:")
post = postfit_error(result.minimizer)
println("total postfit error = $post \n")

println("result:")
println(result)
println(result.minimizer)
