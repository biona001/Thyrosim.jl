using DifferentialEquations
using Thyrosim
using DiffEqCallbacks
using DiffEqParamEstim
using SharedArrays
using Optim

function objective(
    p_being_optimized::Vector,
    fitting_index::Vector,
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
    any(p_being_optimized .< 0.0) && return Inf
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
    schneider_errors = zeros(Threads.nthreads()*8) # spacing by 64 bytes to avoid false sharing
    Threads.@threads for i in 1:num_sample
        #calculate (y - μ)^2 in mcg/kg
        err = abs2(schneider_euthy_dose[i] / schneider_weight[i] -
            mean_passing_dose(p_being_optimized, fitting_index, schneider_height[i],
            schneider_weight[i], schneider_sex[i], schneider_tspan[i], schneider_init_tsh[i],
            schneider_euthy_dose[i], schneider_init_dose[i]))

        id = Threads.threadid()
        schneider_errors[id*8] += err
    end
    schneider_err = sum(schneider_errors) / (2 * abs2(p[64]))
    schneider_err += num_sample * log(2π) / 2 - num_sample * log(p[64])
    verbose && println("schneider's negative loglikelihood = $schneider_err")
    total_neg_logl += schneider_err
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
        tot_loss += n * log(2π) / 2 - n * log(σ)
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
        tot_loss += n * log(2π) / 2 - n * log(σ)
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
        tot_loss += n * log(2π) / 2 - n * log(σ)
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
        tot_loss += n * log(2π) / 2 - n * log(σ)
    end
    return tot_loss
end
function mean_passing_dose( # in mcg/kg
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

    possible_doses = collect(50.0:12.5:325.0)
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
    reverse!(possible_doses)
    for dose in possible_doses
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

    return (maxdose + mindose) / 2 / weight
end
function TSH_within_interval(sol, Vtsh)
    tot_loss = 0
    tsh = sol.u[end][2] * 5.6 / Vtsh
    if 0.5 ≤ tsh ≤ 4.5
        return true
    end
    return false
end

function fit_all()
    fitting_index = [
        1; 11;
        13; 14; 15; 16; 17; 18; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;
        36; 37; 40; 41; 42; 44;
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters
    initial_guess = [
        0.00174155; 0.88;
        0.00998996; 2.85; 6.63*10^-4; 95; 0.00074619; 0.075; 3.3572*10^-4;
        101; 47.64; 0.0; 0.53; 0.226;
        23; 0.118; 0.037; 0.0034; 5; 0.12;
        4.57; 3.90; 11.0; 5.0; 3.5; 8.0;
        1.0; 1.0; 1.0; 1.0]

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()

    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()

    # schneider setup
    train, test, toy = schneider_data();
    train_data = toy
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])

    return optimize(p -> objective(p, fitting_index,
                                   blakesley_time, my400_data, my450_data, my600_data,
                                   jonklaas_time, patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
                                   height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=false), initial_guess,
        NelderMead(), Optim.Options(time_limit = 80000.0, iterations = 10000, g_tol=1e-5))
end

function prefit_error()
    fitting_index = [
        1; 11;
        13; 14; 15; 16; 17; 18; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;
        36; 37; 40; 41; 42; 44;
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters
    initial_guess = [
        0.00174155; 0.88;
        0.00998996; 2.85; 6.63*10^-4; 95; 0.00074619; 0.075; 3.3572*10^-4;
        101; 47.64; 0.0; 0.53; 0.226;
        23; 0.118; 0.037; 0.0034; 5; 0.12;
        4.57; 3.90; 11.0; 5.0; 3.5; 8.0;
        1.0; 1.0; 1.0; 1.0]

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data();
    train_data = toy
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])

    return objective(initial_guess, fitting_index,
              blakesley_time, my400_data, my450_data, my600_data,
              jonklaas_time, patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
              height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=true)
end

function postfit_error(minimizer)
    # need to know fitting index
    fitting_index = [
        1; 11;
        13; 14; 15; 16; 17; 18; 19; # T4 -> T3 conversion
        30; 31; 32; 34; 35;
        36; 37; 40; 41; 42; 44;
        49; 50; 51; 52; 53; 54;  # hill function parameters
        61; 62; 63; 64]          # variance parameters

    # blakesley setup
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data();
    train_data = toy
    height = convert(Vector{Float64}, train_data[!, Symbol("Ht.m")])
    weight = convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")])
    sex    = convert(Vector{Bool}, train_data[!, Symbol("Sex")])
    tspan  = convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")])
    init_tsh   = convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")])
    euthy_dose = convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")])
    init_dose  = convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")])

    return objective(minimizer, fitting_index,
              blakesley_time, my400_data, my450_data, my600_data,
              jonklaas_time, patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
              height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=true)
end

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
