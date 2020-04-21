using Distributed
addprocs(12)

@everywhere begin
    using DifferentialEquations
    using Thyrosim
    using DiffEqCallbacks
    using DiffEqParamEstim
    using SharedArrays
end

using Optim

function objective(
    p_being_optimized::Vector, 
    fitting_index::SharedArray, 
    blakesley_time::Vector,
    blakesley_my400_data::Matrix, 
    blakesley_my450_data::Matrix, 
    blakesley_my600_data::Matrix,
    jonklaas_time::Vector, 
    jonklaas_patient_t4::Matrix, 
    jonklaas_patient_t3::Matrix, 
    jonklaas_patient_tsh::Matrix, 
    jonklaas_patient_param::Matrix, 
    jonklaas_patient_dose::Matrix,
    schneider_height::SharedArray, 
    schneider_weight::SharedArray, 
    schneider_sex::SharedArray, 
    schneider_tspan::SharedArray, 
    schneider_init_tsh::SharedArray, 
    schneider_euthy_dose::SharedArray, 
    schneider_init_dose::SharedArray;      
    verbose::Bool = false #set to true to display intermediate errors
    )
    
    total_scale_error = 0.0
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
    T4_init_400, T3_init_400, TSH_init_400 = blakesley_my400_data[1, :]
    T4_init_450, T3_init_450, TSH_init_450 = blakesley_my450_data[1, :]
    T4_init_600, T3_init_600, TSH_init_600 = blakesley_my600_data[1, :]
    # solve different ODE problems for varying doses
    p_400[55] = 400.0 / 777.0
    p_450[55] = 450.0 / 777.0
    p_600[55] = 600.0 / 777.0
    set_patient_ic!(ic, p_400, T4_init_400, T3_init_400, TSH_init_400, steady_state=true, set_tsh_lag=true)
    set_patient_ic!(ic, p_450, T4_init_450, T3_init_450, TSH_init_450, steady_state=true, set_tsh_lag=true)
    set_patient_ic!(ic, p_600, T4_init_600, T3_init_600, TSH_init_600, steady_state=true, set_tsh_lag=true)
    prob_400 = ODEProblem(thyrosim,ic,tspan,p_400,callback=cbk)
    prob_450 = ODEProblem(thyrosim,ic,tspan,p_450,callback=cbk)
    prob_600 = ODEProblem(thyrosim,ic,tspan,p_600,callback=cbk)
    sol_400 = solve(prob_400, save_idxs=[1, 4, 7])
    sol_450 = solve(prob_450, save_idxs=[1, 4, 7])
    sol_600 = solve(prob_600, save_idxs=[1, 4, 7])
    T4_error = blakesley_t4_error(sol_400, blakesley_time, blakesley_my400_data, p[47]) + 
               blakesley_t4_error(sol_450, blakesley_time, blakesley_my450_data, p[47]) + 
               blakesley_t4_error(sol_600, blakesley_time, blakesley_my600_data, p[47])
    T3_error = blakesley_t3_error(sol_400, blakesley_time, blakesley_my400_data, p[47]) + 
               blakesley_t3_error(sol_450, blakesley_time, blakesley_my450_data, p[47]) + 
               blakesley_t3_error(sol_600, blakesley_time, blakesley_my600_data, p[47])
    TSH_error = blakesley_tsh_error(sol_400, blakesley_time, blakesley_my400_data, p[48]) + 
                blakesley_tsh_error(sol_450, blakesley_time, blakesley_my450_data, p[48]) + 
                blakesley_tsh_error(sol_600, blakesley_time, blakesley_my600_data, p[48])
    blakesley_err = 0.01T4_error + T3_error + TSH_error
    scaled_blakesley_error = blakesley_err / 297 # divide total error by number of data 198 = T4 and TSH, 297 = T4,T3,TSH 
    verbose && println("blakesley error: unscaled = $blakesley_err, scaled = $scaled_blakesley_error")
    total_scale_error += scaled_blakesley_error
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
        jonklaas_err += jonklaas_error(sol, jonklaas_time, jonklaas_patient_t3[i, :], p[47])
    end
    jonklaas_err *= 10.0 # scale by 10 to match schneider&blakesley error range
    scaled_jonklaas_error = jonklaas_err / 135.0
    verbose && println("jonklaas error: unscaled = $jonklaas_err, scaled = $scaled_jonklaas_error")
    total_scale_error += scaled_jonklaas_error
    #
    # Schneider
    #
    num_params = length(p_being_optimized)
    num_sample = length(schneider_height)
    num_params == length(fitting_index) || error("check parameter length bro")
    current_iter = SharedArray{Float64}(p_being_optimized)
    schneider_err = @sync @distributed (+) for i in 1:num_sample
        one_simulation(current_iter, fitting_index, schneider_height[i], 
            schneider_weight[i], schneider_sex[i], schneider_tspan[i], schneider_init_tsh[i], 
            schneider_euthy_dose[i], schneider_init_dose[i])
    end
    scaled_schneider_err = schneider_err / num_sample
    verbose && println("schneider error: unscaled = $schneider_err, scaled = $scaled_schneider_err")
    total_scale_error += scaled_schneider_err
    #
    # Return final error
    #
    return total_scale_error
end

@everywhere function one_simulation(
    current_iter::SharedArray,
    fitting_index::SharedArray, 
    height::Float64, 
    weight::Float64, 
    sex::Bool,
    tspan::Float64,
    initial_tsh::Float64,
    euthyroid_dose::Float64,
    initial_dose::Float64
    )
    
    #initialize simulation parameters
    scale_Vp = true
    dial  = [0.0; 0.88; 0.0; 0.88]
    ic, p = initialize(dial, scale_Vp, height, weight, sex)
    ic[7] = initial_tsh
    tot_loss = zero(Int)
    cbk = PeriodicCallback(add_dose!, 24.0)# function to add dose

    # update parameter for fitting 
    p[fitting_index] .= current_iter
    
    # calculate error for euthyroid dose
    p[55] = euthyroid_dose / 777.0
    prob  = ODEProblem(thyrosim,ic,(0.0, tspan),p,callback=cbk)
    sol   = solve(prob, save_idxs=7)
    
    #increment error
    tot_loss += compute_euthyroid_dose_l2_error(sol, p[48])
    
    # when initial dose != euthyroid dose, calculate error
    if initial_dose != euthyroid_dose
        p[55] = initial_dose / 777.0
        prob = ODEProblem(thyrosim,ic,(0.0, tspan),p,callback=cbk)
        sol = solve(prob, save_idxs=7)
        tot_loss += compute_initial_dose_l2_error(sol, euthyroid_dose, initial_dose, p[48])
    end

    return tot_loss
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
@everywhere function add_dose!(integrator)
    integrator.u[10] += integrator.p[55]
    integrator.u[12] += integrator.p[56]
end
@everywhere function compute_euthyroid_dose_error(sol, Vtsh)
    tot_loss = 0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        total_hours  = sol.t[end]
        TSH_last_day = sol.u[sol.t .>= total_hours - 24] .* 5.6 ./ Vtsh
        if !all(0.5 .≤ TSH_last_day .≤ 4.5)
            tot_loss += 1
        end
    end
    return tot_loss
end
@everywhere function compute_initial_dose_error(sol, Vtsh)
    tot_loss = 0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        total_hours  = sol.t[end]
        TSH_last_day = sol.u[sol.t .>= total_hours - 24] .* 5.6 ./ Vtsh
        if all(0.5 .≤ TSH_last_day .≤ 4.5)
            tot_loss += 1
        end
    end
    return tot_loss
end
# distance to set penalty where the set C = [0.5, 4.5]
# @everywhere function compute_euthyroid_dose_l2_error(sol, Vtsh)
#     tot_loss = 0.0
#     if any((s.retcode != :Success for s in sol))
#         tot_loss = Inf
#     else
#         tsh = sol.u[end] * 5.6 / Vtsh
#         if tsh > 4.5
#             tot_loss += (tsh - 4.5)^2
#         elseif tsh < 0.5
#             tot_loss += (0.5 - tsh)^2
#         end
#     end
#     return tot_loss
# end
# distance to set penalty in log scale
@everywhere function compute_euthyroid_dose_l2_error(sol, Vtsh)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        tsh = sol.u[end] * 5.6 / Vtsh
        if tsh > 4.5
            tot_loss += log(tsh / 4.5)
        elseif tsh < 0.5
            tot_loss += log(0.5 / tsh)
        end
    end
    return tot_loss
end
# distance to set penalty where the set C = [0.0, 0.5] ∪ [4.5, Inf]                       
# @everywhere function compute_initial_dose_l2_error(sol, euthyroid_dose, initial_dose, Vtsh)
#     tot_loss = 0
#     if any((s.retcode != :Success for s in sol))
#         tot_loss = Inf
#     else
#         tsh = sol.u[end] * 5.6 / Vtsh
#         if euthyroid_dose > initial_dose && tsh < 4.5 #original TSH too high
#             tot_loss += (4.5 - tsh)^2
#         elseif euthyroid_dose < initial_dose && tsh > 0.5 #original TSH too low
#             tot_loss += (0.5 - tsh)^2
#         end
#     end
#     return tot_loss
# end
@everywhere function compute_initial_dose_l2_error(sol, euthyroid_dose, initial_dose, Vtsh)
    tot_loss = 0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        tsh = sol.u[end] * 5.6 / Vtsh
        if euthyroid_dose > initial_dose && tsh < 4.5 #original TSH too high
            tot_loss += log(4.5 / tsh)
        elseif euthyroid_dose < initial_dose && tsh > 0.5 #original TSH too low
            tot_loss += log(tsh / 0.5)
        end
    end
    return tot_loss
end
function blakesley_t4_error(sol, time, data, Vp)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        for i in 1:length(time)
            T4_predicted = sol(time[i])[1] * 777.0 / Vp
            tot_loss += (T4_predicted - data[i, 1])^2
        end
    end
    return tot_loss
end
function blakesley_t3_error(sol, time, data, Vp)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        for i in 1:length(time)
            T3_predicted = sol(time[i])[2] * 651.0 / Vp
            tot_loss += (T3_predicted - data[i, 2])^2
        end
    end
    return tot_loss
end
function blakesley_tsh_error(sol, time, data, Vtsh)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        for i in 1:length(time)
            predicted_tsh = sol(time[i])[3] * 5.6 / Vtsh
            tot_loss += (predicted_tsh - data[i, 3])^2
        end
    end
    return tot_loss
end
function jonklaas_error(sol, time, data, Vp)
    tot_loss = 0.0
    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        for i in 1:length(time)
            T3_predicted = sol(time[i]) * 651.0 / Vp
            tot_loss += (T3_predicted - data[i])^2
        end
    end
    return tot_loss
end
function fit_all()
    # initialize initial guess and fitting index
    fitting_index = SharedArray{Int}([12; 29; 28; 45; 30; 31; 49; 50; 51; 52; 53; 54])
    initial_guess = [0.0189; 0.207; 0.8892067744277633;1.6882221360501146;
        69.90379778202167;38.71161774205076; 6.039888256864343; 3.7006563259936747;
        8.748185980217668;6.590694001313398;2.896554559451672;13.013203952637502]

    # blakesley setup 
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    
    # schneider setup
    train, test, toy = schneider_data();
    train_data = toy
    height = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Ht.m")]))
    weight = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")]))
    sex    = SharedArray{Bool}(convert(Vector{Bool}, train_data[!, Symbol("Sex")]))
    tspan  = SharedArray{Float64}(convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")]))
    init_tsh   = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")]))
    euthy_dose = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")]))
    init_dose  = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")]))
    
    return optimize(p -> objective(p, fitting_index, 
                                   blakesley_time, my400_data, my450_data, my600_data,
                                   jonklaas_time, patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
                                   height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=false), initial_guess, 
        NelderMead(), Optim.Options(time_limit = 79200, iterations = 10000, g_tol=1e-5))
end

function prefit_error()
    # initialize initial guess and fitting index
    fitting_index = SharedArray{Int}([12; 29; 28; 45; 30; 31; 49; 50; 51; 52; 53; 54])
    initial_guess = [0.0189; 0.207; 0.8892067744277633;1.6882221360501146;
        69.90379778202167;38.71161774205076; 6.039888256864343; 3.7006563259936747;
        8.748185980217668;6.590694001313398;2.896554559451672;13.013203952637502]

    # blakesley setup 
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data()
    train_data = train
    height = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Ht.m")]))
    weight = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")]))
    sex    = SharedArray{Bool}(convert(Vector{Bool}, train_data[!, Symbol("Sex")]))
    tspan  = SharedArray{Float64}(convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")]))
    init_tsh   = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")]))
    euthy_dose = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")]))
    init_dose  = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")]));
    
    return objective(initial_guess, fitting_index, 
              blakesley_time, my400_data, my450_data, my600_data,
              jonklaas_time, patient_t4, patient_t3, patient_tsh, jonklaas_patient_param, jonklaas_patient_dose,
              height, weight, sex, tspan, init_tsh, euthy_dose, init_dose, verbose=true)
end

function postfit_error(minimizer)
    # need to know fitting index
    fitting_index = SharedArray{Int}([12; 29; 28; 45; 30; 31; 49; 50; 51; 52; 53; 54])

    # blakesley setup 
    blakesley_time, my400_data, my450_data, my600_data = blakesley_data()
    # jonklaas setup
    jonklaas_time = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]
    jonklaas_patient_param, jonklaas_patient_dose, patient_t4, patient_t3, patient_tsh = jonklaas_data()
    # schneider setup
    train, test, toy = schneider_data()
    train_data = train
    height = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Ht.m")]))
    weight = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("Wt.kg")]))
    sex    = SharedArray{Bool}(convert(Vector{Bool}, train_data[!, Symbol("Sex")]))
    tspan  = SharedArray{Float64}(convert(Vector{Float64}, 24train_data[!, Symbol("Days.to.euthyroid")]))
    init_tsh   = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("TSH.preop")]))
    euthy_dose = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.euthyroid.dose")]))
    init_dose  = SharedArray{Float64}(convert(Vector{Float64}, train_data[!, Symbol("LT4.initial.dose")]));
    
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
