using Revise
using ExcelReaders
using Thyrosim

# import data 
my_time, my400_data, my450_data, my600_data = blakesley_data()
patient_param, patient_t4, patient_t3, patient_tsh = jonklaas_data()
train, test, toy = schneider_data()




################# try fitting blakesley data

using Revise
using DifferentialEquations
using Thyrosim
# using Plots
using DiffEqCallbacks
using Optim

my_time, my400_data, my450_data, my600_data = blakesley_data()
ic, p = initialize() # blakesley patients are not hypothyroid 
tspan = (my_time[1], my_time[end])

function condition(u, t, integrator)
    return t - 24.0
end
function affect400!(integrator)
#     T4_dose = 1.6 * 70.0 / 777.0;
    T4_dose = 400.0 / 777.0
    T3_dose = 0.0
    
    ## CHECK: Should probably be adding to pill compartments 10, 12 but do this same way as Simon's code so
    ## we can compare.
    integrator.u[11] += T4_dose
    integrator.u[13] += T3_dose
end
function affect450!(integrator)
#     T4_dose = 1.6 * 70.0 / 777.0;
    T4_dose = 450.0 / 777.0
    T3_dose = 0.0
    
    ## CHECK: Should probably be adding to pill compartments 10, 12 but do this same way as Simon's code so
    ## we can compare.
    integrator.u[11] += T4_dose
    integrator.u[13] += T3_dose
end
function affect600!(integrator)
#     T4_dose = 1.6 * 70.0 / 777.0;
    T4_dose = 600.0 / 777.0
    T3_dose = 0.0
    
    ## CHECK: Should probably be adding to pill compartments 10, 12 but do this same way as Simon's code so
    ## we can compare.
    integrator.u[11] += T4_dose
    integrator.u[13] += T3_dose
end

# gives dose at 24 hour
cbk400 = ContinuousCallback(condition, affect400!);
cbk450 = ContinuousCallback(condition, affect450!);
cbk600 = ContinuousCallback(condition, affect600!);

# define problem and error function
ode = original_thyrosim
prob_400 = ODEProblem(ode,ic,tspan,p,callback=cbk400)
prob_450 = ODEProblem(ode,ic,tspan,p,callback=cbk400)
prob_600 = ODEProblem(ode,ic,tspan,p,callback=cbk400)
cost_function_400 = build_loss_objective(prob,Tsit5(),L2Loss(my_time, my400_data[:, 3]),maxiters=10000,verbose=false, save_idxs=[7])
cost_function_450 = build_loss_objective(prob,Tsit5(),L2Loss(my_time, my450_data[:, 3]),maxiters=10000,verbose=false, save_idxs=[7])
cost_function_600 = build_loss_objective(prob,Tsit5(),L2Loss(my_time, my600_data[:, 3]),maxiters=10000,verbose=false, save_idxs=[7])

#fits blakesley data
result = optimize(cost_function_400, p, BFGS())
result.minimizer


###############################

# fitting_index = [30, 31, 49, 50, 51, 52]
fitting_index = [30, 31]
function loss(p)
	dial = [1.0; 0.88; 1.0; 0.88]

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
	# p[30] = 1166           #Bzero
	# p[31] = 581            #Azero
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

	return cost_function(p)
end


#need cost function which accepts only 1 input vector?
function jonklaas_loss(prob::ODEProblem, data::AbstractVector, times::AbstractVector)
	err = 0.0
	sol = solve(prob, save_idxs=[7]);
	storage = zeros(Float64, 1)
	for i in 1:length(data)
		storage .= sol(times[i])
		err += (data[i] - storage[1])^2
	end
	return err	   
end
jonklaas_loss(prob, my400_data[:, 3], my_time)

L2Loss(my_time, my400_data[:, 3])










# result = optimize(cost_function, p, BFGS())

# result = optimize(cost_function, p, BFGS())



# solve ode before fitting and plot solution
# sol = solve(prob, save_idxs=[1, 4, 7]); # save_idxs allows you to only store values in certain compartments. 
# output_plot(sol)


# try optimizing with BFGS
# lower  = 
# upper  = 
# result = optimize(cost_function, p, BFGS())





####### try fitting only first parameter in only p[1]

#ode
function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - u[1]*u[2]
  du[2] = dy = p[2]*u[2] + u[1]*u[2]
end

#problem & initial condition
u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.5, -3]
prob = ODEProblem(f,u0,tspan,p)

#simulate data
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

#define cost function to calculate error in both compartments
cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t, data),
                                     maxiters=10000,verbose=false, save_idxs=[1, 2])

function closure(p)
	p[2] = -1.5
	err = cost_function(p)
	return err
end

result = optimize(cost_function, [1.0, -3.0])












