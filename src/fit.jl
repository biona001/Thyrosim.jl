function reset_p!(
	p::Vector{Float64}, 
	fitting_index::Vector{Int},
	dial::Vector{Float64} = [1.0; 0.88; 1.0; 0.88]
	)

	ic, p_init = initialize(dial)

	for i in 1:length(p)
		if i in fitting_index 
			continue
		else
			p[i] = p_init[i]
		end
	end

	return nothing
end

"""
Creates a new ODE model where all parameters not being fitted will become constants. 
"""
function create_thyrosim_odes_for_fitting(fitting_index::Vector{Int})
    # maybe the Problem Generator Function will work:
    # http://docs.juliadiffeq.org/v5.0.0/analysis/parameter_estimation.html#The-Problem-Generator-Function-1
    # otherwise not sure how to make this completely general
end

"""
Simulate Schneider patients for a fixed period of time.
"""
function simulate_schneider_patient(ode::Function, tspan::Tuple{Float64, Float64})
	# 1. generate ODEProblem (e.g. prob = ODEProblem(original_thyrosim,ic,tspan,p,callback=cbk))
	# 2. solve the problem (e.g. sol = solve(prob))
	# 3. grab the last timepoint of sol 
end
