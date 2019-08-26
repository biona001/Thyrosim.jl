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

