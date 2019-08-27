__precompile__()

module Thyrosim

    export initialize, initialize_original_thyrosim
    export thyrosim, original_thyrosim
    export blakesley_data, jonklaas_data, schneider_data
    export reset_p!

	using ExcelReaders, CSV

    include("thyrosim_odes.jl")
    include("utilities.jl")
    include("fit.jl")

    # data directory, data not publically availble. 
	datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
