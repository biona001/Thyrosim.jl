__precompile__()

module Thyrosim

    export initialize, original_thyrosim
    export blakesley_data, jonklaas_data, schneider_data

	using ExcelReaders, CSV

    include("thyrosim_odes.jl")
    include("utilities.jl")

    # data directory, data not publically availble. 
	datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
