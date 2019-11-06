__precompile__()

module Thyrosim

    export initialize, initialize_original_thyrosim
    export thyrosim, original_thyrosim
    export blakesley_data, jonklaas_data, schneider_data
    export initialize!, set_patient_ic!

	using ExcelReaders, CSV

    include("thyrosim_odes.jl")
    include("utilities.jl")

    # data directory, data not publically availble. 
	datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
