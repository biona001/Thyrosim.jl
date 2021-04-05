module Thyrosim

    export initialize, initialize_original_thyrosim
    export thyrosim, original_thyrosim, output_plot
    export blakesley_data, jonklaas_data, jonklaas_data_new, schneider_data
    export initialize!, set_patient_ic!, find_patient_ic!
    export plot_jonklaas, plot_blakesley, plot_jonklaas_T3only
    export simulate
    export fat_free_mass, adipose_tissue_free_mass

    using CSV
    using DelimitedFiles
    using DifferentialEquations
    using DataFrames
    
    import Plots:plot, hline!, scatter!

    include("thyrosim_odes.jl")
    include("utilities.jl")

    # data directory, data not publically availble. 
    datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
