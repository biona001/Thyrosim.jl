"""
# Blakesley dataset description:

+ 36 healthy patients, 18 women 18 men, 3 were excluded
+ They were given 400/450/600 μg of oral T4 dose at 24th hour. Measurements were make for up to 5 days
+ All patients have "normal weight for height" 
"""
function blakesley_data()
    datapath = normpath(Thyrosim.datadir())

    my_time    = readdlm(datapath * "/blakesley/blakesley.time.csv", ',')
    my400_data = readdlm(datapath * "/blakesley/blakesley.400data.csv", ',')
    my450_data = readdlm(datapath * "/blakesley/blakesley.450data.csv", ',')
    my600_data = readdlm(datapath * "/blakesley/blakesley.600data.csv", ',')

    return my_time, my400_data, my450_data, my600_data
end

# function parse_blakesley_excel()
#     # data path
#     datapath = normpath(Thyrosim.datadir())
#     data_400_path = datapath * "/blakesley/blakeley_data_400.xlsx"
#     data_450_path = datapath * "/blakesley/blakeley_data_450.xlsx"
#     data_600_path = datapath * "/blakesley/blakeley_data_600.xlsx"

#     # import excel data
#     data_400 = readxl(data_400_path, "Sheet1!A2:D34")
#     data_450 = readxl(data_450_path, "Sheet1!A2:D34")
#     data_600 = readxl(data_600_path, "Sheet1!A2:D34")
#     my_time = data_400[:, 1]
    
#     T4_400 = data_400[:, 2]*10.0
#     T3_400 = data_400[:, 3]
#     TSH_400 = data_400[:, 4]
#     my400_data = [T4_400 T3_400 TSH_400]

#     T4_450 = data_450[:, 2]*10.0
#     T3_450 = data_450[:, 3]
#     TSH_450 = data_450[:, 4]
#     my450_data = [T4_450 T3_450 TSH_450]
    
#     T4_600 = data_600[:, 2]*10.0
#     T3_600 = data_600[:, 3]
#     TSH_600 = data_600[:, 4]
#     my600_data = [T4_600 T3_600 TSH_600]

#     # column order = [T4 T3 TSH]
#     my_time = convert(Vector{Float64}, my_time)
#     my400_data = convert(Matrix{Float64}, my400_data)
#     my450_data = convert(Matrix{Float64}, my450_data)
#     my600_data = convert(Matrix{Float64}, my600_data)

#     return my_time, my400_data, my450_data, my600_data
# end

# make data readable without excelreader 
# my_time, my400_data, my450_data, my600_data = parse_blakesley_excel()
# writedlm("blakesley.time.csv", my_time, ',')
# writedlm("blakesley.400data.csv", my400_data, ',')
# writedlm("blakesley.450data.csv", my450_data, ',')
# writedlm("blakesley.600data.csv", my600_data, ',')

"""
Jonklaas's patient description:

This consists of T4, T3 and TSH measurements of 15 hypothyroid patients 
(for which we know their individual height, weight and sex (mostly women)) 
under T3 replacement therapy. In an initial period, the right dose to give 
the patients was determined by adjusting the dose and measuring the patients 
for 8 weeks. Then the actual timecourse measurements begin for 8 hours after 
the last dose (either 30 or 45 ug T3, depending on the patient), every hour, 
except for the beginning measurement at 30 min.

# Output:
+ patient_param: Matrix of where columns are height (m) weight (kg) and sex (1 = male). Each row is a patient. 
+ patient_dose: T3 oral dose (in μg)
"""
function jonklaas_data()
    datapath = normpath(Thyrosim.datadir())

    patient_dose  = [30 30 45 45 30 30 45 30 30 30 45 45 45 45 30] # mcg of T3
    patient_param = readdlm(datapath * "/jonklaas/patient_param.csv", ',')
    patient_t4    = readdlm(datapath * "/jonklaas/patient_t4.csv", ',')
    patient_t3    = readdlm(datapath * "/jonklaas/patient_t3.csv", ',')
    patient_tsh   = readdlm(datapath * "/jonklaas/patient_tsh.csv", ',')

    return patient_param, patient_dose, patient_t4, patient_t3, patient_tsh
end

"""
New Jonklaas's patient description:

This data contains 4 measurements of FT4, T3 and TSH for 50 patients. These 
patients were never diagnosed with hyper/hypo-thyroidism, and their first 2 
TSH measurements were normal. 

The first 2 measurements were 1 week and 1 day before thyroidectomy, and the
last two measurements were 8 and 16 weeks after surgery. Patient weight were
measured on their first and last visit. 
"""
function jonklaas_data_new()
    datapath = normpath(Thyrosim.datadir())
    full = CSV.read(datapath * "/jonklaas/jonklass_new_data.csv", DataFrame)

    patient_t4 = convert(Matrix{Float64}, full[:, [6, 10, 14, 18]])
    patient_t3 = convert(Matrix{Float64}, full[:, [8, 12, 16, 20]]) ./ 100
    patient_tsh = convert(Matrix{Float64}, full[:, [9, 13, 17, 21]])

    dose_w1 = full[Symbol("1st dose")]
    dose_w4 = full[Symbol("Final dose")]
    patient_dose = convert(Matrix{Float64}, [dose_w1 dose_w4])

    weight_w1 = full[Symbol("Wt 1")] # KG
    weight_w4 = full[Symbol("Wt 4")] # KG
    height = full[Symbol(" Ht (cm)")] ./ 100 # convert to meters
    sex = full[:Sex] .== "M" # 1 is male, 0 is female
    patient_param = convert(Matrix{Float64}, [weight_w1 weight_w4 height sex])

    return patient_param, patient_dose, patient_t4, patient_t3, patient_tsh
end

# function parse_jonklaas_excel()

#     # data path
#     datapath = normpath(Thyrosim.datadir())
#     f = openxl(datapath * "/jonklaas/Daily T3 6-6-17.xlsx")

#     patients      = [1  13 10 18 19 31 27  2  8 20 11  6 22 24  5]
#     patient_dose  = [30 30 45 45 30 30 45 30 30 30 45 45 45 45 30]
#     patient_t4    = zeros(Float64, 15, 10) # each row is a patient
#     patient_t3    = zeros(Float64, 15, 10) # each row is a patient
#     patient_tsh   = zeros(Float64, 15, 10) # each row is a patient
#     patient_param = zeros(Float64, 15, 3)  # height (m), weight (kg), sex (1 = male)

#     counter = 1
#     for i in patients

#         # store patient characterisitcs
#         str = "$i" * "!A1:B14"
#         patient = readxl(f, str)
#         sex = patient[1, 1]
#         patient_param[counter, 1] = patient[13, 2]
#         patient_param[counter, 2] = patient[14, 2]
#         patient_param[counter, 3] = (sex == "Male" ? 1 : 0)

#         # store initial TSH, T4, T3 (note: this is Wk-6-30, NOT baseline)
#         str = "$i" * "!H2:H4"
#         baseline = readxl(f, str)
#         patient_tsh[counter, 1] = baseline[1]
#         patient_t4[counter, 1] = 10baseline[2]   # multiplying by 10 to match units in thyrosim
#         patient_t3[counter, 1] = 0.01baseline[3] # dividing by 100 to match units in thyrosim

#         # store data from 30 min onwards TSH, T4, T3
#         str = "$i" * "!J2:R4"
#         data = readxl(f, str) 
#         patient_tsh[counter, 2:end] .= data[1, :]
#         patient_t4[counter, 2:end] .= 10data[2, :]
#         patient_t3[counter, 2:end] .= 0.01data[3, :]

#         # move to next row
#         counter += 1
#     end

#     patient_param[:, 1] ./= 100 # scale height from cm to meters

#     return patient_param, patient_dose, patient_t4, patient_t3, patient_tsh
# end

# make data readable without excelreader 
# patient_param, patient_t4, patient_t3, patient_tsh = parse_jonklaas_excel()
# writedlm("patient_param.csv", patient_param, ',')
# writedlm("patient_t4.csv", patient_t4, ',')
# writedlm("patient_t3.csv", patient_t3, ',')
# writedlm("patient_tsh.csv", patient_tsh, ',')


"""
Schneider patients are completely thyroidectomized, but have attained 
euthyroidism using T4 oral treatments. Normal TSH level = 0.45-4.5 mIU/ml. 
Here 12.5 μg is the smallest dosing increment.

# Usable data:
- Starting T4 dose
- Euthyroid T4 dose
- Days to euthyroid
- Number of dose changes 
- Initial TSH level
"""
function schneider_data(;exclude_missing=true)
    # data path
    datapath = normpath(Thyrosim.datadir())

    all_data = CSV.read(datapath * "/schneider/merged_schneider.csv", 
        DataFrame, delim=',')
    train_idx = findall(!ismissing, all_data[!, Symbol("6 week TSH")])
    train_data = all_data[train_idx, :]

    if exclude_missing # only some values of TSH.preop is missing
        keep_idx_all = findall(x -> x !== missing, all_data[!, Symbol("TSH.preop")])
        keep_idx_train = findall(x -> x !== missing, train_data[!, Symbol("TSH.preop")])
        all_data = all_data[keep_idx_all, :]
        train_data = train_data[keep_idx_train, :]
    end

    return all_data, train_data
end

function output_plot(sol; title::AbstractString = "Thyrosim simulation", automargins::Bool=true)

    # parameters to adjust figure limits
    p = sol.prob.p 
    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = 777.0 * sol[1, :] / p[47]
    T3 = 651.0 * sol[4, :] / p[47]
    TSH = 5.6 * sol[7, :] / p[48]
    if automargins
        t4lim = max(1.2maximum(T4), 110.0)
        t3lim = max(1.2maximum(T3), 2.5)
        tshlim = max(1.2maximum(TSH), 5.5)
    end

    p1 = plot(sol.t / 24.0, T4, ylim=(0, t4lim), label="",
       ylabel="T4 (mcg/L)", title=title)
    p1 = hline!([45, 120], label= "")
    
    p2 = plot(sol.t / 24.0, T3, ylim=(0, t3lim), label="", 
       ylabel="T3 (mcg/L)")
    p2 = hline!([0.6, 1.8], label= "")
    
    p3 = plot(sol.t / 24.0, TSH, ylim=(0, tshlim), label="",
       ylabel="TSH (mU/L)", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    
    plot(p1, p2, p3, layout=(3, 1))
end

function plot_blakesley(sol, which::AbstractString="400"; 
    title::AbstractString = "Thyrosim simulation (Blakesley data)", automargins::Bool=true)
    markersize = 2
    t_data, data400, data450, data600 = blakesley_data()
    if which == "400"
        data = data400
    elseif which == "450"
        data = data450
    elseif which == "600"
        data = data600
    else
        error("choices for 2nd argument is 400, 450, or 600, must be inputted as strings.")
    end
    
    t_data = t_data / 24.0

    # parameters to adjust figure limits
    p = sol.prob.p 
    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = 777.0 * sol[1, :] / p[47]
    T3 = 651.0 * sol[4, :] / p[47]
    TSH = 5.6 * sol[7, :] / p[48]
    if automargins
        t4lim = max(1.2maximum(T4), 110.0)
        t3lim = max(1.2maximum(T3), 2.5)
        tshlim = max(1.2maximum(TSH), 5.5)
    end

    ## Need to change to pick better y limits!
    p1 = plot(sol.t / 24.0, T4, ylim=(0, t4lim), label="",
       ylabel="T4 (mcg/L)", title=title)
    p1 = hline!([45, 105], label= "")
    p1 = scatter!(t_data, data[:, 1], label="", markersize=markersize)
    
    p2 = plot(sol.t / 24.0, T3, ylim=(0, t3lim), label="", 
       ylabel="T3 (mcg/L)")
    p2 = hline!([0.6, 1.8], label= "")
    p2 = scatter!(t_data, data[:, 2], label="", markersize=markersize)
    
    p3 = plot(sol.t / 24.0, TSH, ylim=(0, tshlim), label="",
       ylabel="TSH (mU/L)", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(t_data, data[:, 3], label="", markersize=markersize)
    
    plot(p1, p2, p3, layout=(3, 1))
end

function plot_jonklaas(sol, T4data::Vector, T3data::Vector, TSHdata::Vector; 
        title::AbstractString="Thyrosim simulation (Jonklaas data)", automargins::Bool=true)

    markersize = 2
    #time in which data are measured
    t_data = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0] 

    # parameters to adjust figure limits
    p = sol.prob.p 
    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = 777.0 * sol[1, :] / p[47]
    T3 = 651.0 * sol[4, :] / p[47]
    TSH = 5.6 * sol[7, :] / p[48]
    if automargins
        t4lim = max(1.2maximum(T4), 110.0)
        t3lim = max(1.2maximum(T3), 2.5)
        tshlim = max(1.2maximum(TSH), 5.5)
    end

    ## Need to change to pick better y limits!
    p1 = plot(sol.t, T4, ylim=(0, t4lim), label="",
       ylabel="T4 (mcg/L)", title=title)
    p1 = hline!([45, 105], label= "")
    p1 = scatter!(t_data, T4data, label="", markersize=markersize)
    
    p2 = plot(sol.t, T3, ylim=(0, t3lim), label="", 
       ylabel="T3 (mcg/L)")
    p2 = hline!([0.6, 1.8], label= "")
    p2 = scatter!(t_data, T3data, label="", markersize=markersize)
    
    p3 = plot(sol.t, TSH, ylim=(0, tshlim), label="",
       ylabel="TSH (mU/L)", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(t_data, TSHdata, label="", markersize=markersize)
    
    plot(p1, p2, p3, layout=(3, 1))
end

function plot_jonklaas_T3only(sol, T3data::Vector; 
    title::AbstractString="Thyrosim simulation (Jonklaas data)", automargins::Bool=true)
    markersize = 2

    #time in which data are measured
    t_data = [0.0; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0] 

    # parameters to adjust figure limits
    p = sol.prob.p 
    t3lim = 4
    T3 = 651.0 * sol[4, :] / p[47]
    if automargins
        t3lim = max(1.2maximum(T3), 2.5)
    end

    ## Need to change to pick better y limits!
    p2 = plot(sol.t, T3, ylim=(0, t3lim), label="", 
       ylabel="T3 (mcg/L)", title=title)
    p2 = hline!([0.6, 1.8], label= "")
    p2 = scatter!(t_data, T3data, label="", markersize=markersize)

    plot(p2)
end

"""
    simulate(h, w, sex, ...)

Simulate a person of known height, weight, and gender for 30 days (default).

If `warmup = true`, will first run the model for 30 days, assuming healthy
thyroid function, to get approximate initial condition. 
"""
function simulate(
    h::Float64, # units meters
    w::Float64, # units kg
    sex::Bool; # true = male, false = female
    days::Int=30, 
    dial=[1.0; 0.88; 1.0; 0.88], 
    T4dose::Float64=0.0, # mcgs
    T3dose::Float64=0.0, # mcgs
    dosing_interval::Float64=24.0, #hours
    warmup::Bool = true,
    fitting_index = Int[],
    parameters = Float64[],
    )
    function add_dose!(integrator)
        integrator.u[10] += integrator.p[55]
        integrator.u[12] += integrator.p[56]
    end
    cbk = PeriodicCallback(add_dose!, dosing_interval) 

    # initialize thyrosim parameters
    ic, p = initialize([1.0; 0.88; 1.0; 0.88], true, h, w, sex)
    p[fitting_index] .= parameters

    # run simulation for 30 days to get approximate steady state conditions
    # this assumes healthy patient without dose
    warmup && find_patient_ic!(ic, p, 30) 

    # setup daily dosing and fitting parameters 
    p[55] = T4dose / 777.0 # daily dose
    p[56] = T3dose / 651.0 # daily dose
    p[57:60] .= dial #set dial
    
    # solve and return ode solution
    prob = ODEProblem(thyrosim,ic,(0.0, 24days),p,callback=cbk)
    return solve(prob)
end
