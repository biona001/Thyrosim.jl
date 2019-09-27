function blakesley_data()
    
    # data path
    datapath = normpath(Thyrosim.datadir())
    data_400_path = datapath * "/blakesley/blakeley_data_400.xlsx"
    data_450_path = datapath * "/blakesley/blakeley_data_450.xlsx"
    data_600_path = datapath * "/blakesley/blakeley_data_600.xlsx"

    # import excel data
    data_400 = readxl(data_400_path, "Sheet1!A2:D34")
    data_450 = readxl(data_450_path, "Sheet1!A2:D34")
    data_600 = readxl(data_600_path, "Sheet1!A2:D34")
    my_time = data_400[:, 1]
    
    T4_400 = data_400[:, 2]*10.0
    T3_400 = data_400[:, 3]
    TSH_400 = data_400[:, 4]
    my400_data = [T4_400 T3_400 TSH_400]

    T4_450 = data_450[:, 2]*10.0
    T3_450 = data_450[:, 3]
    TSH_450 = data_450[:, 4]
    my450_data = [T4_450 T3_450 TSH_450]
    
    T4_600 = data_600[:, 2]*10.0
    T3_600 = data_600[:, 3]
    TSH_600 = data_600[:, 4]
    my600_data = [T4_600 T3_600 TSH_600]

    # column order = [T4 T3 TSH]
    my_time = convert(Vector{Float64}, my_time)
    my400_data = convert(Matrix{Float64}, my400_data)
    my450_data = convert(Matrix{Float64}, my450_data)
    my600_data = convert(Matrix{Float64}, my600_data)

    return my_time, my400_data, my450_data, my600_data
end
# my_time, my400_data, my450_data, my600_data = parse_blakesley()

function jonklaas_data()

    # data path
    datapath = normpath(Thyrosim.datadir())
    f = openxl(datapath * "/jonklaas/Daily T3 6-6-17.xlsx")

    patients      = [1 13 10 18 19 31 27 2 8 20 11 6 22 24 5]
    patient_t4    = zeros(Float64, 15, 10) # each row is a patient
    patient_t3    = zeros(Float64, 15, 10) # each row is a patient
    patient_tsh   = zeros(Float64, 15, 10) # each row is a patient
    patient_param = zeros(Float64, 15, 3)  # height (cm), weight (kg), sex

    counter = 1
    for i in patients

        # store patient characterisitcs
        str = "$i" * "!A1:B14"
        patient = readxl(f, str)
        sex = patient[1, 1]
        patient_param[counter, 1] = patient[13, 2]
        patient_param[counter, 2] = patient[14, 2]
        patient_param[counter, 3] = (sex == "Male" ? 1 : 0)

        # store initial TSH, T4, T3 (note: this is Wk-6-30, NOT baseline)
        str = "$i" * "!H2:H4"
        baseline = readxl(f, str)
        patient_tsh[counter, 1] = baseline[1]
        patient_t4[counter, 1] = 10baseline[2]   # multiplying by 10 to match units in thyrosim
        patient_t3[counter, 1] = 0.01baseline[3] # dividing by 100 to match units in thyrosim

        # store data from 30 min onwards TSH, T4, T3
        str = "$i" * "!J2:R4"
        data = readxl(f, str) 
        patient_tsh[counter, 2:end] .= data[1, :]
        patient_t4[counter, 2:end] .= 10data[2, :]
        patient_t3[counter, 2:end] .= 0.01data[3, :]

        # move to next row
        counter += 1
    end

    return patient_param, patient_t4, patient_t3, patient_tsh
end
# patient_param, patient_t4, patient_t3, patient_tsh = jonklaas_data()

function schneider_data()

    # data path
    datapath = normpath(Thyrosim.datadir())

    # toy = readdlm("schneider_train_15patients.csv", ',', header=true)
    # test = readdlm("schneider_test_154patients.csv", ',', header=true)
    # train = readdlm("schneider_train_400patients.csv", ',', header=true)

    toy   = CSV.read(datapath * "/schneider/schneider_train_15patients.csv", delim=',')
    test  = CSV.read(datapath * "/schneider/schneider_test_154patients.csv", delim=',')
    train = CSV.read(datapath * "/schneider/schneider_train_400patients.csv", delim=',')

    return train, test, toy
end
# train, test, toy = schneider_data()

function output_plot(sol)
    ## Need to change to pick better y limits!
    p = sol.prob.p 
    p1 = plot(sol.t / 24.0, 777.0 * sol[1, :] / p[47], ylim=(0, 115), label="",
       ylabel="T4", title="Thyrosim simulation")
    p1 = hline!([45, 120], label= "")
    
    p2 = plot(sol.t / 24.0, 651.0 * sol[4, :] / p[47], ylim=(0, 4), label="", 
       ylabel="T3")
    p2 = hline!([0.6, 1.8], label= "")
    
    p3 = plot(sol.t / 24.0, 5.6 * sol[7, :] / p[48], ylim=(0, 10), label="",
       ylabel="TSH", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    
    plot(p1, p2, p3, layout=(3, 1))
end

function plot_blakesley(sol, which="400")
    markersize = 2
    t_data, data400, data450, data600 = blakesley_data()
    if which == "400"
        data = data400
    elseif which == "450"
        data = data450
    else
        data = data600
    end
    
    t_data = t_data / 24.0

    ## Need to change to pick better y limits!
    p = sol.prob.p 
    p1 = plot(sol.t / 24.0, 777.0 * sol[1, :] / p[47], ylim=(0, 140), label="",
       ylabel="T4", title="Thyrosim simulation (Blakesley data)")
    p1 = hline!([45, 105], label= "")
    p1 = scatter!(t_data, data[:, 1], label="", markersize=markersize)
    
    p2 = plot(sol.t / 24.0, 651.0 * sol[4, :] / p[47], ylim=(0, 4), label="", 
       ylabel="T3")
    p2 = hline!([0.6, 1.8], label= "")
    p2 = scatter!(t_data, data[:, 2], label="", markersize=markersize)
    
    p3 = plot(sol.t / 24.0, 5.6 * sol[7, :] / p[48], ylim=(0, 10), label="",
       ylabel="TSH", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(t_data, data[:, 3], label="", markersize=markersize)
    
    plot(p1, p2, p3, layout=(3, 1))
end
