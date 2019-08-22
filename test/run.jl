using Revise
using ExcelReaders
using Thyrosim

# import data
my_time, my400_data, my450_data, my600_data = blakesley_data()
patient_param, patient_t4, patient_t3, patient_tsh = jonklaas_data()
train, test, toy = schneider_data()
