from ar_model_fitting import ar_from_csv
arcoef = ar_from_csv(".", "neuromod_output_2024-10-10_12-03-01.csv")
print(arcoef)