function chi2pdof = chi2pdof_pmf(x,emp_counts,pred_counts,n_model_params)

dof = length(x) - n_model_params - 1;
chi2 = sum((emp_counts - pred_counts).^2./pred_counts);
chi2pdof = chi2/dof;