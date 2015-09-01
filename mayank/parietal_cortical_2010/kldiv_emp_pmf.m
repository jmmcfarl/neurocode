function kldiv = kldiv_emp_pmf(emp_pmf,pred_pmf)

used_bins = find(emp_pmf > 0);
kldiv = sum(emp_pmf(used_bins).*log2(emp_pmf(used_bins)./pred_pmf(used_bins)));