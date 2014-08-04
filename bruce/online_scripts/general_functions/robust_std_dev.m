function rstd = robust_std_dev(data)

data = data - nanmedian(data);
rstd = 1.483*nanmedian(abs(data));
