function rstd = robust_std_dev(data)

if size(data,2) > 1
    data = bsxfun(@minus,data,nanmedian(data));
else
    data = data - nanmedian(data);
end
rstd = 1.483*nanmedian(abs(data));
