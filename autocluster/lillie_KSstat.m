function KSstat = lillie_KSstat(x)

outlier_thresh = 3;
used_samps = find(abs(zscore(x)) < outlier_thresh);
[sampleCDF,xCDF] = ecdf(x(used_samps));
xCDF = xCDF(2:end);
nullCDF = normcdf(xCDF, mean(x(used_samps)), std(x(used_samps)));
delta1   = sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
delta2   = sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
deltaCDF = abs([delta1 ; delta2]);
KSstat = max(deltaCDF);
