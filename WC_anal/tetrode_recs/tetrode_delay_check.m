niqf = 2016/2;

[b,a] = butter(2,[0.1/niqf 2/niqf]);

lf12_f = filtfilt(b,a,lf12_data);
lf2_f = filtfilt(b,a,lf2_data);

maxlag = 2016/2;

[datcor,lags] = xcov(lf12_f,lf2_f,maxlag,'coeff');

plot(lags/2016,datcor);grid;shg