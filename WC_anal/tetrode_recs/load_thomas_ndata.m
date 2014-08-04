
lf12_data = reshape(CSC12_Samples,1,numel(CSC12_Samples));
lf1_data = reshape(CSC1_Samples,1,numel(CSC1_Samples));
lf2_data = reshape(CSC2_Samples,1,numel(CSC2_Samples));
lf12_data = detrend(lf12_data);
lf1_data = detrend(lf1_data);
lf2_data = detrend(lf2_data);
time_stamps = interp1q((1:length(CSC12_TimeStamps))',CSC12_TimeStamps',...
    (linspace(1,length(CSC12_TimeStamps),length(lf12_data)))');

time = [1:length(lf12_data)]/2016;

clear CSC*

save lfp_data time* *data