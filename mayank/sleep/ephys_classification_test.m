clear all
close all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs data
% load sleep_dirs_old

% spike_thresh = 10;
%%
dd =29;
% for dd = 1:length(data);
cd(data(dd).dir)
pwd
load procData heka_data heka_time mp_d mp_t ipsi_csc csc_time contra_csc

if length(heka_data) > 0
heka_deriv = [0; diff(heka_data)];
heka_deriv = smooth(heka_deriv,3);
spike_thresh = prctile(heka_deriv,99); %quian-quiroga threshold
% heka_zderiv = zscore(heka_deriv);
spk_times = find([diff(sign([0;diff(heka_deriv)]));0] < 0 & heka_deriv >= spike_thresh);
% spk_times = find([diff(sign([0;diff(heka_zderiv)]));0] < 0);


subplot(3,1,1);
plot(heka_time,heka_deriv); hold on
plot(heka_time(spk_times),heka_deriv(spk_times),'ro');

avg_spk_rate(dd) = length(spk_times)/range(heka_time);

subplot(3,1,2)
plot(heka_time,heka_data); axis tight

heka_Fs = 1/nanmedian(diff(heka_time));
[bb,aa] = butter(4,40/(heka_Fs/2),'low');
heka_filt = filtfilt(bb,aa,heka_data);
subplot(3,1,3)
ksdensity(heka_filt,'bandwidth',0.1,'npoints',1e3);

% end
end
%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
if ~isempty(ipsi_csc)
ctx_lfp = ipsi_csc{3};
else
    ctx_lfp = contra_csc{3};
end


params.Fs = 1/median(diff(csc_time));
params.tapers = [10 19];
[bb,aa] = butter(2,0.05/(params.Fs/2),'high');
ctx_lfp = zscore(filtfilt(bb,aa,ctx_lfp));

win = 25;
[S,f] = mtspectrumsegc(ctx_lfp,win,params,1);

%%
params.Fs = 1/nanmedian(diff(mp_t));
[bb,aa] = butter(2,0.05/(params.Fs/2),'high');
mp_d = zscore(filtfilt(bb,aa,mp_d));
[S_mp,f_mp] = mtspectrumsegc(mp_d,win,params,1);

%%
figure
plot(f,log10(S),f_mp,log10(S_mp),'r')
set(gca,'xscale','log'); 
xlim([0.05 20]);

% pause
% close all
% end