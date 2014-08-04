%% for pyramidal cells
clear all
load pyr_heka_dir

movingwin = [3 3];
params.Fs = 5000;
params.fpass = [2 1000];
params.tapers = [4 7];
params.err = [2 0.05];

niqf = 2500;
lcf = 20/niqf;
hcf = 200/niqf;
[b,a] = butter(2,[lcf hcf]);

thresh_amp = [-0.59 -0.52 -0.58 -0.6 -0.55 -0.58 -0.53 -0.54 -0.54 -0.62 -0.54 -0.58 -0.54 -0.5 -0.54 -0.52 -0.54];

for d = 1:length(f_loc)
 
    load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = []; 
    
    tim_name = [f_loc{d} '_sampletimes'];
    tim_name(1:24) = [];
    
    eval(['data = ' dat_name ';']);
    eval(['clear ' dat_name])
    eval(['times = ' tim_name ';']);
    eval(['clear ' tim_name])
    
sweep_ends = find(diff(times) < 0);
    
sweep_begs = [1;(sweep_ends)+1];
sweep_ends = [sweep_ends;length(times)];
% sMarkers = [sweep_begs sweep_ends];
% close all
overall_filt = [];
overall_down = [];
overall_up = [];
for i = 1:length(sweep_begs)
    
beg_ind = sweep_begs(i);
end_ind = sweep_ends(i);
cur_seg = data(beg_ind:end_ind);
up_times = find(cur_seg > thresh_amp(d));
down_times = find(cur_seg < thresh_amp(d));
% filt_seg1 = wavefilter(data(beg_ind:end_ind)',5);
filt_seg2 = filtfilt(b,a,data(beg_ind:end_ind));
overall_filt = [overall_filt;filt_seg2];
overall_down = [overall_down;filt_seg2(down_times)];
overall_up = [overall_up;filt_seg2(up_times)];
% filt_seg2 = sqrt(jmm_smooth_1d(filt_seg2.^2,80));
% filt_seg2(filt_seg2>0.05) = 0.05;
plot(times(beg_ind:end_ind),data(beg_ind:end_ind))
% hold on
% plot(times(beg_ind:end_ind),filt_seg1-0.6,'r')
% plot(times(beg_ind:end_ind),filt_seg2,'r')

ylim([-0.8 0.2])
pause
clf

end

% [ S{d}, f{d},Serr{d}]= mtspectrumc_unequal_length_trials( data, movingwin, params, sMarkers );
% 
% t_names = ['C:\WC_Germany\Persistent_activity\heka_spectrum\' dat_name];
% plot(f{d},10*log10(S{d}),'linewidth',2)
% hold on
% plot(f{d},10*log10(Serr{d}(1,:)),'--')
% plot(f{d},10*log10(Serr{d}(2,:)),'--')
% print('-dpng',t_names)
% close

%% find all spikes

% figure
% [y,x] = gpkde(overall_filt,-1);
% plot(x,y)
% hold on
% [y,x] = gpkde(overall_up,-1);
% plot(x,y,'r')
% [y,x] = gpkde(overall_down,-1);
% plot(x,y,'g')
% xlim([-0.08 0.08])
% legend('overall','up','down')
% t_names = ['C:\WC_Germany\Persistent_activity\heka_spectrum\hf_dist_20_200hz_' dat_name];
% print('-dpng',t_names)
% close

end





