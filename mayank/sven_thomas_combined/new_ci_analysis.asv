clear all
close all
%% FOR LEC CELLS 
depol_array{1} = 'C:\wc_data\2009-04-21_CWC_LFP\2009-4-21_12s_CI.mat'; %18
depol_array{2} = 'C:\wc_data\2010-05-29_CWC_LFP\2010-05-29_12s_CI.mat'; %19
depol_array{3} = 'C:\wc_data\2010-08-12_CWC_LFP_B\2010-08-12_B_12s_CI.mat'; %20
depol_array{4} = 'C:\wc_data\2012_06_21_A\2012_06_21_1_12s_CI.mat';
depol_array{5} = 'C:\wc_data\2012_06_21_B\2012_06_21_2_12s_CI.mat';
depol_array{6} = 'C:\wc_data\2012_06_23_C\2012_06_23_4_12s_CI.mat';
depol_array{7} = 'C:\wc_data\2012_06_24_A\2012_06_24_1_12s_CI.mat';
depol_array{8} = 'C:\wc_data\2012_06_24_B\2012_06_24_2_12s_CI.mat';
depol_array{9} = 'C:\wc_data\2012_06_24_C\2012_06_24_3_12s_CI.mat';

xi = linspace(-80,-40,100);
bndwdth = 2;
Fs = 2e4;
for d = 1:9
    d
    load(depol_array{d})
    time = (1:length(data))/Fs;
    [mean_sweep(d,:),sweep_mat,t_axis,stim_rate(d),non_stim_rate(d)] = compute_ci_avgs(data(:),time);
    n_stims(d) = size(sweep_mat,1);
    down_t = downsample(t_axis,50);
    down_sweep = downsample(sweep_mat',50)';
    dens_mat(d,:,:) = zeros(length(down_t),100);
    dens_mat_sm(d,:,:) = zeros(length(down_t),100);
    for i = 1:length(down_t)
       dens_mat(d,i,:) = ksdensity(down_sweep(:,i),xi,'width',bndwdth);
    end
    for i = 1:length(xi)
        dens_mat_sm(d,:,i) = jmm_smooth_1d_cor(dens_mat(d,:,i),3);
    end
    % figure
    %     imagesc(sweep_mat)
%     pause
%     close all
clear data
end

% figure
% pcolor(down_t,xi,squeeze(mean(dens_mat))');shading flat
figure
pcolor(down_t,xi,squeeze(mean(dens_mat_sm))');shading flat
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (mV)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

%% NOW COMPUTE BACKGRND SPIKE RATES FOR THESE LEC NEURONS
backgnd_dir{1} = 'C:\wc_data\2009-04-21_CWC_LFP\2009-4-21_22-20-11\';
backgnd_dir{2} = 'C:\wc_data\2010-05-29_CWC_LFP\';
backgnd_dir{3} = 'C:\wc_data\2010-08-12_CWC_LFP_B\2010-8-12_14-55-7\';
backgnd_dir{4} = 'C:\wc_data\2012_06_21_A\2012-6-21_14-33-23\';
backgnd_dir{5} = 'C:\wc_data\2012_06_21_B\2012-6-21_16-1-44\';
backgnd_dir{6} = 'C:\wc_data\2012_06_23_C\2012-6-23_19-39-37\';
backgnd_dir{7} = 'C:\wc_data\2012_06_24_A\2012-6-24_14-49-46\';
backgnd_dir{8} = 'C:\wc_data\2012_06_24_B\2012-6-24_15-46-53\';
backgnd_dir{9} = 'C:\wc_data\2012_06_24_C\2012-6-24_16-24-8\';


clear backgnd_rate 
for d = 1:9
    cd(backgnd_dir{d})
    load ./spike_time_jmm
    load ./sync_times.mat
    if d==4 %only use data after 950 s
        spk_times = round(spkid)/2016;
        t = (1:length(synct))/2016;
        backgnd_rate(d) = length(find(spk_times > 950))/sum(t > 950);        
    else
    backgnd_rate(d) = length(spkid)/length(synct)*2016;  
    end
end
%%
jp = 7;
figure
% shadedErrorBar(t_axis-8,mean(mean_sweep)-jp,std(mean_sweep)/sqrt(9))
shadedErrorBar(t_axis,mean(mean_sweep)-jp,std(mean_sweep)/sqrt(9))
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (mV)','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')
xlim([0 30])

%%
cd C:\WC_Germany\sven_thomas_combined\
save lec_12s_ci_data t_axis mean_sweep backgnd_rate stim_rate non_stim_rate

%% load heka data
% cd C:\wc_data\2010-08-12_CWC_LFP_B\2010-8-12_15-22-58
cd C:\wc_data\2010-05-29_CWC_LFP\
cd C:/wc_data/2009-04-21_CWC_LFP/2009-4-21_22-20-11/
dsf = 8;
Fs = 2016;
Fsd = Fs/dsf;
niqf = Fs/2;
load ./used_data wcv lf8
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

lf8_f = filtfilt(b,a,lf8);
lf8_d = downsample(lf8_f,dsf);
lf8_d = zscore(lf8_d);

nlx_t_d = (1:length(lf8_d))/Fsd;
nlx_t = (1:length(wcv))/Fs;

base_nlxoffset = 0;

%% load Heka data
junction_pot = 7;
% load(depol_array{3})
load(depol_array{1})
Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
time = (1:length(data))/Fs;
data = downsample(data,dsf);
time = downsample(time,dsf);


%% align data
approx_nlx_times = find(nlx_t >= base_nlxoffset & nlx_t < base_nlxoffset + max(time));
[dc_offset,dc_maxcorr] = align_dc_ac_sigs_initial_v2(data,time,wcv(approx_nlx_times));
nlx_d = wcv(approx_nlx_times);
nlx_t = nlx_t(approx_nlx_times)- base_nlxoffset - dc_offset;

nlx_t_d = nlx_t_d - base_nlxoffset - dc_offset;

%% plot data
figure
plot(nlx_t_d,lf8_d,'r')
xlim([90 120])
set(gca,'fontsize',16,'fontname','arial')

figure
plot(time,data-junction_pot)
xlim([90 120])
x = 0:.001:30;

s = zeros(size(x));
s(x >= 9 & x <= 21) = 1;
x = x + 90;
hold on
plot(x,s*20-60,'k')
set(gca,'fontsize',16,'fontname','arial')



%% FOR NEW DISTAL MEC RECS
clear all
close all

Fs = 2e4;
dsf = 40;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 3;

stim_dur = 12;
pause_dur = 18;
first_stim = 9;

depol_array{1} = 'C:\wc_data\2012_06_14\2012_06_14_2_CI.mat'; 
depol_array{2} = 'C:\wc_data\2012_06_17\2012_06_17_1_CI.mat'; 

depol_array{3} = 'C:\wc_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI.mat';
depol_array{4} = 'C:\wc_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_12s_CI.mat';
depol_array{5} = 'C:\wc_data\2009-04-13_B\2009-04-13_CWC_LFP_B_12s_CI_smaller.mat';

xi = linspace(-80,-20,100);
bndwdth = 3;
clear mean_sweep12 stim_rate12 non_stim_rate12
for d = 1:length(depol_array)
    d
    load(depol_array{d})
    time = (1:length(data))/Fs;
    load(depol_array{d})
    time = (1:length(data))/Fs;
    [mean_sweep12(d,:),sweep_mat,t_axis,stim_rate12(d),non_stim_rate12(d)] = compute_ci_avgs(data(:),time);
    n_stims(d) = size(sweep_mat,1);
    down_t = downsample(t_axis,50);
    down_sweep = downsample(sweep_mat',50)';
    dens_mat(d,:,:) = zeros(length(down_t),100);
    dens_mat_sm(d,:,:) = zeros(length(down_t),100);
    for i = 1:length(down_t)
        dens_mat(d,i,:) = ksdensity(down_sweep(:,i),xi,'width',bndwdth);
    end
    for i = 1:length(xi)
        dens_mat_sm(d,:,i) = jmm_smooth_1d_cor(dens_mat(d,:,i),3);
    end    
end

% load ./old_12ci_sweep_mats.mat
% for d = 1:2
%     down_sweep = 100*downsample(sweep_mat_old{d}',50)';
%     odens_mat(d,:,:) = zeros(size(down_sweep,2),100);
%     odens_mat_sm(d,:,:) = zeros(size(down_sweep,2),100);
%     for i = 1:size(down_sweep,2)
%         odens_mat(d,i,:) = ksdensity(down_sweep(:,i),xi,'width',bndwdth);
%     end
%     for i = 1:length(xi)
%         odens_mat_sm(d,:,i) = jmm_smooth_1d_cor(odens_mat(d,:,i),2);
%     end
% end

% figure
% pcolor(down_t,xi,squeeze(mean(dens_mat))');shading flat
figure
pcolor(down_t,xi,squeeze(mean(dens_mat_sm))');shading flat
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (mV)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

%%

% backgnd_dir{1} = 'C:\wc_data\2007-05-31_CWC_LFP\2007-5-31_17-59-31\';
% backgnd_dir{2} = 'C:\wc_data\2007-06-04_CWC_LFP_B\2007-6-4_16-24-12\';
backgnd_dir{1} = 'C:\wc_data\2009-04-07\2009-4-7-19\'; %18
backgnd_dir{2} = 'C:\wc_data\2009-04-13_A\2009-4-13-18\'; %19
backgnd_dir{3} = 'C:\wc_data\2009-04-13_B\'; %20
backgnd_dir{4} = 'C:\wc_data\2012_06_14\2012-6-14_13-40-31\';
backgnd_dir{5} = 'C:\wc_data\2012_06_17\2012-6-17_14-34-10';

clear backgnd_rate 
for d = 1:5
    cd(backgnd_dir{d})
    load ./spike_time_jmm
    load ./sync_times.mat
    backgnd_rate(d) = length(spkid)/length(synct)*2016;  
end

%%
jp = 7;
figure
shadedErrorBar(t_axis,mean(mean_sweep12)-jp,std(mean_sweep12)/sqrt(5))
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (mV)','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')
xlim([0 30])
%%
% cd C:\WC_Germany\persistent_9_27_2010
% load ./all_12s_CI_data.mat
% 
% starttime = find(t_axis > 8, 1,'first')+1;
% endtime = find(t_axis > 24, 1, 'first')-1;
% all_mean_sweeps = [100*mean_sweep; new_mean_sweep(:,starttime:endtime); mean_sweep12(:,starttime:endtime)];
% % all_mean_sweeps = [100*mean_sweep; new_mean_sweep(:,starttime:endtime)];
% stim_rates = [stim_rate(:); stim_rate12(:)];
% 
% cd C:\WC_Germany\sven_thomas_combined\
% save all_12s_mec_CI_new starttime endtime t_axis all_mean_sweeps stim_rates backgnd_rate
%%
% jp = 7;
% figure
% shadedErrorBar(t_axis(starttime:endtime)-8,mean(all_mean_sweeps)-jp,std(all_mean_sweeps)/sqrt(7))
% xlabel('Time (s)','fontsize',18)
% ylabel('Amplitude (mV)','fontsize',18)
% set(gca,'fontsize',14,'fontname','arial')
% xlim([0 16])
