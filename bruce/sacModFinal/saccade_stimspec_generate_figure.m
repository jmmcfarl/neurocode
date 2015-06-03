clear all
close all
type_set = {'worst','avg','best'};
cmap = [1 0 0; 0 0 0; 0 0 1];

% fig_dir = '/Users/james/Analysis/bruce/FINsac_mod/pspec_simulation/';
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/pspec_simulation/';
f1 = figure(); hold on

%%
f2 = figure();

% dname = '~/Analysis/bruce/FINsac_mod/orth_eyetrajectories';
dname = '~/Analysis/bruce/FINsac_mod/orth_eyetrajectories_FIN';
load(dname);

%% load in phosphor response trace, and sample enough repeats to cover the simulated trial
load ~/Data/bruce/misc/ViewsonicFrame.mat
phosphor_Fs = 4e4;
phosphor_lum = -Average3_ViewSonic_FrameRise__Ch1.values;
phosphor_t = (1:length(phosphor_lum))'/phosphor_Fs;
use_seg = find(phosphor_t <= 0.01);
phosphor_t = phosphor_t(use_seg); phosphor_lum = phosphor_lum(use_seg);

phosphor_lum = phosphor_lum/max(phosphor_lum); %scale to have max value of 1

n_traces_per_trial = 2;
phosphor_lum = repmat(phosphor_lum,n_traces_per_trial,1);
phosphor_t = (1:length(phosphor_lum))'/phosphor_Fs;
%%

% H = fspecial('gaussian',10,1); %circular gaussian window for slight frequency smoothing
tfreq_smooth = 5; %smoothing sigma (Hz)
for tt = 1:length(type_set)
    fprintf('Computing for type %s\n',type_set{tt});
    type = type_set{tt};
    
%     sname = [fig_dir 'pspec_calcs_' type];
    sname = [fig_dir 'pspec_calcs_FIN_' type];
    load([sname '.mat']);
    
    %% compute relative modulation of amplitude spectrum
    avg_pow = squeeze(nanmean(all_PP.^2)); %power spectrum during fixation
    avg_pow_eye = squeeze(nanmean(all_PP_eye.^2)); %power spectrum during saccaddes
    
    %     avg_pow = filter2(H,avg_pow);
    %     avg_pow_eye = filter2(H,avg_pow_eye);
    
    %avg together power spectrum for the two eye directions
    avg_pow_eye = 0.5*avg_pow_eye + 0.5*flipud(avg_pow_eye);
    avg_amp_eye = sqrt(avg_pow_eye); %amplitude spectrum during saccades
    avg_amp = sqrt(avg_pow); %amplitude spectrum during fixation
    
    %relative change in amp spectrum during saccades
    rel_amp_diff = (avg_amp_eye - avg_amp)./avg_amp;
    
    %% load the preferred spatial frequencies of all the units used in analysis
    sname = '~/Analysis/bruce/FINsac_mod/used_unit_SFs.mat';
    load(sname);
    
    %% find change in temporal amplitude spectrum at each neuron's preferred SF
    [XX,TT] = meshgrid(fx(fxu),ft(ftu));
    V1 = interp2(XX,TT,rel_amp_diff,used_unit_FSFs,ft(ftu),'bilinear');
    
    %% slightly smooth in temporal frequency
    V1_smoothed = V1;
    tfreq_smoothwin = round(tfreq_smooth/median(diff(ft)));
    for ii = 1:length(used_unit_FSFs)
        V1_smoothed(:,ii) = jmm_smooth_1d_cor(V1_smoothed(:,ii),tfreq_smoothwin);
    end
    %% plot avg (and SEM) across SUs
    figure(f1);
    shadedErrorBar(ft(ftu),mean(V1_smoothed,2),std(V1_smoothed,[],2)/sqrt(length(used_unit_gSFs)),{'color',cmap(tt,:)});
    
    %%
    %pick out chunk of velocity profile to use
    chunk_dur = 0.032;
    beg_offset = 0.001;
    sac_tax = orth_trajects.tax;
    switch type
        case 'avg'
            orth_velprof = orth_trajects.avg_orth_speed;
        case 'best'
            orth_velprof = orth_trajects.avg_acc_speed;
        case 'worst'
            orth_velprof = orth_trajects.avg_inac_speed;
    end
    % orth_velprof = orth_trajects.avg_inac_speed;
    use_chunk = find(sac_tax >= beg_offset & sac_tax <= (beg_offset + chunk_dur));
    sac_tax = sac_tax(use_chunk) - sac_tax(use_chunk(1));
    orth_velprof = orth_velprof(use_chunk);
        
    %% resample both phosphor trace and velocity profile onto sufficiently dense time axis
    trial_Fs = 5e3;
    
    Trial_Tax = 1/trial_Fs:1/trial_Fs:max(sac_tax);
    interp_velprof = interp1(sac_tax,orth_velprof,Trial_Tax);
    interp_phosphor = interp1(phosphor_t,phosphor_lum,Trial_Tax);
    
    Trial_Dur = 0.03;
    ep = find(Trial_Tax >= Trial_Dur,1);
    interp_velprof = interp_velprof(1:ep);
    interp_phosphor = interp_phosphor(1:ep);
    Trial_Tax = Trial_Tax(1:ep);
    
%% plot phosphor response and eye velocity profiles
figure(f2);
subplot(2,1,1)
plot(Trial_Tax*1e3,interp_phosphor)
xlabel('Time (ms)');
ylabel('Relative phosphor intensity');

subplot(2,1,2); hold on
plot(Trial_Tax*1e3,interp_velprof,'color',cmap(tt,:))
xlabel('Time (ms)');
ylabel('Orthogonal eye speed (deg/sec)');

end;


%%
figure(f1);
xlim([-100 100]);
ylim([-0.2 0.2]);

line([-100 100],[0 0],'color','k','linestyle','--');
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Temporal frequency (Hz)');
ylabel('Change in amplitude spectrum (fold-difference)');

fig_width = 6; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'stimulus_temporal_spectra.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
figure(f2);
subplot(2,1,1);
ylim([0 1]);

%PRINT PLOTS
fig_width = 6; rel_height = 1.2;
figufy(f2);
fname = [fig_dir 'eyespeed_phosphor_profiles.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

