clear all
close all

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/pspec_simulation/';

type = 'best'; %avg worst or best, determines which velocity profile to use

%%
dname = '~/Analysis/bruce/FINsac_mod/orth_eyetrajectories';
load(dname);

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

%% load in phosphor response trace, and sample enough repeats to cover the simulated trial
load ~/Data/bruce/misc/ViewsonicFrame.mat
phosphor_Fs = 4e4;
phosphor_lum = -Average3_ViewSonic_FrameRise__Ch1.values;
phosphor_t = (1:length(phosphor_lum))'/phosphor_Fs;
use_seg = find(phosphor_t <= 0.01);
phosphor_t = phosphor_t(use_seg); phosphor_lum = phosphor_lum(use_seg);

phosphor_lum = phosphor_lum/max(phosphor_lum); %scale to have max value of 1

n_traces_per_trial = 5;
phosphor_lum = repmat(phosphor_lum,n_traces_per_trial,1);
phosphor_t = (1:length(phosphor_lum))'/phosphor_Fs;

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
f1 = figure();
subplot(2,1,1)
plot(Trial_Tax*1e3,interp_phosphor)
xlabel('Time (ms)');
ylabel('Relative phosphor intensity');

subplot(2,1,2)
plot(Trial_Tax*1e3,interp_velprof)
xlabel('Time (ms)');
ylabel('Orthogonal eye speed (deg/sec)');

% %PRINT PLOTS
% fig_width = 6; rel_height = 1.2;
% figufy(f1);
% fname = [fig_dir 'eyespeed_phosphor_profiles_' type '.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%

n_rpts = 100; %number of times to repeat simulation for averaging
frame_dur = 0.01; %duration of each stimulus frame
Ntrials = 50;%number of simulated trials to pool for calculating FFT at higher temp freq sampling
Nframes = ceil(Trial_Dur/frame_dur); %number of stimulus frames per trial
wi = 2; %width of stimulus (deg)
bar_width = 0.0565; %width of individual bars (deg)
target_pixsize = 0.5/trial_Fs; %desired spatial sampling
spatial_usfac = round(bar_width/target_pixsize); %required spatial up-sampling factor
nbars = round(wi/bar_width); %number of bars per frame
npix = nbars*spatial_usfac; %number of pixels per frame
pix_size = nbars*bar_width/npix; %pixel size (deg)

dd = 12; %bar density

H = fspecial('gaussian',5,1); %circular gaussian window for slight frequency smoothing
%%
for nn = 1:n_rpts
    fprintf('Sim rpt %d of %d\n',nn,n_rpts);
    
    %create bar white noise stim
    is_zero = rand(Nframes,nbars,Ntrials) > dd/100;
    stim = randi(2,Nframes,nbars,Ntrials);
    stim(stim==2) = -1;
    stim(is_zero) = 0;
    
    %spatial up-sampling
    stim_up = zeros(Nframes,npix,Ntrials);
    for ii = 1:nbars
        for jj = 1:spatial_usfac
            stim_up(:,spatial_usfac*(ii-1)+jj,:) = stim(:,ii,:);
        end
    end
    
    %% temporal up-sampling
    temp_usfac = round(frame_dur*trial_Fs);
    stim_up = repmat(stim_up,[1 1 1 temp_usfac]);
    stim_up = permute(stim_up,[4 1 2 3]);
    % stim_up = reshape(stim_up,[],npix,Ntrials);
    %% simulate frame-by-frame translation, egiven eye velocity
    vel_in_pix = round(interp_velprof/trial_Fs/pix_size);
    stim_up_trans = stim_up;
    for jj = 1:Nframes
        for ii = 2:temp_usfac
            stim_up_trans(ii,jj,:,:) = shift_matrix_Nd(stim_up_trans(ii-1,jj,:,:),vel_in_pix((jj-1)*temp_usfac+ii-1),3);
        end
    end
    
    %these pixels could have zeros introduced from the image translation
    buffer_win = round(0.4/pix_size);
    stim_up_trans(:,:,1:buffer_win,:) = [];
    stim_up(:,:,1:buffer_win,:) = [];
    npix = size(stim_up,3);
    
    stim_up_trans = reshape(stim_up_trans,[],npix,Ntrials);
    stim_up = reshape(stim_up,[],npix,Ntrials);
    
    %factor in phosphor response
    stim_up_trans = bsxfun(@times,stim_up_trans,interp_phosphor');
    stim_up = bsxfun(@times,stim_up,interp_phosphor');
    
    stim_up_trans = permute(stim_up_trans,[1 3 2]);
    stim_up = permute(stim_up,[1 3 2]);
    stim_up_trans = reshape(stim_up_trans,[],npix);
    stim_up = reshape(stim_up,[],npix);
    %% compute and store ffts
    NT = length(Trial_Tax)*Ntrials;
    niqf_x = 1/(2*pix_size);
    niqf_t = trial_Fs/2;
    fx = linspace(-niqf_x,niqf_x,npix);
    ft = linspace(-niqf_t,niqf_t,NT);
    [FX,FT] = meshgrid(fx,ft);
    
    dft = 2*niqf_t/NT;
    dfx = 2*niqf_x/npix;
    PP = abs(fftshift(fft2(stim_up)));
    PP_trans = abs(fftshift(fft2(stim_up_trans)));
    
    fxu = find(abs(fx) <= 200);
    ftu = find(abs(ft) <= 200);
    
    PP = filter2(H,PP);
    PP_trans = filter2(H,PP_trans);
    all_PP(nn,:,:) = PP(ftu,fxu);
    all_PP_eye(nn,:,:) = PP_trans(ftu,fxu);
    
end

%%
sname = [fig_dir 'pspec_calcs_' type];
save(sname,'fx','ft','fxu','ftu','all_PP','all_PP_eye');

%%
avg_PP = squeeze(nanmean(all_PP));
avg_PP_eye = squeeze(nanmean(all_PP_eye));
xr = [-15 15];
tr = [-30 30];

f1 = figure();
subplot(2,1,1);
imagesc(fx(fxu),ft(ftu),avg_PP);
xlim(xr);ylim(tr);
xlabel('Spatial frequency (cyc/deg)');
ylabel('Temporal frequency (Hz)');
set(gca,'ydir','normal');

subplot(2,1,2);
imagesc(fx(fxu),ft(ftu),avg_PP_eye);
xlim(xr);ylim(tr);
xlabel('Spatial frequency (cyc/deg)');
ylabel('Temporal frequency (Hz)');
set(gca,'ydir','normal');

xr = [0 10];
tr = [0 15];
diff_mat = (avg_PP-avg_PP_eye)./avg_PP;
diff_mat_avg = 0.5*diff_mat + 0.5*flipud(diff_mat);
f2 = figure();
imagesc(fx(fxu),ft(ftu),diff_mat_avg);
% imagesc(fx(fxu),ft(ftu),diff_mat);
xlim(xr);ylim(tr);
switch type
    case 'avg'
        caxis([0 0.3]);
    case 'worst'
        caxis([0 0.5]);
end
colorbar
xlabel('Spatial frequency (cyc/deg)');
ylabel('Temporal frequency (Hz)');
set(gca,'ydir','normal');

% %PRINT PLOTS
% fig_width = 6; rel_height = 1.4;
% figufy(f1);
% fname = [fig_dir 'amplitude_spectra_' type '.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 5; rel_height = 0.7;
% figufy(f2);
% fname = [fig_dir 'spectra_diffs_' type '.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
