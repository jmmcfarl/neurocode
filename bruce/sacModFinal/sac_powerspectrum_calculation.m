clear all
close all

%%
% dname = '~/Analysis/bruce/G086/FINsac_mod/gen_trig_avg_data_ori0.mat';
dname = '~/Analysis/bruce/FINsac_mod/orth_eyetrajectories';
load(dname);

% sac_tax = gen_data.raw_eye_lags;
% sac_velprof = gen_data.gsac_rawtavg_eyespeed;
% 
% chunk_dur = 0.032;
% beg_offset = 0.001;
% use_chunk = find(sac_tax >= beg_offset & sac_tax <= (beg_offset + chunk_dur));
% sac_tax = sac_tax(use_chunk) - sac_tax(use_chunk(1));
% 
% avg_delta_Y_frac = mean(abs(gen_data.gsac_delta_Y_frac));
% upper_delta_Y_frac = prctile(abs(gen_data.gsac_delta_Y_frac),75);
% % orth_velprof = sac_velprof(use_chunk)*avg_delta_Y_frac;
% orth_velprof = sac_velprof(use_chunk)*upper_delta_Y_frac;

chunk_dur = 0.032;
beg_offset = 0.001;
sac_tax = orth_trajects.tax;
orth_velprof = orth_trajects.avg_orth_speed;
% orth_velprof = orth_trajects.avg_inac_speed;
use_chunk = find(sac_tax >= beg_offset & sac_tax <= (beg_offset + chunk_dur));
sac_tax = sac_tax(use_chunk) - sac_tax(use_chunk(1));
orth_velprof = orth_velprof(use_chunk);

%%
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

%%
trial_Fs = 5e3;

Trial_Tax = 1/trial_Fs:1/trial_Fs:max(sac_tax);
interp_velprof = interp1(sac_tax,orth_velprof,Trial_Tax);
interp_phosphor = interp1(phosphor_t,phosphor_lum,Trial_Tax);

Trial_Dur = 0.03;
ep = find(Trial_Tax >= Trial_Dur,1);
interp_velprof = interp_velprof(1:ep);
interp_phosphor = interp_phosphor(1:ep);
Trial_Tax = Trial_Tax(1:ep);
%%

n_rpts = 10;
frame_dur = 0.01;
Ntrials = 25;
Nframes = ceil(Trial_Dur/frame_dur);
wi = 2;
bar_width = 0.0565;
target_pixsize = 0.5/trial_Fs;
spatial_usfac = round(bar_width/target_pixsize);
nbars = round(wi/bar_width);
npix = nbars*spatial_usfac;
pix_size = nbars*bar_width/npix;

dd = 12;

H = fspecial('gaussian',5,1);
%%
for nn = 1:n_rpts
    fprintf('Sim rpt %d of %d\n',nn,n_rpts);
    
    is_zero = rand(Nframes,nbars,Ntrials) > dd/100;
    stim = randi(2,Nframes,nbars,Ntrials);
    stim(stim==2) = -1;
    stim(is_zero) = 0;
    
    stim_up = zeros(Nframes,npix,Ntrials);
    for ii = 1:nbars
        for jj = 1:spatial_usfac
            stim_up(:,spatial_usfac*(ii-1)+jj,:) = stim(:,ii,:);
        end
    end
    
    %%
    temp_usfac = round(frame_dur*trial_Fs);
    stim_up = repmat(stim_up,[1 1 1 temp_usfac]);
    stim_up = permute(stim_up,[4 1 2 3]);
    % stim_up = reshape(stim_up,[],npix,Ntrials);
    %%
    vel_in_pix = round(interp_velprof/trial_Fs/pix_size);
    stim_up_trans = stim_up;
    for jj = 1:Nframes
        for ii = 2:temp_usfac
            stim_up_trans(ii,jj,:,:) = shift_matrix_Nd(stim_up_trans(ii-1,jj,:,:),vel_in_pix((jj-1)*temp_usfac+ii-1),3);
        end
    end
    
    stim_up_trans = reshape(stim_up_trans,[],npix,Ntrials);
    stim_up = reshape(stim_up,[],npix,Ntrials);
    
    stim_up_trans = bsxfun(@times,stim_up_trans,interp_phosphor');
    stim_up = bsxfun(@times,stim_up,interp_phosphor');
    
    stim_up_trans = permute(stim_up_trans,[1 3 2]);
    stim_up = permute(stim_up,[1 3 2]);
    stim_up_trans = reshape(stim_up_trans,[],npix);
    stim_up = reshape(stim_up,[],npix);
    %%
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
    
    fxu = find(abs(fx) <= 50);
    ftu = find(abs(ft) <= 100);
    
    PP = filter2(H,PP);
    PP_trans = filter2(H,PP_trans);
    all_PP(nn,:,:) = PP(ftu,fxu);
    all_PP_eye(nn,:,:) = PP_trans(ftu,fxu);
    
end

%%
avg_PP = squeeze(nanmean(all_PP));
avg_PP_eye = squeeze(nanmean(all_PP_eye));
xr = [-15 15];
tr = [-30 30];

f1 = figure();
subplot(2,1,1);
imagesc(fx(fxu),ft(ftu),avg_PP);
xlim(xr);ylim(tr);
subplot(2,1,2);
imagesc(fx(fxu),ft(ftu),avg_PP_eye);
xlim(xr);ylim(tr);

diff_mat = (avg_PP-avg_PP_eye)./avg_PP;
diff_mat_avg = 0.5*diff_mat + 0.5*flipud(diff_mat);
f2 = figure();
imagesc(fx(fxu),ft(ftu),diff_mat_avg);
% imagesc(fx(fxu),ft(ftu),diff_mat);
caxis([0 0.5]);
xlim(xr);ylim(tr);
colorbar

