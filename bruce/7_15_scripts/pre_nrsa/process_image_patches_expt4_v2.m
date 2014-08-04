clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12/G029/

load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
load ./eye_calibration_data

dt = 118/1e4;
frames_per_jump = 44;
ims_per_trial = 4;
dt_per_jump = frames_per_jump*dt;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [-0.4 1.2; -1.22 0.4]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-1 2; -2. 1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

min_trial_dur = 2;
max_sac_amp = 1;

%%
Expt_nu = [14 22 27]; %these are the grating expts
spk_win = 0.15; %in sec
spk_delay = 0.0;
n_allunits = 96;
min_trial_dur = 1.5+spk_win;
all_xo_vec = [];
all_yo_vec = [];
all_image_vec = [];
all_angles = [];
all_jump_sizes = [];
all_expt_vec = [];
all_trial_vec = [];
all_binned_spks = [];
all_stim_num = [];
all_eyepos = [];
all_t = [];

single_units = find(CellList(Expt_nu(1),:,1) > 0);
multi_units = setdiff(1:n_allunits,single_units);
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    n_sus = length(single_units);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     completed_trials = 1:length(Trial_durs); %use all trials
    completed_trials = find(Trial_durs > min_trial_dur);
    Trial_im_nums = [Expts{Expt_nu(ee)}.Trials(:).se];
    Trial_angles = [Expts{Expt_nu(ee)}.Trials(:).Fa];
    jump_sizes = [Expts{Expt_nu(ee)}.Trials(:).sz];
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    % %     correct eye positions (quadratic)
    %     clear corrected*
    %     X = [eye_vals(:,1:2) eye_vals(:,1:2).^2 eye_vals(:,1).*eye_vals(:,2)];
    %     corrected_left(:,1) = bsxfun(@plus,X*bx_left(2:end),bx_left(1));
    %     corrected_left(:,2) = bsxfun(@plus,X*by_left(2:end),by_left(1));
    %     X = [eye_vals(:,3:4) eye_vals(:,3:4).^2 eye_vals(:,3).*eye_vals(:,4)];
    %     corrected_right(:,1) = bsxfun(@plus,X*bx_right(2:end),bx_right(1));
    %     corrected_right(:,2) = bsxfun(@plus,X*by_right(2:end),by_right(1));
    
    eye_vals(:,1:2) = corrected_left;
    eye_vals(:,3:4) = corrected_right;
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts,in_blink);
    
    blink_start_times = [blink_data(:).start_times];
    blink_stop_times = [blink_data(:).stop_times];
    sac_amps = [sac_data(:).amplitude];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
    big_sacs = find(sac_amps > max_sac_amp);
    big_sac_starts = sac_start_times(big_sacs);
    big_sac_ends = sac_end_times(big_sacs);
    
    %     check for blinks or big sacs within each trial
    bad_trials = [];
    for i = 1:length(completed_trials)
        if ~isempty(find(big_sac_starts > Trial_starts(completed_trials(i)) & ...
                big_sac_starts < Trial_ends(completed_trials(i)), 1)) ...
                | ~isempty(find(big_sac_ends >= Trial_starts(completed_trials(i)) & ...
                big_sac_ends <= Trial_ends(completed_trials(i)), 1))
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d sac trials\n',length(bad_trials),length(completed_trials));
    completed_trials(bad_trials) = [];
    
    bad_trials = [];
    for i = 1:length(completed_trials)
        if ~isempty(find(blink_start_times >= Trial_starts(completed_trials(i)) & ...
                blink_start_times <= Trial_ends(completed_trials(i)), 1)) ...
                | ~isempty(find(blink_stop_times >= Trial_starts(completed_trials(i)) & ...
                blink_stop_times <= Trial_ends(completed_trials(i)), 1))
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d blink trials\n',length(bad_trials),length(completed_trials));
    completed_trials(bad_trials) = [];
    
    for i = 1:length(completed_trials)
        cur_images = repmat(Trial_im_nums(completed_trials(i)),ims_per_trial,1);
        
        cur_xo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*cos(degtorad(Trial_angles(completed_trials(i))));
        cur_yo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*sin(degtorad(Trial_angles(completed_trials(i))));
        
        cur_t_starts = Trial_starts(completed_trials(i)) + (0:3)*dt_per_jump + spk_delay;
        cur_t_stops = cur_t_starts + spk_win;
        cur_t_edges = sort([cur_t_starts cur_t_stops]);
        cur_t_cents = 0.5*cur_t_starts + 0.5*cur_t_stops;
    
        interp_eyepos = interp1(eye_ts,eye_vals,cur_t_cents);

        cur_binned_spks = nan(n_allunits,length(cur_t_edges));
        for j = 1:n_allunits
            cur_binned_spks(j,:) = histc(Clusters{(j)}.times,cur_t_edges);
        end
        cur_binned_spks = cur_binned_spks(:,1:2:length(cur_t_edges));
        
        all_eyepos = [all_eyepos; interp_eyepos];
        all_t = [all_t; cur_t_cents'];
        all_image_vec = [all_image_vec; cur_images];
        all_xo_vec = [all_xo_vec; cur_xo'];
        all_yo_vec = [all_yo_vec; cur_xo'];
        all_expt_vec = [all_expt_vec; ones(length(cur_images),1)*Expt_nu(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_images),1)*completed_trials(i)];
        all_binned_spks = [all_binned_spks; cur_binned_spks'];
        all_angles = [all_angles; repmat(Trial_angles(completed_trials(i)),ims_per_trial,1)];
        all_jump_sizes = [all_jump_sizes; repmat(jump_sizes(completed_trials(i)),ims_per_trial,1)];
        all_stim_num = [all_stim_num; (1:ims_per_trial)'];
    end
end

%%
cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(all_image_vec));
used_images = unique(all_image_vec);
tot_samps = length(all_image_vec);
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg);
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
    
    for j = 1:length(cur_samp_set)
        ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
        cur_im = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        all_im_patches(cur_samp_set(j),:,:) = cur_im - mean(cur_im(:));
    end
end

%%
fullX = reshape(all_im_patches,size(all_im_patches,1),length(ypatch_inds)*length(xpatch_inds));
fullX = fullX/std(fullX(:));
[NT,klen] = size(fullX);
sdim = sqrt(klen);
full_binned_spks = all_binned_spks;
%%
cd ~/Data/bruce/7_15_12/G029/
load ./expt1_eyecor_d1p25_nosac_v2 gabor*
load ./temp_gabor_params
tot_out_weights = sum(temp_gabor_weights_cor,2);

gabor_params = gabor_params_f{end};
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
for t = 1:96
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_filts(t,:) = gabor_emp1(:);
    gabor_filts2(t,:) = gabor_emp2(:);
    
end

gabor_out_X1 = fullX*gabor_filts';
gabor_out_X2 = fullX*gabor_filts2';
gabor_out_X = gabor_out_X1.^2 + gabor_out_X2.^2;
gabor_out_X = sqrt(gabor_out_X);

for t = 1:96
        
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    [fitp,grad] = GLMsolve_jmm([gabor_out_X(:,t)], spikebins, [tot_out_weights(t); 0], 1, [], [], [], [], [], [1], 0);
    tot_out_offset(t) = fitp.k(end);
    
end

%% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.15;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);
n_xv_cells = length(xv_set);


%% ESTIMATE LL for each shift in each stimulus frame
max_shift = 25;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

frame_LLs = zeros(NT,n_shifts);
Robs = full_binned_spks(:,tr_set);


gabor_filt_bank1 = reshape(gabor_filts(tr_set,:)',[sdim sdim n_tr_cells]);
gabor_filt_bank2 = reshape(gabor_filts2(tr_set,:)',[sdim sdim n_tr_cells]);
shifted_gabor_bank1 = nan(sdim^2,n_tr_cells);
shifted_gabor_bank2 = nan(sdim^2,n_tr_cells);

% frame_xvLLs = zeros(NT,n_shifts);
% Robsxv = full_binned_spks(:,xv_set);
% xgabor_filt_bank1 = reshape(gabor_filts(xv_set,:)',[sdim sdim n_tr_cells]);
% xgabor_filt_bank2 = reshape(gabor_filts2(xv_set,:)',[sdim sdim n_tr_cells]);
% xshifted_gabor_bank1 = nan(sdim^2,n_tr_cells);
% xshifted_gabor_bank2 = nan(sdim^2,n_tr_cells);

shift_cnt = 1;
for xx = 1:length(x_shifts)
    for yy = 1:length(y_shifts)
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift3d(gabor_filt_bank1,x_shifts(xx),2);
        d2 = dist_shift3d(d2,y_shifts(yy),1);
        shifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
        d2 = dist_shift3d(gabor_filt_bank2,x_shifts(xx),2);
        d2 = dist_shift3d(d2,y_shifts(yy),1);
        shifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
        
        gabor_outs1 = fullX*shifted_gabor_bank1;
        gabor_outs2 = fullX*shifted_gabor_bank2;
        energy_out = sqrt(gabor_outs1.^2 + gabor_outs2.^2);
        
        gfun = bsxfun(@times,energy_out,tot_out_weights(tr_set)');
        gfun = bsxfun(@plus,gfun,tot_out_offset(tr_set));
%         gfun = bsxfun(@times,energy_out,gabor_params(tr_set,7)');
%         gfun = bsxfun(@plus,gfun,gabor_params(tr_set,8)');
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs = Robs.*log(pred_rate) - pred_rate;
        frame_LLs(:,shift_cnt) = sum(LLs,2);
      
%         d2 = dist_shift3d(xgabor_filt_bank1,x_shifts(xx),2);
%         d2 = dist_shift3d(d2,y_shifts(yy),1);
%         xshifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
%         d2 = dist_shift3d(xgabor_filt_bank2,x_shifts(xx),2);
%         d2 = dist_shift3d(d2,y_shifts(yy),1);
%         xshifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
%         
%         gabor_outs1 = fullX*xshifted_gabor_bank1;
%         gabor_outs2 = fullX*xshifted_gabor_bank2;
%         energy_out = sqrt(gabor_outs1.^2 + gabor_outs2.^2);
%         
%         gfun = bsxfun(@times,energy_out,tot_out_weights(xv_set)');
%         gfun = bsxfun(@plus,gfun,tot_out_offset(xv_set));
% %         gfun = bsxfun(@times,energy_out,gabor_params(xv_set,7)');
% %         gfun = bsxfun(@plus,gfun,gabor_params(xv_set,8)');
%         too_large = gfun > 50;
%         pred_rate = log(1+exp(gfun));
%         pred_rate(too_large) = gfun(too_large);
%         pred_rate(pred_rate < 1e-20) = 1e-20;
%         
%         LLs = Robsxv.*log(pred_rate) - pred_rate;
%         frame_xvLLs(:,shift_cnt) = sum(LLs,2);

        shift_cnt = shift_cnt + 1;
    end
end


%%
eps_prior_sigma = 0.15; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize

frame_post = bsxfun(@plus,frame_LLs,leps_prior');

[max_post,max_loc] = max(frame_post,[],2);
x_cor = SH(max_loc,1); 
y_cor = SH(max_loc,2);

% xframe_post = bsxfun(@plus,frame_xvLLs,leps_prior');
% [max_post,max_loc] = max(xframe_post,[],2);
% x_corx = SH(max_loc,1);
% y_corx = SH(max_loc,2);

%%
resh_X = reshape(fullX',[sdim sdim NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,NT)';

%%
for t = 1:96
    cur_gabor_params = gabor_params(t,:);
%     cur_gabor_params(7) = tot_out_weights(t);
%     cur_gabor_params(8) = temp_out_offset(t);
    cur_gabor_params(7) = gabor_params(t,7);
    cur_gabor_params(8) = gabor_params(t,8);
    LL(t) = get_gabor_energy_LL(cur_gabor_params,fullX,full_binned_spks(:,t),XX,YY);
    LL_after(t) = get_gabor_energy_LL(cur_gabor_params,fullX_sh,full_binned_spks(:,t),XX,YY);    
end

%%
Fs_t = 1/.5;
lcf = 1/30;
[fb,fa] = butter(2,lcf/(Fs_t/2),'high');

all_eyepos_f = filtfilt(fb,fa,all_eyepos);

%%
close all
figure(1)
set(gcf,'Position',[500 500 800 800])
set(gca,'ydir','normal');
for i = 1:NT
    imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(frame_LLs(i,:),length(y_shifts),length(x_shifts)));
    pause(0.5)
end
