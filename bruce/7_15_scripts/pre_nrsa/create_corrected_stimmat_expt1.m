clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G029/
load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
load ./eye_calibration_data
% cd ..

dt = 118*2/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
old_dsfrac = 1.25;
old_Fsd = 1/Pix2Deg/old_dsfrac;
dsfrac = 1;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [0 0.8; -0.8 0]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-0.5075 1.3075; -1.3 0.5]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-2 2.7; -2.7 2]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
sdim = length(xpatch_inds);

rf_cent = [0.34 -0.43];

min_trial_dur = 2;
max_sac_amp = 0.5;
%%
Expt_nu = [1 6 16 17 20 25 28];
% Expt_nu = [1 2 17 18 19 20 23 24];
n_allunits = 96;
all_image_vec = [];
all_xo_vec = [];
all_yo_vec = [];
all_expt_vec = [];
all_trial_vec = [];
all_binned_spks = [];
all_trial_durs = [];
all_rep_frames = [];
all_eyepos = [];
all_insac = [];
all_inblink = [];
all_t = [];

for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
    n_sus = length(single_units{ee});
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     used_trials = 1:length(Trial_durs); %use all trials
    used_trials = find(Trial_durs > min_trial_dur);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
%     correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
%     %correct eye positions (quadratic)
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
        
    avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
    
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
    for i = 1:length(used_trials)
        if ~isempty(find(big_sac_starts > Trial_starts(used_trials(i)) & ...
                big_sac_starts < Trial_ends(used_trials(i)), 1)) ...
                 | ~isempty(find(big_sac_ends >= Trial_starts(used_trials(i)) & ...
                big_sac_ends <= Trial_ends(used_trials(i)), 1))
           bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d sac trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    bad_trials = [];
    for i = 1:length(used_trials)
        if ~isempty(find(blink_start_times >= Trial_starts(used_trials(i)) & ...
                blink_start_times <= Trial_ends(used_trials(i)), 1)) ...
                | ~isempty(find(blink_stop_times >= Trial_starts(used_trials(i)) & ...
                blink_stop_times <= Trial_ends(used_trials(i)), 1))
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d blink trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    bad_trials = [];
    for i = 1:length(used_trials)
        if isfield(Expts{Expt_nu(ee)}.Trials(used_trials(i)),'rptframes')
            if ~isempty(Expts{Expt_nu(ee)}.Trials(used_trials(i)).rptframes)
                bad_trials = [bad_trials; i];
            end
        end
    end
    fprintf('Eliminating %d of %d framerep trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    for i = 1:length(used_trials)
        cur_images = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).Seedseq(1:end-1)];
        cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo(1:end-1)];
        cur_yo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).yo(1:2:end-1)];
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*length(cur_images));
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        interp_eyepos = interp1(eye_ts,eye_vals,cur_t_cents);
        interp_insac = round(interp1(eye_ts,in_sac,cur_t_cents));
        interp_inblink = round(interp1(eye_ts,in_blink,cur_t_cents));
        if max(interp_inblink) == 1
            break
        end
        
        cur_binned_spks = nan(n_mus,length(cur_images));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        all_image_vec = [all_image_vec; cur_images];
        all_xo_vec = [all_xo_vec; cur_xo'];
        all_yo_vec = [all_yo_vec; cur_yo'];
        all_expt_vec = [all_expt_vec; ones(length(cur_images),1)*Expt_nu(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_images),1)*used_trials(i)];
        all_binned_spks = [all_binned_spks; cur_binned_spks'];
        %         all_trial_durs = [all_trial_durs; Trial_durs'];
        all_trial_durs = [all_trial_durs; ones(length(cur_images),1)*Trial_durs(used_trials(i))];
        all_eyepos = [all_eyepos; interp_eyepos];
        all_insac = [all_insac; interp_insac'];
        all_inblink = [all_inblink; interp_inblink'];
        all_t = [all_t cur_t_cents];
    end
    
end

%get rid of trials with repeat frames
bad_trials = find(all_rep_frames==1);
fprintf('Eliminating %d of %d frames\n',length(bad_trials),length(all_rep_frames));
all_image_vec(bad_trials) = [];
all_xo_vec(bad_trials) = [];
all_yo_vec(bad_trials) = [];
all_expt_vec(bad_trials) = [];
all_trial_vec(bad_trials) = [];
all_trial_durs(bad_trials) = [];
all_binned_spks(bad_trials,:) = [];
all_eyepos(bad_trials,:) = [];
all_insac(bad_trials) = [];
all_inblink(bad_trials) = [];
all_t(bad_trials) = [];

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
NT = length(all_image_vec);
min_fix_dur = 0.2;
trial_flips = find(diff(all_trial_vec) ~= 0);
used_inds = all_insac==0;
used_inds(trial_flips) = 0;

fix_start_inds = 1+find(used_inds(2:end)==1 & used_inds(1:end-1) == 0);
fix_stop_inds = find(used_inds(2:end) == 0 & used_inds(1:end-1)==1);
fix_start_inds = [1; fix_start_inds];
fix_stop_inds = [fix_stop_inds; NT];

fix_durs = (fix_stop_inds-fix_start_inds)*dt;

used_fixs = find(fix_durs > min_fix_dur);
fix_start_inds = fix_start_inds(used_fixs);
fix_stop_inds = fix_stop_inds(used_fixs);
fix_durs = fix_durs(used_fixs);

%%
cd ~/Data/bruce/7_15_12/G029/
load ./expt1_eyecor_d1p25_nosac_v2 x_cor y_cor
load ./Expt1_newcompiled_data_fixedlag_d1p25_nosac full_t

interp_x_corrs = round(interp1(full_t,x_cor{end},all_t))/old_Fsd;
interp_y_corrs = round(interp1(full_t,y_cor{end},all_t))/old_Fsd;
bad_samps = find(isnan(interp_x_corrs));
fprintf('%d of %d nan samples, ignoring\n',length(bad_samps),length(interp_x_corrs));
interp_x_corrs(bad_samps) = 0;
interp_y_corrs(bad_samps) = 0;
%%
cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(all_image_vec));
used_images = unique(all_image_vec);
tot_samps = length(all_image_vec);
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
all_corim_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
    
    %     cur_eye_x = 0.5*all_eyepos(cur_samp_set,
    
    for j = 1:length(cur_samp_set)
        
        ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
        
        if ~isnan(interp_x_corrs(cur_samp_set(j)))
            ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd + interp_y_corrs(cur_samp_set(j))*Fsd);
            xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd + interp_x_corrs(cur_samp_set(j))*Fsd);
            cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
            
            all_corim_patches(cur_samp_set(j),:,:) = cur_patch;
        end
    end
end

%%
flen = 6;
sdim = length(xpatch_inds);
klen = sdim^2*flen;
used_inds = [];
maxdur = sum(fix_stop_inds-fix_start_inds);

fullX = nan(maxdur,klen);
fullX_cor = nan(maxdur,klen);
full_binned_spks = nan(maxdur,96);

n_fixs = length(fix_start_inds);
used_inds = [];
cnt = 0;
for i = 1:n_fixs
    fprintf('%d of %d\n',i,n_fixs)
    cur_inds = fix_start_inds(i):fix_stop_inds(i);
    temp = makeStimRows(all_im_patches(cur_inds,:,:),flen,1);
    test = all_corim_patches(cur_inds,:,:);
    fullX(cnt + (1:size(temp,1)),:) = temp;
    if max(isnan(test(:))) == 0
        temp2 = makeStimRows(all_corim_patches(cur_inds,:,:),flen,1);
        fullX_cor(cnt + (1:size(temp,1)),:) = temp2;
    end
    full_binned_spks(cnt + (1:size(temp,1)),:) = all_binned_spks(cur_inds(flen:end),:);
    used_inds = [used_inds; cur_inds(flen:end)'];
    cnt = cnt + size(temp,1);
end
fullX(cnt+1:end,:) = [];
fullX_cor(cnt+1:end,:) = [];
full_binned_spks(cnt+1:end,:) = [];


full_image_vec = all_image_vec(used_inds);
full_xo_vec = all_xo_vec(used_inds);
full_yo_vec = all_yo_vec(used_inds);
full_expt_vec = all_expt_vec(used_inds);
full_trial_vec = all_trial_vec(used_inds);
full_t = all_t(used_inds);


%%
cd ~/Data/bruce/7_15_12/G029/
save('corrected_Xmat_flen6_fullres','-v7.3','full*','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
