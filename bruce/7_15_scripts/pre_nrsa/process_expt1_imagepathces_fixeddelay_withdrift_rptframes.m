clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G029/
% load ./jbeG029.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G029Expts.mat
load ./eye_calibration_data

cd ..
cd G034
load ./jbeG034.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G034Expts.mat
% load ./eye_calibration_data

% dt = 118*2/1e4;
dt = 118/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [0 0.8; -0.8 -0.]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-0.5075 1.3075; -1.3 0.5]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
sdim = length(xpatch_inds);

min_trial_dur = 2;
max_sac_amp = 5;
%%
% Expt_nu = [1 6 16 17 20 25 28]; 
Expt_nu = [1 2 17 18 19 20 23 24]; 
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
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    eye_vals(:,1:2) = corrected_left;
    eye_vals(:,3:4) = corrected_right;
    
    avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
    
    
    blink_start_times = [blink_data(:).start_times];
    sac_amps = [sac_data(:).amplitude];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
    big_sacs = find(sac_amps > max_sac_amp);
    big_sac_starts = sac_start_times(big_sacs);
    big_sac_ends = sac_end_times(big_sacs);
    
    %check for blinks or big sacs within each trial
    bad_trials = [];
    for i = 1:length(used_trials)
        if ~isempty(find(big_sac_starts > Trial_starts(used_trials(i)) & ...
                big_sac_starts < Trial_ends(used_trials(i)), 1));
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d sac trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    bad_trials = [];
    for i = 1:length(used_trials)
        if ~isempty(find(blink_start_times > Trial_starts(used_trials(i)) & ...
                blink_start_times < Trial_ends(used_trials(i)), 1));
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d blink trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];    
   
    all_rpt_frames = [];
    for i = 1:length(used_trials)
        cur_images = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).Seedseq];
        cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo];
        cur_yo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).yo];

        cur_xo = cur_xo(floor(1:0.5:86)); %upsample x translations
        
        %add in repeat frames
        cur_repeat_frames = Expts{Expt_nu(ee)}.Trials(used_trials(i)).rptframes;
        cur_repeat_frames = flipud(cur_repeat_frames);
        for cc = 1:length(cur_repeat_frames)
            cur_images = [cur_images(1:cur_repeat_frames(cc)); cur_images(cur_repeat_frames(cc));...
                cur_images(cur_repeat_frames(cc)+1:end)];
            cur_xo = [cur_xo(1:cur_repeat_frames(cc)) cur_xo(cur_repeat_frames(cc)) ...
                cur_xo(cur_repeat_frames(cc)+1:end)];
            cur_yo = [cur_xo(1:cur_repeat_frames(cc)) cur_yo(cur_repeat_frames(cc)) ...
                cur_yo(cur_repeat_frames(cc)+1:end)];
        end
        all_rpt_frames = [all_rpt_frames; cur_repeat_frames];
        max_rpt_frames(ee) = max(all_rpt_frames);
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*length(cur_images));
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        interp_eyepos = interp1(eye_ts,eye_vals,cur_t_cents);
        interp_insac = round(interp1(eye_ts,in_sac,cur_t_cents));
        
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
all_t(bad_trials) = [];

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
NT = length(all_image_vec);
min_fix_dur = 0.05;
trial_flips = 1 + find(diff(all_trial_vec) ~= 0);
sac_inds = find(all_insac==1);
all_break_inds = unique([trial_flips; sac_inds]);
fix_start_inds = [1; all_break_inds+1];
fix_stop_inds = [all_break_inds-1; NT];
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
used_fixs = find(fix_durs > min_fix_dur);
fix_start_inds = fix_start_inds(used_fixs);
fix_stop_inds = fix_stop_inds(used_fixs);
fix_durs = fix_durs(used_fixs);

% avg_x_pos = 0.5*all_eyepos(:,1) + 0.5*all_eyepos(:,3);
% avg_y_pos = 0.5*all_eyepos(:,2) + 0.5*all_eyepos(:,4);

%use only data from left eye
avg_x_pos = all_eyepos(:,1);
avg_y_pos = all_eyepos(:,2);

x_drifts = nan(size(avg_x_pos));
y_drifts = nan(size(avg_y_pos));
for i = 1:length(used_fixs)
    cur_set = fix_start_inds(i):fix_stop_inds(i);
    x_drifts(cur_set) = avg_x_pos(cur_set) - median(avg_x_pos(cur_set));
    y_drifts(cur_set) = avg_y_pos(cur_set) - median(avg_y_pos(cur_set));
end

%%

% all_xo_vec = circshift(all_xo_vec,10); %randomize vector of x shifts
% all_yo_vec = circshift(all_yo_vec,10); %randomize vector of x shifts

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
    
    %     cur_eye_x = 0.5*all_eyepos(cur_samp_set,
    
    for j = 1:length(cur_samp_set)
        if ~isnan(x_drifts(cur_samp_set(j))) && ~isnan(y_drifts(cur_samp_set(j)))
            ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd + y_drifts(cur_samp_set(j))*Fsd);
            xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd + x_drifts(cur_samp_set(j))*Fsd);
%             ypatch_inds_adj = round(ypatch_inds + all_yo_vec(cur_samp_set(j))*Fsd + y_drifts(cur_samp_set(j))*Fsd);
%             xpatch_inds_adj = round(xpatch_inds + all_xo_vec(cur_samp_set(j))*Fsd + x_drifts(cur_samp_set(j))*Fsd);
            cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
            %         cur_patch = cur_patch - mean(cur_patch(:));
            all_im_patches(cur_samp_set(j),:,:) = cur_patch;
        end
    end
end

%%
% flen = 7;
fixed_delay = round(0.05/dt);
trial_changes = find(abs(diff(all_trial_vec)) > 0) + 1;
trial_starts = [1; trial_changes];
trial_stops = [trial_changes-1; length(all_trial_vec)];
fullX = [];
full_binned_spks = [];
used_inds = [];
for i = 1:length(trial_changes)
    fprintf('Compiling trial %d of %d\n',i,length(trial_changes));
    cur_inds = trial_starts(i):trial_stops(i);
    cur_inds = cur_inds((1+fixed_delay):end);
    clen(i) = length(cur_inds);
    fullX = [fullX; reshape(all_im_patches(cur_inds-fixed_delay,:,:),length(cur_inds),sdim^2)];
    full_binned_spks = [full_binned_spks; all_binned_spks(cur_inds,:)];
    used_inds = [used_inds; cur_inds(:) - fixed_delay];
end

%all these are in STIMULUS time not in delayed RESPONSE TIME
full_image_vec = all_image_vec(used_inds);
full_xo_vec = all_xo_vec(used_inds);
full_yo_vec = all_yo_vec(used_inds);
full_expt_vec = all_expt_vec(used_inds);
full_trial_vec = all_trial_vec(used_inds);
full_trail_durs = all_trial_durs(used_inds);
full_eyepos = all_eyepos(used_inds,:);
full_insac = all_insac(used_inds);
full_t = all_t(used_inds);
full_xdrift = x_drifts(used_inds);
full_ydrift = y_drifts(used_inds);
%%
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/
% save Expt1_compiled_data_fixedlag_driftcor fullX full* fixed_delay dt Fsd RF_patch* xax yax *_inds dsfrac
save('Expt1_compiled_data_fixedlag_driftcor_leftonly_rptframes.mat','-v7.3','full*','fixed_delay','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
