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

% cd G034
% load ./jbeG034.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G034Expts.mat
% % load ./eye_calibration_data

dt = 118*2/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.5;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [0 0.8; -0.8 -0.]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

min_trial_dur = 2;
max_sac_amp = 0.5;
%%
% Expt_nu = [1 6 16 17 20 25 28]; %these are the grating expts
Expt_nu = [1 6 16 17 20 25 28]; %these are the grating expts
% Expt_nu = [1 2 17 18 19 20 23 24]; %these are the grating expts
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
%         cur_images = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).Seedseq(1:2:end-1)];
%         cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo(1:2:end-1)];
        cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo(1:end-1)];
        cur_yo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).yo(1:2:end-1)];
                   
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*length(cur_images));
        
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
        all_trial_durs = [all_trial_durs; Trial_durs'];
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

%%

% all_image_vec = all_image_vec(randperm(length(all_image_vec)));

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
        ypatch_inds_adj = round(ypatch_inds + all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds + all_xo_vec(cur_samp_set(j))*Fsd);
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
%         cur_patch = cur_patch - mean(cur_patch(:));
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
    end
end

%%
flen = 7;
trial_changes = find(abs(diff(all_trial_vec)) > 0) + 1;
trial_starts = [1; trial_changes];
trial_stops = [trial_changes-1; length(all_trial_vec)];
fullX = [];
full_binned_spks = [];
for i = 1:length(trial_changes)
    fprintf('Compiling trial %d of %d\n',i,length(trial_changes));
    cur_inds = trial_starts(i):trial_stops(i);
    clen(i) = length(cur_inds);
    fullX = [fullX; makeStimRows(all_im_patches(cur_inds,:,:),flen,1)];
    full_binned_spks = [full_binned_spks; all_binned_spks(cur_inds(flen:end),:)];
end

%%
cd ~/Data/bruce/7_15_12/G029/
% cd ~/Data/bruce/7_15_12/G034/
save Expt1_compiled_data_yinv_xinv full* flen dt Fsd

%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
% 
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
