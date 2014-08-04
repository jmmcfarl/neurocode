clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12
obj_info_dir = '~/James_scripts/data_processing/Images/object_im_info';

% cd G029/
% load ./CellList.mat
% load ./G029Expts.mat
cd G034/
load ./CellList.mat
load ./G034Expts.mat

samp_fac = 0.5;
dt = samp_fac*118/1e4;
fst = 1/dt;
frames_per_jump = 44/samp_fac;

Pix2Deg = 0.018837;
dsfrac = 1.5;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [-0.12 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-0.7 1.5; -1.51 0.7]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [-1.5 2.2; -2.2 1.5]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
sdim = length(xpatch_inds);

rf_cent = [0.34 -0.43];
obj_set = 1142:1369;

min_trial_dur = 1;
max_sac_amp = 0.5;

% load ./all_eyedata_expt3
load ./all_eyedata_expt3_34

% left_sac_amps = sqrt(sum(all_lsac_camps.^2,2));
% right_sac_amps = sqrt(sum(all_rsac_camps.^2,2));
% big_sac_inds = find(left_sac_amps > max_sac_amp | right_sac_amps > max_sac_amp);
big_sac_inds = find(all_sac_amps > max_sac_amp | all_sac_amps > max_sac_amp);
big_sac_start_times = all_t(all_sac_start_inds(big_sac_inds));
big_sac_stop_times = all_t(all_sac_stop_inds(big_sac_inds,1));

blink_start_times = all_t(all_blink_start_inds);
blink_stop_times = all_t(all_blink_stop_inds);

% hf_start_times = all_t(all_hf_start_inds);
% hf_stop_times = all_t(all_hf_stop_inds);

full_image_vec = [];
full_xo_vec = [];
full_yo_vec = [];
full_angles = [];
full_jump_sizes = [];
full_expt_vec = [];
full_trial_vec = [];
full_binned_spks = [];
full_trial_durs = [];
full_t = [];

%%
% Expt_nu = [4 5 9 23 32]; %expt 3
Expt_nu = [7:12]; %expt 3 34
n_allunits = 96;

for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
    n_sus = length(single_units{ee});
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    
    effective_trial_ends = Trial_ends;
    for i = 1:length(Trial_ends)
        cur_blink_times = min(blink_start_times(blink_start_times > Trial_starts(i) & blink_start_times < Trial_ends(i)));
        if isempty(cur_blink_times); cur_blink_times = Inf; end
        cur_big_sac_times = min(big_sac_start_times(big_sac_start_times > Trial_starts(i) & big_sac_start_times < Trial_ends(i)));
        if isempty(cur_big_sac_times); cur_big_sac_times = Inf; end
%         cur_hf_times = min(hf_start_times(hf_start_times > Trial_starts(i) & hf_start_times < Trial_ends(i)));
%         if isempty(cur_hf_times); cur_hf_times = Inf; end
        %        effective_trial_ends(i) = min([cur_hf_times cur_big_sac_times cur_blink_times]);
        effective_trial_ends(i) = min([cur_big_sac_times cur_blink_times]);
    end
    effective_trial_ends(isinf(effective_trial_ends)) = Trial_ends(isinf(effective_trial_ends));
    
    Trial_ends = effective_trial_ends;
    
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    used_trials = find(Trial_durs > min_trial_dur);
    Trial_im_nums = floor([Expts{Expt_nu(ee)}.Trials(:).se]);
    Trial_angles = [Expts{Expt_nu(ee)}.Trials(:).Fa];
    jump_sizes = [Expts{Expt_nu(ee)}.Trials(:).sz];
    
    
    
    tt = nan(length(used_trials),1);
    for i = 1:length(used_trials)
        cur_n_tbins = floor(Trial_durs(used_trials(i))/dt);
        cur_n_shifts = ceil(cur_n_tbins/frames_per_jump);
        tt(i) = cur_n_shifts;
        
        shift_inds = ceil(1/frames_per_jump:1/frames_per_jump:cur_n_tbins/frames_per_jump);
        
        cur_images = repmat(Trial_im_nums(used_trials(i)),cur_n_tbins,1);
        
        cur_xo_vals = (0:(cur_n_shifts-1))*jump_sizes(used_trials(i))*cos(degtorad(Trial_angles(used_trials(i))));
        cur_yo_vals = (0:(cur_n_shifts-1))*jump_sizes(used_trials(i))*sin(degtorad(Trial_angles(used_trials(i))));
        cur_xo = cur_xo_vals(shift_inds);
        cur_yo = cur_yo_vals(shift_inds);
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*cur_n_tbins);
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_binned_spks = nan(n_mus,length(cur_images));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        full_image_vec = [full_image_vec; cur_images];
        full_xo_vec = [full_xo_vec; cur_xo'];
        full_yo_vec = [full_yo_vec; cur_yo'];
        full_expt_vec = [full_expt_vec; ones(length(cur_images),1)*Expt_nu(ee)];
        full_trial_vec = [full_trial_vec; ones(length(cur_images),1)*used_trials(i)];
        full_binned_spks = [full_binned_spks; cur_binned_spks'];
        full_trial_durs = [full_trial_durs; ones(length(cur_images),1)*Trial_durs(used_trials(i))];
        full_t = [full_t cur_t_cents];
    end
end


%%
cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(full_image_vec));
used_images = unique(full_image_vec);
tot_samps = length(full_image_vec);

stim_comp = [full_image_vec full_xo_vec full_yo_vec];
un_stim_vals = unique(stim_comp,'rows');
n_un_stims = size(un_stim_vals,1);

trans_comp = [full_xo_vec full_yo_vec];

all_im_patches = nan(n_un_stims,length(ypatch_inds),length(xpatch_inds));
all_obj_info = nan(n_un_stims,length(ypatch_inds),length(xpatch_inds));
full_stim_ids = nan(size(full_t));

stim_cnt = 1;
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
        
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    if ismember(used_images(i),obj_set)
        cd(obj_info_dir)
        obj_info_id = find(obj_set == used_images(i));
        fname = sprintf('info4%.4d.mat',obj_info_id);
        load(fname);
        
        obj_IMAGE = ov_L;
        obj_IMAGE = flipud(obj_IMAGE);
        obj_IMAGE = round(imresize(obj_IMAGE,1/dsfrac));        
        
    end
    
    cur_samp_set = find(full_image_vec == used_images(i));
    cur_un_samps = unique(trans_comp(cur_samp_set,:),'rows');
    cur_n_samps = size(cur_un_samps,1);
   
    fprintf('Analyzing image %d, %d unique samps\n',i,cur_n_samps);
    for j = 1:cur_n_samps
        
        ypatch_inds_adj = round(ypatch_inds - cur_un_samps(j,2)*Fsd);
        xpatch_inds_adj = round(xpatch_inds - cur_un_samps(j,1)*Fsd);
        
        cur_used_inds = find(full_xo_vec == cur_un_samps(j,1) & full_yo_vec == cur_un_samps(j,2) & ...
            full_image_vec == used_images(i));
        full_stim_ids(cur_used_inds) = stim_cnt;
        
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        all_im_patches(stim_cnt,:,:) = cur_patch;
        
        if ismember(used_images(i),obj_set)
           cur_obj_patch = obj_IMAGE(ypatch_inds_adj,xpatch_inds_adj);
           all_obj_info(stim_cnt,:,:) = cur_obj_patch;
        end
        
        stim_cnt = stim_cnt + 1;
    end
end
stim_cnt = stim_cnt - 1;
resh_all_stims = reshape(all_im_patches,stim_cnt,sdim^2);
resh_all_obj = reshape(all_obj_info,stim_cnt,sdim^2);

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
full_insac = ceil(interp1(all_t,all_insac,full_t));
% full_inhf = ceil(interp1(all_t,all_inhf,full_t));
full_inblink = ceil(interp1(all_t,all_inblink,full_t));
full_insac(1) = 0;
full_eyepos = interp1(all_t,all_eyepos,full_t);

fix_dur_win = 0.15;
NT = length(full_image_vec);
stim_flips = 1+find(diff(full_stim_ids) ~= 0);

trial_start_inds = [1 stim_flips];
trial_stop_inds = [(stim_flips-1) NT];

trial_durs = (trial_stop_inds-trial_start_inds)*dt;
bad_trials = find(trial_durs <= fix_dur_win);
trial_start_inds(bad_trials) = [];
trial_stop_inds(bad_trials) = [];
trial_durs(bad_trials) = [];

trial_stop_winds = trial_start_inds + round(fix_dur_win/dt);


n_trials = length(trial_stop_winds);
trial_eyepos = nan(n_trials,2);
trial_stimnum = nan(n_trials,1);
trial_imnum = nan(n_trials,1);
for i = 1:n_trials
   cur_set = trial_start_inds(i):trial_stop_winds(i);
   trial_eyepos(i,:) = median(full_eyepos(cur_set,1:2));
   uset = ~isnan(full_image_vec(cur_set));
   trial_stimnum(i) = unique(full_stim_ids(cur_set(uset))); 
   trial_imnum(i) = unique(full_image_vec(cur_set(uset))); 
end

%%
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/
save('Expt3_fixbased_data.mat','-v7.3','full_*','trial_*','resh_all*','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
