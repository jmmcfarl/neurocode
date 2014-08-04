clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat
%%
samp_fac = 0.5;
dt = 117.5/1e4;
dt = dt*samp_fac;
fst = 1/dt;
frames_per_jump = 40/samp_fac;
jump_dt = frames_per_jump*dt;

spk_smwin = round(0.015/dt);

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

load ./expt3_eyedata_NS

big_sac_inds = find(all_sac_amps > max_sac_amp | all_sac_amps > max_sac_amp);
big_sac_start_times = all_t(all_sac_start_inds(big_sac_inds));
big_sac_stop_times = all_t(all_sac_stop_inds(big_sac_inds,1));

blink_start_times = all_t(all_blink_start_inds);
blink_stop_times = all_t(all_blink_stop_inds);

full_image_vec = [];
full_expt_vec = [];
full_trial_vec = [];
full_result_vec = [];
full_seof_vec = [];
full_binned_spks = [];
full_smbinned_spks = [];
full_trial_durs = [];
full_t = [];
full_t_since_start = [];
%%
% sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
sim_sac_blocks = [14 18 24 28 37 40 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

skip_set = [1000 2000 3000 401000 402000 403000];

n_allunits = 96;
trial_cnt = 0;
for ee = 1:length(sim_sac_blocks)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(sim_sac_blocks));
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(ee)));
    
    Trial_starts = [Expts{sim_sac_blocks(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{sim_sac_blocks(ee)}.Trials(:).End]/1e4;
    Trial_seof = [Expts{sim_sac_blocks(ee)}.Trials(:).seof];
    Trial_result = [Expts{sim_sac_blocks(ee)}.Trials(:).Result];
    
    effective_trial_ends = Trial_ends;
    for i = 1:length(Trial_ends)
        cur_blink_times = min(blink_start_times(blink_start_times > Trial_starts(i) & blink_start_times < Trial_ends(i)));
        if isempty(cur_blink_times); cur_blink_times = Inf; end
        cur_big_sac_times = min(big_sac_start_times(big_sac_start_times > Trial_starts(i) & big_sac_start_times < Trial_ends(i)));
        if isempty(cur_big_sac_times); cur_big_sac_times = Inf; end
        effective_trial_ends(i) = min([cur_big_sac_times cur_blink_times]);
    end
    effective_trial_ends(isinf(effective_trial_ends)) = Trial_ends(isinf(effective_trial_ends));
  
    
%     Trial_ends = effective_trial_ends; %end trials at blinks/sacs
    
    endEvents = [Expts{sim_sac_blocks(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    used_trials = find(Trial_durs > min_trial_dur);
    
    bad_trials = [];
    shft_cnt = nan(length(used_trials),1);
    for i = 1:length(used_trials)
        cur_images = Expts{sim_sac_blocks(ee)}.Trials(used_trials(i)).Seedseq;
        
        %for some experiments the seedseq gets truncated at 291 samples
        if length(cur_images) == 291 & length(cur_images)*dt*samp_fac < Trial_durs(used_trials(i))
            extra = 320 - 291;
            cur_images = [cur_images; ones(extra,1)*cur_images(end)];
        end
            
          cur_n_tbins = floor(length(cur_images)/samp_fac);
%           cur_n_tbins = ceil(length(cur_images)/samp_fac);

%         cur_n_tbins = floor(Trial_durs(used_trials(i))/dt);
        cur_n_shifts = floor(cur_n_tbins/frames_per_jump);
        shft_cnt(i) = cur_n_shifts;
        tbin_cnt(i) = cur_n_tbins;
        
%         shift_inds = 1:frames_per_jump:length(cur_images);
%         shift_inds = ceil(1/frames_per_jump:1/frames_per_jump:cur_n_tbins/frames_per_jump);
        shift_inds = ceil((1:cur_n_tbins)/2);
        
        cur_images = cur_images(shift_inds);
        
        if max(mod(cur_images,100)) <= 8
%             bad_trials = [bad_trials i];
            cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*cur_n_tbins);
            cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            
            cur_images(length(cur_t_cents)+1:end) = [];
            
%             [length(cur_t_cents) length(cur_images)]
            cur_binned_spks = nan(n_allunits,length(cur_t_cents));
            cur_smbinned_spks = nan(n_allunits,length(cur_t_cents));
            for j = 1:n_allunits
                temp = histc(Clusters{j}.times,cur_t_edges);
                cur_binned_spks(j,:) = temp(1:end-1);
                cur_smbinned_spks(j,:) = jmm_smooth_1d_cor(temp(1:end-1),spk_smwin);
            end
            
            full_image_vec = [full_image_vec; cur_images];
            full_expt_vec = [full_expt_vec; ones(length(cur_images),1)*sim_sac_blocks(ee)];
            full_trial_vec = [full_trial_vec; ones(length(cur_images),1)*(trial_cnt + used_trials(i))];
            full_result_vec = [full_result_vec; ones(length(cur_images),1)*Trial_result(used_trials(i))];
            full_seof_vec = [full_seof_vec; ones(length(cur_images),1)*Trial_seof(used_trials(i))];
            full_binned_spks = [full_binned_spks; cur_binned_spks'];
            full_smbinned_spks = [full_smbinned_spks; cur_smbinned_spks'];
            full_trial_durs = [full_trial_durs; ones(length(cur_images),1)*Trial_durs(used_trials(i))];
            full_t = [full_t cur_t_cents];
            full_t_since_start = [full_t_since_start cur_t_cents - cur_t_cents(1)];
        end
    end
    trial_cnt = trial_cnt + length(Trial_durs);
    n_bad_trials(ee) = length(bad_trials);
    
end



%%
cd ~/Data/bruce/Expt_1_8_13_imfolder/
% cd ~/Data/bruce/Expt_1_9_13_imfolder/
tot_images = length(unique(full_image_vec));
used_images = unique(full_image_vec);
tot_samps = length(full_image_vec);

all_im_patches = nan(tot_images,length(ypatch_inds),length(xpatch_inds));
all_obj_info = nan(tot_images,length(ypatch_inds),length(xpatch_inds));
full_stim_ids = nan(size(full_t));

stim_cnt = 1;
for i = 1:tot_images
    fprintf('Processing image %d of %d\n',i,tot_images);
    if used_images(i) < 1e4
        filename = sprintf('IM100%.4d.png',used_images(i));
    elseif used_images(i) < 1e5
        filename = sprintf('IM10%.5d.png',used_images(i));
    else
        filename = sprintf('IM1%.6d.png',used_images(i));
    end
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
    
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(full_image_vec == used_images(i));
    cur_n_samps = length(cur_samp_set);
    full_stim_ids(cur_samp_set) = stim_cnt;
    
    cur_patch = IMAGE(ypatch_inds,xpatch_inds);
    all_im_patches(stim_cnt,:,:) = cur_patch;
    
    stim_cnt = stim_cnt + 1;
end
stim_cnt = stim_cnt - 1;
resh_all_stims = reshape(all_im_patches,stim_cnt,sdim^2);

%%
full_insac = ceil(interp1(all_t,all_insac,full_t));
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
cd ~/Data/bruce/G075/
save('Expt3_fixbased_data_all.mat','-v7.3','full_*','trial*','resh_all*','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
