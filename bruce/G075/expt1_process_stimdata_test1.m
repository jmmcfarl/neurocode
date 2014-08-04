clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/G075/
load jbeG075Expts.mat
load ./CellList.mat

dt = 2*117.5/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [-0.52 1.3; -1.3 0.5]; %location of RFs in degrees [x1 x2;y1 y2]

RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
sdim = length(xpatch_inds);

rf_cent = [0.34 -0.43];

min_trial_dur = 1;
max_sac_amp = 0.6;

%%
load ~/Data/bruce/grating_params_fin.mat
for i = 1:3
    grating_seq_ori{i} = [grating_params{i}.ori];
    grating_seq_freq{i} = [grating_params{i}.freq];
    grating_seq_phase{i} = [grating_params{i}.phase];
end

%%
load ./expt1_eyedata_new
% load ./expt1_eyedata_NS
% load ./expt1_eyedata_GR

big_sac_inds = find(all_sac_amps > max_sac_amp);
big_sac_start_times = all_t(all_sac_start_inds(big_sac_inds));
big_sac_stop_times = all_t(all_sac_stop_inds(big_sac_inds,1));

blink_start_times = all_t(all_blink_start_inds);
blink_stop_times = all_t(all_blink_stop_inds);

% hf_start_times = all_t(all_hf_start_inds);
% hf_stop_times = all_t(all_hf_stop_inds);

full_image_vec = [];
full_xo_vec = [];
full_yo_vec = [];
full_expt_vec = [];
full_trial_vec = [];
% full_trial_seof = [];
full_binned_spks = [];
full_trial_durs = [];
full_t = [];
full_im_start_inds = [];
%%
% forty_blocks = [1 6 16 17 20 25 28];
forty_blocks = [7 13 17 19 23 26 29 31 36 38 48] - 6;
% forty_blocks = [41 43 44] - 6;
% forty_blocks = [12 34 52] - 6;
n_allunits = 96;
all_rpt_frames = [];

for ee = 1:length(forty_blocks)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(forty_blocks));
    single_units{ee} = find(CellList(forty_blocks(ee),:,1) > 0);
    n_sus = length(single_units{ee});
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',forty_blocks(ee)));
    
    cur_n_trials = length(Expts{forty_blocks(ee)}.Trials);
    bad_t = [];
    for i = 1:cur_n_trials
        if length(Expts{forty_blocks(ee)}.Trials(i).Start) ~= 1
            bad_t = [bad_t i];
        end
    end
    Expts{forty_blocks(ee)}.Trials(bad_t) = [];
    
    Trial_starts = [Expts{forty_blocks(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{forty_blocks(ee)}.Trials(:).End]/1e4;
    Trial_durs = (Trial_ends-Trial_starts);
    Trial_seof = [Expts{forty_blocks(ee)}.Trials(:).seof];
    Trial_imseed = [Expts{forty_blocks(ee)}.Trials(:).imseed];
    
    bad = find(Trial_seof ~= Trial_imseed);
    if ~isempty(bad)
        disp('Bad trials!')
    end
    
    used_trials = find(Trial_durs > min_trial_dur & Trial_seof==0 & Trial_imseed == 0);
%     used_trials = find(Trial_durs > min_trial_dur);
%     used_trials = find(Trial_durs > min_trial_dur);
    
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
    
    %     check for blinks or big sacs within each trial
    bad_trials = [];
    for i = 1:length(used_trials)
        if ~isempty(find(big_sac_start_times > Trial_starts(used_trials(i)) & ...
                big_sac_start_times < Trial_ends(used_trials(i)), 1)) ...
                | ~isempty(find(big_sac_stop_times >= Trial_starts(used_trials(i)) & ...
                big_sac_stop_times <= Trial_ends(used_trials(i)), 1))
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d sac trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    bad_trials = [];
    for i = 1:length(used_trials)
        if isfield(Expts{forty_blocks(ee)}.Trials(used_trials(i)),'rptframes')
            if ~isempty(Expts{forty_blocks(ee)}.Trials(used_trials(i)).rptframes)
                bad_trials = [bad_trials; i];
            end
        end
    end
    fprintf('Eliminating %d of %d framerep trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    for i = 1:length(used_trials)
        cur_seedseq = [Expts{forty_blocks(ee)}.Trials(used_trials(i)).Seedseq];
%         if mod(length(cur_seedseq),2) ~= 0
%             cur_seedseq(end) = [];
%         end
        cur_images = cur_seedseq(1:2:end);
                
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*length(cur_images));
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
                
        cur_binned_spks = nan(n_mus,length(cur_t_cents));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        full_image_vec = [full_image_vec; cur_images];
        full_expt_vec = [full_expt_vec; ones(length(cur_images),1)*forty_blocks(ee)];
        full_trial_vec = [full_trial_vec; ones(length(cur_images),1)*used_trials(i)];
        full_binned_spks = [full_binned_spks; cur_binned_spks'];
        full_trial_durs = [full_trial_durs; ones(length(cur_images),1)*Trial_durs(used_trials(i))];
        full_t = [full_t cur_t_cents];
    end
end

% temp = mod(full_image_vec,1000);
% bad = find(temp > 160);
% full_image_vec(bad) = [];
% full_expt_vec(bad) = [];
% full_trial_vec(bad) = [];
% full_binned_spks(bad,:) = [];
% full_trial_durs(bad) = [];
% full_t(bad) = [];
%%
% all_stim_ori = nan(size(full_image_vec));
% all_stim_freq = nan(size(full_image_vec));
% all_stim_phase = nan(size(full_image_vec));
% 
% stim_offset = 4000;
% cur_set = find(full_image_vec >= stim_offset & full_image_vec < stim_offset+1000); 
% cur_stim_ids = full_image_vec(cur_set)-stim_offset;
% cur_set(cur_stim_ids > 160) = [];
% cur_stim_ids(cur_stim_ids > 160) = [];
% all_stim_ori(cur_set) = grating_seq_ori{1}(cur_stim_ids);
% all_stim_freq(cur_set) = grating_seq_freq{1}(cur_stim_ids);
% all_stim_phase(cur_set) = grating_seq_phase{1}(cur_stim_ids);
% 
% stim_offset = 5000;
% cur_set = find(full_image_vec >= stim_offset & full_image_vec < stim_offset+1000); 
% cur_stim_ids = full_image_vec(cur_set)-stim_offset;
% cur_set(cur_stim_ids > 160) = [];
% cur_stim_ids(cur_stim_ids > 160) = [];
% all_stim_ori(cur_set) = grating_seq_ori{2}(cur_stim_ids);
% all_stim_freq(cur_set) = grating_seq_freq{2}(cur_stim_ids);
% all_stim_phase(cur_set) = grating_seq_phase{2}(cur_stim_ids);
% 
% stim_offset = 6000;
% cur_set = find(full_image_vec >= stim_offset & full_image_vec < stim_offset+1000); 
% cur_stim_ids = full_image_vec(cur_set)-stim_offset;
% cur_set(cur_stim_ids > 160) = [];
% cur_stim_ids(cur_stim_ids > 160) = [];
% all_stim_ori(cur_set) = grating_seq_ori{3}(cur_stim_ids);
% all_stim_freq(cur_set) = grating_seq_freq{3}(cur_stim_ids);
% all_stim_phase(cur_set) = grating_seq_phase{3}(cur_stim_ids);
% 
% % all_stim_ori = mod(all_stim_ori+90,180);
% % all_stim_ori = mod(all_stim_ori-90,180);

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
full_insac = ceil(interp1(all_t,all_insac,full_t));
full_inblink = ceil(interp1(all_t,all_inblink,full_t));

NT = length(full_image_vec);
min_fix_dur = 0.15;
trial_flips = 1+find(diff(full_trial_vec) ~= 0);

used_inds = ones(size(full_t));
% used_inds = full_insac==0;
used_inds(trial_flips) = 0;

fix_start_inds = 1+find(used_inds(2:end)==1 & used_inds(1:end-1) == 0);
fix_stop_inds = find(used_inds(2:end) == 0 & used_inds(1:end-1)==1);
fix_start_inds = [1 fix_start_inds];
fix_stop_inds = [fix_stop_inds NT];

fix_durs = (fix_stop_inds-fix_start_inds)*dt;

used_fixs = find(fix_durs > min_fix_dur);
fix_start_inds = fix_start_inds(used_fixs);
fix_stop_inds = fix_stop_inds(used_fixs);
fix_durs = fix_durs(used_fixs);


%%
cd ~/Data/bruce/Expt_1_8_13_imfolder/
% cd ~/Data/bruce/Expt_1_9_13_imfolder/
tot_images = length(unique(full_image_vec));
used_images = unique(full_image_vec);
tot_samps = length(full_image_vec);
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('IM0000%.3d.png',used_images(i))
%     filename = sprintf('IM000%.4d.pgm',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
    
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(full_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));    
    for j = 1:length(cur_samp_set)
        
        cur_patch = IMAGE(ypatch_inds,xpatch_inds);
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
    end
    
%     cur_patch = shiftdim(IMAGE(ypatch_inds,xpatch_inds),-1);
%     all_im_patches(cur_samp_set,:,:) = repmat(cur_patch,[length(cur_samp_set) 1 1]);
end

%%
fixed_delay = round(0.05/dt);
used_inds = [];
n_fixs = length(fix_start_inds);
for i = 1:n_fixs
    cur_inds = fix_start_inds(i):fix_stop_inds(i);
    cur_inds = cur_inds((1+fixed_delay):end);
    fix_len(i) = length(cur_inds);
    used_inds = [used_inds; cur_inds(:) - fixed_delay];
end
full_binned_spks = full_binned_spks(used_inds+fixed_delay,:);

fullX = nan(length(used_inds),sdim^2);
cnt = 0;
for i = 1:n_fixs
    fprintf('Compiling fixation %d of %d\n',i,n_fixs);
    cur_set = cnt + (1:fix_len(i));
    fullX(cur_set,:) = reshape(all_im_patches(used_inds(cur_set),:,:),fix_len(i),sdim^2);
    cnt = cnt + fix_len(i);
end

%all these are in STIMULUS time not in delayed RESPONSE TIME
full_image_vec = full_image_vec(used_inds);
full_expt_vec = full_expt_vec(used_inds);
full_trial_vec = full_trial_vec(used_inds);
full_trial_durs = full_trial_durs(used_inds);
full_t = full_t(used_inds);
% full_stim_ori = all_stim_ori(used_inds);
% full_stim_freq = all_stim_freq(used_inds);

%%
% for i = 1:100
%     imagesc(reshape(fullX(i,:),sdim,sdim))
%     set(gca,'ydir','normal')
% %     [full_stim_ori(i) full_stim_freq(i)]
%     mod(270-full_stim_ori(i),180)
% 
%     pause
%     clf
% end
%%
cd ~/Data/bruce/G075/
save('Expt12_compiled_data_fixeddelay_d1p25.mat','-v7.3','full*','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
