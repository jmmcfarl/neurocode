clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

% cd G029/
% load ./CellList.mat
% load ./G029Expts.mat

cd G034/
load ./CellList.mat
load ./G034Expts.mat

dt = 118/1e4;
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

min_trial_dur = 2;
max_sac_amp = 0.5;

load ./all_eyedata_expt1_34

big_sac_inds = find(all_sac_amps > max_sac_amp);
big_sac_start_times = all_t(all_sac_start_inds(big_sac_inds));
big_sac_stop_times = all_t(all_sac_stop_inds(big_sac_inds,1));

blink_start_times = all_t(all_blink_start_inds);
blink_stop_times = all_t(all_blink_stop_inds);

hf_start_times = all_t(all_hf_start_inds);
hf_stop_times = all_t(all_hf_stop_inds);

full_image_vec = [];
full_xo_vec = [];
full_yo_vec = [];
full_expt_vec = [];
full_trial_vec = [];
full_binned_spks = [];
full_trial_durs = [];
full_t = [];
% full_im_start_inds = [];
full_im_ids= [];
%%
% Expt_nu = [1 6 16 17 20 25 28];
Expt_nu = [1 2 17 19 23 24]; %expt 1 34
n_allunits = 96;
    all_rpt_frames = [];

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
    used_trials = find(Trial_durs > min_trial_dur);
        
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
    
    for i = 1:length(used_trials)
        cur_images = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).Seedseq(1:2:end-1)];
        cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo(1:end-1)];
        cur_yo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).yo(1:2:end-1)];
                     
        cur_repeat_frames = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).rptframes];
        n_rpt_frames(used_trials(i)) = length(cur_repeat_frames);
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i))+dt*(171+n_rpt_frames(used_trials(i))));
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
%         cur_im_start = Trial_starts(used_trials(i)):2*dt:(Trial_starts(used_trials(i)) + 2*dt*(length(cur_images)-1));
        cur_im_start_inds = 1:2:2*length(cur_images);
        rpt_frame_inds = ceil(cur_repeat_frames/2);
        for cc = 1:length(cur_repeat_frames)
%             cur_im_starts(rpt_frame_inds(cc)+1:end) = cur_im_starts(rpt_frame_inds(cc)+1:end) + dt;
            cur_im_start_inds(rpt_frame_inds(cc)+1:end) = cur_im_start_inds(rpt_frame_inds(cc)+1:end) + 1;
        end
%         all_rpt_frames = [all_rpt_frames; cur_repeat_frames];
        
        cur_binned_spks = nan(n_mus,length(cur_t_cents));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
%         cur_binned_spks(:,del_bins) = [];
        
        cur_im_ids = zeros(length(cur_t_cents),1);
        for j = 1:length(cur_images)
           cur_im_ids(cur_im_start_inds(j):end) = cur_im_ids(cur_im_start_inds(j):end) + 1; 
        end
        
%         full_im_start_inds = [full_im_start_inds; cur_im_start_inds(:)+size(full_binned_spks,1)];
if ~isempty(full_im_ids)
cur_offset = max(full_im_ids);
else
    cur_offset = 0;
end
        full_im_ids = [full_im_ids; cur_im_ids+cur_offset];
        full_image_vec = [full_image_vec; cur_images(cur_im_ids)];
        full_xo_vec = [full_xo_vec; cur_xo(cur_im_ids)'];
        full_yo_vec = [full_yo_vec; cur_yo(cur_im_ids)'];
        full_expt_vec = [full_expt_vec; ones(length(cur_im_ids),1)*Expt_nu(ee)];
        full_trial_vec = [full_trial_vec; ones(length(cur_im_ids),1)*used_trials(i)];
        full_binned_spks = [full_binned_spks; cur_binned_spks'];
        full_trial_durs = [full_trial_durs; ones(length(cur_im_ids),1)*Trial_durs(used_trials(i))];
%         full_t = [full_t cur_t_cents];
        full_t = [full_t cur_t_edges(1:end-1)];
    end
end

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
full_insac = ceil(interp1(all_t,all_insac,full_t));
% full_inhf = ceil(interp1(all_t,all_inhf,full_t));
full_inblink = ceil(interp1(all_t,all_inblink,full_t));

NT = length(full_image_vec);
min_fix_dur = 0.15;
trial_flips = 1+find(diff(full_trial_vec) ~= 0);

used_inds = ones(size(full_t));
% used_inds = full_insac==0;
used_inds(trial_flips) = 0;

fix_start_inds = find(used_inds(2:end)==1 & used_inds(1:end-1) == 0);
fix_stop_inds = find(used_inds(2:end) == 0 & used_inds(1:end-1)==1);
fix_start_inds = [1 fix_start_inds];
fix_stop_inds = [fix_stop_inds NT];

fix_durs = (fix_stop_inds-fix_start_inds+1)*dt;

used_fixs = find(fix_durs > min_fix_dur);
fix_start_inds = fix_start_inds(used_fixs);
fix_stop_inds = fix_stop_inds(used_fixs);
fix_durs = fix_durs(used_fixs);

% from_trial_flip = zeros(length(fix_durs),1);
% from_trial_flip(ismember(fix_start_inds-1,trial_flips)) = 1;
% from_sac = zeros(length(fix_durs),1);
% from_sac(from_trial_flip==0) = 1;
% from_sac(1) = 0;
% sac_inds = find(from_sac==1);
% n_sacs = sum(from_sac);

% corresp_sac_num = nan(size(from_trial_flip));
% sac_amp_left = nan(length(from_trial_flip),2);
% sac_amp_right = nan(length(from_trial_flip),2);
% for i = 1:n_sacs
%     [~,bb] = min(abs(all_t(all_sac_start_inds) - full_t(fix_start_inds(sac_inds(i))-1)));
%     corresp_sac_num(sac_inds(i)) = bb;
% end
% sac_amp_left(sac_inds,:) = all_lsac_camps(corresp_sac_num(sac_inds),:);
% sac_amp_right(sac_inds,:) = all_rsac_camps(corresp_sac_num(sac_inds),:);

%%
cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(full_image_vec));
used_images = unique(full_image_vec);
tot_samps = length(unique(full_im_ids));
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
        
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(full_image_vec == used_images(i));
    [cur_samp_im_ids,cur_samp_full_ids] = unique(full_im_ids(cur_samp_set));
    cur_samp_full_ids = cur_samp_set(cur_samp_full_ids);
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
       
    for j = 1:length(cur_samp_im_ids)
        
        ypatch_inds_adj = round(ypatch_inds - full_yo_vec(cur_samp_full_ids(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - full_xo_vec(cur_samp_full_ids(j))*Fsd);
        
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        all_im_patches(cur_samp_im_ids(j),:,:) = cur_patch;
    end
end

%%
% fixed_delay = round(0.05/dt);
used_inds = [];
n_fixs = length(fix_start_inds);
for i = 1:n_fixs
    cur_inds = fix_start_inds(i):fix_stop_inds(i);
%     cur_inds = cur_inds((1+fixed_delay):end);
    fix_len(i) = length(cur_inds);
    used_inds = [used_inds; cur_inds(:)];
%     used_inds = [used_inds; cur_inds(:) - fixed_delay];
end
% rid = find(full_im_start_inds(used_inds) >= size(full_binned_spks,1)-fixed_delay);
% used_inds(rid) = [];
% fix_len(end) = fix_len(end)-length(rid);
% full_binned_spks = full_binned_spks(used_inds+fixed_delay,:);
% avg_binned_spks = full_binned_spks(full_im_start_inds(used_inds)+fixed_delay,:)+full_binned_spks(full_im_start_inds(used_inds)+fixed_delay+1,:);
% full_binned_spks = avg_binned_spks;
% fullX = nan(length(used_inds),sdim^2);
% cnt = 0;
% for i = 1:n_fixs
%     fprintf('Compiling fixation %d of %d\n',i,n_fixs);
%     cur_set = cnt + (1:fix_len(i));
%     fullX(cur_set,:) = reshape(all_im_patches(used_inds(cur_set),:,:),fix_len(i),sdim^2);
%     cnt = cnt + fix_len(i);
% end
% 
fullX = reshape(all_im_patches,tot_samps,sdim^2);

% %all these are in STIMULUS time not in delayed RESPONSE TIME
% full_image_vec = full_image_vec(used_inds);
% full_xo_vec = full_xo_vec(used_inds);
% full_yo_vec = full_yo_vec(used_inds);
% full_expt_vec = full_expt_vec(used_inds);
% full_trial_vec = full_trial_vec(used_inds);
% full_trial_durs = full_trial_durs(used_inds);
% full_t = full_t(full_im_start_inds(used_inds));

%%
cd ~/Data/bruce/7_15_12/G034/
save('Expt1_newcompiled_data_fixeddelay_d1p25_34_ht.mat','-v7.3','full*','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
