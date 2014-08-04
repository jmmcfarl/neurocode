clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G029/
load ./CellList.mat
load ./G029Expts.mat

dt = 118*2/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [-0.52 1.3; -1.3 0.5]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-0.2 1; -1 0.2]; %location of RFs in degrees [x1 x2;y1 y2]
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

load ./all_eyedata_expt1_v2
left_sac_amps = sqrt(sum(all_lsac_camps.^2,2));
right_sac_amps = sqrt(sum(all_rsac_camps.^2,2));
big_sac_inds = find(left_sac_amps > max_sac_amp | right_sac_amps > max_sac_amp);
big_sac_start_times = all_t(all_sac_start_inds(big_sac_inds));
big_sac_stop_times = all_t(all_lsac_stop_inds(big_sac_inds,1));

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

%%
Expt_nu = [1 6 16 17 20 25 28];
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
    
%     bad_trials = [];
%     for i = 1:length(used_trials)
%         if ~isempty(find(hf_start_times >= Trial_starts(used_trials(i)) & ...
%                 hf_start_times <= Trial_ends(used_trials(i)), 1)) ...
%                 | ~isempty(find(hf_stop_times >= Trial_starts(used_trials(i)) & ...
%                 hf_stop_times <= Trial_ends(used_trials(i)), 1))
%             bad_trials = [bad_trials i];
%         end
%     end
%     fprintf('Eliminating %d of %d hf trials\n',length(bad_trials),length(used_trials));
%     used_trials(bad_trials) = [];
        
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
        if isfield(Expts{Expt_nu(ee)}.Trials(used_trials(i)),'rptframes')
            if ~isempty(Expts{Expt_nu(ee)}.Trials(used_trials(i)).rptframes)
                bad_trials = [bad_trials; i];
            end
        end
    end
    fprintf('Eliminating %d of %d framerep trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    %     used_trials(used_trials==length(Trial_durs)) = []; %to prevent errors when testing misalignment
    %      used_trials(used_trials==1) = []; %to prevent errors when testing misalignment
    
    for i = 1:length(used_trials)
        cur_images = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).Seedseq(1:end-1)];
        cur_xo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).xo(1:end-1)];
        cur_yo = [Expts{Expt_nu(ee)}.Trials(used_trials(i)).yo(1:2:end-1)];
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*length(cur_images));
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

%% PARSE DATA INTO FIXATIONS AND COMPUTE WITHIN FIXATION DRIFT
full_insac = ceil(interp1(all_t,all_insac,full_t));
full_inhf = ceil(interp1(all_t,all_inhf,full_t));
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

from_trial_flip = zeros(length(fix_durs),1);
from_trial_flip(ismember(fix_start_inds-1,trial_flips)) = 1;
from_sac = zeros(length(fix_durs),1);
from_sac(from_trial_flip==0) = 1;
from_sac(1) = 0;
sac_inds = find(from_sac==1);
n_sacs = sum(from_sac);

corresp_sac_num = nan(size(from_trial_flip));
sac_amp_left = nan(length(from_trial_flip),2);
sac_amp_right = nan(length(from_trial_flip),2);
for i = 1:n_sacs
    [~,bb] = min(abs(all_t(all_sac_start_inds) - full_t(fix_start_inds(sac_inds(i))-1)));
    corresp_sac_num(sac_inds(i)) = bb;
end
sac_amp_left(sac_inds,:) = all_lsac_camps(corresp_sac_num(sac_inds),:);
sac_amp_right(sac_inds,:) = all_rsac_camps(corresp_sac_num(sac_inds),:);

%%
cur_full_t = full_t;
load ./gabor_tracking_varmeans 
y_cor_vals = y_cor{end}/Fsd;
x_cor_vals = x_cor{end}/Fsd;
load ./Expt1_newcompiled_data_fixeddelay_d1p25_new.mat full_t %load old t-axis
y_cor_interp = interp1(full_t,y_cor_vals,cur_full_t);
x_cor_interp = interp1(full_t,x_cor_vals,cur_full_t);
x_cor_interp(isnan(x_cor_interp)) = 0;
y_cor_interp(isnan(y_cor_interp)) = 0;

full_t = cur_full_t;

cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(full_image_vec));
used_images = unique(full_image_vec);
tot_samps = length(full_image_vec);
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
        
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(full_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
       
    for j = 1:length(cur_samp_set)
        
        ypatch_inds_adj = round(ypatch_inds - full_yo_vec(cur_samp_set(j))*Fsd + x_cor_interp(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - full_xo_vec(cur_samp_set(j))*Fsd + y_cor_interp(cur_samp_set(j))*Fsd);
        
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
    end
end

%%
flen = 5;
sdim = length(xpatch_inds);
klen = sdim^2*flen;
used_inds = [];
maxdur = sum(fix_stop_inds-fix_start_inds);

fullX = nan(maxdur,klen);
nfull_binned_spks = nan(maxdur,96);

n_fixs = length(fix_start_inds);
used_inds = [];
cnt = 0;
for i = 1:n_fixs
    fprintf('%d of %d\n',i,n_fixs)
    cur_inds = fix_start_inds(i):fix_stop_inds(i);
    temp = makeStimRows(all_im_patches(cur_inds,:,:),flen,1);
    fullX(cnt + (1:size(temp,1)),:) = temp;
    nfull_binned_spks(cnt + (1:size(temp,1)),:) = full_binned_spks(cur_inds(flen:end),:);
    used_inds = [used_inds; cur_inds(flen:end)'];
    cnt = cnt + size(temp,1);
end
fullX(cnt+1:end,:) = [];
nfull_binned_spks(cnt+1:end,:) = [];

full_binned_spks = nfull_binned_spks;

%%
% fixed_delay = round(0.05/dt);
% used_inds = [];
% n_fixs = length(fix_start_inds);
% for i = 1:n_fixs
%     cur_inds = fix_start_inds(i):fix_stop_inds(i);
%     cur_inds = cur_inds((1+fixed_delay):end);
%     fix_len(i) = length(cur_inds);
%     used_inds = [used_inds; cur_inds(:) - fixed_delay];
% end
% full_binned_spks = full_binned_spks(used_inds+fixed_delay,:);
% 
% fullX = nan(length(used_inds),sdim^2);
% cnt = 0;
% for i = 1:n_fixs
%     fprintf('Compiling fixation %d of %d\n',i,n_fixs);
%     cur_set = cnt + (1:fix_len(i));
%     fullX(cur_set,:) = reshape(all_im_patches(used_inds(cur_set),:,:),fix_len(i),sdim^2);
%     cnt = cnt + fix_len(i);
% end

%all these are in STIMULUS time not in delayed RESPONSE TIME
full_image_vec = full_image_vec(used_inds);
full_xo_vec = full_xo_vec(used_inds);
full_yo_vec = full_yo_vec(used_inds);
full_expt_vec = full_expt_vec(used_inds);
full_trial_vec = full_trial_vec(used_inds);
full_trial_durs = full_trial_durs(used_inds);
full_t = full_t(used_inds);

%%
cd ~/Data/bruce/7_15_12/G029/
save('Expt1_corrected_data_flen6_d1p25.mat','-v7.3','full*','sac_*','corresp*','*cor_interp','from*','dt','Fsd','RF_patch*','xax','yax','*_inds','dsfrac');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
