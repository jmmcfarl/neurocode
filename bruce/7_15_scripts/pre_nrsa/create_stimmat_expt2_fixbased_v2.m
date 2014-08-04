clear all
close all
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12
obj_info_dir = '~/James_scripts/data_processing/Images/object_im_info';
obj_set = 1142:1369;

% cd G029/
% load ./CellList.mat
% load ./G029Expts.mat
cd G034/
load ./CellList.mat
load ./G034Expts.mat

Pix2Deg = 0.018837;
dsfrac = 1.5;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [-1.5 2.2; -2.2 1.5]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
sdim = length(xpatch_inds);

% load ./expt2_parsed_g029
load ./expt2_parsed_g034

dt_dsf = 4;
dt = all_t(2)-all_t(1);
dtd = dt*dt_dsf;
min_fix_dur = 0.15;
fix_dur_win = 0.15;

full_image_vec = [];
full_xo_vec = [];
full_yo_vec = [];
full_angles = [];
full_jump_sizes = [];
full_expt_vec = [];
full_trial_vec = [];
full_trial_durs = [];

n_allunits = 96;
single_units = find(CellList(1,:,1) > 0);
n_sus = length(single_units);
multi_units = setdiff(1:n_allunits,single_units);

%%
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_expts = length(Expt_nu);

full_binned_spks = [];
full_t = [];
for i = 1:n_expts
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(i)));
    
    cur_set = find(all_expt_vec==i);
    cur_t_edges = downsample(all_t(cur_set),dt_dsf);
    full_t = [full_t; cur_t_edges];
    
    cur_binned_spks = nan(96,length(cur_t_edges));
    for j = 1:96
        temp = histc(Clusters{j}.times,cur_t_edges);
        cur_binned_spks(j,:) = temp;
    end
    full_binned_spks = [full_binned_spks; cur_binned_spks'];
end

full_image_vec = round(interp1(all_t,all_imnum,full_t));
full_expt_vec = round(interp1(all_t,all_expt_vec,full_t));
full_trial_vec = round(interp1(all_t,all_trial_vec,full_t));
full_eyepos = interp1(all_t,all_eyepos,full_t);
full_blink_inds = round(interp1(full_t,1:length(full_t),all_blink_times));
full_sac_inds = round(interp1(full_t,1:length(full_t),all_sac_times));
full_out_inds = round(interp1(full_t,1:length(full_t),all_out_times));

%% PARSE INTO FIXATIONS
use_data = ones(size(all_expt_vec));
use_data(all_inblink==1) = nan;
use_data(all_inout==1) = 0;
use_data(all_insac==1) = 0;
use_data(isnan(all_imnum)) = 0;
use_imdata = ones(size(all_expt_vec));
use_imdata(isnan(all_imnum)) = 0;

fix_start_inds = 1 + find(use_data(1:end-1) == 0 & use_data(2:end) == 1);
im_start_inds = 1 + find(use_imdata(1:end-1) == 0 & use_imdata(2:end) == 1);
if use_data(1) == 1 
    fix_start_inds = [1; fix_start_inds];
end

poss_ends = 1 + find(use_data(1:end-1) == 1 & use_data(2:end) ~= 1);
fix_end_inds = zeros(size(fix_start_inds));
for i = 1:length(fix_start_inds)
   fix_term = find(poss_ends >= fix_start_inds(i),1,'first');
   fix_end_inds(i) = poss_ends(fix_term);
end

fix_durs = (fix_end_inds-fix_start_inds)*(all_t(2)-all_t(1));
too_short = find(fix_durs < min_fix_dur);
fix_start_inds(too_short) = [];
fix_end_inds(too_short) = [];
fix_start_times = all_t(fix_start_inds);
fix_end_times = all_t(fix_end_inds);
fix_is_imflip = ismember(fix_start_inds,im_start_inds);

full_fix_starts = round(interp1(full_t,1:length(full_t),fix_start_times));
full_fix_ends = round(interp1(full_t,1:length(full_t),fix_end_times));

full_fix_wends = full_fix_starts + round(fix_dur_win/dtd);

n_fixs = length(full_fix_wends);
fix_eyepos = nan(n_fixs,2);
fix_imnum = nan(n_fixs,1);
for i = 1:n_fixs
   cur_set = full_fix_starts(i):full_fix_wends(i);
   fix_eyepos(i,:) = median(full_eyepos(cur_set,1:2));
   uset = ~isnan(full_image_vec(cur_set));
   fix_imnum(i) = unique(full_image_vec(cur_set(uset))); 
end

fix_is_sac = ~fix_is_imflip;
fix_sac_amps = nan(n_fixs,1);
for i = 1:n_fixs
    if fix_is_sac(i)==1
       prev_sac = find(all_sac_times(:,1) < fix_start_times(i),1,'last');
       fix_sac_amps(i) = all_sacamp(prev_sac);
    end
end
%%
cd ~/James_scripts/data_processing/Images/image_set_A
used_images = unique(fix_imnum);
tot_images = length(used_images);
NT = length(full_image_vec);

full_im_patches = nan(n_fixs,length(ypatch_inds),length(xpatch_inds));
full_obj_info = nan(n_fixs,length(ypatch_inds),length(xpatch_inds));
full_fix_ids = nan(size(full_t));

for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
        
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    if ismember(used_images(i),obj_set) & used_images(i) ~= 1167
        cd(obj_info_dir)
        obj_info_id = find(obj_set == used_images(i));
        fname = sprintf('info4%.4d.mat',obj_info_id);
        load(fname);
        
        obj_IMAGE = ov_L;
        obj_IMAGE = flipud(obj_IMAGE);
        obj_IMAGE = round(imresize(obj_IMAGE,1/dsfrac));
    end
    
    cur_fix_set = find(fix_imnum == used_images(i));
    cur_n_fixs = length(cur_fix_set);
   
    fprintf('Analyzing image %d, %d unique fixs\n',i,cur_n_fixs);
    for j = 1:cur_n_fixs
        
        ypatch_inds_adj = round(ypatch_inds - fix_eyepos(cur_fix_set(j),2)*Fsd);
        xpatch_inds_adj = round(xpatch_inds - fix_eyepos(cur_fix_set(j),1)*Fsd);
        
        cur_used_inds = full_fix_starts(cur_fix_set(j)):full_fix_ends(cur_fix_set(j));
        full_fix_ids(cur_used_inds) = cur_fix_set(j);
        
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        full_im_patches(cur_fix_set(j),:,:) = cur_patch;
        
        if ismember(used_images(i),obj_set) & used_images(i) ~= 1167
            cur_obj_patch = obj_IMAGE(ypatch_inds_adj,xpatch_inds_adj);
            full_obj_info(cur_fix_set(j),:,:) = cur_obj_patch;
        end
    end
end

resh_all_stims = reshape(full_im_patches,n_fixs,sdim^2);
resh_all_obj = reshape(full_obj_info,n_fixs,sdim^2);
clear full_im_patches full_obj_info

%%
% cd ~/Data/bruce/7_15_12/G029/
% save('Expt2_compiled_windata_d1p5.mat','-v7.3','full_*','resh_all*','Fsd','RF_patch*','xax','yax','*_inds','dsfrac','dt*','fix*');
cd ~/Data/bruce/7_15_12/G034/
save('Expt2_compiled_windata_d1p5.mat','-v7.3','full_*','resh_all*','Fsd','RF_patch*','xax','yax','*_inds','dsfrac','dt*','fix*');
%%
% fullX = makeStimRows(all_im_patches,flen,0);
% full_binned_spks = all_binned_spks;
%
% cd ~/Data/bruce/7_15_12/G029/
% save Expt1_compiled_data_ms full* flen dt Fsd
