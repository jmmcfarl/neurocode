clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

Expt_name = 'G087';
bar_ori = 0; %bar orientation to use (only for UA recs)

% fig_dir = '/Users/james/Analysis/bruce/variability/figures/';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';
mod_data_name = 'corrected_models2';

%%

micro_thresh = 1; %max amp of microsac (deg)
EP_bounds = 1;%eye position boundary (deg from central FP)
sac_burst_isi = 0.15;
max_gsac_dur = 0.1;

%%

Expt_num = str2num(Expt_name(2:end));
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
        case 296
            bar_ori = 45;
        case 297
            if ~ismember(bar_ori,[0 90])
                error('M297 is either 0 or 90 deg bars');
            end
    end
end

if Expt_num >= 280
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %m270 and 266
    end
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
    use_nPix = 32;
    use_LOOXV = 1;
elseif strcmp(rec_type,'UA')
    rpt_seed = nan;
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
    use_nPix = 16;
    use_LOOXV = 1;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end

et_dir = ['~/Analysis/bruce/' 'G093' '/ET_final_imp/'];
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

% et_mod_data_name = 'full_eyetrack_initmods';
% et_anal_name = 'full_eyetrack';
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_hres_anal_name = [et_anal_name '_hres'];

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
et_hres_anal_name = [et_hres_anal_name sprintf('_ori%d',bar_ori)];

%dont fit stim models using these blocks
ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 289
        ignore_blocks = [27]; %27 is off somehow
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28];
end

%problem with M270 where vd was wrong, need this correction factor to get
%correct units
if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end


is_TBT_expt = false;
if Expt_num >= 275
    is_TBT_expt = true;
end

%%

flen = 15;
spatial_usfac = 2;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
elseif ismember(Expt_num,[296 297])
    use_nPix = 22;
    spatial_usfac = 2;
end

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;
new_dt = 0.0025;
Fr = 1;

backlag = round(0.1/dt);
forlag = round(0.3/dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);

full_nPix=36;
switch Expt_num
    case 270
        full_nPix=32;
    case  287
        full_nPix = 22;
    case 289
        full_nPix = 22;
    case 294
        full_nPix = 20;
        %     case 296
        %         full_nPix = 54;
end

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_sac_amp = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_sac_amp(ii) = Expts{ii}.Stimvals.Fs;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
if length(unique(expt_dds(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dds(cur_block_set));
    cur_block_set(expt_dds(cur_block_set) ~= main_dds) = [];
end


sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

n_blocks = length(cur_block_set);

all_nfs = cellfun(@(x) x.Stimvals.nf,Expts(cur_block_set));
if length(unique(all_nfs)) > 1
    fprintf('Warning, multiple different nfs detected: %.4f\n',all_nfs);
end


if ~isempty(grayback_gs_expts) || ~isempty(imback_gs_expts)
    gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
else
    gsac_amp = unique(expt_sac_amp(cur_block_set));
end
if length(gsac_amp) > 1
    fprintf('Multiple guided sac amps detected!\n');
end
%minimum (parallel) amplitude for a guided saccade to be included in
%analysis
gsac_thresh = mean(gsac_amp)/2;

%%
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/scale_fac; %model dx in deg

%%
full_nPix_us = spatial_usfac*full_nPix;

sparsity = 0.12;
n_rpts = 200;
rpt_dur = round(4/dt);
NT = rpt_dur*n_rpts;

rand_stim = zeros(rpt_dur,full_nPix);
rand_stuff = ones(rpt_dur,full_nPix);
rand_stuff(rand(rpt_dur,full_nPix) > 0.5) = -1;
is_nonzero = rand(rpt_dur,full_nPix) < sparsity;
rand_stim(is_nonzero) = rand_stuff(is_nonzero);

if spatial_usfac > 1
    all_stimmat_up = zeros(size(rand_stim,1),full_nPix_us);
    for ii = 1:size(rand_stim,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = rand_stim(:,ii);
        end
    end
elseif spatial_usfac == 1
    all_stimmat_up = rand_stim;
end

stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

rpt_stim = permute(repmat(all_stimmat_up,[1 1 n_rpts]),[1 3 2]);
rpt_stim = reshape(rpt_stim,NT,[]);

%% select submatrix with central pixels
buffer_pix = floor((full_nPix - use_nPix)/2);

%repeat for up-sampled versions of the Xmatrix
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

%%
fprintf('Loading model fits\n');
cd(mod_data_dir)
load(mod_data_name);

cd(et_dir)
load sim_ep_data
full_ep = hres_ep;

rng(1);
% hres_ep = circshift(hres_ep(:),2e5);
rpt_ep = reshape(hres_ep(1:NT),rpt_dur,n_rpts);
rpt_ep = rpt_ep(:,randperm(n_rpts));
hres_ep = rpt_ep(:);

fin_shift_cor = round(hres_ep/sp_dx);
fin_shift_cor(fin_shift_cor > full_nPix_us) = full_nPix_us;
fin_shift_cor(fin_shift_cor < -full_nPix_us) = -full_nPix_us;
%%
ep_rpt_stim = rpt_stim;
for i=1:NT
    ep_rpt_stim(i,:) = shift_matrix_Nd(rpt_stim(i,:),-fin_shift_cor(i),2);
end

%%
all_Xmat = create_time_embedding(rpt_stim,stim_params_us);
all_Xmat = all_Xmat(:,use_kInds_up);

all_ep_Xmat = create_time_embedding(ep_rpt_stim,stim_params_us);
all_ep_Xmat = all_ep_Xmat(:,use_kInds_up);

%%
cc = 101;
target_mod = ModData(cc).rectGQM;
[~,~,noep_prate] = NMMmodel_eval(target_mod,[],all_Xmat);
[~,~,ep_prate] = NMMmodel_eval(target_mod,[],all_ep_Xmat);

noep_prate = reshape(noep_prate,rpt_dur,n_rpts);
ep_prate = reshape(ep_prate,rpt_dur,n_rpts);

%%
% % spk_jitter = 0.005;
% spk_jitter = 0.005;
% 
% rand_spks = poissrnd(noep_prate);
% rand_ep_spks = poissrnd(ep_prate);
% 
% rand_spk_bins = convert_to_spikebins(rand_spks(:));
% rand_spk_trials = ceil(rand_spk_bins/rpt_dur);
% rand_spk_times = mod(rand_spk_bins,rpt_dur)*dt + randn(size(rand_spk_bins))*spk_jitter;
% 
% rand_ep_spk_bins = convert_to_spikebins(rand_ep_spks(:));
% rand_ep_spk_trials = ceil(rand_ep_spk_bins/rpt_dur);
% rand_ep_spk_times = mod(rand_ep_spk_bins,rpt_dur)*dt + randn(size(rand_ep_spk_bins))*spk_jitter;

%%
new_dt = 0.001;
old_tax = (0:rpt_dur-1)*dt;
new_tax = (0:new_dt:old_tax(end));

new_dur = length(new_tax);

noep_prate_interp = interp1(old_tax,noep_prate,new_tax)*new_dt/dt;
ep_prate_interp = interp1(old_tax,ep_prate,new_tax)*new_dt/dt;

rand_spks = poissrnd(noep_prate_interp);
rand_ep_spks = poissrnd(ep_prate_interp);

rand_spk_bins = convert_to_spikebins(rand_spks(:));
rand_spk_trials = ceil(rand_spk_bins/new_dur);
rand_spk_times = mod(rand_spk_bins,new_dur)*new_dt;

rand_ep_spk_bins = convert_to_spikebins(rand_ep_spks(:));
rand_ep_spk_trials = ceil(rand_ep_spk_bins/new_dur);
rand_ep_spk_times = mod(rand_ep_spk_bins,new_dur)*new_dt;

psth_sm = 0.005/new_dt;
psth = hist(rand_spk_times,new_tax)/n_rpts;
psth_ep = hist(rand_ep_spk_times,new_tax)/n_rpts;
psth = jmm_smooth_1d_cor(psth,psth_sm);
psth_ep = jmm_smooth_1d_cor(psth_ep,psth_sm);

%%
close all
var(psth_ep)/var(psth)

line_height = 0.9;
line_width = 1;

msize = 4;
% xl = [1.5 2.5];
xl = [3 4];
ca = [0 100];

f1 = figure();
subplot(3,2,1)
imagesc(new_tax,1:n_rpts,ep_prate_interp'/new_dt);
xlim(xl);
caxis(ca);
% colorbar
xlabel('Time (s)');

subplot(3,2,2); hold on
% plot(rand_ep_spk_times,rand_ep_spk_trials,'k.','markersize',msize)
for ii = 1:n_rpts
    cur_spks = find(rand_ep_spk_trials == ii);
    for jj = 1:length(cur_spks)
       line(rand_ep_spk_times(cur_spks(jj)) + [0 0],ii + [0 line_height],'color','k','linewidth',line_width);
    end
end
xlim(xl);
ylim([0 n_rpts]);

subplot(3,2,3)
imagesc(new_tax,1:n_rpts,noep_prate_interp'/new_dt);
xlim(xl);
caxis(ca);
% colorbar

subplot(3,2,4); hold on
% plot(rand_spk_times,rand_spk_trials,'k.','markersize',msize)
for ii = 1:n_rpts
    cur_spks = find(rand_spk_trials == ii);
    for jj = 1:length(cur_spks)
       line(rand_spk_times(cur_spks(jj)) + [0 0],ii + [0 line_height],'color','k','linewidth',line_width);
    end
end
xlim(xl);
ylim([0 n_rpts]);

subplot(3,2,5); hold on
plot(new_tax,var(ep_prate_interp')/new_dt^2)
plot(new_tax,var(noep_prate_interp')/new_dt^2,'r')
xlim(xl);
set(gca,'YAxisLocation','right');

subplot(3,2,6); hold on
plot(new_tax,psth_ep/new_dt)
plot(new_tax,psth/new_dt,'r')
xlim(xl);
ylim(ca)

% %PRINT FIGURE
fig_width = 7; rel_height = 1.4;
figufy(f1);
fname = [fig_dir 'Example_illust_rasters3.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
% 
%%
% cd ~/Analysis/bruce/G093/
% sname = 'simtemp';
% save new_tax psth* new_dt rand_spk* *prate* n_rpts rand_* 
%%
% f2 = figure();
% imagesc(old_tax,1:n_rpts,rpt_ep');
% caxis([-0.5 0.5]);

% f3 = figure();
% t1 = 10;
% t2 = 6;
% t3 = 28;
% % for t1 = 9:100
% plot(old_tax,rpt_ep(:,t1));hold on
% plot(old_tax,rpt_ep(:,t2),'r')
% plot(old_tax,rpt_ep(:,t3),'k')
% xlim(xl);
% ylim([-0.4 0.4])
% % pause
% % clf
% % end
% xlabel('Time (s)');
% ylabel('Eye position (deg)');
% 
% f4 = figure();
% hist(full_ep,100);
% xlim([-0.4 0.4])
% xlabel('Eye position (deg)');
% ylabel('Relative frequency');

% % %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.9;
% figufy(f3);
% fname = [fig_dir 'Example_eyetraces.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
% 
% figufy(f4);
% fname = [fig_dir 'Example_eyedist.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);