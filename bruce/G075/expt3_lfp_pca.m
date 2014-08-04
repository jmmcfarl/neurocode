clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 150;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 40]/niqf);
use_lfps = [1:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);


%% COMPILE TRIAL-LEVEL DATA
% sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
sim_sac_blocks = [14 24 37 50] - 6; %sim sac b
% sim_sac_blocks = [14 18 24 28 37 40 50] - 6;

repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
all_expt_id = [];
for bb = 1:length(sim_sac_blocks)
    n_trials(bb) = length(Expts{sim_sac_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id sim_sac_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_seof = cell2mat(trial_seof);
all_trial_completed = cell2mat(trial_completed);
all_trial_durs = cell2mat(trial_durs);

min_trial_dur = 1;
used_trials = find(~ismember(all_trial_seof,repeat_inds) & all_trial_durs > min_trial_dur);
%% COMPILE STIMULUS-LEVEL DATA
rel_start_times = (0:40:280)/stim_fs;
stim_dur = 40/stim_fs;
stims_per_trial = 8;
tot_n_trials = sum(n_trials);
all_stim_start_times = [];
all_stim_rel_num = [];
all_stim_trial_num = [];
all_stim_block_num = [];
all_stim_seof = [];
all_stim_id = [];
for tt = 1:length(used_trials)
    cur_n_stims = floor(all_trial_durs(used_trials(tt))/stim_dur);
    for st = 1:cur_n_stims
        all_stim_start_times = [all_stim_start_times; all_trial_start_times(used_trials(tt))+rel_start_times(st)];
        all_stim_rel_num = [all_stim_rel_num; st];
        all_stim_trial_num = [all_stim_trial_num; used_trials(tt)];
        all_stim_block_num = [all_stim_block_num; all_expt_id(used_trials(tt))];
        all_stim_seof = [all_stim_seof; all_trial_seof(used_trials(tt))];
        all_stim_id = [all_stim_id; all_trial_seof(used_trials(tt)) + st];
    end
end


%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
    all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2);
    %    all_dists(ll) = inf;
    [~,best_loc] = min(all_dists(use_lfps));
    nearest_lfps(ll) = best_loc;
end

%%
stim_dur = 0.47;
trial_cnt = 1;
% desired_dt = 0.003;
desired_dt = 0.006;
all_interp_lfp = [];
all_binned_spks = [];
all_expt_t = [];
all_expt_seof = [];
all_expt_stimid = [];
all_expt_tnum = [];
all_t_since_start = [];
all_trial_id = [];
for bb = 1:length(sim_sac_blocks)
    
    %%
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    cur_stims = find(all_stim_block_num==sim_sac_blocks(bb));
    
    expt_t_axis = [];
    expt_seof = [];
    expt_tnum = [];
    expt_stimid = [];
    t_since_start = [];
    for tt = 1:length(cur_stims)
        cur_t_edges = all_stim_start_times(cur_stims(tt)):desired_dt:(all_stim_start_times(cur_stims(tt)) + stim_dur);
        binned_spks = nan(length(use_units),length(cur_t_edges)-1);
        for cc = 1:length(use_units)
            temp = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            binned_spks(cc,:) = temp(1:end-1);
        end
        all_binned_spks = cat(1,all_binned_spks,binned_spks');
        all_trial_id = [all_trial_id; ones(length(cur_t_edges)-1,1)*trial_cnt];
        expt_seof = [expt_seof; ones(length(cur_t_edges)-1,1)*all_stim_seof(cur_stims(tt))];
        expt_tnum = [expt_tnum; ones(length(cur_t_edges)-1,1)*all_stim_trial_num(cur_stims(tt))];
        expt_stimid = [expt_stimid; ones(length(cur_t_edges)-1,1)*all_stim_id(cur_stims(tt))];
        trial_cnt = trial_cnt + 1;
        
        expt_t_axis = [expt_t_axis cur_t_edges(1:end-1)+desired_dt/2];
        t_since_start = [t_since_start cur_t_edges(1:end-1)-cur_t_edges(1)+desired_dt/2];
    end
    
    interp_lfps = interp1(t_ax,Vmat',expt_t_axis);
    
    all_interp_lfp = cat(1,all_interp_lfp,interp_lfps);
    all_expt_t = cat(1,all_expt_t,expt_t_axis');
    all_expt_seof = [all_expt_seof; expt_seof(:)];
    all_expt_stimid = [all_expt_stimid; expt_stimid(:)];
    all_expt_tnum = [all_expt_tnum; expt_tnum(:)];
    all_t_since_start = cat(1,all_t_since_start,t_since_start');
    clear Vmat
    
end

%%
new_trial_inds = [1; 1 + find(diff(all_expt_tnum) ~= 0)];
new_stim_inds = [1; 1 + find(diff(all_t_since_start) < 0)];

used_stim_inds = new_stim_inds(~ismember(new_stim_inds,new_trial_inds) & ~isnan(all_interp_lfp(new_stim_inds,1)));
% used_stim_inds = new_trial_inds(~isnan(all_interp_lfp(new_trial_inds,1)));

backlag = round(0/desired_dt);
forwardlag = round(0.47/desired_dt);
for cc = 1:length(use_lfps)
    [ev_avg(cc,:),lags] = get_event_trig_avg(all_interp_lfp(:,cc),used_stim_inds,backlag,forwardlag);
end


%%
cur_ch = 1;
[ev_mat,lags] = get_event_trig_mat(all_interp_lfp(:,cur_ch),used_stim_inds,backlag,forwardlag);

trig_cov = cov(ev_mat);

rand_inds = ceil(rand(10000,1)*size(all_interp_lfp,1));
rand_inds(rand_inds < backlag | rand_inds > length(all_expt_t) - forwardlag) = [];
[rand_ev_mat,lags] = get_event_trig_mat(all_interp_lfp(:,cur_ch),rand_inds,backlag,forwardlag);
bad_lines = find(any(isnan(rand_ev_mat),2));
rand_ev_mat(bad_lines,:) = [];
rcov = cov(rand_ev_mat);
[U,S,V] = svd(trig_cov - rcov);

SCORE = ev_mat*U;
nSCORE = zscore(SCORE);
COEFF = U;

% [COEFF, SCORE, LATENT] = princomp(ev_mat);
% SCORE = zscore(SCORE);
% ov_mean = mean(ev_mat);

backlag = round(0/desired_dt);
forwardlag = round(0.15/desired_dt);




%%
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);

x0_avg = 0.35;
y0_avg = -0.4;
cd ~/Data/bruce/Expt_1_8_13_imfolder

new_RF_patch = [-2 3; -3 2]; %location of RFs in degrees [x1 x2;y1 y2]
new_RF_patch2 = [0 1; -1 0]; %location of RFs in degrees [x1 x2;y1 y2]
% new_RF_patch = [-1 2; -2 1]; %location of RFs in degrees [x1 x2;y1 y2]

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));
xpatch_inds_cropped2 = find(xax >= new_RF_patch2(1,1) & xax <= new_RF_patch2(1,2));
ypatch_inds_cropped2 = find(yax >= new_RF_patch2(2,1) & yax <= new_RF_patch2(2,2));

surr_dsf = 5;
per_dim = ceil(length(xpatch_inds_cropped)/surr_dsf);
per_dim2 = length(xpatch_inds_cropped2);

new_crop = find(XAX >= new_RF_patch(1,1) & XAX <= new_RF_patch(1,2) & ...
    YAX >= new_RF_patch(2,1) & YAX <= new_RF_patch(2,2));
new_crop2 = find(XAX >= new_RF_patch2(1,1) & XAX <= new_RF_patch2(1,2) & ...
    YAX >= new_RF_patch2(2,1) & YAX <= new_RF_patch2(2,2));

[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XXc2,YYc2] = meshgrid(xax(xpatch_inds_cropped2),yax(ypatch_inds_cropped2));

cent_gauss_params = [x0_avg y0_avg 0 0 0.25 1];
surr_gauss_params = [x0_avg y0_avg 0 0 1.25 1];

% fft_cout = nan(length(used_stim_inds),per_dim2,per_dim2);
fft_sout = nan(length(used_stim_inds),per_dim,per_dim);
for i = 1:length(used_stim_inds)
    i
    cur_rep = all_expt_stimid(used_stim_inds(i));
%     if cur_rep < 1e4
%         cur_fname = sprintf('IM1%.6d.png',cur_rep);
%     elseif cur_rep < 1e5
%         cur_fname = sprintf('IM1%.5d.png',cur_rep);
%     else 
%         cur_fname = sprintf('IM1%.4d.png',cur_rep);        
%     end
cur_fname = sprintf('IM1%.6d.png',cur_rep);
cur_im = imread(cur_fname);

% cur_im_patch = reshape(double(cur_im(new_crop2)),length(ypatch_inds_cropped2),length(xpatch_inds_cropped2));
% cent_gauss_mask = get_pgauss_mask_v2(XXc2,YYc2,cent_gauss_params);
% cent_im = cur_im_patch.*cent_gauss_mask;

cur_im_patch = reshape(double(cur_im(new_crop)),length(ypatch_inds_cropped),length(xpatch_inds_cropped));
surr_gauss_mask = get_pgauss_mask_v2(XXc,YYc,surr_gauss_params);
surr_im = cur_im_patch.*surr_gauss_mask;

% cent_gauss_mask = get_pgauss_mask_v2(XXc,YYc,cent_gauss_params);
% mask_diff = surr_gauss_mask - cent_gauss_mask;
% mask_diff(mask_diff < 0) = 0;
% mask_diff = mask_diff/sum(mask_diff(:));

% fft_cent = abs(fftshift(fft2(cent_im)));
fft_surr = abs(fftshift(fft2(surr_im)));
fft_surr = smoothn(fft_surr,3);
fft_surr = imresize(fft_surr,1/surr_dsf);
%     fft_cout(i,:,:) = fft_cent;
    fft_sout(i,:,:) = fft_surr;
    
    %     ff = linspace(-Fs/2,Fs/2,length(xpatch_inds_cropped));
   
    
end
%%
ff = linspace(-Fs/2,Fs/2,per_dim);
ff2 = linspace(-Fs/2,Fs/2,per_dim2);

close all
pred_mat = reshape(fft_sout,length(used_stim_inds),per_dim^2);
pred_mat = zscore(pred_mat);
% pred_mat2 = reshape(fft_cout,length(used_stim_inds),per_dim2^2);
% pred_mat2 = zscore(pred_mat2);

% pred_mat2 = pred_mat2 - pred_mat;

% sparse_lambda = 4;
% smooth_lambda = 4000;
sparse_lambda = 1;
smooth_lambda = 2000;
init_beta1 = 0.001*randn(per_dim^2,1);
init_beta2 = 0.001*randn(per_dim2^2,1);
% init_beta = [init_beta1; init_beta2];
init_beta = [init_beta1];
kern_inds = [ones(per_dim^2,1); 2*ones(per_dim2^2,1)];

clear beta
n_used_pcs = 10;
for i = 1:n_used_pcs
    [beta(i,:),best_LL(i)] = smooth_regress_2d(nSCORE(:,i),pred_mat,init_beta,smooth_lambda*per_dim,sparse_lambda*per_dim);
%     [beta2(i,:),best_LL2(i)] = smooth_regress_2d(SCORE(:,i),pred_mat2,init_beta2,smooth_lambda*per_dim2,sparse_lambda*per_dim2);
%     % [beta,best_LL] = smooth_regress_2d(zscore(peak_vals),pred_mat,init_beta,smooth_lambda,sparse_lambda);
%     [beta(i,:),best_LL(i)] = smooth_regress_2d_doubkern(SCORE(:,i),[pred_mat pred_mat2],init_beta,smooth_lambda*per_dim,sparse_lambda*per_dim,kern_inds);
    
%     subplot(2,1,1)
    imagesc(ff,ff,reshape(beta(i,kern_inds==1),per_dim,per_dim));
    cm = caxis();
    cm = max(abs(cm));
    caxis([-cm cm])
%     subplot(2,1,2)
%     imagesc(ff2,ff2,reshape(beta(i,kern_inds==2),per_dim2,per_dim2));
%     cm = caxis();
%     cm = max(abs(cm));
%     caxis([-cm cm])
    pause
    clf
end

%%
LATENT = sqrt(diag(S));
coef_pred = pred_mat*beta';
coef_pred = bsxfun(@times,coef_pred,sqrt(LATENT(1:n_used_pcs))');
model_pred = coef_pred*COEFF(:,1:n_used_pcs)';
model_pred = bsxfun(@plus,model_pred,ov_mean);

ev_mat_ms = bsxfun(@minus,ev_mat,mean(ev_mat));
subplot(2,1,1)
imagesc(ev_mat_ms);colorbar
subplot(2,1,2)
imagesc(model_pred);colorbar

% %%
% wavelet = 'haar';
% % Define wavelet decomposition level
% level = 5;
% % Compute multilevel 2D wavelet decomposition
% [C S] = wavedec2(cur_im_patch,level,wavelet);
% % Define colormap and set rescale value
% map = (gray(256)); 
% rv = length(map);
% % Plot wavelet decomposition using square mode
% plotwavelet2(C,S,level,wavelet,rv,'tree');
% title(['Decomposition at level ',num2str(level)]);
% colormap(gray)