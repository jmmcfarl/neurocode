clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];

% sim_sac_blocks = [14 18 24 28 37 40 50] - 6;
% repeat_inds = [1e3 2e3 3e3 4.01e5 4.02e5 4.03e5];


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

%%
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 100;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
[filt_gb,filt_ga] = butter(2,[30 50]/niqf);
[filt_blf,filt_alf] = butter(2,[1 25]/niqf);
% [filt_blf2,filt_alf2] = butter(2,[2.5 25]/niqf);
gam_smooth = round(0.01*Fsd);

use_lfps = [1:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

scales = logspace(log10(8),log10(75),20);
scales = [scales 85 100 115 130 165];
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%%
% for ss = 1:length(repeat_inds)
%     all_block_stim_lfps{ss} = [];
%     all_block_stim_lfpsf{ss} = [];
%     all_block_stim_phase{ss} = [];
%     all_block_stim_amp{ss} = [];
% end
all_block_stim_lfpslf = [];
all_block_stim_lfpslf2 = [];
all_block_stim_lfpsf = [];
all_block_stim_lfpsg = [];
% all_block_stim_phase = [];
% all_block_stim_amp = [];
all_block_seof = [];
all_block_time = [];
all_block_tsb = [];
all_block_trial = [];
for bb = 1:length(sim_sac_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
    Vmatf = [];
    Vmatlf = [];
    Vmatlf2 = [];
    Vmatg = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        dVf = [];
        dVg = [];
        dVlf = [];
%         dVlf2 = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dec_v  = decimate(V(cur_range),dsf);
            dV = [dV dec_v];
            dVlf = [dVlf filtfilt(filt_blf,filt_alf,dec_v)]; %do some high-pass filtering
%             dVlf2 = [dVlf2 filtfilt(filt_blf2,filt_alf2,dec_v)]; %do some high-pass filtering
            dVf = [dVf filtfilt(filt_b,filt_a,dec_v)]; %do some high-pass filtering
            cur_gamma = filtfilt(filt_gb,filt_ga,dec_v);
            cur_gamma = sqrt(jmm_smooth_1d_cor(cur_gamma.^2,gam_smooth));
            dVg = [dVg cur_gamma];
            %             dVg = [dVg filtfilt(filt_gb,filt_ga,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(:,ll) = dV;
        Vmatf(:,ll) = dVf;
        Vmatlf(:,ll) = dVlf;
%         Vmatlf2(:,ll) = dVlf2;
        Vmatg(:,ll) = dVg;
%         temp = cwt(dV,scales,'cmor1-1');
%         phasegrams(:,:,ll) = angle(temp)';
%         ampgrams(:,:,ll) = abs(temp)';
    end
    %     Vmatg = abs(hilbert(Vmatg));
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    cur_use_trials = find(all_expt_id==sim_sac_blocks(bb) & all_trial_completed==1  & ismember(all_trial_seof,repeat_inds));
    use_trial_start_inds = round(interp1(t_ax,1:length(t_ax),all_trial_start_times(cur_use_trials)));
    use_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),all_trial_stop_times(cur_use_trials)));
    bad_t = find(isnan(use_trial_start_inds) | isnan(use_trial_stop_inds));
    cur_use_trials(bad_t) = [];
    use_trial_start_inds(bad_t) = [];
    use_trial_stop_inds(bad_t) = [];
    
    for tt = 1:length(cur_use_trials)
        cur_set = use_trial_start_inds(tt):use_trial_stop_inds(tt);
        all_block_seof = [all_block_seof; ones(length(cur_set),1)*all_trial_seof(cur_use_trials(tt))];
        all_block_trial = [all_block_trial; ones(length(cur_set),1)*cur_use_trials(tt)];
        all_block_time = [all_block_time; t_ax(cur_set)'];
        all_block_tsb = [all_block_tsb; t_ax(cur_set)'-t_ax(cur_set(1))];
        all_block_stim_lfpslf = [all_block_stim_lfpslf; Vmatlf(cur_set,:)];
%         all_block_stim_lfpslf2 = [all_block_stim_lfpslf2; Vmatlf2(cur_set,:)];
        all_block_stim_lfpsf = [all_block_stim_lfpsf; Vmatf(cur_set,:)];
        all_block_stim_lfpsg = [all_block_stim_lfpsg; Vmatg(cur_set,:)];
%         all_block_stim_phase = cat(1,all_block_stim_phase, phasegrams(cur_set,:,:));
%         all_block_stim_amp = cat(1,all_block_stim_amp,ampgrams(cur_set,:,:));
    end
    
end
%% COMPUTE STIM INFO AND WHITENED LFP AVGS
stim_fs = 1e4/117.5;
stim_inds = 1:40:320;
stim_times = stim_inds/stim_fs;
trial_dur = round(320/stim_fs*Fsd);
stim_dur = trial_dur/8;
lags = (1:trial_dur)/Fsd;

% amp_stds = shiftdim(std(all_block_stim_amp),-1);

% cur_rep_set = [7 8 9];
% cur_rep_set = [1 2 3];
cur_rep_set = [1 2 3 7 8 9];
% cur_rep_set = [4 5 6];
all_stim_cond_mean = [];
all_stim_cond_std = [];
all_stim_avg_out = [];
all_stim_sem_out = [];
all_mat = [];
for ss = 1:length(cur_rep_set)
    ss
    trial_set_start_inds = 1+find(all_block_seof(1:end-1) ~= repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) == repeat_inds(cur_rep_set(ss)));
    if all_block_seof(1) == repeat_inds(cur_rep_set(ss))
        trial_set_start_inds = [1; trial_set_start_inds];
    end
    trial_set_stop_inds = 1+find(all_block_seof(1:end-1) == repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) ~= repeat_inds(cur_rep_set(ss)));
    if all_block_seof(end) == repeat_inds(cur_rep_set(ss))
        trial_set_stop_inds = [trial_set_stop_inds; length(all_block_seof)];
    end
    
    cur_n_trials = length(trial_set_start_inds);
    cur_phase_mat = nan(cur_n_trials,trial_dur,length(wfreqs),length(use_lfps));
    cur_amp_mat = nan(cur_n_trials,trial_dur,length(wfreqs),length(use_lfps));
%     cur_raw_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
    %compile matrix
    for tt = 1:cur_n_trials
        cur_set = trial_set_start_inds(tt):trial_set_stop_inds(tt);
        cur_set(trial_dur+1:end) = [];
        cur_phase_mat(tt,1:length(cur_set),:,:) = all_block_stim_phase(cur_set,:,:);
        cur_amp_mat(tt,1:length(cur_set),:,:) = all_block_stim_amp(cur_set,:,:);
%         cur_raw_mat(tt,1:length(cur_set),:) = all_block_stim_lfpslf(cur_set,:);
    end
    
%     cur_out_mat = bsxfun(@rdivide,cur_amp_mat.*cos(cur_phase_mat),amp_stds);
    cur_out_mat = cur_amp_mat.*cos(cur_phase_mat);
    cur_out_mat(:,1:stim_dur,:,:) = []; %throw out first stim
    
    cur_out_vec = squeeze(nansum(cur_out_mat,3));
    stim_avg_out = squeeze(nanmean(cur_out_vec));
    all_stim_avg_out = cat(1,all_stim_avg_out,stim_avg_out);
    stim_sem_out = squeeze(nanstd(cur_out_vec))/sqrt(cur_n_trials);
    all_stim_sem_out = cat(1,all_stim_sem_out,stim_sem_out);
    
    r_out_mat = permute(cur_out_mat,[2 1 3 4]);
    r_out_mat1 = reshape(r_out_mat,[stim_dur,7,cur_n_trials,length(wfreqs),length(use_lfps)]);
    stim_cond_means = squeeze(nanmean(r_out_mat1,3));
    stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
    all_stim_cond_mean = cat(2,all_stim_cond_mean,stim_cond_means);
    all_stim_cond_std = cat(2,all_stim_cond_std,stim_cond_std);
    
    r_out_mat2 = reshape(r_out_mat,[stim_dur,7*cur_n_trials,length(wfreqs),length(use_lfps)]);
    all_mat = cat(2,all_mat,r_out_mat2);
    
    used_n_trials(ss) = cur_n_trials;
end
ov_marg_mat = reshape(all_mat,[size(all_mat,1)*size(all_mat,2) length(wfreqs) length(use_lfps)]);
ov_marg_means = squeeze(nanmean(ov_marg_mat));
ov_marg_stds = squeeze(nanstd(ov_marg_mat));

sac_marg_means = squeeze(nanmean(all_mat,2));
sac_marg_stds = squeeze(nanstd(all_mat,[],2));
sac_ovmean = repmat(shiftdim(ov_marg_means,-1),[size(sac_marg_means,1) 1 1]);
sac_ovstd = repmat(shiftdim(ov_marg_stds,-1),[size(sac_marg_stds,1) 1 1]);
sac_kl_div = 1./(2*sac_ovstd.^2).*((sac_marg_means - sac_ovmean).^2 + sac_marg_stds.^2 - sac_ovstd.^2) + log(sac_ovstd./sac_marg_stds);

stim_sacmean = repmat(permute(sac_marg_means,[1 4 2 3]),[1 size(all_stim_cond_mean,2) 1 1]);
stim_sacstd = repmat(permute(sac_marg_stds,[1 4 2 3]),[1 size(all_stim_cond_mean,2) 1 1]);
stim_kl_div = 1./(2*stim_sacstd.^2).*((all_stim_cond_mean - stim_sacmean).^2 + ...
    all_stim_cond_std.^2 - stim_sacstd.^2) + log(stim_sacstd./all_stim_cond_std);
avg_stim_kl_div = squeeze(mean(stim_kl_div,2));
elavg_stim_kl_div = squeeze(mean(avg_stim_kl_div,3));

stim_avg_out = reshape(all_stim_avg_out,[stim_dur,size(all_stim_avg_out,1)/stim_dur,length(use_lfps)]);
stim_sem_out = reshape(all_stim_sem_out,[stim_dur,size(all_stim_avg_out,1)/stim_dur,length(use_lfps)]);

%%
trial_seq_avgs = reshape(all_stim_avg_out,[stim_dur*7 6 length(use_lfps)]);
trial_seq_sem = reshape(all_stim_sem_out,[stim_dur*7 6 length(use_lfps)]);
trial_seq_ravgs = reshape(all_astim_avg,[stim_dur*7 6 length(use_lfps)]);
trial_seq_rsem = reshape(all_astim_sem,[stim_dur*7 6 length(use_lfps)]);
trial_seq_gavgs = reshape(all_gstim_avg,[stim_dur*7 6 length(use_lfps)]);
trial_seq_gsem = reshape(all_gstim_sem,[stim_dur*7 6 length(use_lfps)]);

% % close all
% figure
% ch = 20;
% ch2 = 50;
% st_n = 6;
% % shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_ravgs(:,st_n,ch)),squeeze(trial_seq_rsem(:,st_n,ch)),{'color','g'});
% hold on
% shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_ravgs(:,st_n,ch)),squeeze(trial_seq_rsem(:,st_n,ch2)),{'color','k'});
% % shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_avgs(:,st_n,ch)),squeeze(trial_seq_sem(:,st_n,ch)),{'color','g'});
% hold on
% % shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_avgs(:,st_n,ch)),squeeze(trial_seq_sem(:,st_n,ch2)),{'color','k'});
% 
% st_n = 5;
% % shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_avgs(:,st_n,ch)),squeeze(trial_seq_sem(:,st_n,ch)),{'color','r'});
% % shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_avgs(:,st_n,ch2)),squeeze(trial_seq_sem(:,st_n,ch2)),{'color','r'});
% shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_ravgs(:,st_n,ch)),squeeze(trial_seq_rsem(:,st_n,ch)),{'color','b'});
% shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_ravgs(:,st_n,ch2)),squeeze(trial_seq_rsem(:,st_n,ch2)),{'color','r'});
% yl = ylim();
% for i = 1:length(stim_times)
%     line(stim_times([i i]),yl,'color','k')
% end
% xlim([0 0.47*7])

figure
ch = 50;
ch2 = 15;
st_n = 6;
shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_gavgs(:,st_n,ch)),squeeze(trial_seq_gsem(:,st_n,ch)));
st_n = 5;
hold on
shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_gavgs(:,st_n,ch)),squeeze(trial_seq_gsem(:,st_n,ch)),{'color','r'});
% shadedErrorBar(lags(1:stim_dur*7),squeeze(trial_seq_gavgs(:,st_n,ch2)),squeeze(trial_seq_gsem(:,st_n,ch2)),{'color','b'});
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0 0.47*4])
%% compute null distribution of stim info
n_iter = 500;
sz_kl = size(avg_stim_kl_div);
null_stim_kl_div = nan([n_iter,sz_kl]);

for n = 1:n_iter
    fprintf('Iteration %d of %d\n',n,n_iter);
    
    rand_trial_ord = randperm(size(all_mat,2));
    
    all_stim_cond_mean = [];
    all_stim_cond_std = [];
    cur_tind = 0;
    for ss = 1:length(cur_rep_set)
        cur_tind = cur_tind(end) + (1:used_n_trials(ss));
        cur_srange = ((cur_tind(1)-1)*7+1):cur_tind(end)*7;
        r_out_mat1 = reshape(all_mat(:,rand_trial_ord(cur_srange),:,:),[stim_dur 7 used_n_trials(ss) length(wfreqs) length(use_lfps)]);
        
        stim_cond_means = squeeze(nanmean(r_out_mat1,3));
        stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
        all_stim_cond_mean = cat(2,all_stim_cond_mean,stim_cond_means);
        all_stim_cond_std = cat(2,all_stim_cond_std,stim_cond_std);
    end
    
%     stim_sacmean = repmat(permute(sac_marg_means,[1 4 2 3]),[1 size(all_stim_cond_mean,2) 1 1]);
%     stim_sacstd = repmat(permute(sac_marg_stds,[1 4 2 3]),[1 size(all_stim_cond_mean,2) 1 1]);
    stim_kl_div = 1./(2*stim_sacstd.^2).*((all_stim_cond_mean - stim_sacmean).^2 + ...
        all_stim_cond_std.^2 - stim_sacstd.^2) + log(stim_sacstd./all_stim_cond_std);
    null_stim_kl_div(n,:,:,:) = squeeze(mean(stim_kl_div,2));
end
null_stim_avg_kl_div = squeeze(nanmean(null_stim_kl_div));
null_stim_std_kl_div = squeeze(nanstd(null_stim_kl_div));

avg_stim_kl_div_z = (avg_stim_kl_div-null_stim_avg_kl_div)./null_stim_std_kl_div;
avg_stim_kl_div_sig = avg_stim_kl_div;
% avg_stim_kl_div_sig(avg_stim_kl_div_z < 3.09) = nan;
avg_stim_kl_div_sig(avg_stim_kl_div_z < 2.3) = nan;
elavg_stim_kl_div_z = squeeze(mean(avg_stim_kl_div_z,3));
elavg_stim_kl_div_sig = squeeze(mean(avg_stim_kl_div,3));
elavg_stim_kl_div_sig(squeeze(sum(avg_stim_kl_div_z > 3.09,3)) < 96/2) = nan;
sigline = nan(141,1);
for i = 1:141
    [~,sigline(i)] = find(isnan(elavg_stim_kl_div_sig(i,:)),1,'last');
end
%%
ex_el = 6;
figure
pcolor(lags(1:stim_dur),wfreqs,squeeze(avg_stim_kl_div_sig(:,:,ex_el))'/log(2));shading interp
set(gca,'yscale','log')
colorbar

figure
% pcolor(lags(1:stim_dur),wfreqs,elavg_stim_kl_div_z'/log(2));shading interp
pcolor(lags(1:stim_dur),wfreqs,log10(elavg_stim_kl_div'/log(2)));shading interp
% pcolor(lags(1:stim_dur),wfreqs,(elavg_stim_kl_div_sig'/log(2)));shading interp
hold on
plot(lags(1:stim_dur),wfreqs(sigline),'r')
set(gca,'yscale','log')
% caxis([0 20]);colorbar
%% COMPUTE STIM INFO FOR GAMMA AMP and raw AMP
% all_block_stim_lfpsg_norm = zscore(all_block_stim_lfpsg);
% all_block_stim_lfpsg_norm = zscore(sqrt(all_block_stim_lfpsg));
all_block_stim_lfpsg_norm = nan(size(all_block_stim_lfpsg));
for i = 1:96
    all_block_stim_lfpsg_norm(:,i) = zscore(boxcox(all_block_stim_lfpsg(:,i)));
end

cur_rep_set = [1 2 3 7 8 9];
% cur_rep_set = [4 5 6];
all_gstim_cond_mean = [];
all_gstim_cond_std = [];
all_astim_cond_mean = [];
all_astim_cond_std = [];
all_astim_avg = [];
all_astim_sem = [];
all_gstim_avg = [];
all_gstim_sem = [];
all_gmat = [];
all_amat = [];
for ss = 1:length(cur_rep_set)
    ss
    trial_set_start_inds = 1+find(all_block_seof(1:end-1) ~= repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) == repeat_inds(cur_rep_set(ss)));
    if all_block_seof(1) == repeat_inds(cur_rep_set(ss))
        trial_set_start_inds = [1; trial_set_start_inds];
    end
    trial_set_stop_inds = 1+find(all_block_seof(1:end-1) == repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) ~= repeat_inds(cur_rep_set(ss)));
    if all_block_seof(end) == repeat_inds(cur_rep_set(ss))
        trial_set_stop_inds = [trial_set_stop_inds; length(all_block_seof)];
    end
    
    cur_n_trials = length(trial_set_start_inds);
    cur_lfpg_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
    cur_lfp_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
    %compile matrix
    for tt = 1:cur_n_trials
        cur_set = trial_set_start_inds(tt):trial_set_stop_inds(tt);
        cur_set(trial_dur+1:end) = [];
        cur_lfpg_mat(tt,1:length(cur_set),:) = all_block_stim_lfpsg_norm(cur_set,:);
        cur_lfp_mat(tt,1:length(cur_set),:) = all_block_stim_lfpslf(cur_set,:);
    end
    
    cur_lfpg_mat(:,1:stim_dur,:) = []; %throw out first stim of each trial
    cur_lfp_mat(:,1:stim_dur,:) = []; %throw out first stim of each trial
    
    cur_avg = squeeze(nanmean(cur_lfp_mat));
    cur_sem = squeeze(nanstd(cur_lfp_mat))/sqrt(cur_n_trials);
    all_astim_avg = cat(1,all_astim_avg,cur_avg);
    all_astim_sem = cat(1,all_astim_sem,cur_sem);
    cur_avg = squeeze(nanmean(cur_lfpg_mat));
    cur_sem = squeeze(nanstd(cur_lfpg_mat))/sqrt(cur_n_trials);
    all_gstim_avg = cat(1,all_gstim_avg,cur_avg);
    all_gstim_sem = cat(1,all_gstim_sem,cur_sem);
    
    r_out_mat = permute(cur_lfpg_mat,[2 1 3]);
    r_out_mat1 = reshape(r_out_mat,[stim_dur,7,cur_n_trials,length(use_lfps)]);
    stim_cond_means = squeeze(nanmean(r_out_mat1,3));
    stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
    all_gstim_cond_mean = cat(2,all_gstim_cond_mean,stim_cond_means);
    all_gstim_cond_std = cat(2,all_gstim_cond_std,stim_cond_std);
    
    r_out_mat2 = reshape(r_out_mat,[stim_dur,7*cur_n_trials,length(use_lfps)]);
    all_gmat = cat(2,all_gmat,r_out_mat2);

    r_out_mat = permute(cur_lfp_mat,[2 1 3]);
    r_out_mat1 = reshape(r_out_mat,[stim_dur,7,cur_n_trials,length(use_lfps)]);
    stim_cond_means = squeeze(nanmean(r_out_mat1,3));
    stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
    all_astim_cond_mean = cat(2,all_astim_cond_mean,stim_cond_means);
    all_astim_cond_std = cat(2,all_astim_cond_std,stim_cond_std);
    
    r_out_mat2 = reshape(r_out_mat,[stim_dur,7*cur_n_trials,length(use_lfps)]);
    all_amat = cat(2,all_amat,r_out_mat2);
    
    used_n_trials(ss) = cur_n_trials;
end

ov_gmarg_mat = reshape(all_gmat,[size(all_gmat,1)*size(all_gmat,2) length(use_lfps)]);
ov_gmarg_means = squeeze(nanmean(ov_gmarg_mat));
ov_gmarg_stds = squeeze(nanstd(ov_gmarg_mat));
ov_amarg_mat = reshape(all_gmat,[size(all_amat,1)*size(all_amat,2) length(use_lfps)]);
ov_amarg_means = squeeze(nanmean(ov_amarg_mat));
ov_amarg_stds = squeeze(nanstd(ov_amarg_mat));

sac_gmarg_means = squeeze(nanmean(all_gmat,2));
sac_gmarg_stds = squeeze(nanstd(all_gmat,[],2));
sac_govmean = repmat(ov_gmarg_means,[size(sac_gmarg_means,1) 1]);
sac_govstd = repmat(ov_gmarg_stds,[size(sac_gmarg_stds,1) 1]);
sac_gkl_div = 1./(2*sac_govstd.^2).*((sac_gmarg_means - sac_govmean).^2 + sac_gmarg_stds.^2 - sac_govstd.^2) + log(sac_govstd./sac_gmarg_stds);

sac_amarg_means = squeeze(nanmean(all_amat,2));
sac_amarg_stds = squeeze(nanstd(all_amat,[],2));
sac_aovmean = repmat(ov_amarg_means,[size(sac_amarg_means,1) 1]);
sac_aovstd = repmat(ov_amarg_stds,[size(sac_amarg_stds,1) 1]);
sac_akl_div = 1./(2*sac_aovstd.^2).*((sac_amarg_means - sac_aovmean).^2 + sac_amarg_stds.^2 - sac_aovstd.^2) + log(sac_aovstd./sac_amarg_stds);

stim_gsacmean = repmat(permute(sac_gmarg_means,[1 3 2]),[1 size(all_gstim_cond_mean,2) 1]);
stim_gsacstd = repmat(permute(sac_gmarg_stds,[1 3 2]),[1 size(all_gstim_cond_mean,2) 1]);
gstim_kl_div = 1./(2*stim_gsacstd.^2).*((all_gstim_cond_mean - stim_gsacmean).^2 + ...
    all_gstim_cond_std.^2 - stim_gsacstd.^2) + log(stim_gsacstd./all_gstim_cond_std);
avg_gstim_kl_div = squeeze(mean(gstim_kl_div,2));
elavg_gstim_kl_div = squeeze(mean(avg_gstim_kl_div,2));
gamma_stim_avgs = all_gstim_cond_mean;

stim_asacmean = repmat(permute(sac_amarg_means,[1 3 2]),[1 size(all_astim_cond_mean,2) 1]);
stim_asacstd = repmat(permute(sac_amarg_stds,[1 3 2]),[1 size(all_astim_cond_mean,2) 1]);
astim_kl_div = 1./(2*stim_asacstd.^2).*((all_astim_cond_mean - stim_asacmean).^2 + ...
    all_astim_cond_std.^2 - stim_asacstd.^2) + log(stim_asacstd./all_astim_cond_std);
avg_astim_kl_div = squeeze(mean(astim_kl_div,2));
elavg_astim_kl_div = squeeze(mean(avg_astim_kl_div,2));
amp_stim_avgs = all_astim_cond_mean;

astim_avg_out = reshape(all_astim_avg,[stim_dur,size(all_astim_avg,1)/stim_dur,length(use_lfps)]);
astim_sem_out = reshape(all_astim_sem,[stim_dur,size(all_astim_avg,1)/stim_dur,length(use_lfps)]);

gstim_avg_out = reshape(all_gstim_avg,[stim_dur,size(all_gstim_avg,1)/stim_dur,length(use_lfps)]);
gstim_sem_out = reshape(all_gstim_sem,[stim_dur,size(all_gstim_avg,1)/stim_dur,length(use_lfps)]);

% figure
% shadedErrorBar(lags(1:stim_dur),nanmean(avg_gstim_kl_div,2),nanstd(avg_gstim_kl_div,[],2)/sqrt(96));
% hold on
% shadedErrorBar(lags(1:stim_dur),nanmean(avg_astim_kl_div,2),nanstd(avg_astim_kl_div,[],2)/sqrt(96),{'color','r'});
% xlim([0 0.47])


%% compute null distribution of stim info for gamma and amp
n_iter = 500;
sz_kl = size(avg_gstim_kl_div);
gnull_stim_kl_div = nan([n_iter,sz_kl]);
anull_stim_kl_div = nan([n_iter,sz_kl]);

for n = 1:n_iter
    fprintf('Iteration %d of %d\n',n,n_iter);
    
    rand_trial_ord = randperm(size(all_gmat,2));
    
    all_gstim_cond_mean = [];
    all_gstim_cond_std = [];
    all_astim_cond_mean = [];
    all_astim_cond_std = [];
    cur_tind = 0;
    for ss = 1:length(cur_rep_set)
        cur_tind = cur_tind(end) + (1:used_n_trials(ss));
        cur_srange = ((cur_tind(1)-1)*7+1):cur_tind(end)*7;
        r_out_mat1 = reshape(all_gmat(:,rand_trial_ord(cur_srange),:,:),[stim_dur 7 used_n_trials(ss) length(use_lfps)]);
        stim_cond_means = squeeze(nanmean(r_out_mat1,3));
        stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
        all_gstim_cond_mean = cat(2,all_gstim_cond_mean,stim_cond_means);
        all_gstim_cond_std = cat(2,all_gstim_cond_std,stim_cond_std);
        
        r_out_mat1 = reshape(all_amat(:,rand_trial_ord(cur_srange),:,:),[stim_dur 7 used_n_trials(ss) length(use_lfps)]);
        stim_cond_means = squeeze(nanmean(r_out_mat1,3));
        stim_cond_std = squeeze(nanstd(r_out_mat1,[],3));
        all_astim_cond_mean = cat(2,all_astim_cond_mean,stim_cond_means);
        all_astim_cond_std = cat(2,all_astim_cond_std,stim_cond_std);
    end
    
    stim_gkl_div = 1./(2*stim_gsacstd.^2).*((all_gstim_cond_mean - stim_gsacmean).^2 + ...
        all_gstim_cond_std.^2 - stim_gsacstd.^2) + log(stim_gsacstd./all_gstim_cond_std);
    gnull_stim_kl_div(n,:,:) = squeeze(mean(stim_gkl_div,2));
    stim_akl_div = 1./(2*stim_asacstd.^2).*((all_astim_cond_mean - stim_asacmean).^2 + ...
        all_astim_cond_std.^2 - stim_asacstd.^2) + log(stim_asacstd./all_astim_cond_std);
    anull_stim_kl_div(n,:,:) = squeeze(mean(stim_akl_div,2));
end
gnull_stim_avg_kl_div = squeeze(nanmean(gnull_stim_kl_div));
gnull_stim_std_kl_div = squeeze(nanstd(gnull_stim_kl_div));
anull_stim_avg_kl_div = squeeze(nanmean(anull_stim_kl_div));
anull_stim_std_kl_div = squeeze(nanstd(anull_stim_kl_div));

avg_gstim_kl_div_z = (avg_gstim_kl_div-gnull_stim_avg_kl_div)./gnull_stim_std_kl_div;
avg_gstim_kl_div_sig = avg_gstim_kl_div;
avg_gstim_kl_div_sig(avg_gstim_kl_div_z < 2.3) = nan;
avg_astim_kl_div_z = (avg_astim_kl_div-anull_stim_avg_kl_div)./anull_stim_std_kl_div;
avg_astim_kl_div_sig = avg_astim_kl_div;
avg_astim_kl_div_sig(avg_astim_kl_div_z < 2.3) = nan;

%%
figure
shadedErrorBar(lags(1:stim_dur),nanmean(avg_gstim_kl_div_z,2),nanstd(avg_gstim_kl_div_z,[],2)/sqrt(96));
xlim([0 0.47])
ylim([0 7])
figure
hold on
shadedErrorBar(lags(1:stim_dur),nanmean(avg_astim_kl_div_z,2),nanstd(avg_astim_kl_div_z,[],2)/sqrt(96),{'color','r'});
xlim([0 0.47])

figure
pcolor(lags(1:stim_dur),1:96,avg_gstim_kl_div_sig'/log(2));
shading flat;caxis([0.1 0.25]);colorbar

figure
pcolor(lags(1:stim_dur),1:96,avg_astim_kl_div_sig'/log(2));
shading flat;caxis([0.1 0.45]);colorbar

%%
stim_fs = 1e4/117.5;
stim_inds = 1:40:320;
stim_times = stim_inds/stim_fs;
trial_dur = round(320/stim_fs*Fsd);
stim_dur = trial_dur/8;

all_block_stim_lfpsg_norm = zscore(all_block_stim_lfpsg);
all_block_stim_lfps_norm = zscore(all_block_stim_lfpslf);

cur_rep_set = [7 8 9];
all_phase_mat = [];
all_phase_cons = [];
% for ss = 1:length(cur_rep_set)
ss = 2;
trial_set_start_inds = 1+find(all_block_seof(1:end-1) ~= repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) == repeat_inds(cur_rep_set(ss)));
if all_block_seof(1) == repeat_inds(cur_rep_set(ss))
    trial_set_start_inds = [1; trial_set_start_inds];
end
trial_set_stop_inds = 1+find(all_block_seof(1:end-1) == repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) ~= repeat_inds(cur_rep_set(ss)));
if all_block_seof(end) == repeat_inds(cur_rep_set(ss))
    trial_set_stop_inds = [trial_set_stop_inds; length(all_block_seof)];
end
% trial_set_start_inds = 1+find(rand_seof(1:end-1) ~= repeat_inds(cur_rep_set(ss)) & rand_seof(2:end) == repeat_inds(cur_rep_set(ss)));
% if rand_seof(1) == repeat_inds(cur_rep_set(ss))
%     trial_set_start_inds = [1; trial_set_start_inds];
% end
% trial_set_stop_inds = 1+find(rand_seof(1:end-1) == repeat_inds(cur_rep_set(ss)) & rand_seof(2:end) ~= repeat_inds(cur_rep_set(ss)));
% if rand_seof(end) == repeat_inds(cur_rep_set(ss))
%     trial_set_stop_inds = [trial_set_stop_inds; length(all_block_seof)];
% end

cur_n_trials = length(trial_set_start_inds);
cur_lfp_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
cur_lfpg_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
%compile matrix
for tt = 1:cur_n_trials
    cur_set = trial_set_start_inds(tt):trial_set_stop_inds(tt);
    cur_set(trial_dur+1:end) = [];
    cur_lfp_mat(tt,1:length(cur_set),:) = all_block_stim_lfps_norm(cur_set,:);
    cur_lfpg_mat(tt,1:length(cur_set),:) = all_block_stim_lfpsg_norm(cur_set,:);
end
stim_lfp_avg = squeeze(nanmean(cur_lfp_mat));
stim_lfp_sem = squeeze(bsxfun(@rdivide,nanstd(cur_lfp_mat),sqrt(sum(~isnan(cur_lfp_mat)))));
stim_lfpg_avg = squeeze(nanmean(cur_lfpg_mat));
stim_lfpg_sem = squeeze(bsxfun(@rdivide,nanstd(cur_lfpg_mat),sqrt(sum(~isnan(cur_lfpg_mat)))));

istim_lfp_avg = reshape(stim_lfp_avg,stim_dur,8,length(use_lfps));

cur_ch = 1;
shadedErrorBar(lags,stim_lfp_avg(:,cur_ch),stim_lfp_sem(:,cur_ch),{'color','r'})
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%% COMPUTE CORRELATION-VS DISTANCE PARAMETERS FOR STIMULUS-TRIG AVGS

gstim_avg_out_vec = reshape(gstim_avg_out,141*size(gstim_avg_out,2),96);
astim_avg_out_vec = reshape(astim_avg_out,141*size(gstim_avg_out,2),96);



bad_elecs = [16 67 92];
used_elecs = setdiff(use_lfps,bad_elecs);

el_dist = squareform(pdist([X_pos(use_lfps(used_elecs))' Y_pos(use_lfps(used_elecs))']));
el_dist = el_dist(:);
cset = find(el_dist > 0);

cur_corrmat = corr(gstim_avg_out_vec(:,used_elecs));
cur_corrmat = cur_corrmat(:);
cur_corrmat = cur_corrmat(cset);

cur_corrmat2 = corr(astim_avg_out_vec(:,used_elecs));
cur_corrmat2 = cur_corrmat2(:);
cur_corrmat2 = cur_corrmat2(cset);

xx = linspace(0.3,5,100);
figure
plot(el_dist(cset)*0.4,cur_corrmat2,'r.','markersize',6)
hold on
plot(el_dist(cset)*0.4,cur_corrmat,'b.','markersize',6);
% p = polyfit(el_dist(cset)*0.4,cur_corrmat2,1);
% pval = polyval(p,xx);
% plot(xx,pval,'k','linewidth',2)
% p = polyfit(el_dist(cset)*0.4,cur_corrmat,1);
% pval = polyval(p,xx);
% plot(xx,pval,'k','linewidth',2)
[~,ord] = sort(el_dist(cset));
sm_fun = smooth(cur_corrmat(ord),5000,'lowess');
plot(el_dist(cset(ord))*0.4,sm_fun,'k','linewidth',2)
[~,ord] = sort(el_dist(cset));
sm_fun2 = smooth(cur_corrmat2(ord),5000,'lowess');
plot(el_dist(cset(ord))*0.4,sm_fun2,'k','linewidth',2)
xlim([0 5.5])
beta_0 = [0.5 2];
LB = [0 1];
UB = [1 100];
avg_betas_corr2 = lsqnonlin(@(X) cur_corrmat2(:)-expon_decay_fun(X,el_dist(cset)*0.4),beta_0,LB,UB);
avg_betas_corr = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist(cset)*0.4),beta_0,LB,UB);


cur_corrmat = corr(all_ms_lfpg_vec(:,used_elecs));
cur_corrmat = cur_corrmat(:);
cur_corrmat = cur_corrmat(cset);

cur_corrmat2 = corr(all_ms_lfp_vec(:,used_elecs));
cur_corrmat2 = cur_corrmat2(:);
cur_corrmat2 = cur_corrmat2(cset);

xx = linspace(0.3,5,100);
figure
plot(el_dist(cset)*0.4,cur_corrmat2,'r.','markersize',6)
hold on
plot(el_dist(cset)*0.4,cur_corrmat,'b.','markersize',6);
% p = polyfit(el_dist(cset)*0.4,cur_corrmat2,1);
% pval = polyval(p,xx);
% plot(xx,pval,'k','linewidth',2)
% p = polyfit(el_dist(cset)*0.4,cur_corrmat,1);
% pval = polyval(p,xx);
% plot(xx,pval,'k','linewidth',2)
[~,ord] = sort(el_dist(cset));
sm_fun = smooth(cur_corrmat(ord),5000,'lowess');
plot(el_dist(cset(ord))*0.4,sm_fun,'k','linewidth',2)
[~,ord] = sort(el_dist(cset));
sm_fun2 = smooth(cur_corrmat2(ord),5000,'lowess');
plot(el_dist(cset(ord))*0.4,sm_fun2,'k','linewidth',2)
xlim([0 5.5])
spont_betas_corr2 = lsqnonlin(@(X) cur_corrmat2(:)-expon_decay_fun(X,el_dist(cset)*0.4),beta_0,LB,UB);
spont_betas_corr = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist(cset)*0.4),beta_0,LB,UB);

%%
all_ms_lfpg_vec = [];
all_ms_lfp_vec = [];
% cur_rep_set = [1 2 3 7 8 9];
for ss = 1:length(cur_rep_set)
    ss
    trial_set_start_inds = 1+find(all_block_seof(1:end-1) ~= repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) == repeat_inds(cur_rep_set(ss)));
    if all_block_seof(1) == repeat_inds(cur_rep_set(ss))
        trial_set_start_inds = [1; trial_set_start_inds];
    end
    trial_set_stop_inds = 1+find(all_block_seof(1:end-1) == repeat_inds(cur_rep_set(ss)) & all_block_seof(2:end) ~= repeat_inds(cur_rep_set(ss)));
    if all_block_seof(end) == repeat_inds(cur_rep_set(ss))
        trial_set_stop_inds = [trial_set_stop_inds; length(all_block_seof)];
    end
    
    cur_n_trials = length(trial_set_start_inds);
    cur_lfpg_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
    cur_lfp_mat = nan(cur_n_trials,trial_dur,length(use_lfps));
    %compile matrix
    for tt = 1:cur_n_trials
        cur_set = trial_set_start_inds(tt):trial_set_stop_inds(tt);
        cur_set(trial_dur+1:end) = [];
        cur_lfpg_mat(tt,1:length(cur_set),:) = all_block_stim_lfpsg_norm(cur_set,:);
        cur_lfp_mat(tt,1:length(cur_set),:) = all_block_stim_lfpslf2(cur_set,:);
    end
    
    cur_lfpg_mat(:,1:stim_dur,:) = []; %throw out first stim of each trial
    cur_lfp_mat(:,1:stim_dur,:) = []; %throw out first stim of each trial

    cur_ms_lfpg_mat = bsxfun(@minus,cur_lfpg_mat,nanmean(cur_lfpg_mat));
    cur_ms_lfp_mat = bsxfun(@minus,cur_lfp_mat,nanmean(cur_lfp_mat));
    cur_ms_lfpg_mat = permute(cur_ms_lfpg_mat,[2 1 3]);
    cur_ms_lfp_mat = permute(cur_ms_lfp_mat,[2 1 3]);
    
    all_ms_lfpg_vec = cat(1,all_ms_lfpg_vec,reshape(cur_ms_lfpg_mat,[cur_n_trials*size(cur_ms_lfpg_mat,1) 96]));
    all_ms_lfp_vec = cat(1,all_ms_lfp_vec,reshape(cur_ms_lfp_mat,[cur_n_trials*size(cur_ms_lfp_mat,1) 96]));
        
end

all_ms_lfpg_vec(isnan(all_ms_lfpg_vec(:,1)),:) = [];
all_ms_lfp_vec(isnan(all_ms_lfp_vec(:,1)),:) = [];

%%
betas_coh = nan(length(wfreqs),stim_dur);
avg_coh = nan(length(wfreqs),stim_dur);
for ww = 1:length(wfreqs)
    for tt = 1:stim_dur;
        cur_corrmat = corr(squeeze(all_stim_cond_mean(tt,:,ww,used_elecs)));
        cur_corrmat = cur_corrmat(:);
        cur_corrmat = cur_corrmat(cset);
%         betas_coh(ww,tt) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun2(X,el_dist),beta_0,LB,UB);
        betas_coh(ww,tt) = regress(1-cur_corrmat,[el_dist(cset)]);
        avg_coh(ww,tt) =mean(cur_corrmat);
        end
    ww
end

% cur_n_stims = size(stim_avg_out,2);
% betas_coh = nan(cur_n_stims,1);
% avg_coh = nan(cur_n_stims,1);
% all_corrs = [];
% all_dists = [];
% for ww = 1:cur_n_stims
%             cur_corrmat = corr(squeeze(stim_avg_out(:,ww,used_elecs)));
%             cur_corrmat = cur_corrmat(:);
%             cur_corrmat = cur_corrmat(cset);
%         all_corrs = [all_corrs; cur_corrmat];
%         all_dists = [all_dists; el_dist(cset)];
% %         betas_coh(ww,tt) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun2(X,el_dist),beta_0,LB,UB);
%         betas_coh(ww) = regress(1-cur_corrmat,[el_dist(cset)]);
%         avg_coh(ww) =mean(cur_corrmat);
%     ww
% end
% 
% figure;hold on
% use_dists = unique(all_dists);
% for ww = 1:cur_n_stims;
%    plot(use_dists,1-use_dists*betas_coh(ww) ,'r')
%     
% end
%%
% % close all
% figure
lags = (1:(trial_dur))/Fsd;
eps = 1.5;
[~,use_ordx] = sort(X_pos(use_lfps));
[~,use_ordy] = sort(Y_pos(use_lfps));
% for i = 1:size(stim_avg_out,2)
rep_set = [1 2 3 7 8 9];
    cd ~/Data/bruce/Expt_1_8_13_imfolder
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);
x0_avg = 0.35;
y0_avg = -0.4;
xw = 1;yw=1;
xw2 = 0.4;yw2=0.4;

% resh_avg_out = reshape(stim_avg_out,[141*7 6 96]);
resh_avg_out = reshape(all_astim_cond_mean,[141*7 6 96]);

f1 = figure(1);
f2 = figure(2);
stim_cnt = 1;
for ss = 1:length(rep_set)
    cur_rep = repeat_inds(rep_set(ss));
    
    figure(1);clf
    ut_cnt = 1;
    for i = 2:8
        figure(1)
        shadedErrorBar(lags(1:987),squeeze(mean(resh_avg_out(:,ss,:),3)),squeeze(std(resh_avg_out(:,ss,:),[],3)))
        hold on
        up_to = stim_dur*ut_cnt;
        plot(lags(1:up_to),squeeze(mean(resh_avg_out(1:up_to,ss,:),3)),'r','linewidth',2)
        ut_cnt = ut_cnt + 1;
        yl = ylim();
for j = 1:length(stim_times)
    line(stim_times([j j]),yl,'color','k')
end

% %     subplot(2,2,1)
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(stim_avg_out(:,stim_cnt,use_ordx))');shading flat
%     title('Avg LFPs x-sort')
% %     subplot(2,2,3)
% %     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(stim_avg_out(:,stim_cnt,use_ordy))');shading flat
% %     title('Avg LFPs y-sort')
% %     subplot(2,2,2)
% %     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(astim_avg_out(:,stim_cnt,use_ordx))');shading flat
% % %     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(gamma_stim_avgs(:,stim_cnt,use_ordx))');shading flat
% %     title('Avg Gamma x-sort')
% %     subplot(2,2,4)
% %     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(gamma_stim_avgs(:,stim_cnt,use_ordy))');shading flat
% %     title('Avg Gamma y-sort')

    stim_cnt = stim_cnt + 1;
    figure(2)
    if cur_rep < 1e4
        cur_fname = sprintf('IM1%.6d.png',cur_rep+i);
    elseif cur_rep < 1e5
        cur_fname = sprintf('IM1%.5d.png',cur_rep+i);
    else cur_rep < 1e5
        cur_fname = sprintf('IM1%.4d.png',cur_rep+i);
        
    end
    cur_fname
        cur_im = imread(cur_fname);
        imagesc(xax,yax,flipud(cur_im));colormap(gray);set(gca,'ydir','normal');
        hold on
%         rectangle('Position',[x0_avg-xw/2 y0_avg-yw/2 xw,yw],'edgecolor','r');
        rectangle('Position',[x0_avg-xw2/2 y0_avg-yw2/2 xw2,yw2],'edgecolor','r');
        xlim([-2 2]); ylim([-2 2])
                
        pause
    end   
end

%%
% ord = 2;
% cur_sig = all_block_stim_lfpslf(:,1);
% 
% acorr_seq2 = xcov(cur_sig,ord);
% acorr_seq2 = acorr_seq2(ord+1:end);
% A = levinson(acorr_seq2(1:ord),ord);
% est_x2 = filtfilt(-A,1,cur_sig);

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
xw = 1;yw=1;

cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
%fit smoothed retinotopic surface
orientations = linspace(0,pi-pi/12,12);
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

cd ~/Data/bruce/Expt_1_8_13_imfolder


inc = 1:15:stim_dur;


id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));

f1 = figure(1);
f2 = figure(2);
f3 = figure(3);
ss = 6;
cur_rep = repeat_inds(rep_set(ss));
i = 3;
stim_cnt = 7*(ss-1)+i;
for j = 1:length(inc)
    
    figure(1)
    subplot(2,1,1)
    pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(stim_avg_out(:,stim_cnt,use_ordx))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(astim_avg_out(:,stim_cnt,use_ordx))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(all_astim_cond_mean(:,stim_cnt,use_ordx))');shading flat
    title('Avg LFPs x-sort')
    yl = ylim();
    line(lags([inc(j) inc(j)]),yl)
    subplot(2,1,2)
    pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(astim_avg_out(:,stim_cnt,use_ordy))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(all_astim_cond_mean(:,stim_cnt,use_ordy))');shading flat
    title('Avg LFPs y-sort')
    line(lags([inc(j) inc(j)]),yl)
    figure(2)
    if cur_rep < 1e4
        cur_fname = sprintf('IM1%.6d.png',cur_rep+i);
    elseif cur_rep < 1e5
        cur_fname = sprintf('IM1%.5d.png',cur_rep+i);
    else cur_rep < 1e5
        cur_fname = sprintf('IM1%.4d.png',cur_rep+i);
    end
    cur_fname
    cur_im = imread(cur_fname);
    imagesc(xax,yax,flipud(cur_im));colormap(gray);set(gca,'ydir','normal');
    hold on
    rectangle('Position',[x0_avg-xw/2 y0_avg-yw/2 xw,yw],'edgecolor','r');
    xlim([-2 2]); ylim([-2 2])
    plot(interp_x,interp_y,'ro')
    
    figure(3)
    subplot(2,1,1)
    cur_set = squeeze(astim_avg_out(inc(j),stim_cnt,:));
%     cur_set = squeeze(all_astim_cond_mean(inc(j),stim_cnt,:));
    temp = nan(10,10);
    temp(use_ids) = cur_set(id_mat(use_ids));
    imagesc(temp);colorbar
    set(gca,'ydir','normal')
    subplot(2,1,2)
%     Vq = griddata(interp_x,interp_y,cur_set,Xi,Yi);
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
    Vq = F(Xi,Yi);
    pcolor(xi,yi,Vq); shading flat; colorbar;

    pause
    figure(1);clf;figure(3);clf;
end
%%
% close all
figure
lags = (1:(trial_dur))/Fsd;
eps = 1.5;
[~,use_ordx] = sort(X_pos(use_lfps));
[~,use_ordy] = sort(Y_pos(use_lfps));
for i = 1:size(stim_avg_out,2)
    subplot(2,1,1)
    pcolor(lags(1:stim_dur),1:48,squeeze(gamma_stim_avgs(:,i,use_ordx))');shading flat
    subplot(2,1,2)
    pcolor(lags(1:stim_dur),1:48,squeeze(gamma_stim_avgs(:,i,use_ordy))');shading flat
    pause
    clf
end

%%
save('all_data_set_locking_analysis','-v7.3')
save('all_data_set_locking_analysis_lf2','-v7.3','all_block_stim_lfpslf2')

%%
close all
cd ~/Data/bruce/Expt_1_8_13_imfolder
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
[XAX,YAX] = meshgrid(xax,yax);
x0_avg = 0.35;
y0_avg = -0.4;
xw = 1;yw=1;
xw2 = 0.4;yw2=0.4;

xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;

new_RF_patch = [-2 3; -3 2]; %location of RFs in degrees [x1 x2;y1 y2]
xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(XAX >= new_RF_patch(1,1) & XAX <= new_RF_patch(1,2) & ...
    YAX >= new_RF_patch(2,1) & YAX <= new_RF_patch(2,2));

[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));

lags = (1:(trial_dur))/Fsd;
[~,use_ordx] = sort(X_pos(use_lfps));
[~,use_ordy] = sort(Y_pos(use_lfps));

% for i = 1:size(stim_avg_out,2)
rep_set = [1 2 3 7 8 9];

% resh_avg_out = reshape(stim_avg_out,[141*7 6 96]);
resh_avg_out = reshape(all_astim_cond_mean,[141*7 6 96]);

cent_gauss_params = [x0_avg y0_avg 0 0 0.3 1];
surr_gauss_params = [x0_avg y0_avg 0 0 1.5 1];

f1 = figure(1);
f2 = figure(2);
stim_cnt = 1;
for ss = 1:length(rep_set)
    cur_rep = repeat_inds(rep_set(ss));
    
    figure(1);clf
    ut_cnt = 1;
    for i = 2:8
        figure(1)
        shadedErrorBar(lags(1:987),squeeze(mean(resh_avg_out(:,ss,:),3)),squeeze(std(resh_avg_out(:,ss,:),[],3)))
        hold on
        up_to = stim_dur*ut_cnt;
        plot(lags(1:up_to),squeeze(mean(resh_avg_out(1:up_to,ss,:),3)),'r','linewidth',2)
        ut_cnt = ut_cnt + 1;
        yl = ylim();
for j = 1:length(stim_times)
    line(stim_times([j j]),yl,'color','k')
end
    stim_cnt = stim_cnt + 1;
    figure(2)
    if cur_rep < 1e4
        cur_fname = sprintf('IM1%.6d.png',cur_rep+i);
    elseif cur_rep < 1e5
        cur_fname = sprintf('IM1%.5d.png',cur_rep+i);
    else cur_rep < 1e5
        cur_fname = sprintf('IM1%.4d.png',cur_rep+i);
        
    end
    cur_fname
        cur_im = imread(cur_fname);
             
        cur_im_patch = reshape(double(cur_im(new_crop)),length(ypatch_inds_cropped),length(xpatch_inds_cropped));
        cent_gauss_mask = get_pgauss_mask_v2(XXc,YYc,cent_gauss_params);
        surr_gauss_mask = get_pgauss_mask_v2(XXc,YYc,surr_gauss_params);
        mask_diff = surr_gauss_mask - cent_gauss_mask;
        mask_diff(mask_diff < 0) = 0;
        
        cent_im = cur_im_patch.*cent_gauss_mask;
        surr_im = cur_im_patch.*surr_gauss_mask;
%         cent_im = cent_im - mean(cent_im(:));
%         surr_im = surr_im - mean(surr_im(:));
        
        fft_cent = abs(fftshift(fft2(cent_im)));
        fft_surr = abs(fftshift(fft2(surr_im)));
        fft_surr = smoothn(fft_surr,2);
        
        ff = linspace(-Fs/2,Fs/2,length(xpatch_inds_cropped));
        
        subplot(2,2,1)
        imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),flipud(cent_im));colormap(gray);set(gca,'ydir','normal');
        subplot(2,2,2)
        imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),flipud(surr_im));colormap(gray);set(gca,'ydir','normal');
        hold on
         subplot(2,2,3)
       imagesc(ff,ff,(fft_cent)); %caxis([-0.5 1]); %colorbar
       xlim([-10 10]);ylim([-10 10]);caxis([0 15])
        subplot(2,2,4)
        imagesc(ff,ff,(fft_surr)); %caxis([-0.75 0.25]); %colorbar
       xlim([-10 10]);ylim([-10 10]);caxis([0 2.5])
        
        pause
    end   
end

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
xw = 1;yw=1;

cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
%fit smoothed retinotopic surface
orientations = linspace(0,pi-pi/12,12);
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

cd ~/Data/bruce/Expt_1_8_13_imfolder


inc = 1:1:stim_dur;

id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));

%%
close all
ss = 3;
cur_rep = repeat_inds(rep_set(ss));
i = 5;
stim_cnt = 7*(ss-1)+i;
Vq = nan(length(inc),size(Xi,1),size(Xi,2));
for j = 1:length(inc)
    cur_set = squeeze(stim_avg_out(inc(j),stim_cnt,:));
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
    Vq(j,:,:) = F(Xi,Yi);
end
% Vq = abs(Vq);

figure(3)
if cur_rep < 1e4
    cur_fname = sprintf('IM1%.6d.png',cur_rep+i);
elseif cur_rep < 1e5
    cur_fname = sprintf('IM1%.5d.png',cur_rep+i);
else cur_rep < 1e5
    cur_fname = sprintf('IM1%.4d.png',cur_rep+i);
end
cur_fname
cur_im = imread(cur_fname);
imagesc(xax,yax,flipud(cur_im));colormap(gray);set(gca,'ydir','normal');
hold on
rectangle('Position',[x0_avg-xw/2 y0_avg-yw/2 xw,yw],'edgecolor','r');
    plot(interp_x,interp_y,'ro')
xlim([-2 2]); ylim([-2 2])

figure(1)
% pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(stim_avg_out(:,stim_cnt,use_ordx))');shading flat
cur_avg = squeeze(mean(stim_avg_out(:,stim_cnt,:),3));
plot(lags(1:stim_dur),cur_avg,'linewidth',2);shading flat
hold on
% caxis([-7 7]*1e-4)
xlim([0 0.47])

figure(2);
pause;
for j = 1:length(inc)
    j
    figure(2)
    pcolor(xi,yi,squeeze(Vq(j,:,:))); shading flat; colorbar;
    caxis([-5 5]*1e-4);
    figure(1);
    plot(lags(inc(j)),cur_avg(inc(j)),'ro','linewidth',2)
    pause(0.02);
    drawnow();
end

%%
    figure(1)
    subplot(2,1,1)
    pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(stim_avg_out(:,stim_cnt,use_ordx))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(astim_avg_out(:,stim_cnt,use_ordx))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(all_astim_cond_mean(:,stim_cnt,use_ordx))');shading flat
    title('Avg LFPs x-sort')
    yl = ylim();
    line(lags([inc(j) inc(j)]),yl)
    subplot(2,1,2)
    pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(astim_avg_out(:,stim_cnt,use_ordy))');shading flat
%     pcolor(lags(1:stim_dur),1:length(use_lfps),squeeze(all_astim_cond_mean(:,stim_cnt,use_ordy))');shading flat
    title('Avg LFPs y-sort')
    line(lags([inc(j) inc(j)]),yl)
    figure(2)
    if cur_rep < 1e4
        cur_fname = sprintf('IM1%.6d.png',cur_rep+i);
    elseif cur_rep < 1e5
        cur_fname = sprintf('IM1%.5d.png',cur_rep+i);
    else cur_rep < 1e5
        cur_fname = sprintf('IM1%.4d.png',cur_rep+i);
    end
    cur_fname
    cur_im = imread(cur_fname);
    imagesc(xax,yax,flipud(cur_im));colormap(gray);set(gca,'ydir','normal');
    hold on
    rectangle('Position',[x0_avg-xw/2 y0_avg-yw/2 xw,yw],'edgecolor','r');
    xlim([-2 2]); ylim([-2 2])
    plot(interp_x,interp_y,'ro')
    
    figure(3)
    subplot(2,1,1)
    cur_set = squeeze(astim_avg_out(inc(j),stim_cnt,:));
%     cur_set = squeeze(all_astim_cond_mean(inc(j),stim_cnt,:));
    temp = nan(10,10);
    temp(use_ids) = cur_set(id_mat(use_ids));
    imagesc(temp);colorbar
    set(gca,'ydir','normal')
    subplot(2,1,2)
%     Vq = griddata(interp_x,interp_y,cur_set,Xi,Yi);
    pcolor(xi,yi,Vq); shading flat; colorbar;


