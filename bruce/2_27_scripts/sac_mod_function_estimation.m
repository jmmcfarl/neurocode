clear all
close all

cd ~/Data/bruce/2_27_12/M232
% load ./random_bar_eyedata.mat
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians
%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 20;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_trialvec = [];
full_taxis = [];
full_bar_pos = [];
full_bar_Xmat = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_e0 = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
%         cur_bar_e0 = [Expts{bar_expts(ee)}.Trials(tt).e0];
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_binned_spks = nan(length(cur_t_axis),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
%         cur_interp_eye_perp = interp1(all_eye_ts{ee},fix_perp{ee},cur_t_axis);
%         cur_bar_Op_cor = cur_bar_Op - cur_interp_eye_perp;
        
        cur_bar_mat = zeros(length(cur_t_axis),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
%             cur_bar_mat(cur_set,bb) = cur_bar_e0(cur_set);
        end
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
%         all_bar_e0 = [all_bar_e0; cur_bar_e0(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
    cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    
    all_blink_inds = zeros(size(all_t_axis));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
    all_used_inds(all_blink_inds == 1) = 0;
    all_used_inds = logical(all_used_inds);    
    
    full_Xmat = [full_Xmat; all_bar_Xmat(all_used_inds,:)];
    full_spkbinned = [full_spkbinned; all_binned_spikes(all_used_inds,:)];
    full_exptvec = [full_exptvec; ones(sum(all_used_inds),1)*ee];
    full_taxis = [full_taxis; all_t_axis(all_used_inds)];
    full_trialvec = [full_trialvec; all_trial_vec(all_used_inds)];
    full_bar_pos = [full_bar_pos; all_bar_Op(all_used_inds)];
end

%%
lin_X = zeros(length(full_taxis),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_taxis==i);
    lin_X(cur_set,i) = 1;
end

[c,ia,ic] = unique([full_exptvec full_trialvec],'rows');
n_trials = length(ia);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

%%
cd /home/james/Data/bruce/2_27_12/M232
% save random_bar_phasemods_stimmods glm_fit
load ./random_bar_unit_models_ftime glm_fit

%%
all_sac_inds = [];
all_big_inds = [];
all_microsac_inds = [];
all_firstsac_inds = [];
all_secondsac_inds = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    sac_amps = all_sac_peakamps{ee};
    
    cur_sac_set = find(all_expt_num==ee);
    if length(cur_sac_set) ~= length(sac_start_times)
        error('saccade mismatch')
    end
    intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));

    infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    big_sac_set = [infirst_sac_set; insecond_sac_set];
    micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
    micro_sac_set(ismember(micro_sac_set,big_sac_set)) = [];

    sac_times = sac_start_times(intrial_sac_set);
    first_sac_times = sac_start_times(infirst_sac_set);
    second_sac_times = sac_start_times(insecond_sac_set);
    micro_sac_times = sac_start_times(micro_sac_set);
    big_sac_times = sac_start_times(big_sac_set);
    
    cur_e_set = find(full_exptvec == ee);
    
    sac_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),sac_times));
    sac_inds(isnan(sac_inds)) = [];
    micro_sac_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),micro_sac_times));
    micro_sac_inds(isnan(micro_sac_inds)) = [];
    first_sac_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),first_sac_times));
    first_sac_inds(isnan(first_sac_inds)) = [];
    second_sac_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),second_sac_times));
    second_sac_inds(isnan(second_sac_inds)) = [];
    big_sac_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),big_sac_times));
    big_sac_inds(isnan(big_sac_inds)) = [];
    
    all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
    all_big_inds = [all_big_inds; cur_e_set(big_sac_inds(:))];
    all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
    all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
    all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
end

%%
dt = 0.01;
cur_dt = 0.02;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.7];
tent_centers = round(tent_centers/dt);
ntents = length(tent_centers);
% tbmat = construct_tent_bases(tent_centers,1);
tbmat = zeros(ntents,max(tent_centers)+1);
for tt = 1:length(tent_centers)
    cur_beg = (tt-1)*(cur_dt/dt)+1;
    tbmat(tt,cur_beg:(cur_beg+1)) = 1;
end


[ntents,tblen] = size(tbmat);

shift = round(0.2/cur_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;


trial_sac_inds = zeros(size(full_taxis));
trial_sac_inds(all_big_inds) = 1;
trial_sac_mat = zeros(length(full_taxis),ntents);
for i = 1:ntents
    trial_sac_mat(:,i) = conv(trial_sac_inds,tbmat(i,:),'same');
end

trial_msac_inds = zeros(size(full_taxis));
trial_msac_inds(all_microsac_inds) = 1;
trial_msac_mat = zeros(length(full_taxis),ntents);
for i = 1:ntents
    trial_msac_mat(:,i) = conv(trial_msac_inds,tbmat(i,:),'same');
end

trial_fsac_inds = zeros(size(full_taxis));
trial_fsac_inds(all_firstsac_inds) = 1;
trial_fsac_mat = zeros(length(full_taxis),ntents);
for i = 1:ntents
    trial_fsac_mat(:,i) = conv(trial_fsac_inds,tbmat(i,:),'same');
end

trial_ssac_inds = zeros(size(full_taxis));
trial_ssac_inds(all_secondsac_inds) = 1;
trial_ssac_mat = zeros(length(full_taxis),ntents);
for i = 1:ntents
    trial_ssac_mat(:,i) = conv(trial_ssac_inds,tbmat(i,:),'same');
end

%%
n_bar_pos = length(un_bar_pos);
clear temp_kern
for cc = 1:24
    k_mat = reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos);
    temp_kern(cc,:) = sqrt(mean(k_mat.^2,2));
    [~,best_lag(cc)] = max(temp_kern(cc,:));
end
pref_lag = flen - best_lag - 1;
[bps,fls] = meshgrid(1:n_bar_pos,1:flen);
%%
nbins = 15;
bin_edges = linspace(0,100,nbins+1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);

% cur_sac_mat = trial_msac_mat;
cur_sac_mat = trial_sac_mat;
% cur_sac_mat = trial_fsac_mat;

clear cur_stas cur_stas_g g_dist g_x g_bin_cents
for cc = 1:24
    cc
    cur_set = find(fls == best_lag(cc));
    cur_X = full_Xmat(:,cur_set);
    
    Robs = full_spkbinned(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(full_spkbinned(tr_inds,cc));

    ov_sta(cc,:) = Robs'*cur_X./sum(cur_X);
    
    cur_stas(cc,:,:) = nan(length(tent_centers),n_bar_pos);
    for tt = 1:length(tent_centers)
       cur_inds = find(cur_sac_mat(:,tt) > 0);
       cur_stas(cc,tt,:) = (Robs(cur_inds)'*cur_X(cur_inds,:))./sum(cur_X(cur_inds,:));
    end
    
    [~, ~, ~, prate, g] = getLL_GNM(glm_fit(cc),full_Xmat,tr_spkbns,'none');
    pp = prctile(g,bin_edges);
    g_bin_cents(cc,:) = 0.5*pp(1:end-1) + 0.5*pp(2:end);
    cur_stas_g(cc,:,:) = zeros(length(tent_centers),nbins);
    for tt = 1:length(tent_centers)
       cur_inds = find(cur_sac_mat(:,tt) > 0);
       for ii = 1:nbins
           cur_set = cur_inds(g(cur_inds) >= pp(ii) & g(cur_inds) < pp(ii+1));
        cur_stas_g(cc,tt,ii) = mean(Robs(cur_set));
       end
    end
    
    [g_dist(cc,:),g_x(cc,:)] = ksdensity(g);
    
end

%%
sm_sig = 1;
for cc = 1:24
    for tt = 1:length(tent_centers)
        cur_stas(cc,tt,2:end) = jmm_smooth_1d_cor(squeeze(cur_stas(cc,tt,2:end)),sm_sig);
    end
end
%%
sm_sig = 1;
for cc = 1:24
    for tt = 1:length(tent_centers)
        cur_stas_g(cc,tt,:) = jmm_smooth_1d_cor(squeeze(cur_stas_g(cc,tt,:)),sm_sig);
    end
end

%%
zval = find(tent_centers==0);
eval = find(tent_centers*dt > 0.2,1,'first');
sac_mod = squeeze(mean(cur_stas,3));
stim_mod = squeeze(mean(cur_stas,2));
g_trans = squeeze(mean(cur_stas_g,2));
for cc = 1:24
        [~,temp] = min(sac_mod(cc,zval:eval));
    post_sac_min(cc) = zval - 1 + temp;
    [~,temp] = max(sac_mod(cc,post_sac_min(cc):eval));
    post_sac_max(cc) = post_sac_min(cc) - 1 + temp;
end
%%
save gsac_stim_interaction flen dt glm_fit flen n_bar_pos un_bar_pos post_sac_* sac_mod best_lag tent_centers cur_stas* stim_mod g_*
%%
close all
f_tax = dt*(0:-1:(-flen+1));
for cc = 1:24
    cc
    subplot(3,2,1)
    kmat = reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos);
    imagesc(un_bar_pos,f_tax,flipud(kmat))
    xl = xlim();
    line(xl,f_tax(flen - best_lag([cc cc])+1),'color','w')
    xlabel('Bar position (deg)','fontsize',16)
    ylabel('Time lag (s)','fontsize',16)
    title('STRF')
    
    subplot(3,2,2)
    plot(tent_centers*dt,sac_mod(cc,:)/dt)
    hold on
    plot(tent_centers(post_sac_min(cc))*dt,sac_mod(cc,post_sac_min(cc))/dt,'ro','markersize',14,'linewidth',2)
    plot(tent_centers(post_sac_max(cc))*dt,sac_mod(cc,post_sac_max(cc))/dt,'o','color',[0.2 0.8 0.2],'markersize',14,'linewidth',2)
    ylabel('Firing rate (Hz)','fontsize',16)
    xlabel('Time since sac onset (s)','fontsize',16)
    title('Saccade-triggered average rate')
    xlim([-0.1 0.5])
    
    subplot(3,2,3)
    imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(cur_stas(cc,:,2:end))')
    xlabel('Time since sac onset (s)','fontsize',16)
    ylabel('Bar position (deg)','fontsize',16)
    title('Saccade-triggered stimulus tuning (best latency)')
    
    subplot(3,2,4)
    plot(un_bar_pos(2:end),stim_mod(cc,2:end)/dt);
    hold on
    plot(un_bar_pos(2:end),squeeze(cur_stas(cc,post_sac_min(cc),2:end)/dt),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas(cc,post_sac_max(cc),2:end)/dt),'color',[0.2 0.8 0.2])
    xlabel('Bar position (deg)','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)
    title('Stimulus tuning time slices')
      xlim(un_bar_pos([2 end]))
 
    subplot(3,2,5)
    imagesc(tent_centers*dt,g_bin_cents(cc,:),squeeze(cur_stas_g(cc,:,:))')
    set(gca,'ydir','normal')
    xlabel('Time since sac onset (s)','fontsize',16)
    ylabel('Generating signal','fontsize',16)
    title('Saccade modulation of spiking NL')
    
    subplot(3,2,6)
    plot(g_bin_cents(cc,:),g_trans(cc,:)/dt);
    hold on
    plot(g_bin_cents(cc,:),squeeze(cur_stas_g(cc,post_sac_min(cc),:))/dt,'r')
    plot(g_bin_cents(cc,:),squeeze(cur_stas_g(cc,post_sac_max(cc),:))/dt,'color',[0.2 0.8 0.2])
    yl = ylim();
    off = yl(1);
    sc_fac = 0.75*(yl(2) - yl(1));
    xlabel('Generating signal','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)
    xlim([g_bin_cents(cc,1) - 0.05 g_bin_cents(cc,end) + 0.05])
    title('Spiking NL time slices')
    
    plot(g_x(cc,:),g_dist(cc,:)/max(g_dist(cc,:))*sc_fac + off,'k')
    pause
    clf
end






%%
% clear cur_stas
for cc = 1:24
    
    cur_set = find(fls == best_lag(cc));
    cur_X = full_Xmat(:,cur_set);
    
    Robs = full_spkbinned(tr_inds,cc);

    ov_sta(cc,:) = Robs'*cur_X./sum(cur_X);
    
    cur_stas_big(cc,:,:) = nan(length(tent_centers),n_bar_pos);
    for tt = 1:length(tent_centers)
       cur_inds = find(trial_sac_mat(:,tt) > 0);
       cur_stas_big(cc,tt,:) = (Robs(cur_inds)'*cur_X(cur_inds,:))./sum(cur_X(cur_inds,:));
    end
    
end

%%
sm_sig = 1;
for cc = 1:24
    for tt = 1:length(tent_centers)
        cur_stas_big(cc,tt,2:end) = jmm_smooth_1d_cor(squeeze(cur_stas_big(cc,tt,2:end)),sm_sig);
    end
end
norm_stas = bsxfun(@rdivide,cur_stas_big,mean(cur_stas_big,3));
%%
zval = find(tent_centers==0);
eval = find(tent_centers*dt > 0.2,1,'first');
big_sac_mod = squeeze(mean(cur_stas_big,3));
big_stim_mod = squeeze(mean(cur_stas_big,2));
for cc = 1:24
        [~,temp] = min(big_sac_mod(cc,zval:eval));
    post_big_sac_min(cc) = zval - 1 + temp;
    [~,temp] = max(big_sac_mod(cc,post_big_sac_min(cc):eval));
    post_big_sac_max(cc) = post_big_sac_min(cc) - 1 + temp;
end

%%
close all
for cc = 1:24
    cc
    subplot(2,2,1)
    kmat = reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos);
    imagesc(un_bar_pos,1:flen,kmat)
    xl = xlim();
    line(xl,best_lag([cc cc]),'color','w')
    
    subplot(2,2,2)
    plot(tent_centers*dt,big_sac_mod(cc,:))
    hold on
    plot(tent_centers(post_big_sac_min(cc))*dt,big_sac_mod(cc,post_big_sac_min(cc)),'ro')
    plot(tent_centers(post_big_sac_max(cc))*dt,big_sac_mod(cc,post_big_sac_max(cc)),'ko')
    subplot(2,2,3)
    imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(cur_stas_big(cc,:,2:end))')
    subplot(2,2,4)
    plot(un_bar_pos(2:end),big_stim_mod(cc,2:end));
    hold on
    plot(un_bar_pos(2:end),squeeze(cur_stas_big(cc,post_big_sac_min(cc),2:end)),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas_big(cc,post_big_sac_max(cc),2:end)),'k')
%     subplot(3,2,5)
%     imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(norm_stas(cc,:,2:end))')
    pause
    clf
end

%%
for cc = 1:24
   msac_sup_tune(cc,:) = squeeze(cur_stas(cc,post_sac_min(cc),2:end));
   msac_exc_tune(cc,:) = squeeze(cur_stas(cc,post_sac_max(cc),2:end));
   sac_sup_tune(cc,:) = squeeze(cur_stas_big(cc,post_big_sac_min(cc),2:end));
   sac_exc_tune(cc,:) = squeeze(cur_stas_big(cc,post_big_sac_max(cc),2:end));
       
end
%%
close all
for cc = 1:24
    cc
    subplot(2,4,1)
    kmat = reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos);
    imagesc(un_bar_pos,1:flen,kmat)
    xl = xlim();
    line(xl,best_lag([cc cc]),'color','w')
    
    subplot(2,4,2)
    plot(tent_centers*dt,sac_mod(cc,:))
    hold on
    plot(tent_centers(post_sac_min(cc))*dt,sac_mod(cc,post_sac_min(cc)),'ro')
    plot(tent_centers(post_sac_max(cc))*dt,sac_mod(cc,post_sac_max(cc)),'ko')
    subplot(2,4,5)
    imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(cur_stas(cc,:,2:end))')
    subplot(2,4,6)
    plot(un_bar_pos(2:end),stim_mod(cc,2:end));
    hold on
    plot(un_bar_pos(2:end),squeeze(cur_stas(cc,post_sac_min(cc),2:end)),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas(cc,post_sac_max(cc),2:end)),'k')
%     subplot(3,2,5)
%     imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(norm_stas(cc,:,2:end))')

subplot(2,4,3)
    kmat = reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos);
    imagesc(un_bar_pos,1:flen,kmat)
    xl = xlim();
    line(xl,best_lag([cc cc]),'color','w')
    
    subplot(2,4,4)
    plot(tent_centers*dt,big_sac_mod(cc,:))
    hold on
    plot(tent_centers(post_big_sac_min(cc))*dt,big_sac_mod(cc,post_big_sac_min(cc)),'ro')
    plot(tent_centers(post_big_sac_max(cc))*dt,big_sac_mod(cc,post_big_sac_max(cc)),'ko')
    subplot(2,4,7)
    imagesc(tent_centers*dt,un_bar_pos(2:end),squeeze(cur_stas_big(cc,:,2:end))')
    subplot(2,4,8)
    plot(un_bar_pos(2:end),big_stim_mod(cc,2:end));
    hold on
    plot(un_bar_pos(2:end),squeeze(cur_stas_big(cc,post_big_sac_min(cc),2:end)),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas_big(cc,post_big_sac_max(cc),2:end)),'k')
   pause
    clf
end

%%
close all
for i = 1:24
    subplot(2,1,1)
plot(msac_sup_tune(i,:),msac_exc_tune(i,:),'o')
hold on
pp = polyfit(msac_sup_tune(i,:),msac_exc_tune(i,:),1);
xx = linspace(0,max(msac_exc_tune(i,:)),100);
plot(xx,polyval(pp,xx),'b')
line([0 max(msac_exc_tune(i,:))],[0 max(msac_exc_tune(i,:))],'color','k')
% plot(sac_sup_tune(i,:),sac_exc_tune(i,:),'ro')
xl = xlim();
yl = ylim();
xlim([0 max([xl(2) yl(2)])])
ylim([0 max([xl(2) yl(2)])])
    subplot(2,1,2);hold on
   plot(un_bar_pos(2:end),squeeze(cur_stas(i,post_sac_min(i),2:end)),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas(i,post_sac_max(i),2:end)),'k')
pause
clf
end
%%
close all
for i = 1:24
    subplot(2,1,1)
plot(sac_sup_tune(i,:),sac_exc_tune(i,:),'o')
hold on
pp = polyfit(sac_sup_tune(i,:),sac_exc_tune(i,:),1);
xx = linspace(0,max(sac_exc_tune(i,:)),100);
plot(xx,polyval(pp,xx),'b')
line([0 max(sac_exc_tune(i,:))],[0 max(sac_exc_tune(i,:))],'color','k')
xl = xlim();
yl = ylim();
xlim([0 max([xl(2) yl(2)])])
ylim([0 max([xl(2) yl(2)])])
    subplot(2,1,2);hold on
   plot(un_bar_pos(2:end),squeeze(cur_stas_big(i,post_big_sac_min(i),2:end)),'r')
    plot(un_bar_pos(2:end),squeeze(cur_stas_big(i,post_big_sac_max(i),2:end)),'k')
pause
clf
end
