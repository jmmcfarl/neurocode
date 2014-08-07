fname = 'sacStimProc_sta';
anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
cd(anal_dir)
load(fname);

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';

%%
cc = 101;
% fprintf('Starting model fits for unit %d\n',cc);
% loo_cc = find(loo_set == cc); %index within the LOOXV set
% cc_uinds = full_inds(~isnan(Robs_mat(full_inds,cc))); %set of used indices where this unit was isolated
% 
% cur_Robs = Robs_mat(cc_uinds,cc);
% 
% cur_GQM = ModData(cc).rectGQM;
% sacStimProc(cc).ModData = ModData(cc);
% sacStimProc(cc).used = true;
% 
% fprintf('Reconstructing retinal stim for unit %d\n',cc);
% if ismember(cc,loo_set) %if unit is member of LOOXV set, use its unique EP sequence
%     cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
%     cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
%     cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
%     cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
%     [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
%         cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
%     
%     fin_shift_cor = round(fin_tot_corr);
%     
%     %RECOMPUTE XMAT
%     all_shift_stimmat_up = all_stimmat_up;
%     if ~fit_unCor
%         for i=1:NT
%             all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
%         end
%     end
%     all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%     all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
%     
% else %otherwise use overall EP sequence
%     all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
%     all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
% end
% 
% %% FOR GSACS
% cur_Xsac = Xsac(cc_uinds,:); %saccade indicator Xmat
% 
% if is_TBT_expt
%     gs_trials = find(all_trial_Ff > 0);
%     gs_inds = find(ismember(all_trialvec(used_inds),gs_trials));
% else
%     gs_blocks = [imback_gs_expts; grayback_gs_expts];
%     gs_inds = find(ismember(all_blockvec(used_inds),gs_blocks));
% end
% 
% %only use indices during guided saccade expts here
% any_sac_inds = find(ismember(cc_uinds,gs_inds));
% 
% %% Fit spk NL params and refit scale of each filter using target data (within trange of sacs)
% cur_GQM = NMMfit_logexp_spkNL(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
% cur_GQM = NMMfit_scale(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
% 
% stim_mod_signs = [cur_GQM.mods(:).sign];
% [~,~,~,~,filt_outs,fgint] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% fgint = bsxfun(@times,fgint,stim_mod_signs);
% stimG = sum(fgint,2);

%%
cur_gdist = sacStimProc(cc).gsac_TB_gdist;
cur_post_mod = sacStimProc(cc).gsac_post_singmod;
sac_dep_theta = cur_post_mod.spk_NL_params(1) + cur_post_mod.mods(2).filtK;
sac_dep_beta = cur_post_mod.mods(3).filtK + cur_post_mod.mods(1).filtK;

Gtick = sacStimProc(cc).gsac_TB_gX;
Xtick = sacStimProc(cc).gsac_TB_lagX;
sac_dep_gen = bsxfun(@plus,bsxfun(@times,sac_dep_beta,Gtick),sac_dep_theta);
sac_dep_rate = cur_post_mod.spk_NL_params(3)*log(1+exp(cur_post_mod.spk_NL_params(2)*sac_dep_gen));

base_out = cur_GQM.spk_NL_params(3)*log(1+exp(cur_GQM.spk_NL_params(2)*(Gtick + cur_GQM.spk_NL_params(1))));


[mmm,mmmloc] = minmax(sacStimProc(cc).gsac_avg_rate);
mmmloc(2) = mmmloc(2) - 1; %manually set this one bin earlier

h1=figure;
subplot(2,1,1)
plot(slags*dt,sacStimProc(cc).gsac_avg_rate/dt);
hold on
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xl = xlim();
line(xl,[0 0]+sacStimProc(cc).ModData.unit_data.avg_rate,'color','k');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
subplot(2,1,2);
plot(slags*dt,sacStimProc(cc).gsac_post_singmod.mods(3).filtK);
hold on
plot(slags*dt,sacStimProc(cc).gsac_post_singmod.mods(2).filtK,'r');
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xlabel('Time (s)');
ylabel('Filter');

h2=figure;
subplot(2,1,1);
imagesc(slags*dt,Gtick,sac_dep_rate'); set(gca,'ydir','normal'); caxis([0 max(sac_dep_rate(:))]*0.9);
% imagesc(slags*dt,Gtick,log10(sac_dep_rate')); set(gca,'ydir','normal'); caxis(log10([0.01 max(sac_dep_rate(:))]));
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','w');
line(slags(mmmloc([2 2]))*dt,yl,'color','w');
xlabel('Time (s)');
ylabel('Generating signal');
subplot(2,1,2)
imagesc(slags*dt,Gtick,sacStimProc(cc).gsac_TB_rate); set(gca,'ydir','normal'); caxis([0 max(sac_dep_rate(:))]*0.9);
% imagesc(slags*dt,Gtick,log10(sacStimProc(cc).gsac_TB_rate)); set(gca,'ydir','normal');  caxis(log10([0.01 max(sac_dep_rate(:))]));
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','w');
line(slags(mmmloc([2 2]))*dt,yl,'color','w');
xlabel('Time (s)');
ylabel('Generating signal');


h3=figure;
subplot(2,1,1);
plot(Gtick,sac_dep_rate(mmmloc(1),:)/dt); hold on
plot(Gtick,sacStimProc(cc).gsac_TB_rate(:,mmmloc(1))/dt,'r');
plot(Gtick,base_out/dt,'k','linewidth',2);
xlim(Gtick([1 end]));
yl = ylim();
plot(Gtick,cur_gdist/max(cur_gdist)*yl(2)*0.8,'k--');
xlabel('Generating signal');
ylabel('Firing rate (Hz)');
subplot(2,1,2);
plot(Gtick,sac_dep_rate(mmmloc(2),:)/dt); hold on
plot(Gtick,sacStimProc(cc).gsac_TB_rate(:,mmmloc(2))/dt,'r');
plot(Gtick,base_out/dt,'k','linewidth',2);
xlim(Gtick([1 end]));
yl = ylim();
plot(Gtick,cur_gdist/max(cur_gdist)*yl(2)*0.8,'k--');
xlabel('Generating signal');
ylabel('Firing rate (Hz)');

h4=figure;
subplot(2,1,1)
plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo);
hold on
plot(Xtick*dt,sacStimProc(cc).gsac_TB_info,'r');
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xlim(slags([1 end])*dt);
xlabel('Time (s)');
ylabel('Single-spike info (bits/spike)');
subplot(2,1,2)
plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo'.*sacStimProc(cc).gsac_avg_rate);
hold on
plot(Xtick*dt,sacStimProc(cc).gsac_TB_info.*sacStimProc(cc).gsac_TB_avg_rate,'r');
yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xlim(slags([1 end])*dt);
xlabel('Time (s)');
ylabel('Info rate (bits/sec)');


cid = sprintf('E%d_C%d_',Expt_num,cc);

% fig_width = 5; rel_height = 1.6;
% figufy(h1);
% fname = [fig_dir cid 'Filtkerns.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% figufy(h2);
% fname = [fig_dir cid '2drates.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);
% 
% figufy(h3);
% fname = [fig_dir cid 'spkNLexamps.pdf'];
% exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h3);
% 
% figufy(h4);
% fname = [fig_dir cid 'Infokerns.pdf'];
% exportfig(h4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h4);

%%
t_ax = (0:(flen-1))*dt + dt/2;
p_ax = sp_dx*((1:use_nPix_us)-use_nPix_us/2);

yy = [-0.35 0.35];

ov_sta = sacStimProc(cc).ov_phaseDep_sta;
ov_staA = sacStimProc(cc).ov_phaseInd_sta;
% load sacStimProc_sta
% ov_sta = sum(bsxfun(@times,all_Xmat_shift,cur_Robs))/sum(cur_Robs);
% base_a = mean(all_Xmat_shift);
% ov_sta = ov_sta - base_a;
% 
% ov_staA = sum(bsxfun(@times,abs(all_Xmat_shift),cur_Robs))/sum(cur_Robs);
% base_a = mean(abs(all_Xmat_shift));
% ov_staA = ov_staA - base_a;

temp = sacStimProc(cc).gsac_phaseDep_sta;
tempa = sacStimProc(cc).gsac_phaseInd_sta;

h1=figure;
subplot(3,1,1)
imagesc(p_ax,t_ax,reshape(ov_sta,flen,use_nPix_us));
ylim([0 0.11]);
xlim(yy);
xlabel('Rel position (deg)');
ylabel('Lag (s)');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
subplot(3,1,2);
imagesc(p_ax,t_ax,squeeze(temp(mmmloc(1),:,:))); caxis([-cam cam]);
xlim(yy);
ylim([0 0.11]);
xlabel('Rel position (deg)');
ylabel('Lag (s)');
subplot(3,1,3);
imagesc(p_ax,t_ax,squeeze(temp(mmmloc(2),:,:))); caxis([-cam cam]);
ylim([0 0.11]);
xlim(yy);
xlabel('Rel position (deg)');
ylabel('Lag (s)');

% yy = [5 25];
h2=figure;
for jj = 1:4
    subplot(4,1,jj);
    imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(temp(:,3+jj,:))');
    caxis([-cam cam]);ylim(yy); yl = ylim();
    line(slags(mmmloc([1 1]))*dt,yl,'color','k');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k');
    xlabel('Time (s)');
    ylabel('Rel Position (deg)');
end


cid = sprintf('E%d_C%d_',Expt_num,cc);

% fig_width = 5; rel_height = 3;
% figufy(h1);
% fname = [fig_dir cid 'STAexamps.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% fig_width = 5; rel_height = 4;
% figufy(h2);
% fname = [fig_dir cid 'STAslices.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);

%%
cur_Filts = reshape([cur_GQM.mods(:).filtK],[flen use_nPix_us length(cur_GQM.mods)]);
filt_tempkerns = squeeze(std(cur_Filts,[],2));
filt_Signs = [cur_GQM.mods(:).sign];

eKern = mean(filt_tempkerns(:,filt_Signs==1),2);
Ikern = mean(filt_tempkerns(:,filt_Signs==-1),2);

t_ax = (0:(flen-1))*dt + dt/2;
h1=figure;
subplot(3,1,3);
plot(t_ax,eKern,t_ax,Ikern,'r');
xlabel('Lag (s)');
ylabel('Filter env');
legend('Efilts','Ifilts');
xlim([0 0.15]);

subplot(3,1,2);
plot(slags*dt,1+sacStimProc(cc).gsac_post_Egains); hold on
plot(slags*dt,1+sacStimProc(cc).gsac_post_Igains,'r');
xlabel('Time (s)');
ylabel('Gain');
legend('Efilts','Ifilts');
xlim([0 0.15]);

subplot(3,1,1);
plot(slags*dt,1+sacStimProc(cc).gsac_post_Egains); hold on
plot(slags*dt,1+sacStimProc(cc).gsac_post_Igains,'r');
xlabel('Time (s)');
ylabel('Gain');
legend('Efilts','Ifilts');

% fig_width = 5; rel_height = 2.5;
% figufy(h1);
% fname = [fig_dir cid 'EIkerns.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%%

h1=figure;
plot(slags*dt,sacStimProc(cc).gsacGainMod.gain_kernel);
hold on
plot(slags*dt,sacStimProc(cc).gsacGainMod.stim_kernel,'k');
plot(slags*dt,sacStimProc(cc).gsac_post_singmod.mods(3).filtK,'r');
xlabel('Time (s)');
ylabel('Gain kernels');
legend('Post (both)','Pre','Post only');

% fig_width = 5; rel_height = 0.8;
% figufy(h1);
% fname = [fig_dir cid 'Prepost.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%%
mod_ekern = nan(length(slags),flen);
mod_ikern = nan(length(slags),flen);
stim_kern = sacStimProc(cc).gsacGainMod.stim_kernel;
gain_kern = 1+sacStimProc(cc).gsacGainMod.gain_kernel;
for ii = 1:length(slags)
    khist = (ii-flen+1):ii;
    uset = find(khist > 0);
    tkern = ones(flen,1);
    tkern(uset) = tkern(uset) + stim_kern(khist(uset));
    tkern = flipud(tkern);
    
    mod_ekern(ii,:) = (eKern.*tkern)';
    mod_ikern(ii,:) = (Ikern.*tkern)';
end

% mod_gkern = mod_tkern;
% mod_gkern = bsxfun(@times,mod_gkern,gain_kern);

h1=figure;
subplot(2,1,1)
imagesc(slags*dt,t_ax,mod_ekern'); set(gca,'ydir','normal');
yl = ylim();
line([0 0],yl,'color','w');
xlabel('Time (s)');
ylabel('Lag (s)');
line(slags([17 17])*dt,yl,'color','g');
line(slags([19 19])*dt,yl,'color','r');

subplot(2,1,2)
imagesc(slags*dt,t_ax,mod_ikern'); set(gca,'ydir','normal');
line([0 0],yl,'color','w');
xlabel('Time (s)');
ylabel('Lag (s)');

h2=figure;hold on
plot(t_ax,eKern,'k');
plot(t_ax,mod_ekern(17,:),'g');
plot(t_ax,mod_ekern(19,:),'r');
xlabel('Lag (s)');
ylabel('Temp kernel');

fig_width = 5; rel_height = 1.6;
figufy(h1);
fname = [fig_dir cid 'tempkern_mod.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

fig_width = 5; rel_height = 0.8;
figufy(h2);
fname = [fig_dir cid 'tempkern_slices.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

%%
h1 = figure;
subplot(2,1,1);
plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo); hold on
plot(slags*dt,sacStimProc(cc).gsac_modinfo,'r')
plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo,'k')
plot(Xtick*dt,sacStimProc(cc).gsac_TB_info,'g');
xlabel('Time (s)');
ylabel('SS info (bits/spk)');
legend('Post-mod','Pre-post mod','Subspace mod','TB mod');
xlim(slags([1 end])*dt)
subplot(2,1,2);
plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo.*sacStimProc(cc).gsac_avg_rate/dt); hold on
plot(slags*dt,sacStimProc(cc).gsac_modinfo.*sacStimProc(cc).gsac_avg_rate/dt,'r')
plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo.*sacStimProc(cc).gsac_avg_rate/dt,'k')
plot(Xtick*dt,sacStimProc(cc).gsac_TB_info.*sacStimProc(cc).gsac_TB_avg_rate/dt,'g');
xlabel('Time (s)');
ylabel('SS info rate (bits/sec)');
legend('Post-mod','Pre-post mod','Subspace mod','TB mod');
xlim(slags([1 end])*dt)

% fig_width = 5; rel_height = 1.6;
% figufy(h1);
% fname = [fig_dir cid 'allinfo.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%%
tempsub = reshape(sacStimProc(cc).gsac_phaseDep_subfilt,[length(slags) flen use_nPix_us]);
tempasub = reshape(sacStimProc(cc).gsac_phaseInd_subfilt,[length(slags) flen use_nPix_us]);
cam = max(abs(tempsub(:)));
cam2 = max(abs(tempasub(:)));

h1=figure;
subplot(3,2,1)
imagesc(p_ax,t_ax,squeeze(mean(tempsub(1:5,:,:)))); caxis([-cam cam]);
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');
subplot(3,2,3);
imagesc(p_ax,t_ax,squeeze(tempsub(mmmloc(1),:,:))); caxis([-cam cam]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');
subplot(3,2,5);
imagesc(p_ax,t_ax,squeeze(tempsub(mmmloc(2),:,:))); caxis([-cam cam]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');

subplot(3,2,2)
imagesc(p_ax,t_ax,squeeze(mean(tempasub(1:5,:,:)))); caxis([-cam cam]);
ca = caxis(); cam2 = max(abs(ca)); caxis([-cam2 cam2]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');
subplot(3,2,4);
imagesc(p_ax,t_ax,squeeze(tempasub(mmmloc(1),:,:)));  caxis([-cam2 cam2]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');
subplot(3,2,6);
imagesc(p_ax,t_ax,squeeze(tempasub(mmmloc(2),:,:)));  caxis([-cam2 cam2]);
xlabel('Rel Pos (deg)');
ylabel('Lag (s)');


h2=figure;
subplot(2,1,1);
imagesc(slags*dt,p_ax,squeeze(tempsub(:,5,:))');
caxis([-cam cam]);ylim(yy); yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xlabel('Time (s)');
ylabel('Rel Pos (deg)');
subplot(2,1,2);
imagesc(slags*dt,p_ax,squeeze(temp(:,5,:))');
caxis([-cam cam]);ylim(yy); yl = ylim();
line(slags(mmmloc([1 1]))*dt,yl,'color','k');
line(slags(mmmloc([2 2]))*dt,yl,'color','k');
xlabel('Time (s)');
ylabel('Rel Pos (deg)');



% fig_width = 8; rel_height = 1.2;
% figufy(h1);
% fname = [fig_dir cid 'Subspace_filts.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
%
% fig_width = 6; rel_height = 1.5;
% figufy(h2);
% fname = [fig_dir cid 'filtslice_compare.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);

%%

sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
g_tot = stimG;
clear tr_stim
tr_stim{1} = [g_tot];
tr_stim{2} = cur_Xsac;
clear sac_stim_params
sac_stim_params(1) = NMMcreate_stim_params(1);
sac_stim_params(2) = NMMcreate_stim_params([size(cur_Xsac,2)]);
mod_signs = [1 1];
Xtargets = [1 2];
NL_types = {'lin','lin','lin'};
post_gsac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
post_gsac_Smod = NMMfit_filters(post_gsac_Smod,cur_Robs,tr_stim,[],[],silent);


%%
is_X_bs = false(size(all_Xmat_shift));
temp = repmat((1:1:(flen))',1,use_nPix_us);
back_mat = create_time_embedding(temp,NMMcreate_stim_params([flen use_nPix_us]));
% uuu = back_mat > 0;
uuu = back_mat == 0;
for ss = 1:length(big_sacs)
    fprintf('%d/%d\n',ss,length(big_sacs));
    cur_stop_ind = saccade_stop_inds(big_sacs(ss));
    %     cur_stop_ind = saccade_start_inds(big_sacs(ss));
    cur_chunk = cur_stop_ind:(cur_stop_ind + (flen-1));
    is_X_bs(cur_chunk,:) = uuu;
end

shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds,use_kInds_up);
shuf_X(~is_X_bs) = all_Xmat_shift(~is_X_bs);

%%
sacGainMod = sacStimProc(cc).gsacGainMod;
[gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
[shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
% null_prate = ones(size(gain_pred_rate))*mean(cur_Robs);
% null_LLseq = cur_Robs.*log2(null_prate)-null_prate;

% str_stim = tr_stim;
% [~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, cur_Robs, shuf_X);
% fgint = bsxfun(@times,fgint,[cur_GQM.mods(:).sign]);
% sg_tot = sum(fgint,2);
% str_stim{1} = sg_tot;

% % [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% % [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
% [base_LL,~,base_predrate] = NMMmodel_eval(post_gsac_Smod,cur_Robs,tr_stim);
% [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(post_gsac_Smod,cur_Robs,str_stim);
% base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
% baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;

[ta_gain,tlags] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,saccade_stop_inds(big_sacs),14,14);
% [tabase_gain,tlags] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,saccade_stop_inds(big_sacs),14,14);

[nspks,tlags] = get_event_trig_avg_v3(cur_Robs,saccade_stop_inds(big_sacs),14,14);
tot_nspks = nspks*length(big_sacs);
ta_gain = ta_gain./tot_nspks;
% tabase_gain = tabase_gain./tot_nspks;


% [sac_fpost_info,sac_spost_info,sac_subpost_info,sac_info,sac_LL,sac_fpost_LL,sac_spost_LL,sac_subpost_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
% for ii = 1:length(slags)
%     temp = find(cur_Xsac(:,ii) == 1);
%     sac_avgrate(ii) = mean(cur_Robs(temp));
%     cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
%     sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
%     sac_Nspks(ii) = sum(cur_Robs(temp));
%
%     [shuf_LL(ii)] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
%     n_frames(ii) = sum(t_since_sac_start == slags(ii));
% end

%%

% [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% cur_Robs = poissrnd(base_predrate);

msac_times = saccade_start_inds(micro_sacs);
msac_isis = [0;diff(msac_times)];
bad = find(msac_isis < 40);
use_micro_sacs = micro_sacs;
use_micro_sacs(bad) = [];

cur_sac_start_inds = saccade_start_inds(big_sacs);
cur_sac_stop_inds = saccade_stop_inds(big_sacs);
% cur_sac_start_inds = saccade_start_inds(use_micro_sacs);
% cur_sac_stop_inds = saccade_stop_inds(use_micro_sacs);

% rand_jit = randi(100,length(big_sacs),1);
% cur_sac_start_inds = cur_sac_start_inds + rand_jit;
% cur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
% % cur_sac_stop_inds = cur_sac_start_inds;
% bad = find(cur_sac_stop_inds > NT);
% cur_sac_start_inds(bad) = [];
% cur_sac_stop_inds(bad) = [];

lagdist = 15;
all_before = [];
all_after = [];
all_during = [];
for ss = 1:length(cur_sac_start_inds)
    cur_before = (cur_sac_start_inds(ss)-lagdist-10):(cur_sac_start_inds(ss)-1);
    cur_after = cur_sac_stop_inds(ss):(cur_sac_stop_inds(ss)+lagdist);
    cur_during = (cur_sac_start_inds(ss)):(cur_sac_stop_inds(ss)-1);
    %     all_before = cat(1,all_before,find(ismember(cc_uinds,cur_before)));
    %     all_after = cat(1,all_after,find(ismember(cc_uinds,cur_after)));
    %     all_during = cat(1,all_during,find(ismember(cc_uinds,cur_during)));
    
    cur_before_tvec = all_trialvec(used_inds(cur_before));
    cur_before_tvec(cur_before_tvec ~= all_trialvec(used_inds(cur_sac_start_inds(ss)))) = [];
    cur_after_tvec = all_trialvec(used_inds(cur_after));
    cur_after_tvec(cur_after_tvec ~= all_trialvec(used_inds(cur_sac_stop_inds(ss)))) = [];
    cur_during_tvec = all_trialvec(used_inds(cur_during));
    cur_during_tvec(cur_during_tvec ~= all_trialvec(used_inds(cur_sac_start_inds(ss)))) = [];
    
    all_before = cat(2,all_before,cur_before);
    all_after = cat(2,all_after,cur_after);
    all_during = cat(2,all_during,cur_during);
end

%% BEFORE
shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds(all_before),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_before),1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);

sacGainMod = sacStimProc(cc).gsacGainMod;
% sacGainMod = sacStimProc(cc).msacGainMod;
[gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);


[shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);

temp_sparams = NMMcreate_stim_params(1);
tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
tempmod = NMMfit_filters(tempmod,cur_Robs,Gshuff);
[~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs,Gshuff);

gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
shuf_LLseq = cur_Robs.*log2(tempprate)-tempprate;

[base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
[base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;

% cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
cur_stop_inds = find(ismember(cc_uinds,cur_sac_stop_inds));
% [before_ta_gain,tlags] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,[],all_trialvec(used_inds(cc_uinds)));

[before_ta_gain,tlags,~,n_events,before_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
before_ta_gain = nansum(before_ta_mat);

[before_base_gain,tlags,~,~,before_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
before_base_gain = nansum(before_base_mat);
%% AFTER
shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds(all_after),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_after),1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);

sacGainMod = sacStimProc(cc).gsacGainMod;
% sacGainMod = sacStimProc(cc).msacGainMod;
[gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
[shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;

[base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
[base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;

% cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
[after_ta_gain,tlags,~,~,after_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
after_ta_gain = nansum(after_ta_mat);
[after_base_gain,tlags,~,~,after_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
after_base_gain = nansum(after_base_mat);
%% DURING
shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds(all_during),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_during),1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);

sacGainMod = sacStimProc(cc).gsacGainMod;
% sacGainMod = sacStimProc(cc).msacGainMod;
[gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
[shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;

[base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
[base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;

% cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
[during_ta_gain,tlags,~,~,during_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
during_ta_gain = nansum(during_ta_mat);
[during_base_gain,tlags,~,~,during_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
during_base_gain = nansum(during_base_mat);
%%
% cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
[nspks,tlags,~,~,nspks_mat] = get_event_trig_avg_v3(cur_Robs,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
tot_nspks = nansum(nspks_mat);

before_ta_Ngain = before_ta_gain./tot_nspks;
before_base_Ngain = before_base_gain./tot_nspks;
after_ta_Ngain = after_ta_gain./tot_nspks;
after_base_Ngain = after_base_gain./tot_nspks;
during_ta_Ngain = during_ta_gain./tot_nspks;
during_base_Ngain = during_base_gain./tot_nspks;
%
% before_ta_Ngainr = before_ta_Ngain.*length(cur_stop_inds);
% before_base_Ngainr = before_base_Ngain.*length(cur_stop_inds);
% after_ta_Ngainr = after_ta_Ngain.*length(cur_stop_inds);
% after_base_Ngainr = after_base_Ngain.*length(cur_stop_inds);
% during_ta_Ngainr = during_ta_Ngain.*length(cur_stop_inds);
% during_base_Ngainr = during_base_Ngain.*length(cur_stop_inds);

%%
figure();
plot(tlags*dt,before_ta_Ngain); hold on
plot(tlags*dt,after_ta_Ngain,'r');
plot(tlags*dt,during_ta_Ngain,'k');


plot(tlags*dt,before_base_Ngain,'--');hold on
plot(tlags*dt,after_base_Ngain,'r--')
plot(tlags*dt,during_base_Ngain,'k--');

xlabel('Time since fixation onset (s)');
ylabel('Information');
figufy(gcf);
yl = ylim();
line([0 0],yl,'color','m')
axis tight