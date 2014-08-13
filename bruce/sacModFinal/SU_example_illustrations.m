clear all
close all

fit_unCor = false;

Expt_name = 'M296';
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end
Expt_num = str2num(Expt_name(2:end));
bar_ori = 0;

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

fname = 'sacStimProc';
fname = [fname sprintf('_ori%d',bar_ori)];
if fit_unCor
    fname = [fname '_unCor'];
end

anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
cd(anal_dir)
load(fname);

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/summary_figs/';
et_anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];

if strcmp(rec_type,'LP')
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
    use_nPix = 32;
elseif strcmp(rec_type,'UA')
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
    use_nPix = 16;
end

use_measured_pos = 3;
et_anal_name = 'full_eyetrack';

%if using coil initialization
if use_measured_pos == 1
    et_anal_name = [et_anal_name '_Cinit'];
end
%if using coil init with trial-sub
if use_measured_pos == 2
    et_anal_name = [et_anal_name '_CPinit'];
end
if use_measured_pos == 3
    et_anal_name = [et_anal_name '_Rinit'];
end
%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];

load([et_anal_dir et_anal_name],'et_params');


sua_sm = 0.01/dt;

%%
for cc = (n_probes+1):length(sacStimProc);
    if sacStimProc(cc).used && ~isempty(sacStimProc(cc).gsac_post_singmod)
    cur_GQM = sacStimProc(cc).ModData.rectGQM;
    flen = cur_GQM.stim_params(1).stim_dims(1);
    sp_dx = et_params.sp_dx;
    use_nPix = et_params.use_nPix;
    use_nPix_us = use_nPix*et_params.spatial_usfac;
    
    %%
    sac_xr = [-0.1 0.3];
%     sac_xr = [-0.05 0.25];
%     gr = [-1.5 3.5];
    
    close all
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
    % mmmloc(2) = mmmloc(2) - 1; %manually set this one bin earlier
    
    cur_Filts = reshape([cur_GQM.mods(:).filtK],[flen use_nPix_us length(cur_GQM.mods)]);
    filt_tempkerns = squeeze(std(cur_Filts,[],2));
    filt_Signs = [cur_GQM.mods(:).sign];
    
    eKern = mean(filt_tempkerns(:,filt_Signs==1),2);
    Ikern = mean(filt_tempkerns(:,filt_Signs==-1),2);
    t_ax = (0:(flen-1))*dt + dt/2;
    
    h1=figure;
    
    if sua_sm > 0
    sm_arate = jmm_smooth_1d_cor(sacStimProc(cc).gsac_avg_rate/dt,sua_sm);
    else
        sm_arate = sacStimProc(cc).gsac_avg_rate/dt;
    end
    
    %sac-trig avg
    subplot(3,3,1)
    plot(slags*dt,sm_arate);
    hold on; axis tight
    yl = ylim();
    line(slags(mmmloc([1 1]))*dt,yl,'color','k','linestyle','--');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k','linestyle','--');
    xlim(sac_xr);
    line(sac_xr,[0 0]+sacStimProc(cc).ModData.unit_data.avg_rate,'color','k');
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    
    %gain-kernels
    subplot(3,3,2);
    plot(slags*dt,sacStimProc(cc).gsacGainMod.off_kernel);
    hold on;
    plot(slags*dt,sacStimProc(cc).gsacGainMod.gain_kernel,'m');
    plot(slags*dt,sacStimProc(cc).gsacGainMod.stim_kernel,'k');
    plot(slags*dt,sacStimProc(cc).gsac_post_Egains,'g');
    plot(slags*dt,sacStimProc(cc).gsac_post_Igains,'r');
     plot(slags*dt,sacStimProc(cc).gsac_post_singmod.mods(3).filtK,'y');
   axis tight
    yl = ylim();
    xlim(sac_xr)
    line(slags(mmmloc([1 1]))*dt,yl,'color','k','linestyle','--');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k','linestyle','--');
    line(sac_xr,[0 0],'color','k');
    xlabel('Time (s)');
    ylabel('Filter');
    
    %TB rate map
    subplot(3,3,4)
    imagesc(slags*dt,Gtick,sacStimProc(cc).gsac_TB_rate); set(gca,'ydir','normal'); caxis([0 max(sac_dep_rate(:))]*0.65);
    yl = ylim();
% ylim(gr); yl = gr;
    line(slags(mmmloc([1 1]))*dt,yl,'color','w');
    line(slags(mmmloc([2 2]))*dt,yl,'color','w');
    xlim(sac_xr)
    xlabel('Time (s)');
    ylabel('Generating signal');
    
    %response functions
    subplot(3,3,7);
    plot(Gtick,sac_dep_rate(mmmloc(1),:)/dt,'b--'); hold on
    plot(Gtick,sacStimProc(cc).gsac_TB_rate(:,mmmloc(1))/dt,'b');
    plot(Gtick,base_out/dt,'k','linewidth',2);
    xlim(Gtick([1 end]));
% xlim(gr);
    yl = ylim();
    plot(Gtick,cur_gdist/max(cur_gdist)*yl(2)*0.75,'k--');
    xlabel('Generating signal');
    ylabel('Firing rate (Hz)');
    plot(Gtick,sac_dep_rate(mmmloc(2),:)/dt,'r--'); hold on
    plot(Gtick,sacStimProc(cc).gsac_TB_rate(:,mmmloc(2))/dt,'r');
%     ylim([0 200])
   axis tight; 
    xlim(Gtick([1 end]));
xlabel('Generating signal');
    ylabel('Firing rate (Hz)');
    
    %mod info SS
    subplot(3,3,5)
    plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo);
    hold on;
    plot(Xtick*dt,sacStimProc(cc).gsac_TB_info,'r');
    plot(slags*dt,sacStimProc(cc).gsac_modinfo,'k');
    plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo,'m');
    axis tight
    yl = ylim();
    line(slags(mmmloc([1 1]))*dt,yl,'color','k');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k');
    xlim(sac_xr)
    xlabel('Time (s)');
    ylabel('Single-spike info (bits/spike)');
    
    %mod info rate
    subplot(3,3,8)
    plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo'.*sacStimProc(cc).gsac_avg_rate);
    hold on;
    plot(Xtick*dt,sacStimProc(cc).gsac_TB_info.*sacStimProc(cc).gsac_TB_avg_rate,'r');
    plot(slags*dt,sacStimProc(cc).gsac_modinfo'.*sacStimProc(cc).gsac_avg_rate,'k');
    plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo'.*sacStimProc(cc).gsac_avg_rate,'m');
    axis tight
    yl = ylim();
    line(slags(mmmloc([1 1]))*dt,yl,'color','k');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k');
    xlim(sac_xr)
    xlabel('Time (s)');
    ylabel('Info rate (bits/sec)');
    
    subplot(3,3,3)
    plot(t_ax,eKern,t_ax,Ikern,'r');
    xlabel('Lag (s)');
    ylabel('Filter env');
    legend('Efilts','Ifilts');
    xlim([0 0.15]);
    
    
    %STA figs
    ov_sta = sacStimProc(cc).ov_phaseDep_sta;
    ov_staA = sacStimProc(cc).ov_phaseInd_sta;
    [~,sta_peakloc] = max(std(ov_sta,[],2));
    [~,staA_peakloc] = max(std(ov_staA,[],2));
    temp = sacStimProc(cc).gsac_phaseDep_sta;
    tempa = sacStimProc(cc).gsac_phaseInd_sta;
    
    sta_sm = 0.5;
    for iii = 1:size(temp,3)
        temp(:,sta_peakloc,iii) = jmm_smooth_1d_cor(temp(:,sta_peakloc,iii),sta_sm);
    end
    
    stemp = reshape(sacStimProc(cc).gsac_phaseDep_subfilt,length(slags),flen,[]);
%     stempa = reshape(sacStimProc(cc).gsac_phaseInd_subfilt,length(slags),flen,[]);
    subplot(3,3,6)
    imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(temp(:,sta_peakloc,:))');
%     imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(stemp(:,sta_peakloc,:))');
    cam = max(abs(temp(:)));
    caxis([-cam cam]);
    yl = ylim();
    xlim(sac_xr)
    line(slags(mmmloc([1 1]))*dt,yl,'color','k');
    line(slags(mmmloc([2 2]))*dt,yl,'color','k');
    xlabel('Time (s)');
    ylabel('Rel Position (deg)');
    
    subplot(3,3,9)
%     imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(tempa(:,staA_peakloc,:))');
    imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(stemp(:,sta_peakloc,:))');
    cam = max(abs(stemp(:)));
    caxis([-cam cam]);
    yl = ylim();
%     line(slags(mmmloc([1 1]))*dt,yl,'color','k');
%     line(slags(mmmloc([2 2]))*dt,yl,'color','k');
    xlim(sac_xr)
    xlabel('Time (s)');
    ylabel('Rel Position (deg)');
    
    cid = sprintf('E%d_C%d_',Expt_num,cc);
    
    %%
    fig_width = 10; rel_height = 0.8;
    figufy(h1);
    fname = [fig_dir cid sprintf('ori%d_',bar_ori) 'Gsac_mod.pdf'];
    exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(h1);
    
    %%
%     f1 = figure();
%     imagesc(ov_sta); set(gca,'ydir','normal');
%     cam = max(abs(ov_sta(:))); caxis([-cam cam]);
%     xl = xlim();
%     line(xl,[sta_peakloc sta_peakloc],'color','k');
%     figufy(f1);
%     fname = [fig_dir cid sprintf('ori%d_',bar_ori) 'ov_sta.pdf'];
%     exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     close(f1);
    
    %%
    sh = NMMdisplay_model(sacStimProc(cc).ModData.rectGQM);
    Nmods = length(sacStimProc(cc).ModData.rectGQM.mods);
    n_columns = max(round(sqrt(Nmods/2)),1);
    n_rows = ceil(Nmods/n_columns);
    figufy(sh.stim_filts);
    fig_width = 6*n_columns;
    rel_height = 0.8*n_rows/n_columns/2;
    fname = [fig_dir cid sprintf('ori%d_',bar_ori) 'stim_mod.pdf'];
    exportfig(sh.stim_filts,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(sh.stim_filts);
    
    end
end


% %%
% t_ax = (0:(flen-1))*dt + dt/2;
% p_ax = sp_dx*((1:use_nPix_us)-use_nPix_us/2);
% 
% yy = [-0.35 0.35];
% 
% ov_sta = sacStimProc(cc).ov_phaseDep_sta;
% ov_staA = sacStimProc(cc).ov_phaseInd_sta;
% % load sacStimProc_sta
% % ov_sta = sum(bsxfun(@times,all_Xmat_shift,cur_Robs))/sum(cur_Robs);
% % base_a = mean(all_Xmat_shift);
% % ov_sta = ov_sta - base_a;
% %
% % ov_staA = sum(bsxfun(@times,abs(all_Xmat_shift),cur_Robs))/sum(cur_Robs);
% % base_a = mean(abs(all_Xmat_shift));
% % ov_staA = ov_staA - base_a;
% 
% temp = sacStimProc(cc).gsac_phaseDep_sta;
% tempa = sacStimProc(cc).gsac_phaseInd_sta;
% 
% h1=figure;
% subplot(3,1,1)
% imagesc(p_ax,t_ax,reshape(ov_sta,flen,use_nPix_us));
% ylim([0 0.11]);
% xlim(yy);
% xlabel('Rel position (deg)');
% ylabel('Lag (s)');
% ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
% subplot(3,1,2);
% imagesc(p_ax,t_ax,squeeze(temp(mmmloc(1),:,:))); caxis([-cam cam]);
% xlim(yy);
% ylim([0 0.11]);
% xlabel('Rel position (deg)');
% ylabel('Lag (s)');
% subplot(3,1,3);
% imagesc(p_ax,t_ax,squeeze(temp(mmmloc(2),:,:))); caxis([-cam cam]);
% ylim([0 0.11]);
% xlim(yy);
% xlabel('Rel position (deg)');
% ylabel('Lag (s)');
% 
% % yy = [5 25];
% h2=figure;
% for jj = 1:4
%     subplot(4,1,jj);
%     imagesc(slags*dt,(1:use_nPix_us)*sp_dx-use_nPix_us*sp_dx/2,squeeze(temp(:,3+jj,:))');
%     caxis([-cam cam]);ylim(yy); yl = ylim();
%     line(slags(mmmloc([1 1]))*dt,yl,'color','k');
%     line(slags(mmmloc([2 2]))*dt,yl,'color','k');
%     xlabel('Time (s)');
%     ylabel('Rel Position (deg)');
% end
% 
% 
% cid = sprintf('E%d_C%d_',Expt_num,cc);
% 
% % fig_width = 5; rel_height = 3;
% % figufy(h1);
% % fname = [fig_dir cid 'STAexamps.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% %
% % fig_width = 5; rel_height = 4;
% % figufy(h2);
% % fname = [fig_dir cid 'STAslices.pdf'];
% % exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h2);
% 
% %%
% cur_Filts = reshape([cur_GQM.mods(:).filtK],[flen use_nPix_us length(cur_GQM.mods)]);
% filt_tempkerns = squeeze(std(cur_Filts,[],2));
% filt_Signs = [cur_GQM.mods(:).sign];
% 
% eKern = mean(filt_tempkerns(:,filt_Signs==1),2);
% Ikern = mean(filt_tempkerns(:,filt_Signs==-1),2);
% 
% t_ax = (0:(flen-1))*dt + dt/2;
% h1=figure;
% subplot(3,1,3);
% plot(t_ax,eKern,t_ax,Ikern,'r');
% xlabel('Lag (s)');
% ylabel('Filter env');
% legend('Efilts','Ifilts');
% xlim([0 0.15]);
% 
% subplot(3,1,2);
% plot(slags*dt,1+sacStimProc(cc).gsac_post_Egains); hold on
% plot(slags*dt,1+sacStimProc(cc).gsac_post_Igains,'r');
% xlabel('Time (s)');
% ylabel('Gain');
% legend('Efilts','Ifilts');
% xlim([0 0.15]);
% 
% subplot(3,1,1);
% plot(slags*dt,1+sacStimProc(cc).gsac_post_Egains); hold on
% plot(slags*dt,1+sacStimProc(cc).gsac_post_Igains,'r');
% xlabel('Time (s)');
% ylabel('Gain');
% legend('Efilts','Ifilts');
% 
% % fig_width = 5; rel_height = 2.5;
% % figufy(h1);
% % fname = [fig_dir cid 'EIkerns.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% 
% %%
% 
% h1=figure;
% plot(slags*dt,sacStimProc(cc).gsacGainMod.gain_kernel);
% hold on
% plot(slags*dt,sacStimProc(cc).gsacGainMod.stim_kernel,'k');
% plot(slags*dt,sacStimProc(cc).gsac_post_singmod.mods(3).filtK,'r');
% xlabel('Time (s)');
% ylabel('Gain kernels');
% legend('Post (both)','Pre','Post only');
% 
% % fig_width = 5; rel_height = 0.8;
% % figufy(h1);
% % fname = [fig_dir cid 'Prepost.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% 
% %%
% mod_ekern = nan(length(slags),flen);
% mod_ikern = nan(length(slags),flen);
% stim_kern = sacStimProc(cc).gsacGainMod.stim_kernel;
% gain_kern = 1+sacStimProc(cc).gsacGainMod.gain_kernel;
% for ii = 1:length(slags)
%     khist = (ii-flen+1):ii;
%     uset = find(khist > 0);
%     tkern = ones(flen,1);
%     tkern(uset) = tkern(uset) + stim_kern(khist(uset));
%     tkern = flipud(tkern);
%     
%     mod_ekern(ii,:) = (eKern.*tkern)';
%     mod_ikern(ii,:) = (Ikern.*tkern)';
% end
% 
% % mod_gkern = mod_tkern;
% % mod_gkern = bsxfun(@times,mod_gkern,gain_kern);
% 
% h1=figure;
% subplot(2,1,1)
% imagesc(slags*dt,t_ax,mod_ekern'); set(gca,'ydir','normal');
% yl = ylim();
% line([0 0],yl,'color','w');
% xlabel('Time (s)');
% ylabel('Lag (s)');
% line(slags([17 17])*dt,yl,'color','g');
% line(slags([19 19])*dt,yl,'color','r');
% 
% subplot(2,1,2)
% imagesc(slags*dt,t_ax,mod_ikern'); set(gca,'ydir','normal');
% line([0 0],yl,'color','w');
% xlabel('Time (s)');
% ylabel('Lag (s)');
% 
% h2=figure;hold on
% plot(t_ax,eKern,'k');
% plot(t_ax,mod_ekern(17,:),'g');
% plot(t_ax,mod_ekern(19,:),'r');
% xlabel('Lag (s)');
% ylabel('Temp kernel');
% 
% % fig_width = 5; rel_height = 1.6;
% % figufy(h1);
% % fname = [fig_dir cid 'tempkern_mod.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% %
% % fig_width = 5; rel_height = 0.8;
% % figufy(h2);
% % fname = [fig_dir cid 'tempkern_slices.pdf'];
% % exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h2);
% 
% %%
% h1 = figure;
% subplot(2,1,1);
% plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo); hold on
% plot(slags*dt,sacStimProc(cc).gsac_modinfo,'r')
% plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo,'k')
% plot(Xtick*dt,sacStimProc(cc).gsac_TB_info,'g');
% xlabel('Time (s)');
% ylabel('SS info (bits/spk)');
% legend('Post-mod','Pre-post mod','Subspace mod','TB mod');
% xlim(slags([1 end])*dt)
% subplot(2,1,2);
% plot(slags*dt,sacStimProc(cc).gsac_spost_modinfo.*sacStimProc(cc).gsac_avg_rate/dt); hold on
% plot(slags*dt,sacStimProc(cc).gsac_modinfo.*sacStimProc(cc).gsac_avg_rate/dt,'r')
% plot(slags*dt,sacStimProc(cc).gsac_sub_modinfo.*sacStimProc(cc).gsac_avg_rate/dt,'k')
% plot(Xtick*dt,sacStimProc(cc).gsac_TB_info.*sacStimProc(cc).gsac_TB_avg_rate/dt,'g');
% xlabel('Time (s)');
% ylabel('SS info rate (bits/sec)');
% legend('Post-mod','Pre-post mod','Subspace mod','TB mod');
% xlim(slags([1 end])*dt)
% 
% % fig_width = 5; rel_height = 1.6;
% % figufy(h1);
% % fname = [fig_dir cid 'allinfo.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% 
% %%
% tempsub = reshape(sacStimProc(cc).gsac_phaseDep_subfilt,[length(slags) flen use_nPix_us]);
% tempasub = reshape(sacStimProc(cc).gsac_phaseInd_subfilt,[length(slags) flen use_nPix_us]);
% cam = max(abs(tempsub(:)));
% cam2 = max(abs(tempasub(:)));
% 
% h1=figure;
% subplot(3,2,1)
% imagesc(p_ax,t_ax,squeeze(mean(tempsub(1:5,:,:)))); caxis([-cam cam]);
% ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% subplot(3,2,3);
% imagesc(p_ax,t_ax,squeeze(tempsub(mmmloc(1),:,:))); caxis([-cam cam]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% subplot(3,2,5);
% imagesc(p_ax,t_ax,squeeze(tempsub(mmmloc(2),:,:))); caxis([-cam cam]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% 
% subplot(3,2,2)
% imagesc(p_ax,t_ax,squeeze(mean(tempasub(1:5,:,:)))); caxis([-cam cam]);
% ca = caxis(); cam2 = max(abs(ca)); caxis([-cam2 cam2]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% subplot(3,2,4);
% imagesc(p_ax,t_ax,squeeze(tempasub(mmmloc(1),:,:)));  caxis([-cam2 cam2]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% subplot(3,2,6);
% imagesc(p_ax,t_ax,squeeze(tempasub(mmmloc(2),:,:)));  caxis([-cam2 cam2]);
% xlabel('Rel Pos (deg)');
% ylabel('Lag (s)');
% 
% 
% h2=figure;
% subplot(2,1,1);
% imagesc(slags*dt,p_ax,squeeze(tempsub(:,5,:))');
% caxis([-cam cam]);ylim(yy); yl = ylim();
% line(slags(mmmloc([1 1]))*dt,yl,'color','k');
% line(slags(mmmloc([2 2]))*dt,yl,'color','k');
% xlabel('Time (s)');
% ylabel('Rel Pos (deg)');
% subplot(2,1,2);
% imagesc(slags*dt,p_ax,squeeze(temp(:,5,:))');
% caxis([-cam cam]);ylim(yy); yl = ylim();
% line(slags(mmmloc([1 1]))*dt,yl,'color','k');
% line(slags(mmmloc([2 2]))*dt,yl,'color','k');
% xlabel('Time (s)');
% ylabel('Rel Pos (deg)');
% 
% 
% 
% % fig_width = 8; rel_height = 1.2;
% % figufy(h1);
% % fname = [fig_dir cid 'Subspace_filts.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);
% %
% % fig_width = 6; rel_height = 1.5;
% % figufy(h2);
% % fname = [fig_dir cid 'filtslice_compare.pdf'];
% % exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h2);
% 
% %%
% 
% sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
% g_tot = stimG;
% clear tr_stim
% tr_stim{1} = [g_tot];
% tr_stim{2} = cur_Xsac;
% clear sac_stim_params
% sac_stim_params(1) = NMMcreate_stim_params(1);
% sac_stim_params(2) = NMMcreate_stim_params([size(cur_Xsac,2)]);
% mod_signs = [1 1];
% Xtargets = [1 2];
% NL_types = {'lin','lin','lin'};
% post_gsac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
% post_gsac_Smod = NMMfit_filters(post_gsac_Smod,cur_Robs,tr_stim,[],[],silent);
% 
% 
% %%
% is_X_bs = false(size(all_Xmat_shift));
% temp = repmat((1:1:(flen))',1,use_nPix_us);
% back_mat = create_time_embedding(temp,NMMcreate_stim_params([flen use_nPix_us]));
% % uuu = back_mat > 0;
% uuu = back_mat == 0;
% for ss = 1:length(big_sacs)
%     fprintf('%d/%d\n',ss,length(big_sacs));
%     cur_stop_ind = saccade_stop_inds(big_sacs(ss));
%     %     cur_stop_ind = saccade_start_inds(big_sacs(ss));
%     cur_chunk = cur_stop_ind:(cur_stop_ind + (flen-1));
%     is_X_bs(cur_chunk,:) = uuu;
% end
% 
% shuf_stim = all_shift_stimmat_up;
% shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
% shuf_X = create_time_embedding(shuf_stim,stim_params_us);
% shuf_X = shuf_X(used_inds,use_kInds_up);
% shuf_X(~is_X_bs) = all_Xmat_shift(~is_X_bs);
% 
% %%
% sacGainMod = sacStimProc(cc).gsacGainMod;
% [gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
% [shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
% % null_prate = ones(size(gain_pred_rate))*mean(cur_Robs);
% % null_LLseq = cur_Robs.*log2(null_prate)-null_prate;
% 
% % str_stim = tr_stim;
% % [~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, cur_Robs, shuf_X);
% % fgint = bsxfun(@times,fgint,[cur_GQM.mods(:).sign]);
% % sg_tot = sum(fgint,2);
% % str_stim{1} = sg_tot;
% 
% % % [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% % % [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
% % [base_LL,~,base_predrate] = NMMmodel_eval(post_gsac_Smod,cur_Robs,tr_stim);
% % [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(post_gsac_Smod,cur_Robs,str_stim);
% % base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
% % baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;
% 
% [ta_gain,tlags] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,saccade_stop_inds(big_sacs),14,14);
% % [tabase_gain,tlags] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,saccade_stop_inds(big_sacs),14,14);
% 
% [nspks,tlags] = get_event_trig_avg_v3(cur_Robs,saccade_stop_inds(big_sacs),14,14);
% tot_nspks = nspks*length(big_sacs);
% ta_gain = ta_gain./tot_nspks;
% % tabase_gain = tabase_gain./tot_nspks;
% 
% 
% % [sac_fpost_info,sac_spost_info,sac_subpost_info,sac_info,sac_LL,sac_fpost_LL,sac_spost_LL,sac_subpost_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
% % for ii = 1:length(slags)
% %     temp = find(cur_Xsac(:,ii) == 1);
% %     sac_avgrate(ii) = mean(cur_Robs(temp));
% %     cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
% %     sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
% %     sac_Nspks(ii) = sum(cur_Robs(temp));
% %
% %     [shuf_LL(ii)] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% %     n_frames(ii) = sum(t_since_sac_start == slags(ii));
% % end
% 
% %%
% 
% % [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% % cur_Robs = poissrnd(base_predrate);
% 
% msac_times = saccade_start_inds(micro_sacs);
% msac_isis = [0;diff(msac_times)];
% bad = find(msac_isis < 40);
% use_micro_sacs = micro_sacs;
% use_micro_sacs(bad) = [];
% 
% cur_sac_start_inds = saccade_start_inds(big_sacs);
% cur_sac_stop_inds = saccade_stop_inds(big_sacs);
% % cur_sac_start_inds = saccade_start_inds(use_micro_sacs);
% % cur_sac_stop_inds = saccade_stop_inds(use_micro_sacs);
% 
% % rand_jit = randi(100,length(big_sacs),1);
% % cur_sac_start_inds = cur_sac_start_inds + rand_jit;
% % cur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
% % % cur_sac_stop_inds = cur_sac_start_inds;
% % bad = find(cur_sac_stop_inds > NT);
% % cur_sac_start_inds(bad) = [];
% % cur_sac_stop_inds(bad) = [];
% 
% lagdist = 15;
% all_before = [];
% all_after = [];
% all_during = [];
% for ss = 1:length(cur_sac_start_inds)
%     cur_before = (cur_sac_start_inds(ss)-lagdist-10):(cur_sac_start_inds(ss)-1);
%     cur_after = cur_sac_stop_inds(ss):(cur_sac_stop_inds(ss)+lagdist);
%     cur_during = (cur_sac_start_inds(ss)):(cur_sac_stop_inds(ss)-1);
%     %     all_before = cat(1,all_before,find(ismember(cc_uinds,cur_before)));
%     %     all_after = cat(1,all_after,find(ismember(cc_uinds,cur_after)));
%     %     all_during = cat(1,all_during,find(ismember(cc_uinds,cur_during)));
%     
%     cur_before_tvec = all_trialvec(used_inds(cur_before));
%     cur_before_tvec(cur_before_tvec ~= all_trialvec(used_inds(cur_sac_start_inds(ss)))) = [];
%     cur_after_tvec = all_trialvec(used_inds(cur_after));
%     cur_after_tvec(cur_after_tvec ~= all_trialvec(used_inds(cur_sac_stop_inds(ss)))) = [];
%     cur_during_tvec = all_trialvec(used_inds(cur_during));
%     cur_during_tvec(cur_during_tvec ~= all_trialvec(used_inds(cur_sac_start_inds(ss)))) = [];
%     
%     all_before = cat(2,all_before,cur_before);
%     all_after = cat(2,all_after,cur_after);
%     all_during = cat(2,all_during,cur_during);
% end
% 
% %% BEFORE
% shuf_stim = all_shift_stimmat_up;
% shuf_stim(used_inds(all_before),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_before),1)),:);
% shuf_X = create_time_embedding(shuf_stim,stim_params_us);
% shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
% 
% sacGainMod = sacStimProc(cc).gsacGainMod;
% % sacGainMod = sacStimProc(cc).msacGainMod;
% [gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
% 
% 
% [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% 
% temp_sparams = NMMcreate_stim_params(1);
% tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
% tempmod = NMMfit_filters(tempmod,cur_Robs,Gshuff);
% [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs,Gshuff);
% 
% gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% % shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
% shuf_LLseq = cur_Robs.*log2(tempprate)-tempprate;
% 
% [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
% base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
% baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;
% 
% % cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
% cur_stop_inds = find(ismember(cc_uinds,cur_sac_stop_inds));
% % [before_ta_gain,tlags] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,[],all_trialvec(used_inds(cc_uinds)));
% 
% [before_ta_gain,tlags,~,n_events,before_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% before_ta_gain = nansum(before_ta_mat);
% 
% [before_base_gain,tlags,~,~,before_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% before_base_gain = nansum(before_base_mat);
% %% AFTER
% shuf_stim = all_shift_stimmat_up;
% shuf_stim(used_inds(all_after),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_after),1)),:);
% shuf_X = create_time_embedding(shuf_stim,stim_params_us);
% shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
% 
% sacGainMod = sacStimProc(cc).gsacGainMod;
% % sacGainMod = sacStimProc(cc).msacGainMod;
% [gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
% [shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
% 
% [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
% base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
% baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;
% 
% % cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
% [after_ta_gain,tlags,~,~,after_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% after_ta_gain = nansum(after_ta_mat);
% [after_base_gain,tlags,~,~,after_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% after_base_gain = nansum(after_base_mat);
% %% DURING
% shuf_stim = all_shift_stimmat_up;
% shuf_stim(used_inds(all_during),:) = all_shift_stimmat_up(used_inds(randi(NT,length(all_during),1)),:);
% shuf_X = create_time_embedding(shuf_stim,stim_params_us);
% shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
% 
% sacGainMod = sacStimProc(cc).gsacGainMod;
% % sacGainMod = sacStimProc(cc).msacGainMod;
% [gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
% [shuf_LL,shuf_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% shuf_LLseq = cur_Robs.*log2(shuf_pred_rate)-shuf_pred_rate;
% 
% [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% [base_shufLL,~,base_shuff_prate] = NMMmodel_eval(cur_GQM,cur_Robs,shuf_X);
% base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
% baseshuf_LLseq = cur_Robs.*log2(base_shuff_prate)-base_shuff_prate;
% 
% % cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
% [during_ta_gain,tlags,~,~,during_ta_mat] = get_event_trig_avg_v3(gain_LLseq-shuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% during_ta_gain = nansum(during_ta_mat);
% [during_base_gain,tlags,~,~,during_base_mat] = get_event_trig_avg_v3(base_LLseq-baseshuf_LLseq,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% during_base_gain = nansum(during_base_mat);
% %%
% % cur_stop_inds = find(ismember(cc_uinds,saccade_stop_inds(big_sacs)));
% [nspks,tlags,~,~,nspks_mat] = get_event_trig_avg_v3(cur_Robs,cur_stop_inds,lagdist,lagdist,2,all_trialvec(used_inds(cc_uinds)));
% tot_nspks = nansum(nspks_mat);
% 
% before_ta_Ngain = before_ta_gain./tot_nspks;
% before_base_Ngain = before_base_gain./tot_nspks;
% after_ta_Ngain = after_ta_gain./tot_nspks;
% after_base_Ngain = after_base_gain./tot_nspks;
% during_ta_Ngain = during_ta_gain./tot_nspks;
% during_base_Ngain = during_base_gain./tot_nspks;
% %
% % before_ta_Ngainr = before_ta_Ngain.*length(cur_stop_inds);
% % before_base_Ngainr = before_base_Ngain.*length(cur_stop_inds);
% % after_ta_Ngainr = after_ta_Ngain.*length(cur_stop_inds);
% % after_base_Ngainr = after_base_Ngain.*length(cur_stop_inds);
% % during_ta_Ngainr = during_ta_Ngain.*length(cur_stop_inds);
% % during_base_Ngainr = during_base_Ngain.*length(cur_stop_inds);
% 
% %%
% figure();
% plot(tlags*dt,before_ta_Ngain); hold on
% plot(tlags*dt,after_ta_Ngain,'r');
% plot(tlags*dt,during_ta_Ngain,'k');
% 
% 
% plot(tlags*dt,before_base_Ngain,'--');hold on
% plot(tlags*dt,after_base_Ngain,'r--')
% plot(tlags*dt,during_base_Ngain,'k--');
% 
% xlabel('Time since fixation onset (s)');
% ylabel('Information');
% figufy(gcf);
% yl = ylim();
% line([0 0],yl,'color','m')
% axis tight
% 
% 
% 
% 
% %%
% h2=figure;
% subplot(2,1,1);
% imagesc(slags*dt,Gtick,sac_dep_rate'); set(gca,'ydir','normal'); caxis([0 max(sac_dep_rate(:))]*0.9);
% % imagesc(slags*dt,Gtick,log10(sac_dep_rate')); set(gca,'ydir','normal'); caxis(log10([0.01 max(sac_dep_rate(:))]));
% yl = ylim();
% line(slags(mmmloc([1 1]))*dt,yl,'color','w');
% line(slags(mmmloc([2 2]))*dt,yl,'color','w');
% xlabel('Time (s)');
% ylabel('Generating signal');
% subplot(2,1,2)
% imagesc(slags*dt,Gtick,sacStimProc(cc).gsac_TB_rate); set(gca,'ydir','normal'); caxis([0 max(sac_dep_rate(:))]*0.9);
% % imagesc(slags*dt,Gtick,log10(sacStimProc(cc).gsac_TB_rate)); set(gca,'ydir','normal');  caxis(log10([0.01 max(sac_dep_rate(:))]));
% yl = ylim();
% line(slags(mmmloc([1 1]))*dt,yl,'color','w');
% line(slags(mmmloc([2 2]))*dt,yl,'color','w');
% xlabel('Time (s)');
% ylabel('Generating signal');
