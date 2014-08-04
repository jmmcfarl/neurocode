clear all

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_expts = length(Expt_list);

fig_dir =  '/home/james/Analysis/bruce/saccade_modulation/sacmod_summary_figs/';

xr = [-0.2 0.4];
close all

flen = 16;
nPix = 32;
all_SacMod = [];
for ee = 1:n_expts
    Expt_name = Expt_list{ee};
    
    out_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod/'];
    fin_anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
    
    cd(fin_anal_dir);
    fin_mod_name = 'corr_mods_hbar';
    load(fin_mod_name);
    
    cd(out_dir)
%     dname = 'sacmod_stimproc';
    dname = 'sacmod_stimproc';
    load(dname);
    full_n_chs = length(SacMod);
    %%
    uset = 97:full_n_chs;
    uset = uset(arrayfun(@(x)length(x.tot_spks),cor_gqm(uset)) > 0);
    uset = uset([cor_gqm(uset).tot_samps] > 0);
    uset([SacMod(uset).avg_rate] == 0) = [];
    n_units = length(uset);
    fprintf('Using %d units\n',n_units);
    
    all_SacMod = cat(2,all_SacMod,SacMod(uset));
    %%
%     for uu = 1:n_units
%         cur_unit_id = uset(uu);
%         stim_filters = reshape([cor_gqm(cur_unit_id).mod_fit.mods(1:6).filtK],[flen nPix 6]);
%         ford = [1 3 4 2 5 6];
%         slocs = [1 2 3 10 11 12];
%         mg = [0.04 0.06];
%         mg2 = [0.04 0.02];
%         
%         h = figure();
%         for ss = 1:6
%            subplot_tight(6,9,slocs(ss),mg2);
%            imagesc(squeeze(stim_filters(:,:,ford(ss))));
%            ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]); colormap(gray);
%            set(gca,'ydir','normal','xtick',[],'ytick',[]);
%         end
%             
%         subplot_tight(6,9,[19 20 21 28 29 30],mg); hold on
%         plot(sac_bincents*dt,SacMod(cur_unit_id).sac_avg_prate/dt);
%         plot(sac_bincents*dt,SacMod(cur_unit_id).msac_avg_prate/dt,'r');
%         axis tight
%         xlim(xr);
%         line(xr,SacMod(cur_unit_id).avg_rate([1 1])/dt,'color','k','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Average rate (Hz)');
%         
%         subplot_tight(6,9,[37 38 39 46 47 48],mg); hold on;
%         plot(sac_bincents*dt,SacMod(cur_unit_id).sac_offset);
%         plot(sac_bincents*dt,SacMod(cur_unit_id).msac_offset,'r');
%         axis tight
%         xlim(xr);
%         line(xr,[0 0],'color','k','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Sac offset');
%     
%         subplot_tight(6,9,[4 5 6 13 14 15],mg); hold on;
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).sac_Egain,'g');
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).sac_Igain,'k');
%         axis tight
%         xlim(xr);
%         line(xr,[1 1],'color','b','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Sac gain');
% 
%         subplot_tight(6,9,[22 23 24 31 32 33],mg); hold on;
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).msac_Egain,'g');
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).msac_Igain,'k');
%         axis tight
%         xlim(xr);
%         line(xr,[1 1],'color','b','linestyle','--');        
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('MSac gain');
%         
%         subplot_tight(6,9,[40 41 42 49 50 51],mg); hold on;
%         plot(sac_bincents*dt,SacMod(cur_unit_id).avg_sac_info);
%         plot(sac_bincents*dt,SacMod(cur_unit_id).avg_msac_info,'r');
%         axis tight
%         xlim(xr);
%         line(xr,SacMod(cur_unit_id).ov_info([1 1]),'color','k','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Stim info');
%         
%         subplot_tight(6,9,[7 8 9 16 17 18],mg); hold on;
%         plot(SacMod(cur_unit_id).pop_filter,'.-'); axis tight
%         axis tight
%         xl = xlim();
%         line(xl,[0 0],'color','k','linestyle','--');
%         ylabel('Coupling');
%        
%         subplot_tight(6,9,[25 26 27 34 35 36],mg); hold on;
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).pop_sac_gain,'b');
%         plot(sac_bincents*dt,1+SacMod(cur_unit_id).pop_msac_gain,'r');
%         axis tight
%         xlim(xr);
%         line(xr,[1 1],'color','k','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Pop-coupling gain');
%        
%         subplot_tight(6,9,[43 44 45 52 53 54],mg); hold on;
%         plot(sac_bincents*dt,SacMod(cur_unit_id).sac_trig_pop_info,'b');
%         plot(sac_bincents*dt,SacMod(cur_unit_id).msac_trig_pop_info,'r');
%         axis tight
%         xlim(xr);
%         line(xr,SacMod(cur_unit_id).ov_pop_info([1 1]),'color','k','linestyle','--');
%         yl = ylim();
%         line([0 0],yl,'color','k','linestyle','--');
%         ylabel('Pop-coupling info');
% 
%         fig_width = 10;
%         rel_height = 0.6;
%         figufy(h);
% set(gcf, 'PaperSize', [12*0.1 12]);
% % figname = [fig_dir sprintf('%s_Unit%d',Expt_name,cur_unit_id)];
% figname = [fig_dir sprintf('%s_Unit%d_orthmag',Expt_name,cur_unit_id)];
%         exportfig(h,figname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%         print(figname,'-dpng');
%         close(h);
%         
%     end
   
end

%%
close all
min_stim_info = 0.1;
min_avg_rate = 5;

tot_units = length(all_SacMod);
all_stim_info = [all_SacMod(:).ov_info];
all_stim_infor = [all_SacMod(:).ov_info_rate];
all_pop_info = [all_SacMod(:).ov_pop_info];
all_pop_infor = [all_SacMod(:).ov_pop_info_rate];
ov_avg_rate = [all_SacMod(:).avg_rate];

% uset = find(all_stim_info >= min_stim_info);
uset = find(ov_avg_rate/dt >= min_avg_rate);

all_sac_rate = [all_SacMod(:).sac_avg_prate];
all_msac_rate = [all_SacMod(:).msac_avg_prate];
all_sac_rate = bsxfun(@rdivide,all_sac_rate,ov_avg_rate);
all_msac_rate = bsxfun(@rdivide,all_msac_rate,ov_avg_rate);

all_sac_info = [all_SacMod(:).avg_sac_info];
all_msac_info = [all_SacMod(:).avg_msac_info];
all_sac_info = bsxfun(@rdivide,all_sac_info,all_stim_info);
all_msac_info = bsxfun(@rdivide,all_msac_info,all_stim_info);

all_sac_infor = [all_SacMod(:).avg_sac_info_rate];
all_msac_infor = [all_SacMod(:).avg_msac_info_rate];
all_sac_infor = bsxfun(@rdivide,all_sac_infor,all_stim_infor);
all_msac_infor = bsxfun(@rdivide,all_msac_infor,all_stim_infor);

h=figure; 
subplot(3,1,1);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_rate(:,uset),2),std(all_sac_rate(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_rate(:,uset),2),std(all_msac_rate(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
line([-0.2 0.4],[1 1],'color','k');
yl = ylim(); line([0 0],yl,'color','k');
title('Avg rate');

subplot(3,1,2);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_info(:,uset),2),std(all_sac_info(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_info(:,uset),2),std(all_msac_info(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
xlim([-0.2 0.4]);
line([-0.2 0.4],[1 1],'color','k');
ylim([0.5 1.5]);
yl = ylim(); line([0 0],yl,'color','k');
title('Stimulus info');

subplot(3,1,3);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_infor(:,uset),2),std(all_sac_infor(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_infor(:,uset),2),std(all_msac_infor(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
xlim([-0.2 0.4]);
line([-0.2 0.4],[1 1],'color','k');
ylim([0.5 1.5]);
yl = ylim(); line([0 0],yl,'color','k');
title('Stimulus info rate');

% ename = '/home/james/Desktop/lab_meeting_figs/sac_infomod_sizedep.png';
% fig_width = 4; rel_height = 2;
% figufy(h);
% exportfig(h,ename,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
all_sac_pinfo = [all_SacMod(:).sac_trig_pop_info];
all_msac_pinfo = [all_SacMod(:).msac_trig_pop_info];
all_sac_pinfo = bsxfun(@rdivide,all_sac_pinfo,all_pop_info);
all_msac_pinfo = bsxfun(@rdivide,all_msac_pinfo,all_pop_info);

all_sac_pinfor = [all_SacMod(:).sac_trig_pop_info_rate];
all_msac_pinfor = [all_SacMod(:).msac_trig_pop_info_rate];
all_sac_pinfor = bsxfun(@rdivide,all_sac_pinfor,all_pop_infor);
all_msac_pinfor = bsxfun(@rdivide,all_msac_pinfor,all_pop_infor);

h = figure; 
subplot(3,1,1);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_rate(:,uset),2),std(all_sac_rate(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_rate(:,uset),2),std(all_msac_rate(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
line([-0.2 0.4],[1 1],'color','k');
yl = ylim(); line([0 0],yl,'color','k');
title('Avg rate');

subplot(3,1,2);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_pinfo(:,uset),2),std(all_sac_pinfo(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_pinfo(:,uset),2),std(all_msac_pinfo(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
ylim([0.5 2.5])
line([-0.2 0.4],[1 1],'color','k');
yl = ylim(); line([0 0],yl,'color','k');
title('Network info');

subplot(3,1,3);
hold on
shadedErrorBar(sac_bincents*dt,mean(all_sac_pinfor(:,uset),2),std(all_sac_pinfor(:,uset),[],2)/sqrt(length(uset)),{'color','r'});
shadedErrorBar(sac_bincents*dt,mean(all_msac_pinfor(:,uset),2),std(all_msac_pinfor(:,uset),[],2)/sqrt(length(uset)),{'color','b'});
xlim([-0.2 0.4]);
ylim([0.5 2.5])
line([-0.2 0.4],[1 1],'color','k');
yl = ylim(); line([0 0],yl,'color','k');
title('Network info rate');


ename = '/home/james/Desktop/lab_meeting_figs/pop_infomod.pdf';
fig_width = 4; rel_height = 2;
figufy(h);
exportfig(h,ename,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
