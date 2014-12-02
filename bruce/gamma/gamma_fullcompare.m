clear all
close all

all_avgs = [];
all_stds = [];
all_cond_spec = [];
all_cond_lfp = [];
all_cond_names = [];
all_exptnum = [];
%% G091
G091_data = '~/Analysis/bruce/G091/gamma/cond_trig_avg_gamma.mat';
load(G091_data);

all_cond_spec = cat(2,all_cond_spec,cond_trig_spec);
all_cond_lfp = cat(2,all_cond_lfp,cond_trig_lfp);
all_avgs = cat(3,all_avgs,ov_avg_ampgram);
all_stds = cat(3,all_stds,ov_std_ampgram);
all_C{1} = C;
all_cond_names{1} = cond_spec;
all_exptnum = [all_exptnum; 1*ones(size(cond_trig_spec,2),1)];

clear cond_trig_spec ov* cond_spec
%% G091 ori test
G091_data = '~/Analysis/bruce/G091/gamma/ori_cond_trig_avg_gamma.mat';
load(G091_data);
cond_spec = cell(2,1);
cond_spec{1} = 'ori'; cond_spec{2} = 'nph';
all_cond_spec = cat(2,all_cond_spec,cond_trig_spec);
all_cond_lfp = cat(2,all_cond_lfp,cond_trig_lfp);
all_avgs = cat(3,all_avgs,ov_avg_ampgram);
all_stds = cat(3,all_stds,ov_std_ampgram);
all_C{2} = C;
all_cond_names{2} = cond_spec;
all_exptnum = [all_exptnum; 2*ones(size(cond_trig_spec,2),1)];
clear cond_trig_spec ov*

%% G092
G092_data = '~/Analysis/bruce/G092/gamma/cond_trig_avg_gamma.mat';
load(G092_data);

all_cond_spec = cat(2,all_cond_spec,cond_trig_spec);
all_cond_lfp = cat(2,all_cond_lfp,cond_trig_lfp);
all_avgs = cat(3,all_avgs,ov_avg_ampgram);
all_stds = cat(3,all_stds,ov_std_ampgram);
all_C{3} = C;
all_cond_names{3} = cond_spec;
all_exptnum = [all_exptnum; 3*ones(size(cond_trig_spec,2),1)];
clear cond_trig_spec ov*

%% G093 RDS
G093_data = '~/Analysis/bruce/G093/gamma/rdscond_trig_avg_gamma.mat';
load(G093_data);

all_cond_spec = cat(2,all_cond_spec,cond_trig_spec);
all_cond_lfp = cat(2,all_cond_lfp,cond_trig_lfp);
all_avgs = cat(3,all_avgs,ov_avg_ampgram);
all_stds = cat(3,all_stds,ov_std_ampgram);
all_C{4} = C;
all_cond_names{4} = 'sl';
all_exptnum = [all_exptnum; 4*ones(size(cond_trig_spec,2),1)];
clear cond_trig_spec ov*

%% G093 GRATINGS
% G093_data = '~/Analysis/bruce/G093/gamma/grateev_cond_trig_avgs.mat';
% load(G093_data);
% 
% % all_cond_spec = cat(2,all_cond_spec,cond_trig_spec);
% all_avgs = cat(3,all_avgs,ov_avg_ampgram);
% all_stds = cat(3,all_stds,ov_std_ampgraccm);
% all_C{5} = C;
% % all_cond_names{5} = cond_spec;
% % all_exptnum = [all_exptnum; 1*ones(size(cond_trig_spec,2),1)];

%% FIT POWER LAW FUNCTION TO SPECTRA
all_avg_spectrum = squeeze(mean(all_avgs,3));
all_std_spectrum = squeeze(mean(all_stds,3));

% % powfun = @(a,x)(a(1)*x.^a(2));
% powfun = @(a,x)(a(1)*exp(x*a(2)));
% 
% un_wfreqs = linspace(wfreqs(1),wfreqs(end),100)';
% avg_pow_int = interp1(wfreqs,all_avg_spectrum',un_wfreqs);
% std_pow_int = interp1(wfreqs,all_std_spectrum',un_wfreqs);
% 
% for ll = 1:size(all_avgs,1)
%     BETA_0 = [max(all_avg_spectrum(ll,:)) -1];
%     BETA_avg(ll,:) = nlinfit(un_wfreqs,avg_pow_int(:,ll),powfun,BETA_0);
%     powfit = un_wfreqs.^BETA_avg(ll,2)*BETA_avg(ll,1);
%     interp_avgpow(ll,:) = interp1(un_wfreqs,powfit,wfreqs);
%     
%     BETA_std(ll,:) = nlinfit(un_wfreqs,std_pow_int(:,ll),powfun,[max(all_std_spectrum(ll,:)) -1]);
%     powfit = un_wfreqs.^BETA_std(ll,2)*BETA_std(ll,1);
%     interp_stdpow(ll,:) = interp1(un_wfreqs,powfit,wfreqs);
% end

un_wfreqs = linspace(wfreqs(1),wfreqs(end),100)';
avg_pow_int = interp1(wfreqs,log10(all_avg_spectrum'),un_wfreqs);
std_pow_int = interp1(wfreqs,log10(all_std_spectrum'),un_wfreqs);

for ll = 1:size(all_avgs,1)
    BETA_avg(ll,:) = robustfit(un_wfreqs,avg_pow_int(:,ll));
    powfit = un_wfreqs*BETA_avg(ll,2)+BETA_avg(ll,1);
    interp_avgpow(ll,:) = interp1(un_wfreqs,powfit,wfreqs);
    
    BETA_std(ll,:) = robustfit(un_wfreqs,std_pow_int(:,ll));
    powfit = un_wfreqs*BETA_std(ll,2)+BETA_std(ll,1);
    interp_stdpow(ll,:) = interp1(un_wfreqs,powfit,wfreqs);
end
interp_avgpow = 10.^(interp_avgpow);
interp_stdpow = 10.^(interp_stdpow);
norm_cond_spec = bsxfun(@minus,all_cond_spec,reshape(interp_avgpow,[length(use_lfps),1,1,length(wfreqs)]));
norm_cond_spec = bsxfun(@rdivide,norm_cond_spec,reshape(interp_stdpow,[length(use_lfps),1,1,length(wfreqs)]));


znorm_cond_spec = bsxfun(@minus,all_cond_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
znorm_cond_spec = bsxfun(@rdivide,znorm_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));

avg_norm_cond_spec = squeeze(mean(norm_cond_spec));
avg_znorm_cond_spec = squeeze(mean(znorm_cond_spec));
avg_cond_lfp = squeeze(mean(all_cond_lfp));
sem_cond_lfp = squeeze(std(all_cond_lfp))/sqrt(length(use_lfps));
% sem_cond_lfp = squeeze(std(all_cond_lfp));

% cd ~/Analysis/bruce/ 
% save gamma_powfits interp* wfreqs use_lfps all_*_spectrum


%% VISUALIZE G091 ori
close all
cur_conds = find(all_exptnum == 2);
tot_nconds = length(cur_conds);
for cc = 1:tot_nconds
    for nn = 1:length(all_cond_names{2})
        fprintf('%s = %.3f ',all_cond_names{2}{nn},all_C{2}(cc,nn));
    end
    fprintf('\n');
    
    pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cur_conds(cc),:,:))');shading flat
    caxis([-1 1]);
%     colorbar
    pause
    clf
end

%% VISUALIZE G091
close all
cur_conds = find(all_exptnum == 1);
tot_nconds = length(cur_conds);
for cc = 1:tot_nconds
    for nn = 1:length(all_cond_names{1})
        fprintf('%s = %.3f ',all_cond_names{1}{nn},all_C{1}(cc,nn));
    end
    fprintf('\n');
    subplot(2,1,1)
    shadedErrorBar(lags/Fsd,avg_cond_lfp(cur_conds(cc),:),sem_cond_lfp(cur_conds(cc),:));
    xlabel('Time since trial onset (s)','fontsize',14);
    ylabel('LFP amplitude (V)','fontsize',14);
    xlim([-0.1 1]);
    subplot(2,1,2)
    pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cur_conds(cc),:,:))');shading flat
    caxis([-1 1]);
%     colorbar
    ylim([2 100]);
    xlabel('Time since trial onset (s)','fontsize',14);
    ylabel('Frequency (Hz)','fontsize',14);
     xlim([-0.1 1]);
   pause
    clf
end

%% VISUALIZE G092
cd ~/Analysis/bruce/gamma/
xr = [-0.1 1];
close all
figure
cur_conds = find(all_exptnum == 3);
tot_nconds = length(cur_conds);
for cc = 1:tot_nconds
fname = sprintf('G092_specgramZOOM_cond%d',cc);
    for nn = 1:length(all_cond_names{3})
        fprintf('%s = %.3f ',all_cond_names{3}{nn},all_C{3}(cc,nn));
    end
    fprintf('\n');
    cur_title = [];
    cur_title = [cur_title sprintf('%s = %d  ',all_cond_names{3}{1},all_C{3}(cc,1))];
    cur_title = [cur_title sprintf('%s = %d  ',all_cond_names{3}{2},all_C{3}(cc,2))];
    cur_title = [cur_title sprintf('%s = %d  ',all_cond_names{3}{3},all_C{3}(cc,3))];
    cur_title = [cur_title sprintf('%s = %.2f  ',all_cond_names{3}{4},all_C{3}(cc,4))];
    cur_title = [cur_title sprintf('%s = %d  ',all_cond_names{3}{5},all_C{3}(cc,5))];
     cur_title = [cur_title sprintf('%s = %d  ',all_cond_names{3}{6},all_C{3}(cc,6))];

   subplot(2,1,1)
    hold on
    shadedErrorBar(lags/Fsd,avg_cond_lfp(cur_conds(cc),:),sem_cond_lfp(cur_conds(cc),:),{'color','r'});
    xlabel('Time since trial onset (s)','fontsize',12);
    ylabel('LFP amplitude (V)','fontsize',10);
    ylim([-20 10]*1e-5);
    xlim(xr);
   title(cur_title,'fontsize',12,'fontname','arial');
    subplot(2,1,2)
    pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cur_conds(cc),:,:))');shading interp
    caxis([-1 1]);
    %     colorbar
    ylim([2 120]);
    xlabel('Time since trial onset (s)','fontsize',12);
    ylabel('Frequency (Hz)','fontsize',12);
    xlim(xr);
   
    stimes = [0:0.24:4];
    subplot(2,1,1);
    yl = ylim();
    for ii = 1:length(stimes)
        line(stimes([ii ii]),yl,'color','k')
    end
    subplot(2,1,2);
    yl = ylim();
    for ii = 1:length(stimes)
        line(stimes([ii ii]),yl,'color','w')
    end
%     fillPage(gcf,'papersize',[8 8]);
%     print(fname,'-dpng')
%     if cc == 1
%     print(fname,'-dpsc')
%     else
%     print(fname,'-dpsc','-append')
%     end
    pause
    clf
end


%% VISUALIZE G093 RDS
close all
cur_conds = find(all_exptnum == 4);
tot_nconds = length(cur_conds);
for cc = 1:tot_nconds
    fprintf('sl = %.3f\n ',all_C{4}(cc));
    
    subplot(2,1,1)
    hold on
    shadedErrorBar(lags/Fsd,avg_cond_lfp(cur_conds(cc),:),sem_cond_lfp(cur_conds(cc),:),{'color','r'});
    xlabel('Time since trial onset (s)','fontsize',14);
    ylabel('LFP amplitude (V)','fontsize',14);
    xlim([-0.1 1]);
    subplot(2,1,[2])
    pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cur_conds(cc),:,:))');shading flat
    caxis([-0.75 0.75]);
%     colorbar
    ylim([2 100]);
    xlabel('Time since trial onset (s)','fontsize',14);
    ylabel('Frequency (Hz)','fontsize',14);
     xlim([-0.1 1]);
     stimes = [0:0.48:4];
subplot(2,1,1);
    yl = ylim();
for ii = 1:length(stimes)
    line(stimes([ii ii]),yl,'color','r')
end
subplot(2,1,2);
    yl = ylim();
for ii = 1:length(stimes)
    line(stimes([ii ii]),yl,'color','w')
end
    pause
    clf
end

%%
G093_data = '~/Analysis/bruce/G093/gamma/grateev_cond_trig_avgs.mat';
load(G093_data);
znorm_cond_spec = bsxfun(@minus,cond_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
znorm_cond_spec = bsxfun(@rdivide,znorm_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));

avg_norm_cond_spec = squeeze(mean(norm_cond_spec));
avg_znorm_cond_spec = squeeze(mean(znorm_cond_spec));

%%
close all
to_print = 0;
figname = ['~/Analysis/bruce/G093/grating_eta'];
n_events = length(event_names);
avg_trig_sig_cond = squeeze(mean(cond_trig_lfp));
sem_trig_sig_cond = squeeze(std(cond_trig_lfp))/sqrt(length(use_lfps));
for cc = 1:n_events
figname = ['~/Analysis/bruce/G093/grating_eta_' sprintf('Cond %d',cc)];
    
    fprintf('Event type %d\n',cc);
    disp(event_names{cc})
    subplot(3,1,1)
    shadedErrorBar(lags/Fsd,avg_trig_sig_cond(cc,:),sem_trig_sig_cond(cc,:));
    xlim([-0.1 0.6]);
    ylim([-12 4]*1e-5);
    xlabel('Time (s)','fontsize',16); ylabel('Amplitude (V)','fontsize',16);
    xl = xlim();
    line(xl,[0 0],'color','k','linestyle','--');
    subplot(3,1,[2 3])
    pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cc,:,:))');shading flat
    caxis([-1 1])
    xlim([-0.1 0.6]);
    xlabel('Time (s)','fontsize',16); ylabel('Frequency (Hz)','fontsize',16);
    title(sprintf('Event type: %s',event_names{cc}));
    if to_print == 1
        fillPage(gcf,'Papersize',[5 6]);
%         if cc == 1
            print(figname,'-dpng');
%         else
%             print(figname,'-dpsc','-append');
%         end
        close;
    else
        pause
        clf
    end
end
