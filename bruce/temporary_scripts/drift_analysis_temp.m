
ms_eye_vals = corrected_eye_vals_interp;
for ii = 1:n_trials
    temp = find(all_trialvec == ii);
    ms_eye_vals(temp,:) = bsxfun(@minus,ms_eye_vals(temp,:),nanmedian(ms_eye_vals(temp,:)));
end

%%
% fig_dir = '/home/james/Desktop/example_figs_G086/';
fig_dir = '/home/james/Desktop/example_figs_M296/';
close all
min_print_trialdur = 3;
f1 = figure;
for ii = 1:n_trials
    hold on;
    cur_range = find(all_trialvec(used_inds) == ii);
    if length(cur_range)*dt >= min_print_trialdur
        ii
        rt = all_t_axis(used_inds(cur_range));
        sp = rt(1);
        rt = rt - sp;
        shadedErrorBar(rt,fin_tot_corr(cur_range)*drift_stim_dx,fin_tot_std(cur_range)*drift_stim_dx,{'color','r'});
%         shadedErrorBar(rt,vfin_tot_corr(cur_range)*drift_stim_dx,vfin_tot_std(cur_range)*drift_stim_dx,{'color','r'});
%         shadedErrorBar(rt,afin_tot_corr(cur_range)*drift_stim_dx,afin_tot_std(cur_range)*drift_stim_dx,{'color','k'});
        % plot(1:uNT,corrected_eye_vals_interp(used_inds(1:uNT),[2]),'g','linewidth',2)
        % plot(1:uNT,corrected_eye_vals_interp(used_inds(1:uNT),[4]),'b','linewidth',2)
        plot(rt,ms_eye_vals(used_inds(cur_range),[2]),'b','linewidth',1.5)
        plot(rt,ms_eye_vals(used_inds(cur_range),[4]),'m','linewidth',1.5)
        
        
        xlabel('Time (s)');
        ylabel('Eye position (deg)');
        ylim([-0.4 0.4]);
        xlim(rt([1 end]));
        
        yl = ylim();
        cur_fix_starts = find(ismember(fix_start_inds,cur_range));
        for jj = 1:length(cur_fix_starts)
%         for jj = 2
            line(all_t_axis(used_inds(fix_start_inds(cur_fix_starts(jj)))) - [sp sp],yl,'color','k','linestyle','--','linewidth',0.5);
            line(all_t_axis(used_inds(fix_stop_inds(cur_fix_starts(jj)))) - [sp sp],yl,'color','k','linestyle','--','linewidth',0.5);
%               X = [all_t_axis(used_inds(fix_stop_inds(cur_fix_starts(jj-1)))) all_t_axis(used_inds(fix_start_inds(cur_fix_starts(jj))))];
%               X = [X fliplr(X)] - sp;
%               Y = [yl(1) yl(1) yl(2) yl(2)];
%               p = patch(X,Y,'k');
%               set(p,'FaceAlpha',0.25,'linestyle','none');
        end
        line(xlim(),[0 0],'color','k','linestyle','--');
        pause
        clf         

%         figufy(f1);
%         fig_width = 6; rel_height = 0.6;
%         fname = [fig_dir sprintf('examp_trial%d.pdf',ii)];
%         exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%         clf
    end
end
%%
min_fix_dur = round(0.3/dt);
bbuff_win = round(0.05/dt);
ebuff_win = round(0.00/dt);
inf_drift = nan(NT,1);
meas_drift = nan(NT,2);
for ii = 1:n_fixs
    temp = find(fix_ids == ii);
    if length(temp) >= min_fix_dur
        temp(1:bbuff_win) = [];
        temp((end-ebuff_win+1):end) = [];
        inf_drift(temp) = fin_tot_corr(temp) - mean(fin_tot_corr(temp));
        meas_drift(temp,1) = ms_eye_vals(used_inds(temp),2) - mean(ms_eye_vals(used_inds(temp),2));
        meas_drift(temp,2) = ms_eye_vals(used_inds(temp),4) - mean(ms_eye_vals(used_inds(temp),4));
    end
end
inf_drift = inf_drift*drift_stim_dx;
%%
outlier_thresh = 10; %in z-score
ctype = 'spearman';

meas_drift_z = bsxfun(@minus,meas_drift,nanmedian(meas_drift));
meas_drift_z = bsxfun(@rdivide,meas_drift,robust_std_dev(meas_drift_z));

fprintf('\n\n');
uset1 = find(abs(meas_drift_z(:,1)) < outlier_thresh);
corr1 = corr(inf_drift(uset1),meas_drift(uset1,1),'type',ctype)
% fitp1 = regress(inf_drift(uset1),[meas_drift(uset1,1) ones(length(uset1),1)]);
% slope1 = fitp1(1)
fitp1 = robustfit(meas_drift(uset1,1),inf_drift(uset1));
slope1 = fitp1(2)
fitp1G = GMregress(meas_drift(uset1,1),inf_drift(uset1))
uset2 = find(abs(meas_drift_z(:,2)) < outlier_thresh);
corr2 = corr(inf_drift(uset2),meas_drift(uset2,2),'type',ctype)
% fitp2 = regress(inf_drift(uset2),[meas_drift(uset2,2) ones(length(uset2),1)]);
% slope2 = fitp2(1)
fitp2 = robustfit(meas_drift(uset2,2),inf_drift(uset2));
fitp2G = GMregress(meas_drift(uset2,2),inf_drift(uset2))
slope2 = fitp2(2)

uset12 = find(abs(meas_drift_z(:,1)) < outlier_thresh & abs(meas_drift_z(:,2)) < outlier_thresh);
corr12 = corr(meas_drift(uset12,1),meas_drift(uset12,2),'type',ctype)
% fitp12 = regress(meas_drift(uset12,1),[meas_drift(uset12,2) ones(length(uset12),1)]);
% slope12 = fitp12(1)
% fitp12 = robustfit(meas_drift(uset12,1),meas_drift(uset12,2));
fitp12G = GMregress(meas_drift(uset12,1),meas_drift(uset12,2));
slope12 = fitp12(2)

corravg = corr(mean(meas_drift(uset12,:),2),inf_drift(uset12),'type',ctype)
fitpavgG = GMregress(mean(meas_drift(uset12,:),2),inf_drift(uset12));

%%
fig_dir = '/home/james/Desktop/example_figs_M296/';
% close all
f1 = figure();
subplot(2,1,1)
h = DensityPlot_jmm(meas_drift(uset1,1),inf_drift(uset1),'nbins',[100 100],'logsc',...
    'sd',[0.5 0.5],'xrange',[-0.1 0.1],'yrange',[-0.1 0.1]);
set(gca,'ydir','normal');
hold on
xlabel('Measured drift position (deg)');
ylabel('Inferred drift position (deg)');
line([-0.1 0.1],[-0.1 0.1],'color','k','linewidth',2)
xpred = linspace(-0.1,0.1,100);
% ypred = xpred*fitp(1) + fitp(2);
% ypred = xpred*fitp1(2) + fitp1(1);
ypred = xpred*fitp1G;
plot(xpred,ypred,'r','linewidth',2);
title('Comparison to coil 1')

subplot(2,1,2)
h = DensityPlot_jmm(meas_drift(uset2,2),inf_drift(uset2),'nbins',[100 100],'logsc',...
    'sd',[0.5 0.5],'xrange',[-0.1 0.1],'yrange',[-0.1 0.1]);
set(gca,'ydir','normal');
hold on
xlabel('Measured drift position (deg)');
ylabel('Inferred drift position (deg)');
line([-0.1 0.1],[-0.1 0.1],'color','k','linewidth',2)
xpred = linspace(-0.1,0.1,100);
% ypred = xpred*fitp(1) + fitp(2);
% ypred = xpred*fitp2(2) + fitp2(1);
ypred = xpred*fitp2G ;
plot(xpred,ypred,'r','linewidth',2);
title('Comparison to coil 2')

figufy(f1);
fig_width = 7; rel_height = 2;
fname = [fig_dir 'inferred_coil_compare.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
uset = find(abs(meas_drift_z(:,1)) < outlier_thresh & abs(meas_drift_z(:,2)) < outlier_thresh);
fig_dir = '/home/james/Desktop/example_figs_M296/';
% close all
f1 = figure();
h = DensityPlot_jmm(meas_drift(uset,1),meas_drift(uset,2),'nbins',[100 100],'logsc',...
    'sd',[0.5 0.5],'xrange',[-0.1 0.1],'yrange',[-0.1 0.1]);
set(gca,'ydir','normal');
hold on
xlabel('Measured drift position (deg)');
ylabel('Inferred drift position (deg)');
line([-0.1 0.1],[-0.1 0.1],'color','k','linewidth',2)
xpred = linspace(-0.1,0.1,100);
% ypred = xpred*fitp12(1) + fitp12(2);
% ypred = xpred*fitp12(2) + fitp12(1);
ypred = xpred*fitp12G;
plot(xpred,ypred,'r','linewidth',2);
title('Comparison coils1-2')

figufy(f1);
fig_width = 4; rel_height = 1;
fname = [fig_dir 'coil_coil_compare.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
