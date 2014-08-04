
cd ~/Analysis/bruce/summary_analysis/eyetrack_figs/

unit_range = 1:80;
unit_range(9) = [];

target_trial = 11;
data_range = find(all_trialvec(used_inds) == target_trial);

close all
figure
space = 0.1;
ms = 10;
lheight = 0.09;
lwidth = 0.5;
for ii = 1:length(unit_range);
    ii
    cur_spk_inds = find(Robs_mat(data_range,unit_range(ii)) > 0);
%     cur_spk_inds = find(Robs(used_inds(data_range),unit_range(ii)) > 0);
    cur_spk_times = cur_spk_inds*dt;
    jitter = rand(size(cur_spk_times))*dt;
%     plot(cur_spk_times+jitter,ones(size(cur_spk_times))*ii*space,'k.','markersize',ms);
    for jj = 1:length(cur_spk_times)
       line([cur_spk_times(jj) + jitter(jj) cur_spk_times(jj) + jitter(jj)],[ii*space ii*space+lheight],'color','k','linewidth',lwidth); 
    end
hold on
end
set(gca,'ytick',[]);
set(gca,'xtick',[]);
ylabel('Unit number','fontsize',12);
xlabel('Time (s)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
ylim([0 80*space]);
xlim([0.1 0.6])
box off
fillPage(gcf,'papersize',[5 5]);
fname = 'examp_raster5';