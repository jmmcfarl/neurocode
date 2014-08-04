clear all
close all

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/layer_boundaries/';
sname = [fig_dir 'layer_classification'];
load(sname);

xl = [-0.05 0.2];
ca = [0.6 1.4];
% ca = [1 1.4];

    h1 = figure;
    h2 = figure;
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    sname = 'sac_trig_avg_data3';
    sdir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(sdir)
    load(sname);
    
    sname = 'lfp_trig_avgs';
    lfp_dat = load(sname);
    
    figure(h1);
    subplot(2,2,1)
    imagesc(lags*dt,1:24,mua_data.gsac_gray_avg')
    caxis(ca);
    xlim(xl);
    line(xl,boundary_class(ee).lb + [0 0],'color','w');
    line(xl,boundary_class(ee).ub + [0 0],'color','w');
    
    subplot(2,2,2)
    imagesc(lags*dt,1:24,mua_data.msac_gray_avg')
    caxis(ca);
    xlim(xl);
    line(xl,boundary_class(ee).lb + [0 0],'color','w');
    line(xl,boundary_class(ee).ub + [0 0],'color','w');
    
    subplot(2,2,3)
    imagesc(lfp_dat.lags/lfp_dat.Fsd,1:24,lfp_dat.lfp_data.trial_onset_csd)
    caxis([-0.3 0.3]*1e7);
    xlim(xl);
    title('CSD');
    xlabel('Time since trial onset (s)');
    ylabel('Probe');
    line(xl,boundary_class(ee).lb + [0 0],'color','w');
    line(xl,boundary_class(ee).ub + [0 0],'color','w');

    use_w = 5;
    cur_specgram = lfp_dat.lfp_data.onset_specgram;
%     cur_specgram = lfp_dat.lfp_data.sac_specgram;
%     cur_specgram = bsxfun(@minus,cur_specgram,lfp_dat.lfp_data.ampgrams_mean);
%     cur_specgram = bsxfun(@rdivide,cur_specgram,lfp_dat.lfp_data.ampgrams_std);
    
            subplot(2,2,4);
%         imagesc(lfp_dat.lags/lfp_dat.Fsd,1:24,squeeze(cur_specgram(:,use_w,:))');
    imagesc(lfp_dat.lags/lfp_dat.Fsd,1:24,lfp_dat.lfp_data.gsac_csd)
    xlim(xl);
    caxis([-0.2 0.2]*1e7);
%         title(sprintf('Specgram %.1f Hz',lfp_dat.wfreqs(use_w)));
        xlabel('Time since trial onset (s)');
        ylabel('Probe');
    line(xl,boundary_class(ee).lb + [0 0],'color','w');
    line(xl,boundary_class(ee).ub + [0 0],'color','w');

    
    temp = mua_data.gsac_gray_avg;
    [~,minloc] = min(temp);
    [~,maxloc] = max(temp);
    figure(h2);
    subplot(2,1,1)
    plot(lags(minloc)*dt,'o-')
    yl = ylim();
    line(boundary_class(ee).lb + [0 0],yl,'color','k');
    line(boundary_class(ee).ub + [0 0],yl,'color','k');
     subplot(2,1,2)
    plot(lags(maxloc)*dt,'o-')
    yl = ylim();
    line(boundary_class(ee).lb + [0 0],yl,'color','k');
    line(boundary_class(ee).ub + [0 0],yl,'color','k');
   
    pause
    figure(h1);clf;figure(h2);clf;
end







