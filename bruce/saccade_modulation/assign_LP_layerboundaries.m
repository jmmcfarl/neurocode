clear all
close all

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/layer_boundaries/';

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
sname = 'lfp_trig_avgs';
use_w = 5;
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    save_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir);
    load(sname);
    keep = 0;
    
    fig_name = [fig_dir sprintf('%s_layBounds.pdf',Expt_name)];
    
    h = figure();
    while ~keep
        subplot(2,1,1);
        imagesc(lags/Fsd,1:24,lfp_data.trial_onset_csd)
        caxis([-0.3 0.3]*1e7);
        xlim([0 0.2]);
        title('CSD');
        xlabel('Time since trial onset (s)');
        ylabel('Probe');
        
        subplot(2,1,2);
%         pcolor(lags/Fsd,1:24,squeeze(lfp_data.onset_specgram(:,use_w,:))');shading flat
        imagesc(lags/Fsd,1:24,squeeze(lfp_data.onset_specgram(:,use_w,:))');
        xlim([0 0.2]);
%         set(gca,'ydir','reverse');
        title(sprintf('Specgram %.1f Hz',wfreqs(use_w)));
        xlabel('Time since trial onset (s)');
        ylabel('Probe');
        
        
        ub = input('Upper boundary of L4?');
        subplot(2,1,1);
        xl = xlim();
        line(xl,[ub ub]+0.5,'color','w');
        subplot(2,1,2);
        line(xl,[ub ub]+0.5,'color','w');
        
        lb = input('Lower boundary of L4?');
        subplot(2,1,1);
        xl = xlim();
        line(xl,[lb lb]+0.5,'color','w');
        subplot(2,1,2);
        line(xl,[lb lb]+0.5,'color','w');
        
        
        keep = input('Keep boundaries? (0/1)');
        
    end
    
    boundary_class(ee).lb = lb;
    boundary_class(ee).ub = ub;
    boundary_class(ee).Expt_num = str2num(Expt_name(2:end));
    
    fig_width = 5;
    rel_height = 1.6;
    figufy(h);
    exportfig(h,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(h);
end

%%
sname = [fig_dir 'layer_classification'];
save(sname,'boundary_class');