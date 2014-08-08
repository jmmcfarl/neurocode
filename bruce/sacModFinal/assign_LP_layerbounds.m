clear all
close all
%%
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/layer_boundaries/';

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297','M297'};
ori_list = [80 60 135 70 140 90 160 40 45 0 90];
base_sname = 'lfp_trig_avgs';

for ee = 10:length(Expt_list)
    Expt_name = Expt_list{ee};
    bar_ori = ori_list(ee);
    save_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod'];
    cd(save_dir);
    sname = [base_sname sprintf('_ori%d',bar_ori)];
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
        
        
        powSpec = log10(lfp_data.powSpec);
        zpowSpec = (nanzscore(powSpec'));
        f = lfp_data.f;
        logf = logspace(log10(f(2)),log10(f(end)),length(f)-1);
        interp_zpowSpec = interp1(f,zpowSpec',logf)';
        
        subplot(2,1,2);
%         pcolor(lfp_data.f(2:end),1:24,(zpowSpec(:,2:end))); caxis([-2 2]); shading flat; set(gca,'ydir','reverse');
        imagesc(logf,1:24,interp_zpowSpec); caxis([-2 2]); 
        set(gca,'xscale','log');
        xlabel('Frequency (Hz)');
        ylabel('Probe');
        
        
        ub = input('Upper boundary of L4?');
        subplot(2,1,1);
        xl = xlim();
        line(xl,[ub ub]+0.5,'color','w');
        subplot(2,1,2);
        xl = xlim();
        line(xl,[ub ub]+0.5,'color','w');
        
        lb = input('Lower boundary of L4?');
        subplot(2,1,1);
        xl = xlim();
        line(xl,[lb lb]+0.5,'color','w');
        subplot(2,1,2);
        xl = xlim();
        line(xl,[lb lb]+0.5,'color','w');
        
        
        keep = input('Keep boundaries? (0/1)');
        
        subplot(2,1,1)
        title(sprintf('CSD lb%d ub%d',lb,ub));
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