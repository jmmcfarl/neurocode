clear all
close all

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
addpath('G:\WC_Germany\parietal_cortical_2010\')

for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
        
    load used_data lf8

    [desynch_times_lf8,desynch_inds_lf8,P_lf8,f,t] = locate_desynch_times_individual(lf8);
    
    df = f(2)-f(1);
    f_so = find(f > 0.2 & f < 1.5);
    f_hf = find(f > 4);
    max_so_lf8 = max(10*log10(P_lf8(:,f_so)),[],2);
    net_hf_lf8 = zscore(trapz(10*log10(P_lf8(:,f_hf)),2)*df);
        
    %%
    Fig = figure('visible','off');
    subplot(3,1,1)
    pcolor(t,f,10*log10(P_lf8'));shading flat;
    caxis([-25 1]);
    set(gca,'yscale','log')
    subplot(3,1,2)
    plot(t,max_so_lf8), hold on
    xlim([t(1) t(end)])
    ylim([-7 3])
    line([t(1) t(end)],[-6 -6],'Color','k')
    if size(desynch_times_lf8,1) > 0
        plot(desynch_times_lf8(:,1),ones(size(desynch_times_lf8(:,1)))*-5,'go')
        plot(desynch_times_lf8(:,2),ones(size(desynch_times_lf8(:,2)))*-5,'ro')
    end
    subplot(3,1,3)
    plot(t,net_hf_lf8)
    xlim([t(1) t(end)])
    line([t(1) t(end)],[-1.5 -1.5],'Color','k')
    tname = ['G:\WC_Germany\persistent_9_27_2010\desynch_detect_allcort\lf8_' s_name];
    print('-dpng',tname);
    close
    save desynch_times_lf8_v2 desynch_times_lf8 desynch_inds_lf8
    clear desynch_* *_lf8
    
end


