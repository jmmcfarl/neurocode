clear all
close all

load F:\WC_Germany\overall_EC\overall_allcells_dir
addpath('F:\WC_Germany\parietal_cortical_2010\')

for d = 110:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
        
    load used_data lf3
    
%     lf2 = lf2/sess_data(d).gains(2);
%     lf5 = lf5/sess_data(d).gains(5);
%     lf2_r = lf2-lf5;
    
%     [desynch_times_lf8,desynch_inds_lf8,P_lf8,f,t] = locate_desynch_times_individual(lf8);
    [desynch_times_lf3,desynch_inds_lf3,P_lf3,f,t] = locate_desynch_times_individual(lf3);
% if max(isnan(lf2_r)) > 0
%     desynch_times_lf2r = nan;
%     desynch_inds_lf2r = nan;
%     P_lf2r = nan;
%     bad_lf2r(d) = 1;
% else
%     [desynch_times_lf2r,desynch_inds_lf2r,P_lf2r,f,t] = locate_desynch_times_individual(lf2_r);
%     bad_lf2r(d) = 0;
% end

    df = f(2)-f(1);
    f_so = find(f > 0.2 & f < 1.5);
    f_hf = find(f > 4);
%     max_so_lf8 = max(10*log10(P_lf8(:,f_so)),[],2);
%     net_hf_lf8 = zscore(trapz(10*log10(P_lf8(:,f_hf)),2)*df);
    max_so_lf3 = max(10*log10(P_lf3(:,f_so)),[],2);
    net_hf_lf3 = zscore(trapz(10*log10(P_lf3(:,f_hf)),2)*df);
% if bad_lf2r(d) == 0
%     max_so_lf2r = max(10*log10(P_lf2r(:,f_so)),[],2);
%     net_hf_lf2r = zscore(trapz(10*log10(P_lf2r(:,f_hf)),2)*df);
% end
    %%
%     Fig = figure(1);
%     subplot(3,1,1)
%     pcolor(t,f,10*log10(P_lf8'));shading flat;
%     caxis([-25 1]);
%     set(gca,'yscale','log')
%     subplot(3,1,2)
%     plot(t,max_so_lf8), hold on
%     xlim([t(1) t(end)])
%     ylim([-7 3])
%     line([t(1) t(end)],[-6 -6],'Color','k')
%     if size(desynch_times_lf8,1) > 0
%         plot(desynch_times_lf8(:,1),ones(size(desynch_times_lf8(:,1)))*-5,'go')
%         plot(desynch_times_lf8(:,2),ones(size(desynch_times_lf8(:,2)))*-5,'ro')
%     end
%     subplot(3,1,3)
%     plot(t,net_hf_lf8)
%     xlim([t(1) t(end)])
%     line([t(1) t(end)],[-2 -2],'Color','k')
%     tname = ['F:\WC_Germany\overall_EC\desynch_detect\lf8_' s_name];
%     print('-dpng',tname);
%     close
%      save desynch_times_lf8 desynch_times*

Fig = figure(1);
    subplot(3,1,1)
    pcolor(t,f,10*log10(P_lf3'));shading flat;
    caxis([-25 1]);
    set(gca,'yscale','log')
    subplot(3,1,2)
    plot(t,max_so_lf3), hold on
    xlim([t(1) t(end)])
    ylim([-7 3])
    line([t(1) t(end)],[-6 -6],'Color','k')
    if size(desynch_times_lf3,1) > 0
        plot(desynch_times_lf3(:,1),ones(size(desynch_times_lf3(:,1)))*-5,'go')
        plot(desynch_times_lf3(:,2),ones(size(desynch_times_lf3(:,2)))*-5,'ro')
    end
    subplot(3,1,3)
    plot(t,net_hf_lf3)
    xlim([t(1) t(end)])
    line([t(1) t(end)],[-2 -2],'Color','k')
    tname = ['F:\WC_Germany\overall_EC\desynch_detect\lf3_' s_name];
    print('-dpng',tname);
    close
    save desynch_times_lf3 desynch_times*

% if bad_lf2r(d) == 0
% Fig = figure(1);
%     subplot(3,1,1)
%     pcolor(t,f,10*log10(P_lf2r'));shading flat;
%     caxis([-25 1]);
%     set(gca,'yscale','log')
%     subplot(3,1,2)
%     plot(t,max_so_lf2r), hold on
%     xlim([t(1) t(end)])
%     ylim([-7 3])
%     line([t(1) t(end)],[-6 -6],'Color','k')
%     if size(desynch_times_lf2r,1) > 0
%         plot(desynch_times_lf2r(:,1),ones(size(desynch_times_lf2r(:,1)))*-5,'go')
%         plot(desynch_times_lf2r(:,2),ones(size(desynch_times_lf2r(:,2)))*-5,'ro')
%     end
%     subplot(3,1,3)
%     plot(t,net_hf_lf2r)
%     xlim([t(1) t(end)])
%     line([t(1) t(end)],[-2 -2],'Color','k')
%     tname = ['F:\WC_Germany\overall_EC\desynch_detect\lf3_' s_name];
%     print('-dpng',tname);
%     close
%     save desynch_times_lf2r desynch_times*
% end
end

% cd G:\WC_Germany\overall_EC\
% save desynch_times_individual desynch_times*
