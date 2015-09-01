clear all
close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
sess_data = sess_data(thom_pfc);


for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    load used_data lf8 wcv_minus_spike lf4
    [desynch_times_mp{d},desynch_inds_mp{d},P_mp,f,t] = locate_desynch_times_individual(wcv_minus_spike);
    [desynch_times_lf8{d},desynch_inds_lf8{d},P_lf8,f,t] = locate_desynch_times_individual(lf8);
    if sess_data(d).thom_elec
        [desynch_times_lf4{d},desynch_inds_lf4{d},P_lf4,f,t] = locate_desynch_times_individual(lf4);
    end
    
    %compute test statistics from the spectrogram (as used in the
    %locate_desynch_times_individual routine)
    df = f(2)-f(1);
    f_so = find(f > 0.2 & f < 1.5);
    f_hf = find(f > 4);
    max_so_mp = max(10*log10(P_mp(:,f_so)),[],2);
    max_so_lf8 = max(10*log10(P_lf8(:,f_so)),[],2);
    net_hf_mp = zscore(trapz(10*log10(P_mp(:,f_hf)),2)*df);
    net_hf_lf8 = zscore(trapz(10*log10(P_lf8(:,f_hf)),2)*df);
    if sess_data(d).thom_elec
        max_so_lf4 = max(10*log10(P_lf4(:,f_so)),[],2);
        net_hf_lf4 = zscore(trapz(10*log10(P_lf4(:,f_hf)),2)*df);
    end
    
    %%
    Fig = figure(1);
    subplot(3,1,1)
    pcolor(t,f,10*log10(P_lf8'));shading flat;
    caxis([-25 1]);
    set(gca,'yscale','log')
    subplot(3,1,2)
    plot(t,max_so_lf8), hold on
    xlim([t(1) t(end)])
    ylim([-7 3])
    line([t(1) t(end)],[-6 -6],'Color','k')
    if size(desynch_times_lf8{d},1) > 0
        plot(desynch_times_lf8{d}(:,1),ones(size(desynch_times_lf8{d}(:,1)))*-5,'go')
        plot(desynch_times_lf8{d}(:,2),ones(size(desynch_times_lf8{d}(:,2)))*-5,'ro')
    end
    subplot(3,1,3)
    plot(t,net_hf_lf8)
    xlim([t(1) t(end)])
    line([t(1) t(end)],[-2 -2],'Color','k')
    tname = ['F:\WC_Germany\parietal_cortical_2010\desynch_detect\ind_lf8_v2_' sess_data(d).name];
    print('-dpng',tname);
    close
    
    Fig = figure(1);
    subplot(3,1,1)
    pcolor(t,f,10*log10(P_mp'));shading flat;
    caxis([-25 1]);
    set(gca,'yscale','log')
    subplot(3,1,2)
    plot(t,max_so_mp), hold on
    xlim([t(1) t(end)])
    ylim([-7 3])
    line([t(1) t(end)],[-6 -6],'Color','k')
    if size(desynch_times_mp{d},1) > 0
        plot(desynch_times_mp{d}(:,1),ones(size(desynch_times_mp{d}(:,1)))*-5,'go')
        plot(desynch_times_mp{d}(:,2),ones(size(desynch_times_mp{d}(:,2)))*-5,'ro')
    end
    subplot(3,1,3)
    plot(t,net_hf_mp)
    xlim([t(1) t(end)])
    line([t(1) t(end)],[-2 -2],'Color','k')
    tname = ['F:\WC_Germany\parietal_cortical_2010\desynch_detect\ind_mp_v2_' sess_data(d).name];
    print('-dpng',tname);
    close
    
    if sess_data(d).thom_elec
        Fig = figure(1);
        subplot(3,1,1)
        pcolor(t,f,10*log10(P_lf4'));shading flat;
        caxis([-25 1]);
        set(gca,'yscale','log')
        subplot(3,1,2)
        plot(t,max_so_lf4), hold on
        xlim([t(1) t(end)])
        ylim([-7 3])
        line([t(1) t(end)],[-6 -6],'Color','k')
        if size(desynch_times_lf4{d},1) > 0
            plot(desynch_times_lf4{d}(:,1),ones(size(desynch_times_lf4{d}(:,1)))*-5,'go')
            plot(desynch_times_lf4{d}(:,2),ones(size(desynch_times_lf4{d}(:,2)))*-5,'ro')
        end
        subplot(3,1,3)
        plot(t,net_hf_lf4)
        xlim([t(1) t(end)])
        line([t(1) t(end)],[-2 -2],'Color','k')
        tname = ['F:\WC_Germany\parietal_cortical_2010\desynch_detect\ind_lf4_v2_' sess_data(d).name];
        print('-dpng',tname);
        close
    end
    
end

cd F:\WC_Germany\parietal_cortical_2010\
save desynch_times_individual_v2 desynch_times*
