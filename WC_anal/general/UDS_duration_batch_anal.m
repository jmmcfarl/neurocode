clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
g = [0:0.1:20];
for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

[up_state_dur{d},down_state_dur{d},up_trans{d},down_trans{d},sig_f1,sig_f2,threshold{d},effect_size{d}] = UDS_extract_9_5(wcv_minus_spike);
[up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_f1,sig_f2,threshold8{d},effect_size8{d}] = UDS_extract_9_5(lf8);

up_dist{d} = hist(up_state_dur{d},g);
up_dist{d} = up_dist{d}/sum(up_dist{d});
down_dist{d} = hist(down_state_dur{d},g);
down_dist{d} = down_dist{d}/sum(down_dist{d});

up_dist8{d} = hist(up_state_dur8{d},g);
up_dist8{d} = up_dist8{d}/sum(up_dist8{d});
down_dist8{d} = hist(down_state_dur8{d},g);
down_dist8{d} = down_dist8{d}/sum(down_dist8{d});
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    max_up = max([max(up_state_dur{d}) max(up_state_dur8{d})]);
    plot(g,fgsmooth(up_dist{d},3))
    hold on
    plot(g,fgsmooth(up_dist8{d},3),'r')
    xlim([0 max_up+2])
    legend('WCV','LF8')
    title('Up State Duration','FontSize',14)
    subplot(2,1,2)
    plot(g,fgsmooth(down_dist{d},3))
    hold on
    plot(g,fgsmooth(down_dist8{d},3),'r')
    max_down = max([max(down_state_dur{d}) max(down_state_dur8{d})]);
    xlim([0 max_down+2])
    legend('WCV','LF8')
    title('Down State Duration','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\UDS_duration\' f_names{d}];
    print('-dpng',tname)
    close

    save E:\WC_Germany\JMM_Analysis\UDS_dur_data up_state_dur* down_state_dur* up_trans* down_trans* threshold* effect_size*
    
end