cmap = colormap(jet(6));

for i = 18:43

    figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(lagVec/Fs,spkTrigWCV(i,:),'Color','k','linewidth',2)
    hold on
    plot(lagVec/Fs,spkTrigLF2(i,:),'Color',cmap(2,:),'linewidth',2)
    plot(lagVec/Fs,spkTrigLF3(i,:),'Color',cmap(3,:),'linewidth',2)
    plot(lagVec/Fs,spkTrigLF5(i,:),'Color',cmap(4,:),'linewidth',2)
    plot(lagVec/Fs,spkTrigLF7(i,:),'Color',cmap(5,:),'linewidth',2)
    plot(lagVec/Fs,spkTrigLF8(i,:),'Color',cmap(6,:),'linewidth',2)
    legend('WCV','LF2','LF3','LF5','LF7','LF8')
    xlim([-.1 .1])
    grid on
    xlabel('Time (s)','FontSize',14);
    ylabel('Spike Triggered Avg','FontSize',14)
    title(sprintf('STA (z) %d spikes %0.2g rate',length(spkAmp{i}),spkrate(i)),'FontSize',16)
    tname = ['E:\WC_Germany\JMM_Analysis\spk_trig\pyramidal\spk_trig_avg_all_LFP_zoom_' f_names{i}];
    print('-dpng',tname)
    close all

    
end