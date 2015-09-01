%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


% %initialize overall variables
d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd
%     load E:\WC_Germany\JMM_Analysis\wcv_Gains
    load used_data CSC8_SampleFrequencies wcv wcv_minus_spike spkid dt synct lf8 
    load used_data_multLFP_2 lf7 lf3
    load used_data_multLFP lf2 lf5
    
    Fs = mean(CSC8_SampleFrequencies);
    
    %find spk amplitudes
%     for i= 1:length(spkid)
%         begWin = spkid-10;
%         endWin = spkid + 100;
%         [spkAmp{d}(i),peakLoc] = max(hiwcv(begWin:endWin));
%         leftEdge = find(hiwcv(begWin:endWin)>spkAmp{d}(i)/2,1,'first');
%         rightEdge = find(hiwcv(begWin:endWin)>spkAmp{d}(i)/2,1,'last');
%         spkWidth{d}(i) = (rightEdge - leftEdge)/Fs;
%     end
    

%filter    
    nyqf = Fs/2;
    lif = 10/nyqf;
    [b,a] = butter(2,lif,'high');
%     hiwcv = filtfilt(b,a,wcv_minus_spike);
    lf2 = filtfilt(b,a,lf2);
    lf3 = filtfilt(b,a,lf3);
    lf5 = filtfilt(b,a,lf5);
    lf7 = filtfilt(b,a,lf7);
    lf8 = filtfilt(b,a,lf8);


    %zscore    
%     wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
%     wcv_z = wcv_z/std(wcv_z);
    lf8 = lf8-mean(lf8);
    lf8 = lf8/std(lf8);
     lf3 = lf3-mean(lf3);
    lf3 = lf3/std(lf3);

    lf2 = lf2 - mean(lf2);
    lf2 = lf2/std(lf2);
    lf5 = lf5 - mean(lf5);
    lf5 = lf5/std(lf5);
    lf7 = lf7 - mean(lf7);
    lf7 = lf7/std(lf7);
    
% %     downsample signals 
% %     dsf = 3;
% %     Fs_d = round(Fs/dsf);
% %     down_wcv = downsample(wcv_z,dsf);
% %     down_lf8 = downsample(lf8,dsf);
% %     down_lf3 = downsample(lf3,dsf);
% %     maxLag = round(Fs_d*1);
% %     
% %     find spike triggered wcv and lf8
% %     spkVec = zeros(size(down_wcv));
% %     spkVec(round(spkid/dsf))=1;


maxLag = 1;
lagVec = [-round(maxLag*Fs):round(maxLag*Fs)]';
firstSpike = find((synct(spkid)-synct(1))>maxLag,1,'first');
lastSpike = find((synct(end) - synct(spkid))>maxLag,1,'last');
goodSpikes = length(firstSpike:lastSpike);
% spkTrigTot_wcv = zeros(size(lagVec));
spkTrigTot_lf8 = zeros(size(lagVec));
spkTrigTot_lf3 = zeros(size(lagVec));
spkTrigTot_lf2 = zeros(size(lagVec));
spkTrigTot_lf7 = zeros(size(lagVec));
spkTrigTot_lf5 = zeros(size(lagVec));

for i = firstSpike:lastSpike

%    spkTrigTot_wcv = spkTrigTot_wcv + wcv_z(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
   spkTrigTot_lf8 = spkTrigTot_lf8 + lf8(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
   spkTrigTot_lf3 = spkTrigTot_lf3 + lf3(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
    spkTrigTot_lf2 = spkTrigTot_lf2 + lf2(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
    spkTrigTot_lf5 = spkTrigTot_lf5 + lf5(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
    spkTrigTot_lf7 = spkTrigTot_lf7 + lf7(spkid(i)-round(Fs*maxLag):spkid(i)+round(Fs*maxLag));
end

% spkTrigAvg_wcv{d} = spkTrigTot_wcv/goodSpikes;
spkTrigAvg_lf8{d} = spkTrigTot_lf8/goodSpikes;
spkTrigAvg_lf3{d} = spkTrigTot_lf3/goodSpikes;

spkTrigAvg_lf2{d} = spkTrigTot_lf2/goodSpikes;
spkTrigAvg_lf5{d} = spkTrigTot_lf5/goodSpikes;
spkTrigAvg_lf7{d} = spkTrigTot_lf7/goodSpikes;

figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lagVec/Fs,spkTrigAvg_wcv{d},'k','linewidth',2)
hold on
plot(lagVec/Fs,spkTrigAvg_lf2{d},'r','linewidth',2)
plot(lagVec/Fs,spkTrigAvg_lf3{d},'g','linewidth',2)
plot(lagVec/Fs,spkTrigAvg_lf5{d},'c','linewidth',2)
plot(lagVec/Fs,spkTrigAvg_lf8{d},'b','linewidth',2)
legend('LF2','LF3','LF5','LF8')
xlabel('Time (s)','FontSize',14);
    ylabel('Spike Triggered Avg','FontSize',14)
    xlim([-0.1 0.2])
    grid on
    title(sprintf('HP STA %d spikes',length(spkid)),'FontSize',16)
    tname = ['E:\WC_Germany\JMM_Analysis\spk_data\hp_spk_trig_avg_' f_names{d}];
    print('-dpng',tname)
    close all
% 
    
%  spkrate(d) = length(spkid)/max(synct);
    
    d = d +1

%     save E:\WC_Germany\JMM_Analysis\overall_data_jmm_spkdat.mat spkTrigAvg* ...
%         lagVec spkAmp spkWidth d

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_spkdat_hp.mat spkTrigAvg* d lagVec

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_spkdat_hp.mat 


end



% 
% 
% 
