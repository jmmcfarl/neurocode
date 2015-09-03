clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

window = 30; %window size in sec
win_slide = 2; %window slide

dsf = 8;
Fsd = 2016/dsf;
niqf = Fsd/2;
hif = 2/niqf;
[b,a] = butter(2,hif,'low');

for d = 26:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,dsf);
    wcv_down = downsample(wcv_minus_spike,dsf);

    lf8_down = filtfilt(b,a,lf8_down);
    wcv_down = filtfilt(b,a,wcv_down);
    
    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);

record_dur = length(lf8_down_z)/Fsd;
numWins = floor((record_dur-window)/win_slide);

for i = 1:numWins
    
   begPt = (i-1)*win_slide*Fsd+1;
   endPt = begPt+window*Fsd;
   
   corTime(i) = (begPt+endPt)/2/Fsd;
   
   lf8Seg = lf8_down_z(begPt:endPt);
   wcvSeg = wcv_down_z(begPt:endPt);
   
   [u8(i,:),sig8(i,:),t8(i,:),iter] = fit_mix_gaussian(lf8Seg,2);
   [uw(i,:),sigw(i,:),tw(i,:),iter] = fit_mix_gaussian(wcvSeg,2);

   [am,upState] = max(u8(i,:));
   if upState == 2
       downState = 1;
   else
       downState = 2;
       disp('DownState was higher')
   end
   
   dci8(i) = (t8(i,upState)-t8(i,downState))/(t8(i,1)+t8(i,2));
   dciw(i) = (tw(i,upState)-tw(i,downState))/(tw(i,1)+tw(i,2));
 
   clf
   subplot(3,1,1)
   plot(lf8Seg)
   hold on
   plot(wcvSeg,'k')
   subplot(3,1,2)
   plot_mix_gaussian(uw(i,:),sigw(i,:),tw(i,:),wcvSeg)
   subplot(3,1,3)
   plot_mix_gaussian(u8(i,:),sig8(i,:),t8(i,:),lf8Seg)
   
   
   
end
      
gaus_prob8{d} = t8;
gaus_probw{d} = tw;

dci_lf8{d} = dci8;
dci_wcv{d} = dciw;


    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
plot(corTime,dci8,corTime,dciw,'r','linewidth',2)
    legend('Lf8','Wcv')
    xlabel('Time (s)','FontSize',14)
    ylabel('DCI','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\gmm_dci\' f_names{d}];
    print('-dpng',tname);
    
 clear u8 uw sig8 sigw t8 tw numWins lfpSeg wcvSeg corTime lf8_down* wcv_down* dci8 dciw   
    
save E:/WC_Germany/jmM_Analysis/gmm_dci dci_lf8 dci_wcv gaus_prob*

 
end

