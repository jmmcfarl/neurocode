%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables

d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data wcv_minus_spike lf8 CSC8_SampleFrequencies synct

    Fs = mean(CSC8_SampleFrequencies);
    
%get lf8 up transition points
[trans_points] = getUpTransitionPoints(Fs,lf8);
    
    
    %filter data
    nyqf = Fs/2;
    hif = 40/nyqf;
    [b,a] = butter(2,hif,'low');
    lofwcv = filtfilt(b,a,wcv_minus_spike);

    %zscore data
    lofwcv = lofwcv - mean(lofwcv);
    lofwcv = lofwcv/std(lofwcv);


maxLag = 5;
lagVec = [-round(maxLag*Fs):round(maxLag*Fs)]';
firstUp = find((synct(trans_points)-synct(maxLag))>maxLag,1,'first');
lastUp = find((synct(end) - synct(trans_points))>maxLag,1,'last');
goodUps = length(firstUp:lastUp);
upTrigTot_wcv = zeros(size(lagVec));

  for i = firstUp:lastUp

   upTrigTot_wcv = upTrigTot_wcv + lofwcv(trans_points(i)-round(Fs*maxLag):trans_points(i)+round(Fs*maxLag));
  
  end
  
  upTrigAvg_wcv{d} = upTrigTot_wcv/goodUps;

  figure
  set(gcf,'PaperUnits','centimeters');
  set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
  set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
  plot(lagVec/Fs,upTrigAvg_wcv{d},'k','linewidth',2)
  xlabel('Time (s)','FontSize',14)
  ylabel('Avg WCV ZScore','FontSize',14)
  title('LF8 Up Trig Avg Wcv','FontSize',16)
  tname = ['E:\WC_Germany\JMM_Analysis\up_trig_avg\' f_names{d}];
print('-dpng',tname)
close all
  
  figure
  set(gcf,'PaperUnits','centimeters');
  set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
  set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
  plot(lagVec/Fs,upTrigAvg_wcv{d},'k','linewidth',2)
  xlabel('Time (s)','FontSize',14)
  ylabel('Avg WCV ZScore','FontSize',14)
  xlim([-0.5 1])
  grid on
  title('LF8 Up Trig Avg Wcv','FontSize',16)
  tname = ['E:\WC_Germany\JMM_Analysis\up_trig_avg\zoom_' f_names{d}];
print('-dpng',tname)
close all



    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_upTrigWcv d upTrig*

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_upTrigWcv.mat


end
