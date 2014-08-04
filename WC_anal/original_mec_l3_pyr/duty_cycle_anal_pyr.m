%% calculate excess up state

clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_Analysis_pyr\run_hist_thresh_v2\up_per_data_10_28
load('C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28.mat')

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;

Fsd = 2016/dsf;

segDur = 20;
winSlide = 1;

maxlag = 300; %maximum lag in s

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    
 
    load used_data lf8 wcv_minus_spike
    
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    
    lf8_f = downsample(lf8_f,dsf);
    wcv_f = downsample(wcv_f,dsf);
    
    lf8_f = zscore(lf8_f);
    wcv_f = zscore(wcv_f);
    
    %calculate average up state duration in segmented data
    nsamples = length(wcv_f);
    numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
    duty_w{d} = zeros(1,numWins);
    duty_8{d} = zeros(1,numWins);
    excess_duty{d} = zeros(1,numWins);
    
    thresh_w(d) = (max(wcv_f)+min(wcv_f))/2;
    thresh_8(d) = (max(lf8_f)+min(lf8_f))/2;
    disp(sprintf('WCV threshold: %0.3g',thresh_w(d)))
    disp(sprintf('LF8 threshold: %0.3g',thresh_8(d)))
    
    
for w = 1:numWins
    
    begSamp = round(1+(w-1)*winSlide*Fsd);
    endSamp = begSamp + round(segDur*Fsd);
   
    wcv_seg = wcv_f(begSamp:endSamp);
    lf8_seg = lf8_f(begSamp:endSamp);
    
    duty_w{d}(w) = length(find(wcv_seg > thresh_w(d)))/length(wcv_seg);
    duty_8{d}(w) = length(find(lf8_seg > thresh_8(d)))/length(lf8_seg);
    excess_duty{d}(w) = duty_w{d}(w) - duty_8{d}(w);
    
end


[duty_w_acorr(d,:),lags] = xcov(duty_w{d},maxlag,'coeff');
[duty_8_acorr(d,:),lags] = xcov(duty_8{d},maxlag,'coeff');
[excess_acorr(d,:),lags] = xcov(excess_duty{d},maxlag,'coeff');
[avg_acorr(d,:),lags] = xcov(avg_up_dur{d},maxlag,'coeff');
[avg_acorr8(d,:),lags] = xcov(avg_up_dur8{d},maxlag,'coeff');

    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    plot(duty_w{d})
    hold on
    plot(duty_8{d},'r')
    plot(excess_duty{d},'k')
    legend('WCV','LF8','Excess')
    subplot(2,1,2)
    plot(lags,duty_w_acorr(d,:))
    hold on
    plot(lags,duty_8_acorr(d,:),'r')
    plot(lags,excess_acorr(d,:),'k')
    legend('WCV','LF8','Excess')
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\duty_cycle\duty_cycle_' f_names{d}];
    print('-dpng',t_names)
    close all
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    plot(avg_up_dur{d})
    hold on
    plot(avg_up_dur8{d},'r')
    legend('WCV','LFP')
    subplot(2,1,2)
    plot(lags,avg_acorr(d,:))
    hold on
    plot(lags,avg_acorr8(d,:),'r')
    legend('WCV','LFP')
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\duty_cycle\avg_updur_' f_names{d}];
    print('-dpng',t_names)
    close all
    
    
end