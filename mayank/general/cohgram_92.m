clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

params.fpass = [0 2];
% params.err = [2 .05];
params.tapers = [3 5];
dsf = 8;
params.Fs = 2016/dsf;
window = [20 5];

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd


    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,dsf);
    wcv_down = downsample(wcv_minus_spike,dsf);

    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);

    params.err = [0];
        [Sw,t,f] = mtspecgramc(wcv_down_z,window,params);
        
        [dummy,peakfreq] = max(Sw,[],2);
        
params.err = [2 .05];
    [C,phi,S12,S1,S2,t,f,confC,phistd,Cerr]=cohgramc(lf8_down_z,wcv_down_z,window,params);

    
    for i = 1:length(t)
        CoherFun{d}(i) = C(i,peakfreq(i));
        CoherErrFun{d}(:,i) = Cerr(:,i,peakfreq(i));
        phiFun{d}(i) = phi(i,peakfreq(i));
        phiStdFun{d}(i) = phistd(i,peakfreq(i));
    end
    
    Fig = figure(2)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    errorbar(t,CoherFun{d},CoherErrFun{d}(1,:),CoherErrFun{d}(2,:))
    line([0 max(t)],[confC confC],'Color','k')
    subplot(2,1,2)
    errorbar(t,phiFun{d},phiStdFun{d},'r')
    tname = ['E:\WC_Germany\JMM_Analysis\coher_atpeakfreq\' f_names{d}];
    print('-dpng',tname);
%         
emp    clear wcv* lf8* C phi S12 S1 S2 t f confC phistd Cerr Sw t f


end

