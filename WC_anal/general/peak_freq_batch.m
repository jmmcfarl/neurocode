clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

params.fpass = [0 2];
params.err = [0];
params.tapers = [1 1];
params.Fs = 504;
window = [40 10];

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd


    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,4);
    wcv_down = downsample(wcv_minus_spike,4);

    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);


    [S8,t,f] = mtspecgramc(lf8_down_z,window,params);
    [Sw,t,f] = mtspecgramc(wcv_down_z,window,params);

    %set DC term to 0
    S8(:,1) = 0;
    Sw(:,1) = 0;
    
[a8,b8] = max(S8,[],2);
[aw,bw] = max(Sw,[],2);

peakFreq_8{d} = f(b8);
peakFreq_w{d} = f(bw);   
    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     plot(t,peakFreq_8{d},'.-')
%     hold on
%     plot(t,peakFreq_w{d},'r.-')
%     legend('LF8','WCV')
%     tname = ['E:\WC_Germany\JMM_Analysis\peak_freq\' f_names{d}];
%     print('-dpng',tname);
%         
    clear a8 b8 aw bw t f S8 Sw lf8* wcv*


end

peak_pyr = [];
peak_stell = [];
peak_8 = [];

for i = 1:17
    peak_stell = [peak_stell peakFreq_w{i}];
    peak_8 = [peak_8 peakFreq_8{i}];
end

for i = 18:43
    peak_pyr = [peak_pyr peakFreq_w{i}];
    peak_8 = [peak_8 peakFreq_8{i}];
end

freq_grid = [0:.025:1];
np = hist(peak_pyr,freq_grid);
ns = hist(peak_stell,freq_grid);
n8 = hist(peak_8,freq_grid);
np = np/sum(np);
ns = ns/sum(ns);
n8 = n8/sum(n8);

figure
plot(freq_grid,np,'.-')
hold on
plot(freq_grid,ns,'r.-')
plot(freq_grid,n8,'k.-')
legend('Pyr','Stell','Lf8')
