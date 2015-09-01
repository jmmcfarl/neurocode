clear all
close all

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = find_struct_field_vals(sess_data,'layer','5');
deep = [deep find_struct_field_vals(sess_data,'layer','56')];
used_cells = [parietal frontal prefrontal];
n = length(used_cells);

raw_Fs = 2016;
dsf = 8;
lcf = 0.05;
hcf = 100;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
[b,a] = butter(2,[lcf/niqf hcf/niqf]);

params.Fs = Fsd;
params.fpass = [0 80];
params.tapers = [5 9];
params.err = 0;
win = 100;

freqs = linspace(0,80,500);

spectrum_mp = nan(n,length(freqs));
spectrum_lf8 = nan(n,length(freqs));
spectrum_lf7 = nan(n,length(freqs));
spectrum_lf6 = nan(n,length(freqs));
spectrum_lf5 = nan(n,length(freqs));
spectrum_lf4 = nan(n,length(freqs));
spectrum_lf3 = nan(n,length(freqs));
spectrum_lf2 = nan(n,length(freqs));
n8 = 0;
n7 = 0;
n6 = 0;
n5 = 0;
n4 = 0;
n3 = 0;
n2 = 0;

for d = 39:n
    cd(sess_data(used_cells(d)).directory)
    d
    
    load used_data wcv_minus_spike lf8 lf7 lf6 lf5 lf4 lf3 lf2
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
    lf7_f = filtfilt(b,a,lf7);
    lf7_f = downsample(lf7_f,dsf)/sess_data(d).gains(7);
    lf6_f = filtfilt(b,a,lf6);
    lf6_f = downsample(lf6_f,dsf)/sess_data(d).gains(6);
    lf5_f = filtfilt(b,a,lf5);
    lf5_f = downsample(lf5_f,dsf)/sess_data(d).gains(5);
    if exist('lf4','var')
        lf4_f = filtfilt(b,a,lf4);
        lf4_f = downsample(lf4_f,dsf)/sess_data(d).gains(4);
    end
    lf3_f = filtfilt(b,a,lf3);
    lf3_f = downsample(lf3_f,dsf)/sess_data(d).gains(3);
    lf2_f = filtfilt(b,a,lf2);
    lf2_f = downsample(lf2_f,dsf)/sess_data(d).gains(2);
    
    clear wcv_minus_spike lf8 lf7 lf6 lf5 lf4 lf3 lf2
    
    if sess_data(d).gains(1) > 0 & exist('wcv_f','var')
        [S,f]=mtspectrumsegc(wcv_f,win,params);
        spectrum_mp(d,:) = interp1(f,S,freqs);
    end
    if sess_data(d).gains(2) > 0 & exist('wcv_f','var')
        [S,f]=mtspectrumsegc(lf2_f,win,params);
        spectrum_lf2(d,:) = interp1(f,S,freqs);
        n2 = n2+1;
    end
    if sess_data(d).gains(3) > 0 & exist('lf3_f','var')
        [S,f]=mtspectrumsegc(lf3_f,win,params);
        spectrum_lf3(d,:) = interp1(f,S,freqs);
        n3 = n3+1;
    end
    if sess_data(d).gains(4) > 0 & exist('lf4_f','var')
        [S,f]=mtspectrumsegc(lf4_f,win,params);
        spectrum_lf4(d,:) = interp1(f,S,freqs);
        n4 = n4+1;
    end
    if sess_data(d).gains(5) > 0 & exist('lf5_f','var')
        [S,f]=mtspectrumsegc(lf5_f,win,params);
        spectrum_lf5(d,:) = interp1(f,S,freqs);
        n5 = n5+1;
    end
    if sess_data(d).gains(6) > 0 & exist('lf6_f','var')
        [S,f]=mtspectrumsegc(lf6_f,win,params);
        spectrum_lf6(d,:) = interp1(f,S,freqs);
        n6 = n6+1;
    end
    if sess_data(d).gains(7) > 0 & exist('lf7_f','var')
        [S,f]=mtspectrumsegc(lf7_f,win,params);
        spectrum_lf7(d,:) = interp1(f,S,freqs);
        n7 = n7+1;
    end
    if sess_data(d).gains(8) > 0 & exist('lf8_f','var')
        [S,f]=mtspectrumsegc(lf8_f,win,params);
        spectrum_lf8(d,:) = interp1(f,S,freqs);
        n8 = n8+1;
    end
end

%%
cd G:\WC_Germany\parietal_cortical_2010

avg_spectrum8 = nanmean(log10(spectrum_lf8));
sem_spectrum8 = nanstd(log10(spectrum_lf8))/sqrt(n8);
avg_spectrum7 = nanmean(log10(spectrum_lf7));
sem_spectrum7 = nanstd(log10(spectrum_lf7))/sqrt(n7);
avg_spectrum6 = nanmean(log10(spectrum_lf6));
sem_spectrum6 = nanstd(log10(spectrum_lf6))/sqrt(n6);
avg_spectrum5 = nanmean(log10(spectrum_lf5));
sem_spectrum5 = nanstd(log10(spectrum_lf5))/sqrt(n5);
avg_spectrum4 = nanmean(log10(spectrum_lf4));
sem_spectrum4 = nanstd(log10(spectrum_lf4))/sqrt(n4);
avg_spectrum3 = nanmean(log10(spectrum_lf3));
sem_spectrum3 = nanstd(log10(spectrum_lf3))/sqrt(n3);
avg_spectrum2 = nanmean(log10(spectrum_lf2));
sem_spectrum2 = nanstd(log10(spectrum_lf2))/sqrt(n2);

cmap = colormap(jet(7));
errorbar(freqs,avg_spectrum2,sem_spectrum2,'color',cmap(1,:)), hold on
errorbar(freqs,avg_spectrum3,sem_spectrum3,'color',cmap(2,:))
errorbar(freqs,avg_spectrum4,sem_spectrum4,'color',cmap(3,:))
errorbar(freqs,avg_spectrum5,sem_spectrum5,'color',cmap(4,:))
errorbar(freqs,avg_spectrum6,sem_spectrum6,'color',cmap(5,:))
errorbar(freqs,avg_spectrum7,sem_spectrum7,'color',cmap(6,:))
errorbar(freqs,avg_spectrum8,sem_spectrum8,'color',cmap(7,:))
set(gca,'xscale','log')
legend('LF2','LF3','LF4','LF5','LF6','LF7','LF8')
