clear all
load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\pers_times_oversmooth
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
% load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

dsf = 2;
params.Fs = 2016/dsf;
params.err = [1 .01];
params.fpass = [60 200];
params.tapers = [1 1];
bad_f = [100 150 200];

% W = 0.04;
% win = 100;
niqf = 2016/2;
lcf = 10/niqf;
[b,a] = butter(2,lcf,'high');

%for 10 to 400 band
% exclude_f = [13:17 39:43 64:68 90:94 115:119 140:144 166:170 191:193];

%for 60 to 200 band
exclude_f = [19:23 44:48 70:71];

movingwin = [0.5 0.5];
for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_w = zscore(down_w);
    down_8 = zscore(down_8);

    up_trans8{d} = up_trans8{d}*4;
    down_trans8{d} = down_trans8{d}*4;
    up_trans{d} = up_trans{d}*4;
    down_trans{d} = down_trans{d}*4;

    up_markers8 = [up_trans8{d}' down_trans8{d}'];
    up_markers = [up_trans{d}' down_trans{d}'];
    down_markers8 = [down_trans8{d}(1:end-1)' up_trans8{d}(2:end)'];
    down_markers = [down_trans{d}(1:end-1)' down_trans{d}(2:end)'];

    Swup{d} = zeros(length(up_trans{d}),71);
    tot_up_pow{d} = zeros(length(up_trans{d}),1);
    for i = 1:length(up_trans{d})
        if up_state_dur{d}(i) > 0.5
            [Swup{d}(i,:),f]=mtspectrumsegc(down_w(up_trans{d}(i):down_trans{d}(i)),movingwin,params);
            Swup{d}(i,exclude_f) = nan;
            tot_up_pow{d}(i) = nansum(10*log10(Swup{d}(i,:)));
        else
            Swup{d}(i,:) = nan;
            tot_up_pow{d}(i) = nan;
        end
    end
    
    Swdown{d} = zeros(length(up_trans{d})-1,71);
    tot_down_pow{d} = zeros(length(up_trans{d})-1,1);
    for i = 1:length(up_trans{d})-1
        if down_state_dur{d}(i) > 0.5
            [Swdown{d}(i,:),f]=mtspectrumsegc(down_w(down_trans{d}(i):up_trans{d}(i+1)),movingwin,params);
            Swdown{d}(i,exclude_f) = nan;
            tot_down_pow{d}(i) = nansum(10*log10(Swdown{d}(i,:)));
        else
            Swdown{d}(i,:) = nan;
            tot_down_pow{d}(i) = nan;
        end
    end

pers_up_pow{d} = tot_up_pow{d}(pers_up_inds{d});
npers_up_pow{d} = tot_up_pow{d}(npers_up_inds{d});
pers_down_pow{d} = tot_down_pow{d}(pers_down_inds{d}(1:end-1));
npers_down_pow{d} = tot_down_pow{d}(npers_down_inds{d}(1:end-1));
    
z_tot_up_pow{d} = tot_up_pow{d}-nanmean(tot_up_pow{d});
z_tot_up_pow{d} = z_tot_up_pow{d}/nanstd(tot_up_pow{d});
z_tot_down_pow{d} = tot_down_pow{d}-nanmean(tot_down_pow{d});
z_tot_down_pow{d} = z_tot_down_pow{d}/nanstd(tot_down_pow{d});

%get rid of outliers 
z_tot_down_pow{d}(abs(z_tot_down_pow{d}) > 5) = [];
z_tot_up_pow{d} = z_tot_up_pow{d}-nanmean(z_tot_up_pow{d});
z_tot_up_pow{d} = z_tot_up_pow{d}/nanstd(tot_up_pow{d});

% ud_rat = z_tot_up_pow{d}(1:end-1)./z_tot_down_pow{d};
% ud_rat_ff(d) = nanstd(ud_rat)/nanmean(ud_rat);


figure
[up_pow_up_dur_cor(d),up_pow_up_dur_p(d)]= scatter_with_cor(z_tot_up_pow{d},up_state_dur{d})
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\up_dur_v_up_pow_60_200' f_names{d}];
    print('-dpng',t_names)
    close all
    
figure
[down_pow_up_dur_cor(d),down_pow_up_dur_p(d)]=scatter_with_cor(z_tot_down_pow{d},up_state_dur{d}(2:end))
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\up_dur_v_down_pow_60_200' f_names{d}];
    print('-dpng',t_names)
    close all

figure
[down_pow_prev_up_dur_cor(d),down_pow_prev_up_dur_p(d)] = scatter_with_cor(z_tot_down_pow{d},up_state_dur{d}(1:end-1));
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\prev_up_dur_v_down_pow_60_200' f_names{d}];
    print('-dpng',t_names)
    close all
    
    
end

cd C:\WC_Germany\JMM_analysis_pyr
save up_down_spectra_data_individ f S* *dur_cor *dur_p