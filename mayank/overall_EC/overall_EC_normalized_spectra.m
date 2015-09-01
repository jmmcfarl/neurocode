clear all
close all

dsf = 4;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
[b,a] = butter(2,[0.1/niqf 200/niqf]);

load C:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('C:\Code\WC_anal\general\')
addpath('C:\WC_Germany\persistent_2010\')

used_data = [sup_mec sup_lec];

pkfreqs = [];
pkpows = [];
ds_vars = [];

f_i = linspace(0,10,500);

%%
for d = used_data
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    disp(sprintf('session %d',d))
    cd(cdir);
    load ./used_data wcv_minus_spike
    
    wcv = filtfilt(b,a,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    %         wcv = wcv/sess_data(d).gains(1);
    wcv = zscore(wcv);
    df = Fsd/length(wcv);
    t_axis = (1:length(wcv))/Fsd;
    
    clear params
    params.Fs = Fsd;
    params.tapers = [4 7];
    params.err = 0;
    params.fpad = 0;
    win = [1 1];
    W = params.tapers(1)/win(1);
    
    [S,f] = mtspectrumsegc(wcv(:),win,params);
    res_f = 0:df:(Fsd/2);
    one_S = interp1(f,S,res_f);
    lines = [50 150];
    line_f = [];
    for i = 1:length(lines)
        line_f = [line_f find(res_f >= lines(i) - 1.4*W& res_f <= lines(i) + 1.4*W)];
    end
    noline_f = setdiff(1:length(res_f),line_f);
    one_S_r = interp1(res_f(noline_f),one_S(noline_f),res_f);
    two_S = [one_S_r one_S_r(end-1:-1:2)];
    acorr_seq = fftshift(ifft(two_S));
    [~,midpt] = max(acorr_seq);
    acorr_seq = acorr_seq(midpt:end);
    ord = 2;
    A = levinson(acorr_seq(1:ord),ord);
    est_x = filter(-A,1,wcv);
    
    acorr_seq2 = xcov(wcv,ord);
    acorr_seq2 = acorr_seq2(ord+1:end);
    A = levinson(acorr_seq2(1:ord),ord);
    est_x2 = filter(-A,1,wcv);
    
    params.tapers = [3 5];
    params.err = [2 .05];
    params.fpad = 0;
    params.fpass = [0.02 150];
    win = [20 20];
    W = params.tapers(1)/win(1);
    
    
    [S,f,~,~,Serr] = mtspectrumsegc(wcv(:),win,params);
    %         [S_est,f,~,~,S_esterr] = mtspectrumsegc(est_x(:),win,params);
    [S_est2,f,~,~,S_est2err] = mtspectrumsegc(est_x2(:),win,params);
    
    df = f(2)-f(1);
    outside_W = round(3*W/df);
    peakrange = find(f >= 1.5 & f <= 4.5);
    lS = 10*log10(S_est2);
    lS_sm = jmm_smooth_1d_cor(lS,7);
    lS_orig_sm = jmm_smooth_1d_cor(10*log10(S),7);
    
    interp_lS_sm(d,:) = interp1(f,lS_sm,f_i);
    interp_lS_orig_sm(d,:) = interp1(f,lS_orig_sm,f_i);
    
    [peak_theta_pow(d),peak_theta_loc(d)] = max(lS_sm(peakrange));
    peak_theta_freq(d) = f(peakrange(peak_theta_loc(d)));
    dlS_sm = [0 diff(lS_sm)];
    [~,peaklocs] = findpeaks(dlS_sm);
    [~,valleylocs] = findpeaks(-dlS_sm);
    pL = peakrange(peak_theta_loc(d));
    lE = peaklocs(find(peaklocs < pL,1,'last'));
    rE = valleylocs(find(valleylocs > pL,1,'first'));
    peakwidth = (rE-lE)/2;
    
    lp = pL-2*peakwidth;
    rp = pL+2*peakwidth;
    if lp < peakrange(1)
        lp = peakrange(1);
    end
    peak_lci = 10*log10(S_est2err(1,pL));
    above_uci = 10*log10(S_est2err(2,rp));
    below_uci = 10*log10(S_est2err(2,lp));
    
    if peak_lci > above_uci & peak_lci > below_uci
        sig_theta_peak(d) = 1;
    else
        sig_theta_peak(d) = 0;
    end
    
    
    figure('visible','off')
    %         plot_vector(S,f,'l',Serr)
    %         hold on
    %         plot_vector(S_est,f,'l',S_esterr,'r')
    %         shg
    plot_vector(S_est2,f,'l',S_est2err,'k')
    hold on
    line([2 2],[-60 -30])
    line([4.5 4.5],[-60 -30])
    if sig_theta_peak(d) == 1
        plot(peak_theta_freq(d),peak_lci,'ro','linewidth',2)
    else
        plot(peak_theta_freq(d),peak_lci,'ko','linewidth',2)
    end
    plot(f(rp),above_uci,'ko','linewidth',2)
    plot(f(lp),below_uci,'ko','linewidth',2)
    set(gca,'xscale','log')
    xlim([0.02 150])
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer(1),'_',sess_data(d).name);
    f_name = ['C:\WC_Germany\overall_EC\normalized_mp_spectrum\' s_name];
    print(f_name,'-dpng'), close
    
    
end

cd C:\WC_Germany\lec_theta\
save theta_classifications sig_theta_peak peak_* interp* f_i f

%%
l2mec = intersect(layer2,mec);
l2lec = intersect(layer2,lec);

plot_mec_lec_comparison_fig(f_i,interp_lS_orig_sm,sup_mec,sup_lec);
xlim([0 8])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Log Relative Power (AU)','Fontsize',12,'fontname','Arial')
ylim([-31 -1])

plot_mec_lec_comparison_fig(f_i,interp_lS_sm,sup_mec,sup_lec);
xlim([0 8])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Log Relative Power (AU)','Fontsize',12,'fontname','Arial')
ylim([-52 -43])

plot_mec_lec_comparison_fig(f_i,interp_lS_sm,l2lec,l3lec);
xlim([0 8])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Log Relative Power (AU)','Fontsize',12,'fontname','Arial')
ylim([-52 -43])

plot_mec_lec_comparison_fig(f_i,interp_lS_sm,l2lec,l2mec);
xlim([0 8])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Log Relative Power (AU)','Fontsize',12,'fontname','Arial')
ylim([-52 -43])

plot_mec_lec_comparison_fig(f_i,interp_lS_sm,l3lec,l3mec);
xlim([0 8])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Log Relative Power (AU)','Fontsize',12,'fontname','Arial')
ylim([-52 -43])

% figure
% h = errorbar(f_i,nanmean(interp_lS_orig_sm(l2mec,:)),nanstd(interp_lS_orig_sm(l2mec,:))/sqrt(length(l2mec)));
% errorbar_tick(h,.001,'units');
% hold on
% h = errorbar(f_i,nanmean(interp_lS_orig_sm(l2lec,:)),nanstd(interp_lS_orig_sm(l2lec,:))/sqrt(length(l2lec)),'g');
% errorbar_tick(h,.001,'units');


%%
sup_lec_fan = intersect(sup_lec,fan);
sup_lec_pyr = intersect(sup_lec,pyramidal);
sup_lec_mul = intersect(sup_lec,multipolar);
sup_lec_nan = setdiff(sup_lec,[sup_lec_fan sup_lec_pyr sup_lec_mul]);
used_set1 = l2mec;
used_set2 = l3mec;

figure
h = errorbar(f_i,nanmean(interp_lS_orig_sm(used_set1,:)),nanstd(interp_lS_orig_sm(used_set1,:))/sqrt(length(used_set1)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(f_i,nanmean(interp_lS_orig_sm(used_set2,:)),nanstd(interp_lS_orig_sm(used_set2,:))/sqrt(length(used_set2)),'g');
errorbar_tick(h,.001,'units');

figure
h = errorbar(f_i,nanmean(interp_lS_sm(used_set1,:)),nanstd(interp_lS_sm(used_set1,:))/sqrt(length(used_set1)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(f_i,nanmean(interp_lS_sm(used_set2,:)),nanstd(interp_lS_sm(used_set2,:))/sqrt(length(used_set2)),'g');
errorbar_tick(h,.001,'units');