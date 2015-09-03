clear all
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\Code\smoothing\software\')
addpath('G:\Code\general\')
addpath('G:\WC_Germany\hsmm_state_detection\\')

ac_Fs = 2016;
ac_dsf = 16;
ac_Fsd = 2016/ac_dsf;
ac_niqf = 2016/2;
[b,a] = butter(2,40/ac_niqf,'low');

dc_target_fs = 200;

amp_range = linspace(-85,-20,500);
windowSize = 20;
windowSlide = 2;
movingwin = [windowSize windowSlide];

for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    if ~isempty(sess_data(d).heka_dir)
        load used_data wcv_minus_spike
        load aligned_heka
        dc_Fs = 1/(median(diff(dc_time)));
        dc_dsf = round(dc_Fs/dc_target_fs);
        dc_data = downsample(dc_data,dc_dsf);
        if mean(abs(dc_data)) < 1
            dc_data = dc_data*100;
        end

        dc_time = downsample(dc_time,dc_dsf);
        ac_time = downsample(ac_time,ac_dsf);
        
        wcvf = filtfilt(b,a,wcv_minus_spike);
        wcvf = downsample(wcvf,ac_dsf);
        wcvf = zscore(wcvf);
        df = ac_Fsd/length(wcvf);
        
        clear params
        params.Fs = ac_Fsd;
        params.tapers = [4 7];
        params.err = 0;
        params.fpad = 0;
        win = [1 1];
        W = params.tapers(1)/win(1);       
        [S,f] = mtspectrumsegc(wcvf(:),win,params);
        res_f = 0:df:(ac_Fsd/2);
        one_S = interp1(f,S,res_f);
        lines = [50 150];
        line_f = [];
        for i = 1:length(lines)
            line_f = [line_f find(res_f >= lines(i) - 1.4*W& res_f <= lines(i) + 1.4*W)];
        end
        noline_f = setdiff(1:length(res_f),line_f);
        one_S_r = interp1(res_f(noline_f),one_S(noline_f),res_f);
        fit_range = find(res_f >= 0.01 & res_f < 150);
        temp_p = polyfit(log10(res_f(fit_range)),log10(one_S_r(fit_range)),1);
% 
params.Fs = ac_Fsd;
        params.tapers = [2 3];
        params.err = 0;
        ac_movingwin = [10 2];
        params.fpass = [0.01 15];
        [Sg,tg,fg]=mtspecgramc(wcvf(:),ac_movingwin,params);
        logSg = log10(Sg);
        logSp = polyval(temp_p,log10(fg));
        logSp_r = repmat(logSp,length(tg),1);
        logS_r = logSg-logSp_r;
        
        total_dur = ac_time(end);
        numWins = floor((total_dur-windowSize)/windowSlide);
        t_axis = (0:numWins-1)*windowSlide+windowSize/2;
        data_dist = zeros(numWins,length(amp_range));
        for w = 1:numWins
            begT = (w-1)*windowSlide;
            endT = begT + windowSize;
            used_inds = find(dc_time > begT & dc_time < endT);
            data_seg = dc_data(used_inds);
            data_dist(w,:) = gpkde(data_seg(:),1,[-85 -20 500]);
        end


        f1 = figure('visible','off');
        set(f1,'papersize',[2.0 2.5]);
        subplot(2,1,1)
        pcolor(t_axis,amp_range,log(data_dist'));shading flat
        caxis([-10 -2])
        subplot(2,1,2)
        pcolor(tg,fg,log10(Sg)');shading flat
%         caxis([-1.5 1])
%         set(gca,'yscale','log')
        tname = ['G:\WC_Germany\overall_EC\normalized_mp_spectrogram\theta_' s_name];
        print('-dpng',tname);
        close        
        %%
      
    end
end

cd G:\WC_Germany\overall_EC\
% save overall_EC_heka_UDS_data upstate* downstate*