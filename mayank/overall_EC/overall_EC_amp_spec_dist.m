clear all
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\Code\smoothing\software\')
addpath('G:\Code\general\')
addpath('G:\WC_Germany\hsmm_state_detection\\')

ac_Fs = 2016;
ac_dsf = 8;
ac_Fsd = 2016/ac_dsf;
ac_niqf = 2016/2;
[b,a] = butter(2,[0.05/ac_niqf 10/ac_niqf]);
[b2,a2] = butter(2,[0.05/ac_niqf 40/ac_niqf]);

dc_target_fs = 200;

params.Fs = ac_Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 30];
windowSize = 20;
windowSlide = 2;
movingwin = [windowSize windowSlide];

amp_range = linspace(-85,-20,500);
amp_range8 = linspace(-3,5,500);
amp_range3 = linspace(-3,4,500);

% lec = find_struct_field_vals(sess_data,'region','LEC');
for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    if ~isempty(sess_data(d).heka_dir)
        load used_data lf8 lf3 wcv_minus_spike
        load aligned_heka
        dc_Fs = 1/(median(diff(dc_time)));
        dc_dsf = round(dc_Fs/dc_target_fs);
                dc_data = downsample(dc_data,dc_dsf);
        if mean(abs(dc_data)) < 1
            dc_data = dc_data*100;
        end

        dc_time = downsample(dc_time,dc_dsf);
        ac_time = downsample(ac_time,ac_dsf);
        
        lf8f = filtfilt(b,a,lf8);
        lf8f = downsample(lf8f,ac_dsf);
        lf8f = zscore(lf8f);
        lf3f = filtfilt(b,a,lf3);
        lf3f = downsample(lf3f,ac_dsf);
        lf3f = zscore(lf3f);
        wcvf2 = filtfilt(b,a,wcv_minus_spike);
        wcvf2 = downsample(wcvf2,ac_dsf);
        wcvf2 = zscore(wcvf2);
        lf8f2 = filtfilt(b2,a2,lf8);
        lf8f2 = downsample(lf8f2,ac_dsf);
        lf8f2 = zscore(lf8f2);
        lf3f2 = filtfilt(b2,a2,lf3);
        lf3f2 = downsample(lf3f2,ac_dsf);
        lf3f2 = zscore(lf3f2);

        total_dur = ac_time(end);
        numWins = floor((total_dur-windowSize)/windowSlide);
        t_axis = (0:numWins-1)*windowSlide+windowSize/2;
        data_dist = zeros(numWins,length(amp_range));
        data_dist8 = zeros(numWins,length(amp_range));
        data_dist3 = zeros(numWins,length(amp_range));
        for w = 1:numWins
            begT = (w-1)*windowSlide;
            endT = begT + windowSize;
            used_inds = find(dc_time > begT & dc_time < endT);
            data_seg = dc_data(used_inds);
            data_dist(w,:) = gpkde(data_seg(:),1,[-85 -20 500]);
            begInd = round(ac_Fsd*begT+1);
            endInd = begInd + round(ac_Fsd*windowSize);
            data_seg = lf8f(begInd:endInd);
            data_dist8(w,:) = gpkde(data_seg(:),0.06,[-3 5 500]);
            data_seg = lf3f(begInd:endInd);
            data_dist3(w,:) = gpkde(data_seg(:),0.06,[-3 4 500]);
        end

% 
        [Pw,t,f] = mtspecgramc(wcvf2,movingwin,params);
        [P8,t,f] = mtspecgramc(lf8f2,movingwin,params);
        [P3,t,f] = mtspecgramc(lf3f2,movingwin,params);
% 
        f1 = figure('visible','off');
        set(f1,'papersize',[4.0 2.0]);
        subplot(3,2,1)
        pcolor(t_axis,amp_range,log(data_dist'));shading flat
        caxis([-10 -2])
        subplot(3,2,3)
        pcolor(t_axis,amp_range8,log(data_dist8'));shading flat
        caxis([-4 -0.5])
        subplot(3,2,5)
        pcolor(t_axis,amp_range3,log(data_dist3'));shading flat
        caxis([-4 -0.5])
        subplot(3,2,2)
        pcolor(t,f,log(Pw'));shading flat
        caxis([-8 0])
        set(gca,'yscale','log')
        subplot(3,2,4)
        pcolor(t,f,log(P8'));shading flat
        caxis([-8 0])
        set(gca,'yscale','log')
        subplot(3,2,6)
        pcolor(t,f,log(P3'));shading flat
        caxis([-8 0])
        set(gca,'yscale','log')
        tname = ['G:\WC_Germany\overall_EC\amp_spec_comb\' s_name];
        print('-dpng',tname);
        close
%         saveas(f1,tname,'fig')
%         close
        
        %%
        load ec_hmm_state_seq
        dc_data = dc_data(:);
        mp_state_seq = hmm_bbstate_seq;
        dc_upvals = [];
        dc_downvals = [];
        [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,252,length(wcvf2));
        for i = 1:length(mp_state_seq)
            up_trans = find(mp_state_seq{i}(1:end-1)==1 & mp_state_seq{i}(2:end)==2);
            down_trans = find(mp_state_seq{i}(1:end-1)==2 & mp_state_seq{i}(2:end)==1);
            up_trans(up_trans > down_trans(end)) = [];
            down_trans(down_trans < up_trans(1)) = [];
            up_trans = up_trans + new_seg_inds(i,1)-1;
            down_trans = down_trans + new_seg_inds(i,1)-1;
            up_trans = round(up_trans/(ac_dsf/8));
            down_trans = round(down_trans/(ac_dsf/8));
            for j = 1:length(up_trans)-1
               cur_vals = find(dc_time > ac_time(up_trans(j)) & dc_time < ac_time(down_trans(j)));
               dc_upvals = [dc_upvals; dc_data(cur_vals)];
               cur_vals = find(dc_time > ac_time(down_trans(j)) & dc_time < ac_time(up_trans(j+1)));
               dc_downvals = [dc_downvals; dc_data(cur_vals)];
            end
        end
        upstate_mean(d) = mean(dc_upvals);
        downstate_mean(d) = mean(dc_downvals);
        upstate_var(d) = var(dc_upvals);
        downstate_var(d) = var(dc_downvals);
        
    end
end

cd F:\WC_Germany\overall_EC\
save overall_EC_heka_UDS_data upstate* downstate*


%% For eps figure generation
figure
imagesc(t_axis,amp_range,log(data_dist'));shading flat
caxis([-10 -1.5])
ylim([-80 -20])

f_log = logspace(log10(f(5)),log10(f(end)),500);
lPw = 10*log10(Pw);
int_lPw = interp1(f,lPw',f_log);
figure
imagesc(t,f_log,int_lPw);shading flat
ylim([0.05 20])
caxis([-35 5.5])

figure
pcolor(t,f_log,int_lPw);shading flat
ylim([0.12 20])
caxis([-35 5.5])
set(gca,'yscale','log')

xl1 = 240;
xl2 = 258;
figure
plot(dc_time,dc_data)
xlim([xl1 xl2])

wcv_tf = get_lf_features(wcv_minus_spike,2016,252,[1.5 10]);
wcv_lf = get_lf_features(wcv_minus_spike,2016,252,[0.05 2]);
figure
plot(ac_time,zscore(wcv_tf),'k')
hold on
plot(ac_time,zscore(wcv_lf))
plot(ac_time,zscore(lf8f),'r')
xlim([xl1 xl2])

xl1 = 1145;
xl2 = 1180;
figure
plot(dc_time,dc_data)
xlim([xl1 xl2])

figure
plot(ac_time,zscore(lf8f),'r')
xlim([xl1 xl2])