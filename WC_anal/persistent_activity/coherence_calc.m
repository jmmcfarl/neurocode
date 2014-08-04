clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\desynch_extract\desynch_points

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 45];
window = [100 100];
niqf = 2016/2;
hcf1 = 100/niqf;
[b1,a1] = butter(2,hcf1,'low');

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data wcv_minus_spike lf8
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);
    
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
     
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
    if ~isempty(desynch_start{d})
        sMarkers = [[1;desynch_stop{d}'] [desynch_start{d}';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end
    
[Cmn(d,:),Phimn(d,:),Smn,Smm,f,ConfC(d),PhiStd,Cerr(d,:,:)] = coherencyc_unequal_length_trials([down_w down_8],window, params, sMarkers );
    
%     [Sw(d,:),f,varS,C,Swerr(d,:,:)] = mtspectrumsegc(down_w,window,params);
%     [S8(d,:),f,varS,C,S8err(d,:,:)] = mtspectrumsegc(down_8,window,params);
%     
%     plot(f,Cmn(d,:),'linewidth',2)
%     hold on
%     line([0 45],[ConfC(d) ConfC(d)],'Color','k')
%     xlim([0 45])
%     tname = ['C:\WC_Germany\Persistent_activity\coherence\wband_' f_names{d}];
%     print('-dpng',tname);
%     close
%     
%     plot(f,Cmn(d,:),'linewidth',2)
%     hold on
%     plot(f,squeeze(Cerr(d,1,:)),'--')
%     plot(f,squeeze(Cerr(d,2,:)),'--')
%     line([0 1],[ConfC(d) ConfC(d)],'Color','k')
%     xlim([0 1])
%     tname = ['C:\WC_Germany\Persistent_activity\coherence\' f_names{d}];
%     print('-dpng',tname);
%     close

    
    clear down_w down_8 wcv* lf8 
    
end


save C:\WC_Germany\Persistent_activity\coherence\coherence_data Cmn f Phimn Cerr ConfC