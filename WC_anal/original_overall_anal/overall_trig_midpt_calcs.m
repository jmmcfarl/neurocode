clear all
close all
% 
% 
load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data
% beg_lag = find(lags>-2,1,'first');
% end_lag = find(lags>2,1,'first');

for d = 1:length(over_dir)
  
    %mp ups
    [mp_utrig_mp_maxpt(d),maxloc] = max(mp_utrig_mp(d,:));
    mp_utrig_mp_minpt(d) = min(mp_utrig_mp(d,1:maxloc));
    midpt = (mp_utrig_mp_minpt(d)+mp_utrig_mp_maxpt(d))/2;
    mp_utrig_mp_mid(d) = find(mp_utrig_mp(d,1:maxloc)<midpt,1,'last');
    
        [mp_utrig_lf8_maxpt(d),maxloc] = max(mp_utrig_lf8(d,:));
    mp_utrig_lf8_minpt(d) = min(mp_utrig_lf8(d,1:maxloc));
    midpt = (mp_utrig_lf8_minpt(d)+mp_utrig_lf8_maxpt(d))/2;
    mp_utrig_lf8_mid(d) = find(mp_utrig_lf8(d,1:maxloc)<midpt,1,'last');
    
%         [mp_utrig_lf3_maxpt(d),maxloc] = max(mp_utrig_lf3(d,:));
%     mp_utrig_lf3_minpt(d) = min(mp_utrig_lf3(d,1:maxloc));
%     midpt = (mp_utrig_lf3_minpt(d)+mp_utrig_lf3_maxpt(d))/2;
%     mp_utrig_lf3_mid(d) = find(mp_utrig_lf3(d,1:maxloc)<midpt,1,'last');
 
        [mp_utrig_lf2_maxpt(d),maxloc] = max(mp_utrig_lf2(d,:));
    mp_utrig_lf2_minpt(d) = min(mp_utrig_lf2(d,1:maxloc));
    midpt = (mp_utrig_lf2_minpt(d)+mp_utrig_lf2_maxpt(d))/2;
    mp_utrig_lf2_mid(d) = find(mp_utrig_lf2(d,1:maxloc)<midpt,1,'last');
    
    %mp downs
        [mp_dtrig_mp_minpt(d),minloc] = min(mp_dtrig_mp(d,:));
    [mp_dtrig_mp_maxpt(d),maxloc] = max(mp_dtrig_mp(d,1:minloc));
    midpt = (mp_dtrig_mp_minpt(d)+mp_dtrig_mp_maxpt(d))/2;
    mp_dtrig_mp_mid(d) = find(mp_dtrig_mp(d,1:minloc)>midpt,1,'last');
    
    [mp_dtrig_lf8_minpt(d),minloc] = min(mp_dtrig_lf8(d,:));
    [mp_dtrig_lf8_maxpt(d),maxloc] = max(mp_dtrig_lf8(d,1:minloc));
    midpt = (mp_dtrig_lf8_minpt(d)+mp_dtrig_lf8_maxpt(d))/2;
    mp_dtrig_lf8_mid(d) = find(mp_dtrig_lf8(d,1:minloc)>midpt,1,'last');
    
%             [mp_dtrig_lf3_minpt(d),minloc] = min(mp_dtrig_lf3(d,:));
%     [mp_dtrig_lf3_maxpt(d),maxloc] = max(mp_dtrig_lf3(d,1:minloc));
%     midpt = (mp_dtrig_lf3_minpt(d)+mp_dtrig_lf3_maxpt(d))/2;
%     mp_dtrig_lf3_mid(d) = find(mp_dtrig_lf3(d,1:minloc)>midpt,1,'last');
 
            [mp_dtrig_lf2_minpt(d),minloc] = min(mp_dtrig_lf2(d,:));
    [mp_dtrig_lf2_maxpt(d),maxloc] = max(mp_dtrig_lf2(d,1:minloc));
    midpt = (mp_dtrig_lf2_minpt(d)+mp_dtrig_lf2_maxpt(d))/2;
    mp_dtrig_lf2_mid(d) = find(mp_dtrig_lf2(d,1:minloc)>midpt,1,'last');
    
    
     %lfp ups
        [lf8_utrig_mp_maxpt(d),maxloc] = max(lf8_utrig_mp(d,:));
    lf8_utrig_mp_minpt(d) = min(lf8_utrig_mp(d,1:maxloc));
    midpt = (lf8_utrig_mp_minpt(d)+lf8_utrig_mp_maxpt(d))/2;
    lf8_utrig_mp_mid(d) = find(lf8_utrig_mp(d,1:maxloc)<midpt,1,'last');
    
        [lf8_utrig_lf8_maxpt(d),maxloc] = max(lf8_utrig_lf8(d,:));
    lf8_utrig_lf8_minpt(d) = min(lf8_utrig_lf8(d,1:maxloc));
    midpt = (lf8_utrig_lf8_minpt(d)+lf8_utrig_lf8_maxpt(d))/2;
    lf8_utrig_lf8_mid(d) = find(lf8_utrig_lf8(d,1:maxloc)<midpt,1,'last');
    
%         [lf8_utrig_lf3_maxpt(d),maxloc] = max(lf8_utrig_lf3(d,:));
%     lf8_utrig_lf3_minpt(d) = min(lf8_utrig_lf3(d,1:maxloc));
%     midpt = (lf8_utrig_lf3_minpt(d)+lf8_utrig_lf3_maxpt(d))/2;
%     lf8_utrig_lf3_mid(d) = find(lf8_utrig_lf3(d,1:maxloc)<midpt,1,'last');
 
        [lf8_utrig_lf2_maxpt(d),maxloc] = max(lf8_utrig_lf2(d,:));
    lf8_utrig_lf2_minpt(d) = min(lf8_utrig_lf2(d,1:maxloc));
    midpt = (lf8_utrig_lf2_minpt(d)+lf8_utrig_lf2_maxpt(d))/2;
    lf8_utrig_lf2_mid(d) = find(lf8_utrig_lf2(d,1:maxloc)<midpt,1,'last');
    
    %lfp downs
           [lf8_dtrig_mp_minpt(d),minloc] = min(lf8_dtrig_mp(d,:));
    [lf8_dtrig_mp_maxpt(d),maxloc] = max(lf8_dtrig_mp(d,1:minloc));
    midpt = (lf8_dtrig_mp_minpt(d)+lf8_dtrig_mp_maxpt(d))/2;
    lf8_dtrig_mp_mid(d) = find(lf8_dtrig_mp(d,1:minloc)>midpt,1,'last');
    
            [lf8_dtrig_lf8_minpt(d),minloc] = min(lf8_dtrig_lf8(d,:));
    [lf8_dtrig_lf8_maxpt(d),maxloc] = max(lf8_dtrig_lf8(d,1:minloc));
    midpt = (lf8_dtrig_lf8_minpt(d)+lf8_dtrig_lf8_maxpt(d))/2;
    lf8_dtrig_lf8_mid(d) = find(lf8_dtrig_lf8(d,1:minloc)>midpt,1,'last');
    
%             [lf8_dtrig_lf3_minpt(d),minloc] = min(lf8_dtrig_lf3(d,:));
%     [lf8_dtrig_lf3_maxpt(d),maxloc] = max(lf8_dtrig_lf3(d,1:minloc));
%     midpt = (lf8_dtrig_lf3_minpt(d)+lf8_dtrig_lf3_maxpt(d))/2;
%     lf8_dtrig_lf3_mid(d) = find(lf8_dtrig_lf3(d,1:minloc)>midpt,1,'last');
 
            [lf8_dtrig_lf2_minpt(d),minloc] = min(lf8_dtrig_lf2(d,:));
    [lf8_dtrig_lf2_maxpt(d),maxloc] = max(lf8_dtrig_lf2(d,1:minloc));
    midpt = (lf8_dtrig_lf2_minpt(d)+lf8_dtrig_lf2_maxpt(d))/2;
    lf8_dtrig_lf2_mid(d) = find(lf8_dtrig_lf2(d,1:minloc)>midpt,1,'last');
   d
    
end


save C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_midpts *minpt *maxpt *mid lags