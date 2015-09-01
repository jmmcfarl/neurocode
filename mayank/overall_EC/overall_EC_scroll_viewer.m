clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir.mat

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;

lcf = 0.05;
hcf = 20;

lcf2 = 2;
hcf2 = 20;

lcf3 = 2;
hcf3 = 6;

lcf4 = 0.05;
hcf4 = 1;

winSize = 15;

lec = find_struct_field_vals(sess_data,'region','LEC');
sess_data = sess_data(lec);


d = 20;

cdir = sess_data(d).directory;
cdir(1) = 'G';
disp(sprintf('session %d',d))
cd(cdir);

load used_data lf8 wcv_minus_spike

% lf3_f = get_lf_features_caus(lf3,Fs,Fsd,[lcf2 hcf2]);
lf8_f = get_lf_features(lf8,Fs,Fsd,[lcf hcf]);
wcv_f = get_lf_features(wcv_minus_spike,Fs,Fsd,[lcf hcf]);
lf8_f2 = get_lf_features(lf8,Fs,Fsd,[lcf2 hcf2]);
t = (1:length(lf8_f))/Fsd;

% niqf = Fs/2;
% [b,a] = butter(3,[lcf3 hcf3]/niqf);
% lf8_f3 = filtfilt(b,a,lf8);
% lf8_f3 = zscore(downsample(lf8_f3,dsf));
% [b,a] = butter(2,[lcf4 hcf4]/niqf);
% lf8_f4 = filtfilt(b,a,lf8);
% lf8_f4 = zscore(downsample(lf8_f4,dsf));


% lf8_h = hilbert(lf8_f2);
% lf8_a = abs(lf8_h);
% lf8_p = angle(lf8_h);
% lf3_h = hilbert(lf3_f);
% lf3_a = abs(lf3_h);
% lf3_p = angle(lf3_h);

dt = round(winSize*Fsd);
ct = 1;
et = ct + dt;
T = round(t(end)*Fsd);
while ct < T
    plot(t(ct:et),wcv_f(ct:et)), hold on
    plot(t(ct:et),lf8_f(ct:et),'r','linewidth',1), hold on
    plot(t(ct:et),lf8_f2(ct:et)/3-2,'k')
%     plot(t(ct:et),lf7_f(ct:et),'color',[0.6 0.4 0.2],'linewidth',1), 
%         plot(t(ct:et),lf8_f2(ct:et)/3-2,'r','linewidth',1)
% if ~isnan(sess_data(d).gains(6))
%     plot(t(ct:et),lf6_f(ct:et),'k','linewidth',1)
%         plot(t(ct:et),lf6_f2(ct:et)/3-2,'r','linewidth',1)
% end
% plot(t(ct:et),lf5_f(ct:et)/3-2,'k','linewidth',1)
%     plot(t(ct:et),lf8_f2(ct:et)/3-2,'k')
%     plot(t(ct:et),lf5_f2(ct:et)/3-2,'g')
%     plot(t(ct:et),lf3_f(ct:et)/3-4,'k')
%     plot(t(ct:et),lf8_f3(ct:et)/3-4,'k')
%     pause
%     plot(t(ct:et),lf3_f(ct:et)-1,'k')
    xlim([t(ct) t(et)])
%     ylim([-5.5 2.5])
    ct = ct + dt;
    et = ct + dt;
    pause
    clf
end