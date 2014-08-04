clear all
cd C:\WC_Germany\Tetrode_recs\2008-12-14-B\2008-12-14_21-42-30
load lfp_data 

dsf = 4;
niqf = 2016/2;
Fsd = 2016/dsf;
hcf = 200/niqf;
[b,a] = butter(2,hcf,'low');

lf1_f = filtfilt(b,a,lf1_data);
lf12_f = filtfilt(b,a,lf12_data);
lf1_d = downsample(lf1_f,dsf);
lf12_d = downsample(lf12_f,dsf);

win = 100;

params.Fs = Fsd;
params.err = [2 0.05];
params.fpass = [0 30];

[S12,f,varS,C,S12err]=mtspectrumsegc(lf12_d,win,params);
[S1,f,varS,C,S1err]=mtspectrumsegc(lf1_d,win,params);

[C,phi,dummy,dummy,dummy,f]=coherencysegc(lf1_d',lf12_d',win,params);

plot(f,10*log10(S12),'linewidth',2)
hold on
plot(f,10*log10(S1),'r','linewidth',2)
legend('Cortical','MEC')
plot(f,10*log10(S12err(1,:)),'--')
plot(f,10*log10(S12err(2,:)),'--')
plot(f,10*log10(S1err(1,:)),'r--')
plot(f,10*log10(S1err(2,:)),'r--')
xlim([0 1])
print('-dpng','Power_spec')
xlim([0 10])
print('-dpng','Power_spec_wb')
close
figure
plot(f,C)
xlim([0 1])
print('-dpng','Coherence')
close
