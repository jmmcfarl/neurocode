load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
% load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

d = 7;
load(f_loc{d})
dat_name = f_loc{d}(25:end);

eval(['data = ' dat_name '_MP;'])
eval(['time = ' dat_name '_sampletimes;'])
Fs = 5000;
niqf = Fs/2;
[b,a] = butter(2,2/niqf,'high');

sweep_t = find(diff(time) < 0);
low = min(data);
high = (max(data))-0.3

maxlag = 0.5*Fs;
lags = -maxlag:maxlag;
midpt = round(length(lags)/2);
winsize = 1.0*Fs;
winslide = 0.05*Fs;

for i = 1:length(sweep_t)-1
    
    dataseg = data(sweep_t(i)+1:sweep_t(i+1));
    filtseg = filtfilt(b,a,dataseg);
    dfilt = [0;diff(filtseg)];
    dfilt = abs(jmm_smooth_1d(dfilt.^2,50));

    
    %     [Pxx(i,:),f] = pwelch(dataseg,[],[],[],5000);
   datadur = length(dataseg);
    numWins = floor((datadur-winsize)/winslide);
    xcov_mat = zeros(numWins,length(lags));
    for t = 1:numWins
       winbeg = (t-1)*winslide+1;
       winend = winbeg+winsize;
       xcov_mat(t,:) = xcov(filtseg(winbeg:winend),maxlag,'coeff');
    end
[curmaxval,curmaxloc] = max(xcov_mat(:,midpt:end));
curmaxloc = curmaxloc+midpt;
%% optional plotting
subplot(2,1,1)
    plot(time(sweep_t(i)+1:sweep_t(i+1)),dataseg)
%     hold on
%     plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg,'r','linewidth',2)
    ylim([low high])
    grid
    subplot(2,1,2)
    pcolor(linspace(0,time(end),numWins),lags/Fs,xcov_mat');shading flat
    hold on
%     plot(linspace(0,time(end),numWins),lags(curmaxloc)/Fs,'k')
    ylim([0 maxlag/Fs])
    pause
    clf
%%


end





