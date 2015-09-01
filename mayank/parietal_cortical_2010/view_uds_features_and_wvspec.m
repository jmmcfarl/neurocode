function [] = view_uds_features_and_wvspec(sig1,sig2,specgram,freqs,Fs,winsize)

t_axis = (1:length(sig1))/Fs;
T = t_axis(end);
numWins = round(T/winsize);

specgram = log(specgram);
specgram = specgram - repmat(mean(specgram,2),1,length(t_axis));
specgram = specgram./repmat(std(specgram,[],2),1,length(t_axis));



winsize = round(Fs*winsize);
for i = 1:numWins
    inds = ((i-1)*winsize+1):i*winsize;
    inds(inds > length(t_axis)) = [];
    subplot(2,1,1)
    plot(t_axis(inds),sig1(inds),t_axis(inds),sig2(inds),'r')
    subplot(2,1,2)
    pcolor(t_axis(inds),freqs,specgram(:,inds));shading flat
%     set(gca,'yscale','log')
    caxis([-2 2])
    line([t_axis(inds(1)) t_axis(inds(end))],[3 3],'color','w','linewidth',1)
    pause
    clf
end