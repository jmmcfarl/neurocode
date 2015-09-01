function [] = view_uds_features(sigmat,Fs,winsize)

t_axis = (1:size(sigmat,1))/Fs;
T = t_axis(end);
numWins = round(T/winsize);

num_sigs = size(sigmat,2);
cmap = [0.1 0.1 0.1; 0 0.1 0.8; 0.7 0 0.1];

winsize = round(Fs*winsize);
for i = 1:numWins
    inds = ((i-1)*winsize+1):i*winsize;
    inds(inds > length(t_axis)) = [];
    for i = 1:num_sigs
        plot(t_axis(inds),sigmat(inds,i),'color',cmap(i,:))
        hold on
    end
    
    pause
    clf
end