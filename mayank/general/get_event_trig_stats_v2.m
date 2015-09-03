function [sta,U,sing_vals,V,lags] = get_event_trig_stats_v2(events,datamat,backlag,forwardlag,Fs,svnum)

% %INPUTS:
% events = event sample indices
% datamat = TxN matrix of continuous time series.  T is the number of time samples
%     and N is the number of signals
% backlag = number of samples to consider prior to events
% forwardlag = number of samples to consider after event
% Fs = sample frequency
% svnum = number of singular vectors to plot (output only)
% 
% % OUTPUTS:
% sta = matrix of event triggered average signals
% U = matrix of left-singular vectors
% sing_vals = vector of singular values
% V = matrix of right-singular vectors
% lags = vector of time


lags = -backlag:forwardlag;
datalen = size(datamat,1);
num_sigs = size(datamat,2);

%% compute STA
sta_mat = zeros(length(events),length(lags),num_sigs);
for i = 1:length(events)
    if events(i) > length(lags) && datalen-events(i) > length(lags)
        sta_mat(i,:,:) = datamat(events(i)-backlag:events(i)+forwardlag,:);
    else
        sta_mat(i,:,:) = nan;
    end
end

%% plot STA
sta = squeeze(nanmean(sta_mat));
figure(1)
hold on
cmap = colormap(jet(num_sigs));
for i = 1:num_sigs
    plot(lags/Fs,sta(:,i),'Color',cmap(i,:))
end
xlabel('Time (s)')
ylabel('Amplitude')

%% Compute PCA
pcalen = length(lags)*num_sigs;
pca_mat = zeros(length(events),pcalen);
bad_cols = [];
resh_mean = reshape(sta,pcalen,1);
resh_mean_norm = sum(resh_mean.^2);
for i = 1:length(events)
    if events(i) > length(lags) && datalen-events(i) > length(lags)
        curd = reshape(datamat(events(i)-backlag:events(i)+forwardlag,:),pcalen,1);
        pca_mat(i,:) = curd - curd'*resh_mean*resh_mean/resh_mean_norm;
    else
        bad_cols = [bad_cols i];
    end
end

pca_mat(bad_cols,:) = [];

[COEFF,SCORE,latent] = princomp(pca_mat);