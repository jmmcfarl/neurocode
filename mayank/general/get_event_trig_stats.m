function [sta,U,sing_vals,V,lags] = get_event_trig_stats(events,datamat,backlag,forwardlag,Fs,svnum)

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

%% Compute SVD
svdlen = length(lags)*num_sigs;
svd_mat = zeros(length(events),svdlen);
bad_cols = [];
for i = 1:length(events)
    if events(i) > length(lags) && datalen-events(i) > length(lags)
        svd_mat(i,:) = reshape(datamat(events(i)-backlag:events(i)+forwardlag,:),svdlen,1);
    else
        bad_cols = [bad_cols i];
    end
end

svd_mat(bad_cols,:) = [];

%mean subtract data
for i = 1:svdlen
    svd_mat(:,i) = svd_mat(:,i) - mean(svd_mat(:,i));
end

[U,S,V] = svd(svd_mat,0);

%extract singular values
sing_vals = zeros(size(S,1),1);
for i = 1:size(S,1)
    sing_vals(i)=S(i,i);
end

%% plot SVD results
%plot singular values
figure(2)
plot(sing_vals,'.')
xlabel('Singular Value Number')
ylabel('Singular Value')

%plot singular vectors
figure(3)
hold on
for i = 1:svnum
    cur_dat = reshape(V(:,i),length(lags),num_sigs);
    plot(lags/Fs,cur_dat)
    pause
    clf
end

%plot projection of events onto singular vectors
events(bad_cols) = [];
figure(4)
cmap = colormap(jet(svnum));
hold on
for i = 1:svnum
    plot(events/Fs,U(:,i),'.','Color',cmap(i,:))
end
