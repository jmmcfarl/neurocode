function [ev_avg,lags,ev_std,n_events,ev_mat,ev_cis] = get_event_trig_avg_v3(sig,event_inds,backlag,forwardlag,nboot,trial_ids,inc_prev_trial)
%
% [ev_avg,lags] = get_event_trig_avg(sig,event_inds,backlag,forwardlag)
%

%revised may 14 2013 to handle multivariate signals and NAN values

min_nevents = 3; %minimum number of events where we will even compute a triggered avg

if nargin < 5
    nboot = [];
end
if nargin < 6
    trial_ids = [];
end
if isempty(nboot) && nargout > 2
    error('Need to specify number of bootstrap samples for computing error bars');
end
if nargin < 7
    inc_prev_trial = false;
end

orig_size = size(sig);
%if row vector, flip
if orig_size(2) > 1 && orig_size(1) == 1
    sig = sig';
end

%convert signal into 2d array with time as first dimension
orig_size = size(sig);
if length(orig_size) > 2
    sig = reshape(sig,orig_size(1),orig_size(2)*orig_size(3));
end
[NT,p] = size(sig);


lags = (-backlag:forwardlag)';

%get rid of events that happen within the lag-range of the end points
bad_ids = find(event_inds <= backlag);
if ~isempty(bad_ids)
    %     fprintf('Dropping %d early events\n',length(bad_ids));
    event_inds(bad_ids) = [];
end
bad_ids = find(event_inds >= NT - forwardlag);
if ~isempty(bad_ids)
    %     fprintf('Dropping %d late events\n',length(bad_ids));
    event_inds(bad_ids) = [];
end

%if only using within-trial data
if ~isempty(trial_ids)
    event_trial_ids = trial_ids(event_inds);
    event_inds(event_trial_ids == 0 | isinf(event_trial_ids)) = []; %don't use events with 0 or inf trial id
    event_trial_ids = trial_ids(event_inds);
    
    check_trials = 1;
else
    check_trials = 0;
end

n_events = length(event_inds);

%check that we have at least the minimum number of events to work with
if n_events < min_nevents
    ev_avg = nan(length(lags),p);
    ev_std = nan(length(lags),p);
    return
end

if isempty(nboot) || nargout < 3 %if no error bars are requested
    ev_avg = zeros(length(lags),p);
    nan_ind = isnan(sig);
    cnt = zeros(length(lags),p); %need a counter to store how many useable data points for each lag
    for i = 1:n_events
        cur_ids = (event_inds(i)-backlag):(event_inds(i)+forwardlag);
        uset = true(length(cur_ids),p);
        if check_trials
            if ~inc_prev_trial %exclude data from any trial other than the current one
                uset(trial_ids(cur_ids) ~= event_trial_ids(i),:) = false;
            else %exclude data from any trial after the current one
                uset(trial_ids(cur_ids) > event_trial_ids(i),:) = false;
            end
        end
        uset(nan_ind(cur_ids,:)) = false;
        temp_sig = sig(cur_ids,:);
%         temp_sig(~uset) = nan;
%         ev_avg = ev_avg + temp_sig;
ev_avg(uset) = ev_avg(uset) + temp_sig(uset);
        cnt(uset) = cnt(uset) + 1;
    end
    ev_avg = ev_avg./cnt;
    
    if length(orig_size) > 2
        ev_avg = reshape(ev_avg,[length(lags) orig_size(2) orig_size(3)]);
    end
    ev_std = nan;
else %if error bars are requested
    ev_mat = nan(n_events,length(lags),p);
    for i = 1:n_events
        %     fprintf('%d of %d\n',i,n_events);
        cur_ids = (event_inds(i)-backlag):(event_inds(i)+forwardlag);
        if check_trials
            uset = true(size(cur_ids));
            if ~inc_prev_trial %exclude data from any trial other than the current one
                uset(trial_ids(cur_ids(uset)) ~= event_trial_ids(i)) = false;
            else %exclude data from any trial after the current one
                uset(trial_ids(cur_ids(uset)) > event_trial_ids(i)) = false;
            end
            ev_mat(i,uset,:) = sig(cur_ids(uset),:);
        else
            ev_mat(i,:,:) = sig(cur_ids,:);
        end
    end
    ev_mat = reshape(ev_mat,n_events,length(lags)*p);
    boot_avgs = bootstrp(nboot,@nanmean,ev_mat);
    boot_avgs = reshape(boot_avgs,[nboot length(lags) p]);
    ev_avg = squeeze(mean(boot_avgs));
    ev_std = squeeze(std(boot_avgs));
    ev_cis = prctile(boot_avgs,[2.5 97.5]);
end