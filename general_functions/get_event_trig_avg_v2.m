function [ev_avg,lags,ev_std,n_events] = get_event_trig_avg_v2(sig,event_inds,backlag,forwardlag,nboot,trial_ids,inc_prev_trial,use_inds)
%
% [ev_avg,lags] = get_event_trig_avg(sig,event_inds,backlag,forwardlag)
%

%revised may 14 2013 to handle multivariate signals and NAN values

min_nevents = 3;

if nargin < 5
    nboot = [];
end
if nargin < 6
    trial_ids = [];
end
if isempty(nboot) & nargout > 2
    error('Need to specify number of bootstrap samples for computing error bars');
end
if nargin < 7
    inc_prev_trial = 0;
end
if nargin < 8
    use_inds = [];
end

orig_size = size(sig);
if orig_size(2) > 1 && orig_size(1) == 1
    sig = sig';
end
orig_size = size(sig);
if length(orig_size) > 2
    sig = reshape(sig,orig_size(1),orig_size(2)*orig_size(3));
end
[NT,p] = size(sig);

%if using a subset of the data
if ~isempty(use_inds)
    sig = sig(use_inds,:);
    trial_ids = trial_ids(use_inds);
    event_inds = find(ismember(use_inds,event_inds));
    NT = length(use_inds);
end

lags = (-backlag:forwardlag)';

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
    event_inds(event_trial_ids == 0 | isinf(event_trial_ids)) = []; %don't use events with 0 trial id
    event_trial_ids = trial_ids(event_inds);
    
    check_trials = 1;
else
    check_trials = 0;
end

n_events = length(event_inds);

if n_events < min_nevents
    ev_avg = nan(length(lags),p);
    ev_std = nan(length(lags),p);
    return
end

    ev_mat = nan(n_events,length(lags),p);
    for i = 1:n_events
        %     fprintf('%d of %d\n',i,n_events);
        cur_ids = (event_inds(i)-backlag):(event_inds(i)+forwardlag);
        if check_trials
            uset = true(size(cur_ids));
            if ~inc_prev_trial
            uset(trial_ids(cur_ids(uset)) ~= event_trial_ids(i)) = false;
            else
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
