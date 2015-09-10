function Xmat = create_event_Xmat(NT,event_inds,lag_set,trial_vec)
% Xmat = create_event_Xmat(NT,event_inds,lag_set,trial_vec)
% compile set of lagged indicator vectors into predictor matrix
% INPUTS:
%     NT: total number of time points
%     event_inds: index values of the events
%     lag_set: set of lag values (in units of bins)
%     <trial_vec>: [NTx1] vector of trial indices. Used to exclude indicators that get lagged across trial boundaries
% OUTPUTS:
%     Xmat: matrix of indicator vectors
    
if nargin < 4 || isempty(trial_vec)
   trial_vec = zeros(NT,1);  
end
n_lags = length(lag_set);
Xmat = zeros(NT,n_lags);
for ii = 1:n_lags
    cur_target_inds = event_inds + lag_set(ii); %all events offset by current lag
    uu = find(cur_target_inds > 1 & cur_target_inds < NT); %set of in-bound indices
    cur_target_inds = cur_target_inds(uu);
    cur_target_inds(trial_vec(cur_target_inds) ~= trial_vec(event_inds(uu))) = []; %if the lagged value isnt in the same trial, get rid of this
    Xmat(cur_target_inds,ii) = 1;
end
