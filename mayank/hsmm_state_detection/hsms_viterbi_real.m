function [state_seq,max_lik] = hsms_viterbi_real(hmm,emiss)

K=hmm.K;
T=hmm.T;
p=hmm.p;
Fs=hmm.Fs;

%set model parameters to their initialized values
Pi = hmm.Pi;

B = mvgauss_uds_obslike(emiss,hmm);
eps = 1e-300;
B(B<eps) = eps;
psi = cumsum(log(B));
min_state_dur = hmm.min_state_dur;
max_state_dur = hmm.max_state_dur;

%log state dependent duration probabilities
p_d = log(hmm.P');

length_to = -Inf*ones(T,K);
coming_from = zeros(T,K);
length_to(1,:) = Pi;

%for t = 1
edge_set = 1+min_state_dur:1+max_state_dur;
duration_prob = p_d(edge_set-1,:);
obs_prob = psi(edge_set-1,:) - repmat(log(B(1,:)),length(edge_set),1);
edge_lengths = duration_prob+obs_prob;
length_to(edge_set,1) = length_to(1,2) + edge_lengths(:,2);
coming_from(edge_set,2) = 1;
length_to(edge_set,2) = length_to(1,1) + edge_lengths(:,1);
coming_from(edge_set,1) = 1;

for t = 2:T
    edge_set = t+min_state_dur:min(T,t+max_state_dur);
    duration_prob = p_d(edge_set-t,:);
    obs_prob = psi(edge_set-1,:) - repmat(psi(t-1,:),length(edge_set),1);
    edge_lengths = duration_prob+obs_prob;

    %for ups
    longer_path = find(length_to(edge_set,1) <= length_to(t,2) + edge_lengths(:,2));
    length_to(edge_set(longer_path),1) = length_to(t,2) + edge_lengths(longer_path,2);
    coming_from(edge_set(longer_path),1) = t;

    %for down
    longer_path = find(length_to(edge_set,2) <= length_to(t,1) + edge_lengths(:,1));
    length_to(edge_set(longer_path),2) = length_to(t,1) + edge_lengths(longer_path,1);
    coming_from(edge_set(longer_path),2) = t;
end

coming_from(coming_from < 1) = 1;

%%
[max_lik,end_state] = max(length_to(end,:));
cur_trans = T;
cur_state = end_state;
prev_trans = T;
state_seq = ones(T,1);
while cur_trans > 1
    cur_trans = coming_from(cur_trans,cur_state);
    cur_state = setdiff(1:2,cur_state);
    state_seq(cur_trans:prev_trans-1) = cur_state;
    prev_trans = cur_trans;
end

