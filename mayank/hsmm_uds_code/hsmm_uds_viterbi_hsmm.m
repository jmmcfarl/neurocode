function [state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hsmm,emiss)

% Ref: Datta R, Hu J, Ray B (2008) On efficient Viterbi decoding for hidden semi-
% Markov models. 19th International Conference on Pattern Recognition. pp 1–4.

% Compute Viterbi state sequence of emiss given the HSMM specified in structure hsmm
%
% Input: 
%       hmm: structure with the trained parameters of an hsmm
%       emiss: observation sequence
%              
% Output:
%       state_seq: cell array containing the ML state sequence for each UDS
%          segment
%       max_lik: the maximized log likelihood

K=hsmm.K;
T=hsmm.T;
p=hsmm.p;
Fs=hsmm.Fs;

%set model parameters to their initialized values
Pi = hsmm.Pi;

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hsmm.Nsegs
    UDS_seg_inds{i} = hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2)];
end
rel_prob = 0;
for i = 1:hsmm.Nsegs
    rel_prob = rel_prob + sum(hsmm.gamma{i});
end

for i = 1:hsmm.Nsegs
    curT = length(UDS_seg_inds{i});
%     B = mvgauss_uds_obslike_robust(emiss(UDS_seg_inds{i},:),hsmm,i);
%     B = mvgauss_uds_obslike_robust_v2(emiss(UDS_seg_inds{i},:),hsmm,rel_prob,i);
    B = mvgauss_uds_obslike_robust_v3(emiss(UDS_seg_inds{i},:),hsmm,rel_prob,i);
    
    psi = cumsum(log(B));
    min_state_dur = hsmm.min_state_dur;
    max_state_dur = hsmm.max_state_dur;

    %log state dependent duration probabilities
    p_d = log(hsmm.P');

    %initialize
    length_to = -Inf*ones(curT,K);
    coming_from = zeros(curT,K);
    length_to(1,:) = Pi(i,:);

    %for t = 1
    edge_set = 1+min_state_dur:1+max_state_dur; %possible next state transition times
    duration_prob = p_d(edge_set-1,:); %probability of each allowed state duration
    obs_prob = psi(edge_set-1,:) - repmat(log(B(1,:)),length(edge_set),1); %observation likelihoods
    edge_lengths = duration_prob+obs_prob; %total probability of each state transition
    length_to(edge_set,1) = length_to(1,2) + edge_lengths(:,2);
    coming_from(edge_set,2) = 1;
    length_to(edge_set,2) = length_to(1,1) + edge_lengths(:,1);
    coming_from(edge_set,1) = 1;

    %forward recursion
    for t = 2:curT
        edge_set = t+min_state_dur:min(curT,t+max_state_dur);
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

    %% backtrack max like state sequence
    [max_lik(i),end_state] = max(length_to(end,:));
    cur_trans = curT;
    cur_state = end_state;
    prev_trans = curT;
    state_seq{i} = ones(curT,1);
    while cur_trans > 1
        cur_trans = coming_from(cur_trans,cur_state);
        cur_state = setdiff(1:2,cur_state);
        state_seq{i}(cur_trans:prev_trans-1) = cur_state;
        prev_trans = cur_trans;
    end

end

max_lik = sum(max_lik);
