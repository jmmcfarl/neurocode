function [state_seq,llik_best] = hmm_uds_viterbi_seg(hmm,emiss)

%AS in Rabiner 1989

K=hmm.K;
T=hmm.T;
p=hmm.p;

Pi = hmm.Pi;
A = hmm.A;
eps = 1e-200;

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

for i = 1:hmm.Nsegs
    curT = length(UDS_seg_inds{i});
    %initialize variables 
    delta=zeros(curT,K); %delta is the maximized log probability as a function of time
    psi=zeros(curT,K); %Psi stores the most likely preceeding state

    B = mvgauss_uds_obslike_robust(emiss(UDS_seg_inds{i},:),hmm,i);
    
    %initialization
    delta(1,:) = log(Pi(i,:))+log(B(1,:));    % Eq. 105(a) Rabiner

    for t=2:curT
        for k = 1:K
            temp = delta(t-1,:) + log(A(:,k))';
            [delta(t,k),psi(t,k)] = max(temp);
        end
        delta(t,:) = delta(t,:) + log(B(t,:));
    end

    % Backtracking Viterbi
    [llik_best,state_seq{i}(curT)] = max(delta(curT,:));
    for t=curT-1:-1:1,
        state_seq{i}(t) = psi(t+1,state_seq{i}(t+1));
    end
end
