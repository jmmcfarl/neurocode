function [state_seq,llik_best] = hmm_uds_viterbi(hmm,emiss)

K=hmm.K;
T=hmm.T;
p=hmm.p;

Pi = hmm.Pi;
A = hmm.A;

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

for i = 1:hmm.Nsegs
    curT = length(UDS_seg_inds{i});
    %initialize variables (we will work with log probabilities to prevent
    %underflow)
    delta=zeros(curT,K); %rabiner, maximized probability along a particular state sequence
    psi=zeros(curT,K); %rabiner, most likely state that led to a given state in the sequence

    B = mvgauss_uds_obslike(emiss(UDS_seg_inds{i},:),hmm,i);

    %initialization
    delta(1,:) = log(Pi(i,:))+log(B(1,:));    % Eq. 105(a) Rabiner

    for t=2:curT
        for k = 1:K
            temp = delta(t-1,:) + log(A(:,k))';
            [delta(t,k),psi(t,k)] = max(temp);
        end
        delta(t,:) = delta(t,:) + log(B(t,:));
    end

    % Backtracking for Viterbi decoding
    [llik_best,state_seq{i}(curT)] = max(delta(curT,:));
    for t=curT-1:-1:1,
        state_seq{i}(t) = psi(t+1,state_seq{i}(t+1));
    end
end
