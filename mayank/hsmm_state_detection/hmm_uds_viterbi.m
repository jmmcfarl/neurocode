function [state_seq,llik_best] = hmm_uds_viterbi(hmm,emiss)

K=hmm.K;
T=hmm.T;
p=hmm.p;

Pi = hmm.Pi;
A = hmm.A;

%initialize variables (we will work with log probabilities to prevent
%underflow)
delta=zeros(T,K); %rabiner, maximized probability along a particular state sequence
psi=zeros(T,K); %rabiner, most likely state that led to a given state in the sequence

B = mvgauss_uds_obslike(emiss,hmm);

%initialization
delta(1,:) = log(Pi)+log(B(1,:));    % Eq. 105(a) Rabiner

for t=2:T
    for k = 1:K
        temp = delta(t-1,:) + log(A(:,k))';
        [delta(t,k),psi(t,k)] = max(temp);
    end
    delta(t,:) = delta(t,:) + log(B(t,:));
end

% Backtracking for Viterbi decoding
[llik_best,state_seq(T)] = max(delta(T,:));
for t=T-1:-1:1,
    state_seq(t) = psi(t+1,state_seq(t+1));
end

