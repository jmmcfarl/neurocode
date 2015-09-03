function [state_seq,lik_best,delta] = hsms_viterbi_fixmean(hmm,emiss)

K=hmm.K;
T=hmm.T;
p=hmm.p;

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

alpha=zeros(T,K);
beta=zeros(T,K);
gamma=zeros(T,K);

%initialize viterbi bits
delta=zeros(T,K);
psi=zeros(T,K);

B=mvgauss_obslike_fixmean(emiss,hmm);

scale=zeros(T,1);
% Scaling for delta
dscale=zeros(T,1);

tiny=exp(-700);

alpha(1,:)=Pi(:)'.*B(1,:);
scale(1)=sum(alpha(1,:));
alpha(1,:)=alpha(1,:)/(scale(1)+tiny);
% For viterbi decoding
delta(1,:) = alpha(1,:);    % Eq. 32(a) Rabiner (1989)

for t=2:T
    alpha(t,:)=(alpha(t-1,:)*A).*B(t,:);
    scale(t)=sum(alpha(t,:));
    alpha(t,:)=alpha(t,:)/(scale(t)+tiny);

    for l=1:K,
        v=delta(t-1,:).*A(:,l)';
        mv=max(v);
        delta(t,l)=mv*B(t,l);  % Eq 33a Rabiner (1989)
        if length(find(v==mv)) > 1
            % no unique maximum - so pick one at random
            tmp1=find(v==mv);
            tmp2=rand(length(tmp1),1);
            [tmp3,tmp4]=max(tmp2);
            psi(t,l)=tmp4;
        else
            psi(t,l)=find(v==mv);  % ARGMAX; Eq 33b Rabiner (1989)
        end
    end

    % SCALING FOR DELTA ????
    dscale(t)=sum(delta(t,:));
    delta(t,:)=delta(t,:)/(dscale(t)+tiny);
end

% Get beta values for single state decoding
beta(T,:)=ones(1,K)/scale(T);
for t=T-1:-1:1
    beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
end

% Get gamma values for single state decoding
gamma=(alpha.*beta);
gamma=rdiv(gamma,rsum(gamma));
gammasum=sum(gamma);

chi=zeros(T-1,K*K);
for t=1:T-1
    temp=A.*( alpha(t,:)' * (beta(t+1,:).*B(t+1,:)));
    chi(t,:)=temp(:)'/sum(temp(:));
end

likv=sum(log(scale+(scale==0)*tiny));
lik_best=sum(log(dscale+(dscale==0)*tiny));

% Backtracking for Viterbi decoding
state_seq(T) = find(delta(T,:)==max(delta(T,:)),1);% Eq 34b Rabiner;
for t=T-1:-1:1,
    state_seq(t) = psi(t+1,state_seq(t+1));
end

% LP=sum(likv);
% LP_best=sum(lik_best);