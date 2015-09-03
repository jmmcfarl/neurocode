function [gamma] = get_hsmm_gamma(hmm,emiss)
%

K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
windowSize = hmm.windowSize;
windowSlide = hmm.windowSlide;

min_state_dur = hmm.min_state_dur;
dur_range = hmm.dur_range;
D = length(dur_range);

total_dur = T/Fs;
numWins = floor((total_dur-windowSize)/windowSlide);
win_t = (0:numWins-1)*windowSlide+windowSize/2;
t_axis = (1:T)/Fs;
win_t = [0 win_t max(t_axis)];
niqf = Fs/2;
if hmm.hcf > 0
    [filt_b,filt_a] = butter(2,hmm.hcf/niqf,'low');
end

A = [0 1; 1 0]; %can eliminate A entirely, but will leave in for now

ALPHA = zeros(K,D);
bmx=zeros(K,T);
S=zeros(K,T);
E=zeros(K,T);
BETA=ones(K,D);
Ex=ones(K,D);
Sx=ones(K,D);
gamma=zeros(K,T);
Pest=zeros(K,D);
Qest=zeros(T,1);
% Aest = zeros(K,K);


Pi = hmm.Pi';
P = hmm.P;
B = mvgauss_obslike_varmean(emiss,hmm);


%    starttime=clock;
%++++++++++++++++++     Forward     +++++++++++++++++
%---------------    Initialization    ---------------
ALPHA(:)=0;
ALPHA=repmat(Pi,1,D).*P;		%Equation (13)
r=B(1,:)*sum(ALPHA,2);			%Equation (3)
bmx(:,1)=B(1,:)./r;				%Equation (2)
E(:) = 0;
E(:,1)=bmx(:,1).*ALPHA(:,1);		%Equation (5)
S(:) = 0;
S(:,1)=A'*E(:,1);			%Equation (6)

%---------------    Induction    ---------------
for t=2:T
    ALPHA=[repmat(S(:,t-1),1,D-1).*P(:,1:D-1) + ...
        repmat(bmx(:,t-1),1,D-1).*ALPHA(:,2:D) , S(:,t-1).*P(:,D)];		%Equation (12)
    r=(B(t,:)*sum(ALPHA,2));		%Equation (3)
    bmx(:,t)=B(t,:)./r;			%Equation (2)
    E(:,t)=bmx(:,t).*ALPHA(:,1);		%Equation (5)
    S(:,t)=A'*E(:,t);				%Equation (6)
end
%++++++++ To check if the likelihood is increased ++++++++
%     if ir>1
%         %    clock-starttime
%         if (lkh-lkh1)/T<0.001
%             disp('converged')
%             break
%         end
%     end
%++++++++ Backward and Parameter Restimation ++++++++
%---------------    Initialization    ---------------

%     Aest(:)=0;
%     Aest=E(:,T)*ones(1,K);  %Since T_{T|T}(m,n) = E_{T}(m) a_{mn}
Pest(:)=0;
gamma(:)=0;
gamma(:,T)=bmx(:,T).*sum(ALPHA,2);
% Best(:)=0;
% Best(:,O(T))=gamma;
[X,Qest(T)]=max(gamma(:,T));
BETA=repmat(bmx(:,T),1,D);				%Equation (7)
Ex=sum(P.*BETA,2);					%Equation (8)
Sx=A*Ex;						%Equation (9)

%---------------    Induction    ---------------
for t=(T-1):-1:1
    %% for estimate of A
    %         Aest=Aest+E(:,t)*Ex';
    %% for estimate of P
    Pest=Pest+repmat(S(:,t),1,D).*BETA;
    %% for estimate of state at time t
    gamma(:,t)=gamma(:,t+1)+E(:,t).*Sx-S(:,t).*Ex;
    gamma(gamma(:,t)<0,t)=0; % eleminate errors due to inaccurace of the computation.
    [X,Qest(t)]=max(gamma(:,t));
    %     %% for estimate of B
    %     Best(:,O(t))=Best(:,O(t))+gamma;
    BETA=repmat(bmx(:,t),1,D).*[Sx,BETA(:,1:D-1)];	%Equation (14)
    Ex=sum(P.*BETA,2);					%Equation (8)
    Sx=A*Ex;						%Equation (9)
end

