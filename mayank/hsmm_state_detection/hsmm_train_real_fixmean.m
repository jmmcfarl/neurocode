function [hmm,gamma]=hsmm_train_real_fixmean(hmm,emiss,numIt)
%

K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;

min_state_dur = hmm.min_state_dur;
dur_range = hmm.dur_range;
D = length(dur_range);

total_dur = T/Fs;
t_axis = (1:T)/Fs;
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


lkh1 = -Inf;
ir1=numIt;

for ir=1:ir1
    
    Pi = hmm.Pi';
    P = hmm.P;
    B = mvgauss_obslike_fixmean(emiss,hmm);
    
    
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
    lkh=log(r);
    
    %---------------    Induction    ---------------
    for t=2:T
        ALPHA=[repmat(S(:,t-1),1,D-1).*P(:,1:D-1) + ...
            repmat(bmx(:,t-1),1,D-1).*ALPHA(:,2:D) , S(:,t-1).*P(:,D)];		%Equation (12)
        r=(B(t,:)*sum(ALPHA,2));		%Equation (3)
        bmx(:,t)=B(t,:)./r;			%Equation (2)
        E(:,t)=bmx(:,t).*ALPHA(:,1);		%Equation (5)
        S(:,t)=A'*E(:,t);				%Equation (6)
        lkh=lkh+log(r);
    end
    %++++++++ To check if the likelihood is increased ++++++++
    %     if ir>1
    %         %    clock-starttime
    %         if (lkh-lkh1)/T<0.001
    %             disp('converged')
    %             break
    %         end
    %     end
    lkh1=lkh;
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
    
    Pest=Pest+repmat(Pi,1,D).*BETA;    %Since D_{1|T}(m,d) = \PAI(m) P_{m}(d) \Beta_{1}(m,d)
    
    %RESTIMATE PARAMETERS
    Pi=gamma(:,1)./sum(gamma(:,1));
    
    %for estimating state duration dist parameters
    Pest=Pest.*P;
    P_np=Pest./repmat(sum(Pest,2),1,D);
    
    for k = 1:K
        cur_alpha = hmm.state(k).alpha;
        cur_beta = hmm.state(k).beta;
        exp_par = [cur_alpha/cur_beta; psi(cur_alpha)-log(cur_beta)];
        exp_np = [sum(dur_range.*P_np(k,:)); sum(log(dur_range).*P_np(k,:))];
        [est_theta(k,:)] = fminsearch(@(x)gamma_exp_diff(x,exp_np),[cur_alpha cur_beta]);
    end
    

    
    for l = 1:K
        mean_est(l,:) = sum(repmat(gamma(l,:),p,1).*emiss',2)/sum(gamma(l,:),2);
        mdiff = emiss-repmat(mean_est(l,:),T,1);
        var_est(l,:,:) = ((repmat(gamma(l,:),p,1)'.*mdiff)'*mdiff)/sum(gamma(l,:),2);
    end
    
    
    %update all parameters in hmm
    hmm.Pi = Pi';
    for l = 1:K
        new_covar = squeeze(var_est(l,:,:));
        if det(new_covar) > 0
            hmm.state(l).var = squeeze(var_est(l,:,:));
        else
            disp('warning: singular covariance')
        end
        hmm.state(l).mean = mean_est(l,:);
        hmm.P(l,:) = gamma_pmf(dur_range,est_theta(l,1),est_theta(l,2));
        hmm.P(l,1:min_state_dur) = 0;
        hmm.state(l).alpha = est_theta(l,1);
        hmm.state(l).beta = est_theta(l,2);
    end
    
    ir
end
