function [hmm,gamma,hmm_window_post]=hsmm_uds_train(hmm,emiss,numIt)

%% taken from Yu Kobayashi 2006

K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
if strcmp(hmm.meantype,'variable')
    windowSize = hmm.windowSize;
    windowSlide = hmm.windowSlide;
    total_dur = T/Fs;
    numWins = floor((total_dur-windowSize)/windowSlide);
    win_t = (0:numWins-1)*windowSlide+windowSize/2;
    t_axis = (1:T)/Fs;
    win_t = [0 win_t max(t_axis)];
    niqf = Fs/2;
    if hmm.hcf > 0
        [filt_b,filt_a] = butter(2,hmm.hcf/niqf,'low');
    end
end

min_state_dur = hmm.min_state_dur;
dur_range = hmm.dur_range;
D = length(dur_range);

A = [0 1; 1 0]; %can eliminate A entirely, but will leave in for now

ALPHA = zeros(K,D); %forward recursion variable
bmx=zeros(K,T);
S=zeros(K,T);
E=zeros(K,T);
BETA=ones(K,D); %backward recursion variable
Ex=ones(K,D);
Sx=ones(K,D);
gamma=zeros(K,T); %hidden state posterior
Pest=zeros(K,D);
Qest=zeros(T,1);
% Aest = zeros(K,K);

% lkh1 = -Inf;
ir1=numIt;

for ir=1:ir1
    
    Pi = hmm.Pi';
    P = hmm.P;
    B = mvgauss_uds_obslike(emiss,hmm);
     
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
    %++++++++ Backward and Parameter Restimation ++++++++
    %---------------    Initialization    ---------------
    Pest(:)=0;
    gamma(:)=0;
    gamma(:,T)=bmx(:,T).*sum(ALPHA,2);
    [X,Qest(T)]=max(gamma(:,T));
    BETA=repmat(bmx(:,T),1,D);				%Equation (7)
    Ex=sum(P.*BETA,2);					%Equation (8)
    Sx=A*Ex;						%Equation (9)
    
    %---------------    Induction    ---------------
    for t=(T-1):-1:1
        %% for estimate of P
        Pest=Pest+repmat(S(:,t),1,D).*BETA;
        %% for estimate of state at time t
        gamma(:,t)=gamma(:,t+1)+E(:,t).*Sx-S(:,t).*Ex;
        gamma(gamma(:,t)<0,t)=0; % eleminate errors due to inaccurace of the computation.
        [X,Qest(t)]=max(gamma(:,t));
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
        cur_dur_pars = hmm.state(k).dur_pars;
        if strcmp(hmm.state(k).dur_type,'gamma')
            exp_par = [cur_dur_pars(1)/cur_dur_pars(2); psi(cur_dur_pars(1))-log(cur_dur_pars(2))];
            exp_np = [sum(dur_range.*P_np(k,:)); sum(log(dur_range).*P_np(k,:))];
            [est_theta(k,:)] = fminsearch(@(x)gamma_exp_diff(x,exp_np),cur_dur_pars);
        elseif strcmp(hmm.state(k).dur_type,'inv_gauss')
            exp_par(1) = cur_dur_pars(1);
            exp_par(2) = sum(inverse_gaussian_pmf(dur_range,cur_dur_pars(1),cur_dur_pars(2))./dur_range);
            exp_np(1) = sum(P_np(k,:).*dur_range);
            exp_np(2) = sum(P_np(k,:)./dur_range);
            [est_theta(k,:)] = fminsearch(@(x)invgauss_exp_diff(x,exp_np,dur_range),cur_dur_pars);
        else
            error('invalid duration distribution')
        end
    end
    
    if strcmp(hmm.meantype,'variable')
        temp_meanfun = nan(2,numWins);
        tot_prob = nan(2,numWins);
        p_thresh = 0.05;
        for l = 1:K
            for w = 1:numWins
                begT = round((w-1)*windowSlide*Fs)+1;
                endT = begT + round(windowSize*Fs);
                temp_meanfun(l,w) = sum(gamma(l,begT:endT).*emiss(begT:endT,1)',2)/sum(gamma(l,begT:endT),2);
                tot_prob(l,w) = sum(gamma(l,begT:endT),2)/(windowSize*Fs);
            end
            temp_meanfun(l,tot_prob(l,:) < p_thresh) = nan;
        end
        meandiff = temp_meanfun(2,:)-temp_meanfun(1,:);
        avg_meandiff = nanmean(meandiff);
        temp_meanfun(1,isnan(temp_meanfun(1,:))) = temp_meanfun(2,isnan(temp_meanfun(1,:))) - avg_meandiff;
        temp_meanfun(2,isnan(temp_meanfun(2,:))) = temp_meanfun(1,isnan(temp_meanfun(2,:))) + avg_meandiff;
        
        temp_meanfun = [temp_meanfun(:,1) temp_meanfun temp_meanfun(:,end)]';
        meanfun_est = interp1(win_t,temp_meanfun,t_axis);
        for l = 1:K
            meanfun_est(:,l) = filtfilt(filt_b,filt_a,meanfun_est(:,l));
        end
        if p > 1
            for l = 1:K
                fixedmean_est(l,:) = sum(repmat(gamma(l,:),p-1,1).*emiss(:,2:end)',2)/sum(gamma(l,:),2);
            end
        end
    else
        for l = 1:K
            fixedmean = sum(gamma(l,:).*emiss(:,1)',2)/sum(gamma(l,:),2);
        end
    end
    
    for l = 1:K
        if strcmp(hmm.meantype,'variable') && p > 1
            state_mean_est = [meanfun_est(:,l) repmat(fixedmean_est(l,:),T,1)];
        elseif strcmp(hmm.meantype,'variable')
            state_mean_est = meanfun_est(:,l);
        else
            state_mean_est = repmat(fixedmean_est(l,:),T,1);
        end
        mdiff = emiss-state_mean_est;
        var_est(l,:,:) = ((repmat(gamma(l,:),p,1)'.*mdiff)'*mdiff)/sum(gamma(l,:),2);
    end    
    
    %update all parameters in hmm
    hmm.Pi = Pi';
    for l = 1:K
        %for state covariance parameters
        new_covar = squeeze(var_est(l,:,:));
        if det(new_covar) > 0
            hmm.state(l).var = squeeze(var_est(l,:,:));
        else
            disp('warning: singular covariance')
        end
        %for state mean parameters
        if strcmp(hmm.meantype,'variable') 
            if p > 1
                hmm.state(l).fixedmean = fixedmean_est(l,:);
            end
            hmm.state(l).meanfun = meanfun_est(:,l);
        else
            hmm.state(l).fixedmean = fixedmean_est(l,:);
        end
        %for duration model parameters
        if strcmp(hmm.state(l).dur_type,'gamma')
            hmm.P(l,:) = gamma_pmf(dur_range,est_theta(l,1),est_theta(l,2));
        elseif strcmp(hmm.state(l).dur_type,'inv_gauss')
            hmm.P(l,:) = inverse_gaussian_pmf(dur_range,est_theta(l,1),est_theta(l,2));
        end
        if min_state_dur > 1
            hmm.P(l,1:min_state_dur-1) = 0;
        end
        hmm.state(l).dur_pars = est_theta(l,:);
    end
    ir
end

%for outputting the windowed average state posteriors (only for
%nonstationary case
if strcmp(hmm.meantype,'variable')
    hmm_window_post = tot_prob;
else
    hmm_window_post = nan;
end