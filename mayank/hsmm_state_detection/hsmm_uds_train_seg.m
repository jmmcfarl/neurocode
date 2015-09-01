function [hmm,gamma,hmm_window_post,bad_model]=hsmm_uds_train_seg(hmm,emiss)

tolerance = 1e-5;
cur_log_post = -Inf;
delta_log_post = Inf;

%% taken from Yu Kobayashi 2006, modified from source code by JMM
%original function: "hsmm_new.m"
% Author: Shun-Zheng Yu
% Available: http://sist.sysu.edu.cn/~syu/Publications/hsmm_new.m% 
% Ref: Practical Implementation of an Efficient Forward-Backward Algorithm for an Explicit Duration Hidden Markov Model
% by Shun-Zheng Yu and H. Kobayashi
% IEEE Transactions on Signal Processing, Vol. 54, No. 5, MAY 2006, pp. 1947-1951  
%  This program is free software; you can redistribute it and/or
%  modify it under the terms of the GNU General Public License
%  as published by the Free Software Foundation; either version
%  2 of the License, or (at your option) any later version.
%  http://www.gnu.org/licenses/gpl.txt


K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
if strcmp(hmm.meantype,'variable')
    windowSize = hmm.windowSize;
    windowSlide = hmm.windowSlide;
    p_thresh = 0.05;
end

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

min_state_dur = hmm.min_state_dur;
dur_range = hmm.dur_range;
D = length(dur_range);

A = [0 1; 1 0]; %can eliminate A entirely, but will leave in for now
eps = 1e-200;

ir = 0;
LP = [];
while delta_log_post > tolerance
    
    Pi = hmm.Pi;
    P = hmm.P;

    Pest=zeros(K,D);

    for i = 1:hmm.Nsegs

        curT = length(UDS_seg_inds{i});
        ALPHA = zeros(K,D); %forward recursion variable
        bmx=zeros(K,curT);
        S=zeros(K,curT);
        E=zeros(K,curT);
        BETA=ones(K,D); %backward recursion variable
        Ex=ones(K,D);
        Sx=ones(K,D);
        gamma{i}=zeros(K,curT); %hidden state posterior
        Qest=zeros(curT,1);
        
        [B,bad_model] = mvgauss_uds_obslike_robust(emiss(UDS_seg_inds{i},:),hmm,i);
%         B(B < eps) = eps; %prevent likelihoods from rounding to 0
%         %find instances where the emissions model fails and correct the
%         %likelihood
%         false_ups = find(emiss(UDS_seg_inds{i},:) < hmm.state(1).meanfun{i} & ...
%             B(:,1) < B(:,2));
%         B(false_ups,1) = 1;
%         false_downs = find(emiss(UDS_seg_inds{i},:) > hmm.state(2).meanfun{i} & ...
%             B(:,1) > B(:,2));
%         B(false_downs,2) = 1;
%         hmm.model_failure{i} = [false_ups; false_downs];

    
        %++++++++++++++++++     Forward     +++++++++++++++++
        %---------------    Initialization    ---------------
        ALPHA(:)=0;
        ALPHA=repmat(Pi(i,:)',1,D).*P;		%Equation (13)
        r=B(1,:)*sum(ALPHA,2);			%Equation (3)
        bmx(:,1)=B(1,:)./r;				%Equation (2)
        E(:) = 0;
        E(:,1)=bmx(:,1).*ALPHA(:,1);		%Equation (5)
        S(:) = 0;
        S(:,1)=A'*E(:,1);			%Equation (6)
        llik=log(r);

        %---------------    Induction    ---------------
        for t=2:curT
            ALPHA=[repmat(S(:,t-1),1,D-1).*P(:,1:D-1) + ...
                repmat(bmx(:,t-1),1,D-1).*ALPHA(:,2:D) , S(:,t-1).*P(:,D)];		%Equation (12)
            r=(B(t,:)*sum(ALPHA,2));		%Equation (3)
            bmx(:,t)=B(t,:)./r;			%Equation (2)
            E(:,t)=bmx(:,t).*ALPHA(:,1);		%Equation (5)
            S(:,t)=A'*E(:,t);				%Equation (6)
            llik=llik+log(r);
        end
        %++++++++ Backward and Parameter Restimation ++++++++
        %---------------    Initialization    ---------------
%         Pest(:)=0;
%         gamma(:)=0;
        gamma{i}(:,curT)=bmx(:,curT).*sum(ALPHA,2);
        [X,Qest(curT)]=max(gamma{i}(:,curT));
        BETA=repmat(bmx(:,curT),1,D);				%Equation (7)
        Ex=sum(P.*BETA,2);					%Equation (8)
        Sx=A*Ex;						%Equation (9)

        %---------------    Induction    ---------------
        for t=(curT-1):-1:1
            %% for estimate of P
            Pest=Pest+repmat(S(:,t),1,D).*BETA;
            %% for estimate of state at time t
            gamma{i}(:,t)=gamma{i}(:,t+1)+E(:,t).*Sx-S(:,t).*Ex;
            gamma{i}(gamma{i}(:,t)<0,t)=0; % eleminate errors due to inaccurace of the computation.
            [X,Qest(t)]=max(gamma{i}(:,t));
            BETA=repmat(bmx(:,t),1,D).*[Sx,BETA(:,1:D-1)];	%Equation (14)
            Ex=sum(P.*BETA,2);					%Equation (8)
            Sx=A*Ex;						%Equation (9)
        end

        Pest=Pest+repmat(Pi(i,:)',1,D).*BETA;    %Since D_{1|T}(m,d) = \PAI(m) P_{m}(d) \Beta_{1}(m,d)

    end
    
    %RESTIMATE PARAMETERS
    for i = 1:hmm.Nsegs
        Pi(i,:)=gamma{i}(:,1)./sum(gamma{i}(:,1));
    end
    
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
        meandiff = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            total_dur = curT/Fs;
            t_axis{i} = (1:curT)/Fs;
            numWins = floor((total_dur-windowSize)/windowSlide);
            win_t{i} = (0:numWins-1)*windowSlide+windowSize/2;
            win_t{i} = [0 win_t{i} max(t_axis{i})];
            tot_prob{i} = nan(2,numWins);
            temp_meanfun{i} = nan(hmm.K,numWins);

            for l = 1:K
                for w = 1:numWins
                    begT = round((w-1)*windowSlide*Fs)+1;
                    endT = begT + round(windowSize*Fs);
                    temp_meanfun{i}(l,w) = sum(gamma{i}(l,begT:endT).*...
                        emiss(UDS_seg_inds{i}(begT:endT),1)',2)/sum(gamma{i}(l,begT:endT),2);
                    tot_prob{i}(l,w) = sum(gamma{i}(l,begT:endT),2)/(windowSize*Fs);
                end
                temp_meanfun{i}(l,tot_prob{i}(l,:) < p_thresh) = nan;
            end
            meandiff = [meandiff temp_meanfun{i}(2,tot_prob{i}(l,:) >= p_thresh)-...
                temp_meanfun{i}(1,tot_prob{i}(l,:) >= p_thresh)];
        end
        avg_meandiff = nanmean(meandiff);
        for i = 1:hmm.Nsegs
            temp_meanfun{i}(1,isnan(temp_meanfun{i}(1,:))) = temp_meanfun{i}(2,isnan(temp_meanfun{i}(1,:))) - avg_meandiff;
            temp_meanfun{i}(2,isnan(temp_meanfun{i}(2,:))) = temp_meanfun{i}(1,isnan(temp_meanfun{i}(2,:))) + avg_meandiff;
            temp_meanfun{i} = [temp_meanfun{i}(:,1) temp_meanfun{i} temp_meanfun{i}(:,end)]';
            meanfun_est{i} = interp1(win_t{i},temp_meanfun{i},t_axis{i});
        end
        if p > 1
            for l = 1:K
                o_sum = 0;
                g_sum = 0;
                for i = 1:hmm.Nsegs
                    o_sum = o_sum + sum(repmat(gamma{i}(l,:),p-1,1).*emiss(UDS_seg_inds{i},2:end),2);
                    g_sum = g_sum + sum(gamma{i}(l,:),2);
                end
                fixedmean_est{i}(:,l) = o_sum/g_sum;
            end
        end
    else
        for l = 1:K
            o_sum = 0;
            g_sum = 0;
            for i = 1:hmm.Nsegs
               o_sum = o_sum + sum(gamma{i}(l,:).*emiss(UDS_seg_inds{i},1)',2);
               g_sum = g_sum + sum(gamma{i}(l,:),2);
            end
            fixedmean_est(l) = o_sum/g_sum;
        end
    end
    
    for l = 1:K
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            if strcmp(hmm.meantype,'variable') && p > 1
                state_mean_est = [meanfun_est{i}(:,l) repmat(fixedmean_est(l,:),curT,1)];
            elseif strcmp(hmm.meantype,'variable')
                state_mean_est = meanfun_est{i}(:,l);
            else
                state_mean_est = repmat(fixedmean_est(:,l),curT,1);
            end
            mdiff{i} = emiss(UDS_seg_inds{i},:)-state_mean_est;
        end
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            o_sum = o_sum + (repmat(gamma{i}(l,:),p,1)'.*mdiff{i})'*mdiff{i};
            g_sum = g_sum + sum(gamma{i}(l,:),2);
        end
        var_est(l,:,:) = o_sum/g_sum;
    end    
    
    %update all parameters in hmm
    hmm.Pi = Pi;
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
                hmm.state(l).fixedmean = fixedmean_est(:,l);
            end
            for i = 1:hmm.Nsegs
                hmm.state(l).meanfun{i} = meanfun_est{i}(:,l);
            end
        else
            hmm.state(l).fixedmean = fixedmean_est(:,l);
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
    
    LP=[LP; llik];
    delta_log_post = (llik - cur_log_post)/abs(llik);
    cur_log_post = llik;
    fprintf('iteration %i log posterior = %f \n',ir,llik);
end

%for outputting the windowed average state posteriors (only for
%nonstationary case
if strcmp(hmm.meantype,'variable')
    hmm_window_post = tot_prob;
else
    hmm_window_post = nan;
end
hmm.hsmm_LP = LP;