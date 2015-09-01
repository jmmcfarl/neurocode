function [hmm,gamma,hmm_window_post]=hsmm_uds_train_hsmm_bi(hmm,emiss)

%% train parameters of a hidden semi-markov model for inferring UP-DOWN states.
%
% Input:
%       hmm:  structure containing initialized model parameters
%       emiss:  observation sequence

% Output:
%        hmm: structure containing best fit model parameters
%        gamma: sequence of marginal probabilities of being in each state
%        hmm_window_post: total probability of being in each state within
%           each window (for sliding window state-conditional means only)


%% taken from Yu Kobayashi 2006, modified by JMM
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

tolerance = 1e-5;
cur_log_post = -Inf;
delta_log_post = Inf;

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
    
    if ~exist('gamma','var')
        rel_prob = [0.5; 0.5];
    else
        rel_prob = [0; 0];
        for i = 1:hmm.Nsegs
            rel_prob = rel_prob + sum(gamma{i},2);
        end
        rel_prob = rel_prob/sum(rel_prob);
    end
    
    Pi = hmm.Pi;
    P = hmm.P;
    
    Pest=zeros(K,D);
    
    for i = 1:hmm.Nsegs
        
        [B,rob_model_inds{i}] = mvgauss_uds_obslike_robust_bi_v3(emiss(UDS_seg_inds{i},:),hmm,rel_prob,i);
        %                 [B,rob_model_inds{i}] = mvgauss_uds_obslike_robust_v2(emiss(UDS_seg_inds{i},:),hmm,rel_prob,i);
        flex_model_inds{i} = setdiff(1:length(UDS_seg_inds{i}),rob_model_inds{i});
        %         B = mvgauss_uds_obslike(emiss(UDS_seg_inds{i},:),hmm,i);
        %         flex_model_inds{i} = 1:length(UDS_seg_inds{i});
        %         rob_model_inds{i} = [];
        
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
        r = zeros(curT,1);
        
        %++++++++++++++++++     Forward     +++++++++++++++++
        %---------------    Initialization    ---------------
        ALPHA(:)=0;
        ALPHA=repmat(Pi(i,:)',1,D).*P;		%Equation (13)
        r(1)=B(1,:)*sum(ALPHA,2);			%Equation (3)
        bmx(:,1)=B(1,:)./r(1);				%Equation (2)
        E(:) = 0;
        E(:,1)=bmx(:,1).*ALPHA(:,1);		%Equation (5)
        S(:) = 0;
        S(:,1)=A'*E(:,1);			%Equation (6)
        
        %---------------    Induction    ---------------
        for t=2:curT
            ALPHA=[repmat(S(:,t-1),1,D-1).*P(:,1:D-1) + ...
                repmat(bmx(:,t-1),1,D-1).*ALPHA(:,2:D) , S(:,t-1).*P(:,D)];		%Equation (12)
            r(t)=(B(t,:)*sum(ALPHA,2));		%Equation (3)
            bmx(:,t)=B(t,:)./r(t);			%Equation (2)
            E(:,t)=bmx(:,t).*ALPHA(:,1);		%Equation (5)
            S(:,t)=A'*E(:,t);				%Equation (6)
        end
        llik = sum(log(r(flex_model_inds{i})));
        
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
            if ismember(t,flex_model_inds{i})
                Pest=Pest+repmat(S(:,t),1,D).*BETA;
            end
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
            %             exp_par = [cur_dur_pars(1)/cur_dur_pars(2); psi(cur_dur_pars(1))-log(cur_dur_pars(2))];
            exp_np = [sum(dur_range.*P_np(k,:)); sum(log(dur_range).*P_np(k,:))];
            options.Display = 'off';
            est_theta(k,:) = fsolve(@(x)hsmm_gamma_exp_fun(x,exp_np),cur_dur_pars,options);
        elseif strcmp(hmm.state(k).dur_type,'inv_gauss')
            %             exp_par(1) = cur_dur_pars(1);
            %             exp_par(2) = sum(inverse_gaussian_pmf(dur_range,cur_dur_pars(1),cur_dur_pars(2))./dur_range);
            exp_np(1) = sum(P_np(k,:).*dur_range);
            exp_np(2) = sum(P_np(k,:)./dur_range);
            options.Display = 'off';
            est_theta(k,:) = fsolve(@(x)hsmm_ig_exp_fun(x,exp_np,dur_range),cur_dur_pars,options);
        else
            error('invalid duration distribution')
        end
    end
    
    if strcmp(hmm.meantype,'variable')
        meandiff = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            %initialize
            tot_prob{i} = nan(hmm.K,length(UDS_seg_inds{i}));
            meanfun_est_d1{i} = zeros(hmm.K,length(UDS_seg_inds{i}));
            win_kern = ones(round(windowSize*Fs),1); %this is a convolution kernel
            buff_len = round(windowSize*Fs); %this is a buffer window to handle edge effects during convolution
            for l=1:K
                
                %for purposes of interpolation, make sure that we use the
                %first and last samples of the segment (otherwise we
                %attempt extrapolation0
                to_use_inds = flex_model_inds{i};
                if ~ismember(1,to_use_inds)
                    to_use_inds = [1 to_use_inds];
                end
                if ~ismember(length(UDS_seg_inds{i}),to_use_inds)
                    to_use_inds = [to_use_inds length(UDS_seg_inds{i})];
                end
                
                %establish vectors for convolution
                to_conv_num = gamma{i}(l,:)'.*emiss(UDS_seg_inds{i},1);
                to_conv_num = interp1(to_use_inds,to_conv_num(to_use_inds),1:length(to_conv_num));
                to_conv_den = gamma{i}(l,:)';
                to_conv_den = interp1(to_use_inds,to_conv_den(to_use_inds),1:length(to_conv_den));
                
                %add a buffer window onto either side of each vector to
                %minimize convolution edge artifacts
                to_conv_num = [fliplr(to_conv_num(1:buff_len)) to_conv_num fliplr(to_conv_num(end-buff_len+1:end))];
                to_conv_den = [fliplr(to_conv_den(1:buff_len)) to_conv_den fliplr(to_conv_den(end-buff_len+1:end))];
                use_vec = logical(ones(size(to_conv_num)));
                use_vec(1:buff_len) = 0; use_vec(end-buff_len+1:end) = 0;
                
                cnum = conv(to_conv_num,win_kern,'same');
                cden = conv(to_conv_den,win_kern,'same');
                meanfun_est_d1{i}(l,:) = cnum(use_vec)./cden(use_vec);
                tot_prob{i}(l,:) = cden(use_vec)/length(win_kern);
                meanfun_est_d1{i}(l,tot_prob{i}(l,:) < p_thresh) = nan; %set any estimates of meanfunction with insufficient posterior probability to nan
            end
            meandiff = [meandiff meanfun_est_d1{i}(2,tot_prob{i}(l,:) >= p_thresh)-...
                meanfun_est_d1{i}(tot_prob{i}(l,:) >= p_thresh)]; %compute a vector of difference of state meanfunctions
        end
        avg_meandiff = nanmean(meandiff); %average difference of mean functions
        %any instances where there was insufficient posterior probability
        %for one state are given by estimates extrapolated from the other
        %state's mean
        for i = 1:hmm.Nsegs
            meanfun_est_d1{i}(1,isnan(meanfun_est_d1{i}(1,:))) = meanfun_est_d1{i}(2,isnan(meanfun_est_d1{i}(1,:))) - avg_meandiff;
            meanfun_est_d1{i}(2,isnan(meanfun_est_d1{i}(2,:))) = meanfun_est_d1{i}(1,isnan(meanfun_est_d1{i}(2,:))) + avg_meandiff;
            meanfun_est_d1{i} = meanfun_est_d1{i}';
        end
        
        meandiff = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            %initialize
            tot_prob{i} = nan(hmm.K,length(UDS_seg_inds{i}));
            meanfun_est_d2{i} = zeros(hmm.K,length(UDS_seg_inds{i}));
            win_kern = ones(round(windowSize*Fs),1); %this is a convolution kernel
            buff_len = round(windowSize*Fs); %this is a buffer window to handle edge effects during convolution
            for l=1:K
                
                %for purposes of interpolation, make sure that we use the
                %first and last samples of the segment (otherwise we
                %attempt extrapolation0
                to_use_inds = flex_model_inds{i};
                if ~ismember(1,to_use_inds)
                    to_use_inds = [1 to_use_inds];
                end
                if ~ismember(length(UDS_seg_inds{i}),to_use_inds)
                    to_use_inds = [to_use_inds length(UDS_seg_inds{i})];
                end
                
                %establish vectors for convolution
                to_conv_num = gamma{i}(l,:)'.*emiss(UDS_seg_inds{i},2);
                to_conv_num = interp1(to_use_inds,to_conv_num(to_use_inds),1:length(to_conv_num));
                to_conv_den = gamma{i}(l,:)';
                to_conv_den = interp1(to_use_inds,to_conv_den(to_use_inds),1:length(to_conv_den));
                
                %add a buffer window onto either side of each vector to
                %minimize convolution edge artifacts
                to_conv_num = [fliplr(to_conv_num(1:buff_len)) to_conv_num fliplr(to_conv_num(end-buff_len+1:end))];
                to_conv_den = [fliplr(to_conv_den(1:buff_len)) to_conv_den fliplr(to_conv_den(end-buff_len+1:end))];
                use_vec = logical(ones(size(to_conv_num)));
                use_vec(1:buff_len) = 0; use_vec(end-buff_len+1:end) = 0;
                
                cnum = conv(to_conv_num,win_kern,'same');
                cden = conv(to_conv_den,win_kern,'same');
                meanfun_est_d2{i}(l,:) = cnum(use_vec)./cden(use_vec);
                tot_prob{i}(l,:) = cden(use_vec)/length(win_kern);
                meanfun_est_d2{i}(l,tot_prob{i}(l,:) < p_thresh) = nan; %set any estimates of meanfunction with insufficient posterior probability to nan
            end
            meandiff = [meandiff meanfun_est_d2{i}(2,tot_prob{i}(l,:) >= p_thresh)-...
                meanfun_est_d2{i}(tot_prob{i}(l,:) >= p_thresh)]; %compute a vector of difference of state meanfunctions
        end
        avg_meandiff = nanmean(meandiff); %average difference of mean functions
        %any instances where there was insufficient posterior probability
        %for one state are given by estimates extrapolated from the other
        %state's mean
        for i = 1:hmm.Nsegs
            meanfun_est_d2{i}(1,isnan(meanfun_est_d2{i}(1,:))) = meanfun_est_d2{i}(2,isnan(meanfun_est_d2{i}(1,:))) - avg_meandiff;
            meanfun_est_d2{i}(2,isnan(meanfun_est_d2{i}(2,:))) = meanfun_est_d2{i}(1,isnan(meanfun_est_d2{i}(2,:))) + avg_meandiff;
            meanfun_est_d2{i} = meanfun_est_d2{i}';
        end
        
    else
        for l = 1:K
            o_sum = 0;
            g_sum = 0;
            for i = 1:hmm.Nsegs
                %                o_sum = o_sum + sum(gamma{i}(l,:).*emiss(UDS_seg_inds{i},1)',2);
                %                g_sum = g_sum + sum(gamma{i}(l,:),2);
                o_sum = o_sum + sum(gamma{i}(l,flex_model_inds{i}).*emiss(UDS_seg_inds{i}(flex_model_inds{i}),1)',2);
                g_sum = g_sum + sum(gamma{i}(l,flex_model_inds{i}),2);
            end
            fixedmean_est(l) = o_sum/g_sum;
        end
    end
    
    for l = 1:K
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            if strcmp(hmm.meantype,'variable') 
                state_mean_est = [meanfun_est_d1{i}(:,l) meanfun_est_d2{i}(:,l)];
            else
                state_mean_est = repmat(fixedmean_est(:,l),curT,1);
            end
            mdiff{i} = emiss(UDS_seg_inds{i},:)-state_mean_est;
        end
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            %             o_sum = o_sum + (repmat(gamma{i}(l,:),p,1)'.*mdiff{i})'*mdiff{i};
            %             g_sum = g_sum + sum(gamma{i}(l,:),2);
            o_sum = o_sum + (repmat(gamma{i}(l,flex_model_inds{i}),p,1)'.*mdiff{i}(flex_model_inds{i},:))'*mdiff{i}(flex_model_inds{i},:);
            g_sum = g_sum + sum(gamma{i}(l,flex_model_inds{i}),2);
        end
        var_est(l,:,:) = o_sum/g_sum;
        if var_est(l,2,1) ~= var_est(l,1,2)
            new_val = (var_est(l,2,1)+var_est(l,1,2))/2;
            var_est(l,2,1) = new_val;
            var_est(l,1,2) = new_val;
        end
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
            for i = 1:hmm.Nsegs
                hmm.state(l).meanfun{i}(:,1) = meanfun_est_d1{i}(:,l);
                hmm.state(l).meanfun{i}(:,2) = meanfun_est_d2{i}(:,l);
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
        % hmm.P = P_np;
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
hmm.rob_model_inds = rob_model_inds;