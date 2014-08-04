function hmm = hsmm_uds_train_hmm(hmm,emiss)

%Adapted from VARHMMBOX, version 1.1, Iead Rezek, Oxford University, MAR
%2002
% train an HMM
%
% Input:
%       hmm: structure with the initialized parameters of an hmm
%       emiss: observation sequence
%
% Output:
%       hmm: structure of parameters with the following fields
%           p: dimension of emissions
%           T: number of time samples
%           Fs: sample frequency of emissions
%           K: number of hidden states
%           UDS_segs: Nx2 matrix containing the time stamps of the start
%             and stop of each UDS segment
%           Nsegs: number of UDS segments
%           uds_dur: total duration of UDS segments (in s)
%           windowSize: duration of windows for estimating
%             state-conditional means
%           windowSlide: resolution of sliding-window temporal sampling
%           A: KxK transition probability matrix
%           Pi: NxK matrix of initial state probabilities for each of N
%             segments
%           meantype: state-conditional means are either 'variable' or 'fixed'
%           min_mllik: minimum value of the maximal log likelihood for
%             either state before switching to robust estimator
%           hmm_LP: Ix1 vector of the log likelihood at each iteration
%           gamma: cell array containing the sequence of marginal probabilities of each
%              state for each data segment
%           state: 1xK structure array containing the following fields
%                fixedmean: (for fixed state means): state conditional mean
%                meanfun: (for variable state means): state conditional mean function
%                var: state-conditional covariance matrix


tolerance = 1e-5; %convergence tolerance for the likelihood
cur_log_post = -Inf; %will store the current iterations likelihood
delta_log_post = Inf; %will store the current change in likelihood from the previous iteration

%%
K = hmm.K; %number of hidden states
p = hmm.p; %dimensionality of observations
T = hmm.T; %number of time samples
Fs = hmm.Fs; %sampling frequency
if strcmp(hmm.meantype,'variable')
    windowSize = hmm.windowSize;
%     windowSlide = hmm.windowSlide;
    p_thresh = 0.05; %minimum probability of occupying a given state in order to
%     use direct estimation of that state's mean function within a given
%     time window. Default 0.05
end

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

LP = []; %this will store the log likelihoods of the observations given model parameters across iterations of EM

%%

it = 0; %iteration number
while delta_log_post > tolerance
    
    %% compute posterior over hidden states
    %compute likelihood of observations for each model under current
    %parameters
    if ~exist('gamma','var')
        rel_prob = [0.5 0.5];
    else
        rel_prob = [0 0];
        for i = 1:hmm.Nsegs
            rel_prob = rel_prob + sum(gamma{i});
        end
        rel_prob = rel_prob/sum(rel_prob);
    end
    
    %loop over UDS segments
    for i = 1:hmm.Nsegs
        [B,rob_model_inds{i}] = mvgauss_uds_obslike_robust(emiss(UDS_seg_inds{i},:),hmm,rel_prob,i);
        flex_model_inds{i} = setdiff(1:length(UDS_seg_inds{i}),rob_model_inds{i});
        
        %initialize forward and backward messages
        curT = length(UDS_seg_inds{i});
        alpha{i}=zeros(curT,K);
        beta{i}=zeros(curT,K);
        gamma{i}=zeros(curT,K);
        
        scale=zeros(curT,1); %initialize rescaling parameters
        
        %compute rescaled forward messages
        alpha{i}(1,:)=Pi(i,:).*B(1,:);
        scale(1)=sum(alpha{i}(1,:));
        alpha{i}(1,:)=alpha{i}(1,:)/scale(1);
        for t=2:curT
            alpha{i}(t,:)=(alpha{i}(t-1,:)*A).*B(t,:);
            scale(t)=sum(alpha{i}(t,:));
            alpha{i}(t,:)=alpha{i}(t,:)/scale(t);
        end
        
        %compute rescaled backward messages
        beta{i}(curT,:)=ones(1,K)/scale(curT);
        for t=curT-1:-1:1
            beta{i}(t,:)=(beta{i}(t+1,:).*B(t+1,:))*(A')/scale(t);
        end
        
        %compute posteriors over hidden states
        gamma{i}=(alpha{i}.*beta{i});
        gamma{i} = gamma{i}./repmat(sum(gamma{i},2),1,K);
        
        %compute chi (posterior of two consecutive hidden states)
        chi{i}=zeros(curT-1,K*K);
        for t=1:curT-1
            temp=A.*(alpha{i}(t,:)' * (beta{i}(t+1,:).*B(t+1,:)));
            chi{i}(t,:)=temp(:)'/sum(temp(:));
        end
        
        lscale = log(scale);
        %         loglik(i) = sum(lscale); %rabiner eq(103), scale is defined as inverse here though
        loglik(i) = sum(lscale(flex_model_inds{i})); %rabiner eq(103), scale is defined as inverse here though
    end
    
    loglik = sum(loglik); %likelihood over all UDS segments
    LP=[LP; loglik];
    delta_log_post = (loglik - cur_log_post)/abs(loglik);
    cur_log_post = loglik;
    
    %% M STEP
    % transition matrix
    expchi = zeros(1,K^2);
    for i = 1:hmm.Nsegs
        %         expchi = expchi + sum(chi{i},1);
        expchi = expchi + sum(chi{i}(flex_model_inds{i}(1:end-1),:),1);
    end
    expchi=reshape(expchi,K,K);
    A = expchi./repmat(sum(expchi,2),1,K);
    
    % initial state
    for i = 1:hmm.Nsegs
        Pi(i,:)=gamma{i}(1,:)/sum(gamma{i}(1,:));
    end
    
    var_est = zeros(K,p,p);
    if strcmp(hmm.meantype,'fixed')
        fixedmean_est = zeros(K,p);
    elseif p > 1
        fixedmean_est = zeros(K,p-1); %if there is a nonstationary dimension the number of fixed dimensions is p-1
    end
    if strcmp(hmm.meantype,'variable')
        meandiff = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            %initialize
            cur_tot_prob = nan(hmm.K,length(UDS_seg_inds{i}));
            meanfun_est{i} = zeros(hmm.K,length(UDS_seg_inds{i}));
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
                
                to_conv_num = gamma{i}(:,l).*emiss(UDS_seg_inds{i});
                to_conv_num = interp1(to_use_inds,to_conv_num(to_use_inds),1:length(to_conv_num));
                to_conv_den = gamma{i}(:,l);
                to_conv_den = interp1(to_use_inds,to_conv_den(to_use_inds),1:length(to_conv_den));
                
                %add a buffer window onto either side of each vector to
                %minimize convolution edge artifacts
                to_conv_num = [fliplr(to_conv_num(1:buff_len)) to_conv_num fliplr(to_conv_num(end-buff_len+1:end))];
                to_conv_den = [fliplr(to_conv_den(1:buff_len)) to_conv_den fliplr(to_conv_den(end-buff_len+1:end))];
                use_vec = logical(ones(size(to_conv_num)));
                use_vec(1:buff_len) = 0; use_vec(end-buff_len+1:end) = 0;
                
                cnum = conv(to_conv_num,win_kern,'same');
                cden = conv(to_conv_den,win_kern,'same');              
                meanfun_est{i}(l,:) = cnum(use_vec)./cden(use_vec);
                cur_tot_prob(l,:) = cden(use_vec)/length(win_kern);
                meanfun_est{i}(l,cur_tot_prob(l,:) < p_thresh) = nan; %set any estimates of meanfunction with insufficient posterior probability to nan
            end
            meandiff = [meandiff meanfun_est{i}(2,cur_tot_prob(l,:) >= p_thresh)-...
                meanfun_est{i}(1,cur_tot_prob(l,:) >= p_thresh)]; %compute a vector of difference of state meanfunctions
        end
        avg_meandiff = nanmean(meandiff); %average difference of mean functions
        
        %any instances where there was insufficient posterior probability
        %for one state are given by estimates extrapolated from the other
        %state's mean
        for i = 1:hmm.Nsegs
            meanfun_est{i}(1,isnan(meanfun_est{i}(1,:))) = meanfun_est{i}(2,isnan(meanfun_est{i}(1,:))) - avg_meandiff;
            meanfun_est{i}(2,isnan(meanfun_est{i}(2,:))) = meanfun_est{i}(1,isnan(meanfun_est{i}(2,:))) + avg_meandiff;
            meanfun_est{i} = meanfun_est{i}';
        end
     else %for fixed state means
        for l = 1:K
            o_sum = 0;
            g_sum = 0;
            for i = 1:hmm.Nsegs
                %                o_sum = o_sum + sum(repmat(gamma{i}(:,l),1,p).*emiss(UDS_seg_inds{i},:));
                %                g_sum = g_sum + sum(gamma{i}(:,l));
                o_sum = o_sum + sum(repmat(gamma{i}(flex_model_inds{i},l),1,p).*emiss(UDS_seg_inds{i}(flex_model_inds{i}),:));
                g_sum = g_sum + sum(gamma{i}(flex_model_inds{i},l));
            end
            fixedmean_est(l) = o_sum/g_sum;
        end
    end
    
    for l = 1:K %restimate state covariance matrices
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            if strcmp(hmm.meantype,'variable')
                state_mean_est = meanfun_est{i}(:,l);
            else
                state_mean_est = repmat(fixedmean_est(l,:),curT,1);
            end
            mdiff{i} = emiss(UDS_seg_inds{i},:)-state_mean_est;
        end
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            %             o_sum = o_sum + (repmat(gamma{i}(:,l),1,p).*mdiff{i})'*mdiff{i};
            %             g_sum = g_sum + sum(gamma{i}(:,l));
            o_sum = o_sum + (repmat(gamma{i}(flex_model_inds{i},l),1,p).*mdiff{i}(flex_model_inds{i},:))'*mdiff{i}(flex_model_inds{i},:);
            g_sum = g_sum + sum(gamma{i}(flex_model_inds{i},l));
        end
        var_est(l,:,:) = o_sum/g_sum;
    end
    
    %update all parameters in hmm
    hmm.A = A;
    hmm.Pi = Pi;
    for l = 1:K
        new_covar = squeeze(var_est(l,:,:));
        if det(new_covar) > 0
            hmm.state(l).var = squeeze(var_est(l,:,:));
        else
            disp('warning: singular covariance')
        end
        if p > 1 || strcmp(hmm.meantype,'fixed')
            hmm.state(l).fixedmean = fixedmean_est(l,:);
        end
        if strcmp(hmm.meantype,'variable')
            for i = 1:hmm.Nsegs
                hmm.state(l).meanfun{i} = meanfun_est{i}(:,l);
            end
        end
    end
    it = it + 1;
    fprintf('iteration %i log posterior = %f \n',it,loglik);
    %         jmm_print_hmm_params(hmm);
    %     pause
    
end

%store the observation likelihood under the ML parameters
hmm.hmm_LP = LP;
hmm.gamma = gamma;
hmm.rob_model_inds = rob_model_inds;