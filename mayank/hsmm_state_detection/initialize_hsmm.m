function [hmm] = initialize_hsmm(emiss,Fs,num_states,windowSize,windowSlide,state_mean_hcf)

addpath('F:\Code\fullBNT-1.0.4\kPMstats\')
addpath('F:\Code\fullBNT-1.0.4\netlab3.3\')
addpath('F:\WC_Germany\hmm_state_detect')

hmm.p=size(emiss,2);
hmm.T=size(emiss,1);
hmm.Fs = Fs;
hmm.K = num_states;
hmm.windowSize = windowSize;
hmm.windowSlide = windowSlide;
hmm.hcf = state_mean_hcf;
time = (1:hmm.T)/Fs;

%% choose random starting values for A and Pi from dir priors
pi_alpha = 1;
Pi_alpha = ones(1,num_states)*pi_alpha;
Pi = dirichletrnd(Pi_alpha);
A_other = 1/Fs/(num_states-1);
A_self = 1-A_other;
for i = 1:num_states
    for j = 1:num_states
        if i==j
            A(i,j)=A_self;
        else
            A(i,j)=A_other;
        end
    end
end

%%
[state_means,state_t,obs_dist,obsrange,uni_times] = ...
    locate_state_means(emiss(:,1),windowSize,windowSlide,Fs);
state_t  = [0 state_t max(time)];
state_means = [state_means(:,1) state_means state_means(:,end)];
lf_meanfuns = interp1(state_t,state_means',time);

%% initialize emissions model params using GMM fit
mix=gmm(hmm.p,hmm.K,'full');
gmm_options(3) = 1e-10; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 100; %max iterations

temp_emiss = emiss;
temp_emiss(:,1) = emiss(:,1)-lf_meanfuns(:,1);
mix_d = gmmem(mix, temp_emiss, gmm_options);
dstate = find(abs(mix_d.centres(:,1))==min(abs(mix_d.centres(:,1)))); %find mixture component with 0 mean
dstate_var = squeeze(mix_d.covars(:,:,dstate));
if hmm.p > 1
dstate_hf = mix_d.centres(dstate,2:end);
end

temp_emiss = emiss;
temp_emiss(:,1) = emiss(:,1)-lf_meanfuns(:,2);
mix_u = gmmem(mix, temp_emiss, gmm_options);
ustate = find(abs(mix_u.centres(:,1))==min(abs(mix_u.centres(:,1)))); %find mixture component with 0 mean
ustate_var = squeeze(mix_u.covars(:,:,ustate));
if hmm.p > 1
    ustate_hf = mix_u.centres(ustate,2:end);
end

%% create hmm structure
hmm.A = A;
hmm.Pi = Pi;

% state_means = mix.centres(:,1);
% [dummy,state_order] = sort(state_means);
% for i = 1:num_states
%     if hmm.p > 1
%         hmm.state(i).hf_mean=mix.centres(state_order(i),2);
%     end
%     hmm.state(i).lf_meanfun = repmat(mix.centres(state_order(i),1),hmm.T,1);
%     hmm.state(i).var=squeeze(mix.covars(:,:,state_order(i)));
% end

if hmm.p > 1
hmm.state(1).hf_mean = dstate_hf;
hmm.state(2).hf_mean = ustate_hf;
end
hmm.state(1).lf_meanfun = lf_meanfuns(:,1);
hmm.state(2).lf_meanfun = lf_meanfuns(:,2);
hmm.state(1).var = dstate_var;
hmm.state(2).var = ustate_var;

hmm.covtype = 'full'; %default
if state_mean_hcf == 0
    hmm.meantype = 'fixed';
else
    hmm.meantype = 'variable';
end

hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'gamma';
   

