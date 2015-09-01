function [hmm] = initialize_hsmm_fixmeans(emiss,Fs,num_states)

addpath('C:\Code\fullBNT-1.0.4\kPMstats\')
addpath('C:\Code\fullBNT-1.0.4\netlab3.3\')
addpath('G:\WC_Germany\hmm_state_detect')

hmm.p=size(emiss,2);
hmm.T=size(emiss,1);
hmm.Fs = Fs;
hmm.K = num_states;
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

%% initialize emissions model params using GMM fit
mix=gmm(hmm.p,hmm.K,'full');
gmm_options(3) = 1e-10; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 100; %max iterations

mix = gmmem(mix, emiss, gmm_options);

%% create hmm structure
hmm.A = A;
hmm.Pi = Pi;

state_means = mix.centres(:,1);
[dummy,state_order] = sort(state_means);
for i = 1:num_states
    hmm.state(i).mean=mix.centres(state_order(i),:);
    hmm.state(i).var=squeeze(mix.covars(:,:,state_order(i)));
end

hmm.covtype = 'full'; %default

