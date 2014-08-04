%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/G075/
load ./CellList.mat
% single_units = find(CellList(1,:,1) > 0);

load ./Expt12_compiled_data_fixeddelay_d1p25.mat
% load ./NS2_compiled_data_fixeddelay_d1p25.mat
% load ./GR_compiled_data_fixeddelay_d1p25.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);
Nyp = 1024;
Nxp = 1280;


%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));


dt = dt*2;

% SET UP XV CELL SET OF CELLS
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%%
load ./grating_unit_tuning
pref_freq = pref_freq/(Nxp*Pix2Deg);


load ./NS40_gabor_mods
%%
xv_frac = 0.2;
% xv_frac =0;
NT = size(fullX,1);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
% 
% % [~,~,full_expt_inds] = unique(full_expt_vec); 

gab_priors(1).type = 'gauss';
gab_priors(1).theta(2) = 0.5; %prior std

gab_priors(2).type = 'gauss';
gab_priors(2).theta(2) = 0.5;

gab_priors(4).type = 'gam';
gab_priors(4).theta(1) = 8; %shape
gab_priors(4).theta(2) = 0.03; %scale

gab_priors(5).type = 'gam';
gab_priors(5).theta(1) = 8; %shape
gab_priors(5).theta(2) = 0.011; %scale

gab_priors(6).type = 'gam';
gab_priors(6).theta(1) = 2; %shape
gab_priors(6).theta(2) = 2; %scale

mean_x = 0.45;
mean_y = -0.45;
LB = [-0.1 -0.8 0 0.125 0.025 0.2 0 -Inf];
UB = [0.8 0.1 pi 0.4 0.25 6 Inf Inf];
% hold_const = [0 0 1 1 0 1 0 0];
hold_const = [1 1 1 1 1 1 0 0];

for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
%     gab_priors(1).theta(1) = mean_x;
%     gab_priors(2).theta(1) = mean_y;

%     init_params(1) = mean_x; %x0
%     init_params(2) = mean_y; %y0
%     init_params(3) = degtorad(pref_ori(t)); %theta
%     init_params(4) = 1/pref_freq(t); %lambda
%     init_params(5) = 0.5*init_params(4); %sigma
%     init_params(6) = 1.5; %eccentricity of gaussian ellipse
%     init_params(7) = 0; %weight of quad term
%     init_params(8) = 0; %const offset
    
    init_params = gabor_params_f(t,:);
    
%     [gabor_params_f(t,:),LL(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped(tr_samp,:),...
%         full_binned_spks(tr_samp,t),hold_const,LB,UB,gab_priors);
    [gabor_params_f(t,:),LL(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped(tr_samp,:),...
        full_binned_spks(tr_samp,t),hold_const,LB,UB);
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
    
    pout1 = fullX_cropped(xv_samp,:)*gabor_emp1_filt(t,:)';
    pout2 = fullX_cropped(xv_samp,:)*gabor_emp2_filt(t,:)';
    energy_out = gabor_params_f(t,7)*sqrt(pout1.^2+pout2.^2);
    g = energy_out + gabor_params_f(t,8);
    too_large = find(g > 100);
    r = log(1+exp(g));
    r(too_large) = g(too_large);
    gabor_LL(t) = -sum(full_binned_spks(xv_samp,t).*log(r) - r)/sum(full_binned_spks(xv_samp,t));

    
    avg_rate = mean(full_binned_spks(tr_samp,t));
    null_prate = avg_rate*ones(length(xv_samp),1);
    null_xvLL(t) = -sum(full_binned_spks(xv_samp,t).*log(null_prate) - null_prate)/sum(full_binned_spks(xv_samp,t));
    null_prate = avg_rate*ones(length(tr_samp),1);
    null_LL(t) = -sum(full_binned_spks(tr_samp,t).*log(null_prate) - null_prate)/sum(full_binned_spks(tr_samp,t));
    
    
end

%%
% cd ~/Data/bruce/G075/
% save NS40_gabor_mods gabor_params_f init_params gab_priors *_xvLL *LL 
%%
load ./ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
%%
% mean_x = 0.45;
% mean_y = -0.45;
% 
% for t = 1:96
%     t
%     init_params(1) = mean_x; %x0
%     init_params(2) = mean_y; %y0
%     init_params(3) = degtorad(pref_ori(t)); %theta
%     init_params(4) = 1/pref_freq(t); %lambda
% %     init_params(4) = 1/6; %lambda
%     init_params(5) = 0.5*init_params(4); %sigma
%     init_params(6) = 1; %eccentricity of gaussian ellipse
% 
%     gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,init_params,0);
%     gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,init_params,pi/2);
%     
%     gabor_emp1_filt(t,:) = gabor_emp1(:);
%     gabor_emp2_filt(t,:) = gabor_emp2(:);
%     
%     pout1 = fullX_cropped(tr_samp,:)*gabor_emp1_filt(t,:)';
%     pout2 = fullX_cropped(tr_samp,:)*gabor_emp2_filt(t,:)';
%     energy_out = sqrt(pout1.^2+pout2.^2);
% 
%     [B,DEV,STATS] = glmfit(energy_out,full_binned_spks(tr_samp,t),'poisson');
%     
%     weight(t) = B(2);
%     weight_se(t) = STATS.se(2);
%     
%     pout1 = fullX_cropped(xv_samp,:)*gabor_emp1_filt(t,:)';
%     pout2 = fullX_cropped(xv_samp,:)*gabor_emp2_filt(t,:)';
%     energy_out = sqrt(pout1.^2+pout2.^2);
%     pred_r = glmval(B,energy_out,'log');
%         gabor_LL(t) = -sum(full_binned_spks(xv_samp,t).*log(pred_r) - pred_r)/sum(full_binned_spks(xv_samp,t));
% 
%     avg_rate = mean(full_binned_spks(tr_samp,t));
%     null_prate = avg_rate*ones(length(xv_samp),1);
%     null_xvLL(t) = -sum(full_binned_spks(xv_samp,t).*log(null_prate) - null_prate)/sum(full_binned_spks(xv_samp,t));
% 
% end