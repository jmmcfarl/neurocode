%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% cd G034/

% load ./Expt1_compiled_data_ms.mat
% load ./Expt1_compiled_data_yinv.mat
load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt1_compiled_data.mat

fullX = bsxfun(@minus,fullX,mean(fullX));
fullX = bsxfun(@rdivide,fullX,std(fullX));

Pix2Deg = 0.018837;
dsfrac = 1.5;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [0 0.8; -0.8 -0.]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

%%
sdim = length(xpatch_inds);
bandwidth = 0.3;
ar = 1.25;
kern_len = 21-1;
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
orientations = linspace(0,pi,9);
orientations(end) = [];
% orientations = 0.9;
lambdas = [1/6];
% gab_xo = [0.35];
% gab_yo = [-0.48];
gab_xo = [0:0.1:0.8];
gab_yo = [-0.8:0.1:-0];
gabor_bank = zeros(length(orientations)*length(lambdas)*length(gab_xo)*length(gab_yo),size(XX,1),size(XX,2));
cnt = 1;
clear gabor_props
for i = 1:length(orientations)
    for j = 1:length(lambdas)
        for m = 1:length(gab_yo)
            for l = 1:length(gab_xo)
                gabor_bank(cnt,:,:) = get_gabor_template(XX,YY,gab_xo(l),gab_yo(m),orientations(i),lambdas(j),0,bandwidth,ar);
                gabor_props(cnt).xo = gab_xo(l);
                gabor_props(cnt).yo = gab_yo(m);
                gabor_props(cnt).lambda = lambdas(j);
                gabor_props(cnt).ori = orientations(i);
                gabor_props(cnt).bandwidth = bandwidth;
                gabor_props(cnt).ar = ar;
                
                gabor_bank2(cnt,:,:) = get_gabor_template(XX,YY,gab_xo(l),gab_yo(m),orientations(i),lambdas(j),pi/2,bandwidth,ar);
                cnt = cnt + 1;
            end
        end
    end
end
gabor_space_kerns = reshape(gabor_bank,size(gabor_bank,1),sdim^2);
gabor_space_kerns2 = reshape(gabor_bank2,size(gabor_bank,1),sdim^2);
temp_kern = fliplr(gampdf(1:flen,4,1));

gabor_filts = repmat(gabor_space_kerns,[1 1 flen]);
gabor_filts = bsxfun(@times,gabor_filts,reshape(temp_kern,[1 1 flen]));
gabor_filts = reshape(permute(gabor_filts,[1 3 2]),size(gabor_filts,1),flen*sdim^2);

gabor_filts2 = repmat(gabor_space_kerns2,[1 1 flen]);
gabor_filts2 = bsxfun(@times,gabor_filts2,reshape(temp_kern,[1 1 flen]));
gabor_filts2 = reshape(permute(gabor_filts2,[1 3 2]),size(gabor_filts,1),flen*sdim^2);

%%
clear gabor_emp*

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
% load ./onednoise_rffits
% load ./oned_noise_mapping
load ./oned_fixation_fits

for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    %     init_params(1) = rand*0.6+0.1; %x0
    %     init_params(2) = rand*0.6-0.8; %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    %     init_params(3) = degtorad(pref_mu_ori(t))+rand*2*pi; %theta
    init_params(4) = 1/6; %lambda
    for i = 1:flen
        gabor_emp1 = bsxfun(@times,gabor_emp1,temp_kern');
        gabor_emp2 = bsxfun(@times,gabor_emp2,temp_kern');
        
        gabor_emp1_filt(t,:) = gabor_emp1(:);
        gabor_emp2_filt(t,:) = gabor_emp2(:);
    end
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    
    for j = 1:flen
        gabor_emp1 = repmat(gabor_emp1(:),[1 flen])';
        gabor_emp2 = repmat(gabor_emp2(:),[1 flen])';
        
        cur_tkern = zeros(flen,1);
        cur_tkern(j) = 1;
        gabor_emp1 = bsxfun(@times,gabor_emp1,cur_tkern);
        gabor_emp2 = bsxfun(@times,gabor_emp2,cur_tkern);
        
        gabor_emp1_filt(t,j,:) = gabor_emp1(:);
        gabor_emp2_filt(t,j,:) = gabor_emp2(:);
    end
end

%%
% gabor_out_X1 = fullX*gabor_filts';
% gabor_out_X2 = fullX*gabor_filts2';
% gabor_out_X = sqrt(gabor_out_X1.^2 + gabor_out_X2.^2);
% gabor_out_X = zscore(gabor_out_X);
% n_gabor_bases = size(gabor_out_X,2);

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
gabor_emp_X = zscore(gabor_emp_X);


oris = [gabor_props(:).ori];
XO = [gabor_props(:).xo];
YO = [gabor_props(:).yo];


%%
xv_frac = 0.2;
NT = size(fullX,1);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

%%

clear gabor_sta
for t = 96
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    unique_spk_cnts = unique(cur_binned_spks);
    spikebins = [];
    for i = 2:length(unique_spk_cnts)        
        cur_set = find(cur_binned_spks == unique_spk_cnts(i));
        spikebins = [spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
    end
    spikebins = sort(spikebins);
        
    [beta(t,:),dev,stats] = glmfit(gabor_out_X(tr_samp,:),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    
    for or = 1:length(orientations)
        cur_set = find(oris==orientations(or));
        cur_b_im(or,:,:) = reshape(beta(t,cur_set),length(gab_yo),length(gab_xo));
    end
    
    
    pred_r = exp(beta(t,2)*gabor_emp_X(tr_samp,t)+beta(t,1));
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = exp(beta(t,2)*gabor_emp_X(xv_samp,t)+beta(t,1));
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
end

%%
