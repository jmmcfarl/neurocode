clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-7 7;-7 7];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));


addpath(genpath('~/James_scripts'));
%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
% load fixation_stim_righteye
load fixation_data
% load fixation_image_patches
load fixation_image_patches
% X = X_right;
X = X_left;

%% COMPUTE GABOR TEMPLATES
NT = size(X,1);
SDIM = size(X,2);
kern_len = SDIM-1;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

Xmat = reshape(X,NT,SDIM^2);

cellids = [1 2 3 4 5 6 7 8 10];
hor_posd = [4.5 4.72 4.5 4.65 4.65 4.65 4.88 4.57 4.57];
vert_posd = [-2.61 -2 -1.93 -2.08 -1.7 -1.85 -1.7 -1.85 -1.63];
orientations = [0 0.31 0.31 0 0 1.9 pi/2 1.9 0];
lambdas = [1.0 0.35 0.35 0.5 0.67 0.35 0.67 0.35 0.67]*Fsd;
phases = [0 pi/2];
% vert_loc = [12 20 21 19 24 22 24 22 25];
% hor_loc = [14 17 14 16 16 16 19 15 15];
for i = 1:length(cellids)
    [~,vert_ind(i)] = min(abs(yax(ypatch_inds) - vert_posd(i)));
    [~,hor_ind(i)] = min(abs(xax(xpatch_inds) - hor_posd(i)));
end
vert_pos = Yg(vert_ind,1);
hor_pos = Xg(1,hor_ind);
n_used_cells = length(cellids);

for n = 1:n_used_cells;
    fprintf('Fitting cell %d of %d\n',n,n_used_cells);
    init_params(1) = hor_pos(n); %x0
    init_params(2) = vert_pos(n); %y0
    init_params(3) = orientations(n); %theta
    init_params(4) = lambdas(n); %lambda
    init_params(5) = 0.5*init_params(4); %sigma
    init_params(6) = 0; %weight of quad term
    init_params(7) = 0; %weight of lin term 0 phase
    init_params(8) = 0; %weight of lin term pi/2 phase
    init_params(9) = 0; %const offset
    ip(n,:) = init_params(1:5);
    hold_const = [0 0 0 0 0 0 0 0 0];
    [gabor_params(n,:),LL{1}(n)] = fit_gabor_params(init_params,Xmat,spk_cnts(used_fixs,cellids(n)),[SDIM SDIM],hold_const);
end

%%

for blockid = 1:4
    fprintf('Block %d of %d\n',blockid,4);
    cur_data = find(blockids(used_fixs) == blockid);
    
    for n = 1:n_used_cells
        fprintf('Fitting cell %d of %d\n',n,n_used_cells);      
        
        [block_gabor_params(n,blockid,:),block_LL(n,blockid)] = fit_gabor_params(gabor_params(n,:),Xmat(cur_data,:),spk_cnts(used_fixs(cur_data),cellids(n)),[SDIM SDIM],hold_const);
              
    end
end


%%
NT = size(X,1);
SDIM = size(X,2);
kern_len = SDIM-1;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

Xmat = reshape(X,NT,SDIM^2);

muaids = [1 2 4 5 6 9 10 11 12 13 14];
hor_posd = [4.5 4.57 4.65 4.65 2.91 3.51 4.72 4.72 3.0 3.0 4.88];
vert_posd = [-2.0 -1.78 -2.0 -2.0 -1.85 -1.78 -1.78 -1.85 -1.78 -1.78 -1.4];
orientations = [0.39 0.39 0.39 0 1.57 1.57 1.57 0.79 1.57 1.57 0];
lambdas = [0.5 0.5 0.36 0.36 1 0.67 0.36 0.67 0.67 0.67 1]*Fsd;
phases = [0 pi/2];
for i = 1:length(muaids)
   [~,vert_ind(i)] = min(abs(yax(ypatch_inds) - vert_posd(i))); 
   [~,hor_ind(i)] = min(abs(xax(xpatch_inds) - hor_posd(i))); 
end
vert_pos = Yg(vert_ind,1);
hor_pos = Xg(1,hor_ind);
n_used_mus = length(muaids);
gabor_bank_mua = zeros(n_used_mus,length(phases),size(Xg,1),size(Xg,2));

for n = 1:n_used_mus;
    fprintf('Fitting cell %d of %d\n',n,n_used_mus);
    init_params(1) = hor_pos(n); %x0
    init_params(2) = vert_pos(n); %y0
    init_params(3) = orientations(n); %theta
    init_params(4) = lambdas(n); %lambda
    init_params(5) = 0.5*init_params(4); %sigma
    init_params(6) = 0; %weight of quad term
    init_params(7) = 0; %weight of lin term 0 phase
    init_params(8) = 0; %weight of lin term pi/2 phase
    init_params(9) = 0; %const offset
    ip(n,:) = init_params(1:5);
    hold_const = [0 0 0 0 0 0 0 0 0];
    [mgabor_params(n,:),mLL(n)] = fit_gabor_params(init_params,Xmat,mua_cnts(used_fixs,muaids(n)),[SDIM SDIM],hold_const);
end

%%

for blockid = 1:4
    fprintf('Block %d of %d\n',blockid,4);
    cur_data = find(blockids(used_fixs) == blockid);
    
    for n = 1:n_used_mus
        fprintf('Fitting cell %d of %d\n',n,n_used_mus);      
        
        [mblock_gabor_params(n,blockid,:),mblock_LL(n,blockid)] = fit_gabor_params(mgabor_params(n,:),Xmat(cur_data,:),mua_cnts(used_fixs(cur_data),muaids(n)),[SDIM SDIM],hold_const);
              
    end
end




