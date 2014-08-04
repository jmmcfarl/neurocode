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

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_stim_righteye

%% COMPUTE GABOR TEMPLATES
Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

kern_len = 32;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

cellids = [1 2 3 4 5 6 7 8 10];
orientations = [0 0.31 0.31 0 0 1.9 pi/2 1.9 0];
lambdas = [1.0 0.35 0.35 0.5 0.67 0.35 0.67 0.35 0.67]*Fsd;
phases = [0 pi/2];
vert_loc = [12 20 21 19 24 22 24 22 25];
hor_loc = [14 17 14 17 16 16 19 15 15];
vert_pos = Yg(vert_loc,1);
hor_pos = Xg(1,hor_loc);
n_used_cells = length(cellids);

%% CONVOLVE GABOR TEMPLATES WITH IMAGE PATCH
NT = size(X,1);
SDIM = size(X,2);
Xmat = reshape(X,NT,SDIM^2);

for n = 1:n_used_cells;
    init_params(1) = hor_pos(n); %x0
    init_params(2) = vert_pos(n); %y0
    init_params(3) = orientations(n); %theta
    init_params(4) = lambdas(n); %lambda
    init_params(5) = 0.5*init_params(4); %sigma
    init_params(6) = 1/var(Xmat(:)); %weight of quad term
    init_params(7) = 0; %weight of lin term
    init_params(8) = 0; %preferred phase of lin term
    init_params(9) = 0; %const offset
    
    % % hold_const = [1 1 0 0 0 0 0];
    hold_const = [0 0 0 0 0 0 1 1 0];
    
    [gabor_params,LL] = fit_gabor_params(init_params,Xmat,spk_cnts(:,cellids(n)),[SDIM SDIM],hold_const);
    
    hold_const = [0 0 0 0 0 0 0 0 0];
    gabor_params(7) = 1/var(Xmat(:));
    gabor_params(8) = 0;
    [gabor_params2(n,:),LL(n)] = fit_gabor_params(gabor_params,Xmat,spk_cnts(:,cellids(n)),[SDIM SDIM],hold_const);
end