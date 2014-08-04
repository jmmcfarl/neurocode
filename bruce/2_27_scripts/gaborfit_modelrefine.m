clear all
close all

Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_stim_righteye
% load fixation_stim_lefteye

load right_eye_gabortrack_v2
% load left_eye_gabortrack_sm

addpath(genpath('~/James_scripts'));
%%
SDIM = size(X,2);
Xmat_resh = reshape(X,size(X,1),SDIM,SDIM);
for i = 1:size(X,1)
    Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-cur_X_seq(end,i),2,0);
    Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-cur_Y_seq(end,i),1,0);
end

%%
SDIM = 33;
kern_l = SDIM^2;
nmods = 2;
% init_kerns = sta(:);
% init_signs = 1;
init_kerns = randn(kern_l,nmods);
init_signs = ones(nmods,1);

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
defmod.fsdim = SDIM^2;
defmod.pids = 1:defmod.fsdim;
defmod.h = 1;
defmod.SDIM = SDIM;
defmod.locLambda = 0;
defmod.lambda_dX = 1e6;
defmod.lambda_dT = 0;
defmod.lambda_L1x = 1e3;
basis = 'pix';

for cellid = 1:10
    Xmat = reshape(X,[size(X,1) size(X,2)*size(X,3)]);
    fprintf('Fitting cell %d of %d\n',cellid,10);
    glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
    glm0.image_type = '2d';
    
    glm0 = normalizeRFs_full(glm0,Xmat);
    ml_sta_orig(cellid) = fitGLM_lexp(glm0,Xmat,spikebins{cellid},'tots',[],1e-5,1e-6);
    
    Xmat = reshape(Xmat_resh,[size(Xmat_resh,1) size(Xmat_resh,2)*size(Xmat_resh,3)]);
    glm0 = normalizeRFs_full(glm0,Xmat);
    ml_sta_ref(cellid) = fitGLM_lexp(glm0,Xmat,spikebins{cellid},'tots',[],1e-5,1e-6);
end

%%
kern_len = 26;
[XX,YY] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);
orientations = linspace(0,pi,10);
for i = 1:length(orientations)
gabor1 = get_gabor_template(XX,YY,0,0,orientations(i),8.85,0);
nzpad = 30;
zmat = zeros(SDIM+2*nzpad,SDIM+2*nzpad);
cur_mask = zmat;
cur_mask(nzpad+(-13:13)+25,nzpad+(-13:13)+15) = gabor1;
cur_mask = cur_mask(nzpad+(1:SDIM),nzpad+(1:SDIM));

Xmat = reshape(Xmat_resh,[size(X,1) size(X,2)*size(X,3)]);
cellid = 10;
fprintf('Fitting cell %d of %d\n',cellid,10);
glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
glm0.image_type = '2d';
glm0.mods(1) = [];
glm0.mods(1).k = cur_mask(:);
glm0.mods(1).pix = cur_mask(:);
glm0 = normalizeRFs_full(glm0,Xmat);
glm0.const = -0.08;
NT = size(Xmat,1);
[ll0(i),llp,prate] = getLLGLM_lexp(glm0,Xmat,spikebins{cellid},'none');
% Robs = zeros(1,NT);
% ftable = tabulate(spikebins{cellid});
% Robs(ftable(:,1)) = ftable(:,2);
% 
% ml_sta_new(cellid) = fitGLM_lexp(glm0,Xmat,spikebins{cellid},'tots',100,1e-5,1e-8);
% [ll0,llp,prate2] = getLLGLM_lexp(ml_sta_new(cellid),Xmat,spikebins{cellid},'none');
end

    %%
figure

xax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
yax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;

% xl = [-2 2];
% yl = [-2 2];
xl = [-1.2 1.2];
yl = [-1.2 1.2];

for cellid = 1:5
    k_vec = get_k_mat(ml_sta_orig(cellid));
    k_mat = reshape(k_vec,SDIM,SDIM);
    
    k_vec = get_k_mat(ml_sta_new(cellid));
    k_matr = reshape(k_vec,SDIM,SDIM);
    
    minz = min(min(k_mat(:)),min(k_matr(:)));
    maxz = max(max(k_mat(:)),max(k_matr(:)));
    
    subplot(5,4,(cellid-1)*4+1)
    imagesc(xax,yax,k_mat); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
    subplot(5,4,(cellid-1)*4+2)
    imagesc(xax,yax,k_matr); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
end

for cellid = 6:10
    k_vec = get_k_mat(ml_sta_orig(cellid));
    k_mat = reshape(k_vec,SDIM,SDIM);
    
    k_vec = get_k_mat(ml_sta_new(cellid));
    k_matr = reshape(k_vec,SDIM,SDIM);
    
    minz = min(min(k_mat(:)),min(k_matr(:)));
    maxz = max(max(k_mat(:)),max(k_matr(:)));
    
    subplot(5,4,(cellid-6)*4+3)
    imagesc(xax,yax,k_mat); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
    subplot(5,4,(cellid-6)*4+4)
    imagesc(xax,yax,k_matr); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
end
