% reconstruct the stimulus
% created by Yuwei Cui, Feb 28, 2012
% edited by James McFarland

clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/

%%
load Blocks.mat

% original image resolution
% Pix2Deg = 1.1279 / 60;
Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 4;

% block number
blockid = 1;

Ny = 1024;
Nx = 1280;

Nyd = Ny/dsfrac;
Nxd = Nx/dsfrac;

stimtime = Blocks{blockid}.stimtime;
stimID = Blocks{blockid}.stimids;
Nimage = length(stimID);

gaussfilt = fspecial('gaussian',5,0.25);

all_IMage_array = [];

Fsd = 1/Pix2Deg/dsfrac; %new sampling frequency (pix/deg)
sigx = 0.75*Fsd; %sigma for Gaussian
sigy = 0.75*Fsd;

win_x = round(-3*sigx):round(3*sigx);
win_y = round(-3*sigy):round(3*sigy);
[X,Y] = meshgrid(win_x,win_y);
gaussfun = exp(-(X.^2/(2*sigx^2)+Y.^2/(2*sigy^2)));
gaussfun = gaussfun/sum(gaussfun(:));

edge_sig = 0.25*Fsd;
edgeGauss = fspecial('gaussian',round(10*edge_sig),edge_sig);

eps = 1e-10;
bandpass_freqs = [0.3 1;1 2.5];
[freqx,freqy] = meshgrid(linspace(-Fsd/2,Fsd/2,length(win_x)),linspace(-Fsd/2,Fsd/2,length(win_x)));
freqnorm = sqrt(freqx.^2+freqy.^2);
for i = 1:size(bandpass_freqs,1)
    hor_freqrange{i} = find(abs(freqx) >= bandpass_freqs(i,1) & abs(freqx) <= bandpass_freqs(i,2) & freqy == 0);
    ver_freqrange{i} = find(abs(freqy) >= bandpass_freqs(i,1) & abs(freqy) <= bandpass_freqs(i,2) & freqx == 0);
    diagur_freqrange{i} = find(abs(freqy - freqx) < eps & freqnorm >= bandpass_freqs(i,1) & freqnorm <= bandpass_freqs(i,2));
    diagul_freqrange{i} = find(abs(freqy + freqx) < eps & freqnorm >= bandpass_freqs(i,1) & freqnorm <= bandpass_freqs(i,2));
end
freq = linspace(-Fsd/2,Fsd/2,length(win_x));

grid_dx = sigx/2;
grid_dy = sigy/2;
gridX = round(grid_dx/2:grid_dx:(Nxd-grid_dx/2));
gridY = round(grid_dy/2:grid_dy:(Nyd-grid_dy/2));

%%
nzpad = 100;

all_Energy_array = [];

% load images
for i=1:Nimage
    
    fprintf('Image %d of %d\n',i,Nimage);
    if stimID(i)<10
        filename = ['ExptA0000' num2str(stimID(i)) '.png'];
    elseif stimID(i)>=10 && stimID(i)<100
        filename = ['ExptA000' num2str(stimID(i)) '.png'];
    end
    
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    
    [Ny Nx] = size(IMAGEorg);
    
    IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter (anti-aliasing)
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    IMAGE = edgetaper(IMAGE,edgeGauss); %taper image edges
    
    IMAGEp = cat(1,zeros(nzpad,size(IMAGE,2)),IMAGE);
    IMAGEp = cat(1,IMAGEp,zeros(nzpad,size(IMAGEp,2)));
    IMAGEp = cat(2,zeros(size(IMAGEp,1),nzpad),IMAGEp);
    IMAGEp = cat(2,IMAGEp,zeros(size(IMAGEp,1),nzpad));
    
    hor_energy = nan(length(gridY),length(gridX),length(hor_freqrange));
    ver_energy = nan(length(gridY),length(gridX),length(ver_freqrange));
    diagur_energy = nan(length(gridY),length(gridX),length(diagur_freqrange));
    diagul_energy = nan(length(gridY),length(gridX),length(diagul_freqrange));
    for xx = 1:length(gridX)
        for yy = 1:length(gridY)
            IMchunk = IMAGEp(win_y+gridY(yy)+nzpad,win_x+gridX(xx)+nzpad);
            F = fftshift(fft2(IMchunk.*gaussfun));
            for ii = 1:length(hor_freqrange)
                hor_energy(yy,xx,ii) = sum(abs(F(hor_freqrange{ii})));
                ver_energy(yy,xx,ii) = sum(abs(F(ver_freqrange{ii})));
                diagur_energy(yy,xx,ii) = sum(abs(F(diagur_freqrange{ii})));
                diagul_energy(yy,xx,ii) = sum(abs(F(diagul_freqrange{ii})));
            end
        end
    end
    
       subplot(3,2,1)
       imagesc(IMAGE);colormap(gray)
       subplot(3,2,3)
       imagesc(squeeze(hor_energy(:,:,1)));colormap(gray);
       subplot(3,2,4)
       imagesc(squeeze(ver_energy(:,:,1)));colormap(gray);
       subplot(3,2,5);
       imagesc(squeeze(diagur_energy(:,:,1)));colormap(gray);
       subplot(3,2,6);
       imagesc(squeeze(diagul_energy(:,:,1)));colormap(gray);
       pause(2);
       clf
    
    energy_array = cat(4,hor_energy,ver_energy,diagul_energy,diagur_energy);
    
    all_Energy_array = cat(5,all_Energy_array,energy_array);
    
end

%% normalize each set of energy maps to zscore across images/pixels
for i = 1:size(all_Energy_array,3)
    for j = 1:size(all_Energy_array,4)
        cur_set = all_Energy_array(:,:,i,j,:);
        cur_set_norm = (cur_set - mean(cur_set(:)))/std(cur_set(:));
        all_Energy_array(:,:,i,j,:) = cur_set_norm;
    end
end

%%
cd ~/Data/bruce/2_27_12/stimrecon/
save all_image_Energyarrays_block1 all_Energy_array grid* dsfrac Fsd sigx sigy *_freqrange bandpass*