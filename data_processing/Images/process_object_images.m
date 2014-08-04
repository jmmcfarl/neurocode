clear all
close all

outdirectory = '~/James_scripts/data_processing/Images/processed_images2';
obj_directory = '~/James_scripts/data_processing/Images/Object_images/';
cd(obj_directory)
white = 255;
black = 0;
buffer = 10;
max_resc_factor = 2;
resc_bounds = 10;
max_rot = 30;
bound_thick = 1;

dd = dir('./');
dd = dd(4:end);
% for ii = 1:length(dd)
    
filename = 'zebra.jpg';
% filename = dd(ii).name;
disp(filename);
imdata = imread(filename);
imdata = rgb2gray(imdata);
imdata = double(imdata); % convert to double format

%random rotation
offset = 1;
imdata_r = imrotate(imdata+offset,rand*2*max_rot-max_rot);
imdata_r = imdata_r - offset;
imdata_r(imdata_r == -1) = white;
imdata = imdata_r;

inmask = find(imdata < white - buffer);
outmask = find(imdata >= white - buffer);
mask = ones(size(imdata));
mask(outmask) = 0;
[B,L,N,A] = bwboundaries(mask);
enclosed_set = find(A(:,1));
non_enclosed_set = setdiff(1:length(B),enclosed_set);
L(ismember(L,non_enclosed_set)) = 1;
obj_boundary = B{1};

% %%
% % imshow(label2rgb(L, @jet, [.5 .5 .5]))
% % hold on
% % plot(obj_boundary(:,2),obj_boundary(:,1),'r','linewidth',2)
% %%
in_obj = find(L==1);
[sz_y,sz_x] = size(imdata);
[X,Y] = meshgrid(1:sz_x,1:sz_y);
obj_pos = [Y(in_obj) X(in_obj)];
[D,I] = pdist2(obj_boundary,obj_pos,'euclidean','Smallest',1);
boundary_pix = in_obj(D <= bound_thick);
L(L > 1) = 0;
L_old = L;

L(boundary_pix) = 0;
L_old(boundary_pix) = 2;

figure
showIm(L_old);


bounds = prctile(imdata(in_obj),[resc_bounds (100-resc_bounds)]);
resc_factor = (white-black)/diff(bounds);
resc_factor = min(resc_factor,max_resc_factor);
resc_imdata = imdata*resc_factor;
resc_imdata = (resc_imdata - bounds(1) + black);
resc_imdata(resc_imdata < black) = black;
resc_imdata(resc_imdata > white) = white;
% 

% figure
% subplot(2,1,1)
% showIm(imdata)
% hold on
% plot(obj_boundary(:,2),obj_boundary(:,1),'r','linewidth',1.5)
% subplot(2,1,2)
% showIm(resc_imdata)
% hold on
% plot(obj_boundary(:,2),obj_boundary(:,1),'r','linewidth',1.5)

imdata = resc_imdata;

cd(outdirectory)
back_num = ceil(rand*300);
back_filename = sprintf('sim%.4d.png',back_num);
bk_imdata = imread(back_filename);
noise_patt = double(bk_imdata); % convert to double format
[Ny,Nx] = size(bk_imdata);

% Nx = 1280; Ny = 1024;
% % dim = [Ny Ny];
% dim = [Ny Nx];
% noise_patt = spatialPattern(dim,-2);
% Pix2Deg = 0.018837;
% f = highpassfilter(dim,0.2*Pix2Deg,2); %cutoff freq at 2 cycle/degree, assuming 60 pix/degree
% MarginX = (Nx - dim(1))/2;
% noise_patt = real(ifft2(fft2(noise_patt).*f));
% noise_std = 70;
% noise_mean = (white+black)/2;
% noise_patt = noise_patt/std(noise_patt(:))*noise_std + noise_mean;
% noise_patt(noise_patt > white) = white;
% noise_patt(noise_patt < black) = black;

% %hist equalization
% [val indx] = sort(noise_patt(:));
% noise_patt = noise_patt;
% noise_patt(indx) = linspace(black,white,numel(noise_patt));

% bufferX = ones(dim(1),MarginX)*(white+black)/2;
% noise_patt = [bufferX noise_patt bufferX];

[Nyo,Nxo] = size(L);
bk_L = zeros(size(noise_patt));
bk_L(round(Ny/2-Nyo/2):round(Ny/2+Nyo/2)-1,round(Nx/2-Nxo/2):round(Nx/2+Nxo/2)-1) = L;

new_im = noise_patt;
new_im(bk_L==1) = imdata(L==1);

figure
showIm(new_im);

% input('')
% close all
% % end
%%
cd(obj_directory) 
filename = 'cat.jpg';
imdata = imread(filename);
imdata = rgb2gray(imdata);
imdata = double(imdata); % convert to double format

%random rotation
offset = 1;
imdata_r = imrotate(imdata+offset,rand*2*max_rot-max_rot);
imdata_r = imdata_r - offset;
imdata_r(imdata_r == -1) = white;
imdata = imdata_r;

inmask = find(imdata < white - buffer);
outmask = find(imdata >= white - buffer);
mask = ones(size(imdata));
mask(outmask) = 0;
[B,L1,N,A] = bwboundaries(mask);
enclosed_set = find(A(:,1));
non_enclosed_set = setdiff(1:length(B),enclosed_set);
L1(ismember(L1,non_enclosed_set)) = 1;
bounds = prctile(imdata(L1==1),[resc_bounds (100-resc_bounds)]);
resc_factor = (white-black)/diff(bounds);
resc_factor = min(resc_factor,max_resc_factor);
resc_imdata = imdata*resc_factor;
resc_imdata = (resc_imdata - bounds(1) + black);
resc_imdata(resc_imdata < black) = black;
resc_imdata(resc_imdata > white) = white;
imdata1 = resc_imdata;

filename = 'zebra.jpg';
imdata = imread(filename);
imdata = rgb2gray(imdata);
imdata = double(imdata); % convert to double format

%random rotation
offset = 1;
imdata_r = imrotate(imdata+offset,rand*2*max_rot-max_rot);
imdata_r = imdata_r - offset;
imdata_r(imdata_r == -1) = white;
imdata = imdata_r;

inmask = find(imdata < white - buffer);
outmask = find(imdata >= white - buffer);
mask = ones(size(imdata));
mask(outmask) = 0;
[B,L2,N,A] = bwboundaries(mask);
enclosed_set = find(A(:,1));
non_enclosed_set = setdiff(1:length(B),enclosed_set);
L2(ismember(L2,non_enclosed_set)) = 1;
bounds = prctile(imdata(L2==1),[resc_bounds (100-resc_bounds)]);
resc_factor = (white-black)/diff(bounds);
resc_factor = min(resc_factor,max_resc_factor);
resc_imdata = imdata*resc_factor;
resc_imdata = (resc_imdata - bounds(1) + black);
resc_imdata(resc_imdata < black) = black;
resc_imdata(resc_imdata > white) = white;
imdata2 = resc_imdata;

% cd(outdirectory)
% back_filename = 'sim0001.png';
% bk_imdata = imread(back_filename);
% bk_imdata = double(bk_imdata); % convert to double format
% [Ny,Nx] = size(bk_imdata);
Nx = 1280; Ny = 1024;
% dim = [Ny Ny];
dim = [Ny Nx];
noise_patt = spatialPattern(dim,-2);

Pix2Deg = 0.018837;
f = highpassfilter(dim,0.2*Pix2Deg,2); %cutoff freq at 2 cycle/degree, assuming 60 pix/degree
MarginX = (Nx - dim(1))/2;
noise_patt = real(ifft2(fft2(noise_patt).*f));
[val indx] = sort(noise_patt(:));
noise_patt = noise_patt;
noise_patt(indx) = linspace(black,white,numel(noise_patt));

% bufferX = ones(dim(1),MarginX)*(white+black)/2;
% noise_patt = [bufferX noise_patt bufferX];

x_offset = 150;
[Ny1,Nx1] = size(L1);
[Ny2,Nx2] = size(L2);
bk_L1 = zeros(size(noise_patt));
bk_L2 = zeros(size(noise_patt));
bk_L1(round(Ny/2-Ny1/2):round(Ny/2+Ny1/2)-1,round(Nx/2-Nx1/2)-x_offset:round(Nx/2+Nx1/2)-1-x_offset) = L1;
bk_L2(round(Ny/2-Ny2/2):round(Ny/2+Ny2/2)-1,round(Nx/2-Nx2/2)+x_offset:round(Nx/2+Nx2/2)-1+x_offset) = L2;

new_im = noise_patt;
if rand > 0.5
    new_im(bk_L1==1) = imdata1(L1==1);
    new_im(bk_L2==1) = imdata2(L2==1);
else
    new_im(bk_L2==1) = imdata2(L2==1);
    new_im(bk_L1==1) = imdata1(L1==1);
end
showIm(new_im)
