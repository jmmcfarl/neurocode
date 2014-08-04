% reconstruct the stimulus
% created by Yuwei Cui, Feb 28, 2012
% edited by James McFarland

clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/

%%
load Blocks.mat

% original image resolution
Pix2Deg = 1.1279 / 60;

% down-sampling fraction for image
dsfrac = 4;

% block number
blockid = 1;

stimtime = Blocks{blockid}.stimtime;
stimID = Blocks{blockid}.stimids;
Nimage = length(stimID);

gaussfilt = fspecial('gaussian',5,0.25);

all_Image_array = [];

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
    [Nyd Nxd] = size(IMAGE); %down-sampled image size

    %if you want to reflect the image about x or y axes
%     IMAGE = flipud(IMAGE);
%     IMAGE = fliplr(IMAGE);
      
    all_Image_array = cat(3,all_Image_array,IMAGE);

end

cd ~/Data/bruce/2_27_12/stimrecon/
save all_image_mat_ds4 all_Image_array