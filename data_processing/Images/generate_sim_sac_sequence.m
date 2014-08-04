clear all
close all
addpath('~/James_scripts/data_processing/Images/');
addpath(genpath('~/James_scripts/textureSynth/'));
addpath(genpath('~/James_scripts/matlabPyrTools/'));

screenRect = [0 0 1280 1024];
windowRect = [0 0 1280 1024];
black = 0;
white = 255;
gray = (white + black)/2;
% center pixel before scale transformation
XcenterPix = windowRect(3)/2;
YcenterPix = windowRect(4)/2;

Pix2Deg = 0.018837;
min_amp = 0.05;
max_amp = 5;
mu = 2.5;
lambda = 1.4;
max_disp = 6;

sacs_per_im = 4;

%%
imagedirectory = '~/James_scripts/data_processing/Images/processed_images2/';
outdirectory = '~/James_scripts/data_processing/Images/processed_image_seqs';
imagefiles = dir(imagedirectory);

%compile all image file names
imagename = {};
j=1;
for i=1:length(imagefiles)
    if isempty(strfind(imagefiles(i).name, '.png'))
        % if isempty(strfind(imagefiles(i).name, '.imc.rescale'))
        continue;
    else
        if ~exist('imagename')
            imagename = imagefiles(i).name;
        else
            imagename = [imagename; imagefiles(i).name;];
        end
        j=j+1;
    end
end

Nphoto = length(imagename);
imagetex = zeros(Nphoto,1);
im_set = 1:Nphoto;
cur_cnt = 1;

figure
colormap gray
%%
while ~isempty(im_set)
    cd(imagedirectory)
    orig_nimages = length(imagename);
    cur_choice = ceil(rand*length(im_set));
    
    myimgfile = imagename(im_set(cur_choice),:);
    myimgfile = myimgfile{1};
    
    fprintf('Processing image %d of %d\n',cur_choice,orig_nimages);
    
    imdata = imread(myimgfile);
    imdata = double(imdata); % convert to double format
    
    imdataorg = imdata;
    
    %%
    
    if str2num(myimgfile(1)) ~= 4
    %keep trying to generate a sequence of sacs that stays within the window.
    bad_pts = 1; %initialize
    while ~isempty(bad_pts)
        sac_amps = randraw('ig',[mu lambda],round(4*sacs_per_im));
        sac_amps(sac_amps > max_amp) = [];
        
        rand_dirs = rand(sacs_per_im,1)*2*pi;
        X = zeros(sacs_per_im+1,2);
        its = 0;
        for i = 2:sacs_per_im+1
            X(i,1) = X(i-1,1) + sac_amps(i-1)*cos(rand_dirs(i-1));
            X(i,2) = X(i-1,2) + sac_amps(i-1)*sin(rand_dirs(i-1));
        end
        bad_pts = find(abs(X(:,1)) > max_disp | abs(X(:,2)) > max_disp);
        its = its + 1;
    end
    else
        max_rdisp = 3.5;
       X = rand(sacs_per_im+1,2)*2*max_rdisp-max_rdisp; 
    end
        
%     myimgfile(1)
    %if this isn't a natural image, start out after one saccade from the
    %center
%     if str2num(myimgfile(1)) ~= 1
        X(1,:) = []; %get rid of initial point
%     end
  
    X_inds = round(X/Pix2Deg);

    %%
    [Ny,Nx] = size(imdata);
    cd(outdirectory)
    for cur_fix = 1:sacs_per_im
        shift_im = shift_mat2(imdataorg,X_inds(cur_fix,1),X_inds(cur_fix,2),1);
        imdata = shift_im;
        %                 imagesc(imdata)
        %                 hold on
        %                 plot(Nx/2,Ny/2,'b+','linewidth',2,'markersize',8)
        %                 x_cent = round(0.35/Pix2Deg);
        %                 y_cent = round(-0.45/Pix2Deg);
        %                 x_width = round(0.7/Pix2Deg);
        %                 y_width = round(0.7/Pix2Deg);
        %                 x_cent = x_cent - x_width/2;
        %                 y_cent = y_cent - y_width/2;
        %                 rectangle('Position',[Nx/2+x_cent Ny/2+y_cent x_width y_width],'linewidth',1,'edgecolor','r')
        %                 x_cent = 0;
        %                 y_cent = 0;
        %                 x_width = round(5/Pix2Deg);
        %                 y_width = round(5/Pix2Deg);
        %                 x_cent = x_cent - x_width/2;
        %                 y_cent = y_cent - y_width/2;
        %                 rectangle('Position',[Nx/2+x_cent Ny/2+y_cent x_width y_width],'linewidth',1,'edgecolor','w')
        %                 axis off
        %                 pause(0.5)
        
        %print raw image
        imdatapng = imdata./white;
        cur_fname = sprintf('1%.4d',cur_cnt);
        imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
        
        sim_sac_data(cur_cnt).im_name = myimgfile(1:end-4);
        sim_sac_data(cur_cnt).X_trans = X(cur_fix,1);
        sim_sac_data(cur_cnt).Y_trans = X(cur_fix,2);
        
        cur_cnt = cur_cnt + 1;
        
    end
    
    %%
    
    
    im_set(cur_choice) = [];
end

save sim_sac_data sim_sac_data

    