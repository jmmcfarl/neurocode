% Load and pre-process the image files

clear all
addpath('~/James_scripts/data_processing/Images/');
addpath(genpath('~/James_scripts/textureSynth/'));
addpath(genpath('~/James_scripts/matlabPyrTools/'));

imagedirectory = '~/James_scripts/data_processing/Images/vanhateren_imc';
outdirectory = '~/James_scripts/data_processing/Images/processed_images';

Xscale = 1; Yscale=1; % rescale patch

addpath(imagedirectory);
imagefiles = dir(imagedirectory);

%compile all image file names
imagename = {};
j=1;
for i=1:length(imagefiles)
    if isempty(strfind(imagefiles(i).name, '.imc'))
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

screenRect = [0 0 1280 1024];
% windowRect = [0 0 1024 1024];
windowRect = [0 0 1280 1024];
blackColor = 0;
whiteColor = 255;
Nphoto = length(imagename);
imagetex = zeros(Nphoto,1);

% center pixel before scale transformation
XcenterPix = windowRect(3)/2;
YcenterPix = windowRect(4)/2;
grayColor = (whiteColor + blackColor)/2;

%%
% stim_set = 1035:Nphoto;
cd(outdirectory)
dd = dir();
dd(1:3) = [];
stim_set = zeros(length(dd),1);
for i = 1:length(dd)
   stim_set(i) = str2num(dd(i).name(4:8)); 
end

cd(imagedirectory)
for i=1:length(stim_set)
    fprintf('Processing image %d of %d\n',stim_set(i),Nphoto);
    
    myimgfile = imagename(stim_set(i),:);
    myimgfile = myimgfile{1};
    
    % prepare the image
    f1=fopen(myimgfile,'rb','ieee-be');
    w=1536;h=1024;
    imdata = fread(f1,[w,h],'uint16');
    imdata = imdata';
    fclose(f1);
    
    % crop images
    [Ly Lx] = size(imdata);
    MarginX = (windowRect(3) - Lx)/2;
    assert( mod(MarginX,1)==0)
    MarginY = (windowRect(4) - Ly)/2;
    assert( mod(MarginY,1)==0)
    if MarginY>0
        bufferY = ones(MarginY, Lx)*grayColor;
        imdata = [bufferY; imdata; bufferY];
    else
        imdata = imdata((-MarginY+1):(Ly+MarginY),:);
    end
    if MarginX>0
        bufferX = ones(MarginX, size(imdata,1))*grayColor;
        imdata = [bufferX imdata bufferX];
    else
        imdata = imdata(:,(-MarginX+1):(Lx+MarginX));
    end
    
    imdataorg = imdata;
    [Ny Nx] = size(imdata);
    
    %histogram equalization
    [val indx] = sort(imdata(:));
    imdatarescale = imdata;
    imdatarescale(indx) = linspace(blackColor,whiteColor,numel(imdata));
    
    imdata = imdatarescale;
    
    %% add buffers to fill screen
    % adjust margins to fit screen
    [Ly Lx] = size(imdata);
    MarginX = (screenRect(3) - Lx)/2;
    assert( mod(MarginX,1)==0)
    MarginY = (screenRect(4) - Ly)/2;
    assert( mod(MarginY,1)==0)
    bufferY = ones(MarginY, Lx)*grayColor;
    imdata = [bufferY; imdata; bufferY];
    
    bufferX = ones(size(imdata,1),MarginX)*grayColor;
    imdata = [bufferX imdata bufferX];
    
    % cutoff color
    imdata(imdata>whiteColor)=whiteColor;
    imdata(imdata<blackColor)=blackColor;
    
    %print raw image
    imdatapng = imdata./whiteColor;
    cur_fname = myimgfile(1:end-4);
    cur_fname(1:3) = 'raw';
    imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
        
end
