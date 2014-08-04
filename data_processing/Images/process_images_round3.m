% Load and pre-process the image files

clear all
addpath('~/James_scripts/data_processing/Images/');
addpath(genpath('~/James_scripts/textureSynth/'));
addpath(genpath('~/James_scripts/matlabPyrTools/'));

imagedirectory = '~/James_scripts/data_processing/Images/processed_images/';
outdirectory = '~/James_scripts/data_processing/Images/processed_images2';

Xscale = 1; Yscale=1; % rescale patch

addpath(imagedirectory);
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

screenRect = [0 0 1280 1024];
windowRect = [0 0 1024 1024];
blackColor = 0;
whiteColor = 255;
Nphoto = length(imagename);
imagetex = zeros(Nphoto,1);

% center pixel before scale transformation
XcenterPix = windowRect(3)/2;
YcenterPix = windowRect(4)/2;
grayColor = (whiteColor + blackColor)/2;

%% simoncelli parameters
si_Nsc = 6; % Number of scales
si_Nor = 4; % Number of orientations
si_Na = 5;  % Spatial neighborhood is Na x Na coefficients
si_Niter = 10;	% Number of iterations of synthesis loop

%%
cd(imagedirectory)
orig_nimages = length(imagename);
im_set = 1:orig_nimages;
cur_cnt = 1;


while length(im_set) > orig_nimages/2 %while more than half images remain
    cur_choice = ceil(rand*length(im_set));
    
    myimgfile = imagename(im_set(cur_choice),:);
    myimgfile = myimgfile{1};
    
    fprintf('Processing image %d of %d\n',cur_choice,orig_nimages);
    
    imdata = imread(myimgfile);
    imdata = double(imdata); % convert to double format
    
    imdataorg = imdata;
    
%     % crop images
%     [Ly Lx] = size(imdata);
%     MarginX = (windowRect(3) - Lx)/2;
%     assert( mod(MarginX,1)==0)
%     MarginY = (windowRect(4) - Ly)/2;
%     assert( mod(MarginY,1)==0)
%     if MarginY>0
%         bufferY = ones(MarginY, Lx)*grayColor;
%         imdata_cr = [bufferY; imdata; bufferY];
%     else
%         imdata_cr = imdata((-MarginY+1):(Ly+MarginY),:);
%     end
%     if MarginX>0
%         bufferX = ones(MarginX, size(imdata,1))*grayColor;
%         imdata_cr = [bufferX imdata bufferX];
%     else
%         imdata_cr = imdata(:,(-MarginX+1):(Lx+MarginX));
%     end
    
    %         fprintf('Original Size: %d  X %d\n',Nx,Ny);
    [Ny Nx] = size(imdata);

    %% now generate phase-randomized version
    RandomPhase = angle(fft2(rand(Ny, Nx)));
    %generate random phase structure
    ImFourier = fft2(imdata);
    %Fast-Fourier transform
    Amp = abs(ImFourier);
    %amplitude spectrum
    Phase = angle(ImFourier);
    %phase spectrum
    Phase = Phase + RandomPhase;
    %add random phase to original phase
    ImScrambled = ifft2(Amp.*exp(sqrt(-1)*(Phase)));
    %combine Amp and Phase then perform inverse Fourier
    ImScrambled = real(ImScrambled); %get rid of imaginery
    
    %histogram equalization on phase randomized image
    [val indx] = sort(ImScrambled(:));
    imdatarescale = ImScrambled;
    imdatarescale(indx) = linspace(blackColor,whiteColor,numel(imdatarescale));
    
%     tfreqs   = [(Nx/2):-1:-((Nx/2)-1)]/Nx; %spatial frequency axis (nyquist to minus nyquist)
%     [Fx,Fy] = meshgrid(tfreqs);
%     Fabs = sqrt(Fx.^2+Fy.^2);
%     nbins = 100;
%     bin_edges = prctile(Fabs(:),linspace(1,99,nbins+1));
%     bin_centers = 0.5*bin_edges(2:end)+0.5*bin_edges(1:end-1);
%     
%     D1 = fftshift(fft2(imdata)); %compute 2D FFT
%     rabs = abs(D1); %amplitude spectrum
%     D1sc = fftshift(fft2(ImScrambled)); %compute 2D FFT
%     rabssc = abs(D1sc); %amplitude spectrum
%     
%     avg_spec = zeros(size(bin_centers));
%     avg_specsc = zeros(size(bin_centers));
%    for ff = 1:length(bin_centers)
%         temp = find(Fabs >= bin_edges(ff) & Fabs < bin_edges(ff+1));
%         avg_spec(ff) = mean(rabs(temp));
%         avg_specsc(ff) = mean(rabssc(temp));
%     end
%     figure
%     subplot(3,1,1)
%     imagesc(tfreqs,tfreqs,log10(rabs))
%     subplot(3,1,3)
%     plot(bin_centers,log10(avg_spec))
%     subplot(3,1,2)
%     imagesc(tfreqs,tfreqs,log10(rabssc))
%     subplot(3,1,3)
%     hold on
%     plot(bin_centers,log10(avg_specsc),'r')
    
    
    ImScrambled = imdatarescale;
    
    %% Now generate Simoncelli image
%     [Nyc Nxc] = size(imdata);
    
    params = textureAnalysis(imdata, si_Nsc, si_Nor, si_Na);
    si_Im = textureSynthesis(params, [Ny Nx], si_Niter);
    
    %histogram equalization 
    [val indx] = sort(si_Im(:));
    imdatarescale = si_Im;
    imdatarescale(indx) = linspace(blackColor,whiteColor,numel(si_Im));
    si_Im = imdatarescale;
    
%     % adjust margins to fit screen
%     [Ly Lx] = size(si_Im);
%     MarginX = (screenRect(3) - Lx)/2;
%     si_Im = [si_Im(:,end-MarginX+1:end) si_Im si_Im(:,1:MarginX)];
    
%     %% add buffers to fill screen
%     % adjust margins to fit screen
%     [Ly Lx] = size(imdata);
%     MarginX = (screenRect(3) - Lx)/2;
%     assert( mod(MarginX,1)==0)
%     MarginY = (screenRect(4) - Ly)/2;
%     assert( mod(MarginY,1)==0)
%     bufferY = ones(MarginY, Lx)*grayColor;
%     imdata = [bufferY; imdata; bufferY];
%     ImScrambled = [bufferY; ImScrambled; bufferY];
%     si_Im = [bufferY; si_Im; bufferY];
%     
%     bufferX = ones(size(imdata,1),MarginX)*grayColor;
%     imdata = [bufferX imdata bufferX];
%     ImScrambled = [bufferX ImScrambled bufferX];
%     si_Im = [bufferX si_Im bufferX];
    
    % cutoff color
    imdata(imdata>whiteColor)=whiteColor;
    imdata(imdata<blackColor)=blackColor;
    ImScrambled(ImScrambled>whiteColor)=whiteColor;
    ImScrambled(ImScrambled<blackColor)=blackColor;
    si_Im(si_Im>whiteColor)=whiteColor;
    si_Im(si_Im<blackColor)=blackColor;
    
    %print raw image
    imdatapng = imdata./whiteColor;
%     cur_fname = sprintf('raw%.4d',cur_cnt);
    cur_fname = sprintf('1%.4d',cur_cnt);
    imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
    
    ImScrambledpng = ImScrambled./whiteColor;
%     cur_fname = sprintf('phr%.4d',cur_cnt);
    cur_fname = sprintf('2%.4d',cur_cnt);
    imwrite(ImScrambledpng,strcat(outdirectory, '/', cur_fname, '.png'),'png');

    si_Impng = si_Im./whiteColor;
%     cur_fname = sprintf('sim%.4d',cur_cnt);
    cur_fname = sprintf('3%.4d',cur_cnt);
    imwrite(si_Impng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
    
    %%
    im_set(cur_choice) = [];
    cur_cnt = cur_cnt + 1;
    
end

%%
while ~isempty(im_set) %for all remaining images
    cur_choice = ceil(rand*length(im_set));
    
    myimgfile = imagename(im_set(cur_choice),:);
    myimgfile = myimgfile{1};
    
    fprintf('Processing image %d of %d\n',i,Nphoto);
    
    imdata = imread(myimgfile);
    imdata = double(imdata); % convert to double format
    % cutoff color
    imdata(imdata>whiteColor)=whiteColor;
    imdata(imdata<blackColor)=blackColor;
    
    %print raw image
    imdatapng = imdata./whiteColor;
%     cur_fname = sprintf('raw%.4d',cur_cnt);
    cur_fname = sprintf('1%.4d',cur_cnt);
    imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
        
    %%
    im_set(cur_choice) = [];
    cur_cnt = cur_cnt + 1;
    
end