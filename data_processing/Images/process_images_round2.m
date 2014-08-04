% Load and pre-process the image files

clear all
addpath('~/James_scripts/data_processing/Images/');
addpath(genpath('~/James_scripts/textureSynth/'));
addpath(genpath('~/James_scripts/matlabPyrTools/'));

imagedirectory = '/Users/James/Desktop/vanhateren_imc';
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
stim_set = 1:Nphoto;
for i=stim_set
    
    fprintf('Processing image %d of %d\n',i,Nphoto);
    
    myimgfile = imagename(i,:);
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
    %         fprintf('Original Size: %d  X %d\n',Nx,Ny);
    [Ny Nx] = size(imdata);
    
%     prct = prctile(imdata(:),[2.5 97.5]);
%     image_shift(i) = prct(1);
%     image_scale(i) = (prct(2)-prct(1));
%     imdata = imdata - image_shift(i);
%     imdata = imdata/image_scale(i)*whiteColor;
    %     imdata(imdata < blackColor) = blackColor;
    %     imdata(imdata > whiteColor) = whiteColor;
%     image_median(i) = median(imdata(:));
    
    %histogram equalization
    [val indx] = sort(imdata(:));
    imdatarescale = imdata;
    imdatarescale(indx) = linspace(blackColor,whiteColor,numel(imdata));
    
    imdata = imdatarescale;
    
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
    imdatascrescale = ImScrambled;
    imdatascrescale(indx) = linspace(blackColor,whiteColor,numel(imdatascrescale));
    
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
    
    
    ImScrambled = imdatascrescale;
    
    %% Now generate Simoncelli image
    params = textureAnalysis(imdata, si_Nsc, si_Nor, si_Na);
    si_Im = textureSynthesis(params, [Ny Nx], si_Niter);
    
    %histogram equalization 
    [val indx] = sort(si_Im(:));
    imdatascrescale = si_Im;
    imdatascrescale(indx) = linspace(blackColor,whiteColor,numel(si_Im));
    si_Im = imdatarescale;
    
    %% add buffers to fill screen
    % adjust margins to fit screen
    [Ly Lx] = size(imdata);
    MarginX = (screenRect(3) - Lx)/2;
    assert( mod(MarginX,1)==0)
    MarginY = (screenRect(4) - Ly)/2;
    assert( mod(MarginY,1)==0)
    bufferY = ones(MarginY, Lx)*grayColor;
    imdata = [bufferY; imdata; bufferY];
    ImScrambled = [bufferY; ImScrambled; bufferY];
    si_Im = [bufferY; si_Im; bufferY];
    
    bufferX = ones(size(imdata,1),MarginX)*grayColor;
    imdata = [bufferX imdata bufferX];
    ImScrambled = [bufferX ImScrambled bufferX];
    si_Im = [bufferX si_Im bufferX];
    
    % cutoff color
    imdata(imdata>whiteColor)=whiteColor;
    imdata(imdata<blackColor)=blackColor;
    ImScrambled(ImScrambled>whiteColor)=whiteColor;
    ImScrambled(ImScrambled<blackColor)=blackColor;
    si_Im(si_Im>whiteColor)=whiteColor;
    si_Im(si_Im<blackColor)=blackColor;
    
    %print raw image
    imdatapng = imdata./whiteColor;
    cur_fname = myimgfile(1:end-4);
    cur_fname(1:3) = 'raw';
    imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
    
    ImScrambledpng = ImScrambled./whiteColor;
    cur_fname = myimgfile(1:end-4);
    cur_fname(1:3) = 'phr';
    imwrite(ImScrambledpng,strcat(outdirectory, '/', cur_fname, '.png'),'png');

    si_Impng = si_Im./whiteColor;
    cur_fname = myimgfile(1:end-4);
    cur_fname(1:3) = 'sim';
    imwrite(si_Impng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
    
end
