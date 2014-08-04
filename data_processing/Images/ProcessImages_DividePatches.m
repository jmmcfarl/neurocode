% Load and pre-process the image files
% The user need to specify directory of the raw image files ('.imc' format)
% The script would
% (1) load the image files as a matlab matrix
% (2) Crop the image according to the desired screen size
% (3) Break image into patches and rescale it (specified by gridx/gridy;
% Xscale/Yscale
% (4) Perform luminance balance operation for each patch
% (5-a) Save the processed data as '.mat' files in /patches subfolder
% (5-b) save the processed data as '.png' files in /Patches_png subfolder
% (5-c) save the best 60% images as '.png' files in /GoodPatches_png subfolder
% Yuwei Cui, Last Modified Feb 8 2012

clear all
imagedirectory = './imgfile_Example';
imagedirectory_new = [ imagedirectory '/patches'];
imagedirectory_new_png = [imagedirectory '/Patches_png'];
imagedirectory_good_png = [imagedirectory '/GoodPatches_png'];

% gridx = 4; gridy = 4; % divide into gridx X gridy patches
% taskParams.Xscale = 4; taskParams.Yscale=4; % rescale patch

gridx = 1; gridy = 1; % divide into gridx X gridy patches
taskParams.Xscale = 1; taskParams.Yscale=1; % rescale patch

addpath(imagedirectory);
outdirs = {imagedirectory_new,imagedirectory_new_png, imagedirectory_good_png};
for diri = 1:length(outdirs)
    outdir = outdirs{diri};
    if system(['cd ' outdir])~=0
        system(['mkdir ' outdir]);
    end
end
imagefiles = dir(imagedirectory);
mode = 'conversion';
% mode = 'display';
DISPLAYTYPE = 'matlab_figure';
imagename = {};
j=1;
for i=1:length(imagefiles)
    if isempty(strfind(imagefiles(i).name, '.imc')) || ~isempty(strfind(imagefiles(i).name, '.mat'))||~isempty(strfind(imagefiles(i).name, '.jpeg'))
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
stim.imagename = imagename;
if strcmp(DISPLAYTYPE,  'psychophysics_toolbox')
    stim.displayParams = CPLabDisplayParamsWithGammCorr;
else
    stim.displayParams.screenRect = [0 0 1280 1024];
    stim.displayParams.windowRect = [0 0 1280 1024];
    stim.displayParams.blackColor = 0;
    stim.displayParams.whiteColor = 255;
end
stim.Nphoto = length(imagename);
stim.imagetex = zeros(stim.Nphoto,1);


% center pixel before scale transformation
XcenterPix = stim.displayParams.windowRect(3)/2;
YcenterPix = stim.displayParams.windowRect(4)/2;
stim.displayParams.grayColor = ( stim.displayParams.whiteColor + stim.displayParams.blackColor)/2;


cost = zeros(stim.Nphoto, gridx, gridy);

for i=1:stim.Nphoto
    myimgfile =stim.imagename(i,:);
    myimgfile = myimgfile{1};
    
    if strcmp(mode,'conversion')
        % prepare the image
        f1=fopen(myimgfile,'rb','ieee-be');
        %     f1=fopen(myimgfile,'rb');
        w=1536;h=1024;
        imdata = fread(f1,[w,h],'uint16');
        imdata = imdata';
        imdata = imdata-(median(imdata(:))-stim.displayParams.grayColor);
        imdata = imdata/(prctile(imdata(:),95)-prctile(imdata(:),5))*stim.displayParams.whiteColor;
%         imdata(imdata(:)>stim.displayParams.whiteColor) = stim.displayParams.whiteColor;
%         imdata(imdata(:)<stim.displayParams.blackColor) = stim.displayParams.blackColor;
        %         load(strcat(imagedirectory, '/', myimgfile, '.mat'));
        
        % crop images
        [Ly Lx] = size(imdata);
        MarginX = (stim.displayParams.screenRect(3) - Lx)/2;
        assert( mod(MarginX,1)==0)
        MarginY = (stim.displayParams.screenRect(4) - Ly)/2;
        assert( mod(MarginY,1)==0)
        
        if MarginY>0
            bufferY = ones(MarginY, Lx)*stim.displayParams.grayColor;
            imdata = [bufferY; imdata; bufferY];
        else
            imdata = imdata((-MarginY+1):(Ly+MarginY),:);
        end
        
        if MarginX>0
            bufferX = ones(MarginX, size(imdata,1))*stim.displayParams.grayColor;
            imdata = [bufferX imdata bufferX];
        else
            imdata = imdata(:,(-MarginX+1):(Lx+MarginX));
        end
        
        imdataorg = imdata;
        %         fprintf('Original Size: %d  X %d\n',Nx,Ny);
        [Ny Nx] = size(imdata);
        Lx = floor(Nx/gridx);
        Ly = floor(Ny/gridy);
    end
    
    for xi=1:gridx
        for yi=1:gridy
            if strcmp(mode,'conversion')
                                
                imdata = imdataorg(Ly*(yi-1)+(1:Ly), Lx*(xi-1)+(1:Lx));
                imdata = imdata-(median(imdata(:))-stim.displayParams.grayColor);
                imdata = imdata/(prctile(imdata(:),95)-prctile(imdata(:),5))*stim.displayParams.whiteColor;

                f = highpassfilter(size(imdata),0.0167,4); %cutoff freq at 1 cycle/degree, assuming 60 pix/degree
                im_filt = real(ifft2(fft2(imdata).*f));
                imdata = im_filt;
                
                %                 imdata(imdata(:)>stim.displayParams.whiteColor) = stim.displayParams.whiteColor;
%                 imdata(imdata(:)<stim.displayParams.blackColor) = stim.displayParams.blackColor;
                
                [val indx] = sort(imdata(:));
                imdatarescale = zeros(size(imdata));
                binN = 50; pixperbin = floor(length(imdata(:))/binN);
                lumscaleperbin = (stim.displayParams.whiteColor-stim.displayParams.blackColor)/binN;
                for ii=1:binN
                    indxindx = (ii-1)*pixperbin + (1:pixperbin);
                    pixels = imdata(indx(indxindx));
                    pixels = pixels - min(pixels); pixels = pixels/max(pixels);
                    imdatarescale(indx(indxindx)) = pixels * lumscaleperbin + (ii-1)*lumscaleperbin + stim.displayParams.blackColor;
                end
                distorg = histc(imdata(:),0:lumscaleperbin:stim.displayParams.whiteColor); distorg = distorg/sum(distorg(:));
                distnew = histc(imdatarescale(:),0:lumscaleperbin:stim.displayParams.whiteColor); distnew = distnew/sum(distnew(:));
                costnow = sum((distorg-distnew).^2);
                cost(i,xi,yi) = costnow;
                
                fprintf('IMage %d X %d Y %d Cost %3.4f\n',i,xi,yi,costnow);
                imdata = imdatarescale;

                                %                 % adjust mean luminance to gray
                %                 imdata = 0.8*stim.displayParams.grayColor*imdata./mean(mean(imdata));
                %                 % adjust contrast
                %                 imdata = imdata./std(std(imdata))*15;
                % cutoff color
                imdata(imdata>stim.displayParams.whiteColor)=stim.displayParams.whiteColor;
                imdata(imdata<stim.displayParams.blackColor)=stim.displayParams.blackColor;
                
                [Nx, Ny] = size(imdata);
                imdata = imdata(floor((0:Nx*taskParams.Xscale-1)/taskParams.Xscale)+1,floor((0:Ny*taskParams.Yscale-1)/taskParams.Yscale)+1);
                
                figure(1);clf;plot(distnew);hold on;plot(distorg,'r');title(['cost = ',num2str(costnow)]);
                %                 figure(1);clf;hist(imdata(:),lumscaleperbin/2:lumscaleperbin:stim.displayParams.whiteColor); axis tight;pause(.1);
                
                save(strcat(imagedirectory_new, '/', myimgfile, '.x', num2str(xi), 'y', num2str(yi), '.mat'),'imdata');
            else
                load(strcat(imagedirectory_new, '/', myimgfile, '.x', num2str(xi), 'y', num2str(yi), '.mat'),'imdata');
            end
            
            
            switch DISPLAYTYPE
                case 'psychophysics_toolbox'
                    % desired size (after scaling)
                    sx = size(imdata,1);%*taskParams.Xscale;
                    sy = size(imdata,1);%*taskParams.Yscale;
                    stim.destrec =[XcenterPix-sx/2,YcenterPix-sy/2,XcenterPix+sx/2,YcenterPix+sy/2];
                    stim.imagetex(i)=Screen('MakeTexture', stim.displayParams.window, imdata);
                    Screen('DrawTexture', stim.displayParams.window, stim.imagetex(i),[],stim.destrec);
                    Screen(stim.displayParams.window, 'Flip');
                case 'matlab_figure'
                    figure(55);clf;imagesc(imdata);colormap gray; axis tight;
            end
            
            imdatapng = imdata./stim.displayParams.whiteColor;
            imwrite(imdatapng,strcat(imagedirectory_new_png, '/', myimgfile, '.x', num2str(xi), 'y', num2str(yi), '.png'),'png');
            %             figure(1),imagesc(imdata);colormap gray
            %             axis tight;set(gca,'Xtick',[]);set(gca,'Ytick',[]);set(gca,'Ydir','reversed');
            pause(.5)
            %             X=screen('GetImage',stim.displayParams.window);
            %             f=im2frame(X);
            %             M=addframe(M,f);
        end
    end
end

% Find the worst 10% of images and delete them
[costsort costind] = sort(cost(:));
for i=1:round(length(cost(:))*0.6)
    index = costind(i);
    imagei = ceil(index/(gridx*gridy));
    res = index-(imagei-1)*(gridx*gridy);
    gridxi = ceil(res/gridy);
    if gridxi>1
        gridyi = res-(gridxi-1)*gridy;
    else
        gridyi = res;
    end
    myimgfile =stim.imagename(imagei,:);
    myimgfile = myimgfile{1};
    load(strcat(imagedirectory_new, '/', myimgfile, '.x', num2str(gridxi), 'y', num2str(gridyi), '.mat'),'imdata');
    
    % desired field to present the stimuli
    switch DISPLAYTYPE
        case 'psychophysics_toolbox'
            sx = size(imdata,1)*taskParams.Xscale;
            sy = size(imdata,1)*taskParams.Yscale;
            stim.destrec =[XcenterPix-sx/2,YcenterPix-sy/2,XcenterPix+sx/2,YcenterPix+sy/2];
            stim.imagetex(i)=Screen('MakeTexture', stim.displayParams.window, imdata);
            Screen('DrawTexture', stim.displayParams.window, stim.imagetex(i),[],stim.destrec);
            Screen(stim.displayParams.window, 'Flip');
        case 'matlab_figure'
            figure(55);clf;imagesc(imdata);colormap gray; axis tight;
    end
    imdatapng = imdata./stim.displayParams.whiteColor;
    imwrite(imdatapng,strcat(imagedirectory_good_png, '/', myimgfile, '.x', num2str(gridxi), 'y', num2str(gridyi), '.png'),'png');
    fprintf('IMage %d X %d Y %d Cost %3.4f\n',imagei,gridxi,gridyi,costsort(i));
end
clear mex
% M=close(M);