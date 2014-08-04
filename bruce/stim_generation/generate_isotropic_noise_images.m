close all
clear all

save_dir = '~/Data/bruce/isotropic_noise_images';
blackColor = 0;
whiteColor = 255;
% cbounds = [25 75];
cbounds = [5 95];
% cbounds = [10 90];
grayColor = (blackColor+whiteColor)/2;
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
niqf = Fs/2;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);

%%
freq_cent = 4.5; %cyc/deg
freq_a = 9;
orange = [-180 180];
dec_const = 2/niqf;

n_images = 600;
for im_cnt = 1:n_images;
    fprintf('Generating image %d of %d\n',im_cnt,n_images);
    
    [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
    Fx = Fx*niqf;
    Fy = Fy*niqf;
    noise_im(Nyp+1:end,:) = [];
    
    % bounds = prctile(noise_im(:),cbounds);
    % noise_im = (noise_im-bounds(1))/range(bounds);
    % noise_im = noise_im*whiteColor;
    % noise_im(noise_im < blackColor) = blackColor;
    % noise_im(noise_im > whiteColor) = whiteColor;
    
    [val indx] = sort(noise_im(:));
    imdatarescale = noise_im;
    imdatarescale(indx) = linspace(blackColor,whiteColor,numel(noise_im));
    noise_im = imdatarescale;
    
    noise_png = noise_im/whiteColor;
    cur_fname = sprintf('%.4d',im_cnt);
    imwrite(noise_png,strcat(save_dir, '/', cur_fname, '.png'),'png');
end

%% ISOTROPIC NOISE SEQUENCES
freq_cent = 4.5; %cyc/deg
freq_a = 9;
orange = [-180 180];
dec_const = 2/niqf;

n_noise_sequences = 1;
ims_per_sequence = 200;
for seq = 1:n_noise_sequences
    cd ~/Data/bruce/
    fold_name = sprintf('noise_seq_%d',seq);
    if ~exist(fold_name,'dir')
        eval(['!mkdir ' fold_name])
    end
    cd(fold_name)
    
    for im_cnt = 1:ims_per_sequence
        [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
        noise_im(Nyp+1:end,:) = [];
        
        bounds = prctile(noise_im(:),cbounds);
        noise_im = (noise_im-bounds(1))/range(bounds);
        noise_im = noise_im*whiteColor;
        noise_im(noise_im < blackColor) = blackColor;
        noise_im(noise_im > whiteColor) = whiteColor;
        
%         [val indx] = sort(noise_im(:));
%         imdatarescale = noise_im;
%         imdatarescale(indx) = linspace(blackColor,whiteColor,numel(noise_im));
%         noise_im = imdatarescale;
        
        noise_png = noise_im/whiteColor;
        cur_fname = sprintf('%.4d',im_cnt);
        imwrite(noise_png,strcat(cur_fname, '.png'),'png');
    end
end

%% ANISOTROPIC NOISE SEQUENCES
poss_or_cents = -60:30:90;
or_width = 35;
n_poss_or_cents = length(poss_or_cents);

gauss_env_std = 75;
n_gauss = 100;
gfilt = fspecial('gaussian',[gauss_env_std*5 gauss_env_std*5],gauss_env_std);

n_noise_sequences = 1;
ims_per_sequence = 200;
for seq = 1:n_noise_sequences
    rand_seq = ceil(rand(ims_per_sequence,1)*n_poss_or_cents);
    orientation_sequence = poss_or_cents(rand_seq);
    cd ~/Data/bruce/
    fold_name = sprintf('aniso_noise_seq_%d',seq);
    if ~exist(fold_name,'dir')
        eval(['!mkdir ' fold_name])
    end
    cd(fold_name)
    
    for im_cnt = 1:ims_per_sequence
        cur_orange = [-or_width or_width]/2+orientation_sequence(im_cnt);
        [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,cur_orange,dec_const);
        noise_im(Nyp+1:end,:) = [];
        
%         env_freq_cent = 0.00;
%         env_freq_a = 0.25;
% %         [noise_env,w2,Fx,Fy] = create_bandfilt_noise(siz,env_freq_cent/niqf,env_freq_a/niqf,[-180 180]);
%         [noise_env,w2,Fx,Fy] = create_gaussfilt_noise(siz,env_freq_cent/niqf,0,env_freq_a/niqf);
%         noise_env(Nyp+1:end,:) = [];
%         gauss_cent_x = ceil(rand(n_gauss,1)*Nxp);
%         gauss_cent_y = ceil(rand(n_gauss,1)*Nyp);
%         noise_env = zeros(size(noise_im));
%         inds = sub2ind(size(noise_env),gauss_cent_y,gauss_cent_x);
%         noise_env(inds) = 1;
        
        noise_buffer_size = [1500 1500];
        noise_env = poissrnd(n_gauss/numel(noise_im),noise_buffer_size);
        noise_env = conv2(noise_env,gfilt,'same');
        noise_env(:,1:110) = []; noise_env(:,end-109:end) = [];
        noise_env(1:238,:) = []; noise_env(end-237:end,:) = [];
 
        noise_im = noise_im.*noise_env;
        
        
        bounds = prctile(noise_im(:),cbounds);
        noise_im = (noise_im-bounds(1))/range(bounds);
        noise_im = noise_im*whiteColor;
        noise_im(noise_im < blackColor) = blackColor;
        noise_im(noise_im > whiteColor) = whiteColor;
        
%         [val indx] = sort(noise_im(:));
%         imdatarescale = noise_im;
%         imdatarescale(indx) = linspace(blackColor,whiteColor,numel(noise_im));
%         noise_im = imdatarescale;
        
        noise_png = noise_im/whiteColor;
        cur_fname = sprintf('%.4d',im_cnt);
        imwrite(noise_png,strcat(cur_fname, '.png'),'png');
    end
    save orientation_sequence orientation_sequence
end



