function [rng_seed,filter_props] = generate_1d_STnoise(stim_props,num_seqs,rng_seed)
%[rng_seed,filter_props] = generate_1d_STnoise(num_seqs)
%num_seqs: number of noise sequences to generate
%temporal_cutoffs: [lcf hcf] cutoff freqs for temporal filtering (in Hz)
%spatial_cutoffs: [lcf hcf] cutoff freqs for spatial filtering (in cyc/deg)
%
%rng_seed: state of matlab rng
%filter_props: properties of the ST filter used

%% default stimulus properties
if nargin < 2
    num_seqs = 1;
end
if nargin < 3
    rng_seed = [];
end

sz = 6; %diameter of stimulus region (deg)
sub_sz = 2; %diameter of subsampled region
Nf = 400; %number of frames
or = 0; %orientation of bars
dt = 1/100; %frame duration (sec)
temporal_cutoffs = [1 4]; %[lcf hcf] for temporal freqs (Hz)
spatial_cutoffs = [1 15]; %[lcf hcf] for spatial frequencies (cyc/deg)

if ~isempty(rng_seed)
    rng_seed = rng(rng_seed);
else
    rng_seed = rng;
end

if nargin > 0 && ~isempty(stim_props)
    fn = fieldnames(stim_props);
    for ii = 1:length(fn)
        switch fn{ii}
            case 'sz'
                sz = getfield(stim_props,fn{ii});
            case 'Nf'
                Nf = getfield(stim_props,fn{ii});
            case 'or'
                or = getfield(stim_props,fn{ii});
            case 'dt'
                dt = getfield(stim_props,fn{ii});
            case 'sub_sz'
                sub_sz = getfield(stim_props,fn{ii});
            case 'temporal_cutoffs'
                temporal_cutoffs = getfield(stim_props,fn{ii});
            case 'spatial_cutoffs'
                spatial_cutoffs = getfield(stim_props,fn{ii});
            otherwise
                warning('Unusable stim-property listed');
        end
    end
end
%%
save_dir = '~/binoc_expts/gamma_BLnoise/BL_noise_seqs/';
blackColor = 0;
whiteColor = 255;
cbounds = [5 95];
grayColor = (blackColor+whiteColor)/2;
Pix2Deg = 0.018837;

niqf_x = 1/Pix2Deg/2;
niqf_t = 1/dt/2;

Nx = round(sz/Pix2Deg); %diameter in pix
Nx_sub = round(sub_sz/Pix2Deg); %diameter of subsampled region
is_circular = false;
use_radius = Nx/2;

%%
if ~is_circular && or ~= 0 && or ~= 90
    Nx_final = Nx;
    Nx = round(Nx*sqrt(2));
end

%%
[fx,ft] = freqspace([Nf Nx]);
[Fx,Ft] = meshgrid(fx,ft);
fx = fx*niqf_x; Fx = Fx*niqf_x;
ft = ft*niqf_t; Ft = Ft*niqf_t;

%make a gaussian filter to smooth the frequency-domain transfer fnx
win_smooth = 2; %smoothing sigma in units of frequency sampling (1/N)
win = fspecial('gaussian',[Nx Nf],win_smooth);

%create basic transfer fnx
X = [Fx(:) Ft(:)];
Hd = ones(Nf,Nx);
Hd(X(:,1).^2/spatial_cutoffs(1)^2 + X(:,2).^2/temporal_cutoffs(1)^2 < 1) = 0;
Hd(X(:,1).^2/spatial_cutoffs(2)^2 + X(:,2).^2/temporal_cutoffs(2)^2 > 1) = 0;

%for doing anisotropic filtering
% A = atan2(Ft,Fx);
% ori_range = ([-45 45]+45)*pi/180;
% Hdmult = (A > ori_range(1) & A < ori_range(2));
% Hd = Hd.*Hdmult;

%smooth transfer fnx
win = win./max(win(:));
Hd = conv2(Hd,win,'same');

filter_props.transfer_fun = Hd;
filter_props.fx = fx;
filter_props.ft = ft;

for nn = 1:num_seqs
    fprintf('Creating noise sequence %d of %d\n',nn,num_seqs);
    noise = randn(Nf,Nx);
    FT=Hd.*fftshift(fft2(noise));
    noise=real(ifft2(ifftshift(FT)));
    noise = noise';
    
    %saturate luminances
    bounds = prctile(noise(:),cbounds);
    %map luminances to appropriate range
    noise = (noise-bounds(1))/range(bounds);
    noise = noise*whiteColor;
    noise(noise < blackColor) = blackColor;
    noise(noise > whiteColor) = whiteColor;
    noise = noise/whiteColor;
    
    %set up constraints for circular aperture
    xax = (1:Nx); xax = xax - mean(xax);
    [XX,YY] = meshgrid(xax);
    RR = sqrt(XX.^2 + YY.^2);
    use_pix = RR < use_radius;
    
    noise_ims = repmat(noise,[1 1 Nx]);
    if or == 90
        noise_ims = permute(noise_ims,[3 1 2]);
    else
        noise_ims = permute(noise_ims,[1 3 2]);
        %do a rotation of the images of desired
        if or ~= 0
            for ii = 1:Nf
                noise_ims(:,:,ii) = imrotate(noise_ims(:,:,ii),or,'bilinear','crop');
            end
        end
    end
    
    %enforce circular aperture
    if is_circular
        noise_ims = reshape(noise_ims,[],Nf);
        use_pix = use_pix(:);
        noise_ims(~use_pix,:) = grayColor/whiteColor;
        noise_ims = reshape(noise_ims,Nx,Nx,Nf);
    elseif or ~= 0 && or ~= 90
        buff_win = ceil((Nx-Nx_final)/2);
        noise_ims = noise_ims((buff_win+1):(end-buff_win),(buff_win+1):(end-buff_win),:);
    end
        
    % figure;colormap(gray)
    % for ii = 1:Nf
    %    imagesc(squeeze(noise_ims(:,:,ii))');caxis([blackColor whiteColor]/whiteColor)
    %    pause(0.01);
    % end
    %print out sequence of pgm files
    for jj = 1:Nf
        fname = strcat(save_dir,sprintf('%.5d.pgm',nn*1000+jj));
        imwrite(noise_ims(:,:,jj),fname,'pgm');
    end
    
    %if desired, print out a second set of images which are cropped
    %versions. Save with same numbering scheme but with first number == 1
    %rather than 0.
    if ~isnan(sub_sz)
        buff_win = ceil((Nx_final - Nx_sub)/2);
        noise_ims = noise_ims((buff_win+1):(end-buff_win),(buff_win+1):(end-buff_win),:);
        for jj = 1:Nf
            fname = strcat(save_dir,sprintf('%.5d.pgm',1e4+nn*1e3+jj));
            imwrite(noise_ims(:,:,jj),fname,'pgm');
        end
    end
end

%%
% cd ~/Desktop/
% writerObj = VideoWriter('newfile.avi');
% writerObj.FrameRate = 60;
% open(writerObj);
%
% figure;colormap(gray)
% sat_val = prctile(abs(noise(:)),90);
% for ii = 1:Nf
%    imagesc(squeeze(noise_ims(:,:,ii))');caxis([-sat_val sat_val])
% frame = getframe;
% writeVideo(writerObj,frame);
% end
% close(writerObj);