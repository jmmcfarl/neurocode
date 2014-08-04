close all
clear all

save_dir = '~/Data/bruce/Expt_1_9_13_imfolder';
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

freq_cent = 4.5; %cyc/deg
freq_a = 9;
or_width = 35;
dec_const = 2/niqf;

% fore_cent = [-2 0];
fore_width = [3 3];

%%
n_images = 500;
poss_or_cents = -60:30:90;
n_poss_or_cents = length(poss_or_cents);
x0_avg = 0.35;
x0_range = 4.5;
y0_range = 4.5;
y0_avg = 0.43;

DD = sqrt((XAX-x0_avg).^2+(YAX-y0_avg).^2);
[~,RF_loc_ind] = min(DD(:));

clear is_rf_fore

right_edge = [x0_avg + x0_range/fore_width(1) y0_avg];
left_edge = [x0_avg - x0_range/fore_width(1) y0_avg];
% top_edge = [x0_avg y0_avg - y0_range/fore_width(2)];
% bottom_edge = [x0_avg y0_avg + y0_range/fore_width(2)];
off_right = [x0_avg + x0_range/fore_width(1)+1 y0_avg];
off_left = [x0_avg - x0_range/fore_width(1)-1 y0_avg];
% edge_locs_x = [right_edge(1) left_edge(1) top_edge(1) bottom_edge(1)];
% edge_locs_y = [right_edge(2) left_edge(2) top_edge(2) bottom_edge(2)];
edge_locs_x = [right_edge(1) left_edge(1) off_right(1) off_left(1)];
edge_locs_y = [right_edge(2) left_edge(2) off_right(2) off_left(2)];
edge_frac = 1/3;

null_frac = 1/4;

rand_seq = ceil(rand(n_images,1)*n_poss_or_cents);
fore_orientation_sequence = poss_or_cents(rand_seq);
back_orientation_sequence = fore_orientation_sequence + 90;
back_orientation_sequence(back_orientation_sequence > 90) = mod(back_orientation_sequence(back_orientation_sequence > 90),90)-90;
back_orientation_sequence(fore_orientation_sequence==90) = 0;

% x0_sequence = (rand(n_images,1)-0.5)*x0_range +x0_avg;
% y0_sequence = (rand(n_images,1)-0.5)*y0_range +y0_avg;
x0_sequence = ones(n_images,1)*x0_avg;
y0_sequence = ones(n_images,1)*y0_avg;

edge_sequence = ceil(4*rand(n_images,1));
use_edge = rand(n_images,1) < edge_frac;

use_null = rand(n_images,1) < null_frac; %no foreground

x0_sequence(use_edge) = edge_locs_x(edge_sequence(use_edge));
y0_sequence(use_edge) = edge_locs_y(edge_sequence(use_edge));
x0_sequence(use_null) = nan;
y0_sequence(use_null) = nan;


for im_cnt = 1:n_images;
    fprintf('Generating image %d of %d\n',im_cnt,n_images);
    
    or_cent = back_orientation_sequence(im_cnt);
%     [noise2,w,Fx,Fy] = create_bandfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent);
    [noise2,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent,dec_const);
    noise2(Nyp+1:end,:) = [];

    if ~isnan(x0_sequence(im_cnt))
    or_cent = fore_orientation_sequence(im_cnt);
%     [noise,w,Fx,Fy] = create_bandfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent);
    [noise,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent,dec_const);
    noise(Nyp+1:end,:) = [];
      
    foreground = ones(size(noise));
    d = XAX > x0_sequence(im_cnt) + fore_width(1)/2 | XAX < x0_sequence(im_cnt)-fore_width(1)/2 | ...
        YAX > y0_sequence(im_cnt) + fore_width(2)/2 | YAX < y0_sequence(im_cnt) - fore_width(2)/2;
    foreground(d) = 0;
    
    noise_im = noise2;
    noise_im(foreground==1) = noise(foreground==1);
    else
        noise_im = noise2;
    end
    bounds = prctile(noise_im(:),cbounds);
    noise_im = (noise_im-bounds(1))/range(bounds);
    noise_im = noise_im*whiteColor;
    noise_im(noise_im < blackColor) = blackColor;
    noise_im(noise_im > whiteColor) = whiteColor;
    
    noise_png = noise_im/whiteColor;
    
    im_id = 3000000 + im_cnt;
    cur_fname = sprintf('IM%.7d',im_id);
    imwrite(noise_png,strcat(save_dir, '/', cur_fname, '.pgm'),'pgm');
    
%     is_rf_fore(im_cnt) = foreground(RF_loc_ind);
end
cd ~/Data/bruce/
save foreground_data_fin x0_sequence y0_sequence fore_orientation_sequence back_orientation_sequence

%% ISOTROPIC FIGURE GROUND
freq_cent1 = 3.5; %cyc/deg
freq_a1 = 5;
freq_cent2 = 6.5; %cyc/deg
freq_a2 = 5;

n_images = 500;
% poss_or_cents = -60:30:90;
% n_poss_or_cents = length(poss_or_cents);
x0_avg = 0.35;
x0_range = 4.5;
y0_range = 4.5;
y0_avg = 0.43;

DD = sqrt((XAX-x0_avg).^2+(YAX-y0_avg).^2);
[~,RF_loc_ind] = min(DD(:));

clear is_rf_fore

% right_edge = [x0_avg + x0_range/fore_width(1) y0_avg];
% left_edge = [x0_avg - x0_range/fore_width(1) y0_avg];
% top_edge = [x0_avg y0_avg - y0_range/fore_width(2)];
% bottom_edge = [x0_avg y0_avg + y0_range/fore_width(2)];
% edge_locs_x = [right_edge(1) left_edge(1) top_edge(1) bottom_edge(1)];
% edge_locs_y = [right_edge(2) left_edge(2) top_edge(2) bottom_edge(2)];
% edge_frac = 0.3;

fore_freq_sequence = ceil(rand(n_images,1)*2);

% x0_sequence = (rand(n_images,1)-0.5)*x0_range +x0_avg;
% y0_sequence = (rand(n_images,1)-0.5)*y0_range +y0_avg;
x0_sequence = ones(n_images,1)*x0_avg;
y0_sequence = ones(n_images,1)*y0_avg;

edge_sequence = ceil(4*rand(n_images,1));
use_edge = rand(n_images,1) < edge_frac;

use_null = rand(n_images,1) < null_frac; %no foreground


x0_sequence(use_edge) = edge_locs_x(edge_sequence(use_edge));
y0_sequence(use_edge) = edge_locs_y(edge_sequence(use_edge));
x0_sequence(use_null) = nan;
y0_sequence(use_null) = nan;

for im_cnt = 1:n_images;
    fprintf('Generating image %d of %d\n',im_cnt,n_images);
    
    if fore_freq_sequence(im_cnt) == 1
        [noise,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent1/niqf,freq_a1/niqf,[-180 180],dec_const);
        noise(Nyp+1:end,:) = [];
        
        [noise2,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent2/niqf,freq_a2/niqf,[-180 180],dec_const);
        noise2(Nyp+1:end,:) = [];
    else
        [noise,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent2/niqf,freq_a2/niqf,[-180 180],dec_const);
        noise(Nyp+1:end,:) = [];
        
        [noise2,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent1/niqf,freq_a1/niqf,[-180 180],dec_const);
        noise2(Nyp+1:end,:) = [];
    end
    
    if ~isnan(x0_sequence(im_cnt))
    foreground = ones(size(noise));
    d = XAX > x0_sequence(im_cnt) + fore_width(1)/2 | XAX < x0_sequence(im_cnt)-fore_width(1)/2 | ...
        YAX > y0_sequence(im_cnt) + fore_width(2)/2 | YAX < y0_sequence(im_cnt) - fore_width(2)/2;
    foreground(d) = 0;
    
    noise_im = noise2;
    noise_im(foreground==1) = noise(foreground==1);
    else
        noise_im = noise2;
    end
    bounds = prctile(noise_im(:),cbounds);
    noise_im = (noise_im-bounds(1))/range(bounds);
    noise_im = noise_im*whiteColor;
    noise_im(noise_im < blackColor) = blackColor;
    noise_im(noise_im > whiteColor) = whiteColor;
    
    noise_png = noise_im/whiteColor;
    
    im_id = 3000000 + 1000+ im_cnt;
    cur_fname = sprintf('IM%.7d',im_id);
    imwrite(noise_png,strcat(save_dir, '/', cur_fname, '.pgm'),'pgm');
    
%     is_rf_fore(im_cnt) = foreground(RF_loc_ind);
end
cd ~/Data/bruce/
save foreground_data_freq_fin x0_sequence y0_sequence fore_freq_sequence

%% RAPIDLY FLASHED FIGURE GROUND
% N_seqs = 150;
% frames_per_set = 20;
% sets_per_trial = 1;
% 
% poss_or_cents = -60:30:90;
% n_poss_or_cents = length(poss_or_cents);
% x0_avg = 0.35;
% x0_range = 4.5;
% y0_range = 4.5;
% y0_avg = 0.43;
% 
% DD = sqrt((XAX-x0_avg).^2+(YAX-y0_avg).^2);
% [~,RF_loc_ind] = min(DD(:));
% 
% clear is_rf_fore
% 
% right_edge = [x0_avg + x0_range/fore_width(1) y0_avg];
% left_edge = [x0_avg - x0_range/fore_width(1) y0_avg];
% top_edge = [x0_avg y0_avg - y0_range/fore_width(2)];
% bottom_edge = [x0_avg y0_avg + y0_range/fore_width(2)];
% edge_locs_x = [right_edge(1) left_edge(1) top_edge(1) bottom_edge(1)];
% edge_locs_y = [right_edge(2) left_edge(2) top_edge(2) bottom_edge(2)];
% edge_frac = 0.3;
% 
% 
% rand_seq = ceil(rand(N_seqs,sets_per_trial)*n_poss_or_cents);
% fore_orientation_sequence = poss_or_cents(rand_seq);
% back_orientation_sequence = fore_orientation_sequence + 90;
% back_orientation_sequence(back_orientation_sequence > 90) = mod(back_orientation_sequence(back_orientation_sequence > 90),90)-90;
% back_orientation_sequence(fore_orientation_sequence==90) = 0;
% fore_orientation_sequence = fore_orientation_sequence';
% back_orientation_sequence = back_orientation_sequence';
% 
% x0_sequence = (rand(N_seqs,sets_per_trial)-0.5)*x0_range +x0_avg;
% y0_sequence = (rand(N_seqs,sets_per_trial)-0.5)*y0_range +y0_avg;
% 
% edge_sequence = ceil(4*rand(N_seqs,sets_per_trial));
% use_edge = rand(N_seqs,sets_per_trial) < edge_frac;
% x0_sequence(use_edge) = edge_locs_x(edge_sequence(use_edge));
% y0_sequence(use_edge) = edge_locs_y(edge_sequence(use_edge));
% 
% 
% for cur_seq = 1:N_seqs;
%     fprintf('Generating sequence %d of %d\n',cur_seq,N_seqs);
%     for cur_set = 1:sets_per_trial
%         for im = 1:frames_per_set
%             or_cent = fore_orientation_sequence(cur_seq,cur_set);
%             [noise,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent,dec_const);
%             noise(Nyp+1:end,:) = [];
%             
%             or_cent = back_orientation_sequence(cur_seq,cur_set);
%             [noise2,w,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent,dec_const);
%             noise2(Nyp+1:end,:) = [];
%             
%             foreground = ones(size(noise));
%             d = XAX > x0_sequence(cur_seq,cur_set) + fore_width(1)/2 | XAX < x0_sequence(cur_seq,cur_set)-fore_width(1)/2 | ...
%                 YAX > y0_sequence(cur_seq,cur_set) + fore_width(2)/2 | YAX < y0_sequence(cur_seq,cur_set) - fore_width(2)/2;
%             foreground(d) = 0;
%             
%             noise_im = noise2;
%             noise_im(foreground==1) = noise(foreground==1);
%             
%             bounds = prctile(noise_im(:),cbounds);
%             noise_im = (noise_im-bounds(1))/range(bounds);
%             noise_im = noise_im*whiteColor;
%             noise_im(noise_im < blackColor) = blackColor;
%             noise_im(noise_im > whiteColor) = whiteColor;
%             
%             noise_png = noise_im/whiteColor;
%             
%             im_id = 5000000 + 1000*cur_seq + (cur_set-1)*frames_per_set + im;
%             cur_fname = sprintf('IM%.7d',im_id);
%             imwrite(noise_png,strcat(save_dir, '/', cur_fname, '.png'),'png');
%         end
%     end
% end
% 
% cd ~/Data/bruce/
% save foreground_data_flashed x0_sequence y0_sequence fore_orientation_sequence back_orientation_sequence
% 
