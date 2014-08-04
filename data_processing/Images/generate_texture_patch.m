function [texture,comp_gratings] = generate_texture_patch(imSize)

Pix2Deg = 0.018837;
mean_lambda = 0.33/Pix2Deg;
std_lambda = 0.075/Pix2Deg;
min_lambda = 0.05/Pix2Deg;
texture = zeros(imSize);
n_waves = 4;
for i = 1:n_waves
    rand_lambda = randn*std_lambda + mean_lambda;
    rand_lambda = max(min_lambda,rand_lambda);
    rand_theta = 2*pi*rand;
    rand_phase = 2*pi*rand;
    A = gen_square_grating(imSize,rand_lambda,rand_theta,rand_phase);
    comp_gratings(i).lambda = rand_lambda;
    comp_gratings(i).theta = rand_theta;
    comp_gratings(i).phase = rand_phase;
    texture = texture + A;
end


% texture = texture/std(texture(:));
% texture = texture*noise_std + (white+black)/2;
% texture(texture > white) = white;
% texture(texture < black) = black;

% [X,Y] = meshgrid(0.5:imSize);
% X = X - imSize/2;
% Y = Y - imSize/2;
% circ_x = -100; 
% circ_y = 100;
% rad = 50;
% circ_mask = zeros(imSize);
% circ_mask((X-circ_x).^2+(Y-circ_y).^2 < rad.^2) = 1;

