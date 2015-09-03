function F = hsmm_ig_exp_fun(x,np,dur_range)

%for obtaining ML estimates of inverse gaussian state-duration distribution by solving equation 8 in 
% Mitchell C, Jamieson L (1993) Modeling duration in a hidden Markov model
% with the exponential family. IEEE International Conference on Acoustics,
% Speech, and Signal Processing. pp 11331–11334

F = [x(1) - np(1);
    sum(inverse_gaussian_pmf(dur_range,x(1),x(2))./dur_range) - np(2)];