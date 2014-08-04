function F = hsmm_gamma_exp_fun(x,np)

%for obtaining ML estimates of gamma state-duration distribution by solving equation 8 in 
% Mitchell C, Jamieson L (1993) Modeling duration in a hidden Markov model
% with the exponential family. IEEE International Conference on Acoustics,
% Speech, and Signal Processing. pp 11331–11334

F = [x(1)/x(2) - np(1);
    psi(x(1)) - log(x(2)) - np(2)];