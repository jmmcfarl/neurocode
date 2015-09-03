function F = hsmm_gamma_exp_fun(x,np)

F = [x(1)/x(2) - np(1);
    psi(x(1)) - log(x(2)) - np(2)];