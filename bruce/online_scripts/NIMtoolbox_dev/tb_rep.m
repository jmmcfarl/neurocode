function Xout = tb_rep(gin, NLx)
%
% Xout = tb_rep(gin, NLx)
%
% Creates a TxNbfs matrix representing the input vector gin passed through a
% series of tent basis functions specified by edge points in NLx

%%
n_tbs =length(NLx); %number of tent-basis functions
Xout = zeros(length(gin),n_tbs);
Xout(:,1) = get_tentbasis_output(gin,NLx(1),[-Inf NLx(2)]);
Xout(:,end) = get_tentbasis_output(gin,NLx(end),[NLx(end-1) Inf]);
for n = 2:n_tbs-1
    Xout(:,n) = get_tentbasis_output(gin, NLx(n), [NLx(n-1) NLx(n+1)] );
end
