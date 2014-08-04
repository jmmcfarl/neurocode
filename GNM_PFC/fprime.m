function y = fprime ( x, fpr, NLx )
%
% Usage: y = fprime ( x, fpr, NLx )
%
% assigns the values in fpr to x, segregated by dividing function NLx (ends are assumed 0)

y = zeros(length(x),1);

for n = 1:length(NLx)-1
  y((x >= NLx(n)) & (x < NLx(n+1))) = fpr(n);
end

return
    
% %% 2 ways to run it -- but this actually takes a lot longer (!)
% y = x;
% 
% tic
% a = 1:length(x);
% b = find((x < NLx(1)) | (x >= NLx(end)));
% y(b) = 0;
% a = setdiff(a,b);
% 
% for n = 2:length(NLx)
%   b = find(x(a) < NLx(n));
%   y(a(b)) = fpr(n-1);
%   a = setdiff(a,a(b));
% end
% 
% toc

