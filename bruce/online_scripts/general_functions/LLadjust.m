function dLL = LLadjust( k, lr, Nspks, lr2, llist )
%
% Usage: dLL = LLadjust( k, lamrange, Nspks, <lamrange2>, <llist> )
%
% Includes slope-slope now (lamrange2)

if (nargin < 4) || isempty(lr2)
  Ncon2 = 0;
else
  Ncon2 = size(lr2,1);
end

Ncon = size(lr,1);

if (Ncon+Ncon2) == 0
  dLL = 0;
  return;
end

dLL = 0;

if Ncon > 0
  circ = lr(:,2) < 0;
  lr(:,2) = abs(lr(:,2));

  for i = 1:Ncon
    if lr(i,3) > 0
      dLL = dLL + lr(i,1) * sum((k(lr(i,2):(lr(i,3)-1))-(k((lr(i,2)+1):lr(i,3)))).^2);
      %disp(sprintf( '%2d  %0.6f', i, sum((k(lr(i,2):(lr(i,3)-1))-(k((lr(i,2)+1):lr(i,3)))).^2)/Nspks ))
      if circ(i) > 0
        dLL = dLL + lr(i,1) * (k(lr(i,3))-k(lr(i,2))).^2;
      end
    else
      dLL = dLL + lr(i,1) * (k(lr(i,2))-k(-lr(i,3)))^2;
    end
  end
end
if Ncon2 > 0
  for i = 1:Ncon2
      range = lr2(i,2):lr2(i,3);
      chunk = k(range)'; %'
      %chunk2 = [chunk(end) chunk chunk(1)];
	  chunk2 = [2*chunk(1)-chunk(2) chunk 2*chunk(end)-chunk(end-1)];
      dLL = dLL + lr2(i,1) * sum((2*chunk2(2:end-1) - chunk2(1:end-2) - chunk2(3:end)).^2);
      %disp(sprintf( '%2d  %0.6f', i, sum((k(lr(i,2):(lr(i,3)-1))-(k((lr(i,2)+1):lr(i,3)))).^2)/Nspks ))
  end
end
 
if nargin > 4
  if ~isempty(llist)
    dLL = dLL + llist(1) * sum(abs(k(llist(2:end))));
  end
end

dLL = dLL/Nspks;
