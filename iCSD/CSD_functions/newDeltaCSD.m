function [CSD] = newDeltaCSD( data, electrodePos, dt, b0, b1, exCond, topCond, diam )
electrodePos = electrodePos * 1E-3;
switch(nargin)
    case 0
        error('Too few input parameters')
    case 1
        % TODO
        % set default values
end


CSD = F_delta(electrodePos,diam,exCond,topCond)^(-1) * data;
if(b1 ~= 0)
    [n1,n2]=size(CSD);
    CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
    CSD_add(n1+2,:)=zeros(1,n2);
    CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
    CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
end

% plot_CSD(CSD,electrodePos,dt,1,0);


end