function newimg = DownSampleImg(img0, fracx, fracy)
% Function: DownSample 2D image
% Input:
%       img0(NX,NY): original image
%       fracx, fracy: down sample fraction along each dimension
% Output:
%       newimg(NX/fracx, NY/fracy)
% Yuwei Cui, Created by Oct 1, 2011

[NX NY] = size(img0);
if fracx==1 && fracy==1
    newimg = img0;
    return;
end
indx = 1:fracx:NX+1;
indy = 1:fracy:NY+1;

newimg = zeros(length(indx)-1, length(indy)-1);

for x = 1:length(indx)-1
    for y = 1:length(indy)-1
        newimg(x,y) = mean(mean((img0(indx(x):indx(x+1)-1,indy(y):indy(y+1)-1))));
    end
end