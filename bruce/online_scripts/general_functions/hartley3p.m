function [ipat] = hartley3p(tori,tsfreq,tphase,npix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

refcs = linspace(-pi,pi,npix); 
pixcs = zeros(npix^2,2); 
for ix=1:npix; 
	for iy=1:npix;
        pixcs((ix-1)*npix+iy,:) = [refcs(ix),refcs(iy)]; 
	end
end

%jmm: i modified this ad hoc just to get agreement with the STA analysis on
%other stimuli...  Pretty convinced this is now correct
% ang    =  -pi * tori / 180; 
ang    =  pi * (tori-90) / 180; 
rotmat = [[cos(ang),-sin(ang)];[sin(ang),cos(ang)]];  
% rotmat = [[cos(ang),sin(ang)];[-sin(ang),cos(ang)]];  
rotcs  = pixcs *rotmat;
% ipat   = reshape( sin(tsfreq.* rotcs(:,2) + pi*tphase/180),npix,npix); 
ipat   = reshape( sin(tsfreq.* rotcs(:,2) + pi*tphase/180),npix,npix); 

