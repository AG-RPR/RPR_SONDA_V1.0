% Convert degrees of visual angle in pixels
%
% pix = deg2pix(deg, xres, width, view_dist)
%
% pix:          the quantity in pixels
% xres:         horizontal resolution of the screen
% width:        physical width of the screen
% view_dist:    viewing distance

function pix = deg2pix(deg, xres, width, view_dist)

pix=tan(deg2rad(deg/2)) *2*view_dist * xres/width; 

end