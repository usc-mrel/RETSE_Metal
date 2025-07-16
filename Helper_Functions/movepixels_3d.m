function Iout = movepixels_3d(Iin,Tx)
% This function will translate the pixels of an image
%  according to Tx translation images. 
%
% Inputs;
%   Tx, Ty: The transformation images, describing the
%             (backwards) translation of every pixel in x and y direction.
%
% Outputs,
%   Iout : The transformed image
%
  
% Make all x,y,z indices
[x,y,z] = meshgrid(0:size(Iin,2)-1,0:size(Iin,1)-1,0:size(Iin,3)-1);

% Calculate the Transformed coordinates. Note this uses MATLAB's
% convention of y being rows, x being columns. The framework 
% expects the distorted dimension for dual polarity MRI to be in 
% the row dimension. 
Tlocalx = x;
Tlocaly = y + Tx;
Tlocalz = z;

Iout = interp3(x,y,z,Iin,Tlocalx,Tlocaly,Tlocalz,'linear',0);
