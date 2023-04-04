function Iout=affine_transform_2d_double(Iin,M,mode)
% Affine transformation function (Rotation, Translation, Resize)
% This function transforms a volume with a 3x3 transformation matrix 
%
% Iout=affine_transform_2d_double(Iin,Minv,mode)
%
% inputs,
%   Iin: The input image
%   Minv: The (inverse) 3x3 transformation matrix
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            (cubic interpolation only support by compiled mex file)
%            2: cubic interpolation and outsite pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%
% output,
%   Iout: The transformed image
%
% example,
%   % Read image
%   I=im2double(imread('lenag2.png'))
%   % Make a transformation matrix
%   M=make_transformation_matrix([2 3],2,[1.0 1.1]);
%   % Transform the image
%   Iout=affine_transform_2d_double(I,M,0)
%   % Show the image
%   figure, imshow(Iout);
%
% Function is written by D.Kroon University of Twente (February 2009)
  
% Make all x,y indices
[x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);

% Calculate center of the image
mean=size(Iin)/2;

% Make center of the image coordinates 0,0
xd=x-mean(1); 
yd=y-mean(2);

% Calculate the Transformed coordinates
Tlocalx = mean(1) + M(1,1) * xd + M(1,2) *yd + M(1,3) * 1;
Tlocaly = mean(2) + M(2,1) * xd + M(2,2) *yd + M(2,3) * 1;

switch(mode)
	case 0
		Interpolation='bilinear';
		Boundary='replicate';
	case 1
		Interpolation='bilinear';
		Boundary='zero';
	case 2
		Interpolation='bicubic';
		Boundary='replicate';
	otherwise
		Interpolation='bicubic';
		Boundary='zero';
end
Iout=image_interpolation(Iin,Tlocalx,Tlocaly,Interpolation,Boundary);
