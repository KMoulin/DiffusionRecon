% IMAGE REGISTRATION MAIN
% -------------------------------------------------------------------------
close all
clear all
clc

% cd('C:\Documents and Settings\localadmin\My Documents\MATLAB\2DImageRegistrationToolbox')

% Image loading
% load 20101125_T2mappingGood
load e:\20101222_T2mappingBadZoom.mat




h_coreg = figure('name','image coregistration');

% REFERENCE IMG -------------
Ref = squeeze(bigim(2,:,:)); 



% I2 = squeeze(bigim(1,:,:));
I2 = squeeze(bigim(3,:,:));

% croppping of moving image ---------
% tmp = imcrop(I2,[xs ys ws hs]);
% I2(:) = 0;
% I2(ys:ys+hs,xs:xs+ws) = tmp; clear tmp
% ------------------------------------
% I2 = imfilter(I2,fspecial('gaussian', 5, 3));

figure('name','original images')
subplot(121),imshow(Ref,[]),title('reference image')
subplot(122),imshow(I2,[]),title('moving image')

Options = struct('Similarity',{'cc'},...
                      'Registration',{'Both'},...
                      'Penalty',{1e-3},...
                      'MaxRef',{2},...
                      'Grid',{[]},...
                      'Spacing',{[]},...
                      'MaskMoving',{[]},...
                      'MaskStatic',{[]},'Verbose',{0},...
                      'Points1',{[]},'Points2',{[]},'PStrength',{[]},...
                      'Interpolation',{'Linear'},'Scaling',{[1 1]});
                                   
% Similarity Options
% 'd'  : differences between I1 and I2
% 'sd' : squared differences
% 'mi' : normalized (local) mutual information
% 'mip': normalized mutual information with image split in multiple small regions
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
% Registration Options:
%               Rigid    : Translation, Rotation
%               Affine   : Translation, Rotation, Shear, Resize
%               NonRigid : B-spline grid based local registration
%               Both     : Nonrigid and Affine (Default)


[Coregged,O_trans,Spacing,M,param,B,F] = image_registration (I2,Ref,Options);
% M: affine matrix of the registration
% param: parameters of the registration
%        2D - [ deltaY deltaX theta ]
% -------------------------------------------------------------------------

error = abs(Coregged-I2);
% -------------plot-------------------------------------------
figure(h_coreg),drawnow
subplot(221),imshow(Ref,[]),title('reference image')
subplot(222),imshow(I2,[]),title('moving image')
subplot(223),imshow(Coregged,[]),title('coregistered')
subplot(224),imshow(error,[]),title('error = Coreg - Imov')
% -------------plot-------------------------------------------

coregYdisp = param(1);
coregXdisp = param(2);
coregRot = param(3);



