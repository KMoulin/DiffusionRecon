function[P_Endo,P_Epi,LV_Mask,Mask_Depth]= ROI_KM(Dcm)

%  Draw an Endo and Epi ROI for each Slice and apply it to the DWI matrix 
%  
% SYNTAX:  [P_Endo,P_Epi,LV_mask_slc]= ROI_KM(Dcm)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices]
%
%        
% OUTPUTS:  P_Endo - List of Coordinates of the Endocardium ROI
%                 
%           P_Epi - List of Coordinates of the Epicardium ROI
%
%           LV_Mask - Mask matrix 
%                 [y x slices]
%
%           Mask_Depth - Mask of depth based on the Epi/Endo countour 
%                 [y x slices]
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    P_Endo=[];
    P_Epi=[];
    [Xq,Yq] = meshgrid(1:size(Dcm,2),1:size(Dcm,1));
    disp('Create ROI') 
    for cpt_slc=1:1:size(Dcm,3)
       [epicardium, endocardium, LV_tmp] = Spline_Segmentation(Dcm(:,:,cpt_slc,1,1,1,1), 30, 30);
       LV_Mask(:,:,cpt_slc)=LV_tmp;
       P_Endo(:,:,cpt_slc)=endocardium;
       P_Epi(:,:,cpt_slc)=epicardium;
       
       
        Endo_Line = zeros(size(endocardium));
        Epi_Line = ones(size(epicardium));
        PosRoi = cat(1,epicardium,endocardium);
        LineRoi   = cat(1,Epi_Line,Endo_Line);
        Mask_Depth(:,:,cpt_slc) = griddata(PosRoi(:,1),PosRoi(:,2),LineRoi(:,1),Xq,Yq);
       
       
    end
    %%% Create a Depth Mask
    
    
    
end

%
% Cubic Spline Cardiac LV Segmentation  
%
% SYNTAX [epicardium, endocardium] = spline_segmentation(IM, nEpi_max, nEndo_max)
%
% INPUTS: 
%   IM - 2D image matrix to be segmented
%   nEpi_max - MAXIMUM number of spline nodes allowed around epicardium
%   nEndo_max - MAXIMUM number of spline nodes allowed around endocardium
%
% OUTPUTS:
%   epicardium - nx2 matrix of points generated from epicardium spline
%   points now interpolated to 200 points (Modif KM) 
%   
%   endocardium - nx2 matrix of points generated from endocardium spline
%   points now interpolated to 200 points (Modif KM)

% 
%
% Written by Eric Aliotta and Ilya Verzhbinsky, UCLA. 07/13/2016.
% Modified by K�vin Moulin
% Ennis Lab @ UCLA; http://mrrl.ucla.edu


function [epi, endo, LV_mask] = Spline_Segmentation(IM, nEpi_max, nEndo_max)
 
clear 'hepi';
clear 'hendo';
 
figure('units','normalized','position',[0.156770833333333,0.040740740740741,0.577604166666667,0.769444444444444]); 
imagesc(IM); hold on; % here just view the image you want to base your borders on min(min(IM)) max(max(IM)) ,[0 nanmean(nanmean(IM))*5]
title(gca,['Pick up to ' int2str(nEpi_max) ' points around the border around the epicardium']);
 colormap('jet')
% caxis([0 4*nanmedian(nanmedian(IM))])
epi_tmp = zeros(nEpi_max,2);
 
for j = 1:nEpi_max
    
    h = impoint;
    epi_tmp(j,:) = getPosition(h);
    
    x = epi_tmp(1:j,1); y = epi_tmp(1:j,2);
    
    if j > 1
        t = 1:j;
        ts = 1:1/10:j;
        xs = spline(t,x,ts);
        ys = spline(t,y,ts);
        
        if exist('hepi')
            set(hepi,'Visible','off');
        end
        hepi = plot(xs,ys,'r');
    end
    
    dist = pdist([epi_tmp(j,:); epi_tmp(1,:)]);
    
    if j > 2 
        meandist_list = zeros(j, 1);
        for k = 1:j
            if k ~= j
                dist_tmp = pdist([epi_tmp(k,:); epi_tmp(k+1,:)]);
                meandist_list(k) = dist_tmp;
            end
        end
        meandist = mean(meandist_list);
        if dist < meandist
           break 
        end
    end
    
end
nEpi = j;
epi_tmp(nEpi+1:end, :) = [];

% final spline
t = 1:j+1;
ts = 1:1/10:nEpi+1;
x = epi_tmp([1:end,1],1); y = epi_tmp([1:end,1],2);
xs = spline(t,x,ts);
ys = spline(t,y,ts);
 
set(hepi,'Visible','off');
hepi = plot(xs,ys,'m');
 
epicardium = [xs;ys]';
 
endo_tmp = zeros(nEndo_max,2);
 
title(gca,['Pick up to ' int2str(nEndo_max) ' points around the border around the endocardium']);
 
for j = 1:nEndo_max
    
    h = impoint;
    endo_tmp(j,:) = getPosition(h);
    
    x = endo_tmp(1:j,1); y = endo_tmp(1:j,2);
    
    if j > 1
        t = 1:j;
        ts = 1:1/10:j;
        xs = spline(t,x,ts);
        ys = spline(t,y,ts);
        
        if exist('hendo')
            set(hendo,'Visible','off');
        end
        hendo = plot(xs,ys,'g');
    end
    
    dist = pdist([endo_tmp(j,:); endo_tmp(1,:)]);
    
    if j > 2 
        meandist_list = zeros(j, 1);
        for k = 1:j
            if k ~= j
                dist_tmp = pdist([endo_tmp(k,:); endo_tmp(k+1,:)]);
                meandist_list(k) = dist_tmp;
            end
        end
        meandist = mean(meandist_list);
        if dist < meandist
           break 
        end
    end
    
end

nEndo = j;
endo_tmp(nEndo+1:end, :) = [];

% final spline
t = 1:j+1;
ts = 1:1/10:nEndo+1;
x = endo_tmp([1:end,1],1); y = endo_tmp([1:end,1],2);
xs = spline(t,x,ts);
ys = spline(t,y,ts);
 
set(hendo,'Visible','off');
hendo = plot(xs,ys,'m');
 
endocardium = [xs;ys]';

pause;
close;
endo_mask = poly2mask(endocardium(:,1),endocardium(:,2),size(IM,1),size(IM,2));
epi_mask = poly2mask(epicardium(:,1),epicardium(:,2),size(IM,1),size(IM,2));

LV_mask = zeros(size(endo_mask));
LV_mask = LV_mask + epi_mask - (epi_mask & endo_mask);

xq=linspace(1,size(epicardium,1),200);
yq=linspace(1,size(endocardium,1),200);
epi(:,1) = interp1(epicardium(:,1),xq);
epi(:,2) = interp1(epicardium(:,2),xq);
endo(:,1) = interp1(endocardium(:,1),yq);
endo(:,2) = interp1(endocardium(:,2),yq);
return

end