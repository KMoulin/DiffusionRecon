function[Dcm2, Liver_Mask,Liver_ROI]= ROI_Liver_KM(Dcm,n)

%  Draw a ROI for each Slice and apply it to the DWI matrix 
%  
% SYNTAX:  [Dcm2]= ROI_Liver_KM(Dcm)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices b-values directions averages dataset]
%
%        
% OUTPUTS:  Dcm2 - image matrix 
%                 [y x slices b-values directions averages dataset]
% 
%           Liver_Mask - Mask matrix 
%                 [y x slices]
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    Dcm2=[];

    disp('Create ROI') 
    Liver_ROI=zeros(200,2,size(Dcm,3),n);
    for cpt_slc=1:1:size(Dcm,3)
        for cpt_n=1:1:n
        [Liver, Liver_tmp_Mask] = SplineSingleSegmentation_local(Dcm(:,:,cpt_slc,1,1,1,1), 80,Liver_ROI(:,:,cpt_slc,:));
        Liver_ROI(:,:,cpt_slc,cpt_n)=Liver;
        Liver_Mask(:,:,cpt_slc,cpt_n)=Liver_tmp_Mask;
        end
    end
    Liver_Mask(Liver_Mask>0)=1;
    Dcm2=Dcm.*repmat(Liver_Mask,1,1,1,size(Dcm,4),size(Dcm,5),size(Dcm,6),size(Dcm,7));
end

%
% Cubic Spline Liver  Segmentation  
%
% Syntax: [epicardium, endocardium] = spline__single_segmentation(IM, nEpi_max, nEndo_max)
%
% Inputs:
%   IM - 2D image matrix to be segmented
%   nEpi_max - MAXIMUM number of spline nodes allowed around epicardium
%   nEndo_max - MAXIMUM number of spline nodes allowed around endocardium
%
% Output:
%   epicardium - nx2 matrix of points generated from epicardium spline
%   points
%   
%   endocardium - nx2 matrix of points generated from endocardium spline
%   points

% 
%
% Written by Eric Aliotta and Ilya Verzhbinsky, UCLA. 07/13/2016.
% Modified by Kévin Moulin, UCLA. 16/10/16
% Ennis Lab @ UCLA; http://mrrl.ucla.edu


function [Liver_ROI, Liver_mask] = SplineSingleSegmentation_local(IM, nPointSpline,ROI)

 
figure; title(gca,['Pick up to ' int2str(nPointSpline) ' points around the border of the liver']);
imagesc(IM); colormap('gray');hold on; % here just view the image you want to base your borders on
for cpt_n=1:1:size(ROI,3)
 plot(ROI(:,1,cpt_n),ROI(:,2,cpt_n));
end
spline_tmp = zeros(nPointSpline,2);
 
for j = 1:nPointSpline
    
    h = impoint;
    spline_tmp(j,:) = getPosition(h);
    
    x = spline_tmp(1:j,1); y = spline_tmp(1:j,2);
    
    if j > 1
        t = 1:j;
        ts = 1:1/10:j;
        xs = spline(t,x,ts);
        ys = spline(t,y,ts);
        
        if exist('hSpline')
            set(hSpline,'Visible','off');
        end
        hSpline = plot(xs,ys,'r');
    end
    
    dist = pdist([spline_tmp(j,:); spline_tmp(1,:)]);
    
    if j > 2 
        meandist_list = zeros(j, 1);
        for k = 1:j
            if k ~= j
                dist_tmp = pdist([spline_tmp(k,:); spline_tmp(k+1,:)]);
                meandist_list(k) = dist_tmp;
            end
        end
        meandist = mean(meandist_list);
        if dist < meandist
           break 
        end
    end
    
end
nSpline = j;
spline_tmp(nSpline+1:end, :) = [];

% final spline
t = 1:j+1;
ts = 1:1/10:nSpline+1;
x = spline_tmp([1:end,1],1); y = spline_tmp([1:end,1],2);
xs = spline(t,x,ts);
ys = spline(t,y,ts);
 
set(hSpline,'Visible','off');
hSpline = plot(xs,ys,'m');
 
Liver = [xs;ys]';
 
 
pause;
close;
Liver_mask = poly2mask(Liver(:,1),Liver(:,2),size(IM,1),size(IM,2));

xq=linspace(1,size(Liver,1),200);
Liver_ROI(:,1) = interp1(Liver(:,1),xq);
Liver_ROI(:,2) = interp1(Liver(:,2),xq);

return
end