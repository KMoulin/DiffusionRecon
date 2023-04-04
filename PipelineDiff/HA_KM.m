function [HA TRA]= HA_KM( EigVect1, Mask, P_Epi, P_Endo)

% Calculates the angle (degrees) between the primary eigenvector and the SA
% plane; epicardium and endocardium should be a series of points representing the
% boundaries of the myocardium.
%
% SYNTAX:  HA_KM(  EigVect1, Mask, P_Epi, P_Endo)
%
% INPUTS:   EigVect1 - First EigVector image matrix
%                 [y x slices coordinates]
%           
%           Mask -  Mask matrix
%                 [y x slices]
%
%           P_Endo - List of Coordinates of the Endocardium ROI
%                 
%           P_Epi - List of Coordinates of the Endocardium ROI
%
% OUTPUTS:  HA - HA image matrix (units [- pi pi])
%                 [y x slices]
%           
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html



%% using ellipses

yres = size(EigVect1,1);   xres = size(EigVect1,2);
[Xq,Yq] = meshgrid(1:xres,1:yres);

P_Epi=P_Epi(1:end-2,:,:);
P_Endo=P_Endo(1:end-2,:,:);

npts = size(P_Epi,1);
npts2 = size(P_Endo,1);




%%

disp('Generate HA') 
h = waitbar(0,'Generate HA...');
HA = zeros(size(Mask,1),size(Mask,2),size(Mask,3));
TRA = zeros(size(Mask,1),size(Mask,2),size(Mask,3));
for z=1:size(Mask,3)
    
    Vec = zeros(npts,2);
    Vec2 = zeros(npts2,2);
    positions=[];
    vectors=[];
    
    
     Vec(1,:) = P_Epi(1,:,z) - P_Epi(end,:,z);
    for y = 2:npts
        Vec(y,:) = P_Epi(y,:,z) - P_Epi(y-1,:,z);
    end

    Vec2(1,:) = P_Endo(1,:,z) - P_Endo(end,:,z);
    for y = 2:npts2
        Vec2(y,:) = P_Endo(y,:,z) - P_Endo(y-1,:,z);
    end
    
    positions = cat(1,P_Epi(:,:,z),P_Endo(:,:,z));
    vectors   = cat(1,Vec,Vec2);

    

    Vy = griddata(positions(:,1),positions(:,2),vectors(:,2),Xq,Yq);
    Vx = griddata(positions(:,1),positions(:,2),vectors(:,1),Xq,Yq);

    
    for y = 1:yres
        for x = 1:xres
            if Mask(y,x,z) ~= 0
                            
                
%                 if squeeze(EigVect1(y,x,z,3))>0
%                     Fiber_vect = squeeze(EigVect1(y,x,z,:));
%                 else
%                     Fiber_vect = -squeeze(EigVect1(y,x,z,:));
%                 end
                
                Fiber_vect = squeeze(EigVect1(y,x,z,:));
                vect2 = [squeeze(EigVect1(y,x,z,1)) squeeze(EigVect1(y,x,z,2)) 0];
                tang = [Vx(y,x) Vy(y,x) 0];
                tProj = dot(Fiber_vect,tang)/norm(tang)*tang;
                
                TRA(y,x,z)= acos(dot(vect2,tang)/(norm(tang)*norm(vect2)))*180/(pi);
                
                zComp = Fiber_vect(3);
                hyp = sqrt(tProj(1).^2 + tProj(2).^2+ Fiber_vect(3).^2);
                
                
               % TRA(y,x,z) = asin(zComp/hyp)*180/(pi);
                HA(y,x,z) = asin(zComp/hyp)*180/(pi);
                test(y,x,z)=dot(tang,tProj);
                if dot(tang,Fiber_vect) > 0 %&& HA(j,k) > 0  %dot(tang,tProj) > 0
                    HA(y,x,z) = -HA(y,x,z);
                end
                
                if dot(vect2,tang) < 0 %&& HA(j,k) > 0
                    TRA(y,x,z) = TRA(y,x,z)-180;
                end
            end
        end
    end
    waitbar(z/size(Mask,3),h);
end
close(h)


%HA(MASK==0) = nan;

end

