function [HA TRA]= HA_KM2( EigVect1, Mask, P_Epi, P_Endo)

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
% ???? 08.14.2017
% ?????
% Ennis Lab @ UCLA; http://mrrl.ucla.edu





%% using ellipses

yres = size(EigVect1,1);   xres = size(EigVect1,2);
[Xq,Yq] = meshgrid(1:xres,1:yres);

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
                
                
                if(squeeze(EigVect1(y,x,z,3)))>0 || (squeeze(EigVect1(y,x,z,1))<0 & squeeze(EigVect1(y,x,z,2))<0 &squeeze(EigVect1(y,x,z,3))<0)
                    Fiber_vect = squeeze(EigVect1(y,x,z,:));
                else
                    Fiber_vect = -squeeze(EigVect1(y,x,z,:));
                end
               % Fiber_vect = squeeze(EigVect1(y,x,z,:));
                
                SA_vect=[Fiber_vect(1) Fiber_vect(2) 0];
                EPI_vect = [Vx(y,x) Vy(y,x) 0];  
                
                HA(y,x,z)= acos( dot(Fiber_vect,SA_vect)/(norm(Fiber_vect)*norm(SA_vect)) )*180/pi;
                TRA(y,x,z)= acos( dot(Fiber_vect,EPI_vect)/(norm(Fiber_vect)*norm(EPI_vect)) )*180/pi;

                Norm_vect=[ 0 0 1];
                
                Proj_vect = dot(Fiber_vect,EPI_vect)/norm(EPI_vect)*EPI_vect;
                test(y,x)=dot(EPI_vect,Proj_vect);
                test2(y,x)=dot(Fiber_vect,EPI_vect);
                %if dot(EPI_vect,Proj_vect) > 0
                if dot(Fiber_vect,EPI_vect) > 0 %&& HA(j,k) > 0
                    HA(y,x,z) = -abs(HA(y,x,z));
                else
                     HA(y,x,z) = abs(HA(y,x,z));
                end
                
                if dot(Fiber_vect,EPI_vect) < 0 %&& HA(j,k) > 0
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

