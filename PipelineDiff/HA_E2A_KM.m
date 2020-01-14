function [HA TRA E2A]= HA_E2A_KM( EigVect1, EigVect2, Mask, P_Epi, P_Endo)

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

P_Epi=P_Epi(1:end-2,:,:);
P_Endo=P_Endo(1:end-2,:,:);

npts = size(P_Epi,1);
npts2 = size(P_Endo,1);




%%

disp('Generate HA') 
h = waitbar(0,'Generate HA...');
E2A = zeros(size(Mask,1),size(Mask,2),size(Mask,3));
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
                Sheet_vect = squeeze(EigVect2(y,x,z,:));
                
                % Circunferiential Vector definition
                Circ = [Vx(y,x) Vy(y,x) 0];
                
                % Longitudinal Vector definition
                Long= [0 0 1];
                
                % Projection of the Fiber Vector onto the Circunferential
                % direction 
                tProj=dot(Fiber_vect,Circ)*Circ/(norm(Circ)^2);
                Fiber_proj=[tProj(1) tProj(2) Fiber_vect(3)];
                             
                Fiber_proj=Fiber_proj./norm(Fiber_proj);
                
                HA(y,x,z) = asin(Fiber_proj(3)/norm(Fiber_proj))*180/(pi);
                
              
                % Radial Vector definition
                Rad=cross(Circ/norm(Circ),Long/norm(Long));
                
                MidFiber=cross(Fiber_proj/norm(Fiber_proj),Rad/norm(Rad));
                
                
                % Projection of the Sheet Vector onto the Radial
                % direction 
                tProj2=dot(Sheet_vect,Rad)*Rad/(norm(Rad)^2);
                
                % Projection of the Sheet Vector onto the MidFiber
                % direction 
                tProj3=dot(Sheet_vect,MidFiber)*MidFiber/(norm(MidFiber)^2);
                
                             
                E2A(y,x,z) = atan2(norm(tProj2),norm(tProj3))*180/(pi);
                
                
                
                %[vect_proj]= proj_local_KM( Fiber_vect, tang);
                %HA(y,x,z)= proj_local_KM( Fiber_vect, tang);
               vect2 = [squeeze(EigVect1(y,x,z,1)) squeeze(EigVect1(y,x,z,2)) 0];
                TRA(y,x,z)= acos(dot(vect2,Circ)/(norm(Circ)*norm(vect2)))*180/(pi);
                
                
                
               % TRA(y,x,z) = asin(zComp/hyp)*180/(pi);
                
                
                
                
                
%                test(y,x,z)=dot(tang,tProj);
                if dot(Circ,Fiber_vect) > 0 %&& HA(j,k) > 0  %dot(tang,tProj) > 0
                    HA(y,x,z) = -HA(y,x,z);
                end
                
                if dot(vect2,Circ) < 0 %&& HA(j,k) > 0
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


function [HA]= proj_local_KM( vect, plan)

         
         Xscale = dot(vect,plan)/norm(plan); %%%*plan;       
         vect_proj=[Xscale vect(3)];
         HA=asin(Xscale/norm(vect_proj))*180/(pi);

end