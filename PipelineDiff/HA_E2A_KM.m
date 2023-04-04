function [HA TRA E2A RAD_s CIR_s LON_s]= HA_E2A_KM( EigVect1, EigVect2, Mask, P_Epi, P_Endo)

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

RAD_s=zeros(size(Mask,1),size(Mask,2),size(Mask,3),3);
CIR_s=zeros(size(Mask,1),size(Mask,2),size(Mask,3),3);
LON_s=zeros(size(Mask,1),size(Mask,2),size(Mask,3),3);
 FLIP_s= zeros(size(Mask,1),size(Mask,2),size(Mask,3),3);
                

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
            if Mask(y,x,z) ~= 0 & ~isnan(Mask(y,x,z))
                            
                
%                 if squeeze(EigVect1(y,x,z,3))>0
%                     Fiber_vect = squeeze(EigVect1(y,x,z,:));
%                 else
%                     Fiber_vect = -squeeze(EigVect1(y,x,z,:));
%                 end
                
                E1 = squeeze(EigVect1(y,x,z,:));
                E1= E1./norm(E1);
                E2 = squeeze(EigVect2(y,x,z,:));
                
                % Circunferiential Vector definition
                Circ = [Vx(y,x) Vy(y,x) 0];
                Circ= Circ./norm(Circ);
                % Longitudinal Vector definition
                Long= [0 0 1];
                
                % Projection of the Fiber Vector onto the Circunferential
                % direction 
                E1proj=dot(E1,Circ)*Circ/(norm(Circ)^2);
                Fiber_proj=[E1proj(1) E1proj(2) E1(3)];
                             
                Fiber_proj=Fiber_proj./norm(Fiber_proj);
                
                
                %HA(y,x,z) = atan2(Fiber_proj(3), norm([E1proj(1) E1proj(2) 0]))*180/(pi);
               HA(y,x,z) = asin(Fiber_proj(3)/norm(Fiber_proj))*180/(pi);
                
              
                % Radial Vector definition
                Rad=cross(Circ/norm(Circ),Long/norm(Long));
                
                MidFiber=cross(Fiber_proj/norm(Fiber_proj),Rad/norm(Rad));
                %MidFiber=cross(E1/norm(E1),Rad/norm(Rad));
                
                % Projection of the Sheet Vector onto the Radial
                % direction 
                E2proj=dot(E2,Rad)*Rad/(norm(Rad)^2);
                
                % Projection of the Sheet Vector onto the MidFiber
                % direction 
                tProj3=dot(E2,MidFiber)*MidFiber/(norm(MidFiber)^2);
                
                             
                E2A(y,x,z) = atan2(norm(E2proj),norm(tProj3))*180/(pi);
                
                
                
                %[vect_proj]= proj_local_KM( Fiber_vect, tang);
                %HA(y,x,z)= proj_local_KM( Fiber_vect, tang);
               vect2 = [squeeze(EigVect1(y,x,z,1)) squeeze(EigVect1(y,x,z,2)) 0];
                TRA(y,x,z)= acos(dot(vect2,Circ)/(norm(Circ)*norm(vect2)))*180/(pi);
                
                
                
               % TRA(y,x,z) = asin(zComp/hyp)*180/(pi);
                
                RAD_s(y,x,z,:)=Rad;
                CIR_s(y,x,z,:)=Circ;
                LON_s(y,x,z,:)=Long;
                FLIP_s(y,x,z,:)=Fiber_proj;
%                test(y,x,z)=dot(tang,tProj);

                %tang = [Vx(y,x) Vy(y,x) 0];
                 if dot(Fiber_proj  ,Circ) <0 % || dot(Fiber_proj,Long) > 0%dot(Circ,(E1)) < 0   %
                    HA(y,x,z) = -HA(y,x,z);
                   %  FLIP_s(y,x,z)=1;
                 end
                
%                   if abs(Rad(2))>abs(Rad(1))
%                     HA(y,x,z) = -HA(y,x,z);
%                     FLIP_s(y,x,z)=1;
%                   end
                 
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