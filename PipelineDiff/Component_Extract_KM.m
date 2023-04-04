function [Dx, Dy, Dz]= Component_Extract_KM(Sx, Sy, Sz,b, Mat)
 
    %Mat = [abs(abs(Vx)-abs(Vy))';abs(abs(Vx)-abs(Vz))';abs(abs(Vy)-abs(Vz))'];
    
    
     SxV = reshape(Sx, 1, []);
     SyV = reshape(Sy, 1, []);
     SzV = reshape(Sz, 1, []);
     
     S = abs([log((SxV./SyV));log((SxV./SzV));log((SyV./SzV))]); 
     
     SS=(Mat)\S;
     
     Dx = reshape((SS(1,:)), size(Sx,1), size(Sx,2))/-b;
     
     Dy = reshape((SS(2,:)), size(Sy,1), size(Sy,2))/-b;
     
     Dz = reshape((SS(3,:)), size(Sz,1), size(Sz,2))/-b;
     
     Dx(abs(Dx)>1)=1;
     Dy(abs(Dy)>1)=1;
     Dz(abs(Dz)>1)=1;
     
     Dx(isnan(Dx))=0;
     
     Dy(isnan(Dy))=0;
     
     Dz(isnan(Dz))=0;
%     Sxy=log(Sx./Sy);
%     Sxz=log(Sx./Sz);
%     Syz=log(Sy./Sz);
    

end