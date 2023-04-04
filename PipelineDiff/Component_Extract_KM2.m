function [Dx, Dy, Dz]= Component_Extract_KM2(S1, S2, S3, S4, S5, S6, b, Mat)
 
    %Mat = [abs(abs(Vx)-abs(Vy))';abs(abs(Vx)-abs(Vz))';abs(abs(Vy)-abs(Vz))'];
    
    
     S1V = reshape(S1, 1, []);
     S2V = reshape(S2, 1, []);
     S3V = reshape(S3, 1, []);
     S4V = reshape(S4, 1, []);
     S5V = reshape(S5, 1, []);
     S6V = reshape(S6, 1, []);
     
     S = abs([log(S1V./S2V);log(S3V./S4V);log(S5V./S6V)]); 
     
     SS=Mat\S;
     
     Dx = reshape(abs(SS(1,:)), size(S1,1), size(S1,2))*3/b;
     
     Dy = reshape(abs(SS(2,:)), size(S1,1), size(S1,2))*3/b;
     
     Dz = reshape(abs(SS(3,:)), size(S1,1), size(S1,2))*3/b;
     
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