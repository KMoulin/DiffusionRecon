function [Dcm2 Combinaison,ErrorMap]= Error_Map_KM(Dcm, Mask, enum)
 
    Dcm2=[];
    Combinaison=[];
    disp('Error Map') 
   
    for cpt_slc=5:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           k=1;
           ErrorMap=zeros(size(Dcm,1),size(Dcm,2),enum.dataset(cpt_b).dirNum,enum.dataset(cpt_b).dir(1).avgNum);
            h = waitbar(0,'Error Map...');
           for cpt_dir=1:1:(enum.dataset(cpt_b).dirNum-2);          
               tic
               for cpt_dir2=(cpt_dir+1):1:(enum.dataset(cpt_b).dirNum-1); 
                   for cpt_dir3=(cpt_dir2+1):1:enum.dataset(cpt_b).dirNum;  
                       
                        Vx=enum.dataset(cpt_b).dirVector(:,cpt_dir); 
                        Vy=enum.dataset(cpt_b).dirVector(:,cpt_dir2);
                        Vz=enum.dataset(cpt_b).dirVector(:,cpt_dir3);
                       
                         Mat = [abs(abs(Vx)-abs(Vy))';abs(abs(Vx)-abs(Vz))';abs(abs(Vy)-abs(Vz))'];
                         if(cpt_b~=1 &&cpt_dir~=cpt_dir2 && cpt_dir2~=cpt_dir3 && cpt_dir~=cpt_dir3 && mean(mean(isinf(inv(Mat))))~=1 )
                             
                             for cpt_avg=1:1:1%(enum.dataset(cpt_b).dir(cpt_dir).avgNum)
                                for cpt_avg2=1:1:1%(enum.dataset(cpt_b).dir(cpt_dir).avgNum)    
                                    for cpt_avg3=1:1:1%(enum.dataset(cpt_b).dir(cpt_dir).avgNum)
                                       
                                         Dx=[];
                                         Dy=[];
                                         Dz=[];
                                         [Dx,Dy,Dz]=Component_Extract_KM(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,cpt_avg),Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir2,cpt_avg2),Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir3,cpt_avg3),enum.b(cpt_b),Mat);
                                         
                                         Dcm2(:,:,k,1)=Dx;
                                         Dcm2(:,:,k,2)=Dy;
                                         Dcm2(:,:,k,3)=Dz;
                                         Combinaison(k,1,:)=[cpt_dir  cpt_avg];
                                         Combinaison(k,2,:)=[cpt_dir2 cpt_avg2];
                                         Combinaison(k,3,:)=[cpt_dir3 cpt_avg3];
                                            %+(abs(Dx)+abs(Dx)+abs(Dx))/3;%
                                         ErrorMap(:,:,cpt_dir,cpt_avg)  = ErrorMap(:,:,cpt_dir,cpt_avg)  +(abs(Dx.*enum.dataset(cpt_b).dirVector(1,cpt_dir))  + abs(Dy.*enum.dataset(cpt_b).dirVector(2,cpt_dir))  + abs(Dz.*enum.dataset(cpt_b).dirVector(3,cpt_dir)) )  ;
                                         ErrorMap(:,:,cpt_dir2,cpt_avg2)= ErrorMap(:,:,cpt_dir2,cpt_avg2)+(abs(Dx.*enum.dataset(cpt_b).dirVector(1,cpt_dir2)) + abs(Dy.*enum.dataset(cpt_b).dirVector(2,cpt_dir2)) + abs(Dz.*enum.dataset(cpt_b).dirVector(3,cpt_dir2)))  ;
                                         ErrorMap(:,:,cpt_dir3,cpt_avg3)= ErrorMap(:,:,cpt_dir3,cpt_avg3)+(abs(Dx.*enum.dataset(cpt_b).dirVector(1,cpt_dir3)) + abs(Dy.*enum.dataset(cpt_b).dirVector(2,cpt_dir3)) + abs(Dz.*enum.dataset(cpt_b).dirVector(3,cpt_dir3)))  ;
                                        
                                         k=k+1;
                                   end
                               end
                             end
                         end                    
                   end  
               end
               toc
               waitbar(cpt_dir/(enum.dataset(cpt_b).dirNum-2),h); 
            end
           close(h);
       end
       
    end
    
    DxMin=mean(Dcm2(:,:,:,1),3);
    DyMin=mean(Dcm2(:,:,:,2),3);
    DzMin=mean(Dcm2(:,:,:,3),3);

     for cpt_slc=5:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)      
           for cpt_dir=1:1:(enum.dataset(cpt_b).dirNum);                              
                 for cpt_avg=1:1:1%(enum.dataset(cpt_b).dir(cpt_dir).avgNum)
                   if(cpt_b~=1)  
                    num=size(find(((Combinaison(:,:,1)==cpt_dir).*(Combinaison(:,:,2)==cpt_avg))==1));
                    Dmin=(abs(DxMin.*enum.dataset(cpt_b).dirVector(1,cpt_dir))/3  + abs(DyMin.*enum.dataset(cpt_b).dirVector(2,cpt_dir))/3  + abs(DzMin.*enum.dataset(cpt_b).dirVector(3,cpt_dir))/3);
                    ErrorMap(:,:,cpt_dir,cpt_avg)=(ErrorMap(:,:,cpt_dir,cpt_avg)./(num(1))) ;
                   end  
                 end
           end
       end
     end
     
    

end