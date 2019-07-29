function [Dcm2]= VPCA_KM(Dcm, enum,pca_min_ernegy)

% Aply a 24x24 PCA Filter on the DWI and nDWI matrices through
% the average dimension. 
%
% SYNTAX:  [Dcm2]= VPCA_KM(Dcm, enum,pca_min_ernegy);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%         
%           enum - Structure which contains information about the dataset 
%
%           pca_min_ernegy - Rejection limit for the PCA filter 
%                                   (usually 80%)
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages dataset]
%         
%           
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu
   
    Dcm2=[];
    disp('PCA') 
    h = waitbar(0,'PCA filter...');
    
     for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
         for cpt_b=1:1:enum.datasize(cpt_set).b     
           for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir  
               tmpDataDcm=[];                                              
               for cpt_avg=1:1:enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg
                        tmpDataDcm(:,:,cpt_avg)=Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set);                                  
               end 
               %%%%% PCA method %%%%%
%                tmpDataDcm2=tmpDataDcm;
%                for cptx=2:1:size(tmpDataDcm,1)-1
%                     for cpty=2:1:size(tmpDataDcm,2)-1
%                         tmp=[];
%                         tmp=Pca_KM(tmpDataDcm((cptx-1:cptx+1),(cpty-1:cpty+1),:),pca_min_ernegy); % Execute Pca
%                         tmpDataDcm2(cptx,cpty,:)=tmp(2,2);
%                     end
%                end



                for k=24:-1:1
                   if mod(size(tmpDataDcm,1),k)==0
                       break
                   end
                end
                divx=k;
                for k=24:-1:1
                     if mod(size(tmpDataDcm,2),k)==0
                       break
                   end
                end
                divy=k;
                stepx=size(tmpDataDcm,1)/divx;
                stepy=size(tmpDataDcm,2)/divy;
                tmpDataDcm2=zeros(size(tmpDataDcm));
                for cptx=1:1:divx
                    for cpty=1:1:divy
                        tmpDataDcm2(((cptx-1)*stepx)+1:(cptx*stepx),((cpty-1)*stepy)+1:(cpty*stepy),:)=Pca_local_KM(tmpDataDcm(((cptx-1)*stepx)+1:(cptx*stepx),((cpty-1)*stepy)+1:(cpty*stepy),:),pca_min_ernegy); % Execute Pca
                    end
                end



               for cpt_avg=1:1:enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg    
                        Dcm2(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)=tmpDataDcm2(:,:,cpt_avg,cpt_set);                     
               end
           end
         end
        waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
     end
    close(h);  

end

function [matOut] = Pca_local_KM(matIn,pca_min_ernegy)

%%%%%%%%%%%%%%%%%%%%% Sub the [n,p,dim] matrix %%%%%%%%%%%%%%%%%%%%%%

%matIn=tmpDataDcm;
dim=size(matIn,3);

matIn=double(matIn);

[n,p]=size(matIn(:,:,1));

Image=[];
Vector=[];
X=[];
for i=1:1:dim
    Image(i).dat = matIn(:,:,i);
    Vector(i).dat= Image(i).dat(:);
end


for i=1:1:dim
    Vector(i).mean=mean(Vector(i).dat);
    Vector(i).dat=Vector(i).dat-Vector(i).mean;
    X=[X Vector(i).dat];
end


%on crée la matrice ayant 3 individus et n*p réalisations
%X=[vKL_R,vKL_G,vKL_B];

% on calcule la matrice d'inertie, matrice de covariance
%V=(1/(n*p)).*(X'*X);
v=cov(X);
%on calcule les valeurs propres et les vecteurs propres
%[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% full matrix V whose columns are the corresponding eigenvectors so
% that X*V = V*D.
[mat_vect_p,D]=eig(v);

% D est une matrice avec les valeurs propres sur la diagonale,
% on souhaite les avoir dans un vecteur
vpKL=eig(D);

% 'contribution des différents axes'
% 
% Energy_KL1=[num2str((vpKL(1))/(vpKL(1)+vpKL(2)+vpKL(3))*100),' %']
% Energy_KL2=[num2str((vpKL(2))/(vpKL(1)+vpKL(2)+vpKL(3))*100),' %']
% Energy_KL3=[num2str((vpKL(3))/(vpKL(1)+vpKL(2)+vpKL(3))*100),' %']
% 
% Quality=[num2str(((vpKL(2)+vpKL(3)))/(vpKL(1)+vpKL(2)+vpKL(3))*100),' %']
maxEnergy=0;
for i=1:1:dim
    maxEnergy=maxEnergy+vpKL(i);
    Index(i)=false;
end

e=0;
i=0;
for i=1:1:dim
    Y =find(vpKL== max(vpKL));
    if e <= pca_min_ernegy
        Index(Y)=true; % Note the index of the value
        i=i+1;
        e=e+100*vpKL(Y)/maxEnergy;
        vpKL(Y)=[];
    end
end
%Quality=[num2str(e)];
%C projection de l'image sur la nouvelle base
C=X*mat_vect_p;
canal=[];
for i=1:1:dim
    if Index(i)
        canal(i).dat=C(:,i); 
    else
        canal(i).dat=zeros(n*p,1);
    end
end

%===============================================%
% RECONSTITUTION AVEC UN MAX DE VALEURS PROPRES %
%===============================================%
C_rec=[];
for i=1:1:dim
     C_rec=[C_rec canal(i).dat]; %zeros(n*p,1) if don't keep vector
end


X_rec=C_rec*mat_vect_p';

C_rec=[];

for i=1:1:dim
%     if i==1
%         C_rec=[C_rec zeros(n*p,1)];
%     else
%         C_rec=[C_rec X_rec(:,i)]; %zeros(n*p,1) if don't keep vector
%     end
    C_rec=[C_rec X_rec(:,i)]; %zeros(n*p,1) if don't keep vector
end



tmp=[];
for i=1:1:dim
    tmp2=C_rec(:,i);
    tmp(:,:)=reshape(tmp2(:),n,p);
    tmp(:,:)=tmp(:,:)+Vector(i).mean;
    matOut(:,:,i)=tmp(:,:);
end


matOut=uint16(matOut);

end



