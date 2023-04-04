function [matOut] = Pca_KM(matIn,pca_min_ernegy)

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
% 
% figure(1);
% subplot(2, 2, 1);
% imshow(uint8(floor(imgR)));
% 
% % Création du second axe
% subplot(2, 2, 2);
% imshow(uint8(floor(imgG)));
% 
% % Création du troisième axe
% subplot(2, 2, 3);
% imshow(uint8(floor(imgB)));
% 
% 
% figure(2);
% subplot(2, 2, 1);
% imshow(imgR_rec);
% 
% % Création du second axe
% subplot(2, 2, 2);
% imshow(imgG_rec);
% 
% % Création du troisième axe
% subplot(2, 2, 3);
% imshow(imgB_rec);
% 
% final=[];
% for i=1:n
%     for j=1:p
%         tmp=max(imgR_rec(i,j),imgG_rec(i,j));
%         tmp2=max(tmp,imgB_rec(i,j));
%         final(i,j)=tmp2 ;      
%     end
% end
% figure(3);
% imshow(uint8(floor(final)));

%matOut(:,:,1) = imgR_rec;
%matOut(:,:,2) = imgG_rec;
%matOut(:,:,3) = imgB_rec;
%difference_img=double(img)-double(img_rec);
%dif_img=uint8(abs(difference_img));



