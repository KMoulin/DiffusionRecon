function [Tensor, EigValue,EigVector,MD,FA,Trace_DTI]= Calc_Tensor_KM(varargin)

% Generate a tensor from the DWI data for each slice each b-values
% the first b-value is used as the nDWI reference for the tensor
% calculation
%
% SYNTAX:  [Tensor, EigValue,EigVector,MD,FA,Trace_DTI]= Calc_Tensor_KM(Dcm, enum)
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%
%          
% OUTPUTS:  Tensor - Tensor image matrix 
%                 [y x slices b-values [xx xy xz; yx yy yz; zx zy zz]]
%
%           EigValue - Ordered EigValue image matrix 
%                 [y x slices b-values [EigV1 EigV2 EigV3]]
%           
%           EigValue - Ordered EigVector image matrix 
%                 [y x slices b-values [x y z][EigVect1 EigVect2 EigVect3]]
%
%           MD - Mean Diffusivity from Tensor
%                 [y x slices b-values]
%
%           FA - Fraction of anysotropy from Tensor
%                 [y x slices b-values]
%
%           Trace_DTI - Trace image from Tensor
%                 [y x slices b-values]
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    Tensor=[];
    EigValue=[];
    EigVector=[];
    MD=[];
    FA=[];
    Trace_DTI=[];
    
    tmpInput=[];  
    
      if numel(varargin) == 2   
          Dcm=varargin{1};
          enum=varargin{2}(1);
          R=[1 0 0; 0 1 0; 0 0 1];
          
      elseif numel(varargin) == 3
         Dcm=varargin{1};
          enum=varargin{2}(1);
          R=varargin{3};
     else
        
      end
      
    disp('Tensor calculation')
    h = waitbar(0,'Tensor calculation...');
      for cpt_set=1:1:enum.nset
            for cpt_slc=1:1:enum.datasize(cpt_set).slc
                 for cpt_b=2:1:enum.datasize(cpt_set).b     
                          
 
                          tmpInput(:,:,1)=squeeze(Dcm(:,:,cpt_slc,1,1,1,cpt_set));
                          tmpInput(:,:,(2:enum.datasize.dir+1))=squeeze(Dcm(:,:,cpt_slc,cpt_b,:,1,cpt_set));
                          [tmp_Tensor,tmp_EigVector,tmp_EigValue, Mat] = Tensor_local(tmpInput,enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dirVector'*R,enum.b(cpt_b));
                          
                          Tensor(:,:,cpt_slc,cpt_b-1,:)=tmp_Tensor;
                          EigVector(:,:,cpt_slc,cpt_b-1,:,:)=tmp_EigVector;
                          EigValue(:,:,cpt_slc,cpt_b-1,:)=tmp_EigValue;
                          [tmp_MD, tmp_FA, tmp_Trace] = Maps_local(tmp_EigValue); 
                          MD(:,:,cpt_slc,cpt_b-1)=tmp_MD;
                          FA(:,:,cpt_slc,cpt_b-1)=tmp_FA;
                          Trace_DTI(:,:,cpt_slc,cpt_b-1)=tmp_Trace;
                          
                 end
                 waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
            end
       
    end
    close(h);    

end

function [Tensor,EigVector,EigValue, Mat2] = Tensor_local(slice,dir,bvalue)

% slice : x y [B0 B1-dir1 B2-dir2...]
% dir : x1 y1 z1
%       x2 y2 z2
%        .  .  .
%        .  .  .
%
%  Tensor : x y [6]
%  EigVector : x y [3 coor][3 num]
%  EigValue : x y [3]


nb_dir=size(dir,1);
H = zeros(nb_dir,6);

for i=1:nb_dir
   H(i,:)=[dir(i,1)*dir(i,1),dir(i,2)*dir(i,2),dir(i,3)*dir(i,3),2*dir(i,1)*dir(i,2),2*dir(i,1)*dir(i,3),2*dir(i,2)*dir(i,3)];
end
[U,S,V] = svd(H,0);      
H_inv=V*inv(S)*U';                 %H is a 6*30 matrix



Tensor=zeros(size(slice,1),size(slice,2),6);
EigVector=zeros(size(slice,1),size(slice,2),3,3);
EigValue= zeros(size(slice,1),size(slice,2),3);
Mat=zeros(size(slice,1),size(slice,2),1,9);
Mat2=zeros(size(slice,1),size(slice,2),3,3);
for y=1:1:size(slice,1)
    for x=1:1:size(slice,2)
        Y=[];
        
        if slice(y,x,1)~=0
            for z=1:1:nb_dir
                if double(slice(y,x,z+1))~=0
                    Y=[Y;log(double(slice(y,x,1))/double(slice(y,x,z+1)))/bvalue]; % Y = [log(S0/S1), log(S0/S2), log(S0,S3)....]
                else
                    Y=[Y;0]; % Y = [log(S0/S1), log(S0/S2), log(S0,S3)....]
                end
            end
            Tensor(y,x,:)=H_inv*Y;
            Mat=[Tensor(y,x,1),Tensor(y,x,4),Tensor(y,x,5);Tensor(y,x,4),Tensor(y,x,2),Tensor(y,x,6);Tensor(y,x,5),Tensor(y,x,6),Tensor(y,x,3)];
             if mean(mean(isnan(Mat)))>0 || mean(mean(isinf(Mat)))>0
                Mat=zeros(3,3); 
               % Mat2(y,x,1,:)=[0;0;0;0;0;0;0;0;0];
            end
            Mat2(y,x,:,:)=Mat;

            [Vect,Diag]=eig(Mat);
            
           % EigVector(x,y,:,:)=Vect;

            %EigValue(y,x,:)=[abs(Diag(1,1)/bvalue),abs(Diag(2,2)/bvalue),abs(Diag(3,3)/bvalue)];
            EigValue(y,x,:)=[abs(Diag(1,1)),abs(Diag(2,2)),abs(Diag(3,3))];


           % if((EigValue(y,x,1)<0)&&(EigValue(y,x,2)<0)&&(EigValue(y,x,3)<0)), EigValue(y,x,:)=abs(EigValue(y,x,:));end
           % if(EigValue(y,x,1)<=0), EigValue(y,x,1)=eps; end
           % if(EigValue(y,x,2)<=0), EigValue(y,x,2)=eps; end

            [t,index]=sort(EigValue(y,x,:),'descend');
            EigValue(y,x,:)=EigValue(y,x,index);
            EigVector(y,x,:,:)=Vect(:,index);
        else
            Mat2(y,x,:,:)=0;
            Tensor(y,x,:)=0;
            EigValue(y,x,:)=0;
            EigVector(y,x,:,:)=0;
        end
                
    end
end
 %WriteInrTensorData_KM(Mat2, [1 1 1], 'tensor');
end

function [ADC, FA, TRACE, I1, I2, I3] = Maps_local(EigValue)

ADC=zeros(size(EigValue,1),size(EigValue,2));
FA=zeros(size(EigValue,1),size(EigValue,2));
TRACE=zeros(size(EigValue,1),size(EigValue,2));
I1=zeros(size(EigValue,1),size(EigValue,2));
I2=zeros(size(EigValue,1),size(EigValue,2));
I3=zeros(size(EigValue,1),size(EigValue,2));

for y=1:1:size(EigValue,1)
    for x=1:1:size(EigValue,2)
    TRACE(y,x)=EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3);
    ADC(y,x)=(EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3))/3;
    FA(y,x)=sqrt((3*((EigValue(y,x,1)-ADC(y,x))^2+(EigValue(y,x,2)-ADC(y,x))^2+(EigValue(y,x,3)-ADC(y,x))^2))/(2*(EigValue(y,x,1)^2+EigValue(y,x,2)^2+EigValue(y,x,3)^2)));
    
    I1(y,x)=EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3);
    I2(y,x)=EigValue(y,x,1)*EigValue(y,x,2)+EigValue(y,x,2)*EigValue(y,x,3)+EigValue(y,x,1)*EigValue(y,x,3);
    I3(y,x)=EigValue(y,x,1)*EigValue(y,x,2)*EigValue(y,x,3);
    
    end
end
end
