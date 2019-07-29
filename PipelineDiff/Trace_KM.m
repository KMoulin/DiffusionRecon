function [Dcm2 enum2]= Trace_KM(varargin)
    

%  Generate Trace Matrices from the DWI and nDWI matrices by averaging
%  every diffusion directions (practical but incorrect solution) 
%
%  Update the number of directions in enum
%
% SYNTAX:  [Dcm2 enum2]= Trace_KM(Dcm, enum)
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%
%          
% OUTPUTS:  Dcm2 - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum2 - Structure which contains information about the dataset 
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu
    narginchk(2,3);
    if numel(varargin) == 2
          Dcm=varargin{1};
          enum=varargin{2}(1);
          mod=1;

    else
        Dcm=varargin{1};
        enum=varargin{2}(1);
        mod=varargin{3}(1);
    end
          
          
    enum2=enum; 
    Dcm2=[];
    disp('Trace') 
    h = waitbar(0,'Trace...');
     for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
             for cpt_b=1:1:enum.datasize(cpt_set).b             
                       if(cpt_b>1)
                           %Dcm2(:,:,cpt_slc,(cpt_b-1),1,:)=nthroot(tmpDcmB2,enum.dataset(cpt_b).dirNum);
                           if(mod==1)
                               Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=mean(Dcm(:,:,cpt_slc,cpt_b,:,:,cpt_set),5); % For now just mean but later calculate true trace !
                           elseif(mod==2)
                               Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=median(Dcm(:,:,cpt_slc,cpt_b,:,:,cpt_set),5); % For now just mean but later calculate true trace !
                           elseif (mod==3)
                               Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=max(Dcm(:,:,cpt_slc,cpt_b,:,:,cpt_set),[],5); % For now just mean but later calculate true trace !
                           elseif (mod==4)
                               Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=min(Dcm(:,:,cpt_slc,cpt_b,:,:,cpt_set),[],5); % For now just mean but later calculate true trace !
                           else
                               Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=mean(Dcm(:,:,cpt_slc,cpt_b,:,:,cpt_set),5); % For now just mean but later calculate true trace !
                           end
                       else
                           Dcm2(:,:,cpt_slc,cpt_b,1,1,cpt_set )=Dcm(:,:,cpt_slc,cpt_b,1,1 ,cpt_set);
                       end
                       enum2.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir =1;
             end    
           waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
     end
     close(h);    
  
end