function [ADC]= ADC_KM(Dcm, enum)

%  Generate ADC maps: ADC= log(S/S0)/(b0-b) 
%  Use the first b-value as the nDWI 
%   
% SYNTAX:  [ADC]= ADC_KM(Dcm, enum)
%  
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%
%          
% OUTPUTS:  ADC - ADC image matrix (units [mm²/s])
%                 [y x slices b-values directions averages dataset]
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html
   
    
    ADC=[];
    disp('ADC calculation') 
    h = waitbar(0,'ADC calculation...');
    for cpt_set=1:1:enum.nset
         for cpt_slc=1:1:enum.datasize(cpt_set).slc
            for cpt_b=1:1:enum.datasize(cpt_set).b     
              for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                 for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg  
                        if(cpt_b>1)
                            ADC(:,:,cpt_slc,(cpt_b-1),cpt_dir,cpt_avg,cpt_set)=ADCMap_local(squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_set)),squeeze(Dcm(:,:,cpt_slc,1,1,1,cpt_set)),-enum.b(cpt_b));
                        end
                end
             end
           end
           waitbar(cpt_slc/size(enum.slc,2),h);
         end
    end
    close(h);    

end

function [ADC] = ADCMap_local(Vol,VolB0,b_vect)

ADC=zeros(size(Vol,1),size(Vol,2));

for y=1:1:size(Vol,1)
    for x=1:1:size(Vol,2)
        if VolB0(y,x)~=0 
%             variables=[squeeze(VolB0(y,x)) squeeze(Vol(y,x,:))'];
%             [p,S] = polyfit( b_vect,log((variables)),1);
%             if p(1)>0
%                 ADC(y,x) = 0; % D must not be <0 !!
%                %disp('Fits on mean, read: D<0 set to 0')
%             else
%                 ADC(y,x) = -p(1);
%             end
            ADC(y,x)=log(Vol(y,x)/VolB0(y,x))/b_vect;
        else
            ADC(y,x) = 0;
        end
    end
end
end