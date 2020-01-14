function   [Dcm2, enum2]= Demosa_KM(Dcm, enum)

% Decompose the mosa diffusion matrices on slice diffusion matrices based 
% on the informations find in the dicom headers and stored in 'enum'.
% Update enum to count the number of slices properly
%
%
% SYNTAX:  [Dcm2 enum2]= Demosa_KM(Dcm, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%          
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages dataset]
%                
%           enum2 - Structure which contains information about the dataset
% 
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html
   

    Dcm2=[];
    enum2=enum;
    div=find_division_mosa(enum.mosa)
    if enum.mosa>1
        disp('Unmosaic');
        h = waitbar(0,'Unmosaic...');
        for cpt_set=1:1:enum.nset
            for cpt_slc=1:1:enum.datasize(cpt_set).slc
             for cpt_b=1:1:enum.datasize(cpt_set).b     
               for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg 


                               tmpDataDcm=[];
                               tmpDataDcm2=[];
                               tmpDataDcm=Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set);

                               x = size(tmpDataDcm,1)/div;                     %The information of slices dimensions is also store here :tmpInfoDcm.AcquisitionMatrix(2);
                               y = size(tmpDataDcm,2)/div;                     % and here :tmpInfoDcm.AcquisitionMatrix(3); But recently I had problem with this field so I prefer recalcul 
                               j=1;
                               for cpt3=1:1:div
                                    for cpt4=1:1:div         
                                        if j<=(enum.mosa)                  
                                            tmpDataDcm2(:,:,j)=tmpDataDcm(((cpt3-1)*x)+1:(cpt3*x),((cpt4-1)*y)+1:(cpt4*y));                  
                                            if isempty(find(enum2.slc== (enum2.slc(1)+(j-1)*10)))
                                                enum2.slc=[enum2.slc (enum2.slc(1)+(j-1)*10)];        
                                            end
                                            j=j+1;
                                        end
                                    end
                                end

                                Dcm2(:,:,:,cpt_b,cpt_dir,cpt_avg,cpt_set)=tmpDataDcm2;

                       end           
                   end
               end
               waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
            end
        end
        close(h);    

         for cpt_set=1:1:enum.nset
            enum2.datasize(cpt_set).slc=enum.datasize(cpt_set).slc*enum.mosa;
            for cpt=1:1:enum.mosa
                 enum2.dataset(cpt_set).slc(cpt)=enum2.dataset(cpt_set).slc(1);
            end
         end
        
    end
    


end

function [div]=find_division_mosa(slice_per_mosa)
       if slice_per_mosa>2 & slice_per_mosa<=4          % If there're 2,3 or 4 slices in one mosa : the mosa contains 2x2 slices 
            div=2;
        elseif slice_per_mosa>4 & slice_per_mosa<=9      % If there're 4,5,... 9 slices in one mosa : the mosa contains 3x3 slices, ect...
            div=3;
        elseif slice_per_mosa>9 & slice_per_mosa <=16
            div=4;
        elseif slice_per_mosa>16 & slice_per_mosa <=25
            div=5;
        elseif slice_per_mosa>25 & slice_per_mosa <=36
            div=6;
        elseif slice_per_mosa>36 & slice_per_mosa <=49
            div=7;
        elseif slice_per_mosa>49 & slice_per_mosa <=64
            div=8;
       elseif slice_per_mosa>65 & slice_per_mosa <=81
            div=9;
       elseif slice_per_mosa>82 & slice_per_mosa <=100
            div=10;            
       elseif slice_per_mosa>101 & slice_per_mosa <=121
            div=11;
        elseif slice_per_mosa>121 & slice_per_mosa <=144
            div=12;
       else
           div=1;
       end    
end