function [Dcm2 enum2]= Interpolation_KM(Dcm, enum)


% Interpolate the DWI and nDWI matrices by zeros filling
%
%
% SYNTAX:  [Dcm2 DcmB02]= Demosa_KM(Dcm, DcmB0, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - Interpolated DWI image matrix 
%                 [y*2 x*2 slices b-values directions averages dataset]
%
%           enum2 - Structure which contains information about the dataset 
%           
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html
   
    Dcm2=[];
    disp('Zero filling data') 
    h = waitbar(0,'Zero filling data...');
    enum2=enum;
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
         for cpt_b=1:1:enum.datasize(cpt_set).b     
           for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir  
                   for cpt_avg=1:1:enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg
                         tmpZero=[];
                         tmpKspace=[];
                         tmpZero=zeros(2*size(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg),1),2*size(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_set),2));
                         tmpKspace=fftshift(fft2(squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_set))));
                         tmpZero((size(tmpKspace,1)/2):(size(tmpKspace,1)+size(tmpKspace,1)/2-1),(size(tmpKspace,2)/2):(size(tmpKspace,2)+size(tmpKspace,2)/2-1))=tmpKspace;
                         Dcm2(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_set)=abs(ifft2(fftshift(tmpZero)))*4;             
                   end           
               end
           end
           waitbar(cpt_slc/size(enum.slc,2),h);
        end
    end
    enum2.Pixel=enum.Pixel/2;
    close(h);    

end