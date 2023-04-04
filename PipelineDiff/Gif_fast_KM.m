function  Gif_fast_KM(Dcm, FileName)

% Generate an animated Gif for each slice of Dcm
%
%
% SYNTAX:   Gif_KM(Dcm, enum, 'Average')
%          
% INPUTS:   listing - list of files (ex: listing = dir(dcm_dir))
%          
%           enum - existing enumerator
%
%           dataset - specify the number of the given serie 
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

   disp('Gif') 
  
  for cpt_t=1:1:size(Dcm,3)
        tmpDataDcm=squeeze(Dcm(:,:,cpt_t));          
        if cpt_t==1
            imwrite( double( imresize(tmpDataDcm,[4*size(tmpDataDcm,1) 4*size(tmpDataDcm,2)],'nearest') ),[FileName '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);   %%%% First image, delay time = 0.1s         
            j=2;
        else
            imwrite( double( imresize(tmpDataDcm,[4*size(tmpDataDcm,1) 4*size(tmpDataDcm,2)],'nearest') ),[FileName '.gif'],'gif','WriteMode','append','DelayTime',0.1); %%%% Following images

        end
  end
           
 
end