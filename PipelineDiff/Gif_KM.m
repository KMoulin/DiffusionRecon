function  Gif_KM(Dcm, enum, FileName)

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
    h = waitbar(0,['Gif' FileName '...']);
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
            j=1;
            for cpt_b=1:1:enum.datasize(cpt_set).b     
               for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg
                        tmpDataDcm=squeeze(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set));          
                        if j==1
                            imwrite( double( imresize(tmpDataDcm,[4*size(tmpDataDcm,1) 4*size(tmpDataDcm,2)],'nearest') ),[enum.dcm_dir '/Gif/' FileName '_' num2str(cpt_set) '_' num2str(cpt_slc) '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);   %%%% First image, delay time = 0.1s         
                            j=2;
                        else
                            imwrite( double( imresize(tmpDataDcm,[4*size(tmpDataDcm,1) 4*size(tmpDataDcm,2)],'nearest') ),[enum.dcm_dir '/Gif/' FileName '_' num2str(cpt_set) '_' num2str(cpt_slc) '.gif'],'gif','WriteMode','append','DelayTime',0.1); %%%% Following images
                  
                        end
                  end
               end
            end
            waitbar(cpt_slc/size(enum.datasize(cpt_set).slc,2),h);
        end     
    end
     close(h);
end