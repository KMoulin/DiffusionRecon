function [Dcm DcmB0 enum]= AnalyseDataSet2_KM(listing)
    

% Read files from a listing, analyses diffusion dimensions contain in the 
% dataset. Generate an Enumerator "enum" which enumerate the informations 
% find inside the dataset(number of slices, directions, bvalues,..) and
% generate a matrix of diffusion and non diffusion weighted images. 
%
%
% SYNTAX:  [Dcm DcmB0 enum]= AnalyseDataSet_KM(listing);
%  
%
% INPUTS:   listing - list of files (ex: listing = dir(dcm_dir))
%          
% OUTPUTS:  Dcm - DWI image matrix 
%                 [y x slices b-values directions averages]
%
%           DcmB0 - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           enum - Structure which contains information about the dataset 
%                  
%           
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    dicomdict('set','dicom-dict-mosa.txt');
    disp('Analyze the data set') 
    h = waitbar(0,'Analyze the data set...');
    
    k=1;
    
    infoDcm=[]; % infoDcm containt, the b value, dir, position and td of each files !
    FolderName=[];
    
    Dcm=[];
    DcmB0=[];
    
    MaxDir=1;
    MaxAvg=1;
    tmp_avg=[];
    
    enum=[];
    enum.slc=[];     
    enum.b=[];
    enum.dataset=[];
    for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
            
            tmpInfoDcm=dicominfo(listing(cpt).name);
            [FolderName,name,ext] = fileparts(listing(cpt).name);
            
            infoDcm(k).name=[name ext]; 
           
            tmpV=tmpInfoDcm.AcquisitionTime;
            
            infoDcm(k).AcqTime= (str2double(tmpV(1,1:2))*60*60+str2double(tmpV(1,3:4))*60+str2double(tmpV(1,5:6))) + str2double(tmpV(1,8))*10^-1 + str2double(tmpV(1,9))*10^-2 + str2double(tmpV(1,10))*10^-3 + str2double(tmpV(1,11))*10^-4;  %s
            
           if  (strcmp(tmpInfoDcm.Directionality,'NONE') || strcmp(tmpInfoDcm.SequenceName,'ep_b0'))                       
               infoDcm(k).b = tmpInfoDcm.Bvalue;
               infoDcm(k).slc = tmpInfoDcm.SliceLocation; 
               infoDcm(k).dirV=[0 0 0];
           else  
               infoDcm(k).b = tmpInfoDcm.Bvalue;
               infoDcm(k).slc = tmpInfoDcm.SliceLocation; 
               if (strcmp(tmpInfoDcm.Directionality,'DIRECTIONAL') &&isfield(tmpInfoDcm, 'DiffusionDirection'))
                    ImOrien=tmpInfoDcm.ImageOrientationPatient;
                    Im_3=cross(ImOrien(1:3),ImOrien(4:6));
                    Im_1=ImOrien(1:3); % attention: this is because matlab read first the colume and second the row
                    Im_2=ImOrien(4:6);     
                    infoDcm(k).dirV = ([Im_1';Im_2';Im_3']* tmpInfoDcm.DiffusionDirection)';
               else
                    disp('diffusion direction not found')    
                    infoDcm(k).dirV = [0 0 0];
               end
           end
           
            if (isfield(tmpInfoDcm, 'Slice_per_mosa'))
                infoDcm(k).slicePerMosa=tmpInfoDcm.Slice_per_mosa;
                enum.mosa=tmpInfoDcm.Slice_per_mosa;
            else 
                infoDcm(k).slicePerMosa=1;
                enum.mosa=1;
            end
             
            if (isfield(tmpInfoDcm, 'TriggerTime')) infoDcm(k).tt= tmpInfoDcm.TriggerTime;
            else infoDcm(k).tt=0;
            end
            
            if isempty(find(enum.b==infoDcm(k).b)) enum.b=[enum.b infoDcm(k).b];        
            end
            
            if isempty(find(enum.slc==infoDcm(k).slc)) enum.slc_pos=[enum.slc infoDcm(k).slc];        
            end
            
            k=k+1;
        end
        waitbar(cpt/size(listing,1),h);  
    end
    close(h);

    
    
 disp('Generate enumerator') % Count how many there are of Td, b value, direction and slc and store all in the enum !
  h = waitbar(0,'Generate enumerator...');   
  
      
    enum.slc_pos
    
    for cpt_b=1:1:size(enum.b,2)
        enum.dataset(cpt_b).dir=[];
        enum.dataset(cpt_b).dirVector=[];
        enum.dataset(cpt_b).dirNum=0;
    end
    tmp_avg.slc=[];
    for cpt_slc=1:1:size(enum.slc,2)
        tmp_avg.slc(cpt_slc).b=[];
    end
    for cpt=1:1:size(infoDcm,2)    
        
                cpt_b=find(infoDcm(cpt).b==enum.b);
                cpt_slc=find(infoDcm(cpt).slc==enum.slc) ;       
                if enum.dataset(cpt_b).dirNum==0
                    enum.dataset.slc(cpt_slc).b(cpt_b).dirVector=infoDcm(cpt).dirV'; 
                    enum.slc(cpt_slc).b(cpt_b).dir(1).avgNum=1;
                    enum.slc(cpt_slc).b(cpt_b).dir(1).avg=[];
                    enum.slc(cpt_slc).b(cpt_b).dir(1).avg(1).tt=infoDcm(cpt).tt; 
                    enum.slc(cpt_slc).b(cpt_b).dirNum=enum.dataset(cpt_b).dirNum+1;        
                    
                    for cpt_slc2=1:1:size(enum.slc,2)        
                            tmp_avg.slc(cpt_slc2).b(cpt_b).dir(1).avgNum=0;
                            tmp_avg.slc(cpt_slc2).b(cpt_b).dir(1).avg=[];
                    end
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(1).avgNum=1;
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(1).avg(1).filename=infoDcm(cpt).name;
                     
%                     if(cpt_b==1)
%                         DcmB0(:,:,cpt_slc,cpt_b,1,1)=double(dicomread(infoDcm(cpt).name));
%                     else
%                         Dcm(:,:,cpt_slc,(cpt_b-1),1,1)=double(dicomread(infoDcm(cpt).name));
%                     end

               elseif sum(ismember(enum.dataset(cpt_b).dirVector',infoDcm(cpt).dirV,'Rows'))~=1


                    enum.dataset(cpt_b).dirVector=[enum.dataset(cpt_b).dirVector infoDcm(cpt).dirV']; 
                    dirNum=(find(ismember(enum.dataset(cpt_b).dirVector',infoDcm(cpt).dirV,'Rows')==1));    % Number of the current direction

                    enum.dataset(cpt_b).dir(dirNum).avgNum=1;
                    enum.dataset(cpt_b).dir(dirNum).avg=[];
                    enum.dataset(cpt_b).dir(dirNum).avg(1).tt=infoDcm(cpt).tt;
                    enum.dataset(cpt_b).dirNum=enum.dataset(cpt_b).dirNum+1;
                    enum.dataset(cpt_b).dir(dirNum).avg(1).filename=infoDcm(cpt).name;
                     for cpt_slc2=1:1:size(enum.slc,2)
                            tmp_avg.slc(cpt_slc2).b(cpt_b).dir(dirNum).avgNum=0;
                            tmp_avg.slc(cpt_slc2).b(cpt_b).dir(dirNum).avg=[];
                     end
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum=1;
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avg(1).filename=infoDcm(cpt).name;
                    
                    if (dirNum>MaxDir)
                        MaxDir=dirNum;
                    end
%                      if(cpt_b==1)
%                         DcmB0(:,:,cpt_slc,dirNum,1)=double(dicomread(infoDcm(cpt).name));
%                     else
%                         Dcm(:,:,cpt_slc,dirNum,1)=double(dicomread(infoDcm(cpt).name));
%                      end

                     
                     
                else

                    dirNum=(find(ismember(enum.dataset(cpt_b).dirVector',infoDcm(cpt).dirV,'Rows')==1));    % Number of the current direction
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum= tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum+1;
                    
                    tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avg(tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum).filename=infoDcm(cpt).name;
                      
                    enum.dataset(cpt_b).dir(dirNum).avgNum=tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum;
                    enum.dataset(cpt_b).dir(dirNum).avg(enum.dataset(cpt_b).dir(dirNum).avgNum).tt=infoDcm(cpt).tt;
                    
                    if (tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum>MaxAvg)
                        MaxAvg=tmp_avg.slc(cpt_slc).b(cpt_b).dir(dirNum).avgNum;
                    end
%                     if(cpt_b==1)
%                         DcmB0(:,:,cpt_slc,cpt_b,dirNum,enum.dataset(cpt_b).dir(dirNum).avgNum)=double(dicomread(infoDcm(cpt).name));
%                     else
%                         Dcm(:,:,cpt_slc,(cpt_b-1),dirNum,enum.dataset(cpt_b).dir(dirNum).avgNum)=double(dicomread(infoDcm(cpt).name));
%                     end
                end 
                waitbar(cpt/size(listing,1),h);
    end
         
    close(h);
    
    %% Memory pre-alloc
     %Dcm=zeros(size(double(dicomread(infoDcm(1).name)),1),size(double(dicomread(infoDcm(1).name)),2),size(enum.slc,2),size(enum.b,2)-1,MaxDir,MaxAvg);
     %DcmB0=zeros(size(double(dicomread(infoDcm(1).name)),1),size(double(dicomread(infoDcm(1).name)),2),size(enum.slc,2),1,1,MaxAvg);          
  
     
     disp('Create Volumes') % Count how many there are of Td, b value, direction and slc and store all in the enum !
     h = waitbar(0,'Create Volumes...'); 
     for cpt_slc=1:1:size(enum.slc,2)
        for cpt_b=1:1:size(enum.b,2)    
            if cpt_b==1
                for cpt_avg=1:1:tmp_avg.slc(cpt_slc).b(cpt_b).dir(1).avgNum
                        DcmB0(:,:,cpt_slc,1,1,cpt_avg)=double(dicomread(tmp_avg.slc(cpt_slc).b(cpt_b).dir(1).avg(cpt_avg).filename));
                end
            else
                for cpt_dir=1:1:MaxDir
                    for cpt_avg=1:1:tmp_avg.slc(cpt_slc).b(cpt_b).dir(cpt_dir).avgNum
                        Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,cpt_avg)=double(dicomread(tmp_avg.slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename));
                    end 
                end
            end
        end
        waitbar(cpt_slc/size(enum.slc,2),h);
     end
     close(h);
        
     
 for cpt_b=2:1:size(enum.b,2)
    Write_grad_file_KM(enum.dataset(cpt_b).dirVector',cpt_b-1);
 end
 
 

function  Write_grad_file_KM(directions,num)
%
    if ~exist (['grad_direction' num2str(num) '.txt'],'file')
        delete (['grad_direction' num2str(num) '.txt']);
    end
    fid = fopen(['grad_direction' num2str(num) '.txt'],'w');
    xt=size(directions,1);

    %fprintf(fid,[num2str(xt)]);

    for ij=1:1:xt
            fprintf(fid,['\n', num2str(directions(ij,1)), ' ', num2str(directions(ij,2)), ' ', num2str(directions(ij,3))]);
    end

    fclose(fid);

end

end