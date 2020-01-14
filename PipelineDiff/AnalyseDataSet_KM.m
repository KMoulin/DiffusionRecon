function [Dcm enum]= AnalyseDataSet_KM(varargin)
    

% Read files from a listing, analyses diffusion dimensions contain in the 
% dataset. Generate an Enumerator "enum" which enumerate the informations 
% find inside the dataset(number of slices, directions, bvalues,..) and
% generate a matrix of diffusion and non diffusion weighted images. 
%
%
% SYNTAX:   [Dcm DcmB0 enum]= AnalyseDataSet_KM(listing);
%          
%           [Dcm DcmB0 enum]= AnalyseDataSet_KM(listing,enum);    
%           
%           [Dcm DcmB0 enum]= AnalyseDataSet_KM(listing,enum,dataset);  
%
% INPUTS:   listing - list of files (ex: listing = dir(dcm_dir))
%          
%           enum - existing enumerator
%
%           dataset - specify the number of the given serie 
%
% OUTPUTS:  Dcm - DWI image matrix 
%                 [y x slices b-values directions averages dataset]
%
%           
%           enum - Structure which contains information about the dataset                   
%           
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    narginchk(1,3);
    if numel(varargin) == 1
          listing=varargin{1};
          dataset_num=1;
          enum=[];
          enum.slc=[];     
          enum.b=[];
          enum.dataset=[];
          enum.dataset(dataset_num).slc=[];
          

    elseif numel(varargin) == 2
          listing=varargin{1};
          enum=varargin{2}(1);
          dataset_num=1;
    else
          listing=varargin{1};
          enum=varargin{2}(1);
          dataset_num=varargin{3}(1);
    end

    dicomdict('set','dicom-dict-mosa.txt');
    disp('Analyze the data set') 
    h = waitbar(0,'Analyze the data set...');
    
    k=1;
    
    infoDcm=[]; % infoDcm containt, the b value, dir, position and td of each files !
    FolderName=[];
     VectTime=[];
    Dcm=[];

    
%     enum=[];
%     enum.slc=[];     
%     enum.b=[];
%     enum.dataset=[];
%     enum.dataset(dataset_num).slc=[];
    for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
            
            tmpInfoDcm=dicominfo(listing(cpt).name);
            
            if (isfield(tmpInfoDcm, 'Directionality'))  % check if the file is a valid diffusion file
                [FolderName,name,ext] = fileparts(listing(cpt).name);

                infoDcm(k).name=[name ext]; 

                tmpV=tmpInfoDcm.AcquisitionTime;

                infoDcm(k).AcqTime= (str2double(tmpV(1,1:2))*60*60+str2double(tmpV(1,3:4))*60+str2double(tmpV(1,5:6))) + str2double(tmpV(1,8))*10^-1 + str2double(tmpV(1,9))*10^-2 + str2double(tmpV(1,10))*10^-3 + str2double(tmpV(1,11))*10^-4;  %s
                VectTime(k)=infoDcm(k).AcqTime;
           
               if  (strcmp(tmpInfoDcm.Directionality,'NONE') || strcmp(tmpInfoDcm.SequenceName,'ep_b0'))                       
                   infoDcm(k).b = 0;% tmpInfoDcm.Bvalue;
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

                if isempty(find(enum.slc==infoDcm(k).slc)) enum.slc=[enum.slc infoDcm(k).slc];        
                end

                enum.TE(1)=tmpInfoDcm.EchoTime;
                enum.TR=tmpInfoDcm.RepetitionTime;
                enum.Pixel=tmpInfoDcm.PixelSpacing;
                enum.Thickness=tmpInfoDcm.SliceThickness;
                k=k+1;
           end
        end
        waitbar(cpt/size(listing,1),h);  
    end
    close(h);

    
    %%
    disp('Generate enumerator') % Count how many there are of Td, b value, direction and slc and store all in the enum !
    h = waitbar(0,'Generate enumerator...');   
  
    %%% Build Enum Size %%%  
    enum.datasize(dataset_num).b=size(enum.b,2);
    enum.datasize(dataset_num).slc=size(enum.slc,2);
    enum.datasize(dataset_num).dir=1;
    enum.datasize(dataset_num).avg=1;
    
    %%% Build Enum Data %%% 
    for cpt_slc=1:1:enum.datasize.slc
        enum.dataset(dataset_num).slc(cpt_slc).b=[];
        for cpt_b=1:1: enum.datasize.b
             enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir=[];
             enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir=0;
             enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector=[];
            for cpt_dir=1:1: enum.datasize.dir
                 enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg=0;
                 enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg=[];
               for cpt_avg=1:1: enum.datasize.avg  
                    enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).tt=0;
                    enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).acquTime=0;
                    enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename=0;
                    enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).dirVector=[];
                end
            end
        end
    end

    cpt_b=0;
    cpt_slc=0;
    cpt_dir=0;
    cpt_avg=0;
    
    %%% Fill Enum %%%
    for cpt=1:1:size(infoDcm,2)    
        cpt_b=find(infoDcm(cpt).b==enum.b); % Current b
        cpt_slc=find(infoDcm(cpt).slc==enum.slc) ; % Current slc   
        
        %%% Build first Direction found %%%
        if enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir==0
            cpt_dir=1;
            cpt_avg=1;
            
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir=1;
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir=[];
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector=infoDcm(cpt).dirV'; 
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(1).nb_avg=1;
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(1).avg=[];
            
            
        %%% New direction %%%
       elseif sum(ismember(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector',infoDcm(cpt).dirV,'Rows'))~=1
            
            cpt_dir=enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir+1;    % Number of the current direction
            cpt_avg=1;
            
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir=enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir+1;
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector=[enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector infoDcm(cpt).dirV']; 
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg=1;
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg=[];
            

            
        
        %%% Direction already found
        else
            cpt_dir=(find(ismember(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dirVector',infoDcm(cpt).dirV,'Rows')==1));    % Number of the current direction
            cpt_avg=enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg+1; 
            
            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg=enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg+1;
            
%            enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg)=[];
          
        end 
        
        %%% Store Informations
        
        enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).tt=infoDcm(cpt).tt;       
        enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).dirVector=infoDcm(cpt).dirV;        
        enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).acquTime=infoDcm(cpt).AcqTime-min(VectTime);
        enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename=infoDcm(cpt).name;
              
            
         if cpt_dir>   enum.datasize(dataset_num).dir
             enum.datasize(dataset_num).dir=cpt_dir;
         end
         
         if cpt_avg>   enum.datasize(dataset_num).avg
             enum.datasize(dataset_num).avg=cpt_avg;
         end
            
            
        waitbar(cpt/size(listing,1),h);
    end
         
    close(h);
    
    %%
    %%% Create Dcm File %%%
    
      
     disp('Create Volumes') % Count how many there are of Td, b value, direction and slc and store all in the enum !
     h = waitbar(0,'Create Volumes...'); 
     
     %%% Memory pre-alloc
     Dcm=zeros(size(double(dicomread(infoDcm(1).name)),1),size(double(dicomread(infoDcm(1).name)),2), enum.datasize(dataset_num).slc, enum.datasize(dataset_num).b, enum.datasize(dataset_num).dir, enum.datasize(dataset_num).avg);

     for cpt_slc=1:1:enum.datasize(dataset_num).slc
       for cpt_b=1:1:enum.datasize(dataset_num).b     
          for cpt_dir=1:1: enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir           
             for cpt_avg=1:1: enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg
                 
                  tmp_dcm=double(dicomread(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename));
                  if size(tmp_dcm,1)==size(Dcm,1)
                    Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)=double(dicomread(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename));
                  else
                    Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)=double(dicomread(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).filename)'); 
                  end
             end
           end
        end
        waitbar(cpt_slc/size(enum.slc,2),h);
      end
      close(h);
     
     disp('Generate timing tab') % Count how many there are of Td, b value, direction and slc and store all in the enum !
     h = waitbar(0,'Generate timing tab...'); 
       
     VectTime=VectTime-min(VectTime);
     VectTime=sort(VectTime);
     VectRR=diff(VectTime);
     VectRR=[VectRR VectRR(end)];
     enum.Timing=[];
     if enum.mosa==1
         for cpt_time=1:1:size(VectTime,2)
             for cpt_slc=1:1:enum.datasize(dataset_num).slc
               for cpt_b=1:1:enum.datasize(dataset_num).b     
                  for cpt_dir=1:1: enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).nb_dir           
                     for cpt_avg=1:1: enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg 
                             if(enum.dataset(dataset_num).slc(cpt_slc).b(cpt_b).dir(cpt_dir).avg(cpt_avg).acquTime==VectTime(cpt_time))
                                 enum.Timing.time(cpt_time)=VectTime(cpt_time);
                                 enum.Timing.RR(cpt_time)=1000*VectRR(cpt_time);
                                 enum.Timing.slice(cpt_time)=cpt_slc;
                                 enum.Timing.b(cpt_time)=cpt_b;
                                 enum.Timing.dir(cpt_time)=cpt_dir;
                                 enum.Timing.avg(cpt_time)=cpt_avg;
                             end
                     end
                   end
               end
             end
             waitbar(cpt_time/size(VectTime,2),h);
         end
     else
         cpt=1;
         for cpt_time=1:1:size(VectTime,2)            
               for cpt_b=1:1:enum.datasize(dataset_num).b     
                  for cpt_dir=1:1: enum.dataset(dataset_num).slc(1).b(cpt_b).nb_dir           
                     for cpt_avg=1:1: enum.dataset(dataset_num).slc(1).b(cpt_b).dir(cpt_dir).nb_avg 
                             if(enum.dataset(dataset_num).slc(1).b(cpt_b).dir(cpt_dir).avg(cpt_avg).acquTime==VectTime(cpt_time))
                                for cpt_slc=1:1:enum.mosa
                                % slice_order=[1 3 5 2 4];
                                 slice_order=[1:1:enum.mosa];
                                 enum.Timing.time(cpt)=VectTime(cpt_time);
                                 enum.Timing.RR(cpt)=double(1000*VectRR(cpt_time)./enum.mosa);
                                 enum.Timing.slice(cpt)=double(slice_order(cpt_slc));
                                 enum.Timing.b(cpt)=cpt_b;
                                 enum.Timing.dir(cpt)=cpt_dir;
                                 enum.Timing.avg(cpt)=cpt_avg;
                                 cpt=cpt+1;
                                end
                             end
                     end
                   end
               end    
             waitbar(cpt_time/size(VectTime,2),h);
         end    
     end
      close(h);
    
         
 for cpt_b=2:1:enum.datasize.b  
    Write_grad_file_KM(enum.dataset(dataset_num).slc(1).b(cpt_b).dirVector',cpt_b-1);
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