function Moco_on_direction

warning off
FolderName = uigetdir;
cd(FolderName);
listing = dir(FolderName);
InfoDcm=[];
Enum.slc=[];
Enum.b=[];
Enum.dir=[];
norot=  [];
i=0;
dicomdict('set','dicom-dict-mosa.txt');
direct=[1.000000 0.414250 -0.414250;
1.000000 -0.414250 -0.414250;
1.000000 -0.414250 0.414250;
1.000000 0.414250 0.414250;
0.414250 0.414250 1.000000;
0.414250 1.000000 0.414250;
0.414250 1.000000 -0.414250;
0.414250 0.414250 -1.000000;
0.414250 -0.414250 -1.000000;
0.414250 -1.000000 -0.414250;
0.414250 -1.000000 0.414250;
0.414250 -0.414250 1.000000;]; %vecteur de direction pour 12 dir si et seulement si les directions ne sont pas inscrites dans les fichiers dicom
for cpt=3:1:size(listing,1)
    if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
        i=i+1;
        tmpInfoDcm=dicominfo(listing(cpt).name);
        InfoDcm(i).name=tmpInfoDcm.Filename;
      if isempty(find(tmpInfoDcm.SequenceName =='#'))
      %if strcmp(tmpInfoDcm.SequenceName,'ep_b5#1')
           InfoDcm(i).b=0;
           InfoDcm(i).dir=1;
           InfoDcm(i).slc = tmpInfoDcm.SliceLocation; 
           InfoDcm(i).dirV=0;
       else             
         %  InfoDcm(i).b = 700;%
           InfoDcm(i).b = 350; %% str2num(tmpInfoDcm.SequenceName(  (find(tmpInfoDcm.SequenceName=='b')+1) :   (find(tmpInfoDcm.SequenceName=='#')-1)    ));
           InfoDcm(i).dir = str2num(tmpInfoDcm.SequenceName((find(tmpInfoDcm.SequenceName=='#')+1)  :   end)); %+1 for old files 
           InfoDcm(i).slc = tmpInfoDcm.SliceLocation;
           if (isfield(tmpInfoDcm, 'DiffusionDirection'))
            ImOrien=tmpInfoDcm.ImageOrientationPatient;
            Im_3=cross(ImOrien(1:3),ImOrien(4:6));
            Im_1=ImOrien(1:3);%attention: this is because matlab read first the colume and second the row
            Im_2=ImOrien(4:6);     
            InfoDcm(i).dirV = ([Im_1';Im_2';Im_3']* tmpInfoDcm.DiffusionDirection)';
            norot(i,:)=tmpInfoDcm.DiffusionDirection';
           else
            ImOrien=tmpInfoDcm.ImageOrientationPatient;
            Im_3=cross(ImOrien(1:3),ImOrien(4:6));
            Im_1=ImOrien(1:3);%attention: this is because matlab read first the colume and second the row
            Im_2=ImOrien(4:6);     
            InfoDcm(i).dirV = ([Im_1';Im_2';Im_3']* direct(InfoDcm(i).dir,:)')';
            norot(i,:)=tmpInfoDcm.DiffusionDirection';
           end
       end
       
       if isempty(find(Enum.slc==InfoDcm(i).slc))
            Enum.slc=[Enum.slc InfoDcm(i).slc];        
        end
        if isempty(find(Enum.b==InfoDcm(i).b))
            Enum.b=[Enum.b InfoDcm(i).b];        
        end
        if isempty(find(Enum.dir==InfoDcm(i).dir))
            Enum.dir=[Enum.dir InfoDcm(i).dir];        
        end 
    end
end

Enum.slc=sort(Enum.slc);
Enum.b=sort(Enum.b);
Enum.dir=sort(Enum.dir);
Enum
dcm_volume=[];
Volume.slice=[];
for cpt=1:1:size(Enum.slc,2)
    Volume.slice(cpt).b=[];
    for cpt2=1:1:size(Enum.b,2)
        Volume.slice(cpt).b(cpt2).dir=[];
        for cpt3=1:1:size(Enum.dir,2)
                for cpt4=1:1:size(InfoDcm,2)
                     if InfoDcm(cpt4).b==Enum.b(cpt2) && InfoDcm(cpt4).dir==Enum.dir(cpt3) && InfoDcm(cpt4).slc==Enum.slc(cpt)
                         
                         Volume.slice(cpt).b(cpt2).dir(cpt3).name=InfoDcm(cpt4).name;
                        Volume.slice(cpt).b(cpt2).dir(cpt3).nameRecon=InfoDcm(cpt4).name;
                        Volume.slice(cpt).b(cpt2).value=InfoDcm(cpt4).b;
                        Volume.slice(cpt).b(cpt2).dir(cpt3).value=InfoDcm(cpt4).dir;
                        Volume.slice(cpt).b(cpt2).dir(cpt3).vector=InfoDcm(cpt4).dirV;      
                            
                        
                     end
                end
        end
    end
end

for cpt=1:1:size(Enum.slc,2)
     i=1;
     direction=[];
     for cpt2=1:1:size(Enum.b,2)
         if (Volume.slice(cpt).b(cpt2).value==0)
            cpt3_max=1;
         else
             cpt3_max=size(Enum.dir,2);
             
         end
          for cpt3=1:1:cpt3_max 
            tmpDataDcm=dicomread(Volume.slice(cpt).b(cpt2).dir(cpt3).name);
            tmpInfoDcm=dicominfo(Volume.slice(cpt).b(cpt2).dir(cpt3).name);
            dcm_volume(:,:,cpt,i)= scale_image(tmpDataDcm,2);
            
            if Volume.slice(cpt).b(cpt2).value~=0
              direction=[direction ; Volume.slice(cpt).b(cpt2).dir(cpt3).vector];
           end
            i=i+1;
          end
     end
end
 

    disp('Organize the images for MoCo2'); % sort and rename the dicom in a 2D way : b and Td, we create a b for each different slice, b value or direction
    h = waitbar(0,'Organize the images for MoCo2...');
    mkdir([FolderName, '\Moco2\']);
    mkdir([FolderName, '\Moco2\all\']);
    cd(FolderName);
   for cpt=1:1:size(Enum.slc,2)
         mkdir([FolderName, '\Moco2\slc',num2str(cpt),'\']);
         mkdir([FolderName, '\Moco2\slc',num2str(cpt),'\dicom\']);
         j=1;
     for cpt2=1:1:size(Enum.b,2)
         if (Volume.slice(cpt).b(cpt2).value==0)
            cpt3_max=1;
         else
             cpt3_max=size(Enum.dir,2);
             
         end
          for cpt3=1:1:cpt3_max               
                    tmpDataDcm = dicomread(Volume.slice(cpt).b(cpt2).dir(cpt3).nameRecon);
                    tmpInfoDcm = dicominfo(Volume.slice(cpt).b(cpt2).dir(cpt3).name);    
                    tmpInfoDcm.ImageComments = ['b=', Volume.slice(cpt).b(cpt2).value];
                    filenameOut=[[FolderName, '\Moco2\slc',num2str(cpt),'\dicom\'], 'V1_IM_b_',num2str(1),'_TD_', num2str(j),'.dcm'];
                    Volume.slice(cpt).b(cpt2).dir(cpt3).nameRecon=[[FolderName, '\Moco2\slc',num2str(cpt),'\moco\'], 'V1_IM_b_',num2str(1),'_TD_', num2str(j),'_moco.dcm'];
                    dicomwrite(uint16(tmpDataDcm),filenameOut, tmpInfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true); 
                    j=j+1;
          end 
     end
     waitbar(cpt/size(Enum.slc,2),h); 
   end
   close(h); 
   
   disp('Execute MoCo2 please wait...');
   h = waitbar(0,'Execute MoCo2 please wait...');
   for cpt=1:1:size(Enum.slc,2)
         PerformCardiacDTIMoCo_20111025([FolderName, '\Moco2\slc',num2str(cpt),'\']);
         waitbar(cpt/size(Enum.slc,2),h); 
   end
    for cpt=1:1:size(Enum.slc,2)  
     for cpt2=1:1:size(Enum.b,2)
         if (Volume.slice(cpt).b(cpt2).value==0)
            cpt3_max=1;
         else
             cpt3_max=size(Enum.dir,2);
             
         end
          for cpt3=1:1:cpt3_max                     
                    tmpDataDcm = dicomread(Volume.slice(cpt).b(cpt2).dir(cpt3).nameRecon);
                    tmpInfoDcm = dicominfo(Volume.slice(cpt).b(cpt2).dir(cpt3).name);    
                    tmpInfoDcm.ImageComments = ['b=', Volume.slice(cpt).b(cpt2).value];
                    filenameOut=[[FolderName, '\Moco2\all\'], 'V1_IM_slc_',num2str(cpt),'_b_',num2str(Volume.slice(cpt).b(cpt2).value),'_dir_',num2str(Volume.slice(cpt).b(cpt2).dir(cpt3).value),'.dcm'];
                    Volume.slice(cpt).b(cpt2).dir(cpt3).nameRecon=filenameOut;
                    dicomwrite(uint16(tmpDataDcm),filenameOut, tmpInfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true);    
          end 
     end
   end
   close(h); 
  warning on
end