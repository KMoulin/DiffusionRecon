 function [Folder]= Recreate_Dicom_Maps_KM(varargin)

disp('Save Dicom Maps')

narginchk(2,5);
if numel(varargin) == 2
      Map=varargin{1};
      enum=varargin{2}(1);
      Folder = uigetdir;
      Folder_name = 'Map';    
      Serie_Num=1014;
      
elseif numel(varargin) == 3
      Map=varargin{1};
      enum=varargin{2}(1);
      if isempty(varargin{3})
         Folder = uigetdir;
      else
          Folder=varargin{3};
      end
      Folder_name = 'Map';    
      Serie_Num=1014;
      
elseif numel(varargin) == 4
    Map=varargin{1};
      enum=varargin{2}(1);
       if isempty(varargin{3})
         Folder = uigetdir;
      else
          Folder=varargin{3};
      end
      Folder_name = varargin{4};
       Serie_Num=1014;
else
      Map=varargin{1};
      enum=varargin{2}(1);
       if isempty(varargin{3})
         Folder = uigetdir;
      else
          Folder=varargin{3};
      end
      Folder_name = varargin{4};
      Serie_Num=varargin{5}(1);
end


listing = dir([Folder]);


Old_Folder=regexp(enum.dcm_dir,'\','split');

New_Folder=[enum.dcm_dir '\..\' Old_Folder{end} '_' Folder_name '_KM/'];
mkdir(New_Folder);



       k=1;
       slc=[];
       b_val=[];
       for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
    
            tmpinfoDcm=dicominfo([Folder '/' listing(cpt).name]);
            if ~isempty(find(tmpinfoDcm.SequenceName =='b'))
                    
                    infoDcm(k).name=[Folder '/' listing(cpt).name];
                    infoDcm(k).b =str2num(tmpinfoDcm.SequenceName(  (find(tmpinfoDcm.SequenceName=='b')+3:end)));
                    infoDcm(k).slc=tmpinfoDcm.SliceLocation;
                    
                    if isempty(find(slc==infoDcm(k).slc)) slc=[slc infoDcm(k).slc];        
                    end
                    k=k+1;
            end
            
        end
       end


       
       
       for cpt_slc=1:1:size(Map,3)
       	for cpt_b=1:1:size(Map,4)
            for cpt_i=1:1:k-1
                if  infoDcm(cpt_i).slc==slc(cpt_slc)                     
                   tmpinfoDcm=dicominfo(infoDcm(cpt_i).name);
                   tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription '_' Folder_name '_reconstruct_KM'];
                   tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+Serie_Num+cpt_b;
                   tmpdataDcm=dicomread(infoDcm(cpt_i).name);
                   tmpdataDcm=zeros(size(tmpdataDcm));
                   tmpdataDcm=Map(:,:,cpt_slc,cpt_b);
                   [pathstr,name_file,ext] = fileparts(infoDcm(cpt_i).name); 
                   dicomwrite(uint16(tmpdataDcm),[New_Folder '/ ' Folder_name '_slc' num2str(cpt_slc) '_Bval' num2str(cpt_b) '_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true); 
                end
            end
        end
       end
end
