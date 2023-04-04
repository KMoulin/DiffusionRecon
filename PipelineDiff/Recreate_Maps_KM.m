function Recreate_Maps_KM(ADC,Trace,Trace_Norm,T2map,DcmB0,enum,folder)

disp('Save Dicom Maps')

infoDcm=[];
slc=[];

[path name_folder ext]=fileparts(folder);
ListFolder=dir([folder '\..']);
List=struct2cell(ListFolder);
List= cellstr(List(1,:));

IndexOld=(~cellfun(@isempty,strfind(List,'_KM_r')));
ITrace=(~cellfun(@isempty,strfind(List,[name_folder '_TRACEW']))); %
%ITrace=(ITrace-IndexOld);
%ITrace(ITrace<0)=0;
IndexTrace=find(ITrace);
IADC=(~cellfun(@isempty,strfind(List,[name_folder '_ADC']))); %
%IADC=(IADC-IndexOld);
%IADC(IADC<0)=0;
IndexADC=find(IADC);

if(~isempty(IndexADC))
       mkdir([folder '\..\' name_folder '_ADC_KM/']);
       listing = dir([folder '\..\' ListFolder(IndexADC).name]);
       jj=1;
       slc=[];
       for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
    
            tmpinfoDcm=dicominfo([path '\' ListFolder(IndexADC).name '\' listing(cpt).name]);
            if ~isempty(find(tmpinfoDcm.SequenceName =='b'))
                    infoDcm(jj).name=[path '\'  ListFolder(IndexADC).name '\' listing(cpt).name];
                    infoDcm(jj).b =str2num(tmpinfoDcm.SequenceName(  (find(tmpinfoDcm.SequenceName=='b')+3:end)));
                    infoDcm(jj).slice=tmpinfoDcm.SliceLocation;
                    slc=[slc abs(tmpinfoDcm.SliceLocation)];
                    jj=jj+1;
            end
            
        end
       end
       [S I]=sort(slc);
       for cpt_slc=1:1:size(ADC,3)
       	for cpt_b=1:1:size(ADC,4)
           tmpinfoDcm=dicominfo(infoDcm(I(cpt_slc)).name);
           tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription 'B' num2str(enum.b(cpt_b+1)) '_reconstruct_KM'];
           tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+1000+14+cpt_b;
           tmpdataDcm=dicomread(infoDcm(I(cpt_slc)).name);
           tmpdataDcm=zeros(size(tmpdataDcm));
           
           tmpdataDcm=ADC(:,:,I(cpt_slc),cpt_b).*1000000;
           
           
           [pathstr,name_file,ext] = fileparts(infoDcm(I(cpt_slc)).name); 
           
           dicomwrite(uint16(tmpdataDcm),[folder '_ADC_KM/' name_file '_Bval' num2str(enum.b(cpt_b+1)) '_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true); 
        end
       end
end
 
if(~isempty(IndexTrace))
       mkdir([folder '\..\' name_folder '_TRACEW_KM/']);
       mkdir([folder '\..\' name_folder '_TRACEW_NORM_KM/']);
       listing = dir([folder '\..\' ListFolder(IndexTrace).name]);
       jj=1;
       slc=[];
       for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
    
            tmpinfoDcm=dicominfo([path '\' ListFolder(IndexTrace).name '\' listing(cpt).name]);
            if ~isempty(find(tmpinfoDcm.SequenceName =='b'))
                    infoDcm(jj).name=[path '\'  ListFolder(IndexTrace).name '\' listing(cpt).name];
                    infoDcm(jj).b =350;%tmpinfoDcm.Bvalue;% str2num(tmpinfoDcm.SequenceName(  (find(tmpinfoDcm.SequenceName=='b')+3:end)));
                    infoDcm(jj).slice=tmpinfoDcm.SliceLocation;
                    if isempty(find(slc==infoDcm(jj).slice))
                        slc=[slc infoDcm(jj).slice];        
                    end
                    jj=jj+1;
            end
            
        end
       end
       [S I]=sort(slc);
       for cpt=1:1:size(infoDcm,2)
           tmpinfoDcm=dicominfo(infoDcm( cpt).name);
           tmpdataDcm=dicomread(infoDcm( cpt).name);
           tmpdataDcm=zeros(size(tmpdataDcm));   
           
           cpt_slc=find(infoDcm(cpt).slice==slc)
           cpt_b=find(infoDcm(cpt).b==enum.b)
           if cpt_b==1
                tmpdataDcm=DcmB0(:,:,cpt_slc);
               [pathstr,name_file,ext] = fileparts(infoDcm(cpt).name);
               tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription '_Reconstruct_KM'];
               tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+1000+12;
               dicomwrite(uint16(tmpdataDcm),[folder '_TRACEW_KM/' name_file '_Bval' num2str(enum.b(cpt_b)) '_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true);
           else
               tmpdataDcm=Trace(:,:,cpt_slc,cpt_b-1);
               [pathstr,name_file,ext] = fileparts(infoDcm(cpt).name);
               tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription '_Reconstruct_KM'];
               tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+1000+12;
               dicomwrite(uint16(tmpdataDcm),[folder '_TRACEW_KM/' name_file '_Bval' num2str(enum.b(cpt_b)) '_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true);
               
               tmpdataDcm=Trace_Norm(:,:,cpt_slc,cpt_b-1);
               [pathstr,name_file,ext] = fileparts(infoDcm(cpt).name);
               tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription '_NORM_Reconstruct_KM'];
               tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+1000+13;
               dicomwrite(uint16(tmpdataDcm),[folder '_TRACEW_NORM_KM/' name_file '_Bval' num2str(enum.b(cpt_b)) '_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true);
           end
       end
 end

  if(~isempty(IndexTrace))
       mkdir([folder '\..\' name_folder '_T2_KM/']);
       listing = dir([folder '\..\' ListFolder(IndexTrace).name]);

       for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
    
               tmpinfoDcm=dicominfo([path '\' ListFolder(IndexTrace).name '\' listing(cpt).name]);      
               tmpdataDcm=T2map.*10;
               [pathstr,name_file,ext] = fileparts(listing(cpt).name);
               tmpinfoDcm.SeriesDescription=[tmpinfoDcm.SeriesDescription '_NORM_Reconstruct_KM'];
               tmpinfoDcm.SeriesNumber= tmpinfoDcm.SeriesNumber+1000+14;
               dicomwrite(uint16(tmpdataDcm),[folder '_T2_KM/' name_file '_T2_KM' ext], tmpinfoDcm, 'CreateMode', 'copy', 'WritePrivate' ,true);
           
        end
       end
  end

end