function [struct ListNames ListSeries] =Analyze_Dcm_Content_KM()
    dcm_dir = uigetdir;
    cd(dcm_dir);
    listing = dir(dcm_dir);
    
    
    Dcm=zeros(500,500);
    ListNames={};
    ListSeries=[];
    k=1;
    for cpt=1:1:size(listing,1)
        if length(listing(cpt).name)>3 
            if ~listing(cpt).isdir
                
                
                
                tmpInfo=dicominfo(listing(cpt).name);
                tmpDcm=dicomread(listing(cpt).name);
                %Dcm(:,:,k)=double(dicomread(listing(cpt).name));
                %Dcm(k)=tmpInfo.SeriesNumber;
                k=k+1;
                
                
                if isempty(strcmp(ListNames,tmpInfo.PatientID)) 
                    ListNames{end+1}= tmpInfo.PatientID;      
                end
                cpt_name=find(strcmp(ListNames,tmpInfo.PatientID));
                
                if isempty(find(ListSeries==tmpInfo.SeriesNumber)) 
                    ListSeries(end+1)= tmpInfo.SeriesNumber;     
                    cpt_series=find(ListSeries==tmpInfo.SeriesNumber);
                    struct(cpt_name).series(cpt_series).Dcm=[];
                     struct(cpt_name).series(cpt_series).Listing=[];
                end
                cpt_series=find(ListSeries==tmpInfo.SeriesNumber);
                
                % Build the dicom structure
                bdiff=~isempty(strfind(tmpInfo.ImageType,'DIFFUSION'));
                btr=~isempty(strfind(tmpInfo.ImageType,'TRACEW'));
                btr0=~isempty(strfind(tmpInfo.ImageType,'TENSOR_B0'));
                brgb=~isempty(strfind(tmpInfo.ImageType,'RGB'));
                badc=~isempty(strfind(tmpInfo.ImageType,'ADC'));;
                bfa=~isempty(strfind(tmpInfo.ImageType,'FA'));
                bten=~isempty(strfind(tmpInfo.ImageType,'TENSOR'))&isempty(strfind(tmpInfo.ImageType,'TENSOR_B0'));
                struct(cpt_name).name=tmpInfo.PatientID;
                struct(cpt_name).series(cpt_series).name=tmpInfo.SeriesNumber;
                struct(cpt_name).series(cpt_series).Dcm(:,:,end+1)=tmpDcm(:,:,1);
                struct(cpt_name).series(cpt_series).Listing(end+1).name=listing(cpt).name;
                struct(cpt_name).series(cpt_series).Listing(end).isdir=0;
                struct(cpt_name).series(cpt_series).Listing(end).isdiff=bdiff;
                struct(cpt_name).series(cpt_series).Listing(end).istracew=btr;
                struct(cpt_name).series(cpt_series).Listing(end).istraceb0=btr0;
                struct(cpt_name).series(cpt_series).Listing(end).RGB=brgb;
                struct(cpt_name).series(cpt_series).Listing(end).isadc=badc;
                struct(cpt_name).series(cpt_series).Listing(end).isfa=bfa;
                struct(cpt_name).series(cpt_series).Listing(end).istensor=bten;
                
                struct(cpt_name).series(cpt_series).Listing(end).isdwi= bdiff & ~btr & ~btr0 & ~badc & ~bfa & ~bten;
                
                
            end
        end
    end
end