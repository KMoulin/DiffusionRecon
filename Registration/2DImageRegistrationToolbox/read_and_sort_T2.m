function Matrix = read_and_sort_T2(FolderName)

    if (~exist('FolderName', 'var'))
            FolderName = uigetdir;        
    end
    cd(FolderName);
    listing = dir(FolderName);

    InfoDcm=[];
    Enum.slc=[];
    Enum.time=[];
    tmpInfoDcm=[];
    i=0;
    for cpt=3:1:size(listing,1)
        if listing(cpt).name(end-2:end) == 'dcm'| (listing(cpt).name(end-2:end) == 'IMA' )
            i=i+1;
            tmpInfoDcm=dicominfo(listing(cpt).name);
            InfoDcm(i).name=tmpInfoDcm.Filename;
            InfoDcm(i).time=str2num(tmpInfoDcm.AcquisitionTime);
            InfoDcm(i).slc = tmpInfoDcm.SliceLocation; 

            if isempty(find(Enum.slc==InfoDcm(i).slc))
                Enum.slc=[Enum.slc InfoDcm(i).slc];        
            end
            if isempty(find(Enum.time==InfoDcm(i).time))
                Enum.time=[Enum.time InfoDcm(i).time]; 
            end
        end
    end

    Enum.slc=sort(Enum.slc);
    Enum.time=sort(Enum.time);
    Enum
    
    for cpt=1:1:size(Enum.slc,2)
        i=1;
        for cpt2=1:1:size(Enum.time,2)
             for cpt3=1:1:size(InfoDcm,2)
                if InfoDcm(cpt3).slc==Enum.slc(cpt) && InfoDcm(cpt3).time==Enum.time(cpt2) 
                    Matrix(i,:,:,cpt)=double(dicomread(InfoDcm(cpt3).name)); 
                    i=i+1;
                end
            end
        end
    end
    

end