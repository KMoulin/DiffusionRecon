function Dcm =Read_n_Load_KM()
    dcm_dir = uigetdir;
    cd(dcm_dir);
    listing = dir(dcm_dir);

    k=1;
    for cpt=1:1:size(listing,1)
        cpt
         [FolderName, name, fExt] = fileparts([listing(cpt).folder '\' listing(cpt).name]);
        if (strcmp(fExt, '.dcm') | strcmp(fExt, '.IMA') | isempty(fExt)) & ~listing(cpt).isdir

                Dcm(:,:,k)=double(dicomread(listing(cpt).name));
                k=k+1;

        end
        
    end
end