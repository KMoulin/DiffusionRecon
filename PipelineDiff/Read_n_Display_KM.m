
function Dcm =Read_n_Display_KM()
    dcm_dir = uigetdir;
    cd(dcm_dir);
    listing = dir(dcm_dir);

    k=1;
    for cpt=1:1:size(listing,1)
        if length(listing(cpt).name)>3 
            if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'

                figure(1),imagescn(double(dicomread(listing(cpt).name)));
                pause(1);
            end
        end
    end
end
 