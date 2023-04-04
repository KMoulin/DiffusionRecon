function [ADC]= ADC_2B_KM(Dcm, enum)

    ADC=[];
    disp('ADC calculation') 
    h = waitbar(0,'ADC calculation...');
    if(size(enum.b,2)>2)
        for cpt_slc=1:1:size(enum.slc,2)
           for cpt_b=3:1:(size(enum.b,2))
             for cpt_dir=1:1:enum.dataset(cpt_b).dirNum;
                for cpt_avg=1:1:enum.dataset(cpt_b).dir(cpt_dir).avgNum              
                            ADC(:,:,cpt_slc,(cpt_b-2),cpt_dir,cpt_avg)=ADCMap2B_KM(squeeze(Dcm(:,:,cpt_slc,(cpt_b-2),cpt_dir,cpt_avg)),squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,cpt_avg)),-enum.b(cpt_b-1),-enum.b(cpt_b));
                end
             end
           end
           waitbar(cpt_slc/size(enum.slc,2),h);
        end
    end
    close(h);    

end