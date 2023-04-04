function [Dcm Dcm_ph] =Read_n_Load_MS_KM()
    dcm_dir = uigetdir;
    cd(dcm_dir);
    listing = dir(dcm_dir);

    list_rep=[];
    list_cha=[];
    list_rep_ph=[];
    list_cha_ph=[];
    cpt_shot=1;
    cpt_shot_ph=1;
    for cpt=1:1:size(listing,1)
        if length(listing(cpt).name)>3 
            if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA' | listing(cpt).name(end-2:end) == 'ima'

                
                 
                info=dicominfo(listing(cpt).name);
                sseries=(split(info.SeriesInstanceUID,"."));
                
                if (str2num(sseries{1})==str2num(sseries{7}))
                
                    cpt_cha= str2num(sseries{2})+1;
                    cpt_rep= str2num(sseries{3})+1;
                    if str2num(sseries{4})==0 %Mag
                        
                         trig=0;
                         if ~ismember(cpt_cha,list_cha)
                            list_cha=[list_cha cpt_cha];
                            trig=trig+1;
                            
                        end
                        if ~ismember(cpt_rep,list_rep)
                            list_rep=[list_rep cpt_rep];
                            trig=trig+1;
                        end
                        
                        if trig==2
                           cpt_shot= cpt_shot+1;
                        end
                        
                        Dcm(:,:,cpt_cha,cpt_rep,cpt_shot) =double(dicomread(listing(cpt).name));
                        
                       
                           
                    else
                        
                        
                         if ismember(find(list_ph(:,1)==cpt_cha),find(list_ph(:,2)==cpt_rep)) % the pair is found in the table
                            cpt_shot_ph=cpt_shot_ph+1;    
                         else
                             
                         end
                             trig=0;
                         if ~ismember(cpt_cha,list_cha_ph)
                            list_ph=[list_cha_ph cpt_cha];
                            trig=trig+1;
                            
                        end
                        if ~ismember(cpt_rep,list_rep_ph)
                            list_rep_ph=[list_rep_ph cpt_rep];
                            trig=trig+1;
                        end
                        
                        if trig==2
                           cpt_shot_ph= cpt_shot_ph+1;
                        end
                         Dcm_ph(:,:,cpt_cha,cpt_rep,cpt_shot_ph) =double(dicomread(listing(cpt).name));
                    end
                end
            end
        end
    end
end