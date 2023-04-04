function [fMap,VMap,S0Map,ErrorMap]= fIVIM_fixed_KM2(Dcm, enum,Db,mask)

    fMap=[];
    VMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    for cpt_slc=1:1:size(enum.slc,2)         
                        h = waitbar(0,['fIVIM fit slc :' num2str(cpt_slc) '...']);
                        for y=1:1:size(Dcm,1)                          
                            for x=1:1:size(Dcm,2)
                                if min(squeeze(mask(y,x,cpt_slc,:,1,1))~=0)~=0
                                    tmp=squeeze(Dcm(y,x,cpt_slc,:));
                                    [f, V, ResM0]= fitIVIMFlowPerf_multib([0;tmp]',enum.alpha.*enum.alpha,enum.b,Db);
                                    fMap(y,x,cpt_slc,1,1,1)=f;
                                    VMap(y,x,cpt_slc,1,1,1)=V;
                                    ErrorMap(y,x,cpt_slc,1,1,1,:)=ResM0;
                                else
                                    fMap(y,x,cpt_slc,1,1,1)=0;
                                    VMap(y,x,cpt_slc,1,1,1)=0;
                                    ErrorMap(y,x,cpt_slc,1,1,1,1:(size(Dcm,4)+1))=0;
                                end
                            end
                            waitbar(y/size(Dcm,1),h);
                        end
                        close(h); 
    end
       

end