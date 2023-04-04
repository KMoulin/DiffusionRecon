function [fMap,VMap,S0Map,ErrorMap]= fIVIM_fixed_KM(Dcm, enum,Db,Vb,mask)

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
                                    [f, S0 ResM0]= fitIVIMFlowPerf_multib_fixedV_S0(squeeze(Dcm(y,x,cpt_slc,:,1,1)),enum.alpha(2:end).*enum.alpha(2:end),enum.b(2:end),Db,Vb);
                                    fMap(y,x,cpt_slc,1,1,1)=f;
                                    S0Map(y,x,cpt_slc,1,1,1)=S0;
                                    ErrorMap(y,x,cpt_slc,1,1,1,:)=ResM0;
                                else
                                    fMap(y,x,cpt_slc,1,1,1)=0;
                                    S0Map(y,x,cpt_slc,1,1,1)=0;
                                    ErrorMap(y,x,cpt_slc,1,1,1,1:size(Dcm,4))=0;
                                end
                            end
                            waitbar(y/size(Dcm,1),h);
                        end
                        close(h); 
    end
       

end