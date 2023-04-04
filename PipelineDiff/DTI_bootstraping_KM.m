function [bMD, bFA, bE2A, bHA]=DTI_bootstraping_KM(Dcm,enum,LV_Mask, P_Epi, P_Endo,Mask_Depth,nboot)

    tic

    bMD=[];
    bFA=[];
    bE2A=[];
    bHA=[];
    
    b0=mean(Dcm(:,:,1,1,1,:),6);
    for cpt=1:1:nboot   
        
       [order]=floor(rand(5,1)*5+1);
      
       Dcm2=[];
       Dcm2(:,:,1,1,1)=b0;
       Dcm2(:,:,1,2,1:6)=mean(Dcm(:,:,1,2,:,order),6);
       enum2=enum;
       enum2.dataset(1).slc(1).b(1).dir(1).nb_avg=1;
       for cpt_dir=1:1:6
           enum2.dataset(1).slc(1).b(2).dir(cpt_dir).nb_avg=1;
       end
       
       [Dcm2 enum2]= Interpolation_KM(Dcm2, enum2);
       
       
       [Tensor,EigValue,EigVector,MD,FA,Trace_DTI] = Calc_Tensor_KM(Dcm2, enum2);
    
    
       %%%%%%%%%%%%%% Extract HA %%%%%%%%%%%%%
       EigVect1(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,:,1));
       EigVect2(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,:,2));
  
        
        %[HA TRA]= HA_KM( EigVect1, Dcm, P_Epi, P_Endo );
        [HA2 TRA2 E2A RAD_s CIR_s LON_s]= HA_E2A_KM(EigVect1, EigVect2, LV_Mask, P_Epi, P_Endo);
        %[HA_filter]= HA_Filter_KM(HA,LV_Mask ,Mask_Depth,0);

        [HA_filter2]= HA_Filter_KM(HA2,LV_Mask ,Mask_Depth,0);
       % save([enum.dcm_dir '/Maps/DTI.mat'],'Tensor','EigValue','EigVector','MD','FA','Trace_DTI'); 
       
       bHA(:,:,cpt)=HA_filter2;
       bE2A(:,:,cpt)=E2A;
       bMD(:,:,cpt)=MD;
       bFA(:,:,cpt)=FA;
       
    end
    
    toc
end