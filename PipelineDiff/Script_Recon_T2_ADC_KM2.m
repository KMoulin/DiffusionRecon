
clear all
close all
warning off;

enum=[];
Dcm=[];
DcmB0TEShort=[];

%%%%%%%%%%%%%%% UI Management %%%%%%%%%%%%%%%%%%%%%%
UI=UIDiffRecon_KM(true);
disp('Select Folder');
dcm_dir = uigetdir;
cd(dcm_dir);
mkdir([dcm_dir '/Maps'])
listing = dir(dcm_dir);

%%
%%%%%%%%%%%%%%% Create Enum and Vol %%%%%%%%%%%%%%%%
[Dcm enum]= AnalyseDataSet_KM(listing);
enum.dcm_dir=dcm_dir;
enum.nset=1;

enum.TE(2)=31; % ms
enum.nShortTE=enum.dataset(1).slc(1).b(1).dir(1).nb_avg-enum.dataset(1).slc(1).b(end).dir(1).nb_avg;

save([dcm_dir '/Maps/RAW.mat'],'Dcm','enum');
if UI.gif_mode
    Gif_KM(Dcm, enum, 'Raw')
end

%%
%%%%%%%%%%%%%%% Unmosaic %%%%%%%%%
if UI.mosa_mode && enum.mosa>1
    [Dcm, enum]= Demosa_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/Demosa.mat'],'Dcm','enum');
    if UI.gif_mode
        Gif_KM(Dcm, enum, 'Unmosaic')
    end
end

%%
%%%%%%%%%%%%%%% Registration %%%%%%%%%
if UI.rigid_mode
    [Dcm]= RigidRegistration_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/Rigid.mat'],'Dcm','enum');
    if UI.gif_mode
        Gif_KM(Dcm, enum, 'RigidReg')
    end
end

if UI.Nrigid_mode
    [Dcm]= NonRigidRegistration_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/NonRigid.mat'],'Dcm','enum');
    if UI.gif_mode
        Gif_KM(Dcm, enum, 'NonRigidReg')
    end
end

%%
%%%%%%%%%%%%%%% TE short Separation %%%%%%%%%%%
for cpt_slc=1:1:enum.datasize.slc       
        enum.dataset.slc(cpt_slc).b(1).dir(1).nb_avg=enum.dataset.slc(cpt_slc).b(1).dir(1).nb_avg-enum.nShortTE;    
        DcmB0TEShort(:,:,cpt_slc,1,1,:)=Dcm(:,:,cpt_slc,1,1,enum.dataset.slc(cpt_slc).b(1).dir(1).nb_avg+1:end); 
end

save([enum.dcm_dir '/Maps/DcmB0.mat'],'Dcm','DcmB0TEShort','enum');

%%
%%%%%%%%%%%%%%% Calculate SNR %%%%%%%%%
 if UI.ADC_mode  
    SNR=squeeze(squeeze(mean(Dcm(:,:,cpt_slc,1,1,1:enum.dataset.slc(cpt_slc).b(1).dir(1).nb_avg),6))./std(Dcm(:,:,cpt_slc,1,1,1:enum.dataset.slc(cpt_slc).b(1).dir(1).nb_avg)),[],6);
    save([enum.dcm_dir '/Maps/SNR.mat'],'SNR');
 end 
 
%%
%%%%%%%%%%%%%%%% Average %%%%%%%%%
if UI.avg_mode   
    [Dcm enum]= Average_KM(Dcm, enum); 
    save([enum.dcm_dir '/Maps/Average.mat'],'Dcm','enum');
    if UI.gif_mode
        Gif_KM(Dcm, enum, 'Average')
    end
end

%%
%%%%%%%%%%%%%%%% Average and reject %%%%%%%%%
if UI.avg2_mode 
    [Dcm  enum]= Average_and_Reject_KM(Dcm, enum,3e-3);
    save([enum.dcm_dir '/Maps/Average_Reject.mat'],'Dcm','enum');
     if UI.gif_mode
        Gif_KM(Dcm, enum, 'Average_Reject')
    end
end

%%
%%%%%%%%%%%%%%% PCA %%%%%%%%%
if UI.pca_mode
    [Dcm]= VPCA_KM(Dcm, enum,80);
    save([enum.dcm_dir '/Maps/PCA.mat'],'Dcm','enum');
end

%%
%%%%%%%%%%%%%%% tMIP %%%%%%%%%
if UI.tmip_mode
    [Dcm enum]= tMIP_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/tMIP.mat'],'Dcm','enum');
     if UI.gif_mode
        Gif_KM(Dcm, enum, 'tMIP')
    end
end
%%
%%%%%%%%%%%%%%% Zero filling interpolation %%%%%%%%%
if UI.inter_mode
    [Dcm]= Interpolation_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/Interpolation.mat'],'Dcm','enum');
     if UI.gif_mode
        Gif_KM(Dcm, enum, 'Interpolation')
    end
end 

%%
%%%%%%%%%%%%%%% Create ROI %%%%%%%%%
if UI.roi_mode
    %% For Cardiac ROI
    %[DcmB0,P_Endo,P_Epi,LV_mask]= ROI_KM( DcmB0);
    %save([enum.dcm_dir '/Maps/ROI.mat'],'Dcm','DcmB0','enum','P_Endo','P_Epi','LV_mask');
    
    %% For any Circular ROI
    [Dcm,Dcm_Roi,Mask_ROI,~]= ROI_Ellispe_KM(Dcm,1);
    save([enum.dcm_dir '/Maps/ROI.mat'],'Dcm','enum','Mask_ROI');
end
%%
%%%%%%%%%%%%%%% Create Mask %%%%%%%%%
if UI.mask_mode
    [Mask]= Mask_KM(Dcm(:,:,:,1),100,6000);
    Mask(Mask>0)=1;
    Dcm=Apply_Mask_KM(Dcm,Mask);
    save([enum.dcm_dir '/Maps/Mask.mat'],'Mask');
end
%%
%%%%%%%%%%%%%%% Calculate Tensor %%%%%%%%%
if UI.DTI_mode
    [Tensor,EigValue,EigVector,MD,FA,Trace_DTI] = Calc_Tensor_KM(Dcm, enum);
    save([enum.dcm_dir '/Maps/DTI.mat'],'Tensor','EigValue','EigVector','MD','FA','Trace_DTI');  

%% For Cardiac ROI only    
%    %%%%%%%%%%%%%% Extract HA %%%%%%%%%%%%%
%     if UI.roi_mode
%         EigVect1(:,:,1:size(Dcm,3),:)=squeeze(EigVector(:,:,:,1,:,1));
%         EigVect2(:,:,1:size(Dcm,3),:)=squeeze(EigVector(:,:,:,1,:,2));
%         EigVect3(:,:,1:size(Dcm,3),:)=squeeze(EigVector(:,:,:,1,:,3));
%         Elevation(:,:,1:size(Dcm,3),:)=squeeze(EigVector(:,:,:,1,3,1));
%         HA = HA_KM( EigVect1, Dcm, P_Epi, P_Endo );
%         save([enum.dcm_dir '/Maps/HA.mat'],'EigVect1','EigVect2','EigVect3','Elevation','HA');  
%     end

end

%%
%%%%%%%%%%%%%%% Create Trace %%%%%%%%%
if UI.trace_mode
    [Trace enum]= Trace_KM(Dcm, enum);
    [Trace_Norm]= Norm_KM(Dcm, enum);
    if UI.gif_mode
        Gif_KM(Trace, enum, 'Trace')
    end
    save([enum.dcm_dir '/Maps/Trace.mat'],'Trace','Trace_Norm','enum');
end

%%
%%%%%%%%%%%%%%% Calculate ADC %%%%%%%%%
if UI.ADC_mode   
    [ADC]= ADC_KM(Trace, enum);
    save([enum.dcm_dir '/Maps/ADC.mat'],'ADC');
end

%%
%%%%%%%%%%%%%%% Calculate T2 %%%%%%%%%
if UI.ADC_mode && length(enum.TE)>1 
    [T2Map S0MapT2]= T2fit_KM( mean(DcmB0TEShort(:,:,:,1,1,:),6), mean(Dcm(:,:,:,1,1,:),6), enum); 
    save([enum.dcm_dir '/Maps/T2Map.mat'],'T2Map','S0MapT2');
end
 
%%
%%%%%%%%%%%%%%% Recreate Dicom Maps %%%%%%%%%
%% Required the original ADC/TRACE dicom from Siemens
%if UI.ADC_mode  && UI.trace_mode 
%    [Folder]= Recreate_Dicom_Maps_KM(ADC*1e6,enum,[],'ADCMap',1015);
%	Recreate_Dicom_Maps_KM(T2Map*1e1,enum,Folder,'T2Map',1016);
%end

warning on;