
function [UI]=UIDiffRecon_KM(options_mode)
     

%  Generate an user interface with checkable otpions that returns a matrix
%  'UI' which contains the options chosen as boolean. Options can be save and
%  will be automatically loaded the next time the function is called (saved
%  as a matrix named 'defaultUI.mat' located in the same folder as this
%  function)
%
% SYNTAX:  [UI]=UIDiffRecon_KM(options_mode);
%  
%
% INPUTS:   options_mode:
%               - true   Display the interface and return the values
%                           when the user click on the button OK
%               - false:  Return the default options
%          
% OUTPUTS:  UI - Struct of booleans for each options   
%
%   
% Kevin Moulin 07/18/2016 
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    if (~exist('options_mode', 'var'))
        options_mode = false;
    else
        if options_mode ~= true && options_mode ~=false 
            options_mode = false;
        end
    end
    
    [UI]=load_default();

if options_mode==true
    disp('Select Options');
    h_fig =   figure( 'BackingStore', 'on','DoubleBuffer','on','Render', 'zbuffer','Name', 'Pcatmip Options','NumberTitle','off','MenuBar','none','DockControls', 'off','Toolbar','none','units', 'characters','Position',[10 10 80 22],'Units','characters');
    hp_config = uipanel( 'Title','Reconstruction parameter','units', 'characters','Position',[1 1 78 20]);
    
    %%%% Ligne 1 %%%%
    hu_inter =           uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[1 15 25 3], 'String', 'Interpolation','Value', UI.inter_mode);
    hu_avg2 =           uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[1 13 25 3], 'String', 'Avg and reject','Value', UI.avg2_mode);
    hu_pca =            uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[1 11 25 3], 'String', 'Pca','Value', UI.pca_mode);
    hu_tmip =           uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[1 9 25 3], 'String', 'Tmip','Value', UI.tmip_mode);
    hu_avg =            uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[1 7 25 3], 'String', 'Avg','Value', UI.avg_mode);
      
    %%%% Ligne 2 %%%%
    hu_rigid =          uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[25 15 25 3], 'String', 'Rigid Reg','Value', UI.rigid_mode);
    hu_Nrigid =         uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[25 13 25 3], 'String', 'Non Rigid Reg','Value', UI.Nrigid_mode);
    hu_mask =           uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[25 11 25 3], 'String', 'Mask','Value', UI.mask_mode);
    hu_roi =            uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[25 9 25 3], 'String', 'ROI','Value', UI.roi_mode);
    hu_trace =          uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[25 7 25 3], 'String', 'Trace','Value', UI.trace_mode);
    
    %%%% Ligne 3 %%%%
    hu_dti=             uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[50 15 25 3], 'String', 'DTI','Value', UI.DTI_mode);
    hu_adc=             uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[50 13 25 3], 'String', 'ADC','Value', UI.ADC_mode);
    hu_ivim=            uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[50 11 25 3], 'String', 'IVIM','Value', UI.IVIM_mode);
    hu_mosa =           uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[50 9 25 3], 'String', 'Unmosaic','Value', UI.mosa_mode);
    hu_gif =            uicontrol('parent', hp_config, 'Style','checkbox', 'units', 'characters', 'Position',[50 7 20 3], 'String', 'Gif','Value', UI.gif_mode);
      
    hu_valid =          uicontrol('parent', hp_config, 'style','pushbutton','units','characters', 'TooltipString', 'Run the function','tag', 'hu_info','String', 'GO','Position',           [10 1 25 3],'Callback', @buttonCallback);
    hu_save =           uicontrol('parent', hp_config, 'style','pushbutton','units','characters', 'TooltipString', 'Run the function','tag', 'hu_info','String', 'save default','Position', [40 1 25 3],'Callback', @buttonCallback);
    waitfor(h_fig);

end
function buttonCallback(src,evt)
      UI.inter_mode=get(hu_inter, 'Value'); 
      UI.avg2_mode=get(hu_avg2, 'Value');       
      UI.pca_mode=get(hu_pca, 'Value');
      UI.tmip_mode=get(hu_tmip, 'Value');
      UI.avg_mode=get(hu_avg, 'Value');
      
      UI.rigid_mode=get(hu_rigid, 'Value');
      UI.Nrigid_mode=get(hu_Nrigid, 'Value');
      UI.mask_mode=get(hu_mask, 'Value');
      UI.roi_mode=get(hu_roi, 'Value');
      UI.trace_mode=get(hu_trace, 'Value');
      
      UI.DTI_mode=get(hu_dti, 'Value');
      UI.ADC_mode=get(hu_adc, 'Value');
      UI.IVIM_mode=get(hu_ivim, 'Value');
      UI.gif_mode=get(hu_gif, 'Value');
      UI.mosa_mode=get(hu_mosa, 'Value');
       
        if src==hu_valid
           close(h_fig);
        end
        if src==hu_save                    
            save_default(UI);
        end
end
function  save_default(UI)
   [folder, name, ext] = fileparts(which('UIDiffRecon_KM'));
   save ([folder '/' 'defaultUI.mat'],'UI'); 
end
function [UI]=load_default()

         UI.inter_mode = true;
         UI.avg2_mode = false; 
         UI.pca_mode = false;
         UI.tmip_mode = false;
         UI.avg_mode = true;
         UI.pca_min_ernegy = 80; 
         
         UI.rigid_mode = true;
         UI.Nrigid_mode = false;
         UI.mask_mode = true;
         UI.roi_mode = false;
         UI.trace_mode = true;
             
         UI.DTI_mode = false;
         UI.ADC_mode = false;
         UI.IVIM_mode = false;
         UI.mosa_mode = true;
         UI.gif_mode = true;
          
         
    [folder, name, ext] = fileparts(which('UIDiffRecon_KM'));
    if exist([folder '/' 'defaultUI.mat'], 'file') == 2
       d=load([folder '/' 'defaultUI.mat']);
       UI=d.UI;
    else
        save_default(UI); 
    end
end
end