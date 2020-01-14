function [Dcm2]= Apply_Mask_KM(varargin)

        narginchk(2,3);
        if numel(varargin) == 2
              Dcm=varargin{1};
              Mask=varargin{2};
              Mask2=Mask;
        else
              Dcm=varargin{1};
              Mask=varargin{2};
              Mask2=varargin{3};
        end
        
        SMask=Mask+Mask2;
        SMask(SMask==1)=nan;
        SMask(SMask==2)=1;

        Dcm2=Dcm;
       for cpt4=1:1:size(Dcm,4)
           for cpt5=1:1:size(Dcm,5)
               for cpt6=1:1:size(Dcm,6)
                   for cpt7=1:1:size(Dcm,7)
                       for cpt8=1:1:size(Mask,4)
                        Dcm2(:,:,:,cpt4,cpt5,cpt6,cpt7,cpt8)=squeeze(Dcm(:,:,:,cpt4,cpt5,cpt6,cpt7)).*SMask(:,:,:,cpt8);
                       end
                   end
               end
           end
       end

end