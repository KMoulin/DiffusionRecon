function [T1]= T1fit_KM2(DcmB0Rest, enum ,TI)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [Dcm2 DcmB02 enum2]= Average_KM(Dcm, DcmB0, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages]
%
%           DcmB0 - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages]
%
%           DcmB02 - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           enum2 - Structure which contains information about the dataset 
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

   
    T1=[];
    disp('T1 fit') 
    h = waitbar(0,'T1 fit...');
    
    p1 = -8.7189e-17;
    p2 = 3.8419e-13;
    p3 = -4.2203e-10;
    p4 = -1.6477e-07;
    p5 = 0.00010834;
    p6 = 0.99244;
 
    %x=[1:1:2000];  
    %y = p1*x.^5 + p2*x.^4 +  p3*x.^3 + p4*x.^2 + p5*x + p6 ;
    
    DcmB0Rest(isnan(DcmB0Rest))=0;
    DcmB0Rest(isinf(DcmB0Rest))=0;

    for cpt_slc=1:1:size(enum.slc,2)
     for cpty=1:1:size(DcmB0Rest,1)
         for cptx=1:1:size(DcmB0Rest,2)

             if DcmB0Rest(cpty,cptx,cpt_slc)>0          
              lb = [0];    
              ub = [2000];     
              val0 = [800];   
              options = optimset('Display','off');

              T1_relax = @(x,xdata)(xdata - (p1*x(1).^5 + p2*x(1).^4 +  p3*x(1).^3 + p4*x(1).^2 + p5*x(1) + p6)) ;

              [val,resnorm,residual,exitflag] = lsqcurvefit(T1_relax, val0, squeeze(DcmB0Rest(cpty,cptx,cpt_slc)), 0, lb, ub, options); % lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)

              T1(cpty,cptx,cpt_slc)=val(1);
              T1_residual(cpty,cptx,cpt_slc)=residual;
             else
                 T1(cpty,cptx,cpt_slc)=0;
                 T1_residual(cpty,cptx,cpt_slc)=0;
             end

         end
     end
     waitbar(cpt_slc/size(enum.slc,2),h);
    end
 

    close(h);    
    
    
  
      


 
 

end