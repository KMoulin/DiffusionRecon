function [T1 S0]= T1fit_KM3(Img,TIs,TR)

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


%
%
% S2=exp(-TE2/T2) 
% 
% S1=exp(-TE1/T2)
% 
% S2/S1=exp(-TE2/T2 + TE1/2)
% 
% log(S2/S1) = (TE1-TE2)/T2 
% 
% T2=(TE1-TE2)/log(S2/S1)
    
T2=[];
S0=[];
disp('T1 fit') 
h = waitbar(0,'T1 fit...');
for y=1:1:size(Img,1)
    for x=1:1:size(Img,2)
        if Img(y,x,:)~=0
            [tmpS0 tmpT1]=local_T2_fit(squeeze(Img(y,x,:))',TIs,TR);
            S0(y,x)=tmpS0;
            T1(y,x)=tmpT1;      
        else
            S0(y,x)=0;
            T1(y,x)=0;      
        end
    end
    waitbar(y/size(Img,1),h);
end
    
    %%
%     for cpt_slc=1:1:size(enum.slc,2)                 
%        T2(:,:,cpt_slc)=(enum.TE(1)-enum.TE(2))./log(squeeze(S2(:,:,cpt_slc,1,1,1))./squeeze(S1(:,:,cpt_slc,1,1,1)));
%        S0(:,:,cpt_slc,1)=squeeze(S1(:,:,cpt_slc,1,1,1))./squeeze(exp(-enum.TE(1)./T2(:,:,cpt_slc)));
%        S0(:,:,cpt_slc,2)=squeeze(S2(:,:,cpt_slc,1,1,1))./squeeze(exp(-enum.TE(2)./T2(:,:,cpt_slc)));
%        waitbar(cpt_slc/size(enum.slc,2),h);
%     end
%     
%     T2(isinf(T2))=0;
%     T2(isnan(T2))=0;

    close(h);    
end

function [S0 T1]=local_T2_fit(S,TI,TR)

  warning off;
  lb =   [1  1];    
  ub =   [6000   4000];     
  val0 = [1000   1000];   
  options=optimset('lsqcurvefit');
  options.Algorithm='levenberg-marquardt';
  options.TolFun=1e-8;
  options.TolX=1e-8;
  options.Display='off';
  %  options.MaxIter=1;
  %options = optimset('Display','off','DiffMinChange',1e-4); %'Algorithm','levenberg-marquardt',
  % OPTIONS = optimoptions(), 
  
 
 T1_relax = @(x,ydata)( abs(x(1).*( 1 - 2 .* exp(-ydata./x(2) ) + exp(-TR./x(2)))) ); %+ exp( -TR ./ x(2)) )
 
  [val,resnorm,residual,exitflag,output] = lsqcurvefit(T1_relax, val0, TI, S, lb, ub, options); % lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
  %exitflag
  %output.message
  
  
  % exitflag  
% 1 Function converged to a solution x.
% 
% 2 Change in x was less than the specified tolerance.
% 
% 3 Change in the residual was less than the specified tolerance.
% 
% 4 Magnitude of search direction was smaller than the specified tolerance.
% 
% 0 Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.
% 
% -1 Output function terminated the algorithm.
% 
% -2 Problem is infeasible: the bounds lb and ub are inconsistent. 

  S0=val(1);
  T1=val(2);
  
  
  warning on;
end