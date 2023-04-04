function [fMap,DstarMap,DMap,S0Map, ErrorMap]= IVIM_KM(Dcm,  b_vect, Mask)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [fMap,DstarMap,DMap,S0Map, ErrorMap]= IVIM_KM(Dcm,  b_vect, Mask)
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values]
%
%           b_vect - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           Mask - Mask matrix
%                  [y x slices]
%          
% OUTPUTS:  fMap - fraction of perfusion image matrix 
%                (units [%]) [y x slices]
%
%           DstarMap - Coefficient of pseudo perfusion image matrix 
%               (units [mm²/s])  [y x slices]
%           
%           DMap - Coefficient of diffusion image matrix 
%               (units [mm²/s])  [y x slices]
%
%           S0Map - Mean Diffusivity from Tensor
%                 [y x slices]
%
%           ErrorMap - Fraction of anysotropy from Tensor
%                 [y x slices]
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    fMap=zeros(size(Dcm,1),size(Dcm,2),size(Dcm,3));
    DstarMap=zeros(size(Dcm,1),size(Dcm,2),size(Dcm,3));
    S0Map=zeros(size(Dcm,1),size(Dcm,2),size(Dcm,3));
    ErrorMap=zeros(size(Dcm,1),size(Dcm,2),size(Dcm,3));
    disp('IVIM fit') 
   
    
    warning off;
    for slc=1:1:(size(Dcm,3))
        tic
        h = waitbar(0,['IVIM fit slc' num2str(slc) ]);
        for y=2:1:(size(Dcm,1)-1)                          
            for x=2:1:(size(Dcm,2)-1)
                        if Mask(y,x,slc)>0
                            Dat=squeeze(((Dcm(y,x,slc,:))));
                           % Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,slc,:))) );
                            [D, Dstar, f, ResM0]=  fitIVIM_local(Dat,b_vect');                      
                            S0Map(y,x,slc)=0;
                            fMap(y,x,slc)=f;
                            DstarMap(y,x,slc)=Dstar;
                            DMap(y,x,slc)=D;  
                            ErrorMap(y,x,slc,:)=0;     
                        else
                            S0Map(y,x,slc)=0;
                            fMap(y,x,slc)=0;
                            DstarMap(y,x,slc)=0;
                            DMap(y,x,slc)=0;
                            ErrorMap(y,x,slc,:)=0;  
                        end        
            end
            waitbar(y/size(Dcm,1),h);
        end
        close(h);
        toc
    end
     
    warning on;
end

function [D, Dstar, f, ResM0]=  fitIVIM_local(variables,b_vect) % [S0 S1 S2 S3 ... ] [b-val0 b-val1 b-val2 .... ]
   
  S0=0;
  f=0;
  D=0;
  ResM0=0;
  
  lb = [0 0 0];               % Lower bounds
  ub = [1 1 1];            % Uper bounds
  val0 = [0.1 0.001 0.01];  % Init values  
  
  options = optimset('Display','off');
  
  IVIM_M1_relax = @(x,xdata)( ( (1-x(1))*exp(-xdata*x(2)) + x(1)*exp(-xdata*x(3)) ) ); % S0 * ( (1-f)*exp(-b*D) + f*exp(-b*Dstar) )

  [val,resnorm,residual,exitflag] = lsqcurvefit(IVIM_M1_relax, val0,b_vect, variables, lb, ub, options); % lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
  

  f=val(1);
  D=val(2);
  Dstar=val(3);
  ResM0=residual;

end