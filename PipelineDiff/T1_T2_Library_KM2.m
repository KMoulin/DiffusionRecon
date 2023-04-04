
function [Index Library]= T1_T2_Library_KM2(varargin)

% Bloch Equation Simulation, Excercise B-1b
% -----------------------------------------
% 
%clear all;
%close all;

%%



narginchk(1,3);
if numel(varargin) == 1
      enum=varargin{1};

% elseif numel(varargin) == 2
%       listing=varargin{1};
%       enum=varargin{2}(1);
%       dataset_num=1;
% else
%       listing=varargin{1};
%       enum=varargin{2}(1);
%       dataset_num=varargin{3}(1);
end
warning off
GammaH=42.576e6;

Seq.f = [-200:10:200];		    % Hz off-resonance.

Seq.dT=100e-3;               % ms.
Seq.TE = enum.TE(1);		% ms.
Seq.TE2 = enum.TE(2);		% ms.
Seq.TR = enum.TR;	    % ms.
Seq.TRestore=158;   % ms


Seq.TRFill= Seq.TR - enum.TE(2) - Seq.TRestore;
Seq.TRestore_T1=Seq.TRestore+[0:enum.nRestore-1]*Seq.TRFill/enum.nRestore;

enum.TI= Seq.TR - enum.TE(2) - Seq.TRestore_T1

Seq.slc=enum.datasize.slc;
%Seq.nTR=size(enum.Restore,2);

Seq.nAcq=enum.datasize.slc*((enum.datasize.dir)*(enum.datasize.b-1)+1)*(enum.datasize.avg-enum.nRestore);
Seq.nRest=enum.nRestore*enum.datasize.slc;



Seq.nT1=200;
Seq.nT2=1;



RF.Exc   = 90;	    % degree.
RF.Refoc = 180;     % degree.
RF.Rest  = 180;     % degree.

RF.Duration=2000e-3;                                              % us -> ms
RF.Points=ceil(RF.Duration/Seq.dT);                          % Number of point of the RF pulse                                                 % degree
RF.Pbwt=5200;                                                     % Bwt per product 5200 Hz %DFz=Gs*Thickness*GammaH; % Bandwith (Check)



RF.DurationRefoc=4300e-3;                                         % us -> s
RF.PointsRefoc=ceil(RF.DurationRefoc/Seq.dT);                % Number of point of the RF pulse                                            
RF.PbwtRefoc=5200;                                                 % Bwt per product 5200 Hz %DFz=Gs*Thickness*GammaH; % Bandwith (Check)


%%1gauss =100 ?T

RF.Shape90=[];
RF.Shape90=RF_generator_local(RF.Points,Seq.dT,RF.Pbwt,1,0,0,0)   ;
VecT=[1:1:length(RF.Shape90)].*Seq.dT*1e-3;
DD=GammaH*trapz(VecT,abs(RF.Shape90)) %%%% T/S *   Hz/T
RF.Shape90=RF.Exc*RF.Shape90*40/DD;

RF.Shape180=[];
RF.Shape180=RF_generator_local(RF.PointsRefoc,Seq.dT,RF.Pbwt,1,0,0,0)  ;
VecT=[1:1:length(RF.Shape180)].*Seq.dT*1e-3;
DD=GammaH*trapz(VecT,abs(RF.Shape180)) %%%% T/S *   Hz/T
RF.Shape180=RF.Refoc*RF.Shape180*50/DD;

RF.Shape180Restore=[];
RF.Shape180Restore=RF_generator_local(RF.PointsRefoc,Seq.dT,RF.Pbwt,1,0,0,0)  ;
VecT=[1:1:length(RF.Shape180Restore)].*Seq.dT*1e-3;
DD=GammaH*trapz(VecT,abs(RF.Shape180Restore)) %%%% T/S *   Hz/T
RF.Shape180Restore=RF.Rest*RF.Shape180Restore*50/DD;





RR=enum.Timing.RR;
tmp=[];
tmp.cptT1=0;
Mxy=[];
Index=[];

tic

disp('Generate library') 
h = waitbar(0,['Generate library with ' num2str(Seq.nT1) ' T1s']);
 for cptT1=1:1:Seq.nT1
     for cptT2=1:1:min(Seq.nT2,cptT1)    
        
        
        tmp.cptT1= tmp.cptT1+1; % New entry in the library
        
        Seq.T1 = (10*cptT1+500)*1e-3;	% ms.
        Seq.T2 = 50*1e-3;%5*cptT2+25;  % ms
     
        Index(tmp.cptT1,1)= Seq.T1;
        Index(tmp.cptT1,2)= Seq.T2;

        %%% Initialization %%%
        M=zeros(3,length(Seq.f),Seq.slc);
        M(3,:,:)=1;

        tmp.cpt_slc=1;
        tmp.cpt_rest=1;

         Seq.TE = enum.TE(1);
         Seq.TRestore=Seq.TRestore_T1(tmp.cpt_rest);
        %%% Dummy Scans: Steady state %%%
        for cpt=1:1:(Seq.slc*2)
            tmp.cpt=cpt;

            %[M tmp_Mxy]=SeqRestore(Seq,RF,M,RR(cpt),tmp);
            [M tmp_Mxy]=SeqRestore_Local2(Seq,RF,M,Seq.TR,tmp);
           % Mxy(1,tmp.cpt_slc,tmp.cptT1)=tmp_Mxy; %% Store the initial magnetization 

            tmp.cpt_slc=tmp.cpt_slc+1;
            if(tmp.cpt_slc>Seq.slc)
                tmp.cpt_slc=1;
            end

        end    
       
        %%% Real Experiment: B0 + Diff %%%
        for cpt=1:1:Seq.nAcq           
                tmp.cpt=cpt;

                [M tmp_Mxy]=SeqRestore_Local2(Seq,RF,M,RR(cpt),tmp);
                %[M tmp_Mxy]=SeqRestore(Seq,RF,M,Seq.TR,tmp);

              
                Mxy(enum.Timing.slice(cpt),enum.Timing.b(cpt),enum.Timing.dir(cpt),enum.Timing.avg(cpt),tmp.cptT1,:)=tmp_Mxy; %% Store the initial magnetization 
               
                tmp.cpt_slc=tmp.cpt_slc+1;
                if(tmp.cpt_slc>Seq.slc)
                    tmp.cpt_slc=1;                  
                end
        end    

        %%% Real Experiment:Restore FingerPrinting %%%
        cpt_end=tmp.cpt;
        %tmp.cpt_slc=1;
        Seq.TE = enum.TE(2);
        
        for cpt=cpt_end+1:1:Seq.nRest+ cpt_end
            tmp.cpt=cpt;  
            
            Seq.TRestore=Seq.TRestore_T1(tmp.cpt_rest);
            [M tmp_Mxy]=SeqRestore_Local2(Seq,RF,M,RR(cpt),tmp);
            %[M tmp_Mxy]=SeqRestore(Seq,RF,M,Seq.TR,tmp);
            Mxy(enum.Timing.slice(cpt),enum.Timing.b(cpt),enum.Timing.dir(cpt),enum.Timing.avg(cpt),tmp.cptT1,:)=tmp_Mxy; %% Store the initial magnetization            
            tmp.cpt_slc=tmp.cpt_slc+1;
            if(tmp.cpt_slc>Seq.slc)
                tmp.cpt_slc=1;
                tmp.cpt_rest=tmp.cpt_rest+1;
            end

        end
     end
  waitbar(cptT1/Seq.nT1,h);  
 end
close(h);
warning on
toc
Library=Mxy;

end

function [M Msxy]=SeqRestore_Local2(Seq,RF,M,RR,tmp)
     
    PointsTE=Seq.TE/Seq.dT;
    
    TotalPointsTE=PointsTE+RF.Points/2;
    TotalPointsTR=ceil((RR)/Seq.dT-TotalPointsTE);
    
    StartRefoc=RF.Points/2+PointsTE/2-RF.PointsRefoc/2;
    StartRestore=ceil(Seq.TRestore/Seq.dT);
    
    
    bTE = zeros(int32(TotalPointsTE),1);
    bTR = zeros(int32(TotalPointsTR),1);
    gTE = zeros(int32(TotalPointsTE),1);
    gTR = zeros(int32(TotalPointsTR),1);
        
   % 0 -> TE
   bTE = zeros(int32(TotalPointsTE),1);   
   for cpt1=1:1:Seq.slc
       if(cpt1==tmp.cpt_slc)            
           % 90 +180
            bTE(1:1+RF.Points)=complex(real(RF.Shape90),imag(RF.Shape90));
            bTE(StartRefoc:StartRefoc+RF.PointsRefoc)=complex(real(RF.Shape180),imag(RF.Shape180));
       else
           % 180
           bTE(StartRefoc:StartRefoc+RF.PointsRefoc)=complex(real(RF.Shape180),imag(RF.Shape180));
       end
       [mx,my,mz] = bloch(bTE,gTE,Seq.dT*1e-3,Seq.T1,Seq.T2,Seq.f,[0 0 0],0,M(1,:,tmp.cpt_slc),M(2,:,tmp.cpt_slc),M(3,:,tmp.cpt_slc)); 
       M(:,:,tmp.cpt_slc)=[mx;my;mz];  
   end 
   Msxy   = complex(M(1,:,tmp.cpt_slc),M(2,:,tmp.cpt_slc));  
   
   % TE -> RR
   bTR = zeros(int32(TotalPointsTR),1);
    for cpt1=1:1:Seq.slc
        bTR(StartRestore:StartRestore+RF.PointsRefoc)=complex(real(RF.Shape180Restore),imag(RF.Shape180Restore));
       [mx,my,mz] = bloch(bTR,gTR,Seq.dT*1e-3,Seq.T1,Seq.T2,Seq.f,[0 0 0],0,M(1,:,tmp.cpt_slc),M(2,:,tmp.cpt_slc),M(3,:,tmp.cpt_slc)); 
       M(:,:,tmp.cpt_slc)=[mx;my;mz];
   end
     
end

function RF=RF_generator_local(sample_size,dT,DFz,RFAmpl,Shift,grad,Phase)   


   
    RF=[];
    j=1;
    dT=dT*1e-3;
    for k= (-sample_size/2):1:(sample_size/2)
           % h=0.54-0.46*cos(2*pi*(k*dT)/(sample_size*dT));         % hamming window
            h=1/2 * (1 + cos(2*pi* (k*dT) / (sample_size*dT)));     % han window
            u=DFz*k*dT*pi;                                          % U= Bandwith*t*pi
            if u ~= 0
                RF(j)=complex(h*sin(u)/u,0)*RFAmpl;       
            else
                RF(j)=complex(h,0)*RFAmpl;
            end
            
          %  RF(j)=complex(abs(RF(j))*cos(pi*(Shift*grad)*267e6),abs(RF(j))*cos(pi*(Shift*grad)*267e6));
            j=j+1;
    end

    
    %b1 = [zeros(1,500) msinc(sample_size-1000,2) zeros(1,501)];
    %RF=b1;
    %Write_PTA_Pulse('KPulse.pta',TURK,Ph);

end