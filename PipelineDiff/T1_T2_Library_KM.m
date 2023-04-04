
function [Index Library]= T1_T2_Library_KM(varargin)

narginchk(1,7);
if numel(varargin) == 1
      enum=varargin{1};
       T1Range=[200 4000];
      T2Range=[20 200];   
      T1Step=100;
      T2Step=10;
      
      OffRange=[-200 200];
      OffStep=50;
      
elseif numel(varargin) == 2
      enum=varargin{1};
      T1Range=varargin{2};
      T2Range=[20 200];
      T1Step=100;
      T2Step=10;
      
      OffRange=[-200 200];
      OffStep=50;
      
elseif numel(varargin) == 3
     enum=varargin{1};
      T1Range=varargin{2};
      T2Range=varargin{3};
      T1Step=100;
      T2Step=10;
      
      OffRange=[-200 200];
      OffStep=50;
 elseif numel(varargin) == 4
     enum=varargin{1};
      T1Range=varargin{2};
      T2Range=varargin{3};
      T1Step=varargin{4};
      T2Step=10;
      
      OffRange=[-200 200];
      OffStep=50;
      
 elseif numel(varargin) == 5
     enum=varargin{1};
      T1Range=varargin{2};
      T2Range=varargin{3};
      T1Step=varargin{4};
      T2Step=varargin{5};
      
      OffRange=[-200 200];
      OffStep=50;
elseif numel(varargin) == 6
     enum=varargin{1};
      T1Range=varargin{2};
      T2Range=varargin{3};
      T1Step=varargin{4};
      T2Step=varargin{5};
      
      OffRange=varargin{6};
      OffStep=50;
elseif numel(varargin) == 7
     enum=varargin{1};
      T1Range=varargin{2};
      T2Range=varargin{3};
      T1Step=varargin{4};
      T2Step=varargin{5};  
      
      OffRange=varargin{6};
      OffStep=varargin{7};
else
    Index=[];
    Library=[];
    return;  
end


Seq.dT=10e-3;               % ms.
Seq.TE = enum.TE(1);		% ms.
Seq.TE2 = enum.TE(2);		% ms.
Seq.TR = enum.TR;	    % ms.
Seq.TRestore=enum.TRestore;   % ms


Seq.TRFill= Seq.TR - enum.TE(2) - Seq.TRestore;
Seq.TRestore_T1=Seq.TRestore+[0:enum.nRestore-1]*Seq.TRFill/enum.nRestore; %% Check here not correct !

enum.TI= Seq.TR - enum.TE(2) - Seq.TRestore_T1;

Seq.slc=enum.datasize.slc;
%Seq.nTR=size(enum.Restore,2);

Seq.nAcq=enum.datasize.slc*((enum.datasize.dir)*(enum.datasize.b-1)+1)*(enum.datasize.avg-enum.nRestore);
Seq.nRest=enum.nRestore*enum.datasize.slc;



Seq.nT1=200;
Seq.nT2=20;

Seq.T1 = T1Range(1);	% ms.
Seq.T2 = T2Range(1);  % ms
Seq.df = OffRange(1);		    % Hz off-resonance

RF.Exc   = pi/2;	    % radians.
RF.Refoc = pi;          % radians.
RF.Rest  = pi; % radians.


RR=enum.Timing.RR;
tmp=[];
tmp.cptE=0;
Mxy=[];
Index=[];

tic

disp('Generate library')

while Seq.df<=OffRange(2) 
    while Seq.T1<=T1Range(2)
         while Seq.T2<=T2Range(2) && Seq.T2<=Seq.T1                 
            tmp.cptE= tmp.cptE+1; % New entry in the library

            Index(tmp.cptE,1)= Seq.T1;
            Index(tmp.cptE,2)= Seq.T2;
            Index(tmp.cptE,3)= Seq.df;
            %%% Initialization %%%
            M=zeros(3,Seq.slc);
            M(3,:)=1;

            tmp.cpt_slc=1;
            tmp.cpt_rest=1;

             Seq.TE = enum.TE(1);
             Seq.TRestore=Seq.TRestore_T1(tmp.cpt_rest);
            %%% Dummy Scans: Steady state %%%
            for cpt=1:1:(Seq.slc*2)
                tmp.cpt=cpt;

                %[M tmp_Mxy]=SeqRestore(Seq,RF,M,RR(cpt),tmp);
                [M tmp_Mxy]=SeqRestore_local(Seq,RF,M,Seq.TR,tmp);
               % Mxy(1,tmp.cpt_slc,tmp.cptT1)=tmp_Mxy; %% Store the initial magnetization 

                tmp.cpt_slc=tmp.cpt_slc+1;
                if(tmp.cpt_slc>Seq.slc)
                    tmp.cpt_slc=1;
                end

            end    

            %%% Real Experiment: B0 + Diff %%%
            for cpt=1:1:Seq.nAcq           
                    tmp.cpt=cpt;

                    [M tmp_Mxy]=SeqRestore_local(Seq,RF,M,RR(cpt),tmp);
                    %[M tmp_Mxy]=SeqRestore(Seq,RF,M,Seq.TR,tmp);


                    Mxy(enum.Timing.slice(cpt),enum.Timing.b(cpt),enum.Timing.dir(cpt),enum.Timing.avg(cpt),tmp.cptE)=tmp_Mxy; %% Store the initial magnetization 

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
                [M tmp_Mxy]=SeqRestore_local(Seq,RF,M,RR(cpt),tmp);
                %[M tmp_Mxy]=SeqRestore(Seq,RF,M,Seq.TR,tmp);
                Mxy(enum.Timing.slice(cpt),enum.Timing.b(cpt),enum.Timing.dir(cpt),enum.Timing.avg(cpt),tmp.cptE)=tmp_Mxy; %% Store the initial magnetization            
                tmp.cpt_slc=tmp.cpt_slc+1;
                if(tmp.cpt_slc>Seq.slc)
                    tmp.cpt_slc=1;
                    tmp.cpt_rest=tmp.cpt_rest+1;
                end

            end
             Seq.T2=Seq.T2+T2Step;
         end
          Seq.T1=Seq.T1+T1Step;
          Seq.T2=T2Range(1);
    end
    Seq.T1=T1Range(1);
    Seq.T2=T2Range(1);
    Seq.df=Seq.df+OffStep;
end
toc
Library=Mxy;

end

function [M Msxy]=SeqRestore_local(Seq,RF,M,RR,tmp)

    Rflip         = throt(RF.Exc,0);
    Rrefoc        = throt(RF.Refoc,0);
    Rrestore      = throt(RF.Rest,0);

    [Ate,Bte]     = freeprecess(Seq.TE/2,Seq.T1,Seq.T2,Seq.df);
    [Atres,Btres] = freeprecess(Seq.TRestore,Seq.T1,Seq.T2,Seq.df); %  Seq.df
    [Atr,Btr]     = freeprecess(RR-(Seq.TE+Seq.TRestore),Seq.T1,Seq.T2,Seq.df); % RR
   % [Atr,Btr]     = freeprecess(Seq.TR-(Seq.TE+Seq.TRestore),Seq.T1,Seq.T2,Seq.df);  % TR
    
    
%    M = [0;0;1];
%     dT = 2.3;	% ms between pulses.
%     M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.
% 
%     for k = 2:length(rf)
% 	[A,B] = freeprecess(dT,T1,T2,freq(f));
% 	M = A*M+B;				% Propagate to next pulse.
%     	M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
%     end;
%     sig(f) = M(1)+i*M(2);
%    
%    
   
   
    %%% Slice Selective Excitation %%%
    M(:,tmp.cpt_slc) = Rflip*M(:,tmp.cpt_slc);	% Magnetization after tip.
    
    %%% Refocussing applyed to all the slices %%%
    for cpt1=1:1:Seq.slc
        M(:,cpt1) = Ate*M(:,cpt1)+Bte;        % Magnetization at TE/2.
        M(:,cpt1) = Rrefoc* M(:,cpt1);        % Magnetization after Refoc.
        M(:,cpt1) = Ate*M(:,cpt1)+Bte;        % Magnetization at TE.
        
         if(cpt1==tmp.cpt_slc)  % Image acquisition: Save Magnetization at TE
             %Msxy(tmp.cpt)     =squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             %Msz(tmp.cpt_acqu,tmp.cpt_slc,tmp.cptT1)=squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             %Msz2(tmp.cpt,tmp.cptT1)=squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             Msxy   = squeeze(complex(M(1,tmp.cpt_slc),M(2,tmp.cpt_slc)));
         end
         
        M(:,cpt1) = Atres*M(:,cpt1)+Btres;    % Magnetization at TRestore.
        M(:,cpt1) = Rrestore* M(:,cpt1);      % Magnetization after Restore.
        M(:,cpt1) = Atr*M(:,cpt1)+Btr;	      % Magnetization at TR
    end
     
        
end


function [M Msxy]=SeqRestore_local2(Seq,RF,M,RR,tmp)

    Rflip         = throt(RF.Exc,0);
    Rrefoc        = throt(RF.Refoc,0);
    Rrestore      = throt(RF.Rest,0);

    [Ate,Bte]     = freeprecess(Seq.TE/2,Seq.T1,Seq.T2,Seq.df);
    [Atres,Btres] = freeprecess(Seq.TRestore,Seq.T1,Seq.T2,Seq.df); %  Seq.df
    [Atr,Btr]     = freeprecess(RR-(Seq.TE+Seq.TRestore),Seq.T1,Seq.T2,Seq.df); % RR
    
    
    
  
   
    %%% Slice Selective Excitation %%%
    M(:,tmp.cpt_slc)=RotPulse_KM(M(:,tmp.cpt_slc),Seq.df,2);
    %M(:,tmp.cpt_slc) = Rflip*M(:,tmp.cpt_slc);	% Magnetization after tip.
    
    %%% Refocussing applyed to all the slices %%%
    
    for cpt1=1:1:Seq.slc
        M(:,cpt1) = Ate*M(:,cpt1)+Bte;        % Magnetization at TE/2.
        
        M(:,cpt1)=RotPulse_KM(M(:,cpt1),Seq.df,4);
        %M(:,cpt1) = Rrefoc* M(:,cpt1);        % Magnetization after Refoc.
        M(:,cpt1) = Ate*M(:,cpt1)+Bte;        % Magnetization at TE.
        
         if(cpt1==tmp.cpt_slc)  % Image acquisition: Save Magnetization at TE
             %Msxy(tmp.cpt)     =squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             %Msz(tmp.cpt_acqu,tmp.cpt_slc,tmp.cptT1)=squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             %Msz2(tmp.cpt,tmp.cptT1)=squeeze(M(1,tmp.cpt_slc)+i*M(2,tmp.cpt_slc));
             Msxy   = squeeze(complex(M(1,tmp.cpt_slc),M(2,tmp.cpt_slc)));
         end
         
        M(:,cpt1) = Atres*M(:,cpt1)+Btres;    % Magnetization at TRestore.
        M(:,cpt1)=RotPulse_KM(M(:,cpt1),Seq.df,4);
       % M(:,cpt1) = Rrestore* M(:,cpt1);      % Magnetization after Restore.
        M(:,cpt1) = Atr*M(:,cpt1)+Btr;	      % Magnetization at TR
    end
     
        
end
