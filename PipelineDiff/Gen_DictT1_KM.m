


function [T1Dict]=Gen_DictT1_KM(enum)



%% Init %%
M=zeros(3,enum.slc);
M(3,:)=1;


%% Dummy scan to reach steady state
cpt_slc=1;
for cpt=1:1:enum.slc*2
     
    cpt_slc=cpt_slc+1;
    if(cpt_slc>enum.slc)cpt_slc=1;end
end


%% Diffusion Dataset
df = 0;		% Hz off-resonance.
T1 = 1000;	% ms.
T2 = 10*cptT1;	% ms.
TE = 56;		% ms.
TRestore=150; % ms
TR = 2000;	% ms.
nTR=20;

flip = pi/2;	% radians.
flipRefoc = pi;% radians.
flipRestore = pi;%pi;% radians.
slc=5;

M=zeros(3,slc);
M(3,:)=1;

Rflip = yrot(flip);
Rrefoc =xrot(flipRefoc);
Rrestore =xrot(flipRestore);

[Ate,Bte] = freeprecess(TE/2,T1,T2,df);
[Atres,Btres] = freeprecess(TRestore,T1,T2,df); % 
[Atr,Btr] = freeprecess(TR-(TE+TRestore),T1,T2,df); % 


cpt_slc=1;


%%





end
