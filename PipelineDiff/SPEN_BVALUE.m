%% B-value %%

z=[-100:1:100].*1e-3; % m

Gamma=42.5e6; % Hz/T
Gamma2=2*pi*42.5e6; % Hz/T
Ge=8e-3; % T/m
De=8e-3; % ms
dO=120e3; % Hz;
DELTA=15e-3; % ms

B=[];
T902=[];

for cpt=1:1:4
    
dO=cpt*3e4;
R=(dO)/De;

T90=(Gamma*Ge*z+dO/2)/R;
B(cpt,:)=(Gamma2^2)*(Ge^2)*(    (T90.^3)/3 + (T90.^2)*DELTA + (T90-De).*((2.*T90 + DELTA).^2) ...
-(2.*T90 + DELTA).*( (DELTA+De)^2 - (DELTA+T90).^2)+ ((DELTA+De)^3)/3 + ((DELTA+T90).^3)./3);



T902(cpt,:)=(Gamma*Ge*z+dO/2)/R;
end

B2=(Gamma2^2)*(Ge^2)*(    (T90.^3)/3 + (T90.^2)*DELTA );

B(cpt+1,:)=repmat((Gamma2^2)*(Ge^2)*(De^2)*(DELTA-De/3),1,length(z));

B=B.*1e-6;
figure,plot(z,B);


Ge=19e-3; % T/m
De=2e-3; % ms
DELTA=6e-3;
(Gamma2^2)*(Ge^2)*(De^2)*(DELTA-De/3)*1e-4

(Gamma2^2)*(Ge^2)*(De^2)*(DELTA-De)*1e-4

