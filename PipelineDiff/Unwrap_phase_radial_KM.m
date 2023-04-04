function HA_final = Unwrap_phase_radial_KM(HA,Depth)

HA_final=HA.*0;

for y=2:1:size(HA,1)-1
    for x=2:1:size(HA,2)-1
        Filter=[1 1 1;
                0 0 0;
                -1 -1 -1;];
        Filter=[Depth(y+1,x-1) Depth(y+1,x) Depth(y+1,x+1);...
                Depth(y,x-1)  Depth(y,x)  Depth(y,x+1);...
                Depth(y-1,x-1) Depth(y-1,x) Depth(y-1,x+1)];
        Data_Mask=[HA(y+1,x-1) HA(y+1,x) HA(y+1,x+1);...
                   HA(y,x-1)  HA(y,x)  HA(y,x+1);...
                   HA(y-1,x-1) HA(y-1,x) HA(y-1,x+1)];
                  
        HA_final(y,x)=sum(sum(Data_Mask.*Filter));
    end
end






end