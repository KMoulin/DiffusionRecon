function [ADC] = ADCMap2B_KM(Vol,Vol2,b_vect,b_vect2)

    % Ref : 60, 130

    ADC=zeros(size(Vol,1),size(Vol,2),size(Vol,3));

    for y=1:1:size(Vol,1)
        for x=1:1:size(Vol,2)
            for z=1:1:size(Vol,3)
                variables=[squeeze(VolB0(y,x)) squeeze(Vol(y,x,:))'];
                [p,S] = polyfit( b_vect,log((variables)),1);
                if p(1)>0
                    ADC(y,x) = 0; % D must not be <0 !!
                   %disp('Fits on mean, read: D<0 set to 0')
                else
                    ADC(y,x) = -p(1);
                end
                ADC(y,x,z)=log(Vol(y,x,z)/Vol2(y,x,z))/(b_vect2-b_vect);
            end
        end
    end
end