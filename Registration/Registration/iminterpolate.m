%% Interpolate image
function I = iminterpolate(I,sx,sy)

   
    [x,y] = ndgrid(0:(size(I,1)-1), 0:(size(I,2)-1)); 
    x_prime = x + sx; 
    y_prime = y + sy; 
    
 
    I = interpn(x,y,I,x_prime,y_prime,'linear',0); 

end
