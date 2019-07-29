
%% Find update between two images
function [ux,uy] = findupdate(F,M,sx,sy,sigma_i,sigma_x)

    % Interpolate updated image
    M_prime = iminterpolate(M,sx,sy); % intensities at updated points
    
    % image difference
    diff = F - M_prime;
    
    % moving image gradient
    [gy,gx] = gradient(M_prime);   % image gradient
    normg2  = gx.^2 + gy.^2;       % squared norm of gradient
    area    = size(M,1)*size(M,2); % area of moving image
    
 
    scale = diff ./ (normg2 + diff.^2*sigma_i^2/sigma_x^2);
    scale(normg2==0) = 0;
    scale(diff  ==0) = 0;
    ux = gx .* scale;
    uy = gy .* scale;
    
    % Zero non overlapping areas
    ux(F==0)       = 0; 
    uy(F==0)       = 0;
    ux(M_prime==0) = 0; 
    uy(M_prime==0) = 0;

end