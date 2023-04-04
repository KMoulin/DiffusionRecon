%% Get energy
function e = energy(F,M,sx,sy,cx,cy,sigma_i,sigma_x)

    % Intensity difference
    Mp     = iminterpolate(M,sx,sy); % KEVIN USE that
    diff2  = (F-Mp).^2;
    area   = size(M,1)*size(M,2);
    
    % Transformation Gradient
    [gx_sx,gy_sx] = gradient(sx);
    [gx_sy,gy_sy] = gradient(sy);
    
    % Three energy components
    e_sim  = sum(diff2(:)) / area;
    %e_dist = sum((cx(:)-sx(:)).^2 + (cy(:)-sy(:)).^2) / area;
    e_reg  = sum(gx_sx(:).^2 + gy_sy(:).^2) / area;
    
    % Total energy
    e      = e_sim + (sigma_i^2/sigma_x^2) * e_reg;

end