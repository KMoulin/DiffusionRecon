%% Register two images
function [Mp,sx,sy,ux,uy,e] = registerLinear_KM(F,M,opt)

   
    [F,lim] = piggyback(F,1.2);
    [M,lim] = piggyback(M,1.2);
    
  
    sx = zeros(size(M)); % deformation field
    sy = zeros(size(M));
    if isfield(opt,'sx') && isfield(opt,'sy')
        sx = piggyback(opt.sx,1.2);
        sy = piggyback(opt.sy,1.2);
    end
    e_min = 1e+100;     
    
  
    for iter=1:opt.niter

       
        [ux,uy] = findupdate(F,M,sx,sy,opt.sigma_i,opt.sigma_x);

      
        ux    = imgaussian(ux,opt.sigma_fluid);
        uy    = imgaussian(uy,opt.sigma_fluid);


        if opt.diffeomorphic
            [ux,uy] = expfield(ux,uy);
        end

        step  = 1;
        
        [cx,cy] = compose(sx,sy,step*ux,step*uy);

        sx = imgaussian(cx,opt.sigma_diffusion);
        sy = imgaussian(cy,opt.sigma_diffusion);

        e      = energy(F,M,sx,sy,cx,cy,opt.sigma_i,opt.sigma_x);

        if e<e_min
            ux_min = ux; uy_min = uy; % update best fields
            sx_min = sx; sy_min = sy;
            e_min  = e;
        else
            break;
        end


    end
    
    
    Mp = iminterpolateLinear(M,sx_min,sy_min);
    
    
    Mp = Mp(lim(1):lim(2),lim(3):lim(4));
    ux = ux(lim(1):lim(2),lim(3):lim(4));
    uy = uy(lim(1):lim(2),lim(3):lim(4));
    sx = sx(lim(1):lim(2),lim(3):lim(4));
    sy = sy(lim(1):lim(2),lim(3):lim(4));

end