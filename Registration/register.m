%% Register two images
function [Mp,sx,sy,ux,uy,e] = register(F,M,opt)

    %% Piggyback image
    [F,lim] = piggyback(F,1.2);
    [M,lim] = piggyback(M,1.2);
    
    %% T is the deformation from M to F
    sx = zeros(size(M)); % deformation field
    sy = zeros(size(M));
    if isfield(opt,'sx') && isfield(opt,'sy')
        sx = piggyback(opt.sx,1.2);
        sy = piggyback(opt.sy,1.2);
    end
    e_min = 1e+100;      % Minimal energy
    
    %% Iterate update fields
    for iter=1:opt.niter

        % Find update
        [ux,uy] = findupdate(F,M,sx,sy,opt.sigma_i,opt.sigma_x);

        % Regularize update
        ux    = imgaussian(ux,opt.sigma_fluid);
        uy    = imgaussian(uy,opt.sigma_fluid);

        % Get diffeomorphic update
        if opt.diffeomorphic
            [ux,uy] = expfield(ux,uy);
        end

        % Compute step (max half a pixel)
        step  = 1;
        
        % Update correspondence (demons) - additive
        %cx = sx + step*ux;
        %cy = sy + step*uy;

        % Update correspondence (demons) - composition
        [cx,cy] = compose(sx,sy,step*ux,step*uy);
        
        % Regularize deformation
        sx = imgaussian(cx,opt.sigma_diffusion);
        sy = imgaussian(cy,opt.sigma_diffusion);

        % Compute energy
        e      = energy(F,M,sx,sy,cx,cy,opt.sigma_i,opt.sigma_x);
%         disp(['Iteration: ' num2str(iter) ' - ' 'energy: ' num2str(e)]);
        if e<e_min
            ux_min = ux; uy_min = uy; % update best fields
            sx_min = sx; sy_min = sy;
            e_min  = e;
        else
            break;
        end

%         if opt.do_display
%             figure(1);
%             
%             % display deformation
%             subplot(2,4,7); showvector(ux,uy,4,3,lim); title('Update');
%             subplot(2,4,8); showvector(sx,sy,4,3,lim); title('Transformation');
%             drawnow;
%             
%             % Display registration
%             Mp     = iminterpolate(M,sx,sy);
%             diff   = (F-Mp).^2;
%             showimage(F,'Fixed', M,'Moving', Mp,'Warped', diff,'Diff', 'lim',lim,'nbrows',2); drawnow;
% 
%             % Plot energy
%             if opt.do_plotenergy
%                 plot_e(iter) = e;
%                 subplot(2,2,3)
%                 hold on;
%                 plot(1:iter,plot_e,'r-'); xlim([0 opt.niter]);
%                 xlabel('Iteration'); ylabel('Energy');
%                 hold off;
%                 drawnow
%             end
%         end

    end
    
    %% Transform moving image
    Mp = iminterpolate(M,sx_min,sy_min);
    
    %% Unpiggyback
    Mp = Mp(lim(1):lim(2),lim(3):lim(4));
    ux = ux(lim(1):lim(2),lim(3):lim(4));
    uy = uy(lim(1):lim(2),lim(3):lim(4));
    sx = sx(lim(1):lim(2),lim(3):lim(4));
    sy = sy(lim(1):lim(2),lim(3):lim(4));

end