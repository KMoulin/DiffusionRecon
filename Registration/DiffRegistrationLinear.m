% Demons Registration

function [Mp,eng] = DiffRegistrationLinear(F,M)

%  figure(1); clf; colormap gray;

%% Parameters
niter           = 250;
sigma_fluid     = 3.0; % regularize update      field
sigma_diffusion = 1.0; % regularize deformation field
sigma_i         = 1.0; % weight on similarity term
sigma_x         = 1.0; % weight on spatial uncertainties (maximal step)
diffeomorphic   = 1;   % use exp(u)
nlevel          = 2;   % multiresolution
do_display      = 1;   % display iterations

%% Load fixed image
% % load('F:\Matlab_HJ\fastMBD\registration\demons\demons2d\data\MOCO.mat');
% % F=moco(:,:,2,1);

if nlevel == 1
    
    %% Register
    disp(['Register...']);
    opt = struct('niter',niter, 'sigma_fluid',sigma_fluid, 'sigma_diffusion',sigma_diffusion, 'sigma_i',sigma_i, 'sigma_x',sigma_x, 'diffeomorphic',diffeomorphic, 'do_display',do_display, 'do_plotenergy',1);
    [Mp,sx,sy,ux,uy] = registerLinear(F,M,opt);

else
    
    %% Multiresolution
    sx = zeros(size(M)); % deformation field
    sy = zeros(size(M));
    for k=nlevel:-1:1
%         disp(['Register level: ' num2str(k) '...']);

        % downsample
        scale = 2^-(k-1);
        Fl = imresize(F,scale);
        Ml = imresize(M,scale);
        sxl = imresize(sx*scale,scale);
        syl = imresize(sy*scale,scale);

        % register
        opt = struct('niter',niter,...
            'sigma_fluid',sigma_fluid,...
            'sigma_diffusion',sigma_diffusion,...
            'sigma_i',sigma_i,...
            'sigma_x',sigma_x,...
            'diffeomorphic',diffeomorphic,...
            'sx',sxl, 'sy',syl,...
            'do_display',do_display, 'do_plotenergy',1);
        [Mp,sxl,syl,uxl,uyl,eng] = registerLinear(Fl,Ml,opt);
   
        % upsample
        sx = imresize(sxl/scale,size(M));
        sy = imresize(syl/scale,size(M));
    end
    
end

end


