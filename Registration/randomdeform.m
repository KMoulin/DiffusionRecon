%% Random deformation
function [I,sx,sy] = randomdeform(I,maxdeform,sigma)

    if nargin<2; maxdeform = 3; end;
    if nargin<3; sigma     = 3; end; 
    

    s  = RandStream.create('mt19937ar','seed',1); RandStream.setDefaultStream(s); % always same random numbers
    sx = maxdeform * 2*(rand(size(I))-0.5);
    sy = maxdeform * 2*(rand(size(I))-0.5);
    

    sx = imgaussian(sx,sigma);
    sy = imgaussian(sy,sigma);
    

    I = iminterpolate(I,sx,sy);

end