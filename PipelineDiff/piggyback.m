%% Piggyback image
function [I,lim] = piggyback(I,scale)

    if nargin<2; scale = 2; end; 
    
    Ip  = zeros(ceil(size(I)*scale));
    lim = bsxfun(@plus, floor(size(I)*(scale-1)/2), [[1 1];size(I)]); 
    Ip(lim(1):lim(2),lim(3):lim(4)) = I;                              
    I = Ip;

end
