function [Dist] = Hist_KM (varargin)
   

     narginchk(1,3);
    if numel(varargin) == 1
          Matrix=varargin{1};
          edges = [-90:10:90];
          Matrix2=[];
    elseif numel(varargin) == 1
         Matrix=varargin{1};
         edges = varargin{2};
         Matrix2=[];
    else
         Matrix=varargin{1};
         edges = varargin{2};
         Matrix2=varargin{3};
    end
        
    mat2=[];
    Matrix(isnan(Matrix))=0;
    Matrix(isinf(Matrix))=0;
    Dist=0;
    for cpt1=1:1:size(Matrix,3)
       figure (cpt1)
       for cpt2=1:1:size(Matrix,4)
              tmp=squeeze(Matrix(:,:,cpt1));
              ttt=tmp(tmp~=0);
              %Dist(cpt1,cpt2,:)=tmp(tmp~=0);
              subplot(size(Matrix,4),1,cpt2);
              
              h1 = histogram(ttt,edges,'FaceColor',[1,0,0], 'EdgeColor',[1.0,1.0,1.0]);
              
              if isempty(Matrix2)
              else
                hold on    
                tmp=squeeze(Matrix2(:,:,cpt1));
                 ttt=tmp(tmp~=0);
                %Dist(cpt1,cpt2,:)=tmp(tmp~=0);
                 subplot(size(Matrix2,4),1,cpt2);
              
                 h1 = histogram(ttt,edges,'FaceColor',[0,0,1], 'EdgeColor',[1.0,1.0,1.0]);
              end
                  
             
       end
    end        

end