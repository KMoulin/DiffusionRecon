function [Mean, STD] = Stats_KM (Matrix, Mask)
    mat2=[];
    Matrix(isnan(Matrix))=0;
    Matrix(isinf(Matrix))=0;
    for cpt1=1:1:size(Mask,3)
       for cpt2=1:1:size(Mask,4)       
              tmp=squeeze(Matrix(:,:,cpt1)).*squeeze(Mask(:,:,cpt1,cpt2));        
              Mean(cpt1,cpt2)=mean(tmp(tmp~=0),1);
              STD(cpt1,cpt2)=std(tmp(tmp~=0),1);   
       end
    end        

end