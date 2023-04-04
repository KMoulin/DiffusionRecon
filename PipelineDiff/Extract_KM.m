function mat2=Extract_KM(mat)
    mat(isnan(mat))=0;
    mat(isinf(mat))=0;
    mat2=[];
    for cpt1=1:1:size(mat,3)
       for cpt2=1:1:size(mat,4)
          for cpt3=1:1:size(mat,5)
              tmp=squeeze(mat(:,:,cpt1,cpt2,cpt3));
              mat2(cpt1,cpt2,cpt3)=mean(tmp(tmp~=0),1);
          end
       end
    end        

end