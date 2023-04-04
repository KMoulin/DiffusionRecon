function mat2=Extract_SD_KM(mat)

    mat2=[];
    for cpt1=1:1:size(mat,3)
       for cpt2=1:1:size(mat,4)
          for cpt3=1:1:size(mat,5)
              tmp=squeeze(mat(:,:,cpt1,cpt2,cpt3));
              mat2(cpt1,cpt2,cpt3)=std(tmp(tmp~=0),1);
          end
       end
    end        

end