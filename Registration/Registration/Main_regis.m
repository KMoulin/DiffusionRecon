clc;clear;close all;
FolderName = uigetdir; % enter the folder where the native images are stored
cd(FolderName)
listing = dir(FolderName);

for cpt=3:1:size(listing,1)
    if listing(cpt).name(end-2:end) == 'dcm'
        
       filename = listing(cpt).name;
       DCMinfo = dicominfo(filename);
       tt = dicomread(filename);

       % retrieve bval/direction
       bval_str = DCMinfo.SequenceName;
       if isempty(find(bval_str =='#'))
           b_val = 0;
       else
           b_val = str2num(bval_str((find(bval_str=='#')+1):end));
       end
       %% careful : works with acquisition containing b=0, dir = 1-12 ->
       %% b=1->13
       b_val = b_val + 1;       
       DCMinfo.ImageComments = ['b=', num2str(b_val)];

       % retrive TD 
       % BE CAREFUL : implies that the TD is indicated in the sequence name as
       % last caracters i.e. TD630
       TD_str = DCMinfo.ProtocolName;

       ind1 = find(TD_str == 'T');
       TD = TD_str(ind1(end)+2:end);
           
%        % create a table with the TD values and corresponding index is put
%        % into the filename
       if exist('TDindex')
           if isempty(find(str2num(TD)==TDindex))
             TDindex(end+1) = str2num(TD);  
           end
       else
           TDindex(1) = str2num(TD);
       end
%        
       index = find(str2num(TD)==TDindex);     
       DCMinfo.SeriesNumber = 1000+index;


           G0{b_val,index} = tt;

           G1(:,:,b_val,index) = tt;
    end  
end
 Grep =  G1;

for cpt1 = 2:1:size(Grep,3)

     for cpt = 1:1:size(Grep,4)
        imm = Grep(:,:,cpt1,cpt);
        ValPos(cpt)= immean(imm);
    end
      FP = find(ValPos==max(ValPos));
      Fixed = Grep(:,:,cpt1,FP);
%       figure;imagescn(Fixed);
    for cpt2 = 1:1:size(Grep,4)
        Moving = Grep(:,:,cpt1,cpt2);
        [Mp,eng] = DiffRegistration(Fixed,Moving);
%         eng
        if eng<550
         DiffReg(:,:,cpt1,cpt2) =Grep(:,:,cpt1,cpt2);  
        else
         DiffReg(:,:,cpt1,cpt2) =zeros(size(Grep,1),size(Grep,2)); 
        end
    end
        DiffReg(:,:,1,:)= Grep(:,:,1,:);% b0
end
  save('DiffReg_s22.mat','DiffReg');
  

DiffReg =DiffReg(:,:,:,1:9);

for cpt1 = 1:1:size(DiffReg,3)
    image1 = DiffReg(:,:,cpt1,1);
    image2 = DiffReg(:,:,cpt1,2);
    for cpt2 = 1:1:size(DiffReg,4)-2

        imF = SWTimfuseTest_meanLowF(image1,image2);
        image1 = double(imF);
        image2 =DiffReg(:,:,cpt1,cpt2+2);
    end

     [ imF  ]   =  nonPCA_denoising( imF,imF,20,'normal');

    MCres(:,:,cpt1) = uint16(imF);
end
save('MCres_s2.mat','MCres');




