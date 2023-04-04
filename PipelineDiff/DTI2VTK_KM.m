function DTI2VTK_KM(Dcm, Mask,enum,EigNum,ColorInput,OutputFileName,VectSign)

%  Export the tensor as a VTK format 
%
% SYNTAX:  DTI2VTK_KM(Dcm, enum,EigNum,OutputFileName)
%
% INPUTS:   Dcm - EigVector image matrix
%                 [y x slices b-values coordinates eigVector]
%           
%           enum - Structure which contains information about the dataset 
%
%           EigNum - Number of the desired EigVector to be exported
%
%           OutputFileName - Name of the .vtk file
%
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    disp('Generate colormap VTK') 
    
%     colorVect=[];
%     k=1;
%     for x=1:1:11
%         for y=1:1:11
%             for z=1:1:11
%                    k=k+1;
%                    colorMap(k,:)=[x,y,z]/10-1;
%                    colorVect(x,y,z)=k;
%             end
%         end
%     end
%     ColormapToVTK_local(colorMap);
    
    disp('Export to VTK')
    h = waitbar(0,'Export to VTK...');
    
    ColorInput(isnan(ColorInput))=0;
    for cpt_b=1:1:(size(enum.b,2)-1)
        Nodes=[];
        Vectors=[];
        Colors=[];
        for cpt_slc=1:1:size(Dcm,3)
            for y=1:1:size(Dcm,1)
                for x=1:1:size(Dcm,2)    
                    if(squeeze(Mask(y,x,cpt_slc))~=0)
                        Tmp_Vect=VectSign.*[squeeze(Dcm(y,x,cpt_slc,cpt_b,1,EigNum)) squeeze(Dcm(y,x,cpt_slc,cpt_b,2,EigNum)) squeeze(Dcm(y,x,cpt_slc,cpt_b,3,EigNum))];
                        Vectors=[Vectors; Tmp_Vect];
                        Nodes=[Nodes; (y-size(Dcm,1)) (x-size(Dcm,2)) 5*cpt_slc]; %cpt_slc
                        if isempty(ColorInput)
                            Colors=[Colors;abs(Tmp_Vect)];
                        else
                            Colors=[Colors;ColorInput(y,x,cpt_slc)];
                        end
                        
                    end                    
                end
            end
            waitbar(cpt_slc/size(enum.slc,2),h);
        end
        VectorsToVTK_local([Nodes], [Colors], 'Color',[Vectors],'fiber', ['b' num2str(enum.b(cpt_b+1)) '_' OutputFileName]);
    end
    
    close(h);    

end

function [] = VectorsToVTK_local(Nodes, Scalars, ScalarName, Vectors, VectorName, OutputFileName)



fileID = fopen([OutputFileName,'.vtk'],'w');

fprintf(fileID, '%s\n', '# vtk DataFile Version 4.0');
fprintf(fileID, '%s\n', 'vtk output');
fprintf(fileID, '%s\n', 'ASCII');
fprintf(fileID, '%s\n', 'DATASET POLYDATA');
fprintf(fileID, '%s ' , 'POINTS');
fprintf(fileID, '%d ' , size(Nodes,1));
fprintf(fileID, '%s\n', 'float');

for k=1:size(Nodes,1)
    fprintf(fileID, '%1.6f %1.6f %1.6f \n', Nodes(k,:));
end


% if (size(Scalars, 1) == 0)
% end
fprintf(fileID, '\n%s ' , 'POINT_DATA');
fprintf(fileID, '%d \n' , size(Nodes,1) );
if (size(Vectors, 1) > 0)
%     fprintf(fileID, '\n%s ' , 'POINT_DATA');
    fprintf(fileID, '%s\n' , ['VECTORS   ',VectorName,'   double']);
    for k=1:size(Vectors,1)
        fprintf(fileID, '%1.6f %1.6f %1.6f \n', Vectors(k,:));
    end
end
[LookUpTable LookUpValue]=ColorLookupTable_local(256,Scalars);

if (size(Scalars, 1) > 0)
    fprintf(fileID, '%s \n' , ['SCALARS   ',ScalarName,'   double']);
    fprintf(fileID, '%s\n', 'LOOKUP_TABLE default');
    for k=1:size(LookUpValue,1)
        fprintf(fileID, '%1.6f \n', LookUpValue(k,:));
    end
end

ColormapToVTK_local(LookUpTable);

fclose(fileID);

end

function ColormapToVTK_local(colorVect)
  scheme='Fiber_color';
  fid = fopen([scheme '.xml'], 'w');
  fwrite(fid, '<ColorMaps>\n');
  fprintf(fid, '<ColorMap name="%s" space="RGB">\n', scheme);
   
  for k=1:size(colorVect,1)
      fprintf(fid, '  <Point x="%f" o="1" r="%f" g="%f" b="%f"/>\n', [colorVect(k,4) colorVect(k,1) colorVect(k,2) colorVect(k,3)]);
  end
  fwrite(fid, '</ColorMap>\n');
  fwrite(fid, '<ColorMaps>\n');
  fclose(fid);
end

function [LookUpTable LookUpValue]=ColorLookupTable_local(bit,color)
        LookUpTable=[];

        LookUpTable=jet(bit);
         if size(color,2)==1
             LookUpValue=abs(color./max(color))*bit;
         else
            [LookUpValue,dist] = dsearchn(LookUpTable,abs(color./max(color)));
            
         end
         LookUpTable(:,4)=[1:1:bit];
end