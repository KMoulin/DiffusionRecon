function [Image,Infoimage]=opendicomsequence(filename,path)

% Function that open a specified DICOM file sequence, giving a 3D matrix as output. (filename not mandatory)
% 
%         [I,Info]=opendicom;
%         [I,Info]=opendicom(filename,path)
%
%  filename must be a string, with the path of the image. Info is the
%  dicominfo referring to the last image of the sequence

currdir = cd;

% Input check
if nargin == 1
    error('Input must be Filename AND Path, or nothing');
elseif nargin == 0   %no input, UI interface for file selection needed
    [filename,path] = uigetfile('*.*','load first image of the sequence');  % Image 1 name and path loading
elseif nargin > 2
    error('Too many input');
end

%prepare the files names and counting the files
nrooth=sprintf('%s*',filename(1:end-3));      % rooth of the image name sequence
cd (path);                                    % mooving on the image 1 folder
lsfiles=ls(nrooth);
nfile= size(lsfiles,1);                       % number of the files present in the folder


%reeding the first image
filenamelong = sprintf('%s%s',path,filename);
Infoimage=dicominfo(filenamelong);
II=dicomread(Infoimage);
Image=zeros(size(II,1),size(II,2),nfile);
Image(:,:,1)=II;
clear II; 

% DICOM sequence reading loop
for i=2:nfile
    clear Infoimage;clear filenamelong;
    filenamelong = sprintf('%s%s',path,lsfiles(i,:));
    Infoimage=dicominfo(filenamelong);
    Image(:,:,i)=dicomread(Infoimage);
    
end

cd(currdir);
end
