% T2 mapping script (Ruud van Heeswijk 03-07-2014) for separate
% acquisitions run a 'cell' by pressing ctrl+enter after setting its data
% change easily between cells with ctrl+up/down

close all, clear all, clc
%% 1: check file names of the dicoms and find the relevant ones

% % paste the root path of the dicoms here (should be the only thing to set
% % each time in this cell)
% FilePath='E:\Ruud\Data\20140702 T2mapInfarctStEtienne\COEUR_PC_CARDIAC_DOT_ENGINE_20140702_205553_076000\';
% FilePath=uigetdir();
% % check that the path of the registeration is correct
% %addpath('C:\Users\Ruud\Dropbox\Matlab\2DImageRegistrationToolbox');
% 
% % create variables
% dirlist=dir(FilePath);
% dirfileinfo=cell(length(dirlist),1);
% dirflagtime=zeros(length(dirlist),3);
% 
% % check which folders have 1 file and CV in their name
% for counter=1:length(dirlist)
%     dirtempname=[FilePath dirlist(counter,1).name '\'];    
%     if (length(dir(dirtempname))==3 && (~isempty(strfind(dirlist(counter,1).name,'CV'))==1))
%         dirflagtime(counter,1)=counter;
%         
%         filetempname=dir(dirtempname); 
%         filetempname=filetempname(3,1).name;
%         dirfileinfo{counter,1}=filetempname;
%         
%         tempdicominfo=dicominfo([dirtempname dirfileinfo{counter,1}]);
%         dirflagtime(counter,2)=str2num(tempdicominfo.SeriesTime);
%         dirflagtime(counter,3)=1;
%     end
% end
% 
% %reorder along display time and display 
% filereorder=sortrows(dirflagtime,2);
% disp ('The numbers on the left need to be filled out in pairs of 3 in the cell below')
% for counter= 1:length(filereorder)
%     if filereorder(counter,3)==1
%         line=strcat(num2str(filereorder(counter,1)),'......', dirlist(filereorder(counter,1),1).name,1 );
%         disp(line)
%     end
% end

%% 2. read in files, with the series directory numbers in the proper order

% PUT FILE NUMBERS (from previous cell) HERE
% series=['191207'  % in order 60-30-0
%        '241510'
%        '231409'];
% [slices,~]=size(series);
% 
% % Organized in bigimage as: (TET2prepNr,Y(-row),X(+column),Slice)
% for slicevar=1:slices
%     for k=1:3  %number of images in the series
%         dirnumber=str2num(series(slicevar,2*k-1:2*k));
%         filename=dirfileinfo{dirnumber,1}; 
%         path=strcat(FilePath,dirlist(dirnumber,1).name, '\');
%         [image,infoimage]=opendicomsequence(filename,path);
%         bigimage(k,:,:,slicevar)=squeeze(image);
%     end
% end


bigimage=read_and_sort_T2();
% Reference image for later coregistration (middle of the 3)
bigRef = squeeze(bigimage(2,:,:,:)); 
% moving images
bigI1 = squeeze(bigimage(1,:,:,:));
bigI3 = squeeze(bigimage(3,:,:,:));

%% interactive crop
% crop the images to the heart only
nslices=size(bigimage,4);
figure(1)
clims=[0 max(max(bigRef(:,:,ceil(nslices/2))))/2];
imagesc(squeeze(bigRef(:,:,ceil(nslices/2))),clims), axis equal tight, colormap('gray')
% crop by drawing a rectangle around the heart
bb=imrect;
pos=round(getPosition(bb));
close 1
Ref=bigRef(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),:);
I1=bigI1(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),:);
I3=bigI3(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),:);
%check heart and phantom are in images
for kk=1:nslices
    [0 max(max(bigRef(:,:,kk)))/2];
    imagesc(squeeze(Ref(:,:,kk)),clims), axis equal tight off, colormap('gray')
    pause
end
close
%% Coregistration

% set the options of the coregistration
% (make sure the coregistration toolkit is installed!)
Options = struct('Similarity',{'mi'},...
                      'Registration',{'Affine'},...
                      'Penalty',{1e-3},...
                      'MaxRef',{2},...
                      'Grid',{[]},...
                      'Spacing',{[]},...
                      'MaskMoving',{[]},...
                      'MaskStatic',{[]},'Verbose',{0},...
                      'Points1',{[]},'Points2',{[]},'PStrength',{[]},...
                      'Interpolation',{'Linear'},'Scaling',{[1 1]});
                                   
% Similarity Options
% 'd'  : differences between I1 and I2
% 'sd' : squared differences
% 'mi' : normalized (local) mutual information      this one
% 'mip': normalized mutual information with image split in multiple small regions
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
% Registration Options:
%               Rigid    : Translation, Rotation
%               Affine   : Translation, Rotation, Shear, Resize this one
%               NonRigid : B-spline grid based local registration
%               Both     : Nonrigid and Affine (Default)
[cropY,cropX,~]=size(Ref);
Coregged1=zeros(cropY,cropX,nslices);
Coregged3=zeros(cropY,cropX,nslices);
for kk=1:nslices
    [Coregged1(:,:,kk),O_trans1,Spacing1,M1,param1,B1,F1] = image_registration (I1(:,:,kk),Ref(:,:,kk),Options);
    [Coregged3(:,:,kk),O_trans3,Spacing2,M3,param2,B3,F3] = image_registration (I3(:,:,kk),Ref(:,:,kk),Options);
end

%% bypass coregistration if wanted
Coregged1=I1;
Coregged3=I3;

%% plot the images before and after coregistration if wanted
figure(1)
set(gcf,'name','after')
clims=[0 150];
for kk=1:nslices
    clims=[0 max(max(max(bigimage(:,:,:,kk))))/2];
    subplot(231),imagesc(squeeze(Coregged1(:,:,kk)),clims),title('T2prepafter=60'), axis equal tight off, colormap('gray') 
    subplot(232),imagesc(squeeze(Ref(:,:,kk)),clims),title('T2prepafter=30'), axis equal tight off, colormap('gray')
    subplot(233),imagesc(squeeze(Coregged3(:,:,kk)),clims),title('T2prepafter=0'), axis equal tight off, colormap('gray')

    %clims=[0 30];
    subplot(234),imagesc(squeeze(I1(:,:,kk)),clims),title('T2prepbefore=60'), axis equal tight off, colormap('gray') 
    if exist('RefOr')~=0
        subplot(235),imagesc(squeeze(RefOr(:,:,kk)),clims),title('T2prepbefore=30'), axis equal tight off, colormap('gray')
    else
        subplot(235),imagesc(squeeze(Ref(:,:,kk)),clims),title('T2prepbefore=30'), axis equal tight off, colormap('gray')
    end
    subplot(236),imagesc(squeeze(I3(:,:,kk)),clims),title('T2prepbefore=0'), axis equal tight off, colormap('gray')
    pause
end

%% optional: display fit in single pixel
%clims2=[0 350];
close all
slicesel=1;  %choose the slice here
figure(9)
clims=[0 max(max(bigimage(2,:,:,slicesel)))/2];
imagesc(squeeze(Coregged1(:,:,slicesel)),clims), axis equal tight, colormap('gray')
[xp,yp]=ginput(1);
close(9)

[xs ys ~]=size(Ref);
coregset=zeros(3,xs,ys);
coregset(1,:,:)=Coregged1(:,:,slicesel);
coregset(2,:,:)=Ref(:,:,slicesel);
coregset(3,:,:)=Coregged3(:,:,slicesel);

% Fitting parameters
s = fitoptions('Method','NonlinearLeastSquares',...
                'Startpoint',[50 50]);
%                'Lower',[50,5,0],...
%                'Upper',[1000,500,500],...
%                'Startpoint',[200 50 0]);

f = fittype('a*exp(-x/b)+a*0.11','options',s);
  
time=[60 31 0]';    %time=[0 31 60]';
sig=squeeze(coregset(:,round(yp),round(xp)));
figure(2)
plot(time,sig,'o')
xlim([-1 1.1*max(time)]),ylim([0 1.1*max(sig)])
[pixt2fit,pixparms]=fit(time,sig,f)
hold on
plot(pixt2fit,'r')
xlabel('time [ms]')
ylabel('image intensity [a.u.]')

text(30,20,...
 ['Fitted T2 value = ',num2str(pixt2fit.b),' ms & R^2 = ', num2str(pixparms.adjrsquare)],...
 'HorizontalAlignment','center',... 
 'BackgroundColor',[1 1 1]);

%% calculate the T2map for the slices
%set fitting parameters
MappingThreshold=0; %set a threshold intensity value below which the T2 is set to zero 
MinumumIntensity=50; %set a minimum intensity for the fit - the higher this is the better the lungs/noise look

h = waitbar(0,'Waiting for T2mapping ...');

% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Startpoint',[200 50]);
% %                'Lower',[50,5],...
% %                'Upper',[1000,500],...
% %                'Startpoint',[200 50]);

% playing with fit boundary conditions saves time
% first variable is zero crossing (amplitude), second is T2
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[MinumumIntensity,25],...
               'Upper',[1200,200],...
               'Startpoint',[75 70]);


f = fittype('a*exp(-x/b)+0.12*a','options',s);  %offset for bSSFP fitting
%f = fittype('a*exp(-x/b)','options',s);
myomask=ones(size(Ref));

time=[60 31 0]';  % IMPORTANT: double-check
%time=[0 31 60]';  % IMPORTANT: double-check

% make matrices for fitting
[ys xs ~]=size(Ref);

t2map=zeros(ys,xs,nslices);
t2stdevmap=zeros(ys,xs,nslices);
r2map=zeros(ys,xs,nslices);
S0map=zeros(ys,xs,nslices);

warning off
tic %start time keeping
for sn=1:nslices
    coregset=zeros(3,ys,xs);
    coregset(1,:,:)=Coregged1(:,:,sn);
    coregset(2,:,:)=Ref(:,:,sn);
    coregset(3,:,:)=Coregged3(:,:,sn);
    for l=1:ys 
        %replace the 'for' below with 'parfor' to try parallel computing
        for m=1:xs
           if(myomask(l,m)==1&&Coregged3(l,m,sn)>MappingThreshold)     
                [t2fit,parms]=fit(time,squeeze(coregset(:,l,m)),f);
                t2map(l,m,sn)=(t2fit.b); 
                S0map(l,m,sn)=(t2fit.a); 
                %can also saved map of standard deviation and R^2 if desired            
                gg=confint(t2fit);
                t2stdevmap(l,m,sn)=t2fit.b-gg(1,2);
                R2map(l,m,sn)=parms.adjrsquare;
           end
        end
        waitbar(l/ys,h,sprintf('Mapping slice %d...',sn));
     end
end
warning on
fittime=toc; %report time keeping
sprintf('Fitting took %2.0f hours and %2.0f minutes',floor(fittime/3600),rem(fittime,3600)/60) %report time keeping
close(h);
%% plot T2maps depending on ROI
displayall=1;  % toggle whether to display all images or a single
displayimagenumber=1; %only relevant if displayall=0
chooseplotothervars=1; %toggle to plot S0, R2 and the standard deviation
clims=[20 160]; % T2 plot range

plotrows=ceil(nslices/3);

figure(111), set(111,'Color',[1 1 1])
 if nslices>1 && displayall==1
     for kk=1:nslices
         subplot(plotrows,3,kk), imagesc(squeeze(t2map(:,:,kk)),clims),colorbar; 
         axis equal tight off
         title('T2 map')
     end
     if chooseplotothervars==1
         figure(112), set(112,'Color',[1 1 1]), colorbar
         for kk=1:nslices
             clims2=[0 max(max(S0map(:,:,kk)))];
             subplot(plotrows,3,kk),imagesc(squeeze(S0map(:,:,kk)),clims2); 
             axis equal tight off
             title('S0 map')
         end
         figure(113), set(113,'Color',[1 1 1]), colorbar
         for kk=1:nslices
             clims3=[0 150];
             subplot(plotrows,3,kk),imagesc(squeeze(t2stdevmap(:,:,kk)),clims3); 
             axis equal tight off
             title('T2 stdev')
         end
        figure(114), set(114,'Color',[1 1 1]), colorbar
         for kk=1:nslices
             clims4=[0.97 1];
             subplot(plotrows,3,kk),imagesc(squeeze(r2map(:,:,kk)),clims4); 
             axis equal tight off
             title('R^2')
         end
     end
 else
    imagesc(squeeze(t2map(:,:,displayimagenumber)),clims); colorbar
    axis equal tight off
    figure(112), set(112,'Color',[1 1 1])
    imagesc(squeeze(S0map(:,:,1)),clims2)
    axis equal tight off, colorbar
    
end
   

%% select ROI in displayed colormap and calculate mean/stdev
figure(111)
mapmyomask=roipoly;
slicenr=displayimagenumber;
ctr=1;
t2valarr=0;
 for l=1:ys 
    for m=1:xs
        if(mapmyomask(m,l)==1)
           t2valarr(ctr)=t2map(m,l,slicenr);
           ctr=ctr+1;
    end
    end
 end
 myot2=mean(t2valarr)
 myostd=std(t2valarr)
 %% interpolate and display single map if wanted
iptimes=4;
 [dimy,dimx,~]=size(t2map);
[X,Y]=meshgrid(1:dimx,1:dimy);
[XI,YI]=meshgrid(1:1/iptimes:dimx,1:1/iptimes:dimy);
lin=interp2(X,Y,t2map(:,:,displayimagenumber),XI,YI,'cubic'); 
figure(116),set(116,'Color',[1 1 1])
imagesc(lin,clims), axis equal tight off