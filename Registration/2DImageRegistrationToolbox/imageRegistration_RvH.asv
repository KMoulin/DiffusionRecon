% IMAGE REGISTRATION MAIN
% -------------------------------------------------------------------------
% load images
close all
clear all
clc

%read in files, with the series directory numbers in the proper order
% bigimage will have the dimensions: [series,dimx,dimy,slice]
series=['13' '09' '05'];
for k=1:length(series)/2
    filename='MR000002';
    path=strcat('C:\Users\Ruud\Documents\Work CHUV\Trio\20110120 T2map Eline\DICOM\ST000000\SE0000',series(2*k-1:2*k),'\');
    [image,infoimage]=opendicomsequence(filename,path);
    bigimage(k,:,:,:)=image;
    biginfoimage(k)=infoimage;
end

%% interactive crop

% REFERENCE IMG -------------
bigRef = squeeze(bigimage(2,:,:)); 
% moving images
bigI2 = squeeze(bigimage(1,:,:));
bigI3 = squeeze(bigimage(3,:,:));

figure(1)
imagesc(bigRef), axis equal tight, colormap('gray')
bb=imrect;
pos=round(getPosition(bb));
close 1
Ref=bigRef(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
I2=bigI2(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
I3=bigI3(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
%imagesc(I3), axis equal tight off, colormap('gray')
%%
% figure('name','original images')
% subplot(121),imshow(Ref,[]),title('reference image')
% subplot(122),imshow(I2,[]),title('moving image')

Options = struct('Similarity',{'cc'},...
                      'Registration',{'NonRigid'},...
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
% 'mi' : normalized (local) mutual information
% 'mip': normalized mutual information with image split in multiple small regions
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
% Registration Options:
%               Rigid    : Translation, Rotation
%               Affine   : Translation, Rotation, Shear, Resize
%               NonRigid : B-spline grid based local registration
%               Both     : Nonrigid and Affine (Default)


[Coregged1,O_trans1,Spacing1,M1,param1,B1,F1] = image_registration (I2,Ref,Options);
[Coregged2,O_trans2,Spacing2,M2,param2,B2,F2] = image_registration (I3,Ref,Options);
% M: affine matrix of the registration
% param: parameters of the registration
%        2D - [ deltaY deltaX theta ]
% -------------------------------------------------------------------------

% extra coregistration parameters:
%error = abs(Coregged1-I2);
% coregYdisp = param(1);
% coregXdisp = param(2);
% coregRot = param(3);
%% bypass coregistration
Coregged1=I2;
Coregged2=I3;
%% -------------plot-------------------------------------------
figure(1)
clims=[0 600];
subplot(131),imagesc(Coregged1,clims),title('T2prep=0'), axis equal tight off, colormap('gray') 
subplot(132),imagesc(Ref,clims),title('T2prep=31'), axis equal tight off, colormap('gray')
subplot(133),imagesc(Coregged2,clims),title('T2prep=60'), axis equal tight off, colormap('gray')

figure(2)
clims=[0 400];
subplot(131),imagesc(I2,clims),title('T2prep=0'), axis equal tight off, colormap('gray') 
subplot(132),imagesc(Ref,clims),title('T2prep=31'), axis equal tight off, colormap('gray')
subplot(133),imagesc(I3,clims),title('T2prep=60'), axis equal tight off, colormap('gray')

figure(3)
subplot(121),imagesc(abs(I3-Ref),clims/4),title('before'), axis equal tight off, colormap('gray') 
subplot(122),imagesc(abs(Coregged2-Ref),clims/4),title('after'), axis equal tight off, colormap('gray')

%% select 2 ROIs: outer myocardium first, then inner
figure(4)
imagesc(Ref), axis equal tight, colormap('gray')
outer=roipoly;
imagesc(Ref), axis equal tight, colormap('gray')
inner=roipoly;
myomask=outer-inner;
close(4)
%% display fit in single pixel
figure(9)
imagesc(Ref), axis equal tight off, colormap('gray')
[xp,yp]=ginput(1);
close(9)

[xs ys]=size(Ref);
coregset=zeros(3,xs,ys);
coregset(1,:,:)=Coregged1;
coregset(2,:,:)=Ref;
coregset(3,:,:)=Coregged2;


s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[50,5],...
               'Upper',[1000,500],...
               'Startpoint',[200 50]);
f = fittype('a*exp(-x/b)','options',s);
time=[0 31 60]';
sig=squeeze(coregset(:,round(xp),round(yp)));
figure(2)
plot(time,sig,'o')
xlim([0 1.1*max(time)]),ylim([0 1.1*max(sig)])
[pixt2fit,pixparms]=fit(time,sig,f)
hold on
plot(pixt2fit,'r')
xlabel('time [ms]')
ylabel('image intensity [a.u.]')
%% calculate the T2map for a slice
%set fitting parameters
h = waitbar(0,'Waiting for T2mapping ...');
time=[0 31 60]';
fitvars=2;
if fitvars==3
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[20,10,1],...
                   'Upper',[1000,400,500],...
                   'Startpoint',[200 50 50]);
    f = fittype('a*exp(-x/b)+c','options',s);
elseif fitvars==2
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[30,10],...
                   'Upper',[700,400],...
                   'Startpoint',[150 50]);
    f = fittype('a*exp(-x/b)','options',s);
end

% make matrices for fitting
[xs ys]=size(Ref);
coregset=zeros(3,xs,ys);
coregset(1,:,:)=Coregged1;
coregset(2,:,:)=Ref;
coregset(3,:,:)=Coregged2;

t2map=zeros(xs,ys);
t2stdevmap=zeros(xs,ys);
r2map=zeros(xs,ys);


tic %start time keeping

 for l=1:xs 
    for m=1:ys
       if(myomask(l,m)==1)     
            [t2fit,parms]=fit(time,squeeze(coregset(:,l,m)),f);
            t2map(l,m)=t2fit.b;
            if fitvars<length(time)
                gg=confint(t2fit);
                t2stdevmap(l,m)=t2fit.b-gg(1,2);
                r2map(l,m)=parms.adjrsquare;
            end
       end
    end
    waitbar(l/ys);
 end
 fittime=toc %report time keeping
close(h);
%% plot overlaid image
leftshift=0; % can be used in case of a slight image mismatch (due to different animal positioning)
upshift=0;
figure(10)
clims=[0 100];

bgrgb=repmat((Ref-min(min(Ref)))/(max(max(Ref))-min(min(Ref))),[1 1 3]);
imbg=imshow(bgrgb,[0 100],'InitialMagnification','fit');
hold on
colormap('jet')
imoverlay=imagesc(t2map,'XData',[1-leftshift ys-leftshift],'YData',[1-upshift xs-upshift],clims);
set(imoverlay,'AlphaData',myomask);

axis equal tight off

%% calculate mean and stdev in map
ctr=1;
 for l=1:xs 
    for m=1:ys
        if(myomask(l,m)==1)
           t2valarr(ctr)=t2map(l,m);
           ctr=ctr+1;
        end
    end
 end
 myot2=mean(t2valarr)
 myostd=std(t2valarr)
 
 