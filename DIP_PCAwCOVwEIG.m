% DIP Principal Component Analysis using EigenDecomposition

%% Loading the images

close all;clear all;clc
% we go to our working directory
cd('C:\Users\Michael\Desktop\MATLAB\DIP');
wodi=pwd;
% one directory down to read from the file 'Database'
cd([wodi,'\Database']);
% we find all the files in 'Database'
files=ls;
% 1st element in FILES is home directory, 2nd is the above directory, 3rd
% is thumb - all of these are unnecessary
files=files(4:size(files,1),:);
% we are storing all the images in a cell array
Ids=cell(1,size(files,1));

for i=1:size(files,1)
    
    Ids{i}=imresize(double(rgb2gray(imread(files(i,:)))),[100,100]);
    
end
% and we are going to back to the working directory
cd ..

%% Principal Component Analysis with Singular Value Decomposition

Ieig=zeros([size(Ids{1},1)*size(Ids{1},2),numel(Ids)]);
for i=1:numel(Ids)
    temp=Ids{i};
    Ieig(:,i)=temp(:);
end

meanimage=mean(Ieig,2)*ones(1,size(Ieig,2));
Ieig=Ieig-meanimage;

Icov=Ieig*Ieig';
Isvd=Ieig'*Ieig;
[Ivec,Ival]=eig(Isvd);
EigenImages=Ieig*Ivec;
EigenImages=cumsum(fliplr(EigenImages),2);
% EigenImages=fliplr(EigenImages);

% I=Ids{1};
% ww=EigenImages'*(I(:)-meanimage(:,1));
% index=find(ww==min(ww));




dimensionality=0;
perc=0;
while perc<.5
    dimensionality=dimensionality+1;
    perc=0;
    
    for i=1:numel(Ids)
        I=Ids{i};
        temp=Ids{i};
        wk=EigenImages(:,1:dimensionality)'*(temp(:)-meanimage(:,1));
        ww=EigenImages(:,1:dimensionality)'*(I(:)-meanimage(:,1));
        index=find((ww-wk)==min(ww-wk));
        
        if index==i
            prec=perc+.01;
        end
        
    end
    
end
        
% % % 
% % % %% Dimensionality Reduction and Identifying each Image with referencing
% % % 
% % % dimensionality=2;
% % % Sds=zeros(size(Isvd{1},1),dimensionality,numel(Isvd));
% % % for i=1:numel(Isvd)
% % %     temp=Isvd{i};
% % %     Sds(:,:,i)=temp(:,1:dimensionality);
% % % end
% % % 
% % % % % % %% Test Efficiency
% % % % % % 
% % % % % % MSE=zeros(1,numel(Ids));
% % % % % % for i=1:numel(Ids)
% % % % % %     Io=Ids{i};
% % % % % % %     figure,subplot(121),imshow(Io),title('I want this one')
% % % % % %     I=rgb2gray(Io);
% % % % % %     I=double(I);
% % % % % %     [~,S,~]=svd(I);
% % % % % %     sv=diag(S);
% % % % % %     id=sum(sv(1:dimensionality));
% % % % % %     id=repmat(id,[1 1 numel(Ids)]);
% % % % % %     index=find(squeeze(abs(gsv-id))==min(squeeze(min(abs(gsv-id),[],3))));
% % % % % % %     subplot(122),imshow(Ids{index}),title('I found this one')
% % % % % %     MSE(i)=norm(double(rgb2gray(Ids{index}))-I);
% % % % % % end
% % % % % % figure,plot(1:numel(MSE),MSE,'b-*')
% % % 
% % % %% Search
% % % 
% % % cd([wodi,'\Database']);
% % % files=ls;
% % % files=files(4:size(files,1),:);
% % % from=round(size(files,1)*rand);
% % % I=imread(files(from,:));
% % % cd ..
% % % 
% % % % % % I=imread('a.jpg');
% % % % % % I=rgb2gray(I);
% % % % % % I=imresize(I,[100,100]);
% % % 
% % % set(0,'defaulttextinterpreter','none')
% % % 
% % % figure,subplot(121),imshow(I),title('I want this one'),xlabel(['Image ',files(from,:)])
% % % 
% % % I=rgb2gray(I);
% % % I=double(I);
% % % 
% % % [q,~]=eig(cov(I));
% % % p=q(:,1:dimensionality);
% % % g=repmat(p,[1 1 100]);
% % % 
% % % 
% % % distances=abs(Sds-g);
% % % variance=sum(sum(distances,2),1);
% % % 
% % % index=find(variance(:)==min(variance(:)));
% % % subplot(122),imshow(Ids{index}),title('I found this one'),xlabel(['Image ',files(index,:)])
% % % 
% % % set(0,'defaulttextinterpreter','tex')
% % % 
% % % mtit('Mining an image from the dataset')