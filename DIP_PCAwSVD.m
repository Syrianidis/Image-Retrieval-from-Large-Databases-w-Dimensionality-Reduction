% DIP PCA with SVD

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
Isvd=Ids;
for i=1:size(files,1)
    
    Ids{i}=imread(files(i,:));
    Isvd{i}=double(rgb2gray(Ids{i}));
    
end
% and we are going to back to the working directory
cd ..

%% Principal Component Analysis with Singular Value Decomposition

Sds=cell(size(Ids));
for i=1:numel(Sds)
    [~,Sds{i},~]=svd(Isvd{i});
end

%% Dimensionality Reduction and Identifying each Image with referencing

dimensionality=2;
gsv=zeros(1,dimensionality,numel(Sds));
for i=1:size(gsv,3)
    sv=diag(Sds{i});
    gsv(:,:,i)=sv(1:dimensionality);
end

% % % %% Test Efficiency
% % % 
% % % MSE=zeros(1,numel(Ids));
% % % for i=1:numel(Ids)
% % %     Io=Ids{i};
% % % %     figure,subplot(121),imshow(Io),title('I want this one')
% % %     I=rgb2gray(Io);
% % %     I=double(I);
% % %     [~,S,~]=svd(I);
% % %     sv=diag(S);
% % %     id=sum(sv(1:dimensionality));
% % %     id=repmat(id,[1 1 numel(Ids)]);
% % %     index=find(squeeze(abs(gsv-id))==min(squeeze(min(abs(gsv-id),[],3))));
% % % %     subplot(122),imshow(Ids{index}),title('I found this one')
% % %     MSE(i)=norm(double(rgb2gray(Ids{index}))-I);
% % % end
% % % figure,plot(1:numel(MSE),MSE,'b-*')

%% Search

cd([wodi,'\Database']);
files=ls;
files=files(4:size(files,1),:);
from=round(size(files,1)*rand);
I=imread(files(from,:));
cd ..

% % % I=imread('a.jpg');
% % % I=rgb2gray(I);
% % % I=imresize(I,[100,100]);

set(0,'defaulttextinterpreter','none')

figure,subplot(121),imshow(I),title('I want this one'),xlabel(['Image ',files(from,:)])

I=rgb2gray(I);
I=double(I);
[~,S,~]=svd(I);
sv=diag(S);
id=sv(1:dimensionality)';
id=repmat(id,[1 1 numel(Ids)]);

distances=abs(gsv-id);
variance=sum(distances,2);

index=find(variance(:)==min(variance(:)));
subplot(122),imshow(Ids{index}),title('I found this one'),xlabel(['Image ',files(index,:)])

set(0,'defaulttextinterpreter','tex')

mtit('Mining an image from the dataset')

figure
hold on
for i=1:100
    subplot(10,10,i),imagesc(cov(double(rgb2gray(Ids{i})))),axis off
end
hold off