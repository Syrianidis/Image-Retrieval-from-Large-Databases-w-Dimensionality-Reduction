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
    
end
% and we are going to back to the working directory
cd ..

%% Adding White Gaussian Noise to our Dataset

SNR=0;
for i=1:numel(Ids)
    
    Isvd{i}=double(rgb2gray(Ids{i}));
    Isvd{i}=255*addWGN(Isvd{i},0,SNR);
    
end

%% Principal Component Analysis with Singular Value Decomposition

Sds=cell(size(Ids));
for i=1:numel(Sds)
    [~,Sds{i},~]=svd(Isvd{i});
end

%% Dimensionality Reduction and Identifying each Image with referencing


dimensionality=0;
perc=0;
while perc<.95

    dimensionality=dimensionality+1;
    perc=0;

    gsv=zeros(1,dimensionality,numel(Sds));
    for i=1:size(gsv,3)
        sv=diag(Sds{i});
        gsv(:,:,i)=sv(1:dimensionality);
    end

%% Search

    for i=1:numel(Ids)

        cd([wodi,'\Database']);
        files=ls;
        files=files(4:size(files,1),:);
        from=i;
        I=imread(files(from,:));
        cd ..

        I=rgb2gray(I);
        I=double(I);
        [~,S,~]=svd(I);
        sv=diag(S);
        id=sv(1:dimensionality)';
        id=repmat(id,[1 1 numel(Ids)]);
        distances=abs(gsv-id);
        variance=sum(distances,2);

        index=find(variance(:)==min(variance(:)));
%         Iq=double(rgb2gray(Ids{index}));
%         if norm(Iq-I)==0
        if index==i
            perc=perc+.01;
        end

    end
    
    if dimensionality==100
        break
    end
end

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

figure,subplot(121),imshow(I),title('I want this one'),xlabel(['Image ',files(from,:)]),ylabel(['Dimensionality ',num2str(dimensionality)])

I=rgb2gray(I);
I=double(I);
[~,S,~]=svd(I);
sv=diag(S);
id=sv(1:dimensionality)';
id=repmat(id,[1 1 numel(Ids)]);
distances=abs(gsv-id);
variance=sum(distances,2);

index=zeros(1,5);
for i=1:numel(index)
    index(i)=find(variance(:)==min(variance(:)));
    variance(index(i))=Inf;
end
for i=1:5
    subplot(5,2,2*i),imshow(Ids{index(i)}),title(['Image ',files(index(i),:)])
end
mtit(['SNR ',num2str(SNR)])

set(0,'defaulttextinterpreter','tex')