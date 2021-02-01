function [AllData, TT, imStat]=MainGerry(tiffile,f)

% This function is used to identify single cells that will be fed into the 
% analysis algorithm GerryDist() to organize the data into exportable
% tables
% INPUT:
%
%           -tiffile: the name of the .tif file. The tif file is a stack of
%                     grayscale z-sections of a single channel (GFP)

%           -f: multiplying factor of the background value of the image to
%               determine the threshold used for binarization fo the volume

% OUTPUT:
%
%           -AllData: n-by-2 cell array containing a data table for each yeast cell
%                     in the first column and the size of the bounding box
%                     of each cell volume in the second column
%
%           -TT: The final exportable table of p-by-6 size. Each row
%                corresponds to a single "element" (object), and each column
%                corresponds to a measured variable.
%
%           -imStat: a 2-by-1 vector containing useful statistics of the
%           whole image. imStat(1) corresponds to the background value of
%           the image, imStat(2) corresponds to the 50th brightest pixel in
%           the image.
%
% Example:
%
% [AllData, TT, imStat] = MainGerry('9506_01.tif',3);

%% Load and pre-process the stack

IS=bfopen(tiffile);
ISC=cat(3,IS{1,1}{:,1});
nzstacks=size(ISC,3);

ISP=max(ISC,[],3); %max projection

T=minminT(ISP); %find background value
SS=sort(ISP(:),'descend');
H=SS(50); %find 50th brightest pixel
imStat=[T H];
ISPT=ISP-f*T;
ISCT=ISC-f*T;
figure, imshow(ISPT,[0 1000]);

%% Single yeast cell detection and data initialization

IW=imbinarize(ISPT);
IW=bwareaopen(IW,5);
IWD=imdilate(IW,strel('disk',10));
IWD=imclearborder(IWD);
Stats=regionprops(IWD,'Area','BoundingBox');
A=[Stats.Area]';
tf=A<2000 | A>20000;
Stats(tf)=[];

for i=1:numel(Stats)
    rectangle(gca,'Position',Stats(i).BoundingBox,'EdgeColor','r','LineWidth',1)
    text(Stats(i).BoundingBox(1)+10,Stats(i).BoundingBox(2)+10,int2str(i),'Color','r','FontSize',12)
    %idx=idx+1;
end

%ncell=size(bb,1);
ncell=numel(Stats);
AllData=cell(ncell,2);

%% Single yeast analysis

for i=1:ncell

    BB=Stats(i).BoundingBox;
    %BB=bb(i,:);
    CROP=imcrop(ISP,BB);
    figure, imshow(CROP,[],'InitialMagnification',300);
    [a, b]=size(CROP);
    singleplombe=zeros(a,b,nzstacks);
    for j=1:nzstacks
        singleplombe(:,:,j)=imcrop(ISCT(:,:,j),BB);
    end
    
    [spStats, spsize]=GerryDist(singleplombe); % analysis function
    AllData{i,1}=spStats;
    AllData{i,2}=spsize;
end




%% Create final table

TT=[];
sname=[tiffile(1:end-4) '_'];
for i=1:size(AllData,1)
T=AllData{i,1};
ne=size(T,1);
cid=cellstr([repmat(sname,ne,1) int2str(i*ones(ne,1))]);
CID=table(cid,'VariableNames',{'CellID'});
T=[CID T];
TT=[TT;T];
%idx=idx+1;
end
 %idd=idx;
% fname=[tiffile(1:end-4) '_P.tif'];
% export_fig(fname,'-tif');
% close;
end

function thresvalue=minminT(I)

thresvalue = max([min(max(I,[],1)) min(max(I,[],2))])	 ;

end

    
