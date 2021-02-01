function [spStats, spSize]=GerryDist(singleplombe)

% This function is called by the helper function MainGerry() and it
% extracts metrics used to analyze GFP distribution
%
% INPUT:
%
%           -singleplombe: cropped stack from the original tif image
%           containing a single yeast cell
%
% OUTPUT:
%
%           -spStats: n-by-5 table containing the extracted metrics for an
%           individual yeast cell with n "elements" (objects)
%           -spSize: 3-by-1 vector containing the size (in pixels) of the
%           bounding box.

%% Initialize data

spSize=size(singleplombe);

%% Segment GFP signal and extract metrics

SW=imbinarize(singleplombe);
%SW=singleplombe>0;

if any(SW(:))
    StatsW=regionprops3(SW,singleplombe,'Volume', 'VoxelIdxList', 'VoxelValues', 'MeanIntensity');
    
    
    figure,View3D(SW,'r');
    
    
    %% Get maximum length using geodesic distance between endpoints of skeletonized volume
    SWK=bwskel(SW);
    SWKm=max(SWK,[],3);
    figure,View3D(SWK,'b');
    figure, imshow(SWKm,'InitialMagnification',300);
    spStats=regionprops3(SWK,'VoxelIdxList');
    maxLength=zeros(size(spStats,1),1);
    Volume=zeros(size(spStats,1),1);
    MI=zeros(size(spStats,1),1);
    IM=zeros(size(spStats,1),1);
    Vox=cell(size(spStats,1),1);
    
    for i=1:size(spStats,1)
        MASK=false(spSize);
        if iscell(spStats(i,1).VoxelIdxList)
            Pix=spStats.VoxelIdxList{i,1};
        else
            Pix=spStats(i,1).VoxelIdxList;
        end
        MASK(Pix)=true;
        epidx=find(bwmorph3(MASK,'endpoints'));
        maxlength=zeros(numel(epidx),1);
        for j=1:numel(epidx)
            dd=bwdistgeodesic(MASK,epidx(j));
            maxlength(j)=max(dd(:));
        end
        if ~isempty(maxlength)
            maxLength(i)=max(maxlength);
        else
            maxLength(i)=1;
        end
        D=arrayfun(@(x) any(intersect(Pix,x.VoxelIdxList)),table2struct(StatsW));
        k=find(D);
        Volume(i)=StatsW.Volume(k);
        MI(i)=StatsW.MeanIntensity(k);
        IM(i)=sum(StatsW.VoxelValues{k,1});
        Vox{i}=StatsW.VoxelIdxList{k,1};
    end
    
    spStats.maxLength=maxLength;
    spStats.Volume=Volume;
    spStats.MeanIntensity=MI;
    spStats.IntegralIntensity=IM;
    spStats.VoxelIdxList=[];
    spStats.VoxelIdxList=Vox;
    
else
    spStats.maxLength=0;
    spStats.Volume=0;
    spStats.MeanIntensity=0;
    spStats.IntegralIntensity=0;
    spStats.VoxelIdxList=[];
    spStats=struct2table(spStats);
    
end

end

function View3D(ZBW,Col,varargin)


[d,e,f]=size(ZBW);                   
[x, y, z]=meshgrid(1:e,1:d,1:f);       
p=patch(isosurface(x,y,z,ZBW,0.5));
p.FaceColor=Col;
p.EdgeColor='none';
camlight;
lighting gouraud
axis off;
set(gca,'ydir','reverse','zdir','reverse')
set(gcf,'Color','w')


if isempty(varargin)
    d=[2 2 1];
else
    scale=varargin{1};
    ss=scale(2)/scale(1);
    d=[ss ss 1];
end
daspect(d);

end
