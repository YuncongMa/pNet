function h=fVoxel_Map_3View(Anatomy,Voxel_Map,View_Setting,Color_Function,varargin)
% By Yuncong Ma, 11/14/2022
% h=fVoxel_Map_3View(Anatomy,Voxel_Map,View_Setting,Color_Function)
% View_Setting is {[Center_X,Center_Y,Center_Z],[Rotation_XY,Rotation_XZ,Rotation_YZ],Organization}
% Center_X should be integer, Rotation_XY should be integer
% Organization is [1,3;2,0], [2,3;1,0] for displaying XY, YZ in the first row, while XZ
% and black background in the second row
% Organization is [2;1;3] for displaying XZ,YZ,XY in one column
% Color_Function could be set by fColor_Theme or a matrix which each row
% consists of [Value,R,G,B]
% h=fVoxel_Map_3View(_,_,_,_,Options),
% {'Background','w'} the default background color is 'k'
% {'Interval',[Interval_X,Interval_Y]} to separate each subfigure


Options.Background='k';
Options.Interval=[2,2];
Options=fOption('fVoxel_Map_3View',Options,varargin);
if isempty(Options)
    return;
end
if strcmp(Options.Background,'k')
    BackgroundColor=0;
else
    BackgroundColor=1;
end

if ~isequal(size(Anatomy),size(Voxel_Map))
    fprintf('Error in fVoxel_Map_3View: Anatomy and Voxel_Map should have the same size\n');
    return;
end

Dimension=size(Anatomy);

View_Center=View_Setting{1};
Rotation=View_Setting{2};
Organization=View_Setting{3};

Anatomy2D={rot90(Anatomy(:,:,View_Center(3)),Rotation(1)),...
    rot90(squeeze(Anatomy(:,View_Center(2),:)),Rotation(2)),...
    rot90(squeeze(Anatomy(View_Center(1),:,:)),Rotation(3))};

for i=1:3
    range=prctile(Anatomy2D{i}(Anatomy2D{i}>0),[1,99]);
    Anatomy2D{i}=(Anatomy2D{i}-range(1))/diff(range)*0.8;
    Anatomy2D{i}=repmat(Anatomy2D{i},[1,1,3]);
end

Voxel_Map2D={rot90(Voxel_Map(:,:,View_Center(3)),Rotation(1)),...
    rot90(squeeze(Voxel_Map(:,View_Center(2),:)),Rotation(2)),...
    rot90(squeeze(Voxel_Map(View_Center(1),:,:)),Rotation(3))};

Mask=cell(1,3);
for i=1:3
    Voxel_Map2D{i}=fColorize(Voxel_Map2D{i},Color_Function);
    Mask{i}=repmat(sum(Voxel_Map2D{i},3)==0,[1,1,3]);
end

if isequal(Organization,[1,3;2,0]) || isequal(Organization,[2,3;1,0])
    Figure=cell(2,2);
    for x=1:2
        for y=1:2
            i=Organization(x,y);
            if i==0
                Figure{x,y}=BackgroundColor+zeros(Dimension(3),Dimension(3),3);
            else
                Figure{x,y}=Anatomy2D{i}.*Mask{i}+Voxel_Map2D{i}.*(1-Mask{i});
            end
        end
    end
    h=fAssemble_Image(Figure,{'Background',Options.Background},{'Interval',Options.Interval});

elseif isequal(size(Organization),[3,1])
    Figure=cell(3,1);
    for i=1:3
        Figure{i,1}=Anatomy2D{Organization(i)}.*Mask{Organization(i)}+Voxel_Map2D{Organization(i)}.*(1-Mask{Organization(i)});
    end
    h=fAssemble_Image(Figure,{'Background',Options.Background},{'Interval',Options.Interval});

else
    fprintf('Error in fVoxel_Map_3View: unsupported Organization settings\n');
    return;
end

end



