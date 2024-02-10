function Figure = fVoxel_Map(Anatomy,Voxel_Map,Panel_Organization,Color_Function,varargin)
% By Yuncong Ma, Jul. 11, 2021
% Figure=fVoxel_Map(Anatomy,Voxel_Map,Panel_Organization,Color_Function)
% fVoxel_Map is similar to fVoxel_Map, but redesigned
% It supports clustering, add ROI edge map, rotation of Anatomy and
% Voxel_Map, View
% Figure is a [X,Y,3] matrix for imshow
% Panel_Organization could be 0, 2,[3,0],[0,10],[4,5]
% 0 in Color_Function will be regarded as background
% h=fVoxel_Map(_,_,_,_,{"Cluster_Size",Minimum_Cluster_Size})
% Perform 2D cluster size control in default
% Minimum_Cluster_Size is 0 in default
% h=fVoxel_Map(_,_,_,_,{"Cluster_Dimension",3})
% Perform 3D cluster size control
% h=fVoxel_Map(_,_,_,_,{"Rotation",[Rotation_Anatomy,Rotation_Voxel_Map]})
% [Rotation_Anatomy,Rotation_Voxel_Map] is [0,0] in default
% h=fVoxel_Map(_,_,_,_,{"View",'xy'})
% Support 'yz' and 'xz', also support 'XY', 'XZ' and 'YZ'
% h=fVoxel_Map(_,_,_,_,{"Edge",Matrix})
% Matrix is a binarized matrix, in which 1 will be regarded as white
% edges in default
% h=fVoxel_Map(_,_,_,_,{"Edge_Color",RGB_Color})
% RGB_Color is [1 1 1] in default
% Matrix should has the same size as Anatomy
% Rotation of Matrix will be the same as that of Anatomy
% h=fVoxel_Map(_,_,_,_,{"Anatomy_Color_Function",Color_Function})
% Color_Function is [0,0,0,0;1,.8,.8,.8] in default
% h=fVoxel_Map(_,_,_,_,{"Anatomy_Rescale",1})
% Scale the data range of Anatomy for colorization
% The default setting will rescale the data by its 1 and 99 percentiles
% into [0,1]

Options.Cluster_Size=0;
Options.Cluster_Dimension=2;
Options.Edge=[];
Options.Edge_Color=[1 1 1];
Options.Rotation=[0,0];
Options.View='xy';
Options.Anatomy_Color_Function=[0,0,0,0;1,.8,.8,.8];
Options.Anatomy_Rescale=1;
Options=fOption('fVoxel_Map',Options,varargin);
if isempty(Options)
    return;
end

% Anatomy=rot90(Anatomy,Options.Rotation(1));
% Voxel_Map=rot90(Voxel_Map,Options.Rotation(2));
Dimension=size(Voxel_Map);

if ~isequal(size(Anatomy),size(Voxel_Map))
    error('Error in fVoxel_Map: Anatomy and Voxel_Map should have the same size\n');
end

Anatomy=fSlice_Display_3D_Matrix(Anatomy,Options.Rotation(1),Panel_Organization,Options.View);
if Options.Anatomy_Rescale==1
    Anatomy=(Anatomy-prctile(Anatomy(Anatomy>0),1))/(prctile(Anatomy(Anatomy>0),99)-prctile(Anatomy(Anatomy>0),1));
    Anatomy(Anatomy<0)=0;
    Anatomy(Anatomy>1)=1;
    Anatomy=fColorize(Anatomy,Options.Anatomy_Color_Function);
end

Voxel_Map=fColorize(fSlice_Display_3D_Matrix(Voxel_Map,Options.Rotation(2),Panel_Organization,Options.View),Color_Function);

if Options.Cluster_Size>0
    if Options.Cluster_Dimension==3
        Voxel_Map_Mask=reshape(sum(Voxel_Map,3)>0,Dimension);
        mask_map = bwlabeln(Voxel_Map_Mask);
        for i = 1:max(mask_map(:))
            if sum(mask_map(:)==i)<Options.Cluster_Size
                mask_map(mask_map==i)=0;
            end
        end
        mask_map(mask_map>0)=1;
        Voxel_Map_Mask=mask_map;
        Voxel_Map_Mask=reshape(Voxel_Map_Mask,size(Voxel_Map(:,:,1)));
    elseif Options.Cluster_Dimension==2
        Voxel_Map_Mask=sum(Voxel_Map,3)>0;
        mask_map = bwlabel(Voxel_Map_Mask);
        for i = 1:max(mask_map(:))
            if sum(mask_map(:)==i)<Options.Cluster_Size
                mask_map(mask_map==i)=0;
            end
        end
        mask_map(mask_map>0)=1;
        Voxel_Map_Mask=mask_map;
    else
        error('Error in fVoxel_Map: unsupported Cluster_Dimension %f\n',Options.Cluster_Dimension);
    end
else
    Voxel_Map_Mask=sum(Voxel_Map,3)>0;
end

Figure=Anatomy.*(Voxel_Map_Mask==0)+Voxel_Map.*Voxel_Map_Mask;

if ~isempty(Options.Edge)
    Edge=rot90(Options.Edge,Options.Rotation(1));
    if ~isequal(Dimension,size(Options.Edge))
        error('Error in fVoxel_Map: the dimension of Atlas and Edge should be the same\n');
    end
    Edge=fSlice_Display_3D_Matrix(Edge,0,Panel_Organization,Options.View);
    Edge=fColorize(Edge,[0,0,0,0;1,Options.Edge_Color]);
    Figure=Figure.*repmat(sum(Edge,3)==0,[1 1 3])+Edge;
end

end