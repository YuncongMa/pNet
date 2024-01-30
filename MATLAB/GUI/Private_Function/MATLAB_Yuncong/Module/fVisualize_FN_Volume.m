function fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure,varargin)
% Yuncong Ma, 1/29/2024
% Visualize FN in volume format for NMF APP
% fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure)
% Work for both group-level and individualized FN

Options.Color_Range_Style='NMF_1';
Options=fOption('fVisualize_FN_Volume',Options,varargin);
if isempty(Options)
    return;
end

load(File_FN,'FN');
Brain_Mask=fLoad_MATLAB_Single_Variable(File_Brain_Mask);
Overlay_Image=fLoad_MATLAB_Single_Variable(File_Overlay_Image);

if isequal(size(Brain_Mask),size(Overlay_Image))
    Upsampling=1;
else
    Upsampling=round(size(Overlay_Image,1)/size(Brain_Mask,1));
    if ~isequal(size(imresize3(Brain_Mask,Upsampling)),size(Overlay_Image))
        error('Error in fVisualize_FN_Volume: the overlay image does not have an integer times of higher resolution compared to the brain mask');
    end
end

K=size(FN,4);

clear Image

% rescale
switch Options.Color_Range_Style
    case 'NMF_1'
        FN=FN*100;
end

for i=1:K
    
    Map=imresize3(FN(:,:,:,i),Upsampling,'nearest');
    [Brain_Mask_2,~,Crop_Parameter]=fTruncate_Image_3D_4D(imresize3(Brain_Mask,Upsampling,'nearest'),[1,1,1],[2,2,2]);
    Overlay_Image_2=fApply_Cropped_FOV(Overlay_Image,Crop_Parameter);
    Map_2=fApply_Cropped_FOV(Map,Crop_Parameter);
    Max_Dim=max(size(Map_2));
    Crop_Parameter.FOV_Old=[1,Max_Dim;1,Max_Dim;1,Max_Dim];
    Crop_Parameter.FOV=[1,size(Map_2,1);1,size(Map_2,2);1,size(Map_2,3)]+round(([1,1,1]*Max_Dim-size(Map_2))/2)';
    Map_2=fInverse_Crop_EPI_Image_3D_4D(Map_2,Crop_Parameter);
    Overlay_Image_2=fInverse_Crop_EPI_Image_3D_4D(Overlay_Image_2,Crop_Parameter);

   
    threshold=double(round(prctile(Map_2(Brain_Mask_2>0),99.5,"all")));
    Map_Label=bwlabeln(Map_2>threshold*0.9);
    Size=[];
    for j=1:max(Map_Label(:))
        Size(j)=sum(Map_Label(:)==i);
    end
    [~,ps]=max(Map_Label);
    Center=round(fMass_Center(Map_Label==ps(1)));
    switch Options.Color_Range_Style
        case 'NMF_1'
            Color_Range=round([1/2,1]*threshold);
            Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
        case 'ICA'
            if threshold>1
                Color_Range=round([0.5,1]*threshold*10)/10;
                Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
            else
                Color_Range=[threshold/2,threshold];
                Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
            end
    end
    Image=fVoxel_Map_3View(Overlay_Image_2,Map_2,{Center,[1,1,1],[2;1;3]},Color_Function,{'Interval',[0,0]});

    fFigure(1,1,1,'',[200,700]);
    set(gcf,'Visible', 'off');
    fSet_Background([200,700],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
    imshow(Image);
    title({['FN ',num2str(i)]},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.6,.15],[0,.3]);
    fColor_Bar(Color_Function,{500,50,80,30,'w','Arial','bold'},Color_Range,'Horizontal Below Inside');

    saveas(gcf,fullfile(Dir_Figure,[num2str(i),'.jpg']));
    close(gcf);

    Image=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
    Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,1458,1,17],[1,1458,400,417]}});
    imwrite(Image,fullfile(Dir_Figure,[num2str(i),'.jpg']));
end

clear Image
for i=1:K
    Image{i}=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
end
N_Column=10;
N_Row=ceil(K/N_Column);
for i=K+1:(N_Row*N_Column)
    Image{i}=Image{1}*0;
end

imwrite(fAssemble_Image(reshape(Image,[N_Column,N_Row])',{'Interval',[50,5]}),fullfile(Dir_Figure,'All.jpg'));
