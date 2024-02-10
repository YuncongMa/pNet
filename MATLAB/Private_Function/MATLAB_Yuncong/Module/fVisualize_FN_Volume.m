function fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure,varargin)
% Yuncong Ma, 2/1/2024
% Visualize FN in volume format for NMF APP
% fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure)
% Work for both group-level and individualized FN

Options.Color_Range_Style='SR-NMF';
Options.Output_Setting=0; % 1 for gFNs to sychonize center view
Options.Setting_Folder=''; % pFNs to sychonize center view
Options.File_Center=[];
Options=fOption('fVisualize_FN_Volume',Options,varargin);
if isempty(Options)
    return;
end
Color_Range_Style=Options.Color_Range_Style;

load(File_FN,'FN');

[Brain_Template,Flag,Message]=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Data_Input','Brain_Template.mat'));
if Flag
    return
end
Brain_Mask=Brain_Template.Brain_Mask;
Overlay_Image=Brain_Template.Overlay_Image;

if isequal(size(Brain_Mask),size(Overlay_Image))
    Upsampling=1;
else
    Upsampling=round(size(Overlay_Image,1)/size(Brain_Mask,1));
    if ~isequal(size(imresize3(Brain_Mask,Upsampling)),size(Overlay_Image))
        error('Error in fVisualize_FN_Volume: the overlay image does not have an integer times of higher resolution compared to the brain mask');
    end
end

K=size(FN,4);

% setting for colorbar
Color_Bar.Options.Height=50;
Color_Bar.Options.Width=500;
Color_Bar.Options.Ticks=[];
Color_Bar.Options.Tick_Font_Size=20;
Color_Bar.Options.Tick_Font_Color='w';
Color_Bar.Options.Tick_Font_Name='Arial';
Color_Bar.Options.Tick_Font_Weight='bold';
Color_Bar.Options.Tick_Interval=40;
Color_Bar.Options.Style='Horizontal Below Inside';
Color_Bar.Options.Label='';
Color_Bar.Options.Label_Font_Size=20;
Color_Bar.Options.Label_Font_Color='w';
Color_Bar.Options.Label_Font_Name='Arial';
Color_Bar.Options.Label_Font_Weight='bold';
Color_Bar.Options.Label_Interval=120;
switch Color_Range_Style
    case 'SR-NMF'
        Color_Bar.Options.Label='Loading (%)';
    case 'GIG-ICA'
        Color_Bar.Options.Label='Z';
end

% rescale
switch Options.Color_Range_Style
    case 'SR-NMF'
        FN=FN*100;
    case 'GIG-ICA'
        % no change
end

for i=1:K
    
    % crop map and overlay image
    Map=imresize3(FN(:,:,:,i),Upsampling,'nearest');
    [Brain_Mask_2,~,Crop_Parameter]=fTruncate_Image_3D_4D(imresize3(Brain_Mask,Upsampling,'nearest'),[1,1,1],[2,2,2]);
    Overlay_Image_2=fApply_Cropped_FOV(Overlay_Image,Crop_Parameter);
    Map_2=fApply_Cropped_FOV(Map,Crop_Parameter);
    Max_Dim=max(size(Map_2));
    Crop_Parameter.FOV_Old=[1,Max_Dim;1,Max_Dim;1,Max_Dim];
    Crop_Parameter.FOV=[1,size(Map_2,1);1,size(Map_2,2);1,size(Map_2,3)]+round(([1,1,1]*Max_Dim-size(Map_2))/2)';
    Map_2=fInverse_Crop_EPI_Image_3D_4D(Map_2,Crop_Parameter);
    Overlay_Image_2=fInverse_Crop_EPI_Image_3D_4D(Overlay_Image_2,Crop_Parameter);

    % set color range
    switch Color_Range_Style
        case 'SR-NMF'
            threshold=double(round(prctile(Map_2(Brain_Mask_2>0),99.5,"all")));
            Color_Range=round([1/2,1]*threshold);
            Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
        case 'GIG-ICA'
            threshold=double(prctile(Map_2(Brain_Mask_2>0),99.8,"all"));
            if threshold>1
                Color_Range=round([0.5,1]*threshold*10)/10;
            elseif threshold>0.01
                Color_Range=round([0.5,1]*threshold*100)/100;
            else
                Color_Range=[0.5,1]*threshold;
            end
            Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
    end

    % find center
    if ~isempty(Options.Setting_Folder)
        File_Setting=fullfile(Options.Setting_Folder,['FN_',num2str(i),'.json']);
        Setting=load_json_setting(File_Setting);
        Center=Setting.Center;
    else
        Center=large_3view_center(Map_2);
    end

    % get main figure data
    Image=fVoxel_Map_3View(Overlay_Image_2,Map_2,{Center,[1,1,1],[2;1;3]},Color_Function,{'Interval',[0,0]});

    % organize figure
    fFigure(1,1,1,'',[200,700]);
    set(gcf,'Visible', 'off');
    fSet_Background([200,700],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
    imshow(Image);
    title({['FN ',num2str(i)]},'Fontsize',30,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.8,.2],[0,.4]);
    fColor_Bar(Color_Function,Color_Bar.Options);

    % output figure
    saveas(gcf,fullfile(Dir_Figure,[num2str(i),'.jpg']));
    close(gcf);

    % remove white lines due to issues in MATLAB
    Image=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
    Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,1458,1,17],[1,1458,400,417]}});
    imwrite(Image,fullfile(Dir_Figure,[num2str(i),'.jpg']));

    % output settings
    if Options.Output_Setting
        fMake_Folder(fullfile(Dir_Figure,'Figure_Setting'));
        Setting.Color_Function=Color_Function;
        Setting.Color_Range=Color_Range;
        Setting.Center=Center;
        File_Setting=fullfile(Dir_Figure,'Figure_Setting',['FN_',num2str(i),'.json']);
        write_json_setting(Setting, File_Setting);
    end
end

% assemble images
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
