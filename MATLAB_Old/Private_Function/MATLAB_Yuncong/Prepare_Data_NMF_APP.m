% Yuncong Ma, 11/10/2022
% Prepare .nii data in Zheng Zhou's multi site NMF result


%% Zheng Zhou's multi-site result

Path_Raw_Result='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/NMF/Result/Multi_Site/Volume/FN_17_Raw';

Path_APP_Result='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/NMF/Result/Multi_Site/Volume/FN_17';

% Load_Data
fMake_Folder(fullfile(Path_APP_Result,'Load_Data'));
Brain_Mask=fLoad_NII(fullfile(Path_Raw_Result,'atlas','mask_thr0p5_wmparc_2_cc.nii.gz'));
save(fullfile(Path_APP_Result,'Load_Data','Brain_Mask.mat'),'Brain_Mask');
Overlay_Image=single(fLoad_NII(fullfile(Path_Raw_Result,'atlas','old','MNI152_T1_2mm_brain.nii.gz')));
save(fullfile(Path_APP_Result,'Load_Data','Overlay_Image.mat'),'Overlay_Image');

% Compute FN
fMake_Folder(fullfile(Path_APP_Result,'Compute_FN'));

% Group FN
fMake_Folder(fullfile(Path_APP_Result,'Group_FN'));
K=17;
initV=[];
for i=1:K
    Label=['00',num2str(i)];
    initV(:,:,:,i)=fLoad_NII(fullfile(Path_Raw_Result,'fig',['icn_',Label(end-2:end),'.nii.gz']));
end
save(fullfile(Path_APP_Result,'Group_FN','init.mat'),'initV');

% Individualized FN
fMake_Folder(fullfile(Path_APP_Result,'Individualized_FN'));
Dir=dir(fullfile(Path_Raw_Result,'idv_nii','**','fc_4d.nii.gz'));
for i=1:length(Dir)
    [~,temp]=fileparts(Dir(i).folder);
    fMake_Folder(fullfile(Path_APP_Result,'Individualized_FN',temp,'result'));
    clear V
    V{1}=fLoad_NII(fullfile(Path_Raw_Result,'idv_nii',temp,'fc_4d.nii.gz'));
    save(fullfile(Path_APP_Result,'Individualized_FN',temp,'result','final_UV.mat'),'V');
end

%% Visualizaation for gFN and iFN

for i=1:K
    Map=initV(:,:,:,i)*100;
    [Brain_Mask_2,~,Crop_Parameter]=fTruncate_Image_3D_4D(Brain_Mask,[1,1,1],[2,2,2]);
    Overlay_Image_2=fApply_Cropped_FOV(Overlay_Image,Crop_Parameter);
    Map_2=fApply_Cropped_FOV(Map,Crop_Parameter);
    Max_Dim=max(size(Map_2));
    Crop_Parameter.FOV_Old=[1,Max_Dim;1,Max_Dim;1,Max_Dim];
    Crop_Parameter.FOV=[1,size(Map_2,1);1,size(Map_2,2);1,size(Map_2,3)]+round(([1,1,1]*Max_Dim-size(Map_2))/2)';
    Map_2=fInverse_Crop_EPI_Image_3D_4D(Map_2,Crop_Parameter);
    Overlay_Image_2=fInverse_Crop_EPI_Image_3D_4D(Overlay_Image_2,Crop_Parameter);
    threshold=double(round(prctile(Map_2(Brain_Mask_2>0),99.5,"all")));
    Map_Label=bwlabeln(Map_2>=threshold);
    Size=[];
    for j=1:max(Map_Label(:))
        Size(j)=sum(Map_Label(:)==j);
    end
    [~,ps]=max(Map_Label);
    Center=round(fMass_Center(Map_Label==ps(1)));
    Color_Function=fColor_Theme('Seed_Map_3_Positive',round([1/2,1]*threshold));
    Image=fVoxel_Map_3View(Overlay_Image_2,Map_2,{Center,[1,1,1],[2;1;3]},Color_Function,{'Interval',[0,0]});

    fFigure(1,1,1,'',[200,700]);
    set(gcf,'Visible', 'off');
    fSet_Background([200,700],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
    imshow(Image);
    title({['FN ',num2str(i)]},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.6,.15],[0,.3]);
    fColor_Bar(Color_Function,{500,50,80,30,'w','Arial','bold'},round([1/2,1]*threshold),'Horizontal Below Inside');

    Folder=fullfile(Path_APP_Result,'Group_FN');
    saveas(gcf,fullfile(Folder,[num2str(i),'.jpg']));
    close(gcf);

    Image=imread(fullfile(Folder,[num2str(i),'.jpg']));
    Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,1458,1,17],[1,1458,400,417]}});
    imwrite(Image,fullfile(Folder,[num2str(i),'.jpg']));
end

clear Image
for i=1:K
    Image{i}=imread(fullfile(Folder,[num2str(i),'.jpg']));
end
N_Column=10;
N_Row=ceil(K/N_Column);
for i=K+1:(N_Row*N_Column)
    Image{i}=Image{1}*0;
end

imwrite(fAssemble_Image(reshape(Image,[N_Column,N_Row])',{'Interval',[50,5]}),fullfile(Folder,'All.jpg'));



Dir=dir(fullfile(Path_Raw_Result,'idv_nii','**','fc_4d.nii.gz'));
for i0=1:length(Dir)
    [~,temp]=fileparts(Dir(i0).folder);
    fMake_Folder(fullfile(Path_APP_Result,'Individualized_FN',temp,'result'));
    clear V
    V=fLoad_NII(fullfile(Path_Raw_Result,'idv_nii',temp,'fc_4d.nii.gz'));

    for i=1:K
        Map=initV(:,:,:,i)*100;
        [Brain_Mask_2,~,Crop_Parameter]=fTruncate_Image_3D_4D(Brain_Mask,[1,1,1],[2,2,2]);
        Overlay_Image_2=fApply_Cropped_FOV(Overlay_Image,Crop_Parameter);
        Map_2=fApply_Cropped_FOV(Map,Crop_Parameter);
        Max_Dim=max(size(Map_2));
        Crop_Parameter.FOV_Old=[1,Max_Dim;1,Max_Dim;1,Max_Dim];
        Crop_Parameter.FOV=[1,size(Map_2,1);1,size(Map_2,2);1,size(Map_2,3)]+round(([1,1,1]*Max_Dim-size(Map_2))/2)';
        Map_2=fInverse_Crop_EPI_Image_3D_4D(Map_2,Crop_Parameter);
        Overlay_Image_2=fInverse_Crop_EPI_Image_3D_4D(Overlay_Image_2,Crop_Parameter);
        threshold=double(round(prctile(Map_2(Brain_Mask_2>0),99.5,"all")));
        Map_Label=bwlabeln(Map_2>=threshold);
        Size=[];
        for j=1:max(Map_Label(:))
            Size(j)=sum(Map_Label(:)==j);
        end
        [~,ps]=max(Map_Label);
        Center=round(fMass_Center(Map_Label==ps(1)));
        Color_Function=fColor_Theme('Seed_Map_3_Positive',round([1/2,1]*threshold));
        Image=fVoxel_Map_3View(Overlay_Image_2,Map_2,{Center,[1,1,1],[2;1;3]},Color_Function,{'Interval',[0,0]});

        fFigure(1,1,1,'',[200,700]);
        set(gcf,'Visible', 'off');
        fSet_Background([200,700],[0,0,0]);
        fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
        imshow(Image);
        title({['FN ',num2str(i)]},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
        fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.6,.15],[0,.3]);
        fColor_Bar(Color_Function,{500,50,80,30,'w','Arial','bold'},round([1/2,1]*threshold),'Horizontal Below Inside');

        Folder=fullfile(Path_APP_Result,'Individualized_FN',temp);
        saveas(gcf,fullfile(Folder,[num2str(i),'.jpg']));
        close(gcf);

        Image=imread(fullfile(Folder,[num2str(i),'.jpg']));
        Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,1458,1,17],[1,1458,400,417]}});
        imwrite(Image,fullfile(Folder,[num2str(i),'.jpg']));
    end

    clear Image
    for i=1:K
        Image{i}=imread(fullfile(Folder,[num2str(i),'.jpg']));
    end
    N_Column=10;
    N_Row=ceil(K/N_Column);
    for i=K+1:(N_Row*N_Column)
        Image{i}=Image{1}*0;
    end

    imwrite(fAssemble_Image(reshape(Image,[N_Column,N_Row])',{'Interval',[50,5]}),fullfile(Folder,'All.jpg'));
end












