function fVisualize_Statistics(Folder,varargin)
% Yuncong Ma, 5/23/2023
% Visualize statistical tests for FN
% fVisualize_Statistics(Folder)
% Options support 'Surface' and 'HCP_Surface'

Options.Data_Type='Surface';
Options.Data_Format='HCP_Surface';
Options.Folder_Result='';
Options=fOption('fVisualize_Statistics',Options,varargin);
if isempty(Options)
    return;
end

Path_LoadData=fullfile(Options.Folder_Result,'Load_Data');

% T Value for t-test
% Z Value for Wilcoxon tests
load(fullfile(Folder,'Statistics.mat'),'Statistics');
Result=Statistics.Result;
if isfield(Result,'T_Value') && ~isempty(Result.T_Value)
    Value=Result.T_Value;
elseif isfield(Result,'Z_Value') && ~isempty(Result.Z_Value)
    Value=Result.Z_Value;
else
    return;
end

switch Options.Data_Type
    case 'Surface'
        Path_Brain_Surface=fullfile(Path_LoadData,'Brain_Surface.mat');
        if ~exist(Path_Brain_Surface,'file')
            msgbox('Cannot find the brain surface file','Warning',"Warn");
            return;
        end
        Brain_Surface=Load_Mathematica_Mat(Path_Brain_Surface);
        
        Result.P_Value(isnan(Result.P_Value))=1;
        Value(isnan(Value))=0;
        Value=Value.*(Result.P_Value<=Statistics.PValue);
        K=size(Result.P_Value,2);
        clear Image

        fig=waitbar(0,'Visualizing ...');
        for i=1:K
            if ~isvalid(fig)
                return
            end
            waitbar((i-1)/K,fig,'Visualizing ...');

            Max=double(round(10*prctile(abs(Value(:,i)),99))/10);
            if Max==0
                Max=double(round(10*max(abs(Value(:,i))))/10);
                if Max==0
                    Max=1;
                end
            end
            Color_Function=fColor_Theme('Seed_Map_3',[0.00001,Max]);
            Image=fHCP_FS_6View_Column(Brain_Surface.Shape_Inflated,fHCP_To_32K(Value(:,i),Brain_Surface.MW.LR),Color_Function,...
                {'Material','dull'});

            fFigure(1,1,1,'',[200,1100]);
            set(gcf,'Visible', 'off');
            fSet_Background([200,1100],[0,0,0]);
            fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],1,0.15,[.7,.1],[0,.4]);
            imshow(Image);
            title({['FN ',num2str(i)],''},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
            fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.15,[.7,.1],[0,.4]);
            fColor_Bar(Color_Function,{500,50,100,30,'w','Arial','bold'},[-1,1]*Max,'Horizontal Below Inside');

            saveas(gcf,fullfile(Folder,[num2str(i),'.jpg']));
            close(gcf);

            Image=imread(fullfile(Folder,[num2str(i),'.jpg']));
            Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,2292,1,25],[1,2292,380,417]}});
            imwrite(Image,fullfile(Folder,[num2str(i),'.jpg']));

            fCrop_Image_File(fullfile(Folder,[num2str(i),'.jpg']),[0,0.08,1,.93]);
            
        end
        clear Center Index
        if isvalid(fig)
            waitbar(1,fig,'Visualizing ...');
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
        if isvalid(fig)
            close(fig);
        end


    case 'Volume'

        Path_Brain_Mask=fullfile(Path_LoadData,'Brain_Mask.mat');
        if ~exist(Path_Brain_Mask,'file')
            msgbox('Cannot find the brain mask file','Warning',"Warn");
            return;
        end
        Brain_Mask=Load_Mathematica_Mat(Path_Brain_Mask);
        Path_Overlay_Image=fullfile(Path_LoadData,'Overlay_Image.mat');
        if ~exist(Path_Overlay_Image,'file')
            msgbox('Cannot find the overlay image file','Warning',"Warn");
            return;
        end
        Overlay_Image=Load_Mathematica_Mat(Path_Overlay_Image);
        Upsampling=round(size(Overlay_Image,1)/size(Brain_Mask,1));
        load(fullfile(Folder,'Statistics.mat'),'Statistics');
        Result=Statistics.Result;
        Value(isnan(Value))=0;
        Result.P_Value(isnan(Result.P_Value))=1;
        Value=Value.*(Result.P_Value<=Statistics.PValue);
        K=size(Result.P_Value,4);
        clear Image

        fig=waitbar(0,'Visualizing ...');
        for i=1:K
            if ~isvalid(fig)
                return
            end
            waitbar((i-1)/K,fig,'Visualizing ...');

            Map=imresize3(Value(:,:,:,i),Upsampling,'nearest');
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
            Max=double(round(10*prctile(abs(Value(:,:,:,i)),99,'all'))/10);
            if Max==0
                Max=double(round(10*max(abs(Value(:,:,:,i)),[],'all'))/10);
                if Max==0
                    Max=1;
                end
            end
            Color_Function=fColor_Theme('Seed_Map_3',[0.00001,Max]);
            Image=fVoxel_Map_3View(Overlay_Image_2,Map_2,{Center,[1,1,1],[2;1;3]},Color_Function,{'Interval',[0,0]});

            fFigure(1,1,1,'',[200,700]);
            set(gcf,'Visible', 'off');
            fSet_Background([200,700],[0,0,0]);
            fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
            imshow(Image);
            title({['FN ',num2str(i)]},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
            fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.6,.15],[0,.3]);
            fColor_Bar(Color_Function,{500,50,80,30,'w','Arial','bold'},round([-1,1]*threshold),'Horizontal Below Inside');

            saveas(gcf,fullfile(Folder,[num2str(i),'.jpg']));
            close(gcf);

            Image=imread(fullfile(Folder,[num2str(i),'.jpg']));
            Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,1458,1,17],[1,1458,400,417]}});
            imwrite(Image,fullfile(Folder,[num2str(i),'.jpg']));
        end
        if isvalid(fig)
            waitbar(1,fig,'Visualizing ...');
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

        if isvalid(fig)
            close(fig);
        end
end


end