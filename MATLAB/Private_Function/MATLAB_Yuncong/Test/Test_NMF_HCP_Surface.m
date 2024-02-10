% Yuncong Ma, 6/6/2023
% Test computation parts in NMF APP

%% Use HCP data in mat format to test
App_Dir='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/MATLAB_APP/NMF';
Work_Dir='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/NMF/Result/HCP/Surface/Test_FN17';

[Flag,Message]=fCompute_gFN(App_Dir,Work_Dir);

[Flag,Message]=fVisualize_gFN(App_Dir,Work_Dir);

[Flag,Message]=fCompute_pFN(App_Dir,Work_Dir);

[Flag,Message]=fVisualize_pFN(App_Dir,Work_Dir);


% [FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,zeros(59412,10));

%% Test_FN17_200_Combined
load('/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/Template/HCP/My_Template/HCP_FS.mat','HCP_FS');
App_Dir='/Volumes/Scratch_0/pNet';
Work_Dir='/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined';

% Re-run the pFN
[Flag,Message]=fCompute_pFN(App_Dir,Work_Dir);
[Flag,Message]=fVisualize_pFN(App_Dir,Work_Dir);

Data_Behavior=fLoad_MATLAB_Single_Variable('/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/Behavior_Data/HCP_1200/Gender_Binary.mat');
% Data_Behavior=Data_Behavior(1:200*4);
Data_Behavior=Data_Behavior(1:4:200*4);

FID=fopen('/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/Behavior_Data/HCP_1200/Gender_200sub_combined.txt','w');
for i=1:200
    fprintf(FID,'%d\n',num2str(Data_Behavior(i)));
end
fclose(FID);

% [Flag,Message]=fCompute_gFN(App_Dir,Work_Dir);
% 
% [Flag,Message]=fVisualize_gFN(App_Dir,Work_Dir);
% 
% [Flag,Message]=fCompute_pFN(App_Dir,Work_Dir);
% 
% [Flag,Message]=fVisualize_pFN(App_Dir,Work_Dir);

File_Info=Search_And_Sort_Files(fullfile(Work_Dir,'Personalized_FN'),'FN\.mat',{{'Subject',-2},{'Scan',-1}});

File_List={};
for i=1:200
    File_List{i}=File_Info(i).Dir;
end
Result=fAPP_Statistics(File_List,{'Behavior',Data_Behavior},{'Method','Wilcoxon rank sum test'},{'FDR',1});
% Result=fAPP_Statistics(File_List,{'Behavior',Data_Behavior},{'Method','2 sample ttest'},{'FDR',1});


Z_Value=zeros(K,N_Dim);
for i=1:K
    for j=1:N_Dim
        [P_Value(i,j),~,stats]=ranksum(squeeze(Data1(:,i,j)),squeeze(Data2(:,i,j)));
        Z_Value(i,j)=stats.zval;
    end
end

for k=1:17
    temp=Result.P_Value(:,k);
    temp(isnan(temp))=1;
    Map=Result.Z_Value(:,k);
    Map(isnan(Map))=0;
    Map(temp>0.1)=0;
    Max=max([1,prctile(abs(Map(:)),98)]);
    Color_Function=fColor_Theme('Seed_Map_3',[0.1,Max]);
    Image_Data=fHCP_FS_6View_Column(HCP_FS.Shape,fHCP_To_32K(Map,HCP_FS.MW.LR),Color_Function,...
        {'Material','dull'});
    fFigure(1,1,1,'',[200,700]);
%     set(gcf,'Visible', 'off');
    fSet_Background([200,700],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
    imshow(Image_Data);
    title({'Significance'},'Fontsize',30,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
%     fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.8,.15],[0,.3]);
%     fColor_Bar(Color_Function(2:end,:),{500,50,100,30,'w','Arial','bold'},[0,1],'Horizontal Below');
    saveas(gcf,fullfile(Work_Dir,'Statistics',['Gender_',num2str(k),'.jpg']));
    close(gcf);
end


%%
fVisualize_Statistics('/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined/Statistics/Test_3_Gender_200sub_combined',{'Folder_Result','/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined'});


%% Visualization with functional atlas border
App_Dir='/Volumes/Scratch_0/pNet';
Work_Dir='/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined';



fVisualize_Statistics('/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined/Statistics/Test_3_Gender_200sub_combined',...
    {'Folder_Result','/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined'},{'Flag_Atlas',1});

%% Visualize statistics with FN boundary
App_Dir='/Volumes/Scratch_0/pNet';
Work_Dir='/Volumes/Scratch_0/Example/HCP_Surface/Test_FN17_200_Combined';

load(fullfile(Work_Dir,'Statistics/Test_4_Significance/Statistics.mat'));
Statistics_FN=Statistics;
load(fullfile(Work_Dir,'Statistics/Test_3_Gender_200sub_combined/Statistics.mat'));
load(fullfile(Work_Dir,'Load_Data/Brain_Surface.mat'));

Threshold_p=0.05;
Threshold_t=30;

for k=1:17
    Atlas_k=Statistics_FN.Result.T_Value(:,k)>Threshold_t;
    Atlas_k=fHCP_To_32K(Atlas_k,Brain_Surface.MW.LR);
    Border=1-Atlas_k*0;
    Face=Brain_Surface.Shape.LR.faces;
    for k2=1:length(Face)
        ps=Face(k2,:);
        ps2=Atlas_k(ps);
        if length(unique(ps2))>1
            Border(ps)=0;
        end
    end

    temp=Statistics.Result.P_Value(:,k);
    temp(isnan(temp))=1;
    Map=Statistics.Result.Z_Value(:,k).*(Statistics_FN.Result.T_Value(:,k)>Threshold_t);
    Map(isnan(Map))=0;
    Map(temp>Threshold_p)=0;
    Max=max([1,prctile(abs(Map(:)),98)]);
    Color_Function=fColor_Theme('Seed_Map_3',[0.1,Max]);
    Image_Data=fHCP_FS_6View_Column(HCP_FS.Shape,fHCP_To_32K(Map,HCP_FS.MW.LR),Color_Function,...
        {'Material','dull'},{'Atlas',Border},{'Atlas_Color',[1,1,1]/4});
    fFigure(1,1,1,'',[200,700]);
%     set(gcf,'Visible', 'off');
    fSet_Background([200,700],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[1,.9],[0,-.05],1,0.1,[.7,.1],[0,.3]);
    imshow(Image_Data);
    title(['FN ',num2str(k)],'Fontsize',30,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
%     fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.8,.15],[0,.3]);
%     fColor_Bar(Color_Function(2:end,:),{500,50,100,30,'w','Arial','bold'},[0,1],'Horizontal Below');
    saveas(gcf,fullfile(Work_Dir,'Statistics',['Gender_',num2str(k),'.jpg']));
    close(gcf);
end























