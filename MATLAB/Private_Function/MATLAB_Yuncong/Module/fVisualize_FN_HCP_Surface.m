function fVisualize_FN_HCP_Surface(File_FN,File_Brain_Surface,Dir_Figure)
% Yuncong Ma, 2/22/2023
% Visualize FN in HCP Surface format for NMF APP
% fVisualize_FN_HCP_Surface(File_FN,File_Brain_Surface,Dir_Figure)
% Work for both group-level and individualized FN

load(File_FN,'FN');
Brain_Surface=fLoad_MATLAB_Single_Variable(File_Brain_Surface);

k=size(FN,2);
Center=FN';

clear Image

for i=1:k
    threshold=double(round(100*prctile(abs(reshape(Center(i,:),1,[])),99))/100);
    Image=fHCP_FS_6View_Column(Brain_Surface.Shape,fHCP_To_32K(squeeze(Center(i,:))',Brain_Surface.MW.LR),fColor_Theme('Seed_Map_3_Positive',round(100*[threshold/2,threshold])/100),...
        {'Material','dull'});

    % fFigure(i,1,1,'',[200,1100],[]);
    figure('Visible','off','NumberTitle','off','Position',[1,1,200,1100]);
    fSet_Background([200,1100],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],1,0.15,[.7,.1],[0,.4]);
    imshow(Image);
    title({['FN ',num2str(i)],''},'Fontsize',40,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.15,[.7,.1],[0,.4]);
    fColor_Bar(fColor_Theme('Seed_Map_3_Positive',round([threshold/2,threshold]*100)),{500,50,100,30,'w','Arial','bold'},round([-1,1]*threshold*100),'Horizontal Below Inside');

    saveas(gcf,fullfile(Dir_Figure,[num2str(i),'.jpg']));
    close(gcf);

    Image=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
    Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,2292,1,25],[1,2292,380,417]}});
    imwrite(Image,fullfile(Dir_Figure,[num2str(i),'.jpg']));

    fCrop_Image_File(fullfile(Dir_Figure,[num2str(i),'.jpg']),[0,0.08,1,.93]);
end
clear Center Index

clear Image
for i=1:k
    Image{i}=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
end
N_Column=10;
N_Row=ceil(k/N_Column);
for i=k+1:(N_Row*N_Column)
    Image{i}=Image{1}*0;
end

imwrite(fAssemble_Image(reshape(Image,[N_Column,N_Row])',{'Interval',[50,5]}),fullfile(Dir_Figure,'All.jpg'));
end
