function fVisualize_FN_HCP_Surface(File_FN,File_Brain_Surface,Dir_Figure)
% Yuncong Ma, 2/1/2024
% Visualize FN in HCP Surface format for NMF APP
% fVisualize_FN_HCP_Surface(File_FN,File_Brain_Surface,Dir_Figure)
% Work for both group-level and individualized FN

Options.Color_Range_Style='SR-NMF';
Options=fOption('fVisualize_FN_HCP_Surface',Options,varargin);
if isempty(Options)
    return;
end
Color_Range_Style=Options.Color_Range_Style;

load(File_FN,'FN');
Brain_Surface=fLoad_MATLAB_Single_Variable(File_Brain_Surface);

%setting for colorbar
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

k=size(FN,2);
Center=FN';

clear Image

for i=1:k
    % set color range
    switch Color_Range_Style
        case 'SR-NMF'
            threshold=double(round(100*prctile(abs(reshape(Center(i,:),1,[])),99))/100);
            Color_Range=round(100*[threshold/2,threshold])/100;
            Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
        case 'GIG-ICA'
            threshold=double(prctile(Center(i,:),99.8,"all"));
            if threshold>1
                Color_Range=round([0.5,1]*threshold*10)/10;
            elseif threshold>0.01
                Color_Range=round([0.5,1]*threshold*100)/100;
            else
                Color_Range=[0.5,1]*threshold;
            end
            Color_Function=fColor_Theme('Seed_Map_3_Positive',Color_Range);
    end

    % main figure
    Image=fHCP_FS_6View_Column(Brain_Surface.Shape,fHCP_To_32K(squeeze(Center(i,:))',Brain_Surface.MW.LR),fColor_Theme('Seed_Map_3_Positive',round(100*[threshold/2,threshold])/100),...
        {'Material','dull'});

    % fFigure(i,1,1,'',[200,1100],[]);
    figure('Visible','off','NumberTitle','off','Position',[1,1,200,1100]);
    fSet_Background([200,1100],[0,0,0]);
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],1,0.15,[.7,.1],[0,.4]);
    imshow(Image);
    title({['FN ',num2str(i)]},'Fontsize',30,'FontName','Arial','Color',[1 1 1]*0.9,'Fontweight','bold');
    fAxes_Colorbar(1,1,1,'Below',[.9,.8],[0,-.1],2,0.1,[.8,.2],[0,.4]);
    fColor_Bar(Color_Function,Color_Bar.Options);

    saveas(gcf,fullfile(Dir_Figure,[num2str(i),'.jpg']));
    close(gcf);

    Image=imread(fullfile(Dir_Figure,[num2str(i),'.jpg']));
    Image=fAdd_Block_To_Image(Image,{{[0,0,0],[1,2292,1,25],[1,2292,380,417]}});
    imwrite(Image,fullfile(Dir_Figure,[num2str(i),'.jpg']));

    fCrop_Image_File(fullfile(Dir_Figure,[num2str(i),'.jpg']),[0,0.08,1,.93]);
end

% assemble images
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
