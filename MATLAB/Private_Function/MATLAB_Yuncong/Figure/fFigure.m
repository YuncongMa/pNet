function fFigure(Figure_Number,Total_Figure_Number,Organization_Matrix,Figure_Title,Figure_Size,varargin)
% By Yuncong Ma, 1/18/2022
% fFigure(Figure_Number,Total_Figure_Number,Organization_Matrix,Figure_Title,Figure_Size)
% Organization_Matrix=0 do automatic organization
% Organization_Matrix=5 means 5 rows
% Organization_Matrix=[0 5] means 5 columns
% Organization_Matrix=[5 5] means 5 by 5 organization
% Figure_Title could be '', then Figure_Number will be used as the title
% Figure_Size is the pixel size of each figure
% fFigure disable menubar and toolbar by default
% fFigure(~,~,~,~,~,Options) supports Monitor, MenuBar, ToolBar, Visible
% Sidecar
% Monitor is the choice of Monitor for displaying figure
% Monitor is set to 1 in default, and use integer to use the n-th
% monitor
% MenuBar and TooBar are set to 'off' in default
% Visible is set to 'on' in default
% fFigure(1,1,_,Figure_Title,Figure_Size) will place the figure in the
% center of the screen

% Assume figure title bar takes 20 pixels in height

Options.Monitor=1;
Options.MenuBar='off';
Options.ToolBar='off';
Options.Visible='on';
Options=fOption('fFigure',Options,varargin);
if isempty(Options)
    return;
end

Monitor=get(0,'MonitorPositions');
if size(Monitor,1)<Options.Monitor || Options.Monitor<1
    error('fFigure: cannot find %d-th monitor\n',Options.Monitor);
else
    Screen_Size=Monitor(Options.Monitor,3:4);
    Monitor_Offset=Monitor(Options.Monitor,1:2)-[1,1];
end

if Figure_Number==1 && Total_Figure_Number==1
    Figure_Position=Screen_Size/2-Figure_Size/2+Monitor_Offset;
    figure('Name',Figure_Title,'NumberTitle','off','Position',[Figure_Position,Figure_Size]);
    return
end

% N_Row=0;
% N_Column=0;
if length(Organization_Matrix)==1 && Organization_Matrix==0
    N_Row=floor(sqrt(Total_Figure_Number));
    N_Column=ceil(Total_Figure_Number/N_Row);
elseif length(Organization_Matrix)==1 && Organization_Matrix>0
    N_Row=round(Organization_Matrix);
    N_Column=ceil(Total_Figure_Number/N_Row);
elseif length(Organization_Matrix)==2 && Organization_Matrix(1)==0
    N_Column=round(Organization_Matrix(2));
    N_Row=ceil(Total_Figure_Number/N_Column);
elseif length(Organization_Matrix)==2 && Organization_Matrix(1)>0 && Organization_Matrix(2)>0
    N_Row=round(Organization_Matrix(1));
    N_Column=round(Organization_Matrix(2));
else
    fprintf('fFigure: Wrong Setups in Organization_Matrix: %f\n',Organization_Matrix);
    return;
end
if N_Row<=0||N_Column<=0
    fprintf('fFigure: Wrong Setups in Organization_Matrix: %f\n',Organization_Matrix);
    return;
end


if isempty(Figure_Title)
    Figure_Title=['Figure ',num2str(Figure_Number)];
end
if isempty(Figure_Size)
    Figure_Size=round(Screen_Size./[N_Column,N_Row])-[0 20];
end
Step_Size=round((Screen_Size-Figure_Size)./[max([1,N_Column-1]),max([1,N_Row-1])]);
% Step_Size=[min([Figure_Size(1),Step_Size(1)]),min([Figure_Size(2),Step_Size(2)])];

Column_Number=mod(Figure_Number-1,N_Column)+1;
Row_Number=ceil(Figure_Number/N_Column);

Figure_Position=[Column_Number-1,N_Row-Row_Number].*Step_Size+Monitor_Offset;

figure('Name',Figure_Title,'NumberTitle','off','Position',[Figure_Position,Figure_Size]);

if strcmpi(Options.MenuBar,'off')
    set(gcf, 'MenuBar', 'none');
end
if strcmpi(Options.ToolBar,'off')
    set(gcf, 'ToolBar', 'none');
end
if strcmpi(Options.Visible,'off')
    set(gcf,'visible','off');
end

end