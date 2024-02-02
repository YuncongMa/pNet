function fColor_Bar(Color_Function,varargin)
% Designed by Yuncong Ma, 2/1/2024
% fColor_Bar(Color_Function)
% Color_Function is a matrix, each row contains [value, R, G, B]
% value in Color_Function should be sorted from small to large
% Options support 

% fColor_Bar(_,{Width,Height,Interval,Font_Size,Font_Color,Font_Name,Font_Weight},Ticks)
% {Width,Height,Interval,Font_Size,Font_Color,Font_Name,Font_Weight} is {500,100,10,20,'k','Arial','bold'} in default
% {Width,Height,Interval,Font_Size,Font_Color,Font_Name,Font_Weight} could
% be {}, or part of the settings, such as {Width,Height,Interval}, or
% {Width,[],Interval}
% fColor_Bar(_,_,Ticks)
% Ticks is to display, it's derived from Color_Function in default
% fColor_Bar(_,_,_,Style)
% Style is 'Vertical Right' in default, where the Color bar is vertical and the numbers are on the
% right
% Style could be 'Vertical Outside', where the Color bar is vertical and the numbers are on the
% top or bottom of the Color bar, and the Ticks will be set according to
% Style could also be 'Horizontal Below', 'Horizontal Below Inside' or 'Horizontal Outside'
% fColor_Bar(_,_,_,_,Label)
% Label is {Name,Interval,Font_Size,Font_Color,Font_Name,Font_Weight}
% Set to {'',10,20,'k','Arial','bold'} in default

% default setting
Options.Height=50;
Options.Width=500;
Options.Ticks=[];
Options.Tick_Font_Size=20;
Options.Tick_Font_Color='w';
Options.Tick_Font_Name='Arial';
Options.Tick_Font_Weight='bold';
Options.Tick_Interval=10;
Options.Style='Horizontal Below Inside';
Options.Label='';
Options.Label_Font_Size=20;
Options.Label_Font_Color='w';
Options.Label_Font_Name='Arial';
Options.Label_Font_Weight='bold';
Options.Label_Interval=30;

Options=fOption('fColor_Bar',Options,varargin);
if isempty(Options)
    return;
end
    
% color bar
Height=Options.Height;
Width=Options.Width;
% ticks
Style=Options.Style;
Tick_Font_Size=Options.Tick_Font_Size;
Tick_Font_Color=Options.Tick_Font_Color;
Tick_Font_Name=Options.Tick_Font_Name;
Tick_Font_Weight=Options.Tick_Font_Weight;
Tick_Interval=Options.Tick_Interval;
% label
Label=Options.Label;
Label_Font_Size=Options.Label_Font_Size;
Label_Font_Color=Options.Label_Font_Color;
Label_Font_Name=Options.Label_Font_Name;
Label_Font_Weight=Options.Label_Font_Weight;
Label_Interval=Options.Label_Interval;


if isempty(Options.Ticks)
    Ticks=Color_Function([1, size(Color_Function,1)],1);
else
    Ticks=Options.Ticks;
end


if size(Color_Function,2)~=4 || size(Color_Function,1)<2
    error('Error in fColor_Bar: wrong input for Color_Function\n');
end

Color_Range=Color_Function([1,size(Color_Function,1)],1);

switch Style
    
    case {'Vertical Right','Vertical Outside'}
        Step=(Color_Range(2)-Color_Range(1))/Height;
        
        Matrix=zeros(Height+1,Width,3);
        
        for i=0:Height
            temp=Color_Range(1)+i*Step;
            for j=1:size(Color_Function,1)-1
                if Color_Function(j,1)<=temp && Color_Function(j+1,1)>=temp
                    break
                end
            end
            a=temp-Color_Function(j,1);
            b=Color_Function(j+1,1)-temp;
            if a==0
                temp2=Color_Function(j,2:4);
            elseif b==0
                temp2=Color_Function(j+1,2:4);
            else
                temp2=(b*Color_Function(j,2:4)+a*Color_Function(j+1,2:4))/(a+b);
            end
            Matrix(i+1,:,:)=repmat(temp2,[Width,1]);
        end
        Matrix=Matrix(end:-1:1,:,:);
        
        imshow(Matrix);
        axis on
        
        xticks([]);
        yticks([]);
        
        temp=inf;
        if strcmp(Style,'Vertical Outside')
            Ticks=[Color_Function(1,1);Color_Function(end,1)];
        end
        for i=1:length(Ticks)
            if temp==Ticks(i)
                continue;
            end
            temp=double(Ticks(i));
            if strcmp(Style,'Vertical Right')
                text(Width+Tick_Interval,Height-(temp-Color_Range(1))/Step+1,num2str(temp),...
                    'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight);
            else % 'Vertical Outside'
                if i==1 % bottom
                    text(Width/2,Height-(temp-Color_Range(1))/Step+1+Tick_Interval,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','center');
                else % top
                    text(Width/2,Height-(temp-Color_Range(1))/Step+1-Tick_Interval,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','center');
                end
            end
        end
        if ~isempty(Label)
            text(Width+Label_Interval,Height/2,Label,...
                    'Fontsize',Label_Font_Size,'FontName',Label_Font_Name,'Color',Label_Font_Color,'FontWeight',Label_Font_Weight);
        end
        
    case {'Horizontal Below','Horizontal Outside','Horizontal Below Inside'}
        Step=(Color_Range(2)-Color_Range(1))/Width;
        
        Matrix=zeros(Height,Width+1,3);
        
        for i=0:Width
            temp=Color_Range(1)+i*Step;
            for j=1:size(Color_Function,1)-1
                if Color_Function(j,1)<=temp && Color_Function(j+1,1)>=temp
                    break
                end
            end
            a=temp-Color_Function(j,1);
            b=Color_Function(j+1,1)-temp;
            if a==0
                temp2=Color_Function(j,2:4);
            elseif b==0
                temp2=Color_Function(j+1,2:4);
            else
                temp2=(b*Color_Function(j,2:4)+a*Color_Function(j+1,2:4))/(a+b);
            end
            Matrix(:,i+1,:)=reshape(repmat(temp2,[Height,1]),[Height,1,3]);
        end
%         Matrix=Matrix(end:-1:1,:,:);
        
        imshow(Matrix);
        axis on
        
        xticks([]);
        yticks([]);
        
        temp=inf;
        if strcmp(Style,'Horizontal Outside') || strcmp(Style,'Horizontal Below Inside')
            Ticks=[Color_Function(1,1);Color_Function(end,1)];
        end
        for i=1:length(Ticks)
            if temp==Ticks(i)
                continue;
            end
            temp=double(Ticks(i));
            if strcmp(Style,'Horizontal Below')
                text((temp-Color_Range(1))/Step+1,Height+Tick_Interval,num2str(temp),...
                    'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','center');
            elseif strcmp(Style,'Horizontal Below Inside')
                if i==1 % left
                    text((temp-Color_Range(1))/Step+1,Height+Tick_Interval,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','left');
                else % right
                    text((temp-Color_Range(1))/Step+1,Height+Tick_Interval,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','right');
                end
            else % 'Horizontal Outside'
                if i==1 % left
                    text((temp-Color_Range(1))/Step+1-Tick_Interval,Height/2,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','right');
                else % right
                    text((temp-Color_Range(1))/Step+1+Tick_Interval,Height/2,num2str(temp),...
                        'Fontsize',Tick_Font_Size,'FontName',Tick_Font_Name,'Color',Tick_Font_Color,'FontWeight',Tick_Font_Weight,'HorizontalAlignment','left');
                end
            end
        end
        if ~isempty(Label)
            text(Width/2,Height+Label_Interval,Label,...
                    'Fontsize',Label_Font_Size,'FontName',Label_Font_Name,'Color',Label_Font_Color,'FontWeight',Label_Font_Weight,'HorizontalAlignment','center',...
                    'Interpreter','none');
        end
    
    otherwise
        error('Error in fColor_Bar: wrong settings for Style.\n');
end
% set(gcf, 'MenuBar', 'none');
% set(gcf, 'ToolBar', 'none');
end







