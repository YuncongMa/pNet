function fColor_Bar(Color_Function,varargin)
% Designed by Yuncong Ma, 1/19/2022
% fColor_Bar(Color_Function)
% Color_Function is a matrix, each row contains [value, R, G, B]
% value in Color_Function should be sorted from small to large
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


Height=500;
Width=100;
Interval=10;
Font_Size=20;
Font_Color='k';
Font_Name='Arial';
Font_Weight='bold';
Style='Vertical Right';
    
if ~isempty(varargin)
    if length(varargin{1})>7
        error('Error in fColor_Bar: wrong settings for {Width,Height,Interval,Font_Size,Font_Name,Font_Weight}.\n');
    end
    if ~isempty(varargin{1}) && ~isempty(varargin{1}{1})
        Width=varargin{1}{1};
    end
    if length(varargin{1})>1 && ~isempty(varargin{1}{2})
        Height=varargin{1}{2};
    end
    if length(varargin{1})>2 && ~isempty(varargin{1}{3})
        Interval=varargin{1}{3};
    end
    if length(varargin{1})>3 && ~isempty(varargin{1}{4})
        Font_Size=varargin{1}{4};
    end
    if length(varargin{1})>4 && ~isempty(varargin{1}{5})
        Font_Color=varargin{1}{5};
    end
    if length(varargin{1})>5 && ~isempty(varargin{1}{6})
        Font_Name=varargin{1}{6};
    end
    if length(varargin{1})>6 && ~isempty(varargin{1}{7})
        Font_Weight=varargin{1}{7};
    end
    
    if length(varargin)>1
        if ~isempty(varargin{2})
            Ticks=varargin{2};
        else
            Ticks=Color_Function(:,1);
        end
        if length(varargin)>2
            Style=varargin{3};
        end
    else
        Ticks=Color_Function(:,1);
    end
end

if size(Color_Function,2)~=4 || size(Color_Function,1)<2
    error('Error in fColor_Bar: wrong input for Color_Function\n');
end

Color_Range=Color_Function([1,size(Color_Function,1)],1);

Label.Name='';
Label.Interval=10;
Label.Font_Size=20;
Label.Font_Color='k';
Label.Font_Name='Arial';
Label.Font_Weight='bold';

if length(varargin)==4
    if length(varargin{4})>6
        error('Error in fColor_Bar: wrong settings for {Name,Interval,Font_Size,Font_Color,Font_Name,Font_Weight}.\n');
    end
    if ~isempty(varargin{4}) && ~isempty(varargin{4}{1})
        Label.Name=varargin{4}{1};
    end
    if length(varargin{4})>1 && ~isempty(varargin{4}{2})
        Label.Interval=varargin{4}{2};
    end
    if length(varargin{4})>2 && ~isempty(varargin{4}{3})
        Label.Font_Size=varargin{4}{3};
    end
    if length(varargin{4})>3 && ~isempty(varargin{4}{4})
        Label.Font_Color=varargin{4}{4};
    end
    if length(varargin{4})>4 && ~isempty(varargin{4}{5})
        Label.Font_Name=varargin{4}{5};
    end
    if length(varargin{4})>5 && ~isempty(varargin{4}{6})
        Label.Font_Weight=varargin{4}{6};
    end
end

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
                text(Width+Interval,Height-(temp-Color_Range(1))/Step+1,num2str(temp),...
                    'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight);
            else % 'Vertical Outside'
                if i==1 % bottom
                    text(Width/2,Height-(temp-Color_Range(1))/Step+1+Interval,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','center');
                else % top
                    text(Width/2,Height-(temp-Color_Range(1))/Step+1-Interval,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','center');
                end
            end
        end
        if ~isempty(Label.Name)
            text(Width+Label.Interval,Height/2,Label.Name,...
                    'Fontsize',Label.Font_Size,'FontName',Label.Font_Name,'Color',Label.Font_Color,'FontWeight',Label.Font_Weight);
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
                text((temp-Color_Range(1))/Step+1,Height+Interval,num2str(temp),...
                    'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','center');
            elseif strcmp(Style,'Horizontal Below Inside')
                if i==1 % left
                    text((temp-Color_Range(1))/Step+1,Height+Interval,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','left');
                else % right
                    text((temp-Color_Range(1))/Step+1,Height+Interval,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','right');
                end
            else % 'Horizontal Outside'
                if i==1 % left
                    text((temp-Color_Range(1))/Step+1-Interval,Height/2,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','right');
                else % right
                    text((temp-Color_Range(1))/Step+1+Interval,Height/2,num2str(temp),...
                        'Fontsize',Font_Size,'FontName',Font_Name,'Color',Font_Color,'FontWeight',Font_Weight,'HorizontalAlignment','left');
                end
            end
        end
        if ~isempty(Label.Name)
            text(Width/2,Height+Label.Interval,Label.Name,...
                    'Fontsize',Label.Font_Size,'FontName',Label.Font_Name,'Color',Label.Font_Color,'FontWeight',Label.Font_Weight,'HorizontalAlignment','center',...
                    'Interpreter','none');
        end
    
    otherwise
        error('Error in fColor_Bar: wrong settings for Style.\n');
end
% set(gcf, 'MenuBar', 'none');
% set(gcf, 'ToolBar', 'none');
end







