function Image=fAssemble_Image(Image_Cell,varargin)
% By Yuncong Ma, Sep. 11, 2018
% Image=fAssemble_Image(Image_Cell)
% Image_Cell is a cell type matrix
% images in Image_Cell is RGB colored
% h=fAssemble_Image(_,Options)
% {'Background','w'} the default background color is 'k'
% {'Interval',[Interval_X,Interval_Y]} to separate each subfigure

Options.Background='k';
Options.Interval=[2,2];
Options=fOption('fAssemble_Image',Options,varargin);
if isempty(Options)
    return;
end
if strcmp(Options.Background,'k')
    BackgroundColor=0;
else
    BackgroundColor=1;
end

Interval=Options.Interval;

Image=[];
ps=[0,0];
Size=size(Image_Cell);
for x=1:Size(1)
    for y=1:Size(2)
        temp=Image_Cell{x,y};
        if x==1 && y==1
            Image=temp;
        else
            Image(ps(1)+1:ps(1)+size(temp,1),ps(2)+1:ps(2)+size(temp,2),:)=temp;
        end
        
        
        if x<Size(1)
            Image(ps(1)+1+size(temp,1):ps(1)+size(temp,1)+Interval(1),...
                ps(2)+1:ps(2)+size(temp,2),:)=BackgroundColor;
        end
        if y<Size(2)
            Image(ps(1)+1:ps(1)+size(temp,1),...
                ps(2)+1+size(temp,2):ps(2)+size(temp,2)+Interval(2),:)=BackgroundColor;
        end
        if x<Size(1)&&y<Size(2)
            Image(ps(1)+1+size(temp,1):ps(1)+size(temp,1)+Interval(1),...
                ps(2)+1+size(temp,2):ps(2)+size(temp,2)+Interval(2),:)=BackgroundColor;
        end
        
        ps(2)=ps(2)+Interval(2)+size(temp,2);
    end
    ps(1)=ps(1)+Interval(1)+size(temp,1);
    ps(2)=0;
end


end