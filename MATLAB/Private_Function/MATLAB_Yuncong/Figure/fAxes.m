function varargout=fAxes(Panel_Index,N_Total,Matrix_Size,varargin)
% Designed by Yuncong Ma, Jun. 7, 2020
% fAxes(Panel_Index,N_Total,Matrix_Size)
% Panel_Index is the index of current panel to plot
% Panel_Index could be either the linear index or [Row_Index,Column_Index]
% N_Total is the total number of panels
% if Panel_Index is one number, then a new figure panel will be generated
% when Panel_Index>1 && mod(Panel_Index-1,N_Total)==0.
% 
% Matrix_Size is the way to organize panels
% it could be 0, N_Row, [0,N_Column] [N_Row,N_Column]
% eg. fAxes(2,20,[0,4]),fAxes([3,2],20,5), fAxes(4,20,[5,4]), fAxes(5,20,[4,4])
% fAxes(Panel_Index,N_Total,Organization_Size,Panel_Size)
% Panel_Size is for smallen the panel to add some margin, it contains two
% parameters [Size_X,Size_Y], both of them should be ranging (0 1]
% eg. fAxes(Panel_Index,N_Total,Organization_Size,[.8,.6])
% fAxes(Panel_Index,N_Total,Organization_Size,Panel_Size,Offset)
% Offset is to move each subfigure offset from the center position
% Offset contain [Offset_X,Offset_Y]
% if Offset_X is negative, the subfigure panel will be moved to left
% it should range from [-(1-Size_X), (1-Size_X)/2]
% eg. fAxes(Panel_Index,N_Total,Organization_Size,[.8,.6],[0,-0.15])

% N_Row=0;
% N_Column=0;
if length(Matrix_Size)==1 && Matrix_Size==0
    N_Row=floor(sqrt(N_Total));
    N_Column=ceil(N_Total/N_Row);
elseif length(Matrix_Size)==1 && Matrix_Size>0
    N_Row=round(Matrix_Size);
    N_Column=ceil(N_Total/N_Row);
elseif length(Matrix_Size)==2 && Matrix_Size(1)==0
    N_Column=round(Matrix_Size(2));
    N_Row=ceil(N_Total/N_Column);
elseif length(Matrix_Size)==2 && Matrix_Size(1)>0 && Matrix_Size(2)>0
    N_Row=round(Matrix_Size(1));
    N_Column=round(Matrix_Size(2));
else
    error('Error in fAxes: Wrong Setups in Matrix_Size: %f\n',Matrix_Size);
end
if N_Row<=0||N_Column<=0
    error('Error in fAxes: Wrong Setups in Matrix_Size: %f\n',Matrix_Size);
end

% new figure
if length(Panel_Index)==2
    Panel_Index=(Panel_Index(1)-1)*N_Column+Panel_Index(2);
end
if Panel_Index>1
    if mod(Panel_Index-1,N_Total)==0
        figure;
    end
    Panel_Index=mod(Panel_Index-1,N_Total)+1;
end


Panel_Size=[1/N_Column,1/N_Row];
if length(varargin)==1 && length(varargin{1})==2 && sum(0<varargin{1} & varargin{1}<=1)==2
    Shrink_Size=varargin{1};
    Offset=[0 0];
elseif length(varargin)==2 && length(varargin{1})==2 && sum(0<varargin{1} & varargin{1}<=1)==2 ...
    && length(varargin{2})==2 && sum(0<varargin{1} & varargin{1}<=1)==2
    Shrink_Size=varargin{1};
    Offset=varargin{2};
    if sum((abs(Offset)-(1-Shrink_Size)/2)<0.0001)~=2% in case of calculation error in double
        error('Errior in fAxes: Inappropriate settings for Offset\n');
    end
elseif ~isempty(varargin)
    error('Errior in fAxes: Wrong settings for optional settings\n');
else
    Shrink_Size=[1 1];
end

Panel_Position=[mod(Panel_Index-1,N_Column)/N_Column,1-1/N_Row*(1+floor((Panel_Index-1)/N_Column))];

if nargout==0
    axes('Position',[Panel_Position+((1-Shrink_Size)./2+Offset).*Panel_Size,Panel_Size.*Shrink_Size]);
else
    varargout={[Panel_Position+((1-Shrink_Size)./2+Offset).*Panel_Size,Panel_Size.*Shrink_Size]};
end

end




