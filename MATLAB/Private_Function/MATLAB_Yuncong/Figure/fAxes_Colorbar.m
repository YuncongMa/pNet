function varargout=fAxes_Colorbar(Panel_Index,N_Total,Matrix_Size,Style,Shrink_Size,Offset,Figure_Colorbar,Ratio,Shrink_Colorbar,Offset_Colorbar)
% Designed by Yuncong Ma, Jun. 7, 2020
% fAxes_Colorbar is to show paired color bars with corresponding main figures
% fAxes_Colorbar(Panel_Index,N_Total,Matrix_Size,Style,Shrink_Size,Offset,...
% Figure_Colorbar,Ratio,Shrink_Colorbar,Offset_Colorbar)
% Panel_Index is the index of current panel to plot
% Panel_Index could be either the linear index or [Row_Index,Column_Index]
% N_Total is the total number of panels
% if Panel_Index is one number, then a new figure panel will be generated
% when Panel_Index>1 && mod(Panel_Index-1,N_Total)==0.
% Matrix_Size is the way to organize panels
% it could be 0, N_Row, [0,N_Column] [N_Row,N_Column]
% Style is 'Right' or 'Below'
% Shrink_Size is for smallen the panel to add some margin, it contains two
% parameters [Size_X,Size_Y], both of them should be ranging (0 1]
% Offset is to move each subfigure offset from the center position
% Offset contain [Offset_X,Offset_Y]
% if Offset_X is negative, the subfigure panel will be moved to left
% it should range from [-(1-Size_X), (1-Size_X)/2]
% Figure_Colorbar is 1 or 2, main figure or the color bar
% Ratio is the porportion of width used for color bar
% Shrink_Colorbar is the shrink size of color bar
% Offset_Colorbar is the offset to place the color bar
% eg. fAxes_Colorbar(3,20,[5 4],'Right',[.9,.8],[0,-.05],1,0.15,[.9,.6],[-.04,0])

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
    error('Error in fAxes_Colorbar: Wrong Setups in Matrix_Size: %f\n',Matrix_Size);
end
if N_Row<=0||N_Column<=0
    error('Error in fAxes_Colorbar: Wrong Setups in Matrix_Size: %f\n',Matrix_Size);
end
if length(Panel_Index)==2
    Panel_Index=(Panel_Index(1)-1)*N_Column+Panel_Index(2);
end

% new figure
if Panel_Index>1
    if mod(Panel_Index-1,N_Total)==0
        figure;
    end
    Panel_Index=mod(Panel_Index-1,N_Total)+1;
end


if sum((abs(Offset)-(1-Shrink_Size)/2)<0.0001)~=2% in case of calculation error in double
    error('Errior in fAxes_Colorbar: Inappropriate settings for Offset\n');
end

switch Style
    case 'Right'
        Panel_Position=[mod(Panel_Index-1,N_Column)/N_Column,1-1/N_Row*(1+floor((Panel_Index-1)/N_Column))];
        if Figure_Colorbar==1% main figure
            Panel_Size=[1/N_Column*(1-Ratio),1/N_Row];
            Panel_Axes=[Panel_Position+((1-Shrink_Size)./2+Offset).*Panel_Size,Panel_Size.*Shrink_Size];
        else% color bar
            Panel_Size=[1/N_Column*Ratio,1/N_Row];
            Panel_Position=Panel_Position+[1/N_Column*(1-Ratio),0];
            Panel_Axes=[Panel_Position+((1-Shrink_Colorbar)./2+Offset_Colorbar).*Panel_Size,Panel_Size.*Shrink_Colorbar];
        end
    
    case 'Below'
        Panel_Position=[mod(Panel_Index-1,N_Column)/N_Column,1-1/N_Row*(1+floor((Panel_Index-1)/N_Column))]+[0,1/N_Row*Ratio];
        if Figure_Colorbar==1% main figure
            Panel_Size=[1/N_Column,1/N_Row*(1-Ratio)];
            Panel_Axes=[Panel_Position+((1-Shrink_Size)./2+Offset).*Panel_Size,Panel_Size.*Shrink_Size];
        else% color bar
            Panel_Size=[1/N_Column,1/N_Row*Ratio];
            Panel_Position=Panel_Position+[0,-1/N_Row*Ratio];
            Panel_Axes=[Panel_Position+((1-Shrink_Colorbar)./2+Offset_Colorbar).*Panel_Size,Panel_Size.*Shrink_Colorbar];
        end
    
    otherwise
        error('Error in fAxes_Colorbar: unsupported style %s\n',Style);
end

if nargout==0
    axes('Position',Panel_Axes);
else
    varargout={Panel_Axes};
end

end




