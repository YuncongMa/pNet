function [y,Mask]=fApply_Mask(Mask,x,varargin)
% By Yuncong Ma, Aug. 8, 2020
% Apply Mask to x, making it into 2D matrix or structure
% y=fApply_Mask(Mask,x)
% Size of Mask and x are the same
% Mask could be [], use [y,Mask]=fApply_Mask() to update Mask
% Mask could be 'Upper' or 'Lower' for ROI-wise FC matrices
% while x should contain square matrix
% y=fApply_Mask(Mask,x,Dimension_Observation)
% Dimension_Observation=1 means the first dimension is the observation
% Dimension_Observation=-1 means the last dimension is the observation
% Dimension_Observation is set to 0 in default, meaning that x is one
% observation

if ~isempty(varargin)
    Dimension_Observation=varargin{1};
    if ~(Dimension_Observation==1 || Dimension_Observation==-1 || Dimension_Observation==0)
        error('Error in fApply_Mask: Dimension_Observation should be either 0, 1 or -1');
    end
else
    Dimension_Observation=0;
end

Size_x=size(x);

if ischar(Mask) && (length(Size_x)<2 || length(Size_x)>3)
    error('Error in fInverse_Mask: dimension of x should be either 2 or 3 when using Upper or Lower');
end
if ischar(Mask) && strcmp(Mask,'Upper')
    if Dimension_Observation==0
        if length(Size_x)~=2 || Size_x(1)~=Size_x(2)
            error('Error in fInverse_Mask: data should be in 2D square matrix, when using %s',Mask);
        end
        Mask=triu(ones(Size_x),1);
    elseif Dimension_Observation==1
        if length(Size_x)~=3 || Size_x(2)~=Size_x(3)
            error('Error in fInverse_Mask: data should be in 3D (square matrix in [2,3]), when using %s',Mask);
        end
        Mask=triu(ones(Size_x(2:3)),1);
    else
        if length(Size_x)~=3 || Size_x(1)~=Size_x(2)
            error('Error in fInverse_Mask: data should be in 3D (square matrix in [1,2]), when using %s',Mask);
        end
        Mask=triu(ones(Size_x(1:2)),1);
    end
elseif ischar(Mask) && strcmp(Mask,'Lower')
    if Dimension_Observation==0
        if length(Size_x)~=2 || Size_x(1)~=Size_x(2)
            error('Error in fInverse_Mask: data should be in 2D square matrix, when using %s',Mask);
        end
        Mask=tril(ones(Size_x),1);
    elseif Dimension_Observation==1
        if length(Size_x)~=3 || Size_x(2)~=Size_x(3)
            error('Error in fInverse_Mask: data should be in 3D (square matrix in [2,3]), when using %s',Mask);
        end
        Mask=tril(ones(Size_x(2:3)),1);
    else
        if length(Size_x)~=3 || Size_x(1)~=Size_x(2)
            error('Error in fInverse_Mask: data should be in 3D (square matrix in [1,2]), when using %s',Mask);
        end
        Mask=tril(ones(Size_x(1:2)),1);
    end
elseif ischar(Mask)
    error('Error in fInverse_Mask: unknow settings for Mask: %s',Mask);
end

if isempty(Mask)
    Size_Mask=Size_x;
    switch Dimension_Observation
        case 0
        case 1
            Size_Mask(1)=1;
        case -1
            Size_Mask(end)=1;
    end
    if length(Size_Mask)>2
        switch Dimension_Observation
            case 0
            case 1
                Size_Mask(1)=[];
            case -1
                Size_Mask(end)=[];
        end
    end
    Mask=ones(Size_Mask);
end

% no change cases
if ~any(Mask(:)==0) && length(size(Mask))==2 && length(size(x))==2 && ...
        ((Dimension_Observation==-1 && size(Mask,2)==1) || (Dimension_Observation==1 && size(Mask,1)==1))
    y=x;
    return
end

Size_Mask=size(Mask);
if Dimension_Observation~=0
    if length(Size_Mask)==2 && Size_Mask(1)==1
        Size_Mask=Size_Mask(2);
    elseif length(Size_Mask)==2 && Size_Mask(2)==1
        Size_Mask=Size_Mask(1);
    end
end

if Dimension_Observation==0
    if ~isequal(Size_x,Size_Mask)
        error('Error in fApply_Mask: unequal size of Mask and x');
    end
    y(1,:)=x(Mask(:)>0);
    return
end

if Dimension_Observation==1
    Dimension=2:length(Size_x);
else
    Dimension=1:length(Size_x)-1;
end

if ~isequal(Size_x(Dimension),Size_Mask)
    error('Error in fApply_Mask: unequal size of Mask(Dimension) and x');
end

if Dimension_Observation==1
    x=reshape(x,Size_x(1),prod(Size_x(2:end)));
    y=x(:,Mask(:)>0);
else
    x=reshape(x,prod(Size_x(1:end-1)),Size_x(end));
    y=x(Mask(:)>0,:);
end

end