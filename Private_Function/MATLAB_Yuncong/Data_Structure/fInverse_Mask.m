function y=fInverse_Mask(Mask,x,varargin)
% By Yuncong Ma, Aug. 8, 2020
% Inverse function of fApply_Mask
% Use the same parameter settings of Mask, Dimension, Dimension_Observation
% in fApply_Mask
% Support matrix and structure
% y=fInverse_Mask(Mask,x)
% Size of Mask and y are the same
% x is a 2D matrix
% the number of 1 in Mask is as the same as numel(x)
% Mask could be 'Upper' or 'Lower' for ROI-wise FC matrices
% y=fInverse_Mask(Mask,x,Dimension_Observation)
% Dimension_Observation=1 means the first dimension is the observation
% Dimension_Observation=-1 means the last dimension is the observation
% Dimension_Observation is set to 0 in default, meaning that x is one
% observation
% Support different data types

Class=class(x);

if length(varargin)==1
    Dimension_Observation=varargin{1};
    if ~(Dimension_Observation==1 || Dimension_Observation==-1 || Dimension_Observation==0)
        error('Error in fInverse_Mask: Dimension_Observation should be either 0, 1 or -1');
    end
else
    Dimension_Observation=0;
end

Size_x=size(x);

if ischar(Mask)
    if Dimension_Observation==0
        Size_Mask=ceil(sqrt(2*Size_x(2)));
    elseif Dimension_Observation==1
        Size_Mask=ceil(sqrt(2*Size_x(2)));
    elseif Dimension_Observation==-1
        Size_Mask=ceil(sqrt(2*Size_x(1)));
    end
end
Mask_Char=0;
if ischar(Mask) && strcmp(Mask,'Upper')
    Mask=triu(ones(Size_Mask),1);
    Mask_Char=1;
elseif ischar(Mask) && strcmp(Mask,'Lower')
    Mask=tril(ones(Size_Mask),1);
    Mask_Char=2;
elseif ischar(Mask)
    error('Error in fInverse_Mask: unknow settings for Mask: %s',Mask);
end

Size_Mask=size(Mask);
Numel_Mask=sum(Mask(:)>0);

% no change cases
if ~any(Mask(:)==0) && length(size(x))==2 && length(Size_Mask)==2 &&...
        ((Dimension_Observation==-1 && Size_Mask(2)==1) || (Dimension_Observation==1 && Size_Mask(1)==1))
    y=x;
    return
end

if length(Size_Mask)==2 && Size_Mask(1)==1
    Size_Mask=Size_Mask(2);
elseif length(Size_Mask)==2 && Size_Mask(2)==1
    Size_Mask=Size_Mask(1);
end

if Dimension_Observation==0
    if Size_x(2)~=Numel_Mask
        error('Error in fInverse_Mask: unequal 1 in Mask and numel(x)');
    end
    if isnumeric(x)
        y=zeros(size(Mask),Class);
        y(Mask(:)>0)=x;
        if Mask_Char>0
            y=y+y';
        end
    else
        temp=fStructure(fieldnames(x));
        String=['y(1:',num2str(size(Mask,1))];
        for i=2:length(size(Mask))
            String=[String,',1:',num2str(size(Mask,i))];
        end
        eval([String,')=temp;']);
        y(Mask(:)>0)=x;
        if Mask_Char>0
            y=fTranspose(y,Dimension_Observation,Mask_Char);
        end
    end
    return
end

if (Dimension_Observation==1 && ~isequal(Size_x(2),Numel_Mask)) || (Dimension_Observation==-1 && ~isequal(Size_x(1),Numel_Mask))
    error('Error in fInverse_Mask: unequal size of Mask(Dimension) and x');
end

if Dimension_Observation==1
    if isnumeric(x)
        y=zeros([Size_x(1),prod(Size_Mask)],Class);
        y(:,Mask(:)>0)=x;
        y=reshape(y,[Size_x(1),Size_Mask]);
        if Mask_Char>0
            y=y+permute(y,[1,3,2]);
        end
    else
        y(1:Size_x(1),1:prod(Size_Mask))=fStructure(fieldnames(x));
        y(:,Mask(:)>0)=x;
        y=reshape(y,[Size_x(1),Size_Mask]);
        if Mask_Char>0
            y=fTranspose(y,Dimension_Observation,Mask_Char);
        end
    end
else
    if isnumeric(x)
        y=zeros([prod(Size_Mask),Size_x(2)],Class);
        y(Mask(:)>0,:)=x;
        y=reshape(y,[Size_Mask,Size_x(2)]);
        if Mask_Char>0
            y=y+permute(y,[2,1,3]);
        end
    else
        y(1:prod(Size_Mask),1:Size_x(2))=fStructure(fieldnames(x));
        y(Mask(:)>0,:)=x;
        y=reshape(y,[Size_Mask,Size_x(2)]);
        if Mask_Char>0
            y=fTranspose(y,Dimension_Observation,Mask_Char);
        end
    end
end

if ~isa(y,Class)
    eval(['y=',Class,'(y);']);
end

end

function Structure=fStructure(Field_Value_Cell)

for i=1:length(Field_Value_Cell)
    Structure.(Field_Value_Cell{i})=[];
end
end

function y=fTranspose(y,Dimension_Observation,Mask_Char)
if Dimension_Observation==0 && Mask_Char==1
    for i=1:size(y,1)-1
        for j=i+1:size(y,1)
            y(j,i)=y(i,j);
        end
    end
elseif Dimension_Observation==0 && Mask_Char==2
    for i=1:size(y,1)-1
        for j=i+1:size(y,1)
            y(i,j)=y(j,i);
        end
    end
elseif Dimension_Observation==1 && Mask_Char==1
    for i=1:size(y,2)-1
        for j=i+1:size(y,2)
            for k=1:size(y,1)
                y(k,j,i)=y(k,i,j);
            end
        end
    end
elseif Dimension_Observation==1 && Mask_Char==2
    for i=1:size(y,2)-1
        for j=i+1:size(y,2)
            for k=1:size(y,1)
                y(k,i,j)=y(k,j,i);
            end
        end
    end
elseif Dimension_Observation==-1 && Mask_Char==1
    for i=1:size(y,2)-1
        for j=i+1:size(y,2)
            for k=1:size(y,3)
                y(j,i,k)=y(i,j,k);
            end
        end
    end
elseif Dimension_Observation==-1 && Mask_Char==2
    for i=1:size(y,2)-1
        for j=i+1:size(y,2)
            for k=1:size(y,3)
                y(i,j,k)=y(j,i,k);
            end
        end
    end
end

end



