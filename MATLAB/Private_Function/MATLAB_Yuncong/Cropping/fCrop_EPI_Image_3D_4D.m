function [EPI_Image_3D_4D,Crop_Parameter]=fCrop_EPI_Image_3D_4D(EPI_Image_3D_4D,Crop_Size,varargin)
% Yuncong Ma, 11/10/2022
% fCrop_EPI_Image_3D_4D is based on the mass center
% [EPI_Image_3D_4D,Crop_Parameter]=fCrop_EPI_Image_3D_4D(EPI_Image_3D_4D,Crop_Size)
% Crop_Size=[48 32 20] or [48 32 inf]
% Crop_Parameter.FOV_Old=[1,Nx;1,Ny;1,Nz];
% Crop_Parameter.FOV=[X_Min,X_Max;Y_Min,Y_Max;Z_Min,Z_Max];
% [EPI_Image_3D_4D,Crop_Parameter]=fCrop_EPI_Image_3D_4D(EPI_Image_3D_4D,Crop_Size,'Center',center)


Size=size(EPI_Image_3D_4D);

for dim=1:3
    if isfinite(Crop_Size(dim)) && Crop_Size(dim)>Size(dim)
        fprintf('Error in fCrop_EPI_Image_3D_4D: in the crop size in the %d-th dimension is larger than the EPI\n',dim);
        return;
    end
end

if length(Size)==4
    EPI_Average=squeeze(mean(EPI_Image_3D_4D,4));
else
    EPI_Average=EPI_Image_3D_4D;
end

if length(varargin)==2
    if strcmp(varargin{1},'Center') && length(varargin{2})==3
        Center=varargin{2};
    else
        fprintf('Error in fCrop_EPI_Image_3D_4D: wrong input for options\n');
    end
elseif ~isempty(varargin)
    fprintf('Error in fCrop_EPI_Image_3D_4D: wrong input for options\n');
else
    Center=fMass_Center(EPI_Average);
end

Crop_Parameter.FOV_Old=[1,Size(1);1,Size(2);1,Size(3)];
Crop_Parameter.FOV=zeros([3,2]);
for dim=1:3
    if isfinite(Crop_Size(dim))
        temp=[max([1,round(Center(dim)-Crop_Size(dim)/2)]),min([Size(dim),round(Center(dim)+Crop_Size(dim)/2)])];
        switch diff(temp)-Crop_Size(dim)+1
            case -2
                switch temp(1)
                    case Size(dim)-Crop_Size(dim)
                        temp=temp+[-2,0];
                    case 1
                        temp=temp+[0,2];
                    otherwise
                        temp=temp+[-1,1];
                end
            case -1
                switch diff(temp)-Crop_Size(dim)
                    case Size(dim)-Crop_Size(dim)
                        temp=temp+[-1,0];
                    case 1
                        temp=temp+[0,1];
                    otherwise
                        temp=temp+[0,1];
                end
            case 2
                temp=temp+[1,-1];
            case 1
                temp=temp+[1,0];
        end
        Crop_Parameter.FOV(dim,:)=temp;
    else
        Crop_Parameter.FOV(dim,:)=[1,Size(dim)];
    end
end

if length(Size)==3
   EPI_Image_3D_4D=EPI_Image_3D_4D(Crop_Parameter.FOV(1,1):Crop_Parameter.FOV(1,2),...
       Crop_Parameter.FOV(2,1):Crop_Parameter.FOV(2,2),...
       Crop_Parameter.FOV(3,1):Crop_Parameter.FOV(3,2)); 
else
   EPI_Image_3D_4D=EPI_Image_3D_4D(Crop_Parameter.FOV(1,1):Crop_Parameter.FOV(1,2),...
       Crop_Parameter.FOV(2,1):Crop_Parameter.FOV(2,2),...
       Crop_Parameter.FOV(3,1):Crop_Parameter.FOV(3,2),:); 
end


end