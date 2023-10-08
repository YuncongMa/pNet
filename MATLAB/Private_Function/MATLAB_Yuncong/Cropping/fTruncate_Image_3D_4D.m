function [Image_3D_4D_Truncated,Center,Crop_Parameter]=fTruncate_Image_3D_4D(Image_3D_4D,Voxel_Size,Extend)
% By Yuncong Ma, Feb. 16, 2018
% fTruncate_Image_3D_4D is based on the mask applied into the Image_3D_4D
% Truncated Data according to the mask, making the image file smaller
% [Image_3D_4D_Truncated,Center,Crop_Parameter]=fTruncate_Image_3D_4D(Image_3D_4D,Voxel_Size,Extend)
% Voxel_Size=[Vx, Vy, Vz];
% Extend=[Ex,Ey,Ez]; could be postive or inf, Ez could be negative
% Center is used for alignment or motion correction in imregister
% Crop_Parameter contains the original data size in FOV_Old and new size in
% FOV
% fApply_Cropped_FOV could be used to apply the same cropping parameters to
% other data
% fInverse_Crop_EPI_Image_3D_4D could be used to reconstruct the original
% data

Size=size(Image_3D_4D);
Crop_Parameter.FOV_Old=[1,Size(1);1,Size(2);1,Size(3)];

if length(size(Image_3D_4D))==4
    Mask=Image_3D_4D(:,:,:,1)>0;
else
    Mask=Image_3D_4D>0;
end

FOV=zeros(3,2);

temp=squeeze(sum(sum(Mask,2),3))>0;
FOV(1,:)=[find(temp,1),Size(1)-find(temp(end:-1:1),1)+1];

temp=squeeze(sum(sum(Mask,1),3))>0;
FOV(2,:)=[find(temp,1),Size(2)-find(temp(end:-1:1),1)+1];

temp=squeeze(sum(sum(Mask,1),2))>0;
FOV(3,:)=[find(temp,1),Size(3)-find(temp(end:-1:1),1)+1];

for dim=1:3
    if ~isfinite(Extend(dim))
        FOV(dim,1:2)=[1,Size(dim)];
    else
       FOV(dim,:)=[max([1,FOV(dim,1)-Extend(dim)]),min([Size(dim),FOV(dim,2)+Extend(dim)])]; 
    end
end

if Extend(3)<0 % Re-evaluate
    Mask(:,:,setdiff(1:Size(3),FOV(3,1):FOV(3,2)))=0;
    
    temp=squeeze(sum(sum(Mask,2),3))>0;
    FOV(1,:)=[find(temp,1),Size(1)-find(temp(end:-1:1),1)+1];
    
    temp=squeeze(sum(sum(Mask,1),3))>0;
    FOV(2,:)=[find(temp,1),Size(2)-find(temp(end:-1:1),1)+1];
    
    temp=squeeze(sum(sum(Mask,1),2))>0;
    FOV(3,:)=[find(temp,1),Size(3)-find(temp(end:-1:1),1)+1];
    
    for dim=1:3
        if ~isfinite(Extend(dim))
            FOV(dim,1:2)=[1,Size(dim)];
        else
            FOV(dim,:)=[max([1,FOV(dim,1)-Extend(dim)]),min([Size(dim),FOV(dim,2)+Extend(dim)])];
        end
    end
end

Crop_Parameter.FOV=FOV;

Center=Voxel_Size.*(mean(FOV,2)'-(1+Size(1:3))/2);

if length(Size)==4
    Image_3D_4D_Truncated=Image_3D_4D(FOV(1,1):FOV(1,2),FOV(2,1):FOV(2,2),FOV(3,1):FOV(3,2),:);
else
    Image_3D_4D_Truncated=Image_3D_4D(FOV(1,1):FOV(1,2),FOV(2,1):FOV(2,2),FOV(3,1):FOV(3,2));
end

end