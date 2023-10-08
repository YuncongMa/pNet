function Image_3D_4D=fApply_Cropped_FOV(Image_3D_4D,Crop_Parameter)
% Yuncong Ma, 11/10/2022
% Use predefined Crop_Parameter to crops a 3D or 4D image
% Image_3D_4D=fApply_Cropped_FOV(Image_3D_4D,Crop_Parameter)
% works for data have the same FOV and same resolution
% Crop_Parameter.FOV_Old is the orginal FOV: [x_min,x_max;y_min,y_max;z_min,z_max]
% Crop_Parameter.FOV is for the new FOV

Size=size(Image_3D_4D);
if Crop_Parameter.FOV_Old(1,2)-Crop_Parameter.FOV_Old(1,1)+1~=Size(1) ||...
        Crop_Parameter.FOV_Old(2,2)-Crop_Parameter.FOV_Old(2,1)+1~=Size(2) ||...
        Crop_Parameter.FOV_Old(3,2)-Crop_Parameter.FOV_Old(3,1)+1~=Size(3)
    fprintf('Error in fApply_Cropped_FOV: the old FOV does not match to the data\n');
end

if length(Size)==3
   Image_3D_4D=Image_3D_4D(Crop_Parameter.FOV(1,1):Crop_Parameter.FOV(1,2),...
       Crop_Parameter.FOV(2,1):Crop_Parameter.FOV(2,2),...
       Crop_Parameter.FOV(3,1):Crop_Parameter.FOV(3,2)); 
else
   Image_3D_4D=Image_3D_4D(Crop_Parameter.FOV(1,1):Crop_Parameter.FOV(1,2),...
       Crop_Parameter.FOV(2,1):Crop_Parameter.FOV(2,2),...
       Crop_Parameter.FOV(3,1):Crop_Parameter.FOV(3,2),:); 
end

end
