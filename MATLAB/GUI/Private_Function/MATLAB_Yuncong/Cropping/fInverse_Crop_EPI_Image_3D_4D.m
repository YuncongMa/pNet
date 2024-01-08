function EPI_Image_3D_4D=fInverse_Crop_EPI_Image_3D_4D(EPI_Image_3D_4D,Crop_Parameter)
% By Yuncong Ma, Jan. 3, 2021
% An inverse function for fCrop_EPI_Image_3D_4D
% EPI_Image_3D_4D=fInverse_Crop_EPI_Image_3D_4D(EPI_Image_3D_4D,Crop_Parameter)
% Crop_Parameter contains FOV and FOV_Old

if ~isfield(Crop_Parameter,'FOV_Old') || ~isfield(Crop_Parameter,'FOV')
    fprintf('Error in Crop_Parameter: invalid settings of Crop_Parameter\n');
    return
end

Data=zeros([Crop_Parameter.FOV_Old(:,2)',size(EPI_Image_3D_4D,4)]);
Data(Crop_Parameter.FOV(1,1):Crop_Parameter.FOV(1,2),Crop_Parameter.FOV(2,1):Crop_Parameter.FOV(2,2),...
    Crop_Parameter.FOV(3,1):Crop_Parameter.FOV(3,2),:)=EPI_Image_3D_4D;
EPI_Image_3D_4D=Data;

end