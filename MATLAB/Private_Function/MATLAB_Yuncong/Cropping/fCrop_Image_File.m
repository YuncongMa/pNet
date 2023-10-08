function fCrop_Image_File(File_Path,Crop_Parameter)
% By Yuncong Ma, Sep. 9, 2018
% fCrop_Image_File(File_Path,Crop_Parameter)
% Crop a saved image file
% Crop_Parameter is [width_min,height_min,width_max,height_max]
% All parameters should withing [0,1]

Data=imread(File_Path);

Size=size(Data);

ps1=round(Size(1:2).*Crop_Parameter([2,1]));
ps2=round(Size(1:2).*Crop_Parameter([4,3]));

ps1(ps1<1)=1;
ps1(ps1>Size(1:2))=Size(ps1>Size(1:2));

ps2(ps2<1)=1;
ps2(ps2>Size(1:2))=Size(ps2>Size(1:2));

Data=Data(ps1(1):ps2(1),ps1(2):ps2(2),:);

imwrite(Data,File_Path);


end