function Center=fMass_Center(Image_2D_3D)
% By Yuncong Ma, Mar. 1, 2018
% Center=fMass_Center(Image_2D_3D)
% To find the center of the largest connected region
% Image_2D_3D will be binarized by 0

% find largest connected component
Image_2D_3D=Image_2D_3D>0;
if length(size(Image_2D_3D))<2
    fprintf('Error in fMass_Center: wrong size of Image_2D_3D\n');
    Center=[];
    return;
end
CC=bwconncomp(Image_2D_3D);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,ps] = max(numPixels);
Image_2D_3D=0*Image_2D_3D;
Image_2D_3D(CC.PixelIdxList{ps})=1;


Nz=size(Image_2D_3D,3);
Nx=size(Image_2D_3D,1);
Ny=size(Image_2D_3D,2);
Center_List=zeros([Nz,3]);
for z=1:Nz
    Center_X=(1:Nx)*sum(Image_2D_3D(:,:,z),2);
    Center_Y=(1:Ny)*sum(Image_2D_3D(:,:,z),1)';
    Center_List(z,:)=[Center_X,Center_Y,sum(sum(Image_2D_3D(:,:,z)))];
end

Center(1)=sum(Center_List(:,1))/sum(Center_List(:,3));
Center(2)=sum(Center_List(:,2))/sum(Center_List(:,3));

if length(size(Image_2D_3D))==3
   Center(3)=(1:Nz)*Center_List(:,3)/sum(Center_List(:,3)); 
end

end