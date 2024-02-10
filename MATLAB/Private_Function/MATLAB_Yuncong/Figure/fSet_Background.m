function fSet_Background(Figure_Size,Color)
% By Yuncong Ma, Sep. 8, 2018
% fSet_Background(Figure_Size,Color)
% Designed for MATLAB running in macOS
% eg. fSet_Background([900,100],[0,0,0]);

axes('Position',[0,0,1,1]);
Color_3D=zeros(1,1,3);
Color_3D(1,1,:)=Color;
Image=repmat(Color_3D,[Figure_Size([2,1]),1]);
imshow(Image)
ax=gca;
axis(ax,'fill');



