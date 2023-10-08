function Image_Data=fAdd_Block_To_Image(Image_Data,Block)
%
% By Yuncong Ma, Dec. 12, 2018
% Image_Data=fAdd_Block_To_Image(Image_Data,Block)
% Image_Data is a 3D image matrix
% Block is a cell containing {[R,G,B],[x_min,x_max,y_min,y_max],...}
% Coordination system is switched in imshow

for i=1:length(Block)
    block2=Block{i};
    Color=reshape(block2{1},[1,1,3]);
    for j=2:length(block2)
        ps=block2{j};
        Image_Data(ps(1):ps(2),ps(3):ps(4),:)=repmat(Color,[ps(2)-ps(1)+1,ps(4)-ps(3)+1,1]);
    end
end

end