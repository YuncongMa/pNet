function Image_RGB=fColorize(Matrix,Color_Function)
% By Yuncong Ma, Mar. 1, 2018
% Image_RGB=fColorize(Matrix,Color_Function)
% Color_Function=[0,0,0,0;.25,1,0,0;.5,1,1,0];
% for each row [value, R, G, B]
% value should be from the smallest to the largest

if size(Color_Function,2)~=4 || size(Color_Function,1)<2 || ~isequal(Color_Function(:,1),sort(Color_Function(:,1)))
    fprintf('Error in fColorize: wrong settings in Color_Function.\n');
    return;
end

Dim=size(Matrix);
Matrix=reshape(Matrix,[numel(Matrix),1]);
Image_RGB=zeros(numel(Matrix),3);
ps=Matrix<Color_Function(1,1);
Image_RGB(ps,:)=repmat(Color_Function(1,2:4),[sum(ps),1]);
ps=Matrix>=Color_Function(end,1);
Image_RGB(ps,:)=repmat(Color_Function(end,2:4),[sum(ps),1]);

for p=2:size(Color_Function,1)
    ps=Color_Function(p-1,1)<=Matrix & Matrix<Color_Function(p,1);
    RGB=((Matrix-Color_Function(p-1,1)).*ps*Color_Function(p,2:4)+...
        (Color_Function(p,1)-Matrix).*ps*Color_Function(p-1,2:4))/(Color_Function(p,1)-Color_Function(p-1,1));
    Image_RGB(ps,:)=RGB(ps,:);
end

Image_RGB=reshape(Image_RGB,[Dim,3]);

end