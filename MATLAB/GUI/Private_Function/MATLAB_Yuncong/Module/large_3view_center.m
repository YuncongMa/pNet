function Center=large_3view_center(Map)
% Yuncong Ma, 2/1/2024
% Find the 3D center of a 3D map
% mimic large_3view_center in pNet python
% Center=large_3view_center(Map)

Map(Map<0)=0;

Center=zeros(1,3);

temp=squeeze(sum(sum(Map,2),1));
[~,ps]=max(temp);
Center(3)=ps(1);

temp=squeeze(sum(sum(Map,3),1));
[~,ps]=max(temp);
Center(2)=ps(1);

temp=squeeze(sum(sum(Map,3),2));
[~,ps]=max(temp);
Center(1)=ps(1);


end