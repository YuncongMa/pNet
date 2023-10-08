function [Unique_Cell_String,Count]=fUnique_Cell_String(Cell_String)
% Yuncong Ma, 10/31/2022
% [Unique_Cell_String,Count]=fUnique_Cell_String(Cell_String)
% Get unique cell string from a cell string, and count the occurance of
% each string
% eg. Cell_String={'a','b','a','c'}, Unique_Cell_String={'a','b','c'},
% Count=[2,1,1]

Unique_Cell_String=Cell_String(1);
Count=1;
for i=2:length(Cell_String)
    ps=find(strcmp(Unique_Cell_String,Cell_String{i}));
    if isempty(ps)
        Unique_Cell_String{end+1}=Cell_String{i};
        Count(end+1)=1;
    else
        Count(ps)=Count(ps)+1;
    end
end

end