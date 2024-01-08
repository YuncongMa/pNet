function Options=fOption(Function_Name,Default_Options,Option)
% By Yuncong Ma, Mar. 18, 2018
% Options=fOption(Function_Name,Default_Options,Option)
% To read cell or structure type settings into structure like settings
% Default_Options is the default settings in structure type
% Option should be {{'Setting1',Parameters}} or a {Structure}
% Options will be Options.Setting1=Parameters
% if Option is {{'Setting1',P1,P2}}, Options will be
% Options.Setting1={P1,P2};

Options=Default_Options;
if length(Option)==1 && isstruct(Option{1})
    Option_Structure=Option{1};
    Field_Name=fieldnames(Option_Structure);
    for j=1:length(Field_Name)
        if ~isfield(Options,Field_Name{j})
            error('Error in %s: unknown option setting for %s\n',Function_Name,Field_Name{j});
        end
        Options.(Field_Name{j})=Option_Structure.(Field_Name{j});
    end
    return
end 

Option_Cell=Option;
for i=1:length(Option_Cell)
    if length(Option_Cell{i})<2
        error('Error in %s: wrong option settings\n',Function_Name);
    end
    if ~isfield(Options,Option_Cell{i}{1})
        error('Error in %s: unknown option setting for %s\n',Function_Name,Option_Cell{i}{1});
    end
    
    if length(Option_Cell{i})==2
        Options.(Option_Cell{i}{1})=Option_Cell{i}{2};
    else
        Options.(Option_Cell{i}{1})=Option_Cell{i}(2:end);   
    end
end

end