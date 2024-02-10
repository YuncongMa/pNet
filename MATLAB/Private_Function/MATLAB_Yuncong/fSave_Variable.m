function fSave_Variable(Directory_Variable,Variable)
% By Yuncong Ma, Aug. 23, 2018
% 
% fSave_Variable(Directory_Variable,Variable)
% Variable is a cell containing names and variables
% eg. {{'Center',Center},{'Index',Index}}
% all variables will be saved into .mat files

fMake_Folder(Directory_Variable);
N_Variable=length(Variable);

for i=1:N_Variable
    if ~iscell(Variable{i}) || length(Variable{i})~=2
        error('Error in fSave_Variable: wrong format of input for Variable.');
    elseif ~ischar(Variable{i}{1})
        error('Error in fSave_Variable: the %d-th variable has no string name.',i);
    elseif sum(strcmp(Variable{i}{1},{'i','Directory_Variable','N_Variable'}))
        error('Error in fSave_Variable: cannot use variable names including i, Directory_Variable, Variable and N_Variable.');
    end
end

for i=1:N_Variable
    eval([Variable{i}{1},'=Variable{',num2str(i),'}{2};']);
    save(fullfile(Directory_Variable,[Variable{i}{1},'.mat']),Variable{i}{1});
end

end



