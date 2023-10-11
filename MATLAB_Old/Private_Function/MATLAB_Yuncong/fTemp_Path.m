function Path_Temp=fTemp_Path(varargin)
% By Yuncong Ma, Aug. 29, 2018
% Path_Temp=fTemp_Path()
% Path_Temp=fTemp_Path(Path_Unzip)
% Create a temporary path for storing data


if nargin==0
    Path_Temp=tempname;
    fMake_Folder(Path_Temp);
elseif nargin==1
    [~,folder]=fileparts(tempname);
    Path_Temp=fullfile(varargin{1},folder);
    fMake_Folder(Path_Temp);
else
    error('Error in fTemp_Path: wrong inputs for options');
end

end

