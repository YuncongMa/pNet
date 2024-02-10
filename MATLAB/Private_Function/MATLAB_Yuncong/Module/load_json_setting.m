function Setting_Structure=load_json_setting(File_Json)
% Yuncong Ma, 2/1/204
% Mimic load_json_setting in pNet python version
% Setting_Structure=load_json_setting(File_Json)


fid = fopen(File_Json); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
Setting_Structure = jsondecode(str); % Using the jsondecode function to parse JSON from string

end