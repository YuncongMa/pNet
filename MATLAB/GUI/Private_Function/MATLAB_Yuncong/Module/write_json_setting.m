function write_json_setting(Setting_Structure, File_Json)
% Yuncong Ma, 2/1/204
% Mimic write_json_setting in pNet python version
% write_json_setting(Setting_Structure, File_Json)


Setting_Json=jsonencode(Setting_Structure,PrettyPrint=true);

fid = fopen(File_Json, 'w');
fprintf(fid, '%s', Setting_Json);
fclose(fid);
end