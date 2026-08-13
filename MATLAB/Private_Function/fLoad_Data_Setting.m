function fLoad_Data_Setting(action, filePath, data)
% FLOAD_DATA_SETTING loads or saves data settings to a specified file.
%
%   ACTION   : 'Read' or 'Write'
%   FILEPATH : path to the setting file
%   DATA     : data structure to save or load

    if strcmpi(action, 'Write')
        save(filePath, 'data');
    elseif strcmpi(action, 'Read')
        if exist(filePath, 'file')
            load(filePath, 'data');
        else
            data = struct();
        end
    end
end
