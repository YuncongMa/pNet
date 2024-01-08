function Setting=fData_Input_Setting(Read_Write,File,Data_Input)
% Yuncong Ma, 6/19/2023
% To Read mat or write txt files about settings in Data_Input for NMF APP
% Setting=fData_Input_Setting(Read_Write,File,Data_Input)



switch Read_Write
    case 'Read'
        Setting=fLoad_MATLAB_Single_Variable(File);

    case 'Write'
        FID=fopen(File,'w');
        fprintf(FID,'[Data_Type] = %s\n',Data_Input.Data_Type);
        fprintf(FID,'[Data_Format] = %s\n',Data_Input.Data_Format);
        fclose(FID);
    otherwise
        error('Unknown settings for Read or Write : %s\n',Read_Write);
end

end