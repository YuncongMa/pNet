function Setting=fLoad_Data_Setting(Read_Write,File,Load_Data)
% Yuncong Ma, 2/20/2023
% To Read mat or write txt files about settings in Load_Data for NMF APP
% Setting=fLoad_Data_Setting(Read_Write,File,Load_Data)



switch Read_Write
    case 'Read'
        Setting=fLoad_MATLAB_Single_Variable(File);

    case 'Write'
        FID=fopen(File,'w');
        fprintf(FID,'[Data_Type] = %s\n',Load_Data.Data_Type);
        fprintf(FID,'[Data_Format] = %s\n',Load_Data.Data_Format);
        fclose(FID);
    otherwise
        error('Unknown settings for Read or Write : %s\n',Read_Write);
end

end