function Setting=fCompute_FN_Setting(Read_Write,File,Compute_FN)
% Yuncong Ma, 3/17/2023
% To Read mat or write txt files about settings in Compute_FN for NMF APP
% Setting=fCompute_FN_Setting(Read_Write,File,Compute_FN)



switch Read_Write
    case 'Read'
        Setting=fLoad_MATLAB_Single_Variable(File);
        
    case 'Write'
        FID=fopen(File,'w');
        fprintf(FID,'[Method] = %s\n',Compute_FN.Method);
        fprintf(FID,'[K] = %d\n',Compute_FN.K); 
        fprintf(FID,'[GroupFN.Compute_Flag] = %d\n',Compute_FN.GroupFN.Compute_Flag);
        fprintf(FID,'[GroupFN.BootStrap.File_Selection] = %d\n',Compute_FN.GroupFN.BootStrap.File_Selection);
        fprintf(FID,'[GroupFN.BootStrap.Repetition] = %d\n',Compute_FN.GroupFN.BootStrap.Repetition);
        fprintf(FID,'[GroupFN.spaR] = %d\n',Compute_FN.GroupFN.spaR);
        fprintf(FID,'[GroupFN.vxI] = %d\n',Compute_FN.GroupFN.vxI);
        fprintf(FID,'[GroupFN.iterNum] = %d\n',Compute_FN.GroupFN.iterNum);
        fprintf(FID,'[GroupFN.Alpha] = %d\n',Compute_FN.GroupFN.Alpha);
        fprintf(FID,'[GroupFN.Beta] = %d\n',Compute_FN.GroupFN.Beta);
        fprintf(FID,'[PersonalizedFN.Combine_Flag] = %d\n',Compute_FN.PersonalizedFN.Combine_Flag);
        fprintf(FID,'[PersonalizedFN.spaR] = %d\n',Compute_FN.PersonalizedFN.spaR);
        fprintf(FID,'[PersonalizedFN.vxI] = %d\n',Compute_FN.PersonalizedFN.vxI);
        fprintf(FID,'[PersonalizedFN.ard] = %d\n',Compute_FN.PersonalizedFN.ard);
        fprintf(FID,'[PersonalizedFN.eta] = %d\n',Compute_FN.PersonalizedFN.eta);
        fprintf(FID,'[PersonalizedFN.iterNum] = %d\n',Compute_FN.PersonalizedFN.iterNum);
        fprintf(FID,'[PersonalizedFN.alphaS21] = %d\n',Compute_FN.PersonalizedFN.alphaS21);
        fprintf(FID,'[PersonalizedFN.alphaL] = %d\n',Compute_FN.PersonalizedFN.alphaL);
        fprintf(FID,'[Parallel.Flag] = %d\n',Compute_FN.Parallel.Flag);
        fprintf(FID,'[Parallel.N_Thread] = %d\n',Compute_FN.Parallel.N_Thread);
        fclose(FID);
    otherwise
        error('Unknown settings for Read or Write : %s\n',Read_Write);
end

end