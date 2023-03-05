function Setting=fCompute_FN_Setting(Read_Write,File,Compute_FN)
% Yuncong Ma, 2/21/2023
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
        fprintf(FID,'[IndividualizedFN.Combine_Flag] = %d\n',Compute_FN.IndividualizedFN.Combine_Flag);
        fprintf(FID,'[IndividualizedFN.spaR] = %d\n',Compute_FN.IndividualizedFN.spaR);
        fprintf(FID,'[IndividualizedFN.vxI] = %d\n',Compute_FN.IndividualizedFN.vxI);
        fprintf(FID,'[IndividualizedFN.ard] = %d\n',Compute_FN.IndividualizedFN.ard);
        fprintf(FID,'[IndividualizedFN.eta] = %d\n',Compute_FN.IndividualizedFN.eta);
        fprintf(FID,'[IndividualizedFN.iterNum] = %d\n',Compute_FN.IndividualizedFN.iterNum);
        fprintf(FID,'[IndividualizedFN.alphaS21] = %d\n',Compute_FN.IndividualizedFN.alphaS21);
        fprintf(FID,'[IndividualizedFN.alphaL] = %d\n',Compute_FN.IndividualizedFN.alphaL);
        fprintf(FID,'[Parallel.Flag] = %d\n',Compute_FN.Parallel.Flag);
        fprintf(FID,'[Parallel.N_Thread] = %d\n',Compute_FN.Parallel.N_Thread);
        fclose(FID);
    otherwise
        error('Unknown settings for Read or Write : %s\n',Read_Write);
end

end