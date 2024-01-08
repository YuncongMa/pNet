function Setting=fFN_Computation_Setting(Read_Write,File,FN_Computation)
% Yuncong Ma, 6/19/2023
% To Read mat or write txt files about settings in FN_Computation for NMF APP
% Setting=fFN_Computation_Setting(Read_Write,File,FN_Computation)



switch Read_Write
    case 'Read'
        Setting=fLoad_MATLAB_Single_Variable(File);
        
    case 'Write'
        FID=fopen(File,'w');
        fprintf(FID,'[Method] = %s\n',FN_Computation.Method);
        fprintf(FID,'[K] = %d\n',FN_Computation.K); 
        fprintf(FID,'[GroupFN.Compute_Flag] = %d\n',FN_Computation.GroupFN.Compute_Flag);
        fprintf(FID,'[GroupFN.BootStrap.File_Selection] = %d\n',FN_Computation.GroupFN.BootStrap.File_Selection);
        fprintf(FID,'[GroupFN.BootStrap.Repetition] = %d\n',FN_Computation.GroupFN.BootStrap.Repetition);
        fprintf(FID,'[GroupFN.spaR] = %d\n',FN_Computation.GroupFN.spaR);
        fprintf(FID,'[GroupFN.vxI] = %d\n',FN_Computation.GroupFN.vxI);
        fprintf(FID,'[GroupFN.iterNum] = %d\n',FN_Computation.GroupFN.iterNum);
        fprintf(FID,'[GroupFN.Alpha] = %d\n',FN_Computation.GroupFN.Alpha);
        fprintf(FID,'[GroupFN.Beta] = %d\n',FN_Computation.GroupFN.Beta);
        fprintf(FID,'[PersonalizedFN.Combine_Flag] = %d\n',FN_Computation.PersonalizedFN.Combine_Flag);
        fprintf(FID,'[PersonalizedFN.spaR] = %d\n',FN_Computation.PersonalizedFN.spaR);
        fprintf(FID,'[PersonalizedFN.vxI] = %d\n',FN_Computation.PersonalizedFN.vxI);
        fprintf(FID,'[PersonalizedFN.ard] = %d\n',FN_Computation.PersonalizedFN.ard);
        fprintf(FID,'[PersonalizedFN.eta] = %d\n',FN_Computation.PersonalizedFN.eta);
        fprintf(FID,'[PersonalizedFN.iterNum] = %d\n',FN_Computation.PersonalizedFN.iterNum);
        fprintf(FID,'[PersonalizedFN.alphaS21] = %d\n',FN_Computation.PersonalizedFN.alphaS21);
        fprintf(FID,'[PersonalizedFN.alphaL] = %d\n',FN_Computation.PersonalizedFN.alphaL);
        fprintf(FID,'[Parallel.Flag] = %d\n',FN_Computation.Parallel.Flag);
        fprintf(FID,'[Parallel.N_Thread] = %d\n',FN_Computation.Parallel.N_Thread);
        fclose(FID);
    otherwise
        error('Unknown settings for Read or Write : %s\n',Read_Write);
end

end