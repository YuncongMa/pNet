function [FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,FN)
% Yuncong Ma, 2/22/2023
% Reformat FN for both surface and volume formats
% [FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,FN)
% Reformat FN for output when it is a 2D matrix
% Reformat FN for input when it is a file directory

Flag=0;
Message='';

[Setting.Load_Data,Flag,Message] =fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Setting.mat'));
if Flag
    return
end

if ischar(FN)
    [FN,Flag,Message] =fLoad_MATLAB_Single_Variable(FN);
    if Flag
        FN=[];
        return
    end
    switch Setting.Load_Data.Data_Type
        case 'Surface'
            return

        case 'Volume'
            [Brain_Mask,Flag,Message] =fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Brain_Mask.mat'));
            if Flag
                FN=[];
                return
            end
            if ~isequal(size(Brain_Mask),size(FN(:,:,:,1)))
                FN=[];
                Flag=1;
                Message='The brain mask has a different size compared to the FN';
                return
            end
            FN=fApply_Mask(Brain_Mask>0,FN,-1);
        otherwise
            FN=[];
            Flag=1;
            Message=['Unsupported data type: ',Setting.Load_Data.Data_Type];
            return
    end

elseif isnumeric(FN)

    if length(size(FN))~=2
        FN=[];
        Flag=1;
        Message='FN needs to be in 2D matrix when in the NMF APP computation';
        return
    end

    switch Setting.Load_Data.Data_Type
        case 'Surface'
            return

        case 'Volume'
            [Brain_Mask,Flag,Message] =fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Brain_Mask.mat'));
            if Flag
                FN=[];
                return
            end
            if sum(Brain_Mask(:)>0)~=size(FN,1)
                FN=[];
                Flag=1;
                Message='The brain mask has different numbers of voxels compared to the FN';
                return
            end
            FN=fInverse_Mask(Brain_Mask>0,FN,-1);
        otherwise
            FN=[];
            Flag=1;
            Message=['Unsupported data type: ',Setting.Load_Data.Data_Type];
            return
    end
end

end