
# setup all sub-folders in the pNet result folder
dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = pNet.setup_result_folder(dir_pnet_result)

# ============== Setup ============== #
# ============== Data Input
# setup dataInput
pNet.setup_scan_info(
    dir_pnet_dataInput=dir_pnet_dataInput,
    dataType=dataType, dataFormat=dataFormat,
    file_scan=file_scan, file_subject_ID=file_subject_ID,
    file_subject_folder=file_subject_folder, file_group_ID=file_group_ID,
    Combine_Scan=Combine_Scan
)
# setup brain template
# Volume and surface data types require different inputs to compute the brain template
if file_Brain_Template is None:
    if dataType == 'Volume':
        pNet.setup_brain_template(
            dir_pnet_dataInput=dir_pnet_dataInput,
            dataType=dataType, dataFormat=dataFormat,
            templateFormat=templateFormat,
            file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
            maskValue=maskValue
        )
    elif dataType == 'Surface':
        pNet.setup_brain_template(
            dir_pnet_dataInput=dir_pnet_dataInput,
            dataType=dataType, dataFormat=dataFormat,
            templateFormat=templateFormat,
            file_surfL=file_surfL, file_surfR=file_surfR,
            file_maskL=file_maskL, file_maskR=file_maskR,
            maskValue=maskValue,
            file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated
        )
    elif dataType == 'Surface-Volume':
        pNet.setup_brain_template(
            dir_pnet_dataInput=dir_pnet_dataInput,
            dataType=dataType, dataFormat=dataFormat,
            templateFormat=templateFormat,
            file_surfL=file_surfL, file_surfR=file_surfR,
            file_maskL=file_maskL, file_maskR=file_maskR,
            file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
            maskValue=maskValue,
            file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated
        )

else:
    pNet.setup_brain_template(dir_pnet_dataInput, file_Brain_Template)

# ============== FN Computation
pNet.setup_NMF_setting(
    dir_pnet_result=dir_pnet_result,
    K=K,
    Combine_Scan=Combine_Scan,
    file_gFN=file_gFN,
    samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS,
    maxIter=maxIter, minIter=minIter, meanFitRatio=meanFitRatio, error=error, normW=normW,
    Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL,
    vxI=vxI, ard=ard, eta=eta,
    nRepeat=nRepeat,
    Computation_Mode=Computation_Mode,
    dataPrecision=dataPrecision,
    outputFormat=outputFormat
)

# =============== Visualization
pNet.setup_Visualization(
    dir_pnet_result=dir_pnet_result,
    synchronized_view=synchronized_view,
    synchronized_colorbar=synchronized_colorbar
)

# =============== Server
pNet.setup_server(
    dir_pnet=dir_pnet,
    dir_pnet_result=dir_pnet_result,
    dir_python=dir_python,
    submit_command=submit_command,
    thread_command=thread_command,
    memory_command=memory_command,
    log_command=log_command,
    computation_resource=computation_resource
)

print('All setups are finished\n', flush=True)

# ============== Run ============== #
print('Start to run\n', flush=True)
# FN Computation
pNet.run_FN_Computation_torch_server(dir_pnet_result)

# Quality Control
pNet.run_quality_control_torch_server(dir_pnet_result)

# Visualization
pNet.run_Visualization_server(dir_pnet_result)

print('All run s are finished\n', flush=True)
