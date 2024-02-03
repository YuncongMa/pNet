# Yuncong Ma, 12/18/2023
# Make a web page based report for fast visual examination

from Module.Visualization import *


dir_python = os.path.dirname(os.path.abspath(__file__))


def run_web_report(dir_pnet_result: str):
    """
    generate HTML based web report

    :param dir_pnet_result:
    :return:

    Yuncong Ma, 12/20/2023
    """

    # get directories of sub-folders
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, _ = setup_result_folder(dir_pnet_result)

    # log file
    logFile = os.path.join(dir_pnet_result, 'Log_Report.log')
    logFile = open(logFile, 'w')
    print_log('\nStart generating HTML based web report at ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())) + '\n',
              logFile=logFile, stop=False)

    # load settings for data input and FN computation
    if not os.path.isfile(os.path.join(dir_pnet_dataInput, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder Data_Input')
    if not os.path.isfile(os.path.join(dir_pnet_FNC, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder FN_Computation')
    settingDataInput = load_json_setting(os.path.join(dir_pnet_dataInput, 'Setting.json'))
    settingFNC = load_json_setting(os.path.join(dir_pnet_FNC, 'Setting.json'))
    setting = {'Data_Input': settingDataInput, 'FN_Computation': settingFNC}
    print_log('Settings are loaded from folder Data_Input and FN_Computation', logFile=logFile, stop=False)

    # load basic settings
    dataType = setting['Data_Input']['Data_Type']
    dataFormat = setting['Data_Input']['Data_Format']

    # info about the fMRI dataset
    file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
    file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
    file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')

    list_scan = load_txt_list(file_scan)
    nScan = len(list_scan)
    list_subject_ID = load_txt_list(file_subject_ID)
    nSubject = len(np.unique(list_subject_ID))
    list_subject, subject_index = np.unique(load_txt_list(file_subject_ID), return_index=True)
    list_subject_ID_unqiue = list_subject_ID[subject_index]
    list_subject_folder = load_txt_list(file_subject_folder)
    list_subject_folder_unique, folder_index = np.unique(list_subject_folder, return_index=True)
    nFolder = len(list_subject_folder_unique)

    # template for web page
    template_summary = os.path.join(dir_python, 'Web_Template_Summary.html')

    # Generate the summary web page
    file_summary = os.path.join(dir_pnet_result, 'Summary.html')
    pnet_FN_method = setting['FN_Computation']['Method']
    K = setting['FN_Computation']['K']

    with open(template_summary, 'r') as file:
        html_as_string = file.read()
    # report title
    report_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    html_as_string = html_as_string.replace('{$report_time$}', str(report_time))
    # setting
    html_as_string = html_as_string.replace('{$pnet_FN_method$}', str(pnet_FN_method))
    html_as_string = html_as_string.replace('{$K$}', str(K))
    html_as_string = html_as_string.replace('{$dataType$}', str(dataType))
    html_as_string = html_as_string.replace('{$dataFormat$}', str(dataFormat))
    html_as_string = html_as_string.replace('{$nScan$}', str(nScan))
    html_as_string = html_as_string.replace('{$nSubject$}', str(nSubject))
    # Compress images if not done yet
    image_size = (2000, 10000)
    if not os.path.isfile(os.path.join(dir_pnet_gFN, 'All(Compressed).jpg')):
        compress_image(os.path.join(dir_pnet_gFN, 'All(Compressed).jpg'),
                       os.path.join(dir_pnet_gFN, 'All.jpg'),
                       image_size=image_size)
    for i in range(nFolder):
        if os.path.isfile(os.path.join(dir_pnet_pFN, list_subject_folder_unique[i], 'All.jpg')) and not os.path.isfile(os.path.join(dir_pnet_pFN, list_subject_folder_unique[i], 'All(Compressed).jpg')):
            compress_image(os.path.join(dir_pnet_pFN, list_subject_folder_unique[i], 'All(Compressed).jpg'),
                           os.path.join(dir_pnet_pFN, list_subject_folder_unique[i], 'All.jpg'),
                           image_size=image_size)

    # gFN
    if setting['FN_Computation']['Group_FN']['file_gFN'] is None:
        text_gFN = 'The group FNs are derived using the whole fMRI dataset'
    else:
        text_gFN = 'The group FNs are loaded from precomputed results at ' + setting['FN_Computation']['Group_FN']['file_gFN']
    html_as_string = html_as_string.replace('{$text_gFN$}', str(text_gFN))
    figure_gFN = './Group_FN/All(Compressed).jpg'
    html_as_string = html_as_string.replace('{$figure_gFN$}', str(figure_gFN))
    # pFN examples
    nMax = 10
    frames = [Image.open(os.path.join(dir_pnet_pFN, list_subject_folder[subject_index[i]], 'All(Compressed).jpg')) for i in range(np.minimum(nMax, nSubject))]
    frame_one = frames[0]
    frame_one.save(os.path.join(dir_pnet_pFN, 'Example.gif'), format="GIF", append_images=frames, save_all=True, duration=1000, loop=0, optimize=False)
    pFN_example = './' + os.path.join('Personalized_FN', 'Example.gif')
    html_as_string = html_as_string.replace('{$pFN_example$}', str(pFN_example))
    # pFN links
    link_pFN = ''
    pre_sub = list_subject_ID_unqiue[0]
    for i in range(nFolder):
        if list_subject_ID[folder_index[i]] != pre_sub:
            link_pFN = link_pFN + "<br />"
            pre_sub = list_subject_ID[folder_index[i]]
        file_pFN_indv = './' + os.path.join('Personalized_FN', list_subject_folder_unique[i], 'All(Compressed).jpg')
        link_pFN = link_pFN + f" <a href='{file_pFN_indv}' target='_blank' title='{list_subject_folder_unique[i]}'>({list_subject_folder_unique[i]})</a>\n"
    html_as_string = html_as_string.replace('{$link_pFN$}', str(link_pFN))
    # QC
    Result = load_matlab_single_variable(os.path.join(dir_pnet_QC, 'Result.mat'))
    n_missmatch = np.sum(np.min(Result['Delta_Spatial_Correspondence'][0, 0], axis=1) < 0)
    ps_missmatch = np.where(np.min(Result['Delta_Spatial_Correspondence'][0, 0], axis=1) < 0)
    if n_missmatch == 0:
        text_qc = 'pFNs of all scans passed QC. <br />'
    else:
        text_qc = 'There are ' + str(n_missmatch) + ' scans failed the QC, meaning that they have at least one pFN miss matched. <br />'
        text_qc = ''
        for i in range(n_missmatch):
            text_qc = text_qc + list_subject_folder_unique[ps_missmatch[i]] + ' <br />'
    html_as_string = html_as_string.replace('{$text_qc$}', str(text_qc))

    file_summary = open(file_summary, 'w')
    print(html_as_string, file=file_summary)






