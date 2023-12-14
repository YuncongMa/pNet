# Yuncong Ma, 12/14/2023
# Make a web page based report for fast visual examination

import os
from Data_Input import *

dir_python = os.path.dirname(os.path.abspath(__file__))


def run_web_report(dir_pnet_result: str):
    """
    generate HTML based web report

    :param dir_pnet_result:
    :return:

    Yuncong Ma, 12/14/2023
    """

    # get directories of sub-folders
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, _, _ = setup_result_folder(dir_pnet_result)

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
    print('Settings are loaded from folder Data_Input and FN_Computation', file=logFile_FNC, flush=True)

    # load basic settings
    dataType = setting['Data_Input']['Data_Type']
    dataFormat = setting['Data_Input']['Data_Format']

    file_template = os.path.join(dir_python, 'Report_Template.html')
