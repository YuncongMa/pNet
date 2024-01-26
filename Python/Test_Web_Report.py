import pNet
import os


dir_pnet_result = os.path.join(pNet.dir_pNet, 'Test/Test_FN17_Server_2')

pNet.visualize_quality_control(dir_pnet_result)

# pNet.run_web_report(dir_pnet_result)

