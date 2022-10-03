import InSAR_Pipeline as isp
import os

file_prefix = "infiles"
output_prefix = "outfiles"
scene = "kumamoto_2016"
IW =  "IW1"

firstBurstIndex = 1
lastBurstIndex = 9
ML_nRgLooks = 6

output_dir = scene + "_" + IW + "_" + "wsl"
output_path = os.path.join(output_prefix, output_dir)

os.makedirs(output_path, exist_ok=True)

# Nepal 2015 nepal_2015
if scene == "nepal_2015":
    filename_1 = os.path.join(file_prefix,'S1A_IW_SLC__1SSV_20150417T001852_20150417T001922_005516_0070C1_460B.zip')
    filename_2 = os.path.join(file_prefix,'S1A_IW_SLC__1SDV_20150429T001907_20150429T001935_005691_0074DC_7332.zip')

# Rockridge 2019 rockridge_2019
elif scene == "rockridge_2019":
    filename_1 = os.path.join(file_prefix,'S1B_IW_SLC__1SDV_20190626T020647_20190626T020714_016861_01FB9D_6003.zip')
    filename_2 = os.path.join(file_prefix,'S1B_IW_SLC__1SDV_20190708T020648_20190708T020715_017036_0200CC_D043.zip')

# Kumamoto Earthquake 2016 kumamoto_2016
elif scene == "kumamoto_2016":
    filename_1 = os.path.join(file_prefix,'S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip')
    filename_2 = os.path.join(file_prefix,'S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip')

else:
    raise RuntimeError("Specify a valid scene")


# Define input/output files
out_filename_i =os.path.join(output_path, 'InSAR_pipeline_I')

in_filename_ii =os.path.join(output_path, 'InSAR_pipeline_I.dim')
out_filename_ii =os.path.join(output_path, 'InSAR_pipeline_II')

in_filename_iii = os.path.join(output_path, 'InSAR_pipeline_II.dim')
out_filename_iii =os.path.join(output_path, 'InSAR_pipeline_III')

in_filename_iv = os.path.join(output_path, 'InSAR_pipeline_II.dim')
out_filename_iv =os.path.join(output_path, 'InSAR_pipeline_IV_snaphu_test')


# Call the processing pipeline step by step
# isp.InSAR_pipeline_I(filename_1, filename_2, IW, firstBurstIndex, lastBurstIndex, out_filename_i)

#isp.InSAR_pipeline_II(in_filename_ii, ML_nRgLooks, out_filename_ii)

#isp.InSAR_pipeline_III(in_filename_iii, out_filename_iii)

isp.InSAR_pipeline_IV(in_filename_iv, out_filename_iv)