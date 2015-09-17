import os
import petl as etl
import subprocess

PIPELINES_FN='/Users/rpearson/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/pipeline_parameters.xlsx'
SCRIPT_FN=''

tbl_pipelines = etl.fromxlsx(PIPELINES_FN)
var_names = tbl_pipelines.header()
for rec in tbl_pipelines.data():
    for i, var_name in enumerate(var_names):
        print(var_name, rec[i])
        os.environ[var_name] = rec[i]
    subprocess.call(SCRIPT_FN)
l