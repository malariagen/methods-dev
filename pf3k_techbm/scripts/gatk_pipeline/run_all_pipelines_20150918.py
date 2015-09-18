# To run (ideally using screen):
# python $HOME/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/run_all_pipelines_20150918.py >> $HOME/src/github/malariagen/methods-dev/pf3k_techbm/log/gatk_pipeline_20150918.log
import os
import petl as etl
import subprocess

PIPELINES_FN='%s/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/pipeline_parameters_20150918.xlsx' % os.environ['HOME']
SCRIPT_FN='%s/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/gatk_pipeline_20150918.sh' % os.environ['HOME']

tbl_pipelines = etl.fromxlsx(PIPELINES_FN)
var_names = tbl_pipelines.header()
for rec in tbl_pipelines.data():
    print(rec)
    for i, var_name in enumerate(var_names):
        print(var_name, rec[i])
        os.environ[var_name] = str(rec[i])
    subprocess.call(SCRIPT_FN)
