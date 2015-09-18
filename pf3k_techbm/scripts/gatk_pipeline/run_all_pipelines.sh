#!/bin/bash
export PIPELINES_FN='~/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/pipeline_parameters.txt'
export PIPELINES_FN='/Users/rpearson/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/pipeline_parameters.txt'
line=$(head -n 1 $PIPELINES_FN)
IFS='\t' read -a var_names <<< "$line"
echo "${var_names[3]}"

IFS=$'\r\n' :; vars=(<$PIPELINES_FN)

IFS=$'\r\n' read -d '' -r -a lines=<$PIPELINES_FN
${lines[0]}

readarray lines < $PIPELINES_FN

IFS=$'\r\n' read -d '' -r -a lines < $PIPELINES_FN