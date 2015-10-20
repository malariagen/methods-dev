import petl as etl
import subprocess
import os

etl.repr_index_header = True

plate_name = "W36000_QC417090_418525_417089_409377_20151001"
oxford_plate_name = "DK_KR_iPLEX_DK_20151020"
imeta_cmd = "imeta qu -z seq -d sequenom_plate = %s" % plate_name
irods_data_dir = "/nfs/team112_internal/rp7/recon/sanger_sequenom/irods_data/%s" % plate_name
combined_data_fn = "/nfs/team112_internal/rp7/recon/sanger_sequenom/combined_data/%s.txt" % plate_name
processed_data_fn = "/nfs/team112_internal/rp7/recon/sanger_sequenom/processed_data/%s_heights.xls" % oxford_plate_name
sample_mappings_fn = "/nfs/team112_internal/rp7/recon/sanger_sequenom/sample_mappings/Chris_Jacob_SQNM_various_Oct2015.xlsx"
sample_mappings_sheet_name = "sample_sets"
sample_set_code = "DK1092"
rextract_irods = False

subprocess.call("mkdir -p %s" % irods_data_dir, shell=True)

p = subprocess.Popen(imeta_cmd, shell=True, stdout=subprocess.PIPE)
imeta_out = p.stdout.read().decode("utf-8").split("\n")
for i in range(int(len(imeta_out)/3)):
    print(".", end="")
    irods_dir = imeta_out[i*3].replace("collection: ", "")
    irods_name = imeta_out[(i*3)+1].replace("dataObj: ", "")
    irods_path = "%s/%s" % (irods_dir, irods_name)
    local_path = "%s/%s" % (irods_data_dir, irods_name)
    if not os.path.exists(local_path) or rextract_irods:
        icp_cmd = "iget -f %s %s" % (irods_path, local_path)
        _ = subprocess.call(icp_cmd, shell=True)
    if i == 0:
        combined_file_cmd = "cat %s > %s" % (local_path, combined_data_fn)
    else:
        combined_file_cmd = "tail -n +2 %s >> %s" % (local_path, combined_data_fn)
    _ = subprocess.call(combined_file_cmd, shell=True)

tbl_sample_mappings = etl.fromxlsx(sample_mappings_fn, sample_mappings_sheet_name).selecteq('384_sample_set_code', sample_set_code)

tbl_processed = (etl.fromtsv(combined_data_fn)
    .convertnumbers()
    .update('PLATE', oxford_plate_name)
    .convert('ASSAY_ID', 'replace', 'W36000-', '')
    .leftjoin(tbl_sample_mappings, lkey='WELL_POSITION', rkey='384_position_code')
    .cutout('SAMPLE_ID')
    .rename('sample_code', 'SAMPLE_ID')
    .rename('STATUS', 'DESCRIPTION')
    .cut(['CUSTOMER', 'PROJECT', 'PLATE', 'EXPERIMENT', 'CHIP', 'WELL_POSITION',
    'ASSAY_ID', 'GENOTYPE_ID', 'DESCRIPTION', 'SAMPLE_ID', 'ALLELE', 'MASS',
    'HEIGHT'])
    .sort(['PLATE', 'WELL_POSITION', 'ASSAY_ID', 'MASS'])
)

tbl_processed.totsv(processed_data_fn)
