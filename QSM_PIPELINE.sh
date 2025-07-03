############################
#### Eleonora Patitucci #####
############################
#### 03.07.2025 ############
############################

# patituccie@cardiff.ac.uk




#### 1. Download DICOMs from XNAT
# Use the script 'DOWNLOADDATA_FROM_XNAT' (not shown here)




#### 2. Process DICOMs to generate magnitude/phase NIfTI files
# Define paths

wd=your_directory_with_dicoms                # Directory containing raw DICOM files
subjects=your_list_of_participants           # File with list of participant IDs (one per line)
output_dir=your_output_directory             # Where processed outputs will be stored

##### Tools required: ROMEO (for unwrapping), MEDI (for QSM reconstruction) 
# (Make sure the functions below know where these tools are installed and acquisition params)

#Processing 3T data: loop over subjects, process DICOMs with a MATLAB function
for ID in $(cat ${subjects});
	do matlab -nodesktop -r "process_3T_mag_phase_DICOM('${wd}/${ID}','${output_dir}/${ID}');";
done 

#Processing 7T data: loop over subjects, process separately magnitude + phase DICOMs

for ID in $(cat ${subjects});
	do matlab -nodesktop -r "process_ASPIRE_7T_mag_phase_DICOM('${mag_dicoms_folder}','${phase_dicoms_folder}','${output_dir}/${ID}');";
done




###3. Phase unwrapping using ROMEO
cd ${output_dir}                   # Go to output directory
chmod 777 romeo_cmd.txt             # Make ROMEO command file executable
./romeo_cmd.txt                     # Run ROMEO phase unwrapping (wrapped in this text file)




####4. Generate the local field (RDF = residual dipole field / local field map)

# 3T
matlab -nodesktop -r "local_field_gen_3T('${output_dir}')"
# 7T
matlab -nodesktop -r "local_field_gen_7T('${output_dir}')"




####5. Get susceptibility maps using NDI or MEDI algorithms
# 3T NDI
matlab -nodesktop -r "calc_qsm_ndi_3T('${output_dir}');"
# 3T MEDI
matlab -nodesktop -r "calc_qsm_medi_3T('${output_dir}');"
# 7T NDI
matlab -nodesktop -r "calc_qsm_ndi_7T('${output_dir}');"
# 7T MEDI
matlab -nodesktop -r "calc_qsm_medi_7T('${output_dir}');"




#6. Apply CSF shift correction to NDI QSM
matlab -nodesktop -r "csf_shift_internal_check('${output_dir}');quit()"  
# The above generates mean.txt containing mean CSF value

# Read the mean value and apply shift
mean=`cat ${output_dir}/mean.txt`;
fslmaths ${output_dir}/qsm_ndi_${ref}.nii -sub $mean  ${output_dir}/qsm_ndi_${ref}_csf_shift.nii; 
# This subtracts the mean CSF susceptibility to normalize QSM map

