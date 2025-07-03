function process_3T_mag_phase_DICOM(DICOM_folder, output_folder)

%Michael Germuska and Eleonora Patitucci - March 2023
% create nifti images from DICOM data and write out DICOM parameters to a
% text file for further processing


%update tools location 
romeo_path ='~/matlab/mritools_Linux_3.3.5/bin/';
addpath('~/matlab/mritools_Linux_3.3.5/matlab/NIfTI_20140122');
run('~/matlab/MEDI_toolbox/MEDI_set_path.m');


%modified to run on matlab 2015
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_DICOM(DICOM_folder);

% Estimate the frequency offset in each of the voxel using a complex
% fitting (even echo spacing)
[iFreq_raw, N_std] = Fit_ppm_complex(iField);


% Compute magnitude image
iMag = sqrt(sum(abs(iField).^2,4));

 
% write out magnitude, phase, and complex nifti files
fn_phase = [output_folder '/phase.nii'];
fn_mag = [output_folder '/mag.nii'];
fn_complex = [output_folder '/complex.nii'];
fn_unwrapped = [output_folder '/unwrapped.nii'];


disp('writing nifti data');

phase_nii = make_nii(iFreq_raw);
phase_nii.hdr.dime.pixdim(2:4) = voxel_size;
phase_nii.img = flip(phase_nii.img,2)
save_nii(phase_nii, fn_phase);

mag_nii = make_nii(iMag);
mag_nii.hdr.dime.pixdim(2:4) = voxel_size;
mag_nii.img = flip(mag_nii.img,2)
save_nii(mag_nii, fn_mag);

complex_nii = make_nii(iField);
complex_nii.hdr.dime.pixdim(2:4) = voxel_size;
complex_nii.img = flip(complex_nii.img,2)
save_nii(complex_nii, fn_complex);


%create string for romeo command 

%check number of echoes 
TE_str=[num2str(1000*TE(1)), ' ',   num2str(1000*TE(2)) , ' ', num2str(1000*TE(3)), ' ' , num2str(1000*TE(4)), ' ' , num2str(1000*TE(5)), ' ' , num2str(1000*TE(6)), ' ' , num2str(1000*TE(7)) ];
romeo_cmd= [romeo_path, 'romeo -p ', fn_phase , ' -o ' , fn_unwrapped, ' -m ', fn_mag, ' -t [', TE_str , '] -v'];

romeo_fileID = fopen([output_folder '/romeo_cmd.txt'],'w');
fprintf(romeo_fileID, '%s', romeo_cmd);
fprintf(romeo_fileID, ['\ngzip ' output_folder '/*.nii -f']);
fclose(romeo_fileID);

disp('write DICOM parameters to text file');

fileID = fopen([output_folder '/DICOM_par.txt'],'w');
fprintf(fileID, '%s\n', 'voxel size');
fprintf(fileID, '%.3f %s %.3f %s %.3f', voxel_size(1), ' ', voxel_size(2) ,' ', voxel_size(3));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'matrix size');
fprintf(fileID,  '%i %s %i %s %i', matrix_size(1), ' ', matrix_size(2) ,' ', matrix_size(3));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'centre frequency');
fprintf(fileID,  '%i\n', CF);

fprintf(fileID, '%s\n', 'delta TE');
fprintf(fileID, '%.4f', delta_TE);
fprintf(fileID, '\n');

%check number of echoes 
fprintf(fileID, '%s\n', 'TE');
fprintf(fileID, '%.4f %s %.4f %s %.4f %s %.4f %s %.4f %s %.4f %s %.4f %s', TE(1), ' ', TE(2) ,' ', TE(3) ,' ', TE(4),' ', TE(5),' ', TE(6) ,' ', TE(7));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'B0 dir');
fprintf(fileID, '%.4f %s %.4f %s %.4f', B0_dir(1), ' ', B0_dir(2) ,' ', B0_dir(3));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'matrix size');
fprintf(fileID,  '%i %s %i %s %i', matrix_size(1), ' ', matrix_size(2) ,' ', matrix_size(3));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'centre frequency');
fprintf(fileID,  '%i\n', CF);

fprintf(fileID, '%s\n', 'delta TE');
fprintf(fileID, '%.4f', delta_TE);

fprintf(fileID, '\n');
%check number of echoes 
fprintf(fileID, '%s\n', 'TE');
fprintf(fileID, '%.4f %s %.4f %s %.4f %s %.4f %s %.4f %s %.4f %s %.4f %s', TE(1), ' ', TE(2) ,' ', TE(3) ,' ', TE(4),' ', TE(5),' ', TE(6) ,' ', TE(7));
fprintf(fileID, '\n');

fprintf(fileID, '%s\n', 'B0 dir');
fprintf(fileID, '%.4f %s %.4f %s %.4f', B0_dir(1), ' ', B0_dir(2) ,' ', B0_dir(3));
fprintf(fileID, '\n');

fclose(fileID);

exit


end
