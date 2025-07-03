function local_field_gen_3T(output_dir)

%Michael Germuska 
%Eleonora Patitucci

% tools location
run('~/matlab/MEDI_toolbox/MEDI_set_path.m');
addpath('~/matlab/mritools_Linux_3.3.5/matlab/NIfTI_20140122');

% Input Files
fn_field = fullfile(output_dir, 'complex.nii.gz');
fn_phase = fullfile(output_dir, 'unwrapped.nii.gz');

field_data=load_untouch_nii(fn_field);
phase_data=load_untouch_nii(fn_phase);

iField=field_data.img;
iFreq=phase_data.img;

% load parameters for background field removal
fileID = fopen([output_dir '/DICOM_par.txt'],'r');
tline = fgetl(fileID);
voxel_size = str2num(fgetl(fileID));
tline = fgetl(fileID);
matrix_size = str2num(fgetl(fileID));
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);
B0_dir = str2num(fgetl(fileID));

fclose(fileID);

% compute N_std .... need to see if there is a better way of doing this
disp('estimate noise');
[iFreq_raw, N_std] = Fit_ppm_complex(iField);

iMag = sqrt(sum(abs(iField).^2,4));
% Use FSL BET to extract brain mask
% But first, check if any manual adjustment of the masks has been done
if exist(fullfile(output_dir, 'mask.nii.gz'), 'file') ~= 2
    disp('create mask');
    Mask = BET(iMag,matrix_size,voxel_size, 0.25); % set fractional threshold to 0.3 (try 0.25 ... 0.3 is even too harsh sometimes!)
else
    disp('use manually adjusted mask');
    mask_data=load_untouch_nii(fullfile(output_dir, 'mask.nii.gz'));
    Mask=double(mask_data.img);
end

% Background field removal using Projection onto Dipole Fields
disp('background field removal');
%tol = str2num('0.01');
%n_CG = str2num('1000');

%RDF = PDF(iFreq, N_std, Mask, matrix_size,voxel_size, B0_dir, tol, n_CG);
RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);

disp('save RDF');mask_data

RDF_nii = make_nii(RDF);
RDF_nii.hdr.dime.pixdim(2:4) = voxel_size;
fn_rdf = fullfile(output_dir, 'RDF.nii.gz');
save_nii(RDF_nii, fn_rdf);

disp('save Mask');

mask_nii = make_nii(Mask);
mask_nii.hdr.dime.pixdim(2:4) = voxel_size;
fn_mask = fullfile(output_dir, 'mask.nii.gz');
save_nii(mask_nii, fn_mask);

exit

end
