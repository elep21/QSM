function calc_qsm_medi_3T(output_dir)


%Michael Germuska 
%Eleonora Patitucci

% tools location
run('~/matlab/MEDI_toolbox/MEDI_set_path.m');
addpath('~/matlab/mritools_Linux_3.3.5/matlab/NIfTI_20140122');
addpath('~/matlab/NDI_Toolbox');

%--------------------------------------------------------------------------
%% load data 
%--------------------------------------------------------------------------

% Input Files
fn_rdf = fullfile(output_dir, 'RDF.nii.gz');
fn_phase = fullfile(output_dir, 'unwrapped.nii.gz');
fn_complex = fullfile(output_dir, 'complex.nii.gz');
fn_mask = fullfile(output_dir, 'mask.nii.gz');
fn_mask_old = fullfile(output_dir, 'mask.nii.gz');

rdf_data=load_nii(fn_rdf);
phase_data=load_untouch_nii(fn_phase);
complex_data=load_nii(fn_complex);
mask_data=load_nii(fn_mask);
mask_old=load_nii(fn_mask_old);

RDF=rdf_data.img;
Mask=mask_data.img;
iFreq=phase_data.img;
iField=complex_data.img;
Mask_old=mask_old.img;

iMag = sqrt(sum(abs(iField).^2,4));


% load parameters for qsm calculation
fileID = fopen([output_dir '/DICOM_par.txt'],'r');
tline = fgetl(fileID);
voxel_size = str2num(fgetl(fileID));
tline = fgetl(fileID);
matrix_size = str2num(fgetl(fileID));
tline = fgetl(fileID);
CF = str2num(fgetl(fileID));
tline = fgetl(fileID);
delta_TE = str2num(fgetl(fileID));
tline = fgetl(fileID);
TE = str2num(fgetl(fileID));
tline = fgetl(fileID);
B0_dir = str2num(fgetl(fileID));




fclose(fileID);

%--------------------------------------------------------------------------
%% MEDI
%--------------------------------------------------------------------------

% Background field removal using Laplacian Boundary Value
% RDF = LBV(iFreq,Mask,matrix_size,voxel_size,0.005);
%%%% NMR Biomed 2014;27(3):312-319

% R2* map needed for ventricular CSF mask
R2s = arlo(TE, abs(iField));

% Ventricular CSF mask for zero referencing 
%	Requirement:
%		R2s:	R2* map

%figure, imshow(Mask(:,:,30))

Mask_CSF = extract_CSF(R2s, Mask_old, voxel_size);





[iFreq_raw, N_std] = Fit_ppm_complex(iField);

cd (output_dir)
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir Mask_CSF;

% Morphology enabled dipole inversion with zero reference using CSF (MEDI+0)
QSM = MEDI_L1('lambda',1000,'lambda_CSF',100,'merit','smv',5);

qsm_nii = make_nii(QSM);
qsm_nii.hdr.dime.pixdim(2:4) = voxel_size;
fn_qsm = fullfile(output_dir, 'qsm_medi.nii');
save_nii(qsm_nii, fn_qsm);

end
