function calc_qsm_ndi_3T(output_dir)

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
fn_complex = fullfile(output_dir, 'complex.nii.gz');
fn_mask = fullfile(output_dir, 'mask.nii.gz');
rdf_data=load_nii(fn_rdf);
complex_data=load_nii(fn_complex);
mask_data=load_nii(fn_mask);
RDF=rdf_data.img;
Mask=mask_data.img;
iField=complex_data.img;

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
tline = fgetl(fileID);
tline = fgetl(fileID);
B0_dir = str2num(fgetl(fileID));

fclose(fileID);


N = size(Mask);
kernel=dipole_kernel(matrix_size, voxel_size, B0_dir);


%--------------------------------------------------------------------------
%% TKD
%--------------------------------------------------------------------------

%kthre = 0.15;       % truncation threshold for TKD recon

%krnl = kernel;
%kernel_inv = zeros(N);
%kernel_inv( abs(krnl) > kthre ) = 1 ./ krnl(abs(krnl) > kthre);

%tic
%    chi_tkd = real( ifftn( fftn(RDF(:,:,:,1)) .* kernel_inv ) ) .* Mask; 
%toc

%convert Chi to ppm
%QSM = chi_tkd/(2*pi*delta_TE*CF)*1e6.*Mask;

% save qsm results to nifti
%qsm_nii = make_nii(QSM);
%qsm_nii.hdr.dime.pixdim(2:4) = voxel_size;
%fn_qsm = fullfile(output_dir, 'qsm_tkd.nii.gz');
%save_nii(qsm_nii, fn_qsm);



%--------------------------------------------------------------------------
%% NDI
%--------------------------------------------------------------------------

step_size = 1;      % gradient descent step size
num_iter = 1000;  %1000  

nd = 1;             % number of head directions to use in the NDI recon, can be set between 1 and 5

phs_use = RDF(:,:,:,1:nd);

%scale magnitude data between 0 and 1
mgn_use=iMag.*Mask;
mgn_use=mgn_use./max(mgn_use(:));

M2 = repmat(mean(mgn_use(:,:,:,1:nd),4).^2, [1,1,1,nd]);        % magnitude weighting


if nd == 1
    % 1-direction may require more iterations: decrease convergence tolerance
    tol = 0.25; % 0.15 checking impact 0.25
else
    tol = 1;
end

Chi = zeross(N);
grad_prev = 0;

tic
for t = 1:num_iter
    temp = M2 .* sin(ifft(ifft(ifft(kernel(:,:,:,1:nd) .* repmat(fftn(Chi),[1,1,1,nd]), [], 1), [], 2), [], 3) - phs_use);

    grad_f = 2 * sum(ifft(ifft(ifft(kernel(:,:,:,1:nd) .* fft(fft(fft(temp, [], 1), [], 2), [], 3), [], 1), [], 2), [], 3), 4);

    Chi = Chi - step_size * real(grad_f);

    update_grad = rmse(grad_prev, grad_f);

    disp(['iter: ', num2str(t), '   grad update:', num2str(update_grad)])

    if update_grad < tol
        break
    end

    grad_prev = grad_f;
end
toc

%convert Chi to ppm
QSM = Chi/(2*pi*delta_TE*CF)*1e6.*Mask;

% save qsm results to nifti
qsm_nii = make_nii(QSM);
qsm_nii.hdr.dime.pixdim(2:4) = voxel_size;
fn_qsm = fullfile(output_dir, 'qsm_ndi.nii.gz');
save_nii(qsm_nii, fn_qsm);

exit

end
