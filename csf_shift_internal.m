function QSM_0ref=csf_shift_internal(output_dir)
	%Stefano Zappalà	
	%Michael Germuska


	% tools location
    run('~/matlab/MEDI_toolbox/MEDI_set_path.m');
    addpath('~/matlab/mritools_Linux_3.3.5/matlab/NIfTI_20140122');
    addpath('~/matlab/NDI_Toolbox');

    %--------------------------------------------------------------------------
    %% load data 
    %--------------------------------------------------------------------------

    % Input Files
    fn_qsm = fullfile(output_dir, 'qsm_ndi.nii.gz');
    fn_complex = fullfile(output_dir, 'complex.nii.gz');
    fn_mask = fullfile(output_dir, 'mask.nii.gz');
    qsm_data=load_nii(fn_qsm);
    complex_data=load_nii(fn_complex);
    mask_data=load_nii(fn_mask);
    QSM=qsm_data.img;
    Mask=mask_data.img;
    iField=complex_data.img;


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


    % R2* map needed for ventricular CSF mask
    % TE  array containing te values (in s)
    R2s = arlo(TE, abs(iField));

    % Ventricular CSF mask for zero referencing 
    %	Requirement:
    %		R2s:	R2* map
    Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
    

    Mask_CSF_nii = make_nii(Mask_CSF);
    Mask_CSF_nii.hdr.dime.pixdim(2:4) = voxel_size;
    fn_rdf = fullfile(output_dir, ['Mask_CSF_HR.nii.gz']);
    save_nii(Mask_CSF_nii, fn_rdf);

    title_str=split(output_dir,'/');
    figure;isosurface(Mask_CSF);title(title_str(end-1))
    
    % Zero reference using CSF
    %mean in CSF voxel 
    QSM_0ref = mean(QSM(Mask_CSF==1),'all');
    
    fileID = fopen([output_dir '/mean.txt'],'w');
    fprintf(fileID, '%.3f', QSM_0ref);
    fclose(fileID);
    
    drawnow
end
