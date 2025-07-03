%written by Michael Germuska
%for SIEMENS DICOM images consisting of magnitude and phase parts
%   output
%   iField - the multi-echo complex MR image
%   voxel_size - size of the voxel in mm
%   matrix_size - dimension of the field of view
%   CF - central frequency in Hz
%   delta_TE - TE spacing in s
%   TE - echo times in s
%   B0_dir - direction of the B0 field


function [iMag,iPhase,iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_7T_MG(path_mag,path_ph)

addpath('natsortfiles_test/');

disp('--> get magnitude DICOM list ...');

% read in DICOMs of both magnitude and raw unfiltered phase images
% read in magnitude DICOMs
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));

% germuska@cubric sort files so that lack of leading zeros is taken care of
mag_list = natsortfiles(mag_list);


% get the sequence parameters
dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);

CF = dicom_info.ImagingFrequency *1e6;
voxel_size(1,1) = single(dicom_info.PixelSpacing(1));
voxel_size(2,1) = single(dicom_info.PixelSpacing(2));
voxel_size(3,1) = single(dicom_info.SliceThickness);

matrix_size(1) = single(dicom_info.Width);
matrix_size(2) = single(dicom_info.Height);

disp('finding min and max slice locations...');

minSlice = 1e10;
maxSlice = -1e10;
for i = 1:length(mag_list)
    dicom_info = dicominfo([path_mag,filesep,mag_list(i).name]);
        if dicom_info.SliceLocation<minSlice
            minSlice = dicom_info.SliceLocation;
            minLoc = dicom_info.ImagePositionPatient;
        end
        if dicom_info.SliceLocation>maxSlice
            maxSlice = dicom_info.SliceLocation;
            maxLoc = dicom_info.ImagePositionPatient;
        end            
end
matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;    
Affine2D = reshape(dicom_info.ImageOrientationPatient,[3 2]);
Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
B0_dir = Affine3D\[0 0 1]';

calc_EchoTrainLength = length(mag_list) / matrix_size(3);

dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);
EchoTrainLength = dicom_info.EchoTrainLength;

if calc_EchoTrainLength ~= EchoTrainLength
    disp(['Warning! Apparent number of echoes is ' int2str(calc_EchoTrainLength) ' number of echoes in DICOM header is ' int2str(EchoTrainLength)]);
    disp(['Setting EchoTrainLength to ', int2str(calc_EchoTrainLength)])
    EchoTrainLength = calc_EchoTrainLength;
end


for i = 1:EchoTrainLength % read in TEs    
    dicom_info = dicominfo([path_mag,filesep,mag_list(1+(i-1)*(length(mag_list))./EchoTrainLength).name]);
   % TE(dicom_info.EchoNumbers) = dicom_info.EchoTime*1e-3; 5 matlab before
   % 2021a does not read dicom info correctly for 7T qsm data... fudge by
   % using 'i' as echo number (should work as dicoms sorted before reading in)
    TE(i) = dicom_info.EchoTime*1e-3;
    
end

if length(TE)==1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end

disp('--> read magnitude DICOMS ...');

for i = 1:length(mag_list)
    [NS,NE] = ind2sub([length(mag_list)./EchoTrainLength,EchoTrainLength],i);
    iMag(:,:,NS,NE) = permute(single(dicomread([path_mag,filesep,mag_list(i).name])),[2,1]); 
       
end

disp('--> get phase DICOM list ...');

% read in phase DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% germuska@cubric sort files so that lack of leading zeros is taken care of
ph_list = natsortfiles(ph_list);

disp('--> read phase DICOMS ...');

for i = 1:length(ph_list)
    [NS,NE] = ind2sub([length(ph_list)./EchoTrainLength,EchoTrainLength],i);
    iPhase(:,:,NS,NE) = permute(single(dicomread([path_ph,filesep,ph_list(i).name])),[2,1]);     % covert to [-pi pi] range
    iPhase(:,:,NS,NE) = iPhase(:,:,NS,NE)/4095*2*pi - pi;
    
end

size_mag = size(iMag)
size_Phase = size(iPhase)

iField = iMag.*exp(-1i*iPhase);

end
