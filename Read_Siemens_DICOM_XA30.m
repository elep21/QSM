function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,files]=Read_Siemens_DICOM_XA30(DicomFolder)

filelist = dir(DicomFolder);
i=1;
while i<=length(filelist)
    if filelist(i).isdir==1
        filelist = filelist([1:i-1 i+1:end]);   % eliminate folders
    else
        i=i+1;
    end
end

fnTemp = [DicomFolder '/' filelist(1).name];

info = dicominfo(fnTemp);
matrix_size(1) = single(info.Width);
matrix_size(2) = single(info.Height);

matrix_size(3) = info.NumberOfFrames;

matrix_size(4) = 1;


NumEcho = single(info.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1.EchoTrainLength);
voxel_size(1,1) = single(info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1));
voxel_size(2,1) = single(info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2));
voxel_size(3,1) = single(info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness);

CF = info.SharedFunctionalGroupsSequence.Item_1.MRImagingModifierSequence.Item_1.TransmitterFrequency *1e6;
isz=[matrix_size(1) matrix_size(2) matrix_size(3) matrix_size(4)*NumEcho];
iPhase = single(zeros(isz));
iMag = single(zeros(isz));
filesM={};
filesP={};

TE = single(zeros([NumEcho 1]));
p_ImagePositionPatient = zeros(3, matrix_size(4)*NumEcho);
m_ImagePositionPatient = zeros(3, matrix_size(4)*NumEcho);
r_EchoNumber = zeros(matrix_size(4)*NumEcho,1);
minSlice = 1e10;
maxSlice = -1e10;

% for (deidentified?) VIDA dicoms there is no longer an EchoTime field
% detecting here
filename=[DicomFolder '/' filelist(1).name];
info2 = struct();
info2.EchoNumbers = info.PerFrameFunctionalGroupsSequence.Item_1.FrameContentSequence.Item_1.DimensionIndexValues(2);

%info2 = get_dcm_tags(filename, {'EchoNumbers'});

progress='';

if ~isfield(info2, 'EchoNumber')
    TE=[];
    % need to loop through all images to find the 
    % EchoTime->EchoNumber mapping
    for i = 1:length(filelist)
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Collecting echo times %d', i);
        fprintf(progress);
        filename=[DicomFolder '/' filelist(i).name];
        
        fnTemp = [DicomFolder '/' filelist(i).name];

        info = dicominfo(fnTemp);
        info2.EchoNumbers = info.PerFrameFunctionalGroupsSequence.Item_1.FrameContentSequence.Item_1.DimensionIndexValues(2);

        if isempty(find(TE==info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime * 1e-3))
            TE = sort([TE; info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime*1e-3]);
        end
    end
    fprintf('\n')
    if numel(TE) ~= NumEcho
        error(['Number of echo times found (' num2str(numel(TE)) ...
               ') is different from EchoTrainLength (' num2str(NumEcho) ').']);
    end
end

rctr = 0; 
ictr=0;
progress='';
for i = 1:length(filelist)
    filename=[DicomFolder '/' filelist(i).name];
    %info2 = get_dcm_tags(filename);
    fnTemp = [DicomFolder '/' filelist(i).name];

        
    info = dicominfo(fnTemp);
    if info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient(3)<minSlice
        minSlice = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient(3);
        minLoc = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
    end
    if info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient(3)>maxSlice
        maxSlice = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient(3);
        maxLoc = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
    end
    if isfield(info2, 'EchoNumber')
        if info2.EchoNumbers>NumEcho
            TE = [TE; zeros([info2.EchoNumbers - NumEcho 1])];
            NumEcho = info2.EchoNumbers;
        end
        if TE(info2.EchoNumbers)==0
            TE(info2.EchoNumbers)=info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime*1e-3;
        end
    else
        info2.EchoNumbers=find(TE==info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime*1e-3);
    end
    if (info.PerFrameFunctionalGroupsSequence.Item_1.MRImageFrameTypeSequence.Item_1.FrameType(18)=='P')||(info.PerFrameFunctionalGroupsSequence.Item_1.MRImageFrameTypeSequence.Item_1.FrameType(18)=='p')
        rctr = rctr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        p_ImagePositionPatient(:,rctr) = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
        r_EchoNumber(rctr) = info2.EchoNumbers;
        ph = pagetranspose(single(dicomread(filename)));
        iPhase(:,:,:,rctr)  = (ph*info.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleSlope + info.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleIntercept)/single(max(ph(:)))*pi;%phase
        filesP{rctr} = filename;
    elseif (info.PerFrameFunctionalGroupsSequence.Item_1.MRImageFrameTypeSequence.Item_1.FrameType(18)=='M')||(info.PerFrameFunctionalGroupsSequence.Item_1.MRImageFrameTypeSequence.Item_1.FrameType(18)=='m')
        ictr = ictr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        m_ImagePositionPatient(:,ictr) = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
        i_EchoNumber(ictr) = info2.EchoNumbers;
        iMag(:,:,:,ictr)  = pagetranspose(single(dicomread(filename)));%magnitude
        filesM{ictr} = filename;
    end
end
fprintf('\n');

matrix_size(4) = round(norm(maxLoc - minLoc)/voxel_size(3))+1 ;

Affine2D = reshape(info.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient,[3 2]);
Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
%B0_dir = Affine3D\[0 0 1]';

iop = info.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient;
row = iop(1:3);
col = iop(4:6);

B0_dir = cross(row, col);
B0_dir = round(B0_dir / norm(B0_dir));

files.Affine3D=Affine3D;
files.minLoc=minLoc;
files.maxLoc=maxLoc;
n_exp=matrix_size(4)*NumEcho;

p_sz=size(p_ImagePositionPatient); 
p_sz(1)=1;
p_minLoc=repmat(minLoc, p_sz);
if n_exp ~= p_sz(2)
    warning(['Number of phase images (' ... 
        num2str(p_sz(2)) ...
        ') is different from the expected number (' ...
        num2str(n_exp) ...
        ').'])
end

p_slice = int32(round(sqrt(sum((p_ImagePositionPatient-p_minLoc).^2,1))/voxel_size(3)) +1);
p_ind = sub2ind([matrix_size(4) NumEcho], p_slice(:), int32(r_EchoNumber(:)));
iPhase(:,:,:,p_ind) = iPhase;
files.P=cell(matrix_size(4), NumEcho);
files.P(p_ind)=filesP;
if n_exp ~= size(iPhase,4)
    iPhase = padarray(iPhase, double([0 0  0 n_exp-size(iPhase,4)]), 'post');
    warning(['Some phase images are missing']); 
end

m_sz=size(m_ImagePositionPatient); m_sz(1)=1;
if n_exp ~= m_sz(2)
    warning(['Number of magnitude images (' ... 
        num2str(m_sz(2)) ...
        ') is different from the expected number (' ...
        num2str(n_exp) ...
        ').'])
end
if m_sz(2) ~= p_sz(2)
    warning(['Number of phase images (' ... 
        num2str(p_sz(2)) ...
        ') is different from number of magnitude images (' ...
        num2str(m_sz(2)) ...
        ').'])
end
m_minLoc=repmat(minLoc, m_sz);

m_slice = int32(round(sqrt(sum((m_ImagePositionPatient-m_minLoc).^2,1))/voxel_size(3)) +1);
m_ind = sub2ind([matrix_size(4) NumEcho], m_slice(:), int32(i_EchoNumber(:)));
iMag(:,:,:,m_ind) = iMag;
files.M=cell(matrix_size(4), NumEcho);
files.M(m_ind)=filesM;
if n_exp ~= size(iMag,4)
    iMag = padarray(iMag, double([0 0 0 n_exp-size(iMag,4)]), 'post');
    warning(['Some magnitdue images are missing']); 
end

iField = reshape(iMag.*exp(-1i*iPhase), ...
    [matrix_size(1) matrix_size(2) matrix_size(3) NumEcho]);
% iField = permute(iField,[2 1 3 4 5]); %This is because the first dimension is row in DICOM but COLUMN in MATLAB
% iField(:,:,1:2:end,:) = -iField(:,:,1:2:end,:);
if length(TE)==1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end

end

function info = get_dcm_tags(filename, tags)
if nargin<2
    tags={'SliceLocation','ImagePositionPatient',...
        'EchoTime','EchoNumbers','ImageType','RescaleSlope','RescaleIntercept'};
end
attrs=dicomattrs(filename);
info=struct;
for t=tags
    t=char(t);
    [gr, el] = dicomlookup(t);
    if ~strcmp(t,'ImageType')
        for i=1:length(attrs);
            if (attrs(i).Group==gr)&&(attrs(i).Element==el)
                eval(['info.' t '= sscanf(char(attrs(i).Data), ''%f\\'');'])
                break;
            end
        end
    else
        for i=1:length(attrs);
            if (attrs(i).Group==gr)&&(attrs(i).Element==el)
                info.ImageType = char(attrs(i).Data);
            end
        end
    end
end
end