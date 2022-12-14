function [regParam, projParam] = DefaultProcessingParams(dataSet)
% Registration parameters
regParam.refScan = []; % which scans (in concatenated units) should be used to generate the reference image
regParam.refRun = 1; % which run should be used to generate the reference image (only used if refScan is empty)
regParam.chunkSize = 100; % register this many scans at once, to avoid loading too much data
regParam.histmatch = false; % true; % enable if registration struggles due to bleaching
regParam.avgT = 0; % 15;  temporal Gaussian filter width (not downsampling). scans, not seconds
regParam.avgTsigma = 0; % 5
regParam.binXY = 2; % spatial downsampling factor - bigger number -> faster
regParam.binT = 1; % temporal downsampling - keep at 1
regParam.highpass = 0; % use highpass filtering to clean data during registration
regParam.lowpass = 0; % use lowpass filtering to clean data during registration
regParam.medFilter = [0,0,0]; % dimensions (pixels/scans) to use for median filtering filtering to clean data during registration. set to zeros to skip med filtering
regParam.edges = [60,60,40,40]; % cut this many pixels off from the [left, right, top, bottom] to avoid edges contaminating the registration
regParam.name = ''; % names can be used to distinguish different versions of registration
% Projection parameters
projParam.type = 'max';
switch lower(dataSet)
    case 'afferents'
        regParam.refChan = 'green'; 
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = false;
    case 'astrocyte'
        regParam.refChan = 'red'; 
        regParam.method = 'translation'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'macrophage'
        regParam.refChan = 'red';
        regParam.method = 'translation'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'pollen'
        regParam.refChan = 'green';
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 1; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = true;
    case 'vasculature'
        regParam.refChan = 'green';
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 2; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    otherwise
        regParam.refChan = 'green';
        regParam.method = 'translation'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 1; % Hz
        projParam.color = {'red','green'}; % 
        projParam.vol = false;
end
projParam.umPerPixel_target = 1; % target um/pix for projections. spatial downsampling 
projParam.overwrite = false;
end