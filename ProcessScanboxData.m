clear; clc; close all;
dataDir = 'D:\2photon\'; % 'D:\2photon\Simone\'; %'C:\2photon';
dataSet = 'Afferents'; % 'Pollen'; % 'Vasculature'; %  'Astrocyte'; %  'Anatomy'; %      'Neutrophil_Simone'; %  'NGC'; % 'Neutrophil'; % 
% Parse data table
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx'; %'D:\MATLAB\ImagingDatasets.xlsx'; % 'D:\MATLAB\NGCdata.xlsx';  Simone
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), ...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), 'ref',find(contains(colNames, 'Ref')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt);
% Registration parameters
regParam.refScan = []; % which scans (in concatenated units) should be used to generate the reference image
regParam.refRun = NaN; % which run should be used to generate the reference image (only used if refScan is empty)
regParam.chunkSize = 100; % register this many scans at once, to avoid loading too much data
regParam.histmatch = false; % true; % enable if registration struggles due to bleaching
regParam.avgT = 0; % 15;  temporal Gaussian filter width (not downsampling). scans, not seconds
regParam.avgTsigma = 0; % 5
regParam.binXY = 2; % spatial downsampling factor - bigger number -> faster
regParam.binT = 1; % temporal downsampling - keep at 1
regParam.highpass = 0; % use highpass filtering to clean data during registration
regParam.lowpass = 0; % use lowpass filtering to clean data during registration
regParam.medFilter = [0,0,0]; % dimensions (pixels/scans) to use for median filtering filtering to clean data during registration. set to zeros to skip med filtering
regParam.minInt = 1500; % used for edge selection
regParam.edges = [60,60,40,40]; % cut this many pixels off from the [left, right, top, bottom] to avoid edges contaminating the registration
regParam.name = ''; % names can be used to distinguish different versions of registration
% Projection parameters
projParam.type = 'max';
switch lower(dataSet)
    case 'vasculature'
        regParam.refChan = 'green';
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 2; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'astrocyte'
        regParam.refChan = 'red'; 
        regParam.method = 'translation'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'afferents'
        regParam.refChan = 'green'; 
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = false;
    case 'pollen'
        regParam.refChan = 'green';
        regParam.method = 'affine'; % use affine to get deformation, rigid just to align data
        projParam.rate_target = 1; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = true;
end
projParam.umPerPixel_target = 1; % target um/pix for projections. spatial downsampling 
projParam.overwrite = false;

% Various useful subsets of the data to choose from
%xDone = find([dataTable{:,dataCol.done}] > 0);
%x3D = intersect( xDone, find([dataTable{:,dataCol.volume}] == 1 ));
%x2D = intersect( xDone, find( [dataTable{:,dataCol.volume}] == 0 ));
% Choose which subset to  process
xPresent = 110; %44; %[18,22,24,30:32]; % flip(100:102); %45:47; % [66:69]; %6;  62,64,
Npresent = numel(xPresent);
overwrite = false;
for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Parse data table
    [expt{x}, runInfo{x}, regParam, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir, regParam, projParam);
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x}); 

    % Get locomotion data
    for runs = expt{x}.runs % flip(expt{x}.runs) % 
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',true ); % plot(loco{x}(regParam.refRun).Vdown)
    end
    % Use the longest period of stillness to define the reference image, or define it by hand
    if ~isempty(loco{x}(1).quad)
        if any(cellfun(@isempty, {loco{x}.stateDown})) %isempty(loco{x}.stateDown)
            loco{x} = GetLocoState(expt{x}, loco{x}, 'dir',strcat(dataDir, expt{x}.mouse,'\'), 'name',expt{x}.mouse, 'var','velocity', 'show',false); %
        end
        % Determine reference run and scans (longest pre-CSD epoch of stillness)
        [~,tformPath]= FileFinder(expt{x}.dir, 'contains','regTform');
        if isempty(tformPath)
            [regParam.refRun, regParam.refScan] = DetermineReference(expt{x}, Tscan{x}, loco{x}, 1:2); % 1:4
        else
            a = load(tformPath{1}, 'params');
            regParam = a.params; clearvars a;
        end
    else
        regParam.refRun = 1;  % plot(loco{x}(regParam.refRun).Vdown)
        regParam.refScan = 130:260; %expt{x}.scanLim
    end

    % Register individual runs, starting with the reference run
    runInfo{x}(regParam.refRun) = RegisterRun( runInfo{x}(regParam.refRun), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid', 'interp',false);
    for runs = flip(setdiff(1:expt{x}.Nruns, regParam.refRun)) % expt{x}.runs % 
        runInfo{x}(runs) = RegisterRun( runInfo{x}(runs), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid', 'interp',false); %'rigid' , 'edges',[80,80,40,20]
    end
    
    % Concatenate runs and metadata
    catInfo{x} = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat', 'overwrite',true); % Get concatenated metadata
    interRunShift = ConcatenateExptRuns(expt{x}, runInfo{x}, catInfo{x}, 'refRun',regParam.refRun, 'refChan',regParam.refChan, 'setEdge',regParam.edges, 'sbx','sbxdft'); %1  expt{x}.refChan
    catProj = WriteSbxProjection(expt{x}.sbx.cat, catInfo{x}, 'chan','green', 'type','cat', 'overwrite',overwrite, 'monochrome',true, 'RGB',true); % 
    % Write projections of concatenated data
    projParam.edge = regParam.edges; %[40,150,40,40]; %[60,60,40,20]; % [40,40,40,80];%segParams{x}.edges; %[60,60,20,20]; % [70 40 20 20];
    projParam.z = {4:7, 8:10}; %{17:22, 24:27, 18, 19, 20, 21, 22, 23}; %{3:17, 18:28, 18, 19, 20, 21, 22, 23}; %{17:22, 23:30}; %{29:56, 5:25}; % {7:10, 18:20, 27:30};  % 1; %
    projParam.overwrite = false;
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); % write projections of unregistered data by run

    %
    CatInterpZ(catInfo{x}.path, catInfo{x}, 'refScan',regParam.refScan, 'edges',regParam.edges );
    MakeSbxZ_new(catInfo{x}.path, catInfo{x}); % , shiftPath
    
    % Register the concatenated data (see RegisterCat3D, AlignPlanes and RegisterSBX for more info)   
    fprintf('\n   Performing planar registration... '); % (reference averaged over scans %i - %i)
    %regParam.edges = GetEdges3D( catProj, varargin )
    regParam.refScan = regParam.refScan + expt{x}.scanLims(regParam.refRun);
    %repairStruct = struct('z',1:expt{x}.Nplane, 'scan',5001:expt{x}.totScan);
    regParam = AlignPlanes( expt{x}.sbx.cat, catInfo{x}, regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg ); %, 'repair',repairStruct
    [~,deform{x}] = GetDeformCat3D(catInfo{x}, 'show',true, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last'));  
    regProj = WriteSbxProjection(expt{x}.sbx.reg, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','reg', 'overwrite',overwrite);
 
    % Generate downsampled, possibly z-projected, movies for each run from the concatenated data
    %[projParam.edge, x_result, y_result] = GetEdges3D( regProj(:,:,:,2), 'show',true );
    % {
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); %   
    %}
    % {
    if strcmpi(dataSet, 'afferents')
        % segement the calcium imaging data
        try
            segParams = GetSegParams( catInfo{x} );
        catch
            zProj{1} = projParam.z; % {4:10, 12:20, 22:28};
            segParams = struct('name','', 'z',zProj , 'edges',projParam.edge, 'seg_scan',[], 'cens_scans',1914:2099, 'min_foot',50, 'blur',2, 'corr_thresh_pct',75, 'chunk_size',30, 'hp_cutoff',51, 'min_vol',100, 'xyproj_width',100 );
        end
        segParams = SegmentSbx(catInfo{x}, segParams, 'overwrite',false); 
           
        % Generate mean/max projection images
        maxProjPath = [expt{x}.dir, expt{x}.name, '_maxProj.tif']; %[expt{x}.dir, expt{x}.name, '_affineProj.tif'];
        meanProjPath = [expt{x}.dir, expt{x}.name, '_meanProj.tif'];
        if ~exist(maxProjPath, 'file')
            [~,segProjPath] = FileFinder(expt{x}.dir, 'contains','segProj', 'type','tif'); %segProjPath = segProjPath{1};
            if ~isempty(segProjPath)
                for z = flip(1:numel(segProjPath))
                    tempSegMov{z} = imread_big(segProjPath{z});
                end
                tempSegMov = cat(4, tempSegMov{:});
                expt{x}.maxProj = max(tempSegMov, [] , [3,4]);
                saveastiff(expt{x}.maxProj, maxProjPath);
                expt{x}.meanProj = mean(tempSegMov, [3,4]);
                saveastiff(uint16(expt{x}.meanProj), meanProjPath);
                clearvars tempSegMov;
            else
                fprintf('\nseg projections do not exist - cannot generate mean/max projections!');
            end
        else
            expt{x}.maxProj = loadtiff( maxProjPath );
            expt{x}.meanProj = loadtiff( meanProjPath );
        end
        
        % Make the final ROIs from the segmentation data
        ROI{x} = MakeROI3D(expt{x}, 'overwrite',false, 'corrPrct',75, 'minFoot',50);  % , 'corrPrct',90, [ROI{x}, preROI{x}] overwriteROI
        expt{x}.Nroi = numel(ROI{x}); % Nauto(x)
        [expt{x}.roiProj, ~, expt{x}.roiLabel] = VisualizeSegmentation(expt{x}, ROI{x}, 'overwrite',false);   %  
        % Extract/process fluor signals
        for runs = flip(expt{x}.runs)
            fluor{x}(runs) = GetFluor3D(runInfo{x}(runs), 'overwrite',false);  % expt{x}, r
        end
        fluor{x} = GetROIfluor(expt{x}, catInfo{x}, ROI{x}, fluor{x}, deform{x}, loco{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'lp',0, 'deconvolve',false, 'overwrite',false); 
        % Show the final results
        defVars = {'transAP', 'transML', 'transMag', 'scaleAP', 'scaleML', 'scaleMag', 'stretchAP', 'stretchML', 'shearAP', 'shearML', 'shearMag', 'shiftZ', 'DshiftZ'}; %
        NdefVars = numel( defVars ); %, 'dShiftZ'
        allVars = [{'fluor'}, defVars, {'velocity'}]; NallVars = numel(allVars);
        viewLims = struct('trans',[-Inf,Inf], 'scale',[-Inf,Inf], 'stretch',[-Inf,Inf], 'shear',[-0.03,0.03], 'shift',[-3, 3], 'velocity',[-3, 15], 'fluor',[-5, 5]); 
        ViewResults3D( expt{x}, Tscan{x}, deform{x}, loco{x}, fluor{x}, allVars, ROI{x}, 'cat',true, 'limits',viewLims); %
    end
end
clearvars Expt;