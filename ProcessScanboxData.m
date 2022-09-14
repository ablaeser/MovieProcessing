clear; clc; close all;
dataDir = 'D:\2photon\'; % 'D:\2photon\Simone\'; %'C:\2photon';
dataSet = 'Vasculature'; % 'Astrocyte'; %  'Afferents'; %  'Neutrophil_Simone'; %  'NGC'; % 'Neutrophil'; % 
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
regParam.refChan = ''; % 'green';  %  'red'; %
regParam.refScan = [];
regParam.refRun = NaN; % only used if refScan is empty
regParam.window = 99; 
regParam.chunkSize = 1000; %
regParam.histmatch = true; % false; %
regParam.avgT = 0; % 15;  temporal Gaussian filter width (not downsampling). scans, not seconds
regParam.avgTsigma = 0; % 5
regParam.binXY = 2; % spatial downsampling
regParam.binT = 1; % temporal downsampling
regParam.prereg = false;
regParam.highpass = 0;
regParam.lowpass = 0;
regParam.medFilter = [0,0,0];
regParam.minInt = 1500;
regParam.edges = [80,80,40,20]; % [60,100,60,60]; %[100,120,10,80]; %[120,100,30,30]; % [80,90,20,20]; %[60,30,20,20]; % [80,90,20,20]; % [20,60,20,20];
regParam.name = '';
regParam.method = 'affine'; %'translation'; % 
% Projection parameters
projParam.type = 'max';
projParam.umPerPixel_target = 1; % target um/pix for projections. spatial downsampling 
projParam.rate_target = 2; % target scanrate for projections, in Hz. temporal downsampling
projParam.color = {'red','green'};
projParam.overwrite = true;

% Various useful subsets of the data to choose from
%xDone = find([dataTable{:,dataCol.done}] > 0);
%x3D = intersect( xDone, find([dataTable{:,dataCol.volume}] == 1 ));
%x2D = intersect( xDone, find( [dataTable{:,dataCol.volume}] == 0 ));
% Choose which subset to  process
xPresent = 25; % flip(100:102); %45:47; % [66:69]; %6;  62,64,
Npresent = numel(xPresent);
overwrite = false;
for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Parse data table
    [expt{x}, runInfo{x}, regParam, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir, regParam, projParam);
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x}); 

    % Get locomotion data
    for runs = flip(expt{x}.runs) % expt{x}.runs %
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',true );
    end
    % Determine reference run and scans (longest pre-CSD epoch of stillness)
    [regParam.refRun, regParam.refScan] = DetermineReference(expt{x}, Tscan{x}, loco{x});
    expt{x}.refRun = regParam.refRun; % vestigal at this point?

    % Register individual runs
    %{
    runInfo{x}(regParam.refRun) = RegisterRun( runInfo{x}(regParam.refRun), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid');
    for runs = flip(setdiff(1:expt{x}.Nruns, regParam.refRun)) % expt{x}.runs %
        runInfo{x}(runs) = RegisterRun( runInfo{x}(runs), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid'); %'rigid' , 'edges',[80,80,40,20]
        %WriteSbxPlaneTif(runInfo{x}(runs).path, runInfo{x}(runs), 1, 'chan',regParams.refChan, 'edges',regParams.edges, 'binT',8, 'monochrome',true, 'RGB',false, 'overwrite',false, 'verbose',true);
    end
    %}
    
    % Concatenate runs and metadata
    catInfo(x) = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat', 'overwrite',overwrite); % Get concatenated metadata
    interRunShift = ConcatenateExptRuns(expt{x}, runInfo{x}, catInfo(x), 'refRun',regParam.refRun, 'refChan',expt{x}.refChan, 'setEdge',[70   160   100    80]); %1 
    WriteSbxProjection(expt{x}.sbx.cat, catInfo(x), 'chan','both', 'type','cat', 'overwrite',overwrite, 'monochrome',true, 'RGB',true); % , 'bin',50
    
    % Save metadata to avoid having to rerun all this code of other scripts
    %{
    exptSavePath = sprintf('%s%s_expt.mat', expt{x}.dir, expt{x}.name);
    Expt = expt{x};
    save(exptSavePath, 'Expt');
    %}

    % Use the longest period of stillness to define the reference image
    % Register the concatenated data (see RegisterCat3D, AlignPlanes and RegisterSBX for more info)   
    fprintf('\n   Performing planar registration... '); % (reference averaged over scans %i - %i)
    regParam.refChan = expt{x}.refChan;
    regParam.refScan = regParam.refScan + expt{x}.scanLims(regParam.refRun);
    regParam = AlignPlanes( expt{x}.sbx.cat, catInfo(x), regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg );
    WriteSbxProjection(expt{x}.sbx.reg, catInfo(x), 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','reg', 'overwrite',overwrite);
    [~,deform{x}] = GetDeformCat3D( catInfo(x), 'show',true, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last') );  % [~, deform{x}, regParams, badInd] = 

    % Generate downsampled, possibly z-projected, movies for each run from the concatenated data
    % {
    projParam.edge = [80, 80, 40, 20];
    projParam.z = {6:7, 8:10}; 
    projParam.vol = true;
    projParam.overwrite = true;
    projParam.umPerPixel_target = 1;
    projParam.rate_target = 2;
    projParam = GenerateExptProjections(expt{x}, catInfo(x), Tscan{x}, projParam); %  
    %}
end
clearvars Expt;