function [expt, runInfo, varargout] = ParseDataTable(dataTable, x, dataCol, dataDir, varargin) % dataDir, projParam  dataTablePath, dataSet,
IP = inputParser;
addRequired( IP, 'dataTable', @iscell )
addRequired( IP, 'x', @isnumeric )
addRequired( IP, 'dataCol', @isstruct )
addRequired( IP, 'dataDir', @ischar )
addOptional( IP, 'regParams', struct(), @isstruct)
addOptional( IP, 'projParam', struct(), @isstruct)
parse( IP, dataTable, x, dataCol, dataDir, varargin{:} );  % , ROI
regParams = IP.Results.regParams;
projParam = IP.Results.projParam;

% Experiment metadata
expt.mouse = dataTable{x,dataCol.mouse};
expt.date = dataTable{x,dataCol.date};  if isnumeric(expt.date), expt.date = num2str(expt.date); end
if any(ismissing(dataTable{x,dataCol.FOV})) % check if an FOV # is specified
    expt.fov = 1;
    expt.dir = sprintf('%s%s\\%s\\', dataDir, expt.mouse, expt.date); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt.mouse, expt.date, expt.fov, expt.mouse);
    expt.name = sprintf('%s_%s', expt.mouse, expt.date); 
else
    expt.fov = dataTable{x,dataCol.FOV};
    %if ischar(expt.fov), expt.fov = str2double(expt.fov); end
    if isnumeric(expt.fov), expt.fov = num2str(expt.fov); end
    expt.dir = sprintf('%s%s\\%s_FOV%s\\', dataDir, expt.mouse, expt.date, expt.fov); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt.mouse, expt.date, expt.fov, expt.mouse);
    expt.name = sprintf('%s_%s_FOV%s', expt.mouse, expt.date, expt.fov); 
end

if exist(expt.dir, 'dir')
    [~,exptPath] = FileFinder(expt.dir, 'contains','expt', 'type','mat');
    if ~isempty(exptPath)
        fprintf('\nLoading %s',exptPath{1})
        load(exptPath{1});
        expt = Expt;
    else
        fprintf('\n\nParsing %s\n', expt.name )
        runFolders = FileFinder(expt.dir, 'contains','run', 'type',0);
        [expt.runs, ~] = sort(cellfun(@GetRunNumber, runFolders, 'UniformOutput',true)', 'ascend'); %#ok<UDIM>  sortInd
        expt.Nruns = numel(expt.runs);
        expt.refChan = dataTable{x,dataCol.ref}; % determine which channel to use as reference for registration and concatenation
        % Check for CSD runs
        missingData = cellfun(@all, cellfun(@ismissing, dataTable, 'UniformOutput',false), 'UniformOutput',true );
        if isempty(dataCol.csd) || missingData(x,dataCol.csd)
            expt.csd = NaN;  expt.Ncsd = 0;
        else
            expt.csd = dataTable{x,dataCol.csd};  expt.Ncsd = numel(expt.csd);
        end
        if isnan(expt.csd), expt.preRuns = expt.runs; else, expt.preRuns = 1:expt.csd(1)-1; end
        expt.vasc_chan = 'red'; % dataTable{x,dataCol.vasc};
    end

    % Get run metadata
    runInfo = GetRunInfo(expt, dataDir);

    % Concatenate runs metadata
    expt.Nrow = runInfo(1).sz(1); expt.Ncol = runInfo(1).sz(2);
    if all([runInfo.otlevels]==runInfo(1).otlevels) && all([runInfo.nchan]==runInfo(1).nchan)
        expt.Nplane = runInfo(1).Nplane; % otlevels
        expt.Nchan = runInfo(1).nchan;
    else
        error('All runs must have the same number of planes and channels')
    end
    expt.Nscan = floor([runInfo.nframes]/expt.Nplane); expt.totScan = sum(expt.Nscan); expt.totFrame = sum([runInfo.totFrame]);
    expt.scanLims = [0, cumsum(expt.Nscan)];
    expt.frameRate = runInfo(1).framerate;
    expt.scanRate = runInfo(1).framerate/expt.Nplane;
    expt.sbx.cat = strcat(expt.dir, expt.name, '.sbxcat '); %'.sbx_interp'
    expt.sbx.reg = strcat(expt.dir, expt.name, '.sbxreg ');
    if isfield(runInfo(1).config, 'magnification_list')
        expt.zoom = str2double(runInfo(1).config.magnification_list(runInfo(1).config.magnification,:));
    elseif isfield(runInfo(1).config, 'magnification')
        warning('\nmagnification_list does not exist, using magnification field instead')
        %magList = [1, 1.2, 1.4, 1.7, 2, 2.4, 2.8, 3.4, 4, 4.8, 5.7, 6.7, 8];
        expt.zoom = runInfo(1).config.magnification; %magList()
    else
        expt.zoom = NaN;
    end
    expt.umPerPixel = (1/0.53)/expt.zoom;
    expt.refRun = floor(median(expt.runs));

    % Check if projection planes are specified
    if (expt.Nplane > 1) && ((isfield(dataCol, 'Ztop') && ~isempty(dataCol.Ztop)) || (isfield(dataCol, 'Zbot') && ~isempty(dataCol.Zbot)))
        zRange = sort([dataTable{x,dataCol.Ztop}, dataTable{x,dataCol.Zbot}]);
        projParam.z{1} = zRange(1):zRange(end);
    else
        projParam.z = {};
    end

    % Check if a reference color for registration is specified
    if isempty(dataCol.ref)
        regParams.refChan = 'green';
        warning('\nNo reference channel was specified, using %s by default', regParams.refChan);
    else
        regParams.refChan = dataTable{x,dataCol.ref};
    end

    %
    varargout{1} = regParams;
    varargout{2} = projParam;
else
    error('Experiment directory %s does not exist!', expt.dir)
end
end