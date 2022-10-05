function params = AlignPlanes(sbxInputPath, sbxInfo, params, varargin)
% Use TurboReg to affine register each plane of a volume scan movie
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'params', @isstruct )
addOptional( IP, 'refVol', [], @isnumeric )
addParameter( IP, 'outPath', '', @ischar )
addParameter( IP, 'sbx', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxInputPath, sbxInfo, params, varargin{:} ); 
refVol = IP.Results.refVol;
writeSbx = IP.Results.sbx;
overwrite = IP.Results.overwrite;
sbxOutputPath = IP.Results.outPath;
tformPath = sprintf('%s%s%s_regTforms.mat', sbxInfo.dir, sbxInfo.exptName, params.name); %strcat(sbxInfo.dir,'\',sbxInfo.exptName,'_regTforms.mat');

% If a regTforms file already exists, load it
if exist( tformPath, 'file' ) && ~overwrite
    fprintf('\nLoading %s... ', tformPath );
    load( tformPath ); % 'sbxPath', 'sbxInfo', 'refVol', 'params', 'tforms_all'
    zReg = find( any( cellfun(@isempty, regTform), 2 ) )';
    if ~isempty(zReg),  fprintf('\n   %s already exists - finishing off from z = %i', tformPath, zReg(1) ); end
else
    zReg = 1:sbxInfo.Nplane;
    regTform = cell(sbxInfo.Nplane, sbxInfo.totScan); 
end

if isempty( params.refScan )
    if sbxInfo.Nrun <= 2
        params.refRun = 1;
    else
        params.refRun = 2;
    end
    params.refScan = sbxInfo.scanLim(params.refRun)+2:sbxInfo.scanLim(params.refRun+1)-2;
end

% Check if data has been partially analyzed and pick up where it left off 
if ~isempty(params.name), params.name = strcat('_', params.name); end

%register one z plane at a time
if ~isempty(zReg)
    % Define the reference volume (use middle 50% of data by default)
    if isempty(refVol)
        if isempty(params.refScan) 
            params.refScan = ceil(sbxInfo.Nplane/4):floor(3*sbxInfo.Nplane/4); %   %ceil(sbxInfo.Nplane/2)-25:ceil(sbxInfo.Nplane/2)+25; 
        end
        fprintf('\nAveraging scans %i - %i for reference volume\n', params.refScan(1), params.refScan(end));
        refVol = WriteSbxProjection(sbxInputPath, sbxInfo, 'firstScan',params.refScan(1), 'Nscan',numel(params.refScan), 'type','ref', 'chan',params.refChan, ...
            'verbose',true, 'monochrome',true, 'RGB',false, 'overwrite',overwrite); % 
    end
    refVol = uint16(refVol);
    
    % Define edges
    if isempty(params.edges)
        params.edges = GetEdges( refVol(:,:,end), 'minInt',params.minInt, 'show',true ); 
    else
        ShowEdges( params.edges, refVol(:,:,round(sbxInfo.Nplane/2)) );
    end
    
    % Register plane by plane
    firstScan = 1; %sbxInfo.Nplane = 300;
    if isempty(params.name)
        nameStr = '';
    else
        nameStr = ['_',params.name];
    end
    regDir = strcat(sbxInfo.dir,'\RegTemp\'); mkdir(regDir);
    fix(clock)
    fprintf('\n');
    for z = 1:sbxInfo.Nplane % zReg
        fprintf('Registering plane %d of %d...  ',z, sbxInfo.Nplane);
        regName = sprintf('%s_z%02d%s', sbxInfo.exptName, z, nameStr);
        tic
        regTform(z,:) = RegisterSBX(sbxInputPath, sbxInfo, refVol(:,:,z), params, z, 'firstScan',firstScan, 'Nscan',sbxInfo.totScan, ...
            'dir',regDir, 'name',regName, 'verbose',false, 'intTif',false, 'finalTif',false, 'overwrite',overwrite);  % 
        toc
        fprintf('\nSaving %s ', tformPath);
        save(tformPath, 'sbxInputPath', 'sbxInfo', 'refVol', 'params', 'regTform', '-mat');
        toc
    end
    %{ 
    % Repair bad registration
    tempAff = cell(sbxInfo.Nplane, sbxInfo.Nplane); 
    zRepair = 6;
    repairScan = 350;
    Nrepair = 600;
    
    for z = zRepair
        tempAff(z,:) = AffineTurboReg(sbxInputPath, sbxInfo, refVol(:,:,z), params, z, 'firstScan',repairScan, 'sbxInfo.Nplane',Nrepair, ...
            'dir',regDir, 'name',sprintf('%s_z%02d', sbxInfo.exptName, z), 'verbose',false, 'intTif',true, 'finalTif',true); 
        %regTform(z,repairScan:repairScan+Nrepair-1) = 
    end
    for z = zRepair
        regTform(z,repairScan:repairScan+Nrepair-1) = tempAff(z,repairScan:repairScan+Nrepair-1);
    end
    
    for z = 1
        for s = 1:sbxInfo.totScan
            deformCatStruct.transAP(s,z) = regTform{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
            deformCatStruct.transML(s,z) = regTform{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
            deformCatStruct.scaleAP(s,z) = regTform{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
            deformCatStruct.scaleML(s,z) = regTform{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
            deformCatStruct.shearAP(s,z) = regTform{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
            deformCatStruct.shearML(s,z) = regTform{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
        end
    end
    plot( 
    %}
end

% Apply the transforms and generate sbx_affine
if writeSbx && (~exist(sbxOutputPath,'file') || overwrite)
    Nrow = sbxInfo.sz(1); 
    Ncol = sbxInfo.sz(2); 
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.totScan, 'N',10);
    tic
    w = waitbar(0,'writing .sbxreg');
    rw = SbxWriter(sbxOutputPath, sbxInfo, '.sbxreg', true); % rw = pipe.io.RegWriter(sbxOutputPath, sbxInfo, '.sbxreg', true);
    fprintf('\n     Writing registered sbx file'); 
    imRef = imref2d([Nrow,Ncol]);
    if sbxInfo.nchan == 1
        [pmt, ~] = DeterminePMT(params.refChan, sbxInfo);
        tic
        for chunk = 1:Nchunk
            data_chunk = readSBX(sbxInputPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), pmt, []);
            if sbxInfo.Nplane == 1
                % Single plane, single color
                for s = 1:chunkLength(chunk)
                    data_chunk(:,:,s) = imwarp(data_chunk(:,:,s), regTform{1,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                end
            else
                for s = 1:chunkLength(chunk)
                    for z = 1:sbxInfo.Nplane
                        data_chunk(:,:,z,s) = imwarp(data_chunk(:,:,z,s), regTform{z,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                    end
                end
                data_chunk = reshape(data_chunk, [size(data_chunk,[1,2]), prod(size(data_chunk,[3,4]))]);
            end

            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
            toc
        end
    else
        tic
        for chunk = 1:Nchunk
            data_chunk = readSBX(sbxInputPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), -1, []);
            if sbxInfo.Nplane == 1
                % Single plane, multi color
                data_chunk = permute(data_chunk, [2,3,4,1]);
                for s = 1:chunkLength(chunk)
                    for chan = 1:2
                        data_chunk(:,:,s,chan) = imwarp(data_chunk(:,:,s,chan), regTform{1,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                    end
                end
                data_chunk = permute(data_chunk, [4,1,2,3]);
            else
                % multi plane, multi color
                data_chunk = permute(data_chunk, [2,3,4,5,1]);
                for s = 1:chunkLength(chunk)
                    for chan = 1:2
                        for z = 1:sbxInfo.Nplane
                            data_chunk(:,:,z,s,chan) = imwarp(data_chunk(:,:,z,s,chan), regTform{z,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                        end
                    end
                end
                data_chunk = permute(data_chunk, [5,1,2,3,4]);
                data_chunk = reshape(data_chunk, [size(data_chunk,[1,2,3]), prod(size(data_chunk,[4,5]))]);
            end
            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
            toc
        end
    end
    rw.delete;
    delete(w);
    toc
end
end