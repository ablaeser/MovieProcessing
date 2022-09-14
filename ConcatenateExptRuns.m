function interRunShift = ConcatenateExptRuns(expt, runInfo, catInfo, varargin )
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
%addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'minInt', 1000, @isnumeric )
addParameter( IP, 'maxEdge', [100, 100, 120, 120], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'setEdge', [], @isnumeric )
addParameter( IP, 'refRun', 1, @isnumeric )
addParameter( IP, 'refChan', 'green', @ischar )
addParameter( IP, 'sbx', 'sbxz', @ischar )
%addParameter( IP, 'proj', 'interp', @ischar )
addParameter( IP, 'ext', 'sbxcat', @ischar )
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, expt, runInfo, catInfo, varargin{:} ); % mouse, exptDate,
minIntOrig = IP.Results.minInt;
if expt.Nchan == 2 && expt.Nplane > 1
    minInt = minIntOrig/2^8;
else
    minInt = minIntOrig;
end
edgeMax = IP.Results.maxEdge;
edgeSet = IP.Results.setEdge;
sbxType = IP.Results.sbx;
%projType = IP.Results.proj;
refRun = IP.Results.refRun;
refChan = IP.Results.refChan;
refChanInd = find(strcmpi(refChan, {'red','green'} ));
if isempty(refChanInd),  error('Invalid reference channel'); end
catExt = IP.Results.ext; %'.sbxcat'; % '.sbxcat'
overwrite = IP.Results.overwrite;

interRunShift = zeros(expt.Nruns, 3);

catName = expt.name; % sprintf('%s_FOV%i', , expt.fov);
catPathRoot = sprintf('%s%s', expt.dir, catName);
catSbxPath = strcat(catPathRoot, '.', catExt);
if expt.Nruns > 1
    catDir = strcat(expt.dir,'Concat\'); mkdir(catDir)
    if overwrite || ~exist(catSbxPath,'file')
        Nchan = runInfo(1).nchan; Nscan = [runInfo.Nscan];  %totScan = sum(Nscan);
        sbxPath = cell(1,expt.Nruns);  runProjPath = cell(1,expt.Nruns);  runProj = cell(expt.Nruns, 1);
        cropProj = cell(expt.Nruns, 1); shiftProj = cell(1,expt.Nruns); nullProj = cell(1,expt.Nruns); %hpProj = cell(expt.Nruns, 1);
        if isempty(edgeSet), projEdges = zeros( expt.Nruns, 4 ); else, projEdges = repmat(edgeSet, expt.Nruns, 1 ); end
        if expt.Nplane > 1
            Nplane = runInfo(1).Nplane;   % Nrow = runInfo(1).sz(1); Ncol = runInfo(1).sz(2);
            for runs = 1:expt.Nruns % expt.runs
                sbxPath{runs} = sprintf('%s%s.%s', runInfo(runs).dir, runInfo(runs).fileName, sbxType ); % sbx_affine  sbxz
                runProj{runs} = WriteSbxProjection(sbxPath{runs}, runInfo(runs), 'verbose',true, 'chan',refChan, 'monochrome',true, 'RGB',false, 'type','Z', 'overwrite',false);
                if ndims(runProj{runs}) == 4, runProj{runs} = squeeze(runProj{runs}(:,:,refChanInd,:)); end
                if isempty(edgeSet)
                    projEdges(runs,:) = GetEdges( runProj{runs}(:,:,round(expt.Nplane/2)), 'minInt',minInt, 'show',true );
                else
                    ShowEdges(projEdges(runs,:), runProj{runs}(:,:,round(expt.Nplane/2)));
                    %pause;
                end
                %saveastiff( runProj{r}, sprintf('%sTest\\%s_Proj_run%i.tif', expt.dir, expt.name, r ) );
            end
            close all;
            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);


            %filtSigma = 5;
            for runs = 1:expt.Nruns % expt.runs
                cropProj{runs} = runProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                %hpProj{r} = cropProj{r} - imgaussfilt(cropProj{r}, filtSigma);
                %saveastiff( uint16(cropProj{runs}), sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, runs ) );
                %saveastiff( hpProj{r}, sprintf('%s%s_HP_run%i.tif', catDir, expt.name, r ) );
            end

            % Calculate shifts between runs
            zUse = 1:Nplane;
            %refRun = 2;  %10:Nplane;
            interRunDFTshift = zeros(expt.Nruns,3);  interRunZshift = zeros(numel(zUse), expt.Nruns);
            for runs = 1:expt.Nruns
                if runs ~= refRun
                    if Nplane > 1
                        interRunDFTshift(runs,:) = dftregistration3D(fftn(cropProj{refRun}(:,:,zUse)),  fftn(cropProj{runs}(:,:,zUse)), 4);
                        interRunZshift(:,runs) = InterpolateZshift(cropProj{refRun}(:,:,zUse),  cropProj{runs}(:,:,zUse), 2);
                    else
                        zUse = 1;
                        tempOut = dftregistration(fft2(cropProj{refRun}),  fft2(cropProj{runs})); % [error,diffphase,net_row_shift,net_col_shift]
                        interRunDFTshift(runs,1:2) = tempOut(3:4);
                    end
                end
            end
            %interRunShift = zeros(expt.Nruns,3);
            interRunShift(:,[1,2,3]) = interRunDFTshift(:,[2,1,3]); % imtranslate expects [x,y,z], not [y,x,z]
            interRunShift = round( interRunShift );

            % Apply shifts to each run's projection
            shiftProj = cell(expt.Nruns,expt.Nchan);
            for runs = 1:expt.Nruns
                if Nplane > 1
                    shiftProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,:) );
                else
                    shiftProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,1:2) );
                end
                saveastiff( [cropProj{refRun}, shiftProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, runs ) );
            end
            shiftProjCat = cat( 4, shiftProj{:} ); % projCat = cat( 4, nullProj{:} );
            for z = 1:Nplane
                %saveastiff( squeeze( projCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%sTest\\%s_ind_z%i.tif', expt.dir, expt.name, z ) )
                saveastiff( squeeze( shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_cat_z%i.tif', catDir, expt.name, z ) )
            end

            % Make metadata file
            %catInfo = ConcatenateRunInfo(expt, runInfo, 'suffix','sbxcat'); %SpoofSBXinfo3D(Nx, Ny, Nplane, totScan, Nchan);
            %save( catInfoPath, 'info' );
        else
            for runs = expt.runs
                sbxPath{runs} = runInfo(runs).path; %sprintf('%s%s.sbx', runInfo(r).dir, runInfo(r).fileName );
                runProj{runs} = WriteSbxProjection(runInfo(runs).path, runInfo(runs), 'binT',1, 'verbose',true, 'chan',refChan, 'monochrome',true, 'type','dft', 'overwrite',false);
                if isempty(edgeSet)
                    projEdges(runs,:) = GetEdges( runProj{runs}(:,:,end), 'minInt',minInt, 'show',true );
                else
                    ShowEdges(projEdges(runs,:), runProj{runs});
                end
            end
            close all;

            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);
            % Crop the reference projections
            for runs = expt.runs
                cropProj{runs} = runProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                %saveastiff( cropProj{runs}, sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, runs ) );
            end
            close all;

            % Calculate shifts between runs
            refFT = fft2(cropProj{refRun});
            interRunDFTshift = zeros(expt.Nruns,2);
            for runs = 1:expt.Nruns
                if runs ~= refRun
                    tempOut = dftregistration(refFT,  fft2(cropProj{runs})); % [error,diffphase,net_row_shift,net_col_shift]  hpProj
                    interRunDFTshift(runs,1:2) = tempOut(3:4);
                end
            end
            interRunShift = round( flip(interRunDFTshift, 2) ); % flip shifts to be compatible with imtranslate

            % Apply shifts to each run's projection
            for runs = 1:expt.Nruns
                nullProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,1:2) );
                %saveastiff( nullProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:), sprintf('%s%s_nullProj_run%i.tif', catDir, expt.name, r ) );
                shiftProj{runs} = nullProj{runs};
                %saveastiff( [cropProj{refRun}, shiftProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, r ) );
            end
            projCat = uint16(cat( 4, nullProj{:} )) ;
            shiftProjCat = uint16(cat( 4, shiftProj{:} ));
            for z = 1:expt.Nplane
                saveastiff( squeeze( projCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_ind_z%i.tif', catDir, expt.name, z ) )
                saveastiff( squeeze( shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_cat_z%i.tif', catDir, expt.name, z ) )
            end
        end

        % Write the sbxcat file
        fprintf('\n     Writing %s\n', catSbxPath); tic
        rw = SbxWriter(catSbxPath, catInfo, catExt, true); % pipe.io.RegWriter(catSbxPath, catInfo, catExt, true);
        w = waitbar(0, sprintf('writing %s',catExt));
        for runs = 1:expt.Nruns %expt.runs
            % Load run stack
            fprintf('\n   Loading %s... ', sbxPath{runs}); tic
            if expt.Nplane > 1
                runStack = readSBX(sbxPath{runs}, runInfo(runs), 1, Nscan(runs), -1, []); % [c,x,y,z,t]
                if Nchan == 2
                    runStack = permute(runStack, [2,3,4,5,1]); % [x,y,z,t,c]
                    % Apply shifts to each scan
                    if ~all(interRunShift(runs,:) == 0)
                        for s = 1:Nscan(runs)
                            for c = 1:Nchan
                                runStack(:,:,:,s,c) = imtranslate( runStack(:,:,:,s,c), interRunShift(runs,:)  );
                            end
                        end
                    end
                    runStack = permute(runStack, [5,1,2,3,4]); % [c,x,y,z,t]
                    runStack = reshape(runStack, [size(runStack,[1,2,3]), prod(size(runStack,[4,5]))]);
                    rw.write( runStack ); %
                elseif Nchan == 1
                    % Apply shifts to each scan
                    if ~all(interRunShift(runs,:) == 0)
                        for s = 1:Nscan(runs)
                            runStack(:,:,:,s) = imtranslate( runStack(:,:,:,s), interRunShift(runs,:)  );
                        end
                    end
                    % Write the data to sbxcat
                    runStack = reshape(runStack, [size(runStack,[1,2]), prod(size(runStack,[3,4]))] );
                    rw.write( runStack ); % rw.write(squeeze(uint16(tempScan)));
                end
            else
                if Nchan == 2
                    runStack = readSBX(sbxPath{runs}, runInfo(runs), 1, Nscan(runs), -1, []); % [c,x,y,t]
                    runStack = permute(runStack, [2,3,4,1]); % [x,y,t,c]
                    % Apply shifts to each scan
                    if ~all(interRunShift(runs,:) == 0)
                        for s = 1:Nscan(runs)
                            for c = 1:Nchan
                                runStack(:,:,s,c) = imtranslate( runStack(:,:,s,c), interRunShift(runs,:)  );
                            end
                        end
                    end
                    runStack = permute(runStack, [4,1,2,3]); % [c,x,y,t]
                    rw.write( runStack ); %
                elseif Nchan == 1
                    runStack = readSBX(sbxPath{runs}, runInfo(runs), 1, Nscan(runs), -1, []); % [x,y,t]
                    % Apply shifts to each scan
                    if ~all(interRunShift(runs,:) == 0)
                        for s = 1:Nscan(runs)
                            runStack(:,:,s) = imtranslate( runStack(:,:,s), interRunShift(runs,:)  );
                        end
                    end
                    rw.write( runStack ); % Write the data to sbxcat
                end
            end
            waitbar( runs/expt.Nruns, w );
            toc
        end
        rw.delete;
        delete(w);
    else
        fprintf('\n%s already exists!', catSbxPath);
    end
    %{
    for z = [2,4,6,8,10,14]
        WriteSbxPlaneTif(catSbxPath, catInfo, z, 'dir',catDir, 'name',[catName,'_cat'], 'overwrite',true  );
    end
    %}
    %WriteSbxProjection(catSbxPath, catInfo, catProjPath, 'dir',catDir, 'name',[catName,'_cat'], 'overwrite',overwrite, 'binT',16); % , 'firstScan',100, 'Nscan',500
elseif ~exist(catSbxPath, 'file') || overwrite
    sbxSourcePath = runInfo.path; %strcat(,'z');
    fprintf('\nSingle run experiment: copying %s to %s', sbxSourcePath, catSbxPath);
    copyfile(sbxSourcePath, catSbxPath);
    %{
    % copy-paste the projections too
    if expt.Nplane == 1
        FileFinder(runInfo.dir, 'contains','')
    else
        FileFinder(runInfo.dir, 'contains','_Z_')
    end
    %}
    %projSourcePath = strcat(runInfo.dir, runInfo.fileName, '_interpProj.tif');
    %copyfile(projSourcePath, catProjPath);
end
end