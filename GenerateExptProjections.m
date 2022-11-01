function projParam = GenerateExptProjections(expt, catInfo, Tscan, varargin) % , runProj, regProj
% Generate downsampled, possibly z-projected, movies for each run from the concatenated data
if isempty(varargin)
    projParam.overwrite = false;
else
    projParam = varargin{1};
end
projParam.dir = strcat(expt.dir, 'Projections\'); mkdir(projParam.dir);
projParamPath = sprintf('%s%s_projParam.mat', projParam.dir, expt.name);
chanName = {'red','green'};
if nargin > 3  %~exist(projParamPath, 'file') % || projParam.overwrite
    runsDir = sprintf('%sRuns\\',projParam.dir); mkdir(runsDir);
    if expt.Nplane == 1, projParam.z = {1}; end
    projParam.Ncolor = numel(projParam.color);
    projParam.Nz = numel(projParam.z);
    projParam.bin = max([1, round(expt.scanRate/projParam.rate_target)]);
    projParam.rate_bin = expt.scanRate/projParam.bin;
    projParam.scaleFactor = projParam.umPerPixel_target/expt.umPerPixel; % round(projParam.umPerPixel_target/expt.umPerPixel); %
    if projParam.scaleFactor < 1, projParam.scaleFactor = 1; end
    projParam.umPerPixel_scale = expt.umPerPixel*projParam.scaleFactor;
    disp(projParam);

    % Indivdual run projections
    fprintf('\nGenerating projections:\n');
    Tcat = vertcat(Tscan{:});
    rawVol = cell(expt.Nruns, 1); regVol = cell(expt.Nruns, 1);
    rawRunProj = cell(expt.Nruns, 2, projParam.Nz); regRunProj = cell(expt.Nruns, 2, projParam.Nz); runBinLims = cell(1,expt.Nruns); %Tproj = cell(1,expt.Nruns);
    runCell = cell(expt.Nruns, 3, projParam.Nz); runVolCell = cell(expt.Nruns, 3); catCell = cell(projParam.Nz, 3);
    projParam.path.run.raw.z = runCell; projParam.path.run.reg.z = runCell;  % struct('raw',[], 'reg',[]);
    projParam.path.run.raw.vol = runVolCell; projParam.path.run.reg.vol = runVolCell;
    projParam.path.cat.raw.z = catCell; projParam.path.cat.reg.z = catCell;
    projParam.path.cat.raw.vol = cell(1,3); projParam.path.cat.reg.vol = cell(1,3);
    for runs = 1:expt.Nruns
        tempName = sprintf('%s_run%i',expt.name, runs); % _%s , projParam.color{c}
        if expt.Nplane == 1
            if exist(expt.sbx.cat, 'file')
                %[stackOut, stackChan, binLims, projPaths] = WriteSbxPlaneTif(sbxPath, sbxInfo, z, varargin) 
                [~, rawRunProj(runs,:,1), runBinLims{runs}, tempRawProjPath] = WriteSbxPlaneTif(expt.sbx.cat, catInfo, 1, 'chan','both', 'firstScan',expt.scanLims(runs)+projParam.bin+1, 'Nscan',expt.Nscan(runs)-projParam.bin, ...
                    'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'RGB',false, 'dir',runsDir, 'name',tempName, 'type','raw', 'overwrite',projParam.overwrite );
                %rawRunProj{runs,1,1} = tempRawProj(:,:,:,1);  rawRunProj{runs,2,1} = tempRawProj(:,:,:,2);
                projParam.path.run.raw.z(runs,:,1) = tempRawProjPath;
            end
            if exist(expt.sbx.reg, 'file')
                [~, regRunProj(runs,:,1), runBinLims{runs}, tempRegProjPath] = WriteSbxPlaneTif(expt.sbx.reg, catInfo, 1, 'chan','green', 'firstScan',expt.scanLims(runs)+projParam.bin+1, 'Nscan',expt.Nscan(runs)-projParam.bin, ...
                    'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'RGB',false, 'dir',runsDir, 'name',tempName, 'type','reg', 'overwrite',projParam.overwrite ); % 'both'
                %regRunProj{runs,1,1} = tempRegProj(:,:,:,1);  regRunProj{runs,2,Z} = tempRegProj(:,:,:,2);
                projParam.path.run.reg.z(runs,:,1) = tempRegProjPath;
            end
        else
            % Z projections
            if exist(expt.sbx.reg, 'file')
                [tempRegProj, runBinLims{runs}, tempRegProjPath] = WriteSbxZproj(expt.sbx.reg, catInfo, 'z',projParam.z, 'chan','both', 'dir',runsDir, 'name',tempName, 'sbxType','reg', 'projType',projParam.type, 'monochrome',true,...
                    'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite); % expt.Nscan(runs)
                for Z = 1:projParam.Nz
                    if size(tempRegProj{Z}, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                        regRunProj(runs,2,:) = tempRegProj;
                        projParam.path.run.reg.z(runs,2,:) = tempRegProjPath;
                    elseif size(tempRegProj, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                        regRunProj(runs,1,:) = tempRegProj;
                        projParam.path.run.reg.z(runs,2,:) = tempRegProjPath;
                    else
                        error('2 channel case needs to be fixed')
                    end
                end
            end
            if exist(expt.sbx.cat, 'file')
                [tempRawProj, runBinLims{runs}, tempRawProjPath] = WriteSbxZproj(expt.sbx.cat, catInfo, 'z',projParam.z, 'chan','both', 'dir',runsDir, 'name',tempName, 'sbxType','raw', 'projType',projParam.type, 'monochrome',true,...
                    'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite);
                for Z = 1:projParam.Nz
                    if size(tempRawProj{Z}, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                        rawRunProj(runs,2,:) = tempRawProj;
                        projParam.path.run.raw.z(runs,2,:) = tempRawProjPath;
                    elseif size(tempRawProj, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                        rawRunProj(runs,1,:) = tempRawProj;
                        projParam.path.run.raw.z(runs,2,:) = tempRawProjPath;
                    else
                        error('2 channel case needs to be fixed')
                        %rawRunProj(runs,:,:) = permute( tempRawProj(), [] ); % tempRawProj = 
                        %rawRunProj{runs,2,Z} = tempRawProj(:,:,:,2);
                    end
                end
            end

            %{
            if projParam.Nz > 0
                
                for Z = 1:projParam.Nz
                    if exist(expt.sbx.reg, 'file')
                        [tempRegProj, runBinLims{runs}, tempRegProjPath] = WriteSbxZproj(expt.sbx.reg, catInfo, 'z',projParam.z{Z}, 'chan','both', 'dir',runsDir, 'name',tempName, 'sbxType','reg', 'projType',projParam.type, 'monochrome',true,...
                            'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite); % expt.Nscan(runs)
                        if size(tempRegProj, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                            regRunProj{runs,2,Z} = tempRegProj;
                        elseif size(tempRegProj, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                            regRunProj{runs,1,Z} = tempRegProj;
                        else
                            regRunProj{runs,1,Z} = tempRegProj(:,:,:,1);
                            regRunProj{runs,2,Z} = tempRegProj(:,:,:,2);
                        end
                        projParam.path.run.reg.z(runs,:,Z) = tempRegProjPath;
                    end
                    if exist(expt.sbx.cat, 'file')
                        [tempRawProj, runBinLims{runs}, tempRawProjPath] = WriteSbxZproj(expt.sbx.cat, catInfo, 'z',projParam.z{Z}, 'chan','both', 'dir',runsDir, 'name',tempName, 'sbxType','raw', 'projType',projParam.type, 'monochrome',true,...
                            'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite);
                        if size(tempRawProj, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                            rawRunProj{runs,2,Z} = tempRawProj;
                        elseif size(tempRawProj, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                            rawRunProj{runs,1,Z} = tempRawProj;
                        else
                            rawRunProj{runs,1,Z} = tempRawProj(:,:,:,1);
                            rawRunProj{runs,2,Z} = tempRawProj(:,:,:,2);
                        end
                        projParam.path.run.raw.z(runs,:,Z) = tempRawProjPath;
                    end
                end
            end
            %}
            % Volume tifs
            if projParam.vol
                if exist(expt.sbx.reg, 'file')
                    [regVol{runs}, ~, projParam.path.run.reg.vol(runs,:)] = WriteSbxVolumeTif(expt.sbx.reg, catInfo, 'chan','both', 'dir',runsDir, 'name',tempName, 'type','reg', 'monochrome',true, 'RGB',true,...
                        'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'binT',projParam.bin, 'overwrite',projParam.overwrite);
                end
                if exist(expt.sbx.cat, 'file')
                    [rawVol{runs}, ~, projParam.path.run.raw.vol(runs,:)] = WriteSbxVolumeTif(expt.sbx.cat, catInfo, 'chan','both', 'dir',runsDir, 'name',tempName, 'type','raw', 'monochrome',true, 'RGB',true,...
                        'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'binT',projParam.bin, 'overwrite',projParam.overwrite); %projParam.overwrite
                end
            end
        end
        % Account for temporal binning in time vectors
        if projParam.bin == 1
            projParam.Tproj(runs) = Tscan(runs);
        else
            for b = flip(1:size(runBinLims{runs},1))
                projParam.Tproj{runs}(b) =  mean( Tcat(runBinLims{runs}(b,1):runBinLims{runs}(b,2)) );
            end
        end
    end
    projParam.Nbin = cellfun(@numel, projParam.Tproj);
    projParam.totBin = sum(projParam.Nbin);
    projParam.binLims = [0,cumsum(projParam.Nbin)];

    % Concatenate runs
    if expt.Nruns > 1
        fprintf('\nGenerating concatenated projections')
        if expt.Nplane == 1
            for chan = find(any(~cellfun(@isempty, rawRunProj(:,:,1)))) %find(~cellfun(@isempty, rawRunProj(:,:,1)))
                if exist(expt.sbx.cat, 'file')
                    projParam.path.cat.raw.z{chan,1} = sprintf('%s%s_raw_%s.tif', projParam.dir, expt.name, projParam.color{chan} );
                    if ~exist(projParam.path.cat.raw.z{chan,1}, 'file')
                        WriteTiff(cat(3, rawRunProj{:,chan,1}), projParam.path.cat.raw.z{chan,1} );
                    end
                end
                if exist(expt.sbx.reg, 'file')
                    projParam.path.cat.reg.z{chan,1} = sprintf('%s%s_reg_%s.tif', projParam.dir, expt.name, projParam.color{chan} );
                    if ~exist(projParam.path.cat.reg.z{chan,1}, 'file')
                        WriteTiff(cat(3, regRunProj{:,chan,1}), projParam.path.cat.reg.z{chan,1} );
                    end
                end
            end
        else
            % Concatenated z-projections
            for Z = 1:projParam.Nz
                for chan = find(all(~cellfun(@isempty, rawRunProj(:,:,Z))))
                    projParam.path.cat.raw.z{Z,chan} = sprintf('%s%s_%s_z%i-%i_raw_%s.tif', projParam.dir, expt.name, projParam.type, projParam.z{Z}(1), projParam.z{Z}(end), chanName{chan} ); % projParam.color{chan}
                    if ~exist(projParam.path.cat.raw.z{Z,chan}, 'file') % exist(expt.sbx.cat, 'file') &&
                        WriteTiff(cat(3, rawRunProj{:,chan,Z}), projParam.path.cat.raw.z{Z,chan} );
                    end
                end
                for chan = find(all(~cellfun(@isempty, regRunProj(:,:,Z))))
                    projParam.path.cat.reg.z{Z,chan} = sprintf('%s%s_%s_z%i-%i_reg_%s.tif', projParam.dir, expt.name, projParam.type, projParam.z{Z}(1), projParam.z{Z}(end), chanName{chan} );
                    if ~exist(projParam.path.cat.reg.z{Z,chan}, 'file') % exist(expt.sbx.reg, 'file') && 
                        WriteTiff(cat(3, regRunProj{:,chan,Z}), projParam.path.cat.reg.z{Z,chan} );
                    end
                end
            end
            % Concatenated Volume tifs
            if projParam.vol
                if exist(expt.sbx.cat, 'file')
                    projParam.path.cat.raw.vol{3} = sprintf('%s%s_raw_vol_RGB.tif', projParam.dir, expt.name );
                    if ~exist(projParam.path.cat.raw.vol{3},'file') || projParam.overwrite
                        tempCat = uint16(cat(5, rawVol{:}));
                        for chan = flip(1:size(tempCat,4))
                            projParam.path.cat.raw.vol{chan} = sprintf('%s%s_raw_vol_%s.tif', projParam.dir, expt.name, chanName{chan} );
                            if ~exist(projParam.path.cat.raw.vol{chan},'file') || projParam.overwrite
                                fprintf('\nWriting %s', projParam.path.cat.raw.vol{chan})
                                bfsave( tempCat(:,:,:,chan,:), projParam.path.cat.raw.vol{chan});
                            end

                            chanLower = prctile(reshape(tempCat(:,:,:,chan,:),[],1), 5);
                            chanUpper = max(reshape(tempCat(:,:,:,chan,:),[],1)); %prctile(stackChan{chan}(:), 1);
                            %fprintf('\nRescaling %s channel : [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper);
                            tempCat(:,:,:,chan,:) = rescale(tempCat(:,:,:,chan,:), 0, 2^8-1, 'inputMin',chanLower, 'inputMax',chanUpper);
                        end
                        fprintf('\nWriting %s', projParam.path.cat.raw.vol{3})
                        bfsave( uint8(tempCat), projParam.path.cat.raw.vol{3});
                    end
                end
                if exist(expt.sbx.reg, 'file')
                    projParam.path.cat.reg.vol{3} = sprintf('%s%s_reg_vol_RGB.tif', projParam.dir, expt.name );
                    if ~exist(projParam.path.cat.reg.vol{3},'file') || projParam.overwrite
                        tempCat = uint16(cat(5, regVol{:}));
                        for chan = flip(1:size(tempCat,4))
                            projParam.path.cat.reg.vol{chan} = sprintf('%s%s_reg_vol_%s.tif', projParam.dir, expt.name, chanName{chan} );
                            if ~exist(projParam.path.cat.reg.vol{chan},'file') || projParam.overwrite
                                fprintf('\nWriting %s', projParam.path.cat.reg.vol{chan})
                                bfsave( tempCat(:,:,:,chan,:), projParam.path.cat.reg.vol{chan});
                            end

                            chanLower = prctile(reshape(tempCat(:,:,:,chan,:),[],1), 5);
                            chanUpper = max(reshape(tempCat(:,:,:,chan,:),[],1)); %prctile(stackChan{chan}(:), 1);
                            %fprintf('\nRescaling %s channel : [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper);
                            tempCat(:,:,:,chan,:) = rescale(tempCat(:,:,:,chan,:), 0, 2^8-1, 'inputMin',chanLower, 'inputMax',chanUpper);
                        end
                        fprintf('\nWriting %s', projParam.path.cat.reg.vol{3})
                        bfsave( uint8(tempCat), projParam.path.cat.reg.vol{3});
                    end
                end
            end
        end
    end
    % Save the parameters, including Tproj
    fprintf('\nSaving %s', projParamPath)
    save(projParamPath, 'projParam');
else
    fprintf('\nLoading %s', projParamPath)
    load(projParamPath, 'projParam');
end