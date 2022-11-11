function CatInterpZ(sbxPath, sbxInfo, varargin)
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addParameter(IP, 'chunkSize', 15, @isnumeric )
addParameter(IP, 'refChan', 'green', @ischar ) % for scanbox, 1 = green, 2 = red. -1 = both
addParameter(IP, 'refType', 'mean', @ischar); %options are 'median' or 'mean' for projecting the reference volume
addParameter(IP, 'refScan', [], @isnumeric )
addParameter(IP, 'scale', 4, @isnumeric ) %
addParameter(IP, 'edges',[0,0,0,0], @isnumeric); % [left, right, top, bottom]
addParameter(IP, 'planescorr', 3, @isnumeric); %how many planes up and down to search for Z DFT reg
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep',0.95, @isnumeric); %what proportion of frame to accountf or with shifts
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxPath, sbxInfo, varargin{:} );
refChan = IP.Results.refChan;
[usePMT, ~] = DeterminePMT(refChan, sbxInfo);
refScan = IP.Results.refScan;
edges = IP.Results.edges;
scale = IP.Results.scale;
planescorr = IP.Results.planescorr;
overwrite = IP.Results.overwrite;
% save the results
shiftPath = strcat(sbxInfo.dir, sbxInfo.exptName,'_zinterp.mat'); %IP.Results.shiftpath;
if ~exist(shiftPath,'file') || overwrite
    Ncheck = numel(planescorr+1:sbxInfo.Nplane-planescorr);
    % Generate reference volume
    fprintf('\nAveraging scans %i - %i for reference volume\n', refScan(1), refScan(end));
    ref_vol = WriteSbxProjection(sbxPath, sbxInfo, 'firstScan',refScan(1), 'Nscan',numel(refScan), 'type','zref', 'chan',refChan, 'verbose',true, 'monochrome',true, 'overwrite',false); % , 'edges',edges , 'scale',scale
    ref_vol = ref_vol(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:);
    ref_vol = imresize(ref_vol,1/scale);
    corr_length = 2*planescorr+1;
    % Interpolate each scan to the reference
    fprintf('\nEstimating z-shift relative to reference volume')
    w = waitbar(0,'Estimating z-shift relative to reference volume...');
    RS = zeros(sbxInfo.Nplane, sbxInfo.totScan); CS = zeros(sbxInfo.Nplane, sbxInfo.totScan); ZS = zeros(sbxInfo.Nplane, sbxInfo.totScan); % row shift, column shifts, z shifts
    for scan = 1:sbxInfo.totScan
        % Load the current scan
        %curr_scan = WriteSbxProjection(sbxPath, sbxInfo, 'firstScan',scan, 'Nscan',1, 'chan',refChan, 'verbose',false, 'monochrome',false, 'edges',edges, 'scale',scale, 'overwrite',false); %
        curr_scan = readSBX(sbxPath, sbxInfo, scan, 1, usePMT); % (path, info, firstScan, Nscan, pmt, z)
        curr_scan = curr_scan(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:);
        curr_scan = imresize(curr_scan,1/scale);
        for z = planescorr+1:sbxInfo.Nplane-planescorr % j = considered plane
            corr_z = nan(1,corr_length);
            for i = z-planescorr:z+planescorr % i = reference plane
                corr_z(1,i-z+planescorr+1) = corr2(ref_vol(:,:,i), curr_scan(:,:,z));
            end
            [~,J] = max(corr_z);
            % Set interpolation vectors and degree regarding matrix size
            if  J-5 > 0 && J+5 <= corr_length
                idx = 5;
            else
                idx = min([corr_length-J, J-1]);
            end
            % Find Z shift based on spatial correlations
            x = J-idx:0.01:J+idx;
            FitOrder = idx;
            P = polyfit(J-idx:J+idx, corr_z(J-idx:J+idx),FitOrder);
            CorrelationFit = polyval(P, x);
            [~,x_max_ind] = max(CorrelationFit); % max of the polynomial curve
            % recalcutate x and y shifts
            output = dftregistrationAlex(fft2(ref_vol(:,:,z)), fft2(curr_scan(:,:,z)), 100);
            row_shift(z,1) = output(1);
            column_shift(z,1) = output(2);
            z_shift(z,1) = x(x_max_ind)-planescorr-1;
        end
        z_shift_length = numel(z_shift);
        % add shifts into output matrix
        RS(:,scan) = cat(1,row_shift, zeros(planescorr,1));
        CS(:,scan) = cat(1,column_shift, zeros(planescorr,1));

        % ensuring strict monotony, necessary?
        psteps = (1:z_shift_length)'; %ones(z_shift_length,1).*(1:z_shift_length)';
        zaux = -z_shift + psteps;
        count = 0;
        while ~issorted(zaux) && count < 10
            count = count + 1;
            for plane = 2:z_shift_length
                if zaux(plane) - zaux(plane-1) <= 0
                    if ismember(count, [1,2,5,6,9,10]) % mod(count,2)==1;
                        zaux(plane-1) = NaN;
                    else
                        zaux(plane) = NaN;
                    end
                end
            end
            zaux = naninterp(zaux);
        end
        z_shift = zaux - psteps;
        ZS(:,scan) = cat(1,-z_shift, zeros(planescorr,1));
        waitbar(scan/sbxInfo.totScan, w)
    end
    save(shiftPath, 'RS','CS','ZS','scale', '-mat'); % 'RS1','CS1','ZS1','RS2','CS2','RS3','CS3','RS_chunk','CS_chunk','ZS_chunk','scale'
else
    fprintf('\n%s already exists\n', shiftPath);
    %load(shiftPath)
end