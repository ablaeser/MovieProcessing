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
parse( IP, sbxPath, sbxInfo, varargin{:} );

refChan = IP.Results.refChan;
[usePMT, ~] = DeterminePMT(refChan, sbxInfo);
refScan = IP.Results.refScan;
edges = IP.Results.edges;
scale = IP.Results.scale;
planescorr = IP.Results.planescorr;
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

% save the results 
shiftPath = strcat(sbxInfo.dir, sbxInfo.exptName,'_zinterp.mat'); %IP.Results.shiftpath;
save(shiftPath, 'RS','CS','ZS','scale', '-mat'); % 'RS1','CS1','ZS1','RS2','CS2','RS3','CS3','RS_chunk','CS_chunk','ZS_chunk','scale'



% Get z shift
%[RS_chunk, CS_chunk, ZS_chunk] = EstimateZshift(ref_vol, ref_cat, planescorr);


% Preallocating space for output variables
%{

for t=1:Size(4) % for each time frame
    % pick corresponding reference
    reft = ReferenceVolumes(:,:,:,ceil(t/Chunck));
    % preallocate output vectors
    row_shift = zeros((sbxInfo.Nplane-planescorr),1);
    column_shift = zeros((sbxInfo.Nplane-planescorr),1);
    z_shift = zeros((sbxInfo.Nplane-planescorr),1);

    for z = planescorr+1:sbxInfo.Nplane-planescorr % j = considered plane
        corr_z = ones(1,corr_length)*NaN;

        for i = z-planescorr:z+planescorr % i = reference plane
            ref_slice = reft(:,:,i);
            slice = Volume(:,:,z,t);
            corr_z(1,i-z+planescorr+1) = corr2(ref_slice,slice);
        end
        [~,J] = max(corr_z);
        % Set interpolation vectors and degree regarding matrix size
        if  (J-5 > 0) & (J+5 <= size(corr_z))
            idx = 5;
        else
            idx = min([length(corr_z)-J, J-1]);
        end
        % Find Z shift based on spatial correlations
        x = J-idx:0.01:J+idx;
        FitOrder = idx;
        P = polyfit(J-idx:J+idx, corr_z(J-idx:J+idx),FitOrder);
        CorrelationFit = polyval(P, x);
        [~,I] = max(CorrelationFit); % max of the polynomial curve
        % recalcutate x and y shifts
        output = dftregistrationAlex(fft2(reft(:,:,z)),...
            fft2(Volume(:,:,z,t)),100);
        row_shift(z,1) = output(1);
        column_shift(z,1) = output(2);
        z_shift(z,1) = x(I)-planescorr-1;
    end
    % add shifts into output matrix
    RowShifts(:,t) = cat(1,row_shift, zeros(planescorr,1));
    ColumnShifts(:,t) = cat(1,column_shift, zeros(planescorr,1));

    % ensuring strict monotony, necessary?
    psteps = ones(length(z_shift),1).*(1:length(z_shift))';
    zaux = -z_shift + psteps;
    count = 0;
    while ~issorted(zaux) & count < 10
        count = count + 1;
        for plane = 2:size(z_shift)
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
    ZShifts(:,t) = cat(1,-z_shift, zeros(planescorr,1));
end
end
%}

%{
RS = cell(3,Nchunk); CS = cell(3,Nchunk); ZS = cell(1,Nchunk);
tic
w = waitbar(0,'Z interpolation...');
for c = 1:5 %Nchunk %parfor
    % Load current chunk
    raw_chunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), refPMT ); % readSBX(sbxPath, sbxInfo, chunkFrames*(chunk-1)+1, chunkFrames, refChan, [] );  
    raw_chunk = raw_chunk(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:);
    raw_chunk = imresize(raw_chunk,1/scale);

    % First round of XY DFT registration
    % make single target volume for entire chunk  %ref1 = defineReference(raw_chunk, size(raw_chunk,4), refType);
    if strcmpi(refType, 'mean')
        ref1 = mean(raw_chunk, 4);
    elseif strcmpi(refType, 'median')
        ref1 = median(raw_chunk, 4);
    end
    [RS1, CS1] = DetermineXYShiftsFBS(raw_chunk, blurFactor, keepingFactor, ref1); %calculate DFT shift between each volume and target vol
    reg_chunk1 = ApplyXYShiftsFBS(raw_chunk,RS1,CS1); %apply the shift (need it for subsequent steps)
    %raw_chunk = []; %clear raw chunk to save space

    % First round of Z interpolation registration
    % make new reference volume based on XY DFT registration % ref2 = defineReference(reg_chunk1, size(reg_chunk1,4), refType); 
    if strcmpi(refType, 'mean')
        ref2 = mean(reg_chunk1, 4);
    elseif strcmpi(refType, 'median')
        ref2 = median(reg_chunk1, 4);
    end
    [RS2, CS2, ZS1] = InterpolateZshift(ref2, reg_chunk1, planescorr);
    %[RS2,CS2,ZS1] = ComputeZshiftInterpolateFBS(ref2, reg_chunk1, planescorr, edges); %calculate interpolated zshifts by fitting polynomial
    reg_chunk2 = ApplyZShiftInterpolateFBS(reg_chunk1, ZS1, CS2, RS2); %apply zshift

    %find NaNs and replace with zeros
    nan_idx = find(isnan(reg_chunk2));
    reg_chunk2(nan_idx) = 0;
    %reg_chunk1 = []; %clear registration from step 1

    % Final round of XY DFT registration
    %make new reference volume based on XYZ DFT (4)  %ref3 = defineReference(reg_chunk2, size(reg_chunk2,4), refType);
    if strcmpi(refType, 'mean')
        ref3 = mean(reg_chunk2, 4);
    elseif strcmpi(refType, 'median')
        ref3 = median(reg_chunk2, 4);
    end
    [RS3, CS3] = DetermineXYShiftsFBS(reg_chunk2, blurFactor, keepingFactor, ref3); %calculate DFT shift between each volume and target vol
    %reg_chunk2 = []; %clear to save space

    % Combine Row and Column Shifts from previous steps
    RS(:,c) = {RS1,RS2,RS3};
    CS(:,c) = {CS1,CS2,CS3};
    ZS(c) = {ZS1};

    % Retain final reference volumes for stitching later
    ref_all{c} = ref3;

    waitbar(c/Nchunk, w) 
end
delete(w);
toc
%Fix intra-chunk discontinuities
ref_final = ref_all{1};
%concatenate all of the reference volumes from above
ref_cat = ref_all{1};
for i = 2:Nchunk
    ref_cat = cat(4,ref_cat,ref_all{i});
end

%Compute XYZ shifts from each chunk's reference volume
[RS_chunk,CS_chunk,ZS_chunk] = InterpolateZshift(ref_final, ref_cat, planescorr); %ComputeZshiftInterpolateFBS(ref_final, ref_cat, planescorr, edges);
clearvars ref_cat %ref_cat = [];

%convert local shift correction cells to matrix form
CS1 = [CS{1,:}];
CS2 = [CS{2,:}];
CS3 = [CS{3,:}];
RS1 = [RS{1,:}];
RS2 = [RS{2,:}];
RS3 = [RS{3,:}];
ZS1 = [ZS{:}];

%stretch the intra-chunk corrections to apply to every frame
RS_chunk = imresize(RS_chunk,size(RS1),'nearest');
CS_chunk = imresize(CS_chunk,size(CS1),'nearest');
ZS_chunk = imresize(ZS_chunk,size(ZS1),'nearest');

%save the DFT and optotune transformations
if ~isempty(shiftPath) %saveToggle
    save(shiftPath,'RS1','CS1','ZS1','RS2','CS2','RS3','CS3','RS_chunk','CS_chunk','ZS_chunk','scale', '-mat');
end
%}
