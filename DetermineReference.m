function [refRun, refScan] = DetermineReference(expt, Tscan, loco)

    % Find the longest epoch of stillness, and use that period as the reference point for preregistration
    minStillDur = 30; trim = 5; %minTrim = 0.5;
    Nepoch = 0;
    while Nepoch == 0
        close all;
        [stillEpoch, ~, stillSumm] = BinStillEpochs(expt, Tscan, loco, [], [], [], 'criterion','speed', 'show',false, 'minStillDur',minStillDur, 'trim',trim); % true
        Nepoch =  stillSumm.Nepoch;
        minStillDur = minStillDur-1;
        if trim >= minStillDur/4, trim = minStillDur/10; end
        %trim = max([trim-0.5, minTrim]);
    end

    %[stillEpoch, ~, stillSumm] = BinStillEpochs(expt, Tscan, loco, [], [], [], 'criterion','speed', 'show',true, 'minStillDur',30, 'trim',5);
    if stillSumm.Nepoch > 0
        if ~isnan(expt.csd)
            [~,longestEpRunInd] = max(stillSumm.dur_trim(stillSumm.epoch_run < expt.csd)); % find the longest pre-CSD epoch
        else
            [~,longestEpRunInd] = max(stillSumm.dur_trim); % find the longest pre-CSD epoch
        end
        refRun = stillSumm.epoch_run(longestEpRunInd);
        [~,longestEpInd] = max(stillEpoch(refRun).Nscan_trim);
        refScan = stillEpoch(refRun).scan_trim{longestEpInd};
    end
end