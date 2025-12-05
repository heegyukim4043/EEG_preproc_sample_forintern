
    
    % load data after ICA
%     EEG = pop_loadset('XXXX.set');
%     EEG = pop_biosig
    
    EEG = pop_iclabel(EEG, 'default');
    
    rejIdx=[];
    cutProb=0.5; % 50
    for iICA = 1 : EEG.rank
        [maxProb maxIdx]= max(EEG.etc.ic_classification.ICLabel.classifications(iICA, :));
        % 1: brain / 2: Muscle / 3: Eye / 4: Heart / 5: Line Noise / 6: Channel Noise / 7: Other
        if maxIdx ~= 1 && maxIdx ~= 7 && maxProb > cutProb
            rejIdx = [rejIdx iICA];
        end
    end
%        pop_viewprops( EEG , 0 , 1 : EEG.rank);
    EEG.etc.rejIdx = rejIdx;
    EEG = pop_subcomp( EEG, rejIdx, 0);
    
