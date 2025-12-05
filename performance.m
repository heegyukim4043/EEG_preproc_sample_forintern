%% load and calculate offline acc
clear; clc;
% Note: ft_preproc_bandpassfilter contains filtfilt() for zero-phase
% forward-reverse IIR filter
freq_stop = [59.5 60.5];
nbHarmonics = 4;
% SSVEP_ch = {'O1','O2','PO3','PO4','Oz'};
load('ch_biosemi_32.mat');  % biosemi32_locs
label = re_ch(1:12);
label(13) = {'X'};
label(14:32) = re_ch(13:31);



ch_select = {'P1','P2','PZ','PO3','PO4','PO7','PO8','POZ','O1','O2','OZ'};
% ch_select = {'P4','P3','PZ','O1','O2'};
% FB_sub = [13 70;27 70;41 70];
grid_var = 0:0.25:2;
% frame = [0.1 5.1] ;
stim_freq = [];
for i = 1 : 5
    stim_freq = cat(2,stim_freq, [14:21]+0.2*(i-1) );
end
% FB_sub = [stim_freq(1)-1 nbHarmonics*ceil(stim_freq(end))+1;...
%     2*stim_freq(1)-1 nbHarmonics*ceil(stim_freq(end))+1;...
%     3*stim_freq(1)-1 nbHarmonics*ceil(stim_freq(end))+1;...
%     4*stim_freq(1)-1 nbHarmonics*ceil(stim_freq(end))+1];
for fb_idx =1 : nbHarmonics
    FB_sub(fb_idx,:) =[ fb_idx*stim_freq(1)-1, nbHarmonics*ceil(stim_freq(end))+1];
end
stim_type = 11: 50;
% sub_name = {'hg','ck','sh','sh_l','jy'};

num_trial = 6;
%%  FB -CCA
for nbsub =1 
    
    for i = [1: num_trial]
        fname = sprintf('./data/run_feedback_%d.gdf',i);
        EEG(i) = pop_biosig(fname);
        disp(EEG);
    end
    
    for frame_idx = 1 : 10
                frame = [0.1 0.5*frame_idx+0.1] ;
%         frame = [0.1 2.1] ;
        
        
        %             fname = sprintf('./SSVEP_merged/s%02d.mat', nbsub);
        
        %     EEG.raw = reref(EEG.raw, []);
        %     EEG.raw = EEG.raw(:, 1:750, :);
        %     EEG.frame = [0 2.5];
        
        t= frame(1):1/EEG(1).srate:frame(2);
        t= t(1:end-1);
        
        max_len= 0;
        for run_idx =1:num_trial
            clear temp_data;
            temp_data= EEG(run_idx).data;
            %             disp(size(temp_data,2))
            if size(temp_data,2)>max_len
                max_len = size(temp_data,2);
            end
        end
        
        EEG_ch = label;
        clear temp_SSVEP_EEG;
        for i = 1 :num_trial
            clear Dataout;
            [Dataout, Chanlocs] = reref( EEG(i).data, [2], 'keepref', 'on');
            if size(EEG(i).data,2)< max_len
                temp_SSVEP_EEG(:,:,i) = cat(2,Dataout(:,:),zeros(size(Dataout,1),max_len-size(EEG(i).data,2)));
            else
                temp_SSVEP_EEG(:,:,i) = Dataout(:,:);
            end
        end
        
        temp_dataset=[];
        for run_idx =1:num_trial
            temp_event = [EEG(run_idx).event.type];
            temp_latency = [EEG(run_idx).event.latency];
            for i = 1 : length(stim_type)
                %                         temp_stim_type=stim_type(i);
                temp_stim_type= find(temp_event == stim_type(i) );
                list_type(i,:) = temp_stim_type(1);
            end
            for i = 1 : length(stim_type)
                list_latency(i,:) = temp_latency(list_type(i,:));
            end
            for class_idx = 1 : length(stim_type)
                for trial_idx = 1: size(list_type,2)
                    start_idx = frame(1)*EEG(1).srate+ list_latency(class_idx,trial_idx);
                    end_idx = (frame(2)-frame(1))*EEG(1).srate + start_idx;
                    temp_dataset(:,:,run_idx,class_idx,trial_idx) = temp_SSVEP_EEG(:,start_idx:end_idx-1,run_idx);
                end
            end
        end
        SSVEP_EEG = reshape(temp_dataset,size(temp_dataset,1),size(temp_dataset,2),size(temp_dataset,3)*size(temp_dataset,4));
        
        %             labels = EEG.true_label;
        labels_true=[];
        for i = 1 : size(temp_dataset,4)
            labels_true = cat(2,labels_true,i*ones(1,num_trial));
        end
        EEG_filtered = [];
        
        reference_sig = []; % [2*nb_harmonics x time x class]
        for k_class = 1:size(temp_dataset,4)
            reference_tmp = [];
            for k_harmonics = 1:nbHarmonics
                tmp_sincos = [sin(2*pi*stim_freq(k_class)*k_harmonics*t);cos(2*pi*stim_freq(k_class)*k_harmonics*t)];
                reference_tmp = cat(1, reference_tmp, tmp_sincos);
            end
            reference_sig = cat(3, reference_sig, reference_tmp);
        end
        
        predict_offline = [];
        rho_class = [];
        % FBCCA and standard CCA
        for ind_a= 1
            for ind_b =1 
%                 a = 0 + 0.25*(ind_a-1);
%                 b = 0 + 0.25*(ind_b-1);
                a= 1;
                b= 0;
                for nb_FB = 1:size(FB_sub,1)

                    wn(nb_FB) = nb_FB^(-a)+b;
                    EEG_filtered = [];
                    % step1. pre-processing (filtering)
                    for i=1:size(SSVEP_EEG, 3)
                        tmp_filtered = ft_preproc_bandpassfilter(SSVEP_EEG(:,:,i), EEG(1).srate, FB_sub(nb_FB,:), 4, 'but');
                        tmp_notch = ft_preproc_bandstopfilter(tmp_filtered, EEG(1).srate, freq_stop, 4, 'but');
                        
                        EEG_filtered = cat(3, EEG_filtered, tmp_notch);
                    end
                    
                    for fft_trial_idx =1 : size(EEG_filtered,3)
                        fft(EEG_filtered(:,:,fft_trial_idx))/length(t);
                        
                    end
                    % w/ cross-validation
                    
                    % w/o cross-validation
                    % step2. assign class and perform CCA
                    for j_trial=1:size(EEG_filtered,3)
                        EEG_test = EEG_filtered(ismember(EEG_ch,ch_select), :, j_trial);
                        for k_class=1:size(reference_sig,3)
                            [Wx, Wy, rho] = canoncorr(EEG_test', reference_sig(:,:,k_class)');
                            %                 [Wx_hat, Wy_hat, rho_hat] = cca(EEG_test, reference_sig(:,:,k_class));
                            rho_class(nb_FB, j_trial, k_class) = max(rho);
                        end
                        [~, predict_offline(j_trial)] = max(squeeze(rho_class(nb_FB, j_trial, :)));
                    end
                    % CCA acc for each sub band
                    offline_acc(nb_FB, nbsub) = sum(predict_offline == labels_true) / length(predict_offline);
                end
                % sub-band operation ends
                
                % calculate weighted rho for FBCCA
                weighted_rho = zeros(length(labels_true), size(reference_sig,3));
                for nb_FB = 1:size(FB_sub,1)
                    weighted_rho = weighted_rho + wn(nb_FB) * squeeze(rho_class(nb_FB, :, :)).^2; % [trial x class]
                end
                [~, predict_fbcca] = max(weighted_rho, [], 2);
                offline_acc_fbcca(nbsub,frame_idx) = sum(predict_fbcca'== labels_true) / length(predict_fbcca);
%                  offline_acc_fbcca_grid(nbsub,ind_a,ind_b) = sum(predict_fbcca'== labels_true) / length(predict_fbcca);
                disp('done');
            end
            %         offline_acc_cca(:) = offline_acc(1, :);
        end
    end
end