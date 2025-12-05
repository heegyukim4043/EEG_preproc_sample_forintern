%% path, init
addpath('eeglab2022.0');
eeglab;

clc;
clear;
close all;
%% initial value
clc; clear; close all;

data_dir = './9_speller';
% data_name ='1_jeh_M_filtered_modified';
data_name ='sample_data';
% sample_data
% data_name = 'record';

freq=[1 50]; % filter band
frame =[0 3]; % window, epoch length
m_stimuli = 14:22;
trial_lnum= 5;
%% load data csv
data_name_csv ='1_jeh_M_filtered_modified';
file_name= sprintf('%s/%s.csv',data_dir, data_name_csv);
EEG = pop_importdata('dataformat','ascii','nbchan',0,'data',file_name,'srate',300,'pnts',0,'xmin',0);
chan_list = {'P3','C3','F3','Fz','F4','C4','P4','Cz','Pz','Fp1','Fp2','T3','T5','O1','O2','F7','F8','T6','T4','A1','A2'};

for i = 1 : length(chan_list)
    EEG.chanlocs(i).labels = chan_list{i};
end

% a = pop_chanedit(EEG, 'lookup');
% figure,
% pop_eegplot(EEG);
% spectopo(EEG.data,0,EEG.srate);
%% load data gdf
file_name= sprintf('%s/%s.gdf',data_dir, data_name);

EEG=pop_biosig(file_name);
figure,
pop_eegplot(EEG);
spectopo(EEG.data,0,EEG.srate);

%% channel selection
% EEG = pop_select( EEG, 'channel',{'P3','C3','F3','Fz','F4','C4','P4','Cz','Pz','Fp1','Fp2','T3','T5','O1','O2','F7','F8','T6','T4'});
% EEG_1 = EEG;
% ssvep_ch = find(ismember({EEG.chanlocs.labels}, {'P3', 'P4', 'Pz', 'O1', 'O2'}));
aux_ch = find(ismember({EEG.chanlocs.labels}, {'X3', 'X2', 'X1', 'A2'}));
EEG.data(aux_ch, :) = [];
EEG.chanlocs(aux_ch)= [];
EEG.nbchan = size(EEG.data,1);

pop_eegplot(EEG);
%% reref - 

ref_ch = find(ismember({EEG.chanlocs.labels}, {'A1','A2'}));
EEG_reref = pop_reref(EEG,  ref_ch,'keepref','off');

pop_eegplot(EEG_reref);

%% band pass filter
% EEG_2 = EEG_1;

EEG = pop_eegfiltnew(EEG, 'locutoff',freq(1),'hicutoff',freq(2),'plotfreqz',0);
% EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',50,'plotfreqz',0);
EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1,'plotfreqz',0);

pop_eegplot(EEG);
% EEG.data= ft_preproc_bandpassfilter(EEG.data, EEG.srate, freq, 4, 'but');
% EEG.data= ft_preproc_bandstopfilter(EEG.data, EEG.srate, [59 61], 4, 'but');
% pop_eegplot(EEG);
figure,
spectopo(EEG.data,0,EEG.srate);


%% ICA 
a = pop_chanedit(EEG, 'lookup');
EEG.chaninfo.filename = a.chaninfo.filename;
EEG.chanlocs = a.chanlocs;
% EEG = pop_chanedit(EEG, 'lookup');
%Run the independent component analysis (ICA) 
EEG = pop_runica(EEG, 'runica');
%Use ICLabel to classify the independent components.
EEG = iclabel(EEG);
% Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
EEG = pop_icflag(EEG, [NaN NaN;0.80 1;0.80 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
% EEG = pop_selectcomps(EEG, [1:size(EEG.data,1)] ); % plot
EEG = pop_subcomp( EEG, [], 0);% means removing components flagged for rejection


% pop_eegplot(EEG);
save_name= sprintf('%s/%s.mat',data_dir, data_name);
save(save_name,'EEG','-v7.3');
%%
%% epoch
% load(save_name); % EEG


temp_events={EEG.event.type};
tmep_string= erase(temp_events,'condition ');

for i = 1 : size(tmep_string,2)
    event_marker(i) = str2num(tmep_string{i});
end

temp_latency = [EEG.event.latency];

for marker_index =1 : 9
    list_event(:,marker_index)=find(event_marker == marker_index+ 10) ;  %% marker: 11 - 19 
end

for marker_index =1 : 9
    latency(:,marker_index) = temp_latency(list_event(:,marker_index));
    event(:,marker_index) = event_marker(list_event(:,marker_index));  %% trial, class
end

re_data = [];
frame =[0 3]; % window, epoch length

for class =1: 9  %9
    trial_data =[];
    for trial = 1:5 %5
        start_id = latency(trial, class) + frame(1)*EEG.srate ;
        end_id = start_id + (frame(2)- frame(1))*EEG.srate -1 ;
        tmp_epoch1p = EEG.data(:, start_id:end_id);

        trial_data = cat(3,trial_data,tmp_epoch1p);
    end
    re_data = cat(4,re_data,trial_data);  % ch, time, trial, class
end 



%% frequeny 

srate = EEG.srate  ;               
T = 1/srate;            % Sampling period       
L = size(re_data,2);   % Length of signal
t = (0:L-1)*T; 
f = srate/L*(0:(L/2));

% band_range 
band_range = [4,8;8,12;12,30];
for band_idx =1 : size(band_range,1)
    band_sep_range{band_idx} = find(f>= band_range(band_idx,1) & f<= band_range(band_idx,2));
    
end


for subj = 1
    for class =1: 9  %9
        for trial = 1:5 %5
            
            [spectra,freqs,speccomp,contrib,specstd]= spectopo(re_data(:,:,trial,class),0, srate,'plot','off','winsize',L);
            
            re_spectra(:,:,trial,class,subj) = spectra;
        end
    end
    
end

a=rand(1,1001);
b=sort(a);
plot(a)
hold on
plot(b)
scatter(length(b)*0.9,b(round(end*90/100)));


plot(f,re_spectra(:,:,1,1))

xlim([0 50])

%% 

[ersp,itc,powbase,times,freqs,erspboot,itcboot]=newtimef(re_data(1,:,:,1),size(re_data,2),[0 3000], EEG.srate );

% alpha = rand(1,19);

topoplot(alpha*10,EEG.chanlocs);
colorbar;
caxis([0 10]);
