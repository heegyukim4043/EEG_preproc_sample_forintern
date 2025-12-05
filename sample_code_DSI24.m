%% path, init
addpath('eeglab2022.0');
addpath('fieldtrip-20190618');

eeglab;
ft_defaults;
clc;
clear;
close all;
%% initialize value
clc; clear; close all;

data_dir = './9_speller';
data_name ='sample_data'; % sprintf('sample%d',2)

% sprintf('x%d%d',1,2)

freq=[1 45]; % filter band
frame =[0 3]; % window, epoch length
% m_stimuli = 14:22;
trial_lnum= 5;



%% load data
file_name= sprintf('%s/%s.gdf',data_dir, data_name);
EEG = pop_biosig(file_name);
% size(EEG.data)
% EEG.nbchan


% EEG.data = awgn(EEG.data,10,'measured');
pop_eegplot(EEG);
figure,
spectopo(EEG.data,0,EEG.srate);  


%% channel selection
EEG_1 = EEG;

%  magic(3)
ssvep_ch = find(ismember({EEG_1.chanlocs.labels}, {'P3', 'P4', 'Pz', 'O1', 'O2'}));
aux_ch = find(ismember({EEG_1.chanlocs.labels}, {'X3', 'X2', 'X1', 'A2'}));

% t_mat=EEG.data(ssvep_ch,:);

EEG_1.data(aux_ch, :) = [];
EEG_1.chanlocs(aux_ch, :) = [];

% temp_EEG = EEG_1.data(ssvep_ch, :);

pop_eegplot(EEG_1);
%% band pass filter
EEG_2 = EEG_1;
EEG_2.data= ft_preproc_bandpassfilter(EEG_2.data, EEG_2.srate, freq, 4, 'but');
% EEG_2.data= ft_preproc_bandpassfilter(EEG_2.data, EEG_2.srate, [0.5 55], 4, 'but');

EEG_2.data= ft_preproc_bandstopfilter(EEG_2.data, EEG_2.srate, [59 61], 4, 'brickwall');

pop_eegplot(EEG_2);
figure,
spectopo(EEG_2.data,0,EEG_2.srate);
%% marker
EEG_3 =EEG_2;
temp_events={EEG_3.event.type};
tmep_string= erase(temp_events,'condition ');

for i = 1 : size(tmep_string,2)
    event_marker(i) = str2num(tmep_string{i});
end

temp_latency = [EEG_3.event.latency];

for marker_index =1 : 9
    list_event(:,marker_index)=find(event_marker == marker_index+ 10) ;  %% marker: 11 - 19 
end


for marker_index =1 : 9
    latency(:,marker_index) = temp_latency(list_event(:,marker_index));
    event(:,marker_index) = event_marker(list_event(:,marker_index));  %% trial, class
end
%% data epoching sample
frame =[0 3]; % window, epoch length
trial_data =[];

trial =1;
class =1;

diff(frame) % [0 3 ] > 3-0

start_id = latency(trial, class) + frame(1)*EEG_3.srate ;
end_id = start_id + (frame(2) - frame(1))*EEG_3.srate ;
 

tmp_epoch1p = EEG_3.data(:, start_id:end_id);
trial_data = cat(3,trial_data,tmp_epoch1p);


test_a = 1 ;
test_b = 2 ;

test_a = cat(3,test_a,test_b);



%% data epoching
re_data = [];
frame =[0 3]; % window, epoch length

for class =1: 9  %9
    trial_data =[];
    for trial = 1:5 %5
        start_id = latency(trial, class) + frame(1)*EEG_3.srate ;
        end_id = start_id + (frame(2)- frame(1))*EEG_3.srate;
        tmp_epoch1p = EEG_3.data(:, start_id:end_id);
        
        re_data(:,:,trial,class) = tmp_epoch1p;
%         trial_data = cat(3,trial_data,tmp_epoch1p);
    end
%     re_data = cat(4,re_data,trial_data);  % ch, time, trial, class
end 
%% 

EEG_3.data = awgn(EEG_2.data,10,'measured');
eegplot(EEG_2.data,'data2' ,EEG_3.data);

figure, 
subplot(1,2,1)
spectopo(EEG_2.data,0,EEG_2.srate);
ylim([-40 25]);
subplot(1,2,2)
spectopo(EEG_3.data,0,EEG_2.srate);
ylim([-40 25]);

%%
dat_mat = reshape(re_data,19,[]);

[w,s] = runica(dat_mat,'extended',1);

mixing=w*s*dat_mat;


re_mixed = reshape(mixing,19,901,[]);

eegplot(re_mixed)

