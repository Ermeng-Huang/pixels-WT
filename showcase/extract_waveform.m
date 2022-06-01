% extract waveform of certain neuron   
clc
clear

addpath('D:\code\neuropixel-utils')
addpath(genpath('D:\code\Kilosort2'))
channelMapFile='D:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';

filepath='/neuropixel/neupix/Raw_Data';
folder='191023-DPA-Learning10_30_g0_imec0_cleaned';
cidx=10521;

fl=dir(fullfile(filepath,'*',folder,'cluster_info.tsv'));

waveform=cell(0,5);
ks=Neuropixel.KiloSortDataset(fl.folder,'channelMap',channelMapFile);
ks.load();

try
    snippetSetTet = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 100, 'best_n_channels', 50, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    snippetSetBest = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 500, 'best_n_channels', 1, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    waveform{end+1,1}=rootpath;
    waveform{end,2}=cluster_ids(cidx);
    waveform{end,3}=snippetSetTet.data;
    waveform{end,4}=mean(snippetSetBest.data,3);
    waveform{end,5}=snippetSetTet.channel_ids_by_cluster;
catch
    errors{end+1}=sprintf('%s %d waveform error',onefile.folder,cidx);
end


save(fullfile(fl.folder,sprintf('waveform_%d.mat',cidx)),'waveform');
