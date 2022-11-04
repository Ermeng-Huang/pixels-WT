% extract waveform of certain neuron   
clc
clear

addpath('D:\code\neuropixel-utils')
addpath(genpath('D:\code\Kilosort2'))
channelMapFile='D:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';

filepath='F:\neupix\raw\';
folder='191106-DPA-LearningN_31_gc_imec1_cleaned';
cidx=10273;

fl=dir(fullfile(filepath,folder,'cluster_info.tsv'));


waveform=cell(0,5);
ks=Neuropixel.KiloSortDataset(fl.folder,'channelMap',channelMapFile);
ks.load();

try
    snippetSetTet = ks.getWaveformsFromRawData('cluster_ids',rem(cidx,10000),'num_waveforms', 100, 'best_n_channels', 50, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    snippetSetBest = ks.getWaveformsFromRawData('cluster_ids',rem(cidx,10000),'num_waveforms', 500, 'best_n_channels', 1, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    waveform{end+1,1}=folder;
    waveform{end,2}=cidx;
    waveform{end,3}=snippetSetTet.data;
    waveform{end,4}=mean(snippetSetBest.data,3);
    waveform{end,5}=snippetSetTet.channel_ids_by_cluster;
catch
    errors{end+1}=sprintf('%s %d waveform error',onefile.folder,cidx);
end

savepath=dir(fullfile('F:\neupix\WF\neuropixel','*',folder,'waveform.mat'));

save(fullfile(savepath.folder,sprintf('waveform_%d.mat',cidx)),'waveform');
load(fullfile(savepath.folder,sprintf('waveform_%d.mat',cidx)))
%% plot waveform
channel=waveform{1,5};
waveform=waveform{1,3};
f=figure('Color','w','Position',[100,100,100,1000]);
channel_best=double(channel(1));
channel_plot=channel_best+(-6:1:6);
channel_plot=channel_plot(ismember(channel_plot,channel));
for i=1:9
    subplot(9,1,i)
    plot(mean(waveform(channel==channel_plot(find(channel_plot==channel_best)-5+i),:,:),3),'k-','LineWidth',2)
    set(gca,'XTick',1:30:91,'XTickLabel',{'-1','0','1','2'})
    ylim([-40 25])
    xlabel('time (ms)')
    ylabel('MicroVolt')
    title(sprintf('channel %d',channel(1)+i-5))
    box off
end

exportgraphics(f,fullfile('F:\neupix\WF',sprintf('waveform-%d-%s.pdf',cidx,folder)))