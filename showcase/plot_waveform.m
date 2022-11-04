clear
% filename='324_M34_20191102_g1_imec0_cleaned';
filename='10562_M43_20191203_g0_imec1_cleaned';
homedir='I:\WT\WF\neuropixel';
% id=str2double(regexp(filename,'(\w*)(?=_M)','match','once'));
% folder=regexp(filename,'(?<=_)(\w*)','match','once');
id=10521;
folder='191023-DPA-Learning10_30_g0_imec0_cleaned';
f=dir(fullfile(homedir,'*',folder,'waveform_*.mat'));
fstr=load(fullfile(f.folder,f.name));
 
waveform=fstr.waveform{cell2mat(fstr.waveform(:,2))==id,3};
channel=fstr.waveform{cell2mat(fstr.waveform(:,2))==id,5};
f=figure('Color','w','Position',[100,100,100,1000]);
for i=1:9
    subplot(9,1,10-i)
    plot(mean(waveform(channel==(channel(1)+i-5),:,:),3),'k-','LineWidth',2)
    set(gca,'XTick',1:30:91,'XTickLabel',{'-1','0','1','2'})
    ylim([-40 25])
    xlabel('time (ms)')
    ylabel('MicroVolt')
    title(sprintf('channel %d',channel(1)+i-5))
    box off
end
exportgraphics(f,sprintf('waveform-%s.pdf',filename))
Footer
© 2022 GitHub, Inc.
Footer navigation

    Terms
    Privacy
    Security
    Status
    Docs
    Contact GitHub
    Pricing
    API
    Training
    Blog
    About

