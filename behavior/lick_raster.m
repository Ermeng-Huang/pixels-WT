clear
f1=dir('E:\WT\ser\*.ser');
lick_CR_matrix=cell(0);
lick_Hit_matrix=cell(0);
lick_matrix=cell(0);
for i=1:size(f1,1)
    try
    Data=ser2mat(fullfile('E:\WT\ser\',f1(i,1).name))';
    Sample=Data((Data(:,2)==9&Data(:,3)==1)|(Data(:,2)==10&Data(:,3)==100),:);
    Test=Data((Data(:,2)==9&Data(:,3)==3)|(Data(:,2)==10&Data(:,3)==2),:);
    Response=Data((Data(:,2)==4&Data(:,3)==1)|(Data(:,2)==5&Data(:,3)==1)|(Data(:,2)==6&Data(:,3)==1)|(Data(:,2)==7&Data(:,3)==1),2);    
    Response(Response(:,1)==5|Response(:,1)==7,2)=1; %2-correct,3-WT window,4-delay
    Response(:,3)=0;
    a=40;
    while a<=length(Response)
        goodOff=nnz(Response(a-39:a,:)==5|Response(a-39:a,:)==7);
        if goodOff>=30 %.75 correct rate
            Response(a-39:a,3)=1;
        end
        a=a+1;
    end
    Response(:,4)=(Test(:,1)-Sample(:,1))/1000;
%     Response((Test-Sample)<5000,3)=1; % delay
    
    Sample0=[];    
    Sample0=Sample(all(Response(:,2:4),2),:);   
    Response0=Response(all(Response(:,2:4),2),:);
    for t=1:size(Sample0)
        lick_matrix{end+1,1}=(Data((Data(:,1)>(Sample0(t,1)-3*1000))&(Data(:,1)<(Sample0(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample0(t,1));
        lick_matrix{end,2}=Response0(t,:);
    end
    
    
    if false
        Sample1=[];
        Response1=[];
        Sample1=Sample(Response(:,2)==1&Response(:,1)==5&Response(:,3)==1,:);
        Response1=Response(Response(:,2)==1&Response(:,1)==5&Response(:,3)==1,:);
        %     lick_matrix=zeros(1,17*1000);
        for t=1:size(Sample1)
            lick_CR_matrix{end+1}=(Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1));
        end
        lick=[];
        for t=1:size(Sample1)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick_CR(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
        Sample2=[];
        Response2=[];
        Sample2=Sample(Response(:,2)==1&Response(:,1)==7&Response(:,3)==1,:);
        Response2=Response(Response(:,2)==1&Response(:,1)==7&Response(:,3)==1,:);
        for t=1:size(Sample2)
            lick_Hit_matrix{end+1}=(Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1));
        end
        lick=[];
        for t=1:size(Sample2)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick_Hit(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
    end
    catch
        disp(f1(i,1).name)
        disp(i)
    end
end
% save('lick_raster','lick_matrix')
% save('lick6fromser','lick_CR','lick_Hit','lick_Hit_matrix','lick_CR_matrix')

%% plot lick raster
close all
fh=figure('Color','w','Position',[100,100,400,400]);
hold on
i=1;
while i<=40
    
    if lick_matrix{i+10,2}(:,4)<5 && lick_matrix{i+10,2}(:,1)==7
        fill([-4,14,14,-4]*1000,[40+1-i-0.5,40+1-i-0.5,40+1-i+0.5,40+1-i+0.5],'r','FaceAlpha',0.5,'EdgeColor','none')
    elseif lick_matrix{i+10,2}(:,4)<5 && lick_matrix{i+10,2}(:,1)==5
        fill([-4,14,14,-4]*1000,[40+1-i-0.5,40+1-i-0.5,40+1-i+0.5,40+1-i+0.5],'b','FaceAlpha',0.5,'EdgeColor','none')
    elseif lick_matrix{i+10,2}(:,4)>5 && lick_matrix{i+10,2}(:,1)==7
        fill([-4,14,14,-4]*1000,[40+1-i-0.5,40+1-i-0.5,40+1-i+0.5,40+1-i+0.5],'m','FaceAlpha',0.5,'EdgeColor','none')
    elseif lick_matrix{i+10,2}(:,4)>5 && lick_matrix{i+10,2}(:,1)==5
        fill([-4,14,14,-4]*1000,[40+1-i-0.5,40+1-i-0.5,40+1-i+0.5,40+1-i+0.5],'c','FaceAlpha',0.5,'EdgeColor','none')
    end
    if ~isempty(lick_matrix{i+10,1})
        plot([lick_matrix{i+10,1},lick_matrix{i+10,1}], [40+1-i-0.5 40+1-i+0.5],'k')
    end
    i=i+1;
end
plot([0 0],[0.5,40.5],'k--')
plot([4000 4000],[0.5,40.5],'k--')
plot([7000 7000],[0.5,40.5],'k--')

set(gca,'XTickLabel',0:2:14,'XLim',[-1,11]*1000,'YLim',[0.5,40.5])
xlabel('Time (s)')
ylabel('Trial #')
exportgraphics(fh,'lick_raster.pdf','ContentType','vector');

return

hold on
time=size(lick_Hit,2);
Time = [1:time, fliplr(1:time)];
Highervalue = smooth(mean(lick_Hit,1))' + std(lick_Hit,0,1)/sqrt(size(lick_Hit,1));
Lowervalue = smooth(mean(lick_Hit,1))' - std(lick_Hit,0,1)/sqrt(size(lick_Hit,1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'r', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,smooth(mean(lick_Hit,1))','r-')
hold on
time=size(lick_CR,2);
Time = [1:time, fliplr(1:time)];
Highervalue = smooth(mean(lick_CR,1))' + std(lick_CR,0,1)/sqrt(size(lick_CR,1));
Lowervalue = smooth(mean(lick_CR,1))' - std(lick_CR,0,1)/sqrt(size(lick_CR,1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'b', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,smooth(mean(lick_CR,1))','b-')

plot([15 15],[0 8],'k--')
plot([20 20],[0 8],'k--')
plot([35 35],[0 8],'k--')
plot([40 40],[0 8],'k--')
plot([45 45],[0 8],'k--')
plot([50 50],[0 8],'k--')
plot([55 55],[0 8],'k--')
plot([60 60],[0 8],'k--')
set(gca,'XTick',1:5:85,'XTickLabel',{' ',' ',' ','0',' ',' ',' ',' ',' ','6',' ',' ',' ','','','12'})
xlabel('Time(s)','FontSize',10);
ylabel('Lick(Hz)','FontSize',10);
box off
save('lick3fromser','lick_CR','lick_Hit')
saveas(gcf,'lick3','fig')
% legend('Hit Trials','Correct rejection trials')


clear
load('D:\code\114_sorted_file_path.mat')
addpath('D:\code\fieldtrip-20200320')
ft_defaults
timeLim=[-3,14];
sps=30000;
binSize=0.25;
for i=1:size(sorted_fpath,1)
    rootpath=sorted_fpath{i};
    try
        bnc=h5read(fullfile('D:\neupix\wyt\DataSum',rootpath,'events.hdf5'),'/events');
        lick{i,2}=h5read(fullfile('D:\neupix\wyt\DataSum',rootpath,'events.hdf5'),'/trials')';
    catch
        bnc=h5read(fullfile('D:\neupix\wyt\DataSum\singleProbe',rootpath,'events.hdf5'),'/events');
        lick{i,2}=h5read(fullfile('D:\neupix\wyt\DataSum\singleProbe',rootpath,'events.hdf5'),'/trials')';
    end
    if nnz(bnc(1,:)<0)>0
        disp(rootpath)
    end
    bnc(3,:)=bitand(bnc(2,:),1);
    
    cfg=struct();
    cfg.trl=[lick{i,2}(:,1)+timeLim(1)*sps,lick{i,2}(:,1)+timeLim(2)*sps,zeros(size(lick{i,2},1),1)+timeLim(1)*sps,lick{i,2}];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    FT_SPIKE.label{1}=1;
    FT_SPIKE.timestamp=cell(1,1);
    FT_SPIKE.timestamp{1}=[bnc(1,:)',bnc(3,:)'];
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    for n=1:size(lick{i,2},1)
        lick{i,1}(n,:)=histcounts(FT_SPIKE.time{1}(1,FT_SPIKE.trial{1}==n),FT_SPIKE.trialtime(1,1):binSize:FT_SPIKE.trialtime(1,2))./binSize;
    end
    lick{i,3}=rootpath;
    
end

for i=1:size(lick,1)
    trials=lick{i,2};    
    trials(:,9)=trials(:,5)-trials(:,6)-1;    
    trials(trials(:,9)~=-1,9)=1;
    trials(:,10)=0;
    trials(trials(:,7)==trials(:,9),10)=1;
    trials(:,11)=0;
    a=40;
    while a<=length(trials)
        goodOff=nnz(xor(trials(a-39:a,5)==trials(a-39:a,6),trials(a-39:a,7)>0));
        if goodOff>=30 %.75 correct rate
            trials(a-39:a,11)=1;
        end
        a=a+1;
    end
    
    lick_3{1,1}(i,:)=mean(lick{i,1}(trials(:,8)==3,:),1);
    lick_3{1,2}(i,:)=mean(lick{i,1}(trials(:,8)==3&trials(:,9)==1&trials(:,10)==1&trials(:,11)==1,:),1);
    lick_3{1,3}(i,:)=mean(lick{i,1}(trials(:,8)==3&trials(:,9)==1&trials(:,10)==0&trials(:,11)==1,:),1);
    lick_3{1,4}(i,:)=mean(lick{i,1}(trials(:,8)==3&trials(:,9)==-1&trials(:,10)==0&trials(:,11)==1,:),1);
    lick_3{1,5}(i,:)=mean(lick{i,1}(trials(:,8)==3&trials(:,9)==-1&trials(:,10)==1&trials(:,11)==1,:),1);
    lick_6{1,1}(i,:)=mean(lick{i,1}(trials(:,8)==6,:),1);
    lick_6{1,2}(i,:)=mean(lick{i,1}(trials(:,8)==6&trials(:,9)==1&trials(:,10)==1&trials(:,11)==1,:),1);
    lick_6{1,3}(i,:)=mean(lick{i,1}(trials(:,8)==6&trials(:,9)==1&trials(:,10)==0&trials(:,11)==1,:),1);
    lick_6{1,4}(i,:)=mean(lick{i,1}(trials(:,8)==6&trials(:,9)==-1&trials(:,10)==0&trials(:,11)==1,:),1);
    lick_6{1,5}(i,:)=mean(lick{i,1}(trials(:,8)==6&trials(:,9)==-1&trials(:,10)==1&trials(:,11)==1,:),1);
    clear trials
end

save('D:/code/lick.mat','lick','lick_3','lick_6')

hold on
time=size(lick_6{1,1},2);
Time = [1:time, fliplr(1:time)];
Highervalue = (mean(lick_6{1,2},1)) + std(lick_6{1,2},0,1)/sqrt(size(lick_6{1,2},1));
Lowervalue = (mean(lick_6{1,2},1)) - std(lick_6{1,2},0,1)/sqrt(size(lick_6{1,2},1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'r', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,mean(lick_6{1,2},1),'r-')

% Highervalue = (smooth(mean(lick_6{1,3},1)))' + std(lick_6{1,3},0,1)/sqrt(size(lick_6{1,3},1));
% Lowervalue = (smooth(mean(lick_6{1,3},1)))' - std(lick_6{1,3},0,1)/sqrt(size(lick_6{1,3},1));
% value = [Highervalue, fliplr(Lowervalue)];
% a = fill(Time, value, 'b', 'edgecolor','none');
% alpha(a,0.2);
% plot(1:time,smooth(mean(lick_6{1,3},1))','b-')
% 
% Highervalue = (smooth(mean(lick_6{1,4},1)))' + std(lick_6{1,4},0,1)/sqrt(size(lick_6{1,4},1));
% Lowervalue = (smooth(mean(lick_6{1,4},1)))' - std(lick_6{1,4},0,1)/sqrt(size(lick_6{1,4},1));
% value = [Highervalue, fliplr(Lowervalue)];
% a = fill(Time, value, 'g', 'edgecolor','none');
% alpha(a,0.2);
% plot(1:time,smooth(mean(lick_6{1,4},1))','g-')

Highervalue = (mean(lick_6{1,5},1)) + std(lick_6{1,5},0,1)/sqrt(size(lick_6{1,5},1));
Lowervalue = (mean(lick_6{1,5},1)) - std(lick_6{1,5},0,1)/sqrt(size(lick_6{1,5},1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'b', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,mean(lick_6{1,5},1),'b-')


hold on
time=size(lick_3{1,1},2);
Time = [1:time, fliplr(1:time)];
Highervalue = (mean(lick_3{1,2},1)) + std(lick_3{1,2},0,1)/sqrt(size(lick_3{1,2},1));
Lowervalue = (mean(lick_3{1,2},1)) - std(lick_3{1,2},0,1)/sqrt(size(lick_3{1,2},1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'r', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,mean(lick_3{1,2},1),'r-')

% Highervalue = (smooth(mean(lick_3{1,3},1)))' + std(lick_3{1,3},0,1)/sqrt(size(lick_3{1,3},1));
% Lowervalue = (smooth(mean(lick_3{1,3},1)))' - std(lick_3{1,3},0,1)/sqrt(size(lick_3{1,3},1));
% value = [Highervalue, fliplr(Lowervalue)];
% a = fill(Time, value, 'b', 'edgecolor','none');
% alpha(a,0.2);
% plot(1:time,smooth(mean(lick_3{1,3},1))','b-')
% 
% Highervalue = (smooth(mean(lick_3{1,4},1)))' + std(lick_3{1,4},0,1)/sqrt(size(lick_3{1,4},1));
% Lowervalue = (smooth(mean(lick_3{1,4},1)))' - std(lick_3{1,4},0,1)/sqrt(size(lick_3{1,4},1));
% value = [Highervalue, fliplr(Lowervalue)];
% a = fill(Time, value, 'g', 'edgecolor','none');
% alpha(a,0.2);
% plot(1:time,smooth(mean(lick_3{1,4},1))','g-')

Highervalue = (mean(lick_3{1,5},1)) + std(lick_3{1,5},0,1)/sqrt(size(lick_3{1,5},1));
Lowervalue = (mean(lick_3{1,5},1)) - std(lick_3{1,5},0,1)/sqrt(size(lick_3{1,5},1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'b', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,mean(lick_3{1,5},1),'b-')