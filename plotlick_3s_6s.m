%% 同时画3s和6s在hit和CR trial里面的lick rate
clear
f1=dir('F:\neupix-old\ser\*.ser');
lick_CR_matrix=cell(0);
lick_Hit_matrix=cell(0);
lick_matrix=cell(0);
for i=1:size(f1,1)
    try
        Data=ser2mat(fullfile('F:\neupix-old\ser\',f1(i,1).name))';
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
        Response(:,4)=(Test(:,1)-Sample(:,1))/1000; % delay
        
%         Sample0=[];
%         Sample0=Sample(all(Response(:,2:3),2),:);
%         Response0=Response(all(Response(:,2:3),2),:);
        
        %     for t=1:size(Sample0)
        %         lick_matrix{end+1,1}=(Data((Data(:,1)>(Sample0(t,1)-3*1000))&(Data(:,1)<(Sample0(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample0(t,1));
        %         lick_matrix{end,2}=Response0(t,:);
        %     end
        
        % 3s-CR
        Sample1=Sample(all(Response(:,2:3)==1,2)&Response(:,1)==5&Response(:,4)==4,:);
        Response1=Response(all(Response(:,2:3)==1,2)&Response(:,1)==5&Response(:,4)==4,:);
        lick=[];
        for t=1:size(Sample1)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick3_CR(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
        
        % 3s-Hit
        Sample2=Sample(all(Response(:,2:3)==1,2)&Response(:,1)==7&Response(:,4)==4,:);
        Response2=Response(all(Response(:,2:3)==1,2)&Response(:,1)==7&Response(:,4)==4,:);
        lick=[];
        for t=1:size(Sample2)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick3_Hit(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
        
        % 6s-CR
        Sample1=Sample(all(Response(:,2:3)==1,2)&Response(:,1)==5&Response(:,4)==7,:);
        Response1=Response(all(Response(:,2:3)==1,2)&Response(:,1)==5&Response(:,4)==7,:);
        lick=[];
        for t=1:size(Sample1)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick6_CR(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
        
        % 6s-Hit
        Sample2=Sample(all(Response(:,2:3)==1,2)&Response(:,1)==7&Response(:,4)==7,:);
        Response2=Response(all(Response(:,2:3)==1,2)&Response(:,1)==7&Response(:,4)==7,:);
        lick=[];
        for t=1:size(Sample2)
            lick(t,:)=histcounts((Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1))/1000,-3:0.2:14)/0.2;
        end
        lick6_Hit(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);  
        
    catch
        disp(f1(i,1).name)
        disp(i)
    end
end


%%
fh=figure('Color','w','Position',[100,100,500,500]);
hold on
plotLick(lick3_Hit,'r')
plotLick(lick3_CR,'b')
plotLick(lick6_Hit,'m')
plotLick(lick6_CR,'c')


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
exportgraphics(fh,'lick_rate.pdf','ContentType','vector');
%% function
function plotLick(in,color)
time=size(in,2);
Time = [1:time, fliplr(1:time)];
Highervalue = smooth(mean(in,1))' + std(in,0,1)/sqrt(size(in,1));
Lowervalue = smooth(mean(in,1))' - std(in,0,1)/sqrt(size(in,1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value,color, 'edgecolor','none');
alpha(a,0.2);
plot(1:time,smooth(mean(in,1))','-','color',color)
end
