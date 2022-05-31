clear
f1=dir('E:\WT\ser\*.ser');

[RT3_h,RT6_h,RTswitch_h,RTswitch2_h]=deal(cell(0));
for i=1:size(f1,1)

    Data=ser2mat(fullfile('E:\WT\ser\',f1(i,1).name))';
    Sample=Data((Data(:,2)==9&Data(:,3)==1)|(Data(:,2)==10&Data(:,3)==100),:);
    Test=Data((Data(:,2)==9&Data(:,3)==3)|(Data(:,2)==10&Data(:,3)==2),:);
    Response=Data((Data(:,2)==4&Data(:,3)==1)|(Data(:,2)==5&Data(:,3)==1)|(Data(:,2)==6&Data(:,3)==1)|(Data(:,2)==7&Data(:,3)==1),2);    
    Response(Response(:,1)==5|Response(:,1)==7,2)=1; %2-correct,3-WT window,4-delay
    Response(:,3)=0;
    Delay=round((Test(:,1)-Sample(:,1))/1000)-1;
    a=40;
    while a<=length(Response)
        goodOff=nnz(Response(a-39:a,:)==5|Response(a-39:a,:)==7);
        if goodOff>=30 %.75 correct rate
            Response(a-39:a,3)=1;
        end
        a=a+1;
    end
    
    %3s
    try
    RT3_h{end+1}=arrayfun(@(x)Data(find((Data(:,1)>x)&(Data(:,1)<(x+3*1000))&Data(:,2)==0&Data(:,3)==76,1),1)...       
        ,Test(Response(:,1)==7&Response(:,3)==1&Delay==3,1)) -Test(Response(:,1)==7&Response(:,3)==1&Delay==3,1);
    catch
    end
    %6s    
    try
    RT6_h{end+1}=arrayfun(@(x)Data(find((Data(:,1)>x)&(Data(:,1)<(x+3*1000))&Data(:,2)==0&Data(:,3)==76,1),1)...       
        ,Test(Response(:,1)==7&Response(:,3)==1&Delay==6,1)) -Test(Response(:,1)==7&Response(:,3)==1&Delay==6,1);
        catch
    end
    %switch trial
    try
    RTswitch_h{end+1}=arrayfun(@(x)Data(find((Data(:,1)>x)&(Data(:,1)<(x+3*1000))&Data(:,2)==0&Data(:,3)==76,1),1)...       
        ,Test(Response(:,1)==7&Response(:,3)==1&cat(1,0,diff(Delay))==3,1)) -Test(Response(:,1)==7&Response(:,3)==1&cat(1,0,diff(Delay))==3,1);
        catch
    end
    
    try
    RTswitch2_h{end+1}=arrayfun(@(x)Data(find((Data(:,1)>x)&(Data(:,1)<(x+3*1000))&Data(:,2)==0&Data(:,3)==76,1),1)...       
        ,Test(Response(:,1)==7&Response(:,3)==1&cat(1,0,diff(Delay))==-3,1)) -Test(Response(:,1)==7&Response(:,3)==1&cat(1,0,diff(Delay))==-3,1);
        catch
    end
end

h_3=double(cell2mat(RT3_h'))/1000;
h_6=double(cell2mat(RT6_h'))/1000;
h_switch=double(cell2mat(RTswitch_h'))/1000;
h_switch2=double(cell2mat(RTswitch2_h'))/1000;

%% plot lick raster

m3_h=mean(h_3);
s3_h=std(h_3)/sqrt(numel(h_3));
m6_h=mean(h_6);
s6_h=std(h_6)/sqrt(numel(h_6));
m6_switch=mean(h_switch);
s6_switch=std(h_switch)/sqrt(numel(h_switch));
m6_switch2=mean(h_switch2);
s6_switch2=std(h_switch2)/sqrt(numel(h_switch2));

fh=figure('Color','w','Position',[100,100,250,300]);
hold on
bar(1:4,[m3_h,m6_h,m6_switch,m6_switch2],'w')
errorbar(1:4,[m3_h,m6_h,m6_switch,m6_switch2],[s3_h,s6_h,s6_switch,s6_switch2],'k.')
set(gca,'Xlim',[0.5,4.5],'XTick',1:4,'XTickLabel',{'3s hit trial','6s hit trial','switched trial (3->6)','switched trial (6->3)'},'XTickLabelRotation',30,'Ylim',[0,1])
ylabel('Reaction Time (s)')
exportgraphics(fh,'ReactionTime.pdf','ContentType','vector');

anovan([h_3;h_6;h_switch;h_switch2],{[zeros(numel(h_3),1);ones(numel(h_6),1);2*ones(numel(h_switch),1);3*ones(numel(h_switch2),1)]})
anovan([h_3;h_6],{[zeros(numel(h_3),1);ones(numel(h_6),1)]})
anovan([h_6;h_switch],{[ones(numel(h_6),1);2*ones(numel(h_switch),1)],})

text(0.8,0.2,sprintf('p(correct vs error)=%0.3f',p(1)))
text(0.8,0.18,sprintf('p(delay)=%0.3f',p(2)))

