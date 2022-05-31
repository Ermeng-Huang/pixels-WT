clear
f1=dir('E:\WT\ser\*.ser');

[lick3_count_h,lick3_count_f,lick6_count_h,lick6_count_f]=deal(cell(0));
for i=1:size(f1,1)

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
%     Response(:,4)=(Test(:,1)-Sample(:,1))/1000;
    Response((Test(:,1)-Sample(:,1))<5000,4)=1; % delay    
    lick3_count_h{end+1}=arrayfun(@(x)nnz((Data(:,1)>(x+1*1000))&(Data(:,1)<(x+4*1000))&Data(:,2)==0&Data(:,3)==76)...
        ,Sample(Response(:,1)==7&Response(:,3)==1&Response(:,4),1));
    lick3_count_f{end+1}=arrayfun(@(x)nnz((Data(:,1)>(x+1*1000))&(Data(:,1)<(x+4*1000))&Data(:,2)==0&Data(:,3)==76)...
        ,Sample(Response(:,1)==4&Response(:,4),1));
    
    Response(:,4)=0;
    Response((Test(:,1)-Sample(:,1))>5000,4)=1; % delay
    lick6_count_h{end+1}=arrayfun(@(x)nnz((Data(:,1)>(x+1*1000))&(Data(:,1)<(x+7*1000))&Data(:,2)==0&Data(:,3)==76)...
        ,Sample(Response(:,1)==7&Response(:,3)==1&Response(:,4),1));
    lick6_count_f{end+1}=arrayfun(@(x)nnz((Data(:,1)>(x+1*1000))&(Data(:,1)<(x+7*1000))&Data(:,2)==0&Data(:,3)==76)...
        ,Sample(Response(:,1)==4&Response(:,4),1));
end

h_3=cell2mat(lick3_count_h')/3;
f_3=cell2mat(lick3_count_f')/3;
h_6=cell2mat(lick6_count_h')/6;
f_6=cell2mat(lick6_count_f')/6;
%% plot lick raster

m3_h=mean(h_3);
m3_f=mean(f_3);
s3_h=std(h_3)/sqrt(numel(h_3));
s3_f=std(f_3)/sqrt(numel(f_3));
m6_h=mean(h_6);
m6_f=mean(f_6);
s6_h=std(h_6)/sqrt(numel(h_6));
s6_f=std(f_6)/sqrt(numel(f_6));

fh=figure('Color','w','Position',[100,100,250,300]);
hold on
errorbar(1:4,[m3_h,m3_f,m6_h,m6_f],[s3_h,s3_f,s6_h,s6_f],'ko')
set(gca,'Xlim',[0.5,4.5],'XTick',1:4,'XTickLabel',{'3s hit trial','3s false trial','6s hit trial','6s false trial'},'XTickLabelRotation',30,'Ylim',[0,0.2])
p=anovan([h_3;f_3;h_6;f_6],{[zeros(numel(h_3),1);ones(numel(f_3),1);zeros(numel(h_6),1);ones(numel(f_6),1)]...
    ,[zeros(numel([h_3;f_3]),1);ones(numel([h_6;f_6]),1)]});
text(0.8,0.2,sprintf('p(correct vs error)=%0.3f',p(1)))
text(0.8,0.18,sprintf('p(delay)=%0.3f',p(2)))
ylabel('average lick frequency (Hz)')
exportgraphics(fh,'lick_freq_correctVSerror.pdf','ContentType','vector');

