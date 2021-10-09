clear
clc
close all
level=5;
rootpath=fullfile('F:','neupix','STP','coactivate','0624','20ms');
fl.congru=dir(fullfile(rootpath,'0315_stp_congru_*.mat'));
fl.incongru=dir(fullfile(rootpath,'0315_stp_incongru_*.mat'));
fl.nonmem=dir(fullfile(rootpath,'0315_stp_non-mem_*.mat'));
% fl.both=dir(fullfile(rootpath,'0315_stp_incongru-congru_*.mat')); %%

memtypes=convertCharsToStrings(fieldnames(fl))';

stats=struct();
proability=cell(0);diff_reg=cell(0);same_reg=cell(0);proability_shuffled=cell(0);high_reg=cell(0);low_reg=cell(0);
IsProgressive=cell(0);IsRegressive=cell(0);IsComplex=cell(0);is_OBM1=cell(0);
for memtype=memtypes
    stats.(memtype)=struct();
    stats.(memtype).postspk=[];
    stats.(memtype).postspk_shuffled=[];
    stats.(memtype).sess=[];
    stats.(memtype).sess_suids=[];
    for fidx=1:size(fl.(memtype))
        fstr=load(fullfile(fl.(memtype)(fidx).folder,fl.(memtype)(fidx).name));
        stats.(memtype).postspk=[stats.(memtype).postspk;cellfun(@(x)x,fstr.postspk,'UniformOutput',false)'];
%         stats.(memtype).postspk_shuffled=[stats.(memtype).postspk_shuffled;...
%             cellfun(@(x)x.shuffle100ms,fstr.postspk,'UniformOutput',false)'];        
        stats.(memtype).sess=[stats.(memtype).sess;repmat(unique(fstr.suid_sess),size(fstr.suid,1),1)];
        stats.(memtype).sess_suids = [stats.(memtype).sess_suids;fstr.suid];
    end
    stats.(memtype).reg=bz.hist.tag_hist_reg(stats.(memtype),'type','neupix');
    [is_diff,is_same]=diff_at_level(stats.(memtype).reg,level);
    [is_high,is_low]=high_low(stats.(memtype).reg,level);
%     [IsProgressive{end+1},IsRegressive{end+1},IsComplex{end+1}]=Progressive(stats.(memtype));
    is_OBM1{end+1}=OBM1(stats.(memtype).reg,level);
    diff_reg{end+1}=all(is_diff,2);
    same_reg{end+1}=all(is_same,2);
    high_reg{end+1}=all(is_high,2);
    low_reg{end+1}=all(is_low,2);
    
    proability{end+1}=100*(cell2mat(cellfun(@(x)(x(4,:)),stats.(memtype).postspk,'UniformOutput',false)));
%     proability_shuffled{end+1}=100*(cell2mat(cellfun(@(x)(x(4,:)),stats.(memtype).postspk_shuffled,'UniformOutput',false)));  
%     proability{end+1}=100*(cell2mat(cellfun(@(x)(x(4,:)),stats.(memtype).postspk,'UniformOutput',false))...
%         -cell2mat(cellfun(@(x)(x(4,:)),stats.(memtype).postspk_shuffled,'UniformOutput',false)));  
%     proability{end+1}=cell2mat(cellfun(@(x)(x(4,:)-x(2,:).*x(3,:)),stats.(memtype).postspk,'UniformOutput',false));    
%     proability{end+1}=stats.(memtype).postspk(:,4)- (stats.(memtype).postspk(:,2).*stats.(memtype).postspk(:,2));
  
end

rootpath=fullfile('F:','neupix','STP','coactivate','0624','100ms');
fl.congru=dir(fullfile(rootpath,'0315_stp_congru_*.mat'));
fl.incongru=dir(fullfile(rootpath,'0315_stp_incongru_*.mat'));
fl.nonmem=dir(fullfile(rootpath,'0315_stp_non-mem_*.mat'));
% fl.both=dir(fullfile(rootpath,'0315_stp_incongru-congru_*.mat')); %%

memtypes=convertCharsToStrings(fieldnames(fl))';
stats=struct();
for memtype=memtypes
    stats.(memtype)=struct();
    stats.(memtype).postspk=[];
    stats.(memtype).postspk_shuffled=[];
    stats.(memtype).sess=[];
    stats.(memtype).sess_suids=[];
    for fidx=1:size(fl.(memtype))
        fstr=load(fullfile(fl.(memtype)(fidx).folder,fl.(memtype)(fidx).name));
        stats.(memtype).postspk_shuffled=[stats.(memtype).postspk_shuffled;cellfun(@(x)x,fstr.postspk_shuffled,'UniformOutput',false)'];
%         stats.(memtype).postspk_shuffled=[stats.(memtype).postspk_shuffled;...
%             cellfun(@(x)x.shuffle100ms,fstr.postspk,'UniformOutput',false)'];        
        stats.(memtype).sess=[stats.(memtype).sess;repmat(unique(fstr.suid_sess),size(fstr.suid,1),1)];
        stats.(memtype).sess_suids = [stats.(memtype).sess_suids;fstr.suid];
    end
    
    proability_shuffled{end+1}=100*(cell2mat(cellfun(@(x)(x(4,:)),stats.(memtype).postspk_shuffled,'UniformOutput',false)));

end
for i=1:3
    bad=any([isnan(proability{i}),isnan(proability_shuffled{i})],2);
    diff_reg{i}(bad,:)=[];
    same_reg{i}(bad,:)=[];
    high_reg{i}(bad,:)=[];
    low_reg{i}(bad,:)=[];
    proability{i}(bad,:)=[];
    proability_shuffled{i}(bad,:)=[];
%     IsProgressive{i}(bad,:)=[];
%     IsRegressive{i}(bad,:)=[];
%     IsComplex{i}(bad,:)=[];
    is_OBM1{i}(bad,:)=[];
    

end
%% plot 
c={'r','b','k','m'};
t={'congru-cross region','congru-within region','incongru-cross region','incongru-within region','nonmem-cross region','nonmem-within region'};
L={'-'};
fh=figure('Color','w','Position',[100,100,400,300]);
 
subplot(1,2,1) 
hold on
for i=1:3
    nnz(diff_reg{i}) 
    plotLine(proability{i}(diff_reg{i}&is_OBM1{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i},:),c{i},L{1});
%     plotLine(cell2mat(cellfun(@(x)x(2,:),proability{i}(diff_reg{i},:),'UniformOutput',false)),c{i},L{2});      
end
title('cross region')
ylim([0 5])

set(gca,'XTick',1:2.5:10,'XTickLabel',{'-200' '-150','-100','-50',' ','100','120','140','160','180','200'},'XTickLabelRotation',45)
ylabel('Shuffled-reducted post cell spike increase (%)')
xlabel('time lag (ms)')
grid on
subplot(1,2,2) 
hold on
for i=1:3        
    plotLine(proability{i}(same_reg{i}&is_OBM1{i},:)-proability_shuffled{i}(same_reg{i}&is_OBM1{i},:),c{i},L{1});
%     plotLine(cell2mat(cellfun(@(x)x(2,:),proability{i}(same_reg{i},:),'UniformOutput',false)),c{i},L{2});      
end
title('within region')
ylim([0 5])
set(gca,'XTick',1:2.5:10,'XTickLabel',{'-200' '-150','-100','-50',' ','100','120','140','160','180','200'},'XTickLabelRotation',45)
ylabel('Shuffled-reducted post cell spike increase (%)')
xlabel('time lag (ms)')
grid on


exportgraphics(fh,fullfile(rootpath,'coactiavtion_shuffledtime.pdf'),'ContentType','vector')

%%% P-value
x=[];g1=[];g2=[];
for i=1:3
    x=[x;proability{i}(diff_reg{i}&is_OBM1{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i},:)];
    g1=[g1;i*ones(size(proability{i}(diff_reg{i}&is_OBM1{i},:)))];
    g2=[g2;repmat([1:10],size(proability{i}(diff_reg{i}&is_OBM1{i},:),1),1)];   
end
        
p=anovan(reshape(x,1,[]),{reshape(g1,1,[]),reshape(g2,1,[])},'display','off');

x=[];g1=[];g2=[];
for i=1:3
    x=[x;proability{i}(same_reg{i}&is_OBM1{i},:)-proability_shuffled{i}(same_reg{i}&is_OBM1{i},:)];
    g1=[g1;i*ones(size(proability{i}(same_reg{i}&is_OBM1{i},:)))];
    g2=[g2;repmat([1:10],size(proability{i}(same_reg{i}&is_OBM1{i},:),1),1)];   
end
        
p=anovan(reshape(x,1,[]),{reshape(g1,1,[]),reshape(g2,1,[])},'display','off');


%% plot in different hierarchy
c={'r','b','k','m'};
L={'-','-'};
t={'congru-cross region','congru-within region','incongru-cross region','incongru-within region','nonmem-cross region','nonmem-within region'};
fh=figure('Color','w','Position',[100,100,300,300]);

hold on
for i=1    
%     nnz(diff_reg{i}&low_reg{i})
    plotLine(proability{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:),c{i},'-');  
    plotLine(proability{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:),c{i},'--');
%     plotLine(cell2mat(cellfun(@(x)x(2,:),proability{i}(diff_reg{i},:),'UniformOutput',false)),c{i},L{2});      
end
title('cross-region')
ylim([-1 5])
set(gca,'XTick',1:2.5:10,'XTickLabel',{'-200' '-150','-100','-50',' ','100','120','140','160','180','200'},'XTickLabelRotation',45)
ylabel('Shuffled-reducted pre-cell coactivation increase (%)')
xlabel('time lag (ms)')
exportgraphics(fh,fullfile(rootpath,'coactiavtion_h.pdf'),'ContentType','vector')

% %%% P-value
x=[];g1=[];g2=[];
i=1;
x=[proability{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:)...
    ;proability{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:)-proability_shuffled{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:)];
g1=[1*ones(size(proability{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:)))...
    ;2*ones(size(proability{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:)))];
g2=[repmat([1:10],size(proability{i}(diff_reg{i}&is_OBM1{i}&high_reg{i},:),1),1)...
    ;repmat([1:10],size(proability{i}(diff_reg{i}&is_OBM1{i}&low_reg{i},:),1),1)];
p=anovan(reshape(x,1,[]),{reshape(g1,1,[]),reshape(g2,1,[])},'display','off');

%%
c={'r','b','k','m'};
t={'cross region','within region','incongru-cross region','incongru-within region'};
fh=figure('Color','w','Position',[100,100,600,300]);
subplot(1,2,1)
hold on
for i=1:2
    plotLine(proability{i}(diff_reg{i}&IsProgressive{i},:)-proability_shuffled{i}(diff_reg{i}&IsProgressive{i},:)...
        ,c{i},'-');  
    plotLine(proability{i}(diff_reg{i}&IsRegressive{i},:)-proability_shuffled{i}(diff_reg{i}&IsRegressive{i},:)...
        ,c{i},'--'); 
    plotLine(proability{i}(diff_reg{i}&IsComplex{i},:)-proability_shuffled{i}(diff_reg{i}&IsComplex{i},:)...
        ,c{i},':'); 
end
legend('Progressive-congru','Regressive-congru','Complex-congru','Progressive-incongru','Regressive-incongru','Complex-incongru')
title(t{1})
ylim([-1 5])
set(gca,'XTick',1:2.5:10,'XTickLabel',{'-2000' '-1500','-1000','-500',' ','1000','1200','1400','1600','1800','2000'},'XTickLabelRotation',45)
ylabel('Shuffled-reducted post cell spike increase (%)')

subplot(1,2,2)
hold on
for i=1:2
    plotLine(proability{i}(same_reg{i}&IsProgressive{i},:)-proability_shuffled{i}(same_reg{i}&IsProgressive{i},:)...
        ,c{i},'-');  
    plotLine(proability{i}(same_reg{i}&IsRegressive{i},:)-proability_shuffled{i}(same_reg{i}&IsRegressive{i},:)...
        ,c{i},'--'); 
    plotLine(proability{i}(same_reg{i}&IsComplex{i},:)-proability_shuffled{i}(same_reg{i}&IsComplex{i},:)...
        ,c{i},':'); 
end
legend('Progressive-congru','Regressive-congru','Complex-congru','Progressive-incongru','Regressive-incongru','Complex-incongru')
title(t{2})
ylim([-1 5])
set(gca,'XTick',1:2.5:10,'XTickLabel',{'-2000' '-1500','-1000','-500',' ','1000','1200','1400','1600','1800','2000'},'XTickLabelRotation',45)
ylabel('Shuffled-reducted post cell spike increase (%)')


fh=figure('Color','w','Position',[100,100,450,300]);
subplot(1,2,1)
hold on
plotLine(cell2mat(arrayfun(@(i)proability{i}(diff_reg{i}&IsProgressive{i},:)-proability_shuffled{i}(diff_reg{i}&IsProgressive{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},'-');
plotLine(cell2mat(arrayfun(@(i)proability{i}(diff_reg{i}&IsRegressive{i},:)-proability_shuffled{i}(diff_reg{i}&IsRegressive{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},'--');
plotLine(cell2mat(arrayfun(@(i)proability{i}(diff_reg{i}&IsComplex{i},:)-proability_shuffled{i}(diff_reg{i}&IsComplex{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},':');
title(t{1})
subplot(1,2,2)
hold on
plotLine(cell2mat(arrayfun(@(i)proability{i}(same_reg{i}&IsProgressive{i},:)-proability_shuffled{i}(same_reg{i}&IsProgressive{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},'-');
plotLine(cell2mat(arrayfun(@(i)proability{i}(same_reg{i}&IsRegressive{i},:)-proability_shuffled{i}(same_reg{i}&IsRegressive{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},'--');
plotLine(cell2mat(arrayfun(@(i)proability{i}(same_reg{i}&IsComplex{i},:)-proability_shuffled{i}(same_reg{i}&IsComplex{i},:)...
    ,[1,2],'UniformOutput',false)'),c{1},':');
title(t{2})


%%
function [is_diff,is_same]=diff_at_level(reg,level)
is_diff=false(numel(reg),1);
is_same=false(numel(reg),1);
for i=1:numel(reg)   
    for treedep=level
        if all(cellfun(@(x)ismember(x{2},{'CTX'}),reg{i})) && all(cellfun(@(x)~isempty(x{treedep}),reg{i}))           
            if numel(unique(cellfun(@(x)x{treedep},reg{i},'UniformOutput',false)))==1
                is_same(i,1)=true;
            elseif numel(unique(cellfun(@(x)x{treedep},reg{i},'UniformOutput',false)))==numel(reg{i})
                is_diff(i,1)=true;
            end            
        end
    end
end
end
function is_OBM1=OBM1(reg,level)
is_OBM1=false(numel(reg),1);
load('D:\code\pixels-master\jpsth\+bz\+hist\OBM1Map.mat')
for i=1:numel(reg)   
    for treedep=level
        if all(cellfun(@(x)ismember(x{2},{'CTX'}),reg{i})) && all(cellfun(@(x)~isempty(x{treedep}),reg{i}))  
            try
                if nnz(cellfun(@(x)OBM1map(x{treedep}),reg{i}))==3
                    is_OBM1(i,1)=true;
                end
            catch
                is_OBM1(i,1)=false;
            end
                
        end
    end
end
end

function [is_loop]=IsLoop(in)
load('D:\code\pixels-master\jpsth\+bz\+hist\rings_bz.mat')
is_loop=false(numel(in.sess),3);
for i=1:numel(in.sess)    
    for rsize=3:5       
        if sum((rings{in.sess(i),rsize-2}==in.sess_suids(i,1))+2*(rings{in.sess(i),rsize-2}==in.sess_suids(i,2))...
                +3*(rings{in.sess(i),rsize-2}==in.sess_suids(i,3)),2)==6         
           disp((rings{in.sess(i),rsize-2}==in.sess_suids(i,1))+2*(rings{in.sess(i),rsize-2}==in.sess_suids(i,2))...
                +3*(rings{in.sess(i),rsize-2}==in.sess_suids(i,3)))      
            
        end
    end
end
end

function [is_Progressive,is_Regressive,is_Complex]=Progressive(in)
load('D:\code\pixels-master\jpsth\+bz\+hist\commap.mat')

is_Progressive=false(numel(in.sess),1);
is_Regressive=false(numel(in.sess),1);
is_Complex=false(numel(in.sess),1);
for i=1:numel(in.sess)  
    for cellID=1:3
        try
            try
                com(cellID)=commap.(sprintf('s%d',in.sess(i))).s1(in.sess_suids(i,cellID));
            catch
                com(cellID)=commap.(sprintf('s%d',in.sess(i))).s2(in.sess_suids(i,cellID));
            end
        catch
            com(cellID)=-100;
        end
    end
    if nnz(com==-100)==0
        if com(1)<com(3) && com(2)<com(3)
            is_Progressive(i)=true;
        elseif com(1)>com(3) && com(2)>com(3)
            is_Regressive(i)=true;
        else
            is_Complex(i)=true;
        end
    end
end
end

function plotLine(Data,c,l)
    m=mean(Data,1);
    plot(1:size(m,2),m,c,'LineStyle',l);
%     s=std(Data)/sqrt(size(Data,1));
%     fill([1:size(m,2),fliplr(1:size(m,2))],[m+s,fliplr(m-s)],c,'EdgeColor','none','FaceAlpha',0.2);
    ci=bootci(1000,@(x) mean(x),Data);
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end

function [is_high,is_low]=high_low(reg,level)
load('D:/code/ratiomap.mat')
is_high=false(numel(reg),1);
is_low=false(numel(reg),1);
for i=1:numel(reg)
    try 
    if all(cellfun(@(x)ismember(x{2},{'CTX'}),reg{i})) && all(cellfun(@(x)~isempty(x{level}),reg{i}))    
        ratioCompare=cellfun(@(x)ratiomap(x{5}),reg{i}(1:end-1))<ratiomap(reg{i}{end}{5});
        if all(ratioCompare==1)
            is_high(i,1)=true; %high--low
        elseif ~any(ratioCompare)
            is_low(i,1)=true;
        end
    end
    catch
    end
end
end
%%