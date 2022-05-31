function sample_correct_error(type,opt)
arguments
    type (1,:) char {mustBeMember(type,{'Ordinal(From-3Sample)','Ordinal(From-1Test)','Ordinal(Combined)','Duration-coding'})} = 'Ordinal(Combined)'
    opt.single_bin (1,1) logical = true
    opt.bin (1,:) char {mustBeMember( opt.bin,{'first1','late1'})} = 'first1'
    opt.plot_scatter (1,1) logical = false
    opt.plot_per_su (1,1) logical = true
    opt.plot_per_trial (1,1) logical = false
    opt.ctx_only (1,1) logical = false
    opt.plot_showcase (1,1) logical = false
end
set(groot,'defaultTextFontSize',10);
colors={'r','b'};
phase=["sample","early_delay","late_delay"];
phase_bin=[3,8;0,1;1,2.5;2.5,4;4,5;8,13];
sample=load_bin_sel('sample');
homedir=ephys.util.getHomedir('type','raw');
stats=struct();

for i=1:numel(phase)
    aucc_perneuron=[];
    auce_perneuron=[];
    typesel=find(sample(:,i));
    idx=1;
    stats.type1=nan(0,4);    
    stats.type2=nan(0,4);
    for ii=typesel'        
        if strcmp(phase(i),"presample")||strcmp(phase(i),"postreward")
            fpath=fullfile('F:\neupix\DataSum',meta_6.allpath{ii},'FR_All_ITI_250_-3_13.hdf5');  
        else
            fpath=fullfile('F:\neupix\DataSum',meta_6.allpath{ii},'FR_All_250.hdf5');
        end
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr=h5read(fpath,'/FR_All');
        if nnz(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==3)==0 ...
                ||nnz(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==3)==0 ...
                ||nnz(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6)==0 ...
                ||nnz(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6)==0
            continue
        end
        bins=(phase_bin(i,1)*4+1:phase_bin(i,2)*4)+12;
        
        if strcmp(phase(i),"late_delay")||strcmp(phase(i),"test")
            cs1=mean(fr(trials(:,9)~=0&trials(:,10)~=0&trials(:,5)==4&trials(:,8)==3, suid==meta_6.allcid(ii), bins),3);
            cs2=mean(fr(trials(:,9)~=0&trials(:,10)~=0&trials(:,5)==8&trials(:,8)==3, suid==meta_6.allcid(ii), bins),3);
            es1=mean(fr(trials(:,10)==0&trials(:,5)==4&trials(:,8)==3, suid==meta_6.allcid(ii), bins),3);
            es2=mean(fr(trials(:,10)==0&trials(:,5)==8&trials(:,8)==3, suid==meta_6.allcid(ii), bins),3);
            cs1=[cs1;mean(fr(trials(:,9)~=0&trials(:,10)~=0&trials(:,5)==4&trials(:,8)==6, suid==meta_6.allcid(ii), bins+12),3)];
            cs2=[cs2;mean(fr(trials(:,9)~=0&trials(:,10)~=0&trials(:,5)==8&trials(:,8)==6, suid==meta_6.allcid(ii), bins+12),3)];
            es1=[es1;mean(fr(trials(:,10)==0&trials(:,5)==4&trials(:,8)==6, suid==meta_6.allcid(ii), bins+12),3)];
            es2=[es2;mean(fr(trials(:,10)==0&trials(:,5)==8&trials(:,8)==6, suid==meta_6.allcid(ii), bins+12),3)];
        else
            cs1=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 ,suid==meta_6.allcid(ii),bins),3);
            cs2=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 ,suid==meta_6.allcid(ii),bins),3);
            es1=mean(fr(trials(:,10)==0 &trials(:,5)==4 , suid==meta_6.allcid(ii),bins),3);
            es2=mean(fr(trials(:,10)==0 &trials(:,5)==8 ,suid==meta_6.allcid(ii),bins),3);
        end
        pref=(mean(cs1)- mean(cs2))>0;
        
        delaymm=mean([mean(cs1),mean(cs2)]);
        delaystd=std([cs1;cs2]);
        if delaystd==0, continue;  end
        
        cs1=(cs1-delaymm)./delaystd;
        cs2=(cs2-delaymm)./delaystd;
        es1=(es1-delaymm)./delaystd;
        es2=(es2-delaymm)./delaystd;  
        if pref
            stats.type1(end+1,:)=[mean(cs1),mean(cs2),mean(es1),mean(es2)];
            [~,~,~,aucc_perneuron(end+1,:)]=perfcurve([zeros(numel(cs1),1);ones(numel(cs2),1)],[cs1;cs2],0);
            [~,~,~,auce_perneuron(end+1,:)]=perfcurve([zeros(numel(es1),1);ones(numel(es2),1)],[es1;es2],0);
        else
            stats.type2(end+1,:)=[mean(cs1),mean(cs2),mean(es1),mean(es2)];
            [~,~,~,aucc_perneuron(end+1,:)]=perfcurve([zeros(numel(cs2),1);ones(numel(cs1),1)],[cs2;cs1],0);
            [~,~,~,auce_perneuron(end+1,:)]=perfcurve([zeros(numel(es2),1);ones(numel(es1),1)],[es2;es1],0);
        end
        idx=idx+1;
        
    end
    %auc
    onecolumn=size(stats.type1,1)+size(stats.type2,1);
    [~,~,~,aucc(i)]=perfcurve((1:2*onecolumn)>onecolumn,[stats.type1(:,1);stats.type2(:,2);stats.type1(:,2);stats.type2(:,1)],0);
    [~,~,~,auce(i)]=perfcurve((1:2*onecolumn)>onecolumn,[stats.type1(:,3);stats.type2(:,4);stats.type1(:,4);stats.type2(:,3)],0);
    
    aucc_m(i)=mean(aucc_perneuron);
    aucc_std=std(aucc_perneuron)/sqrt(numel(aucc_perneuron));
    auce_m(i)=mean(auce_perneuron);
    auce_std=std(auce_perneuron)/sqrt(numel(auce_perneuron));
end

%% bar
fh=figure('Color','w','Position',[100,100,750,176]);
bar([aucc([2:6,1]);auce([2:6,1])]')
legend({'correct','error'});
set(gca,'XTick',1:6,'XTickLabel',phase([2:6,1]))
xlabel('Phase')
ylabel('AUC')
exportgraphics(fh,fullfile('F:','neupix','samp_auc_diff_bin.pdf'),'ContentType','vector');

fh=figure('Color','w','Position',[100,100,750,176]);
bar([aucc_m([2:6,1]);auce_m([2:6,1])]')
legend({'correct','error'});
set(gca,'XTick',1:6,'XTickLabel',phase([2:6,1]))
xlabel('Phase')
ylabel('AUC')
exportgraphics(fh,fullfile('F:','neupix','samp_auc_diff_bin.pdf'),'ContentType','vector');

end



%% function
function out=load_bin_sel(type)
out=[];
phase=["presample","sample","early_delay","late_delay","test","postreward"];
for i=1:numel(phase)
fstr=load(fullfile('F:\neupix\per_sec\diff_bin_0316\',sprintf('%s_34regions.mat',phase(i))));
if strcmp(type,'sample')
    out=[out,fstr.samp_any_sel];
elseif strcmp(type,'duration')
    out=[out,fstr.dur_any_sel];
elseif strcmp(type,'ordinal')
    out=[out,fstr.seq_any_sel];    
end
end
end

function out=auc_per_su(data,pos,xq)
out=struct();
labels=(1:numel([data{1};data{2}]))>numel(data{1});
scores=[data{1};data{2}];
[xc,yc,~,out.aucc]=perfcurve(labels,scores,pos);
[G,uxc]=findgroups(xc);
mmy=splitapply(@(x) mean(x), yc, G);
out.yc=interp1(uxc,mmy,xq);

labels=(1:numel([data{3};data{4}]))>numel(data{3});
scores=[data{3};data{4}];
[xe,ye,~,out.auce]=perfcurve(labels,scores,pos);
[G,uxe]=findgroups(xe);
mmy=splitapply(@(x) mean(x), ye, G);
out.ye=interp1(uxe,mmy,xq);

end
