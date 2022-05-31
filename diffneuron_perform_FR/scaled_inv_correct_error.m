
function sust_trans_correct_error(type,opt)
arguments
    type (1,:) char {mustBeMember(type,{'both','only6','only3'})} = 'both'
    opt.single_bin (1,1) logical = true
    opt.bin (1,:) char {mustBeMember( opt.bin,{'first1','late1'})} = 'first1'
    opt.plot_scatter (1,1) logical = false
    opt.plot_per_su (1,1) logical = true
    opt.plot_per_trial (1,1) logical = false
    opt.ctx_only (1,1) logical = false
    opt.plot_showcase (1,1) logical = false
end
set(groot,'defaultTextFontSize',10);
fidx=1;
colors={'r','b'};
% out=plot_wave_3_6();

meta_6=ephys.util.load_meta('delay',6);
meta_3=ephys.util.load_meta('delay',6);

% load('F:\neupix\per_sec\wave_corr.mat','out')
% [scale,invarient]=deal(zeros(size(meta_6.sess)));
% for i=1:numel(out.s1.sess)
%     if nnz(meta_6.sess==out.s1.sess(i)&meta_6.allcid==out.s1.su(i))==0
%         disp(i)
%     end
%     scale(meta_6.sess==out.s1.sess(i)&meta_6.allcid==out.s1.su(i),1)=out.s1.corr(i,1)>=0.75&out.s1.corr(i,2)<0.75;
%     invarient(meta_6.sess==out.s1.sess(i)&meta_6.allcid==out.s1.su(i),1)=out.s1.corr(i,1)<0.75&out.s1.corr(i,2)>=0.75;
%     scale(meta_6.sess==out.s2.sess(i)&meta_6.allcid==out.s2.su(i),2)=out.s2.corr(i,1)>=0.75&out.s2.corr(i,2)<0.75;
%     invarient(meta_6.sess==out.s2.sess(i)&meta_6.allcid==out.s2.su(i),2)=out.s2.corr(i,1)<0.75&out.s2.corr(i,2)>=0.75;
% end
% save('F:\neupix\per_sec\wave_corr.mat','out')
% corr=out;
load('F:\neupix\per_sec\scale_invarient.mat','scale','invarient','corr')

stats=struct();
types=[1,2;3,4]; % don't consider sustained and transient
% out=plot_wave_3_6();

for onetype=types'
    if strcmp(type,'both')        
        if opt.ctx_only
            typesel=find(any(invarient,2)' & ismember(meta_3.mem_type,onetype) & strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),''));
        else
            typesel=find(any(invarient,2)' & ismember(meta_3.mem_type,onetype) & ismember(meta_6.mem_type,onetype));
        end
        if opt.single_bin
            if strcmp(opt.bin,'first1')
                subsel=meta_3.per_bin(1,typesel)~=0 & meta_6.per_bin(1,typesel)~=0;
                typesel=typesel(subsel);
            elseif strcmp(opt.bin,'late1')
                subsel=meta_3.per_bin(end,typesel)~=0 & meta_6.per_bin(end,typesel)~=0;
                typesel=typesel(subsel);
            end
        end        
    elseif strcmp(type,'3s')        
        if opt.ctx_only
            typesel=find(ismember(meta_3.mem_type,onetype) & ~ismember(meta_6.mem_type,onetype) & strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),''));
        else
            typesel=find(ismember(meta_3.mem_type,onetype) & ~ismember(meta_6.mem_type,onetype));
        end
        if opt.single_bin
            if strcmp (opt.bin,'first1')
                subsel=meta_3.per_bin(1,typesel)~=0;
                typesel=typesel(subsel);
            elseif strcmp (opt.bin,'late1')
                subsel=meta_3.per_bin(end,typesel)~=0;
                typesel=typesel(subsel);
            end
        end
        
    elseif strcmp(type,'6s')        
        if opt.ctx_only
            typesel=find(~ismember(meta_3.mem_type,onetype) & ismember(meta_6.mem_type,onetype) & strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),''));
        else
            typesel=find(~ismember(meta_3.mem_type,onetype) & ismember(meta_6.mem_type,onetype));
        end
        if opt.single_bin
            if strcmp (opt.bin,'first1')
                subsel=meta_6.per_bin(1,typesel)~=0;
                typesel=typesel(subsel);
            elseif strcmp (opt.bin,'late1')
                subsel=meta_6.per_bin(end,typesel)~=0;
                typesel=typesel(subsel);
            end
        end
    end   
    homedir=ephys.util.getHomedir('type','raw');
    for delay=[3,6]
        idx=1;
        stats.(sprintf('delay%ds_type%d',delay,(onetype(1)-1)/2+1))=nan(0,4);
        stats.(sprintf('delay%ds_type%d_pertrial',delay,(onetype(1)-1)/2+1))=cell(0,4);
        if opt.single_bin
            if strcmp (opt.bin,'first1')
                bins=5;
            elseif strcmp (opt.bin,'late1')
                bins=delay+4;
            elseif strcmp (opt.bin,'sample')
                bins=4;
            elseif strcmp (opt.bin,'test')
                bins=delay+5;
            elseif strcmp (opt.bin,'ITI')
                bins=delay+8;         
            end
        else
            bins=find(meta.per_bin(:,ii)~=0)+4;
        end
        for ii=typesel
            fpath=fullfile(homedir,meta_3.allpath{ii},'FR_All_1000.hdf5');
            trials=h5read(fpath,'/Trials');
            suid=h5read(fpath,'/SU_id');
            fr=h5read(fpath,'/FR_All');

            if nnz(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==3)==0 ...
                    ||nnz(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==3)==0 ...
                    ||nnz(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6)==0 ...
                    ||nnz(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6)==0
                continue
            end
            cs1=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==delay,suid==meta_3.allcid(ii),bins),3);
            cs2=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==delay,suid==meta_3.allcid(ii),bins),3);
            es1=mean(fr(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==delay,suid==meta_3.allcid(ii),bins),3);
            es2=mean(fr(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==delay,suid==meta_3.allcid(ii),bins),3);
            
            delaymm=mean([mean(cs1),mean(cs2)]);
            delaystd=std([cs1;cs2]);
            if delaystd==0, continue;  end
            
            cs1=(cs1-delaymm)./delaystd;
            cs2=(cs2-delaymm)./delaystd;
            es1=(es1-delaymm)./delaystd;
            es2=(es2-delaymm)./delaystd;
            
            %         if onetype==2 && mean(cs1)<mean(cs2),keyboard;end
            if opt.plot_showcase
                if min(cellfun(@(x) numel(x),{cs1,cs2,es1,es2}))<10
                    continue
                end
                [~,~,~,AUCC]=perfcurve([zeros(size(cs1));ones(size(cs2))],[cs1;cs2],0);
                [~,~,~,AUCE]=perfcurve([zeros(size(es1));ones(size(es2))],[es1;es2],0);
                
                if abs(AUCC-0.5)>0.30 && abs(AUCE-0.5)<0.15
                    figure()
                    subplot(1,3,1)
                    hold on
                    cs1c=squeeze(mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6,suid==meta.allcid(ii),3:10),[1,2]));
                    cs2c=squeeze(mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6,suid==meta.allcid(ii),3:10),[1 2]));
                    es1c=squeeze(mean(fr(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6,suid==meta.allcid(ii),3:10),[1 2]));
                    es2c=squeeze(mean(fr(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6,suid==meta.allcid(ii),3:10),[1 2]));
                    
                    plot(cs1c,'-r')
                    plot(cs2c,'-b')
                    plot(es1c,':r')
                    plot(es2c,':b')
                    xline(1.5,'--k')
                    
                    
                    subplot(1,3,2)
                    hold on
                    histogram(cs1,-3:0.3:3,'FaceColor','r','FaceAlpha',0.4)
                    histogram(cs2,-3:0.3:3,'FaceColor','b','FaceAlpha',0.4)
                    subplot(1,3,3)
                    hold on
                    histogram(es1,-3:0.3:3,'FaceColor','r','FaceAlpha',0.4)
                    histogram(es2,-3:0.3:3,'FaceColor','b','FaceAlpha',0.4)
                    sgtitle(ii);
                    keyboard()
                end
            end
            
            if opt.plot_scatter
                ch(fidx)=scatter(mean(cs1),mean(cs2),36,colors{fidx},'MarkerFaceColor',colors{fidx},'MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
                eh=scatter(mean(es1),mean(es2),16,'k','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
                plot([mean(cs1),mean(es1)],[mean(cs2),mean(es2)],':','LineWidth',0.5,'Color',colors{fidx});
            end
            
            stats.(sprintf('delay%ds_type%d',delay,(onetype(1)-1)/2+1))(idx,:)=[mean(cs1),mean(cs2),mean(es1),mean(es2)];
            stats.(sprintf('delay%ds_type%d_pertrial',delay,(onetype(1)-1)/2+1))(idx,:)={cs1,cs2,es1,es2};
            
            idx=idx+1;
            
        end
    end
    fidx=fidx+1;
end
if opt.plot_scatter
    plot([-10,10],[-10,10],'--k')
    xlim([-2,2]);
    ylim(xlim());
    xlabel('Sample 1 normalized FR (Z-score)')
    ylabel('Sample 2 normalized FR (Z-score)')
    legend([ch(1),ch(2),eh],{'S1 correct trials','S2 correct trials','Error trials'});
end

% SU per session mean
% correct trial
if opt.plot_per_su
    fh=figure('Color','w','Position',[100,100,750,176]);
    for delay=6
        subplot(1,3,1);
        hold on;
        ph=histogram([stats.(sprintf('delay%ds_type1',delay))(:,1);stats.(sprintf('delay%ds_type2',delay))(:,2)],-2:0.1:2,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
        nph=histogram([stats.(sprintf('delay%ds_type1',delay))(:,2);stats.(sprintf('delay%ds_type2',delay))(:,1)],-2:0.1:2,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
        xline(nanmean([stats.(sprintf('delay%ds_type1',delay))(:,1);stats.(sprintf('delay%ds_type2',delay))(:,2)]),'--r','LineWidth',1);
        xline(nanmean([stats.(sprintf('delay%ds_type1',delay))(:,2);stats.(sprintf('delay%ds_type2',delay))(:,1)]),'--k','LineWidth',1)
        xlim([-1.5,1.5]);
        title('Correct trials')
        xlabel('Normalized FR (Z-Score)');
        ylabel('Probability')
        text(max(xlim()),max(ylim()),num2str(numel(ph.Data)),'HorizontalAlignment','right','VerticalAlignment','top');
        
        %error trial
        subplot(1,3,2);
        hold on;
        ph=histogram([stats.(sprintf('delay%ds_type1',delay))(:,3);stats.(sprintf('delay%ds_type2',delay))(:,4)],-2:0.1:2,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
        nph=histogram([stats.(sprintf('delay%ds_type1',delay))(:,4);stats.(sprintf('delay%ds_type2',delay))(:,3)],-2:0.1:2,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
        xline(nanmean([stats.(sprintf('delay%ds_type1',delay))(:,3);stats.(sprintf('delay%ds_type2',delay))(:,4)]),'--r','LineWidth',1);
        xline(nanmean([stats.(sprintf('delay%ds_type1',delay))(:,4);stats.(sprintf('delay%ds_type2',delay))(:,3)]),'--k','LineWidth',1)        
        xlim([-1.5,1.5]);
        title('Error trials');
        legend([ph,nph],{'Prefered','Non-prefered'});
        xlabel('Normalized FR (Z-Score)');
        ylabel('Probability')
        
        %auc
        onecolumn=size(stats.(sprintf('delay%ds_type1',delay)),1)+size(stats.(sprintf('delay%ds_type2',delay)),1);
        [xc,yc,~,aucc]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('delay%ds_type1',delay))(:,1)...
            ,;stats.(sprintf('delay%ds_type2',delay))(:,2);stats.(sprintf('delay%ds_type1',delay))(:,2)...
            ,;stats.(sprintf('delay%ds_type2',delay))(:,1)],0);
        [xe,ye,~,auce]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('delay%ds_type1',delay))(:,3)...
            ,;stats.(sprintf('delay%ds_type2',delay))(:,4);stats.(sprintf('delay%ds_type1',delay))(:,4)...
            ,;stats.(sprintf('delay%ds_type2',delay))(:,3)],0);
        subplot(1,3,3);
        hold on;
        hc=plot(xc,yc,'-r','LineWidth',1);
        he=plot(xe,ye,'-k','LineWidth',1);
        legend([hc,he],{sprintf('Correct trials AUC=%0.3f',aucc),...
            sprintf('Error trials AUC=%0.3f',auce)},...
            'Location','southeast');
        xlabel('False positive rate (fpr)');
        ylabel('True positive rate (tpr)');
        if opt.single_bin
            sgtitle(sprintf('%s averaged cross-trial bin %s',type,opt.bin));
        else
            sgtitle(sprintf('%s averaged cross-trial',type));
        end
    end
    exportgraphics(fh,fullfile('F:','neupix',sprintf('%s_cross_trial_bin%s.pdf',type,opt.bin)),'ContentType','vector');
end

%Per trial
if opt.plot_per_trial
    %correct trials
    figure('Color','w','Position',[100,100,1200,300]);
    subplot(1,3,1);
    hold on;
    xbins=-3:0.15:3;
    prefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,1);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,2)],...
        'UniformOutput',false));
    nonprefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,2);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,1)],...
        'UniformOutput',false));
    
    ph=bar(xbins(1:end-1)+0.075,mean(prefmat),'FaceColor','r','FaceAlpha',0.4);
    nph=bar(xbins(1:end-1)+0.075,mean(nonprefmat),'FaceColor','k','FaceAlpha',0.4);
    
    pcom=sum((xbins(1:end-1)+0.075).*mean(prefmat))./sum(mean(prefmat));
    npcom=sum((xbins(1:end-1)+0.075).*mean(nonprefmat))./sum(mean(nonprefmat));
    xline(pcom,'--r','LineWidth',1);
    xline(npcom,'--k','LineWidth',1);
    title('Correct trials')
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
    
    %error trials
    subplot(1,3,2);
    hold on;
    xbins=-3:0.15:3;
    prefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,3);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,4)],...
        'UniformOutput',false));
    nonprefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,4);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,3)],...
        'UniformOutput',false));
    
    ph=bar(xbins(1:end-1)+0.075,nanmean(prefmat),'FaceColor','r','FaceAlpha',0.4);
    nph=bar(xbins(1:end-1)+0.075,nanmean(nonprefmat),'FaceColor','k','FaceAlpha',0.4);
    
    pcom=sum((xbins(1:end-1)+0.075).*nanmean(prefmat))./sum(nanmean(prefmat));
    npcom=sum((xbins(1:end-1)+0.075).*nanmean(nonprefmat))./sum(nanmean(nonprefmat));
    xline(pcom,'--r','LineWidth',1);
    xline(npcom,'--k','LineWidth',1);
    
    title('Error trials');
    legend([ph,nph],{'Prefered','Non-prefered'});
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
    
    %auc
    
    subplot(1,3,3);
    hold on;
    xq=0.005:0.01:1;
    prefs1=arrayfun(@(x) auc_per_su(stats.(sprintf('type%d_pertrial',types(1)))(x,:),0,xq),...
        1:size(stats.(sprintf('type%d_pertrial',types(1))),1));
    prefs2=arrayfun(@(x) auc_per_su(stats.(sprintf('type%d_pertrial',types(2)))(x,:),1,xq),...
        1:size(stats.(sprintf('type%d_pertrial',types(2))),1));
    
    hc=plot(xq,mean([cell2mat({prefs1.yc}.');cell2mat({prefs2.yc}.')]),'-r','LineWidth',1);
    he=plot(xq,mean([cell2mat({prefs1.ye}.');cell2mat({prefs2.ye}.')]),'-k','LineWidth',1);
    legend([hc,he],{sprintf(...
        'Correct trials AUC=%0.2f',mean([prefs1.aucc,prefs2.aucc]))...
        ,sprintf(...
        'Error trials AUC=%0.2f',mean([prefs1.auce,prefs2.auce]))},...
        'Location','southeast');
    xlabel('False positive rate (fpr)');
    ylabel('True positive rate (tpr)');
    sgtitle(sprintf('%s stats from trials, averaged over neurons',type));
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

%% calculated COM
function out=plot_wave_3_6(opt)
arguments
    opt.sess (1,1) double {mustBeMember(opt.sess,1:116)} = 102
    opt.sortby (1,:) char {mustBeMember(opt.sortby,{'3s','6s'})} = '6s'
    opt.exportgraphics (1,1) logical = false
    opt.plot_global (1,1) logical = false
    opt.plot_session (1,1) logical = false
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1
    opt.bootrpt (1,1) double {mustBeInteger,mustBePositive} = 3;
    opt.plot_2d_corr (1,1) logical = false
    opt.plot_corr_dist (1,1) logical = true
end
out=struct();
persistent com_map_3_sel com_map_6_sel     

if isempty(com_map_3_sel)
    com_map_3_sel=get_com_map('curve',true,'rnd_half',false,'delay',3,'cell_type','all');
    com_map_6_sel=get_com_map('curve',true,'rnd_half',false,'delay',6,'cell_type','all');
end
for alt_comb=opt.comb_set
    com_map_3=com_map_3_sel;
    com_map_6=com_map_6_sel;                            
            
    fss=intersect(fieldnames(com_map_6),fieldnames(com_map_3)).';
    
    for pref=["s1","s2"]
        [immata,immatb,su,sess]=deal([]); 
        for fs1=fss
            fs=char(fs1);            
            curr_key=intersect(cell2mat(com_map_3.(fs).(pref).keys),cell2mat(com_map_6.(fs).(pref).keys));
            heat3=cell2mat(com_map_3.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
            heat6=cell2mat(com_map_6.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
            immata=[immata;heat3];
            immatb=[immatb;heat6];  
            su=[su,curr_key];
            sess=[sess,repmat(str2double(regexp(fs,'(?<=s)\d*','once','match')),size(curr_key))];
        end
        flat_mat=@(x) x(:);
        scale_mat=@(x) (x(:,1:2:size(x,2))+x(:,2:2:size(x,2)))./2;
        for i=1:size(immata,1)
            out.(pref).corr(i,1)=corr(flat_mat(immata(i,:)),flat_mat(scale_mat(immatb(i,:)))); % scaled
            out.(pref).corr(i,2)=corr(flat_mat(immata(i,:)),flat_mat(immatb(i,1:12))); % invariant
            out.(pref).corr(i,3)=corr(flat_mat(immata(i,:)),flat_mat(immatb(i,13:end))); % delay
        end
        out.(pref).su=su;
        out.(pref).sess=sess;
        
    end  
    
end

save('F:\neupix\per_sec\wave_corr.mat','out')
end

function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.per_sec_stats (1,1) logical = false % calculate COM using per-second mean as basis for normalized firing rate, default is coss-delay mean
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.cell_type (1,:) char {mustBeMember(opt.cell_type,{'any_s1','any_s2','any_nonmem','ctx_sel','ctx_trans','all'})} = 'ctx_trans' % select (sub) populations
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
    opt.partial (1,:) char {mustBeMember(opt.partial,{'full','early3in6','late3in6'})}='full' % for TCOM correlation between 3s and 6s trials
    opt.plot_COM_scheme (1,1) logical = false % for TCOM illustration
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
    opt.alt_3_6 (1,1) logical = false
end
persistent com_str opt_

if opt.delay==6
    warning('Delay set to default 6')
end

if true || isempty(com_str) || ~isequaln(opt,opt_)
    meta_str=ephys.util.load_meta('type','neupix','delay',opt.delay);
    if opt.alt_3_6
        meta_str_alt=ephys.util.load_meta('type','neupix','delay',setdiff([3,6],opt.delay));
    end
    homedir=ephys.util.getHomedir('type','raw');
    fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
    com_str=struct();
    for ii=1:size(fl,1)
        if strlength(opt.onepath)==0
            dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
            fpath=fullfile(fl(ii).folder,fl(ii).name);
        else
            dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
            if isempty(dpath)
                dpath=opt.onepath;
            end
            fpath=fullfile(homedir,dpath,'FR_All_ 250.hdf5');
        end
        fpath=replace(fpath,'\',filesep());
        pc_stem=replace(dpath,'/','\');
        sesssel=startsWith(meta_str.allpath,pc_stem);
        if ~any(sesssel), continue;end
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        switch opt.cell_type
            case 'any_s1'
                mcid1=meta_str.allcid(ismember(meta_str.mem_type,1:2) & sesssel.');
                mcid2=[];
            case 'any_s2'
                mcid1=[];
                mcid2=meta_str.allcid(ismember(meta_str.mem_type,3:4) & sesssel.');
            case 'any_nonmem'
                mcid1=meta_str.allcid(meta_str.mem_type==0 & sesssel.');
                mcid2=[];
            case 'ctx_sel'
                mcid1=meta_str.allcid(ismember(meta_str.mem_type,1:2) & sesssel.'& strcmp(meta_str.reg_tree(2,:),'CTX'));
                mcid2=meta_str.allcid(ismember(meta_str.mem_type,3:4) & sesssel.'& strcmp(meta_str.reg_tree(2,:),'CTX'));
            case 'ctx_trans'
                mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
                mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
                if opt.alt_3_6    
                    mcid1alt=meta_str_alt.allcid(meta_str_alt.mem_type==2 & sesssel.' & strcmp(meta_str_alt.reg_tree(2,:),'CTX'));
                    mcid2alt=meta_str_alt.allcid(meta_str_alt.mem_type==4 & sesssel.' & strcmp(meta_str_alt.reg_tree(2,:),'CTX'));
                    mcid1=setdiff(mcid1alt,mcid1);
                    mcid2=setdiff(mcid2alt,mcid2);
                end
             case 'all'
                mcid1=meta_str.allcid(sesssel.');
                mcid2=meta_str.allcid(sesssel.');
        end
        msel1=find(ismember(suid,mcid1));
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2)
            if strlength(opt.onepath)==0
                continue
            else
                break
            end
        end
        sessid=ephys.path2sessid(pc_stem);
        s1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);
        s2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);
        
        e1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,10)==0);
        e2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,10)==0);
        
        sess=['s',num2str(sessid)];
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        if opt.rnd_half
            for ff=["s1a","s2a","s1b","s2b","s1e","s2e"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["s1aheat","s2aheat","s1acurve","s2acurve","s1aanticurve","s2aanticurve",...
                        "s1bheat","s2bheat","s1bcurve","s2bcurve","s1banticurve","s2banticurve",...
                        "s1eheat","s2eheat","s1ecurve","s2ecurve","s1eanticurve","s2eanticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            s1a=randsample(s1sel,floor(numel(s1sel)./2));
%             s1a=s1sel(1:2:end);
            s1b=s1sel(~ismember(s1sel,s1a));

            s2a=randsample(s2sel,floor(numel(s2sel)./2));
%             s2a=s2sel(1:2:end);
            s2b=s2sel(~ismember(s2sel,s2a));
            if nnz(s1a)>2 && nnz(s1b)>2 && nnz(s2a)>2 && nnz(s2b)>2 && nnz(e1sel)>2 && nnz(e2sel)>2
                com_str=per_su_process(sess,suid,msel1,fr,s1a,s2a,com_str,'s1a',opt);
                com_str=per_su_process(sess,suid,msel1,fr,s1b,s2b,com_str,'s1b',opt);
                com_str=per_su_process(sess,suid,msel2,fr,s2a,s1a,com_str,'s2a',opt);
                com_str=per_su_process(sess,suid,msel2,fr,s2b,s1b,com_str,'s2b',opt);
                com_str=per_su_process(sess,suid,msel1,fr,e1sel,e2sel,com_str,'s1e',opt);
                com_str=per_su_process(sess,suid,msel2,fr,e2sel,e1sel,com_str,'s2e',opt);
            end
        else
            for ff=["s1","s2"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["s1heat","s2heat","s1curve","s2curve","s1anticurve","s2anticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            com_str=per_su_process(sess,suid,msel1,fr,s1sel,s2sel,com_str,'s1',opt);
            com_str=per_su_process(sess,suid,msel2,fr,s2sel,s1sel,com_str,'s2',opt);
        end
        if ~strlength(opt.onepath)==0
            break;
        end
    end
end
com_str_=com_str;
opt_=opt;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
%TODO if decision
if opt.decision
    stats_window=(opt.delay*4+17):44;
else
    if opt.delay==6 && strcmp(opt.partial,'early3in6')
        stats_window=17:28;
    elseif opt.delay==6 && strcmp(opt.partial,'late3in6')
        stats_window=29:40;
    else
        stats_window=17:(opt.delay*4+16);
    end
end
for su=reshape(msel,1,[])
    perfmat=squeeze(fr(pref_sel,su,:));
    npmat=squeeze(fr(nonpref_sel,su,:));
    basemm=mean([mean(perfmat(:,stats_window),1);mean(npmat(:,stats_window),1)]);

%% Obsolete code to calculate COM of selectivity indices        
%         sel_idx(sel_idx<0)=0;
%         com=sum((1:numel(stats_window)).*sel_idx(stats_window))./sum(sel_idx(stats_window));
%         if ~isfinite(com)
%             fprintf('Moved sess %d su%d, infinite TCOM\n',sess,su)
% %             keyboard()
%             com=0;
%         end

    
        if ~opt.per_sec_stats
            basemm=mean(basemm);
        end
        itimm=mean(fr([pref_sel;nonpref_sel],su,1:12),'all');
        %TODO compare the effect of smooth
        if strcmp(opt.cell_type,'any_nonmem')
            mm=smooth(squeeze(mean(fr([pref_sel;nonpref_sel],su,:))),5).';
            mm_pref=mm(stats_window)-itimm;
        else
            mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
            mm_pref=mm(stats_window)-basemm;
        end
        if max(mm_pref)<=0,disp('NOPEAK');continue;end % work around 6s paritial
        if opt.selidx
            sel_vec=[mean(perfmat,1);mean(npmat,1)];
            sel_idx=(-diff(sel_vec)./sum(sel_vec));
            sel_idx(all(sel_vec==0))=0;
            curve=sel_idx(stats_window);
        elseif contains(opt.cell_type,'any')
            curve=squeeze(mean(fr(pref_sel,su,:))).'-itimm;
            anticurve=squeeze(mean(fr(nonpref_sel,su,:))).'-itimm;
        elseif opt.delay==6 % full, early or late
            curve=squeeze(mean(fr(pref_sel,su,17:40))).'-basemm;
            anticurve=squeeze(mean(fr(nonpref_sel,su,17:40))).'-basemm;
        else
            curve=mm_pref;
            anticurve=squeeze(mean(fr(nonpref_sel,su,stats_window))).'-basemm;
        end
        mm_pref(mm_pref<0)=0;
        com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);
        if opt.delay==6 && strcmp(opt.partial,'late3in6')
            com=com+12;
        end
        if opt.plot_COM_scheme
%             template=[zeros(1,2),0:0.2:1,1:-0.2:0,zeros(1,10)];
%             if corr(mm_pref.',template.','type','Spearman')>0.2
                fc_scheme(curve,mm_pref,com)
%             end
        end
    com_str.(sess).(samp)(suid(su))=com;
    if opt.curve
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
        if exist('anticurve','var')
            com_str.(sess).([samp,'anticurve'])(suid(su))=anticurve;
        end
        heatcent=squeeze(fr(pref_sel,su,stats_window))-basemm; %centralized norm. firing rate for heatmap plot
        heatnorm=heatcent./max(abs(heatcent));
        heatnorm(heatnorm<0)=0;
        if size(heatnorm,1)>10
            if numel(curve)>numel(stats_window)
                curve=curve(stats_window);
            end
            cc=arrayfun(@(x) min(corrcoef(heatnorm(x,:),curve),[],'all'),1:size(heatnorm,1));
            [~,idx]=sort(cc,'descend');
            heatnorm=heatnorm(idx(1:10),:);
        end
        com_str.(sess).([samp,'heat'])(suid(su))=heatnorm;
    end
end
end

%% for COM scheme illustration
function fc_scheme(curve,mm_pref,com)
curve(curve<0)=0;
close all
fh=figure('Color','w', 'Position',[32,32,275,225]);
bar(mm_pref,'k');
ylabel('Baseline-subtracted firing rate w (Hz)')
set(gca,'XTick',0:4:24,'XTickLabel',0:6)
xlabel('Time t (0.25 to 6 sec in step of 0.25 sec)')
xline(com,'--r');
keyboard()
%             exportgraphics(gcf(),'COM_illustration_L0.pdf','ContentType','vector')
        
end

