for i=18 %7 non-memory
    hist_coeff_mem_nonmem(i,'non-mem',...,
        'prefix','0315',...,
        'tsbin_size',600,...,
        'postspike',true,...,
        'fc_effi',false',...,
        'fc_prob',false',...,
        'type','neupix',...,
        'laser','on',...,
        'assembly',true)
end
function hist_coeff_mem_nonmem(sess,mtype,opt)
arguments
    sess (1,1) double {mustBeInteger,mustBePositive,mustBeNonempty} =2
    mtype (1,:) char {mustBeMember(mtype,{'congru','incongru','non-mem'})} ='congru'
    opt.sess (1,1) double
    opt.prefix (1,:) char = '0315'
    opt.tsbin_size (1,1) double = 600
    opt.postspike (1,1) logical = true
    opt.fc_effi (1,1) logical = false
    opt.fc_prob (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'on'
    opt.assembly (1,1) logical = true
end

if ispc
    load('sig.mat','sig'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
else
    load(fullfile('/home','hem','datashare','sig.mat'),'sig');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
end
switch mtype
    case 'congru'
        typesel=sig.sess==sess & (all(ismember(sig.mem_type,1:2),2) | all(ismember(sig.mem_type,3:4),2));
    case 'incongru'
        typesel=sig.sess==sess ...
            & ((ismember(sig.mem_type(:,1),1:2) & ismember(sig.mem_type(:,2),3:4)) ...
            |(ismember(sig.mem_type(:,1),3:4) & ismember(sig.mem_type(:,2),1:2)));
    case 'non-mem'
        typesel=sig.sess==sess & all(sig.mem_type==0,2);
end
dim=nnz(typesel);
if dim==0
    return
end

idces=find(typesel);
postspk=cell(0); 
sidx=1;
sess_suids=sig.suid(idces,:);
for i=reshape(idces,1,[])    
    fprintf('%d of %d\n',find(idces==i),dim);
    if nnz(sess_suids(:,2)==sig.suid(i,2))<=2 
        continue
    end
    if sidx>1 && nnz(suid(:,end)==sig.suid(i,2))>0         
         continue
    end
    stats.sess=sig.sess(i);
    stats.sess_suids=sig.suid(i,2);
    reg=tag_hist_reg(stats,'type','neupix');
    if ~(cellfun(@(x)ismember(x{1},{'CH','BS'}),reg{1}) && cellfun(@(x)~isempty(x{5}),reg{1}))
        continue
    end
        
    pair=sess_suids(sess_suids(:,2)==sig.suid(i,2),:);
    p=nchoosek(1:size(pair,1),3);
    stats.sess=repmat(sig.sess(i),size(p,1),1);
    stats.sess_suids=[arrayfun(@(x)pair(x,1),p),repmat(sig.suid(i,2),size(p,1),1)];
    reg=tag_hist_reg(stats,'type','neupix');
    [is_diff,is_same]=diff_at_level(reg,5);   
    p(~(is_diff|is_same),:)=[];   
    if isempty(p)
        continue
    end
    
    for pp=1:size(p,1) 
        [postspk{sidx},avail(sidx)]=...
            history_coeff_coactivate_trial(sig.sess(i),pair(p(pp,:),:),...
            'tsbin_size',opt.tsbin_size,...
            'postspike',opt.postspike,...
            'fc_effi',opt.fc_effi,...
            'fc_prob',opt.fc_prob,...
            'type',opt.type,...
            'laser',opt.laser,...
            'assembly',opt.assembly,...
            'sess',sess);
        
        suid(sidx,:)=[pair(p(pp,:),1);unique(pair(p(pp,:),2));]';
        suid_sess(sidx,:)=sig.sess(i);
        sidx=sidx+1;
    end 
    clear pair p
end

%TODO return if whatever-empty

% bad=all(postspk(:,1:11)==0,2);
% postspk(bad,:)=[];
% sess_suids(bad,:)=[];
% maxiter(bad,:)=[];
if ~isempty(postspk)
    postspk(~avail)=[];
    suid(~avail,:)=[];
    suid_sess(~avail,:)=[];
    if ispc
        save(fullfile('F:','neupix','STP','coactivate','0702-3neuron',sprintf('%s_stp_%s_%d_%d.mat',...
            opt.prefix,mtype,sess,opt.tsbin_size)),...
            'postspk','suid','suid_sess','mtype');
    else
        save(fullfile('/home','hem','datashare','0620-time',sprintf('%s_stp_%s_%d_%d.mat',...
            opt.prefix,mtype,sess,opt.tsbin_size)),...
            'postspk','suid','suid_sess','mtype');
    end
end
end

function [out,out_avail]=history_coeff_coactivate_trial(sessid,suid,opt)
arguments
    sessid (1,1) int32 
    suid (:,2) int32        
    opt.sess (1,1) double
    opt.prefix (1,:) char = '0315'
    opt.tsbin_size (1,1) double = 600
    opt.postspike (1,1) logical = false
    opt.fc_effi (1,1) logical = false
    opt.fc_prob (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
    opt.assembly (1,1) logical = false
    % Optional TODO EPOCH
end

persistent spike
if isempty(spike) || (~isempty(spike)&& spike.sess~=sessid)
   [spike.spkID,spike.spkTS,spike.trials,~,~]=getSPKID_TS(sessid); 
   spike.sess=sessid;
end

bitmask=2.^(0:size(suid,1)-1)';
out_avail=false;

[avail(1),tspre1]=pre_process(suid(1,1),spike.spkTS,spike.spkID...
    ,spike.trials(spike.trials(:,8)==6&spike.trials(:,9)==1&spike.trials(:,10)==1,:),opt);
[avail(2),tspre2]=pre_process(suid(2,1),spike.spkTS,spike.spkID...
    ,spike.trials(spike.trials(:,8)==6&spike.trials(:,9)==1&spike.trials(:,10)==1,:),opt);
[avail(3),tspost]=pre_process(suid(1,2),spike.spkTS,spike.spkID...
    ,spike.trials(spike.trials(:,8)==6&spike.trials(:,9)==1&spike.trials(:,10)==1,:),opt);
[avail(4),tspre3]=pre_process(suid(3,1),spike.spkTS,spike.spkID...
    ,spike.trials(spike.trials(:,8)==6&spike.trials(:,9)==1&spike.trials(:,10)==1,:),opt);
if ~any(avail)
    return
end
shuffle=["raw","shuffle100ms"];
for Isshuffle=shuffle
    if strcmp(Isshuffle,'raw')        
        tmin=-3;
        tmax=11;
        bin_size=opt.tsbin_size/30000;
        histpre1=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspre1,'UniformOutput',false));
        histpre2=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspre2,'UniformOutput',false));
        histpost=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspost,'UniformOutput',false));
        histpre3=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspre3,'UniformOutput',false));
        if nnz(histpost)==0
           out.(Isshuffle)='error';
           return            
        end
        temp=zeros(10,0);
        for i=reshape((0:length(tspre1)-1)*(tmax-tmin)/bin_size+(11:(tmax-tmin)/bin_size)',1,[])
            if histpost(i)==0
                continue
            end
            temp(:,end+1)=[histpre1(i-10:i-1)',histpre2(i-10:i-1)',histpre3(i-10:i-1)']*bitmask;
        end        
        for i=1:10
            x(1,i)=nnz(temp(i,:)==0);
            x(2,i)=nnz(temp(i,:)==1);
            x(3,i)=nnz(temp(i,:)==2);
            x(4,i)=nnz(temp(i,:)==7);
        end
        out.(Isshuffle)=x./size(temp,2);
        out_avail=true;
    else
        tmin=-3;
        tmax=11;
        bin_size=opt.tsbin_size/30000;
        out_temp=zeros(4,10,0);
        for shuffle_time=1:100            
            switch Isshuffle
                case 'shuffle20ms'  
                    histpre1=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*1000,1,length(x)),tspre1,'UniformOutput',false),'UniformOutput',false));
                    histpre2=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*1000,1,length(x)),tspre2,'UniformOutput',false),'UniformOutput',false));
                    histpost=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*1000,1,length(x)),tspost,'UniformOutput',false),'UniformOutput',false));
                    histpre3=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*1000,1,length(x)),tspre3,'UniformOutput',false),'UniformOutput',false));
                case 'shuffle100ms'
                    histpre1=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*5000,1,length(x)),tspre1,'UniformOutput',false),'UniformOutput',false));
                    histpre2=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*5000,1,length(x)),tspre2,'UniformOutput',false),'UniformOutput',false));
                    histpost=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*5000,1,length(x)),tspost,'UniformOutput',false),'UniformOutput',false));
                    histpre3=cell2mat(cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,...
                        cellfun(@(x)x+0.001*randi([-1 1]*bin_size*5000,1,length(x)),tspre3,'UniformOutput',false),'UniformOutput',false));
            end
            if nnz(histpost)==0
                continue
            end
            temp=zeros(10,0);            
            for i=reshape((0:length(tspre1)-1)*(tmax-tmin)/bin_size+(11:(tmax-tmin)/bin_size)',1,[])
                if histpost(i)==0
                    continue
                end
                temp(:,end+1)=[histpre1(i-10:i-1)',histpre2(i-10:i-1)',histpre3(i-10:i-1)']*bitmask;
            end
            for i=1:10
                x(1,i)=nnz(temp(i,:)==0);
                x(2,i)=nnz(temp(i,:)==1);
                x(3,i)=nnz(temp(i,:)==2);
                x(4,i)=nnz(temp(i,:)==7);
            end            
             out_temp(:,:,end+1)=x./size(temp,2);
        end
        if exist('out_temp','var')
            out.(Isshuffle)=mean(out_temp,3);            
        else
            out.(Isshuffle)='error';
            out_avail=false;
            return
        end      
    end    
end
end

function [avail,out]=pre_process(id,spkTS,spkId,trials,~)
timlim=[-4,14];
sps=30000;
if ispc
    addpath(fullfile('D:','code','fieldtrip-20200320'))
else
    addpath('fieldtrip-20200320')
end
ft_defaults
FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(id)));
FT_SPIKE.timestamp{1}=cell(1,numel(id));
FT_SPIKE.timestamp{1}=spkTS(spkId==id)';
sps=30000;
cfg=struct();
cfg.trl=[trials(:,1)+timlim(1)*sps,trials(:,1)+timlim(2)*sps,zeros(size(trials,1),1)+timlim(1)*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
for i=1:size(FT_SPIKE.trialtime,1)
    out{i}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1}==i);
end
avail=true;


end

function [spkID_,spkTS_,trials_,SU_id_,folder_]=getSPKID_TS(fidx)
arguments
    fidx (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(fidx,1)}
end

persistent spkID spkTS trials SU_id folder fidx_

if isempty(fidx_) || fidx ~= fidx_
    if ispc
        homedir=fullfile('F:','neupix','SPKINFO'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        homedir=fullfile('/home','zx','neupix','SPKINFO'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    folder=replace(sessid2path(fidx),'\',filesep());
    trials=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/Trials');
    SU_id=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/SU_id');
%     FR_All=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/FR_All');
    spkID=[];spkTS=[];
    if sum(trials(:,9))<40 || numel(SU_id)<2  % apply well-trained criteria
        trials=[];SU_id=[];return;
    end

    cstr=h5info(fullfile(homedir,folder,'spike_info.hdf5')); % probe for available probes
    for prb=1:size(cstr.Groups,1) % concatenate same session data for cross probe function coupling
        prbName=cstr.Groups(prb).Name;
        spkID=cat(1,spkID,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/clusters']));
        spkTS=cat(1,spkTS,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/times']));
    end

    susel=ismember(spkID,SU_id); % data cleaning by FR and contam rate criteria %TODO optional waveform cleaning
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
end

spkID_=spkID;
spkTS_=spkTS;
trials_=trials;
SU_id_=SU_id;
folder_=folder;

end

function [out,homedir]=sessid2path(sessid,opt)
arguments
    sessid (1,1) double {mustBeInteger,mustBePositive}
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent map
if ispc
    homedir=fullfile('E:','hem','bzdata');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
else
    homedir=fullfile('/home','hem','datashare');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
end
if isempty(map)
    
    allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','int32','ValueType','char');
    for i=1:numel(uidx)
        map(uidx(i))=upath{i};
    end
end

out=map(sessid);

end
function reg=tag_hist_reg(per_type_stats,opt)
arguments
    per_type_stats (1,1) struct
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent meta
if isempty(meta)
    if ispc
        homedir=fullfile('E:','hem','bzdata');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
    else
        homedir=fullfile('/home','hem','datashare');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
    end
    meta.trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
    meta.wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p');
    meta.selec=h5read(fullfile(homedir,'transient_6.hdf5'),'/selectivity');
    meta.allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    meta.allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
    meta.reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
    meta.mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');   
    meta.sessid=int32(cellfun(@(x) path2sessid(x,'type',opt.type),meta.allpath));
end
reg_map=containers.Map('KeyType','int32','ValueType','any');
for i=1:numel(meta.sessid)
    reg_map(bitshift(meta.sessid(i),16)+int32(meta.allcid(i)))=meta.reg_tree(:,i);
end
reg=arrayfun(@(y) arrayfun(@(x) ...
    reg_map(bitshift(int32(per_type_stats.sess(y)),16)+int32(per_type_stats.sess_suids(y,x))),...
    1:size(per_type_stats.sess_suids,2),'UniformOutput',false), (1:numel(per_type_stats.sess))', 'UniformOutput',false);
end
function out=path2sessid(path,opt)
arguments
    path (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent map

if isempty(map)
    if ispc
        homedir=fullfile('E:','hem','bzdata');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
    else
        homedir=fullfile('/home','hem','datashare');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% path
    end
    if strcmp(opt.type,'neupix')
        allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    else
        fullpath=deblank(h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/path'));
        allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
    end
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(uidx)
        map(upath{i})=uidx(i);
    end
end
out=map(path);

end

function [is_diff,is_same]=diff_at_level(reg,level)
is_diff=false(numel(reg),1);
is_same=false(numel(reg),1);
for i=1:numel(reg)   
    for treedep=level
        if all(cellfun(@(x)ismember(x{1},{'CH','BS'}),reg{i})) && all(cellfun(@(x)~isempty(x{treedep}),reg{i}))           
            if numel(unique(cellfun(@(x)x{treedep},reg{i},'UniformOutput',false)))==1
                is_same(i,1)=true;
            elseif numel(unique(cellfun(@(x)x{treedep},reg{i},'UniformOutput',false)))==numel(reg{i})
                is_diff(i,1)=true;
            end            
        end
    end
end
end