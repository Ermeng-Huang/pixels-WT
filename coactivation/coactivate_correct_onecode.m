clc
clear
warning('off');
for i=1:19 %7 non-memory
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
    opt.tsbin_size (1,1) double = 6000
    opt.postspike (1,1) logical = true
    opt.fc_effi (1,1) logical = false
    opt.fc_prob (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'on'
    opt.assembly (1,1) logical = true
end
% sig=load_sig_pair(opt);
load('sig.mat','sig')
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
avail=[];
sidx=1;
sess_suids=sig.suid(idces,:);
for i=reshape(idces,1,[])    
    if nnz(sess_suids(:,2)==sig.suid(i,2))<=1 || nnz(sess_suids(:,2)==sig.suid(i,2))>5
        continue
    end
    if sidx>1 && nnz(suid(:,3)==sig.suid(i,2))
         continue
    end
    pair=sess_suids(sess_suids(:,2)==sig.suid(i,2),:);    
    p=nchoosek(1:size(pair,1),2);
    for pp=1:size(p,1)
    fprintf('%d of %d\n',sidx,dim);
    [postspk{sidx},avail(sidx,:)]=history_coeff_coactivate_trial(sig.sess(i),pair(p(pp,:),:),...
        'tsbin_size',opt.tsbin_size,...
        'postspike',opt.postspike,...
        'fc_effi',opt.fc_effi,...
        'fc_prob',opt.fc_prob,...
        'type',opt.type,...
        'laser',opt.laser,...
        'assembly',opt.assembly,...
        'sess',sess);  
 
   
    suid(sidx,:)=[pair(p(pp,:),1);unique(pair(p(pp,:),2))]';
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


postspk(~all(avail,2))=[];
suid(~all(avail,2),:)=[];
suid_sess(~all(avail,2),:)=[];
if ~isempty(postspk)
    
    save(fullfile('F:','neupix','STP','coactivate','0704-trial',sprintf('%s_stp_%s_%d_%d.mat',...
        opt.prefix,mtype,sess,opt.tsbin_size)),...
        'postspk','suid','suid_sess','mtype');
end
end

function [out,out_avail]=history_coeff_coactivate_trial(sessid,suid,opt)
arguments
    sessid (1,1) int32 
    suid (2,2) int32        
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

bitmask=2.^(0:1)';
out_avail=false;

trial.correct=spike.trials(:,8)==6&spike.trials(:,9)==1&spike.trials(:,10)==1;
trial.error=spike.trials(:,8)==6&spike.trials(:,10)==0;
trial.min=min(nnz(trial.correct),nnz(trial.error));
if trial.min==0
   out='error';
   return
end

[avail(1),tspre1]=pre_process(suid(1,1),spike.spkTS,spike.spkID,spike.trials,opt);
[avail(2),tspre2]=pre_process(suid(2,1),spike.spkTS,spike.spkID,spike.trials,opt);
[avail(3),tspost]=pre_process(suid(1,2),spike.spkTS,spike.spkID,spike.trials,opt);
if ~any(avail)
    return
end

tmin=-3;
tmax=11;
bin_size=opt.tsbin_size/30000;
histpre1=cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspre1,'UniformOutput',false);
histpre2=cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspre2,'UniformOutput',false);
histpost=cellfun(@(x)histcounts(x,tmin:bin_size:tmax)>0,tspost,'UniformOutput',false);
if nnz(cell2mat(histpost))== 0            
    out='error';    
    return    
end
histpost1=cellfun(@(x)x(x>10),(cellfun(@(x)find(x~=0),histpost,'UniformOutput',false)),'UniformOutput',false);
bad=find(cell2mat(cellfun(@(x)isempty(x),histpost1,'UniformOutput',false)));
if ~isempty(bad)
    trial.correct(bad)=[];
    trial.error(bad)=[];
    trial.min=min(nnz(trial.correct),nnz(trial.error));
    if trial.min==0
        out='error';
        return
    end
   histpre1(bad)=[];histpre2(bad)=[];histpost1(bad)=[];
end  
    
n=cellfun(@(x)size(x,2),histpost1,'UniformOutput',false); 

for i =1:size(histpost1,2)
    stp{i}=cell2mat(arrayfun(@(x)[histpre1{i}(x-10:x-1)',histpre2{i}(x-10:x-1)']*bitmask,histpost1{i},'UniformOutput',false));
end

shuffle=["correct","error"];
    
for Isshuffle=shuffle
    if nnz(trial.(Isshuffle))==trial.min
        temp.(Isshuffle).stp=stp(trial.(Isshuffle));    
        temp.(Isshuffle).num=cell2mat(arrayfun(@(x)(n{x}),find(trial.(Isshuffle)),'UniformOutput',false));
        for i=1:10            
            out.(Isshuffle)(:,i)=mean(reshape(cell2mat(cellfun(@(x)nnz(x(i,:)==3),temp.(Isshuffle).stp,'UniformOutput',false))./temp.(Isshuffle).num,[],1));
        end 
        out_avail(1)=true;
    else
        trial.good=find(trial.(Isshuffle));  
        id=cell2mat(arrayfun(@(x)trial.good(randperm(length(trial.good),trial.min)),1:500,'UniformOutput',false));
        temp.(Isshuffle).stp=arrayfun(@(x)stp{x},id,'UniformOutput',false); 
        temp.(Isshuffle).num=cell2mat(arrayfun(@(x)(n{x}),id,'UniformOutput',false));
        for i=1:10            
            out.(Isshuffle)(:,i)=mean(reshape(cell2mat(cellfun(@(x)nnz(x(i,:)==3),temp.(Isshuffle).stp,'UniformOutput',false))./temp.(Isshuffle).num,[],1));
        end 
        out_avail(2)=true;
    end    
    
    
end
end

function cal(histpre1,histpre2,histpost)

end
function [avail,out]=pre_process(id,spkTS,spkId,trials,~)
timlim=[-4,14];
sps=30000;
addpath('D:/code/fieldtrip-20200320')
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

%TODO merge with reg_conn_bz script

function [spkID,spkTS,trials,SU_id,folder]=getSPKID_TS(fidx)
arguments
    fidx (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(fidx,1)}
end

homedir=fullfile('F:','neupix','SPKINFO'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [out,homedir]=sessid2path(sessid,opt)
arguments
    sessid (1,1) double {mustBeInteger,mustBePositive}
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent map
homedir=fullfile('E:','hem','bzdata');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(map)
    
    if strcmp(opt.type,'neupix')
        allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    else
        fullpath=deblank(h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/path'));
        allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
    end
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','int32','ValueType','char');
    for i=1:numel(uidx)
        map(uidx(i))=upath{i};
    end
end

out=map(sessid);

end

