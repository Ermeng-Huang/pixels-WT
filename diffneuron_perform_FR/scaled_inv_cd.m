load('F:\neupix\per_sec\scale_invarient.mat','scale','invarient','corr')


typesel=find(any(invarient,2)' & strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),''));
homedir=ephys.util.getHomedir('type','raw');
phase=["early_delay","late_delay"];
for i=1:numel(phase)
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
    fidx=fidx+1;
end
%% CD
cdMat=FR_S1-FR_S2;
cdDelay=mean(cdMat,2);
cdDelay=cdDelay/norm(cdDelay);
[qq,rr]=qr(cdDelay);
cdfree=[FR_S1,FR_S2].'*qq(:,2:end);
[~,score,~]=pca(cdfree,'NumComponents',20);
score=mat2cell(score,ones(6,1)*size(score,1)/6,20);
% [CD_s1off,CD_s2off,CD_s1on,CD_s2on]=deal(score{1},score{2},score{3},score{4});

%% PCA
% [~,score,~]=pca(FR_All','NumComponents',20);
% score=mat2cell(score,ones(4,1)*size(score,1)/4,20);
% [PCA_s1off,PCA_s2off,PCA_s1on,PCA_s2on]=deal(score{1},score{2},score{3},score{4});

%%
fh=figure('Color','w','Position',[100,100,450,300]);
c={'k','r','b'};
% subplot(1,2,1)
view(160,-10)
hold on
for i=1:3
    plot3(score{i}(:,1),score{i}(:,2),FR_S1(:,32*(i-1)+(1:32)).'*cdDelay,'-','Color',c{i},'LineWidth',1.5);
    plot3(score{i+3}(:,1),score{i+3}(:,2),FR_S2(:,32*(i-1)+(1:32)).'*cdDelay,'--','Color',c{i},'LineWidth',1.5);
end

xlabel('PC1')
ylabel('PC2')
zlabel('CD')