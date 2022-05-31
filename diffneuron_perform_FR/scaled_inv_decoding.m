clc
clear
load('F:\neupix\per_sec\scale_invarient.mat','scale','invarient','corr')
meta_6=ephys.util.load_meta('delay',6);
typesel_scale=find(any(scale,2)' & strcmp(meta_6.reg_tree(2,:),'CTX') & ~strcmp(meta_6.reg_tree(5,:),''));
dec_scale=Decoding(data_dec(typesel_scale));
typesel_invariant=find(any(invarient,2)' & strcmp(meta_6.reg_tree(2,:),'CTX') & ~strcmp(meta_6.reg_tree(5,:),''));
dec_invariant=Decoding(data_dec(typesel_invariant(randperm(length(typesel_invariant),length(typesel_scale)))));
%%
fh=figure('Color','w','Position',[100,100,250,200]);

hold on
bar(1:4,cellfun(@mean,[dec_scale.cvcorr(1:2),dec_invariant.cvcorr(1:2)]),'w')
errorbar(1:4,cellfun(@mean,[dec_scale.cvcorr(1:2),dec_invariant.cvcorr(1:2)]),cellfun(@(x)std(x)/sqrt(length(x)),[dec_scale.cvcorr(1:2),dec_invariant.cvcorr(1:2)])...
    ,'lineStyle','none','color','k')
set(gca,'XTick',1:4,'XtickLabel',{'scaled neurons-ED','scaled neurons-LD'...
    ,'invariant neurons-ED','invariant neurons-LD'},'XTickLabelRotation',45)
ylabel('decoding accuracy')
%% function
function out=Decoding(decdata,opt)
%%%
% Input:
%     decdata is structure, including two structures (train and test, each including two 4-D matrix s1 and s2)
%     decdata.train.s1 = (cellID * trlN * bins * rpts) matrix
%%%
arguments
    decdata (1,1) struct
    opt.decoder (1,:) char {mustBeMember(opt.decoder,{'SVM','LDA','NB'})} ='SVM';
end

bins=size(decdata.train.s1,3);
trlN=size(decdata.train.s1,2);
rpts=size(decdata.train.s1,4);


y=[zeros(trlN,1);ones(trlN,1)];
out=struct();
cv=cvpartition(trlN,'KFold',5);

cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),[1:5],'UniformOutput',false));
cv_testID=cell2mat(arrayfun(@(x)test(cv,x),[1:5],'UniformOutput',false));
% s1kf=arrayfun(@(x)decdata.train.s1(:,x,1,1),cv_trainingID,'UniformOutput',false);
% s1kf_1=reshape(s1kf(cv_trainingID),[],10);
switch(nnz(contains(fieldnames(decdata),'test')))
    case 1 
        for rpt=1:rpts
            cvresult=cell(1,bins);
            shufresult=cell(1,bins);
            weightresult=cell(1,bins);
            biasresult=cell(1,bins);
            fprintf('rpts %d of %d\n',rpt,rpts);
            for bin=1:bins                
                for kf=1:cv.NumTestSets
                    s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                    s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
%                     varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
%                     s1kf=s1kf(varsel,:);
%                     s2kf=s2kf(varsel,:);                    
                    Xkf=cat(2,s1kf,s2kf)';
                    ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                    
                    s1Tkf=decdata.test.s1(:,cv_testID(:,kf),bin,rpt);
                    s2Tkf=decdata.test.s2(:,cv_testID(:,kf),bin,rpt);
                    XTkf=cat(2,s1Tkf,s2Tkf)';
                    yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
                                       
                    yshufkf=ykf(randperm(numel(ykf)));
                    if strcmp(opt.decoder,'SVM')
                        SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                        SVMMshuf=fitcsvm(Xkf,yshufkf,'KernelFunction','linear','Standardize',true);
                    elseif strcmp(opt.decoder,'LDA')
                        LDAM=fitcdiscr(Xkf,ykf);
                    elseif strcmp(opt.decoder,'NB')
                        NBM=fitcnb(Xkf,ykf);
                    end
                    
                    cvresult{:,bin}=cat(1,cvresult{:,bin},SVMM.predict(XTkf)==yTkf);
                    shufresult{:,bin}=cat(1,shufresult{:,bin},SVMMshuf.predict(XTkf)==yTkf);  
                    weightresult{:,bin}=cat(2,weightresult{:,bin},SVMM.Beta);
                    biasresult{:,bin}=cat(2,biasresult{:,bin},SVMM.Bias); 
                end                
            end
            cvcorr{rpt}=cellfun(@(x)nnz(x)/length(x),cvresult);
            shufcorr{rpt}=cellfun(@(x)nnz(x)/length(x),shufresult);
            weight{rpt}=cellfun(@(x)mean(x,2),weightresult,'UniformOutput',false);
            bias{rpt}=cellfun(@(x)mean(x,2),biasresult,'UniformOutput',false);
        end
    case 2   
        out.cvcorr2=cell(1,bins);
        for rpt=1:rpts
            fprintf('rpts %d of %d\n',rpt,rpts);
            for bin=1:bins                
                for kf=1:cv.NumTestSets
                    s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                    s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
                    varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
                    s1kf=s1kf(varsel,:);
                    s2kf=s2kf(varsel,:);                    
                    Xkf=cat(2,s1kf,s2kf)';
                    ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                    
                    s1Tkf=decdata.test1.s1(varsel,cv_testID(:,kf),bin,rpt);
                    s2Tkf=decdata.test1.s2(varsel,cv_testID(:,kf),bin,rpt);
                    XTkf1=cat(2,s1Tkf,s2Tkf)';
                    yTkf1=y([cv_testID(:,kf);cv_testID(:,kf)]);
                    
                    
                    XTkf2=[decdata.test2.s1(varsel,:,bin,rpt),decdata.test2.s2(varsel,:,bin,rpt)]';
                    yTkf2=[zeros(size(decdata.test2.s1,2),1);ones(size(decdata.test2.s1,2),1)];
                    
                    if strcmp(opt.decoder,'SVM')
                        SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                        modelPredict=SVMM.predict(XTkf1);
                    elseif strcmp(opt.decoder,'LDA')
                        LDAM=fitcdiscr(Xkf,ykf);
                        modelPredict=LDAM.predict(XTkf1);
                    elseif strcmp(opt.decoder,'NB')
                        NBM=fitcnb(Xkf,ykf);
                        modelPredict=NBM.predict(XTkf1);
                    end
                    
                    out.cvcorr{:,bin}=cat(1,out.cvcorr{:,bin},modelPredict==yTkf1);
                    out.cvcorr2{:,bin}=cat(1,out.cvcorr2{:,bin},SVMM.predict(XTkf2)==yTkf2);
                    
                    for i=1:200
                        yshufTkf=yTkf1(randperm(numel(yTkf1)));
                        cvshufresult=modelPredict==yshufTkf;
                        out.shufcorr{:,bin}=cat(1,out.shufcorr{:,bin},cvshufresult);
                    end                    
                end               
               
            end
        end
end
for bin=1:bins    
    out.cvcorr{1,bin}=cellfun(@(x)x(1,bin),cvcorr);
    out.shufcorr{1,bin}=cellfun(@(x)x(1,bin),shufcorr);
    out.weight{1,bin}=cell2mat(cellfun(@(x)x(1,bin),weight));
    out.bias{1,bin}=cell2mat(cellfun(@(x)x(1,bin),bias));
    
end
out.P=cellfun(@(x,y)ranksum(x,y),out.cvcorr,out.shufcorr);
% out.p_permutation=cellfun(@(x,y)statistics.permutationTest(x,y,100),out.cvcorr,out.shufcorr);

end

function out=data_dec(typesel)
meta_6=ephys.util.load_meta('delay',6);
phase=["early_delay","late_delay"];
phase_bin=[1,2.5;2.5,4];

for i=1:numel(phase)
    idx=1;
    for ii=typesel        
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
                
        delaymm=mean([mean(cs1),mean(cs2)]);
        delaystd=std([cs1;cs2]);
        if delaystd==0, continue;  end
        
        cs1=(cs1-delaymm)./delaystd;
        cs2=(cs2-delaymm)./delaystd;
        es1=(es1-delaymm)./delaystd;
        es2=(es2-delaymm)./delaystd;  
        
        for rpt=1:100
            out.train.s1(idx,:,i,rpt)=cs1(randperm(length(cs1),10));
            out.train.s2(idx,:,i,rpt)=cs2(randperm(length(cs2),10));
            out.train.s1(idx,:,i+2,rpt)=cs1(randperm(length(cs1),10));
            out.train.s2(idx,:,i+2,rpt)=cs2(randperm(length(cs2),10));
            
        end
        idx=idx+1;
    end
    
end
out.test=out.train;
end