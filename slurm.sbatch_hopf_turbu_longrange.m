#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=Turbulence_SClong
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-301
#SBATCH --output=Turbu_SC%A_%a.out
#SBATCH --error=Turbu_SC%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

NPARCELLS=1000;
NR=400;
NRini=20;
NRfin=380;
NFUTURE=10;
NSUBSIM=100;

load empirical_spacorr_rest.mat;
load schaefercog.mat;
load SCFClongrange.mat;
load hpcdata1003_f_diff_fce.mat;
load enstrophy_cascade_rest.mat
load Yeo7vector.mat;

%%% Networks RSN

ind=find(Yeo7vector==1);
fc(ind,ind)=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
Vise=fc(Isubdi);

ind=find(Yeo7vector==2);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
SomMote=fc(Isubdi);

ind=find(Yeo7vector==3);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
DorsAttne=fc(Isubdi);

ind=find(Yeo7vector==4);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
SalVentAttne=fc(Isubdi);

ind=find(Yeo7vector==5);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
Limbice=fc(Isubdi);

ind=find(Yeo7vector==6);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
Conte=fc(Isubdi);

ind=find(Yeo7vector==7);
fc=fce(ind,ind);
Isubdi = find(tril(ones(length(ind)),-1));
Defaulte=fc(Isubdi);

%%%%

lambda=round(lambda,2);

G_range=0.:0.01:3;

G=G_range(s);

empcorrfcn=corrfcn;
empgrandcorrfcn=grandcorrfcn;
rr=zeros(NPARCELLS,NPARCELLS);
IndLong=[];
for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
        if rr(i,j)>40
            IndLong=[IndLong sub2ind([NPARCELLS NPARCELLS],i,j)];
        end
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end

Isubdiag = find(tril(ones(NPARCELLS),-1));

C=zeros(NPARCELLS,NPARCELLS);

LAMBDA=[0.26 0.22 0.18 0.14 0.10 0.06 0.02];

NLAMBDA=length(LAMBDA);
C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
[aux indsca]=min(abs(LAMBDA-lambda));
ilam=1;
for lambda2=LAMBDA
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C1(ilam,i,j)=exp(-lambda2*rr(i,j));
        end
    end
    ilam=ilam+1;
end

%%%
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

% Parameters HOPF
Tmax=1200;
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

linfunc = @(A, x)(A(1)*x+A(2));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');

% fcsimul_sub=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcn=zeros(NPARCELLS,NR);
% fcsimultot_sub=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcntot=zeros(NPARCELLS,NR);
% fcsimultotrnd_sub=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcntotrnd=zeros(NPARCELLS,NR);
% fcsimulSC_sub=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcnSC=zeros(NPARCELLS,NR);
lam_mean_spatime_enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
ensspasub=zeros(NSUBSIM,NPARCELLS);
ensspasub1=zeros(NSUBSIM,NPARCELLS);
ensspasubtot=zeros(NSUBSIM,NPARCELLS);
ensspasub1tot=zeros(NSUBSIM,NPARCELLS);
ensspasubtotrnd=zeros(NSUBSIM,NPARCELLS);
ensspasub1totrnd=zeros(NSUBSIM,NPARCELLS);
ensspasubSC=zeros(NSUBSIM,NPARCELLS);
ensspasub1SC=zeros(NSUBSIM,NPARCELLS);
mutinf1=zeros(1,NLAMBDA);
mutinf1tot=zeros(1,NLAMBDA);
mutinf1totrnd=zeros(1,NLAMBDA);
mutinf1SC=zeros(1,NLAMBDA);
mutinfo_p=zeros(NSUBSIM,NFUTURE);
mutinfo_ptot=zeros(NSUBSIM,NFUTURE);
mutinfo_ptotrnd=zeros(NSUBSIM,NFUTURE);
mutinfo_pSC=zeros(NSUBSIM,NFUTURE);


err_hete=zeros(1,NSUBSIM);
err_hetetot=zeros(1,NSUBSIM);
err_hetetotrnd=zeros(1,NSUBSIM);
err_heteSC=zeros(1,NSUBSIM);
Inflam=zeros(1,NSUBSIM);
Inflamtot=zeros(1,NSUBSIM);
Inflamtotrnd=zeros(1,NSUBSIM);
InflamSC=zeros(1,NSUBSIM);
Err_Rlam=zeros(1,NSUBSIM);
Err_Rlamtot=zeros(1,NSUBSIM);
Err_Rlamtotrnd=zeros(1,NSUBSIM);
Err_RlamSC=zeros(1,NSUBSIM);
Err_Inflam=zeros(1,NSUBSIM);
Err_Inflamtot=zeros(1,NSUBSIM);
Err_Inflamtotrnd=zeros(1,NSUBSIM);
Err_InflamSC=zeros(1,NSUBSIM);
Rmeta=zeros(1,NSUBSIM);
Rmetatot=zeros(1,NSUBSIM);
Rmetatotrnd=zeros(1,NSUBSIM);
RmetaSC=zeros(1,NSUBSIM);
fcfittlong=zeros(1,NSUBSIM);
fcfittlongtot=zeros(1,NSUBSIM);
fcfittlongtotrnd=zeros(1,NSUBSIM);
fcfittlongSC=zeros(1,NSUBSIM);
fclong=zeros(1,NSUBSIM);
fclongtot=zeros(1,NSUBSIM);
fclongtotrnd=zeros(1,NSUBSIM);
fclongSC=zeros(1,NSUBSIM);
% segint=zeros(1,NSUBSIM);
% seginttot=zeros(1,NSUBSIM);
% seginttotrnd=zeros(1,NSUBSIM);
% segintSC=zeros(1,NSUBSIM);
CorrRSN=zeros(7,NSUBSIM);
CorrRSNtot=zeros(7,NSUBSIM);
CorrRSNtotrnd=zeros(7,NSUBSIM);
CorrRSNSC=zeros(7,NSUBSIM);

G

IClong=find(Clong>0);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        C(i,j)=exp(-lambda*rr(i,j));
    end
    C(i,i)=0;
end
Ctot=C;
Ctot(IClong)=Clong(IClong);
Coriginal=C;

Cxhemi=Ctot(1:500,501:1000);
Cxhemi2=Ctot(501:1000,1:500);
Isubdiaghemi = find(tril(ones(500),-1));
Clonghemi=Clong(1:500,1:500);
Coriginalhemi=Coriginal(1:500,1:500);
Clonghemi2=Clong(501:1000,501:1000);
Coriginalhemi2=Coriginal(501:1000,501:1000);

factor=max(max(C));
C=C/factor*0.2;
Ctot=Ctot/factor*0.2;
SC=SC/factor*0.2;

for sub=1:NSUBSIM
    sub
    Clongrnd=zeros(500,500);
    Clongvec=Clonghemi(Isubdiaghemi);
    indexpermuted=randperm(length(Clongvec));
    Clongvecrnd=Clongvec(indexpermuted);
    Clongrnd(Isubdiaghemi)=Clongvecrnd;
    Clongrnd=Clongrnd+Clongrnd';
    IClongrnd=find(Clongrnd>0);  
    Ctotrndhemi=Coriginalhemi;
    Ctotrndhemi(IClongrnd)=Clongrnd(IClongrnd);
    
    Clongrnd=zeros(500,500);
    Clongvec=Clonghemi2(Isubdiaghemi);
    Clongvecrnd=Clongvec(indexpermuted);
    Clongrnd(Isubdiaghemi)=Clongvecrnd;
    Clongrnd=Clongrnd+Clongrnd';
    IClongrnd=find(Clongrnd>0);  
    Ctotrndhemi2=Coriginalhemi2;
    Ctotrndhemi2(IClongrnd)=Clongrnd(IClongrnd);
    
    Ctotrnd=[Ctotrndhemi Cxhemi;Cxhemi2 Ctotrndhemi2];
    Ctotrnd=Ctotrnd/factor*0.2;
    
    wC = G*C;
    wCtot= G*Ctot;
    wCtotrnd= G*Ctotrnd;
    sumC = repmat(sum(wC,2),1,2);
    sumCtot = repmat(sum(wCtot,2),1,2);
    sumCtotrnd = repmat(sum(wCtotrnd,2),1,2);
    wCSC= G*SC;
    sumCSC = repmat(sum(wCSC,2),1,2);
    
    %% Hopf Simulation
    a=-0.02*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    fcsimul=corrcoef(signal_filt');
    fclong(sub)=nanmean(abs(fcsimul(IndLong)));
    
    for i=1:NPARCELLS
        numind=zeros(1,NR);
        corrfcn_1=zeros(1,NR);
        for j=1:NPARCELLS
            r=rr(i,j);
            index=floor(r/delta)+1;
            if index==NR+1
                index=NR;
            end
            mcc=fcsimul(i,j);
            if ~isnan(mcc)
                corrfcn_1(index)=corrfcn_1(index)+mcc;
                numind(index)=numind(index)+1;
            end
        end
        corrfcn(i,:)=corrfcn_1./numind;
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            ilam=ilam+1;
        end
    end
    Rspatime=squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    Rsub=nanstd(Rspatime(:));
    ensspasub(sub,:)=(nanmean(Rspatime,2))';
    
    ifut=1;
    for lam=1:NFUTURE
        [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(indsca,:,ifut:end)))',(squeeze(lam_mean_spatime_enstrophy(indsca,:,1:end-(ifut-1))))');
        mutinfo_p(sub,ifut)=nanmean(cc(find(pp(:)<0.05)));
        ifut=ifut+1;
    end
        
    ilam=1;
    for lam=LAMBDA
        if ilam>1
            [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))');
            mutinf1(ilam)=nanmean(cc(find(pp(:)<0.05)));
        end
        ilam=ilam+1;
    end
    
    %%% Perturbation
    
    a=-0.02+0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn=0;
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    Rspatime1=zeros(NPARCELLS,Tmax);
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    for i=1:NPARCELLS
        %%% enstrophy
        enstrophy=nansum(repmat(squeeze(C1(indsca,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(indsca,i,:));
        Rspatime1(i,:)=abs(enstrophy);
    end
    ensspasub1(sub,:)=(nanmean(Rspatime1,2))';
    
    %%% Long range
    
    a=-0.02*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    for t=0:dt:2000
        suma = wCtot*z - sumCtot.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    for t=0:dt:((Tmax-1)*TR)
        suma = wCtot*z - sumCtot.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    
    fcsimultot=corrcoef(signal_filt');
    fclongtot(sub)=nanmean(abs(fcsimultot(IndLong)));

    
    for i=1:NPARCELLS
        numind=zeros(1,NR);
        corrfcn_1=zeros(1,NR);
        for j=1:NPARCELLS
            r=rr(i,j);
            index=floor(r/delta)+1;
            if index==NR+1
                index=NR;
            end
            mcc=fcsimultot(i,j);
            if ~isnan(mcc)
                corrfcn_1(index)=corrfcn_1(index)+mcc;
                numind(index)=numind(index)+1;
            end
        end
        corrfcntot(i,:)=corrfcn_1./numind;
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            ilam=ilam+1;
        end
    end
    Rspatime=squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    Rsubtot=nanstd(Rspatime(:));
    ensspasubtot(sub,:)=(nanmean(Rspatime,2))';
    
    ifut=1;
    for lam=1:NFUTURE
        [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(indsca,:,ifut:end)))',(squeeze(lam_mean_spatime_enstrophy(indsca,:,1:end-(ifut-1))))');
        mutinfo_ptot(sub,ifut)=nanmean(cc(find(pp(:)<0.05)));
        ifut=ifut+1;
    end
    
    ilam=1;
    for lam=LAMBDA
        if ilam>1
            [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))');
            mutinf1tot(ilam)=nanmean(cc(find(pp(:)<0.05)));
        end
        ilam=ilam+1;
    end
    
    %%% Perturbation
    
    a=-0.02+0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn=0;
    for t=0:dt:((Tmax-1)*TR)
        suma = wCtot*z - sumCtot.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    Rspatime1=zeros(NPARCELLS,Tmax);
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    for i=1:NPARCELLS
        %%% enstrophy
        enstrophy=nansum(repmat(squeeze(C1(indsca,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(indsca,i,:));
        Rspatime1(i,:)=abs(enstrophy);
    end
    ensspasub1tot(sub,:)=(nanmean(Rspatime1,2))';
    
    %%%%%%
    % long range randomized
    
    
    a=-0.02*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    for t=0:dt:2000
        suma = wCtotrnd*z - sumCtotrnd.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    for t=0:dt:((Tmax-1)*TR)
        suma = wCtotrnd*z - sumCtotrnd.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    
    fcsimultotrnd=corrcoef(signal_filt');
    fclongtotrnd(sub)=nanmean(abs(fcsimultotrnd(IndLong)));

    
    for i=1:NPARCELLS
        numind=zeros(1,NR);
        corrfcn_1=zeros(1,NR);
        for j=1:NPARCELLS
            r=rr(i,j);
            index=floor(r/delta)+1;
            if index==NR+1
                index=NR;
            end
            mcc=fcsimultotrnd(i,j);
            if ~isnan(mcc)
                corrfcn_1(index)=corrfcn_1(index)+mcc;
                numind(index)=numind(index)+1;
            end
        end
        corrfcntotrnd(i,:)=corrfcn_1./numind;
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            ilam=ilam+1;
        end
    end
    Rspatime=squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    Rsubtotrnd=nanstd(Rspatime(:));
    ensspasubtotrnd(sub,:)=(nanmean(Rspatime,2))';
    
    ifut=1;
    for lam=1:NFUTURE
        [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(indsca,:,ifut:end)))',(squeeze(lam_mean_spatime_enstrophy(indsca,:,1:end-(ifut-1))))');
        mutinfo_ptotrnd(sub,ifut)=nanmean(cc(find(pp(:)<0.05)));
        ifut=ifut+1;
    end
    
    ilam=1;
    for lam=LAMBDA
        if ilam>1
            [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))');
            mutinf1totrnd(ilam)=nanmean(cc(find(pp(:)<0.05)));
        end
        ilam=ilam+1;
    end
    
    %%% Perturbation
    
    a=-0.02+0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn=0;
    for t=0:dt:((Tmax-1)*TR)
        suma = wCtotrnd*z - sumCtotrnd.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    Rspatime1=zeros(NPARCELLS,Tmax);
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    for i=1:NPARCELLS
        %%% enstrophy
        enstrophy=nansum(repmat(squeeze(C1(indsca,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(indsca,i,:));
        Rspatime1(i,:)=abs(enstrophy);
    end
    ensspasub1totrnd(sub,:)=(nanmean(Rspatime1,2))';    
    
    
    %%%%%%%
    
    %% Hopf Simulation
    a=-0.02*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wCSC*z - sumCSC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wCSC*z - sumCSC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    fcsimulSC=corrcoef(signal_filt');
    fclongSC(sub)=nanmean(abs(fcsimulSC(IndLong)));
    
    for i=1:NPARCELLS
        numind=zeros(1,NR);
        corrfcn_1=zeros(1,NR);
        for j=1:NPARCELLS
            r=rr(i,j);
            index=floor(r/delta)+1;
            if index==NR+1
                index=NR;
            end
            mcc=fcsimulSC(i,j);
            if ~isnan(mcc)
                corrfcn_1(index)=corrfcn_1(index)+mcc;
                numind(index)=numind(index)+1;
            end
        end
        corrfcnSC(i,:)=corrfcn_1./numind;
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            ilam=ilam+1;
        end
    end
    Rspatime=squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    RsubSC=nanstd(Rspatime(:));
    ensspasubSC(sub,:)=(nanmean(Rspatime,2))';
    
    ifut=1;
    for lam=1:NFUTURE
        [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(indsca,:,ifut:end)))',(squeeze(lam_mean_spatime_enstrophy(indsca,:,1:end-(ifut-1))))');
        mutinfo_pSC(sub,ifut)=nanmean(cc(find(pp(:)<0.05)));
        ifut=ifut+1;
    end
    
    ilam=1;
    for lam=LAMBDA
        if ilam>1
            [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))');
            mutinf1SC(ilam)=nanmean(cc(find(pp(:)<0.05)));
        end
        ilam=ilam+1;
    end
    
    %%% Perturbation
    
    a=-0.02+0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn=0;
    for t=0:dt:((Tmax-1)*TR)
        suma = wCSC*z - sumCSC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    Rspatime1=zeros(NPARCELLS,Tmax);
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    for i=1:NPARCELLS
        %%% enstrophy
        enstrophy=nansum(repmat(squeeze(C1(indsca,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(indsca,i,:));
        Rspatime1(i,:)=abs(enstrophy);
    end
    ensspasub1SC(sub,:)=(nanmean(Rspatime1,2))';
    
    %%%%%%%%%%%
    %% Observables
    %%%%%%%%%%%
    
    for i=1:NPARCELLS
        for k=NRini:NRfin
            err11(k)=(corrfcn(i,k)-empcorrfcn(i,k))^2;
        end
        err1(i)=(nanmean(err11(NRini:NRfin)));
    end
    err_hete(sub)=sqrt(nanmean(err1));
    
    for i=1:NPARCELLS
        for k=NRini:NRfin
            err11(k)=(corrfcntot(i,k)-empcorrfcn(i,k))^2;
        end
        err1(i)=(nanmean(err11(NRini:NRfin)));
    end
    err_hetetot(sub)=sqrt(nanmean(err1));

    for i=1:NPARCELLS
        for k=NRini:NRfin
            err11(k)=(corrfcntotrnd(i,k)-empcorrfcn(i,k))^2;
        end
        err1(i)=(nanmean(err11(NRini:NRfin)));
    end
    err_hetetotrnd(sub)=sqrt(nanmean(err1));
    
    for i=1:NPARCELLS
        for k=NRini:NRfin
            err11(k)=(corrfcnSC(i,k)-empcorrfcn(i,k))^2;
        end
        err1(i)=(nanmean(err11(NRini:NRfin)));
    end
    err_heteSC(sub)=sqrt(nanmean(err1));
    
    %%
    Inflam2=mutinf1(2:end);
    InfoCascade(sub,:)=Inflam2;
    Inflam(sub)=nanmean(Inflam2);
    Err_Rlam(sub)=sqrt((Rsub-Rmeta_emp)^2);
    Err_Inflam(sub)=sqrt(nanmean((Inflam2-Inflam_emp).^2));
    
    Inflam2=mutinf1tot(2:end);
    InfoCascadetot(sub,:)=Inflam2;
    Inflamtot(sub)=nanmean(Inflam2);
    Err_Rlamtot(sub)=sqrt((Rsubtot-Rmeta_emp)^2);
    Err_Inflamtot(sub)=sqrt(nanmean((Inflam2-Inflam_emp).^2));

    Inflam2=mutinf1totrnd(2:end);
    InfoCascadetotrnd(sub,:)=Inflam2;
    Inflamtotrnd(sub)=nanmean(Inflam2);
    Err_Rlamtotrnd(sub)=sqrt((Rsubtotrnd-Rmeta_emp)^2);
    Err_Inflamtotrnd(sub)=sqrt(nanmean((Inflam2-Inflam_emp).^2));
    
    Inflam2=mutinf1SC(2:end);
    InfoCascadeSC(sub,:)=Inflam2;
    InflamSC(sub)=nanmean(Inflam2);
    Err_RlamSC(sub)=sqrt((RsubSC-Rmeta_emp)^2);
    Err_InflamSC(sub)=sqrt(nanmean((Inflam2-Inflam_emp).^2));
    
    %%
    Rmeta(sub)=Rsub;
    Rmetatot(sub)=Rsubtot;
    Rmetatotrnd(sub)=Rsubtotrnd;
    RmetaSC(sub)=RsubSC;
    
    %%
    k=1;
    for i=IndLong
        errfcfl(k)=(fce(i)-fcsimul(i))^2;
        k=k+1;
    end
    fcfittlong(sub)=sqrt(nanmean(errfcfl));
    
    k=1;
    for i=IndLong
        errfcfl(k)=(fce(i)-fcsimultot(i))^2;
        k=k+1;
    end
    fcfittlongtot(sub)=sqrt(nanmean(errfcfl));
    
    k=1;
    for i=IndLong
        errfcfl(k)=(fce(i)-fcsimultotrnd(i))^2;
        k=k+1;
    end
    fcfittlongtotrnd(sub)=sqrt(nanmean(errfcfl));
    
    k=1;
    for i=IndLong
        errfcfl(k)=(fce(i)-fcsimulSC(i))^2;
        k=k+1;
    end
    fcfittlongSC(sub)=sqrt(nanmean(errfcfl));
    
    %%
%     
%     [MM QQ]=community_louvain(fcsimul,[],[],'negative_sym');
%     modularity=QQ;
%     integration=nanmean(nanmean(abs(fcsimul)-eye(NPARCELLS)));
%     segint(sub)=modularity*integration;
%     
%     [MM QQ]=community_louvain(fcsimultot,[],[],'negative_sym');
%     modularity=QQ;
%     integration=nanmean(nanmean(abs(fcsimultot)-eye(NPARCELLS)));
%     seginttot(sub)=modularity*integration;
% 
%     [MM QQ]=community_louvain(fcsimultotrnd,[],[],'negative_sym');
%     modularity=QQ;
%     integration=nanmean(nanmean(abs(fcsimultotrnd)-eye(NPARCELLS)));
%     seginttotrnd(sub)=modularity*integration;
%     
%     [MM QQ]=community_louvain(fcsimulSC,[],[],'negative_sym');
%     modularity=QQ;
%     integration=nanmean(nanmean(abs(fcsimulSC)-eye(NPARCELLS)));
%     segintSC(sub)=modularity*integration;
%     
    %%
    
    ind=find(Yeo7vector==1);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    Vis1=fc1(Isubdi);
    Vis2=fc2(Isubdi);
    Vis22=fc22(Isubdi);
    Vis3=fc3(Isubdi);
    CorrRSN(1,sub)=sqrt(nanmean((Vise-Vis1).^2));
    CorrRSNtot(1,sub)=sqrt(nanmean((Vise-Vis2).^2));
    CorrRSNtotrnd(1,sub)=sqrt(nanmean((Vise-Vis22).^2));
    CorrRSNSC(1,sub)=sqrt(nanmean((Vise-Vis3).^2));
    
    ind=find(Yeo7vector==2);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    SomMot1=fc1(Isubdi);
    SomMot2=fc2(Isubdi);
    SomMot22=fc22(Isubdi);
    SomMot3=fc3(Isubdi);
    CorrRSN(2,sub)=sqrt(nanmean((SomMote-SomMot1).^2));
    CorrRSNtot(2,sub)=sqrt(nanmean((SomMote-SomMot2).^2));
    CorrRSNtotrnd(2,sub)=sqrt(nanmean((SomMote-SomMot22).^2));
    CorrRSNSC(2,sub)=sqrt(nanmean((SomMote-SomMot3).^2));
    
    ind=find(Yeo7vector==3);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    DorsAttn1=fc1(Isubdi);
    DorsAttn2=fc2(Isubdi);
    DorsAttn22=fc22(Isubdi);
    DorsAttn3=fc3(Isubdi);
    CorrRSN(3,sub)=sqrt(nanmean((DorsAttne-DorsAttn1).^2));
    CorrRSNtot(3,sub)=sqrt(nanmean((DorsAttne-DorsAttn2).^2));
    CorrRSNtotrnd(3,sub)=sqrt(nanmean((DorsAttne-DorsAttn22).^2));
    CorrRSNSC(3,sub)=sqrt(nanmean((DorsAttne-DorsAttn3).^2));
    
    ind=find(Yeo7vector==4);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    SalVentAttn1=fc1(Isubdi);
    SalVentAttn2=fc2(Isubdi);
    SalVentAttn22=fc22(Isubdi);
    SalVentAttn3=fc3(Isubdi);
    CorrRSN(4,sub)=sqrt(nanmean((SalVentAttne-SalVentAttn1).^2));
    CorrRSNtot(4,sub)=sqrt(nanmean((SalVentAttne-SalVentAttn2).^2));
    CorrRSNtotrnd(4,sub)=sqrt(nanmean((SalVentAttne-SalVentAttn22).^2));
    CorrRSNSC(4,sub)=sqrt(nanmean((SalVentAttne-SalVentAttn3).^2));
    
    ind=find(Yeo7vector==5);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    Limbic1=fc1(Isubdi);
    Limbic2=fc2(Isubdi);
    Limbic22=fc22(Isubdi);
    Limbic3=fc3(Isubdi);
    CorrRSN(5,sub)=sqrt(nanmean((Limbice-Limbic1).^2));
    CorrRSNtot(5,sub)=sqrt(nanmean((Limbice-Limbic2).^2));
    CorrRSNtotrnd(5,sub)=sqrt(nanmean((Limbice-Limbic22).^2));
    CorrRSNSC(5,sub)=sqrt(nanmean((Limbice-Limbic3).^2));
    
    ind=find(Yeo7vector==6);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    Cont1=fc1(Isubdi);
    Cont2=fc2(Isubdi);
    Cont22=fc22(Isubdi);
    Cont3=fc3(Isubdi);
    CorrRSN(6,sub)=sqrt(nanmean((Conte-Cont1).^2));
    CorrRSNtot(6,sub)=sqrt(nanmean((Conte-Cont2).^2));
    CorrRSNtotrnd(6,sub)=sqrt(nanmean((Conte-Cont22).^2));
    CorrRSNSC(6,sub)=sqrt(nanmean((Conte-Cont3).^2));
    
    ind=find(Yeo7vector==7);
    fc1=fcsimul(ind,ind);
    fc2=fcsimultot(ind,ind);
    fc22=fcsimultotrnd(ind,ind);
    fc3=fcsimulSC(ind,ind);
    Isubdi = find(tril(ones(length(ind)),-1));
    Default1=fc1(Isubdi);
    Default2=fc2(Isubdi);
    Default22=fc22(Isubdi);
    Default3=fc3(Isubdi);
    CorrRSN(7,sub)=sqrt(nanmean((Defaulte-Default1).^2));
    CorrRSNtot(7,sub)=sqrt(nanmean((Defaulte-Default2).^2));
    CorrRSNtotrnd(7,sub)=sqrt(nanmean((Defaulte-Default22).^2));
    CorrRSNSC(7,sub)=sqrt(nanmean((Defaulte-Default3).^2));
    %%%%%%%%%%
    
%     fcsimul_sub(sub,:,:)=fcsimul;
%     fcsimultot_sub(sub,:,:)=fcsimultot;
%     fcsimultotrnd_sub(sub,:,:)=fcsimultotrnd;
%     fcsimulSC_sub(sub,:,:)=fcsimulSC;
end

% fcsimul_avg=nanmean(fcsimul_sub,1);
% fcsimultot_avg=nanmean(fcsimultot_sub,1);
% fcsimultotrnd_avg=nanmean(fcsimultotrnd_sub,1);
% fcsimulSC_avg=nanmean(fcsimulSC_sub,1);

infocapacity=nanmean(nanstd(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)));
susceptibility=nanmean(nanmean(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)));

infocapacitytot=nanmean(nanstd(ensspasub1tot-ones(NSUBSIM,1)*nanmean(ensspasubtot)));
susceptibilitytot=nanmean(nanmean(ensspasub1tot-ones(NSUBSIM,1)*nanmean(ensspasubtot)));

infocapacitytotrnd=nanmean(nanstd(ensspasub1totrnd-ones(NSUBSIM,1)*nanmean(ensspasubtotrnd)));
susceptibilitytotrnd=nanmean(nanmean(ensspasub1totrnd-ones(NSUBSIM,1)*nanmean(ensspasubtotrnd)));

infocapacitySC=nanmean(nanstd(ensspasub1SC-ones(NSUBSIM,1)*nanmean(ensspasubSC)));
susceptibilitySC=nanmean(nanmean(ensspasub1SC-ones(NSUBSIM,1)*nanmean(ensspasubSC)));

errinfocapacity=nanstd(nanstd(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)));
errsusceptibility=nanstd(nanmean(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)));

errinfocapacitytot=nanstd(nanstd(ensspasub1tot-ones(NSUBSIM,1)*nanmean(ensspasubtot)));
errsusceptibilitytot=nanstd(nanmean(ensspasub1tot-ones(NSUBSIM,1)*nanmean(ensspasubtot)));

errinfocapacitytotrnd=nanstd(nanstd(ensspasub1totrnd-ones(NSUBSIM,1)*nanmean(ensspasubtotrnd)));
errsusceptibilitytotrnd=nanstd(nanmean(ensspasub1totrnd-ones(NSUBSIM,1)*nanmean(ensspasubtotrnd)));

errinfocapacitySC=nanstd(nanstd(ensspasub1SC-ones(NSUBSIM,1)*nanmean(ensspasubSC)));
errsusceptibilitySC=nanstd(nanmean(ensspasub1SC-ones(NSUBSIM,1)*nanmean(ensspasubSC)));

save(sprintf('WG_%03d.mat',s),'fclong','fclongtot','fclongtotrnd','fclongSC','InfoCascade','InfoCascadetot','InfoCascadeSC','InfoCascadetotrnd','mutinfo_p','mutinfo_ptot','mutinfo_ptotrnd','mutinfo_pSC','CorrRSN','CorrRSNtot','CorrRSNtotrnd','CorrRSNSC','Inflam','Err_Rlam','Err_Inflam','Inflamtot','Err_Rlamtot','Err_Inflamtot','Inflamtotrnd','Err_Rlamtotrnd','Err_Inflamtotrnd','InflamSC','Err_RlamSC','Err_InflamSC','infocapacity','infocapacitytot','infocapacitytotrnd','infocapacitySC','susceptibility','susceptibilitytot','susceptibilitytotrnd','susceptibilitySC','errinfocapacity','errinfocapacitytot','errinfocapacitytotrnd','errinfocapacitySC','errsusceptibility','errsusceptibilitytot','errsusceptibilitytotrnd','errsusceptibilitySC','Rmeta','Rmetatot','Rmetatotrnd','RmetaSC','fcfittlong','fcfittlongtot','fcfittlongtotrnd','fcfittlongSC','err_hete','err_hetetot','err_hetetotrnd','err_heteSC');

EOF