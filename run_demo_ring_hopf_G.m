clear all;

NPARCELLS=1000;

for i=1:NPARCELLS
    for j=1:NPARCELLS
        d1=abs(j-i);
        d2=NPARCELLS-d1;
        rr(i,j)=min(d1,d2)/10;
    end
end

maxrr=max(max(rr));

LAMBDA=[4 2 1 0.5 0.25 0.125 0.0625 0.0312];
NLAMBDA=length(LAMBDA);
C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
C=zeros(NPARCELLS,NPARCELLS);
C2=zeros(NPARCELLS,NPARCELLS);

lambda=1;  %%1 works!!
[aux Nscale]=find(LAMBDA==lambda);
%% repeat for lambda=4;
for i=1:NPARCELLS
    for j=1:NPARCELLS
        C(i,j)=exp(-lambda*rr(i,j));
    end
    C(i,i)=0;
end

ilam=1;
for lambda2=LAMBDA
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C1(ilam,i,j)=exp(-lambda2*rr(i,j));
        end
        C1(ilam,i,i)=1;
    end
    ilam=ilam+1;
end

%%%%

% Parameters HOPF
Tmax=1000;
TR=1;
f_diff=0.025*ones(1,NPARCELLS);
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

G_range=0.:0.001:0.1;

beta=0.05;

lam_mean_spatime_enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
lam_mean_spatime_enstrophy2=zeros(NLAMBDA,NPARCELLS,Tmax);
Phases=zeros(NPARCELLS,Tmax);
Phases2=zeros(NPARCELLS,Tmax);
mutinf1_1=zeros(1,NLAMBDA);
mutinf12_1=zeros(1,NLAMBDA);

C2=C;
for i=1:NPARCELLS
    for j=1:NPARCELLS
        if rand<beta
            out=randperm(NPARCELLS);
            k=out(1);
            l=out(end);
            C2(l,k)=0.25;
        end
    end
end

for i=1:NPARCELLS
    C2(i,i)=0;
end

ii=1;

for G=G_range
    for sub=1:100
        wC = G*squeeze(C);
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        
        %% Hopf Simulation
        a=0.1*ones(NPARCELLS,2);
        xs=zeros(Tmax,NPARCELLS);
        z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 2000 time steps
        for t=0:dt:2000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        end
        % actual modeling
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
            Xanalytic = hilbert(demean(ts(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
        
        ccnn=0;
        cclong=0;
        for i=1:NPARCELLS
            for j=1:NPARCELLS
                if rr(i,j)>0.8*maxrr
                    cclong=cclong+(corr2(ts(i,:),ts(j,:)));
                    ccnn=ccnn+1;
                end
            end
        end
        FClarge(sub)=cclong/ccnn;
        
        %%%%% Small World
        wC = G*squeeze(C2);
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        
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
        % actual modeling
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
            Xanalytic = hilbert(demean(ts(seed,:)));
            Phases2(seed,:) = angle(Xanalytic);
        end
        
        ccnn=0;
        cclong=0;
        for i=1:NPARCELLS
            for j=1:NPARCELLS
                if rr(i,j)>0.8*maxrr
                    cclong=cclong+(corr2(ts(i,:),ts(j,:)));
                    ccnn=ccnn+1;
                end
            end
        end
        FClarge2(sub)=cclong/ccnn;
        
        %%%%%
        ilam=1;
        for lambda=LAMBDA
            for i=1:NPARCELLS
                enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(squeeze(C1(ilam,i,:)));
                lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
                enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases2),sin(Phases2)))/sum(squeeze(C1(ilam,i,:)));
                lam_mean_spatime_enstrophy2(ilam,i,:)=abs(enstrophy);
            end
            if ilam>1
                [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))');
                mutinf1_1(ilam)=nanmean(cc(find(pp(:)<0.05)));
                [cc pp]=corr((squeeze(lam_mean_spatime_enstrophy2(ilam,:,2:end)))',(squeeze(lam_mean_spatime_enstrophy2(ilam-1,:,1:end-1)))');
                mutinf12_1(ilam)=nanmean(cc(find(pp(:)<0.05)));
            end
            ilam=ilam+1;
        end
        
        Rspatime=squeeze(lam_mean_spatime_enstrophy(Nscale,:,:));
        Rspatime2=squeeze(lam_mean_spatime_enstrophy2(Nscale,:,:));
        Turbulence(ii,sub)=nanstd(Rspatime(:));
        Turbulence2(ii,sub)=nanstd(Rspatime2(:));
        Transfer(ii,sub)=nanmean(mutinf1_1(2:end));
        Transfer2(ii,sub)=nanmean(mutinf12_1(2:end));
        TransferInfo(ii,sub,:)=mutinf1_1(2:end);
        TransferInfo2(ii,sub,:)=mutinf12_1(2:end);
    end
    ii=ii+1;
end

save('resultsring.mat','Turbulence','Turbulence2','Transfer','Transfer2','TransferInfo','TransferInfo2','FClarge','FClarge2');

%% Example Figures

figure(1)
%% short-range
plot(G_range,mean(Turbulence,2),'k');
hold on;
%% short-range + long-range exceptions (small world like)
plot(G_range,mean(Turbulence2,2),'r');



