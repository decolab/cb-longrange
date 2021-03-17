clear all;

G_range=0.:0.01:3;
LAMBDA=[0.26 0.22 0.18 0.14 0.10 0.06 0.02];
DIST=1./LAMBDA;
NFUTURE=10;

for s=1:301
    fileName = sprintf('WG_%03d.mat',s);
    load(fileName);
    Iwe=s;
    err_hete_all(Iwe,:)=err_hete;
    err_hetetot_all(Iwe,:)=err_hetetot;
    err_hetetotrnd_all(Iwe,:)=err_hetetotrnd;
    err_heteSC_all(Iwe,:)=err_heteSC;
    fcfittlong_all(Iwe,:)=fcfittlong;
    fcfittlongtot_all(Iwe,:)=fcfittlongtot;
    fcfittlongtotrnd_all(Iwe,:)=fcfittlongtotrnd;
    fcfittlongSC_all(Iwe,:)=fcfittlongSC;
    fclong_all(Iwe,:)=abs(fclong);
    fclongtot_all(Iwe,:)=abs(fclongtot);
    fclongtotrnd_all(Iwe,:)=abs(fclongtotrnd);
    fclongSC_all(Iwe,:)=abs(fclongSC);
    Rmeta_all(Iwe,:)=Rmeta;
    Rmetatot_all(Iwe,:)=Rmetatot;
    Rmetatotrnd_all(Iwe,:)=Rmetatotrnd;
    RmetaSC_all(Iwe,:)=RmetaSC;
    infocapacity_all(Iwe)=infocapacity;
    infocapacitytot_all(Iwe)=infocapacitytot;
    infocapacitytotrnd_all(Iwe)=infocapacitytotrnd;
    infocapacitySC_all(Iwe)=infocapacitySC;
    susceptibility_all(Iwe)=susceptibility;
    susceptibilitytot_all(Iwe)=susceptibilitytot;
    susceptibilitytotrnd_all(Iwe)=susceptibilitytotrnd;
    susceptibilitySC_all(Iwe)=susceptibilitySC;
    errinfocapacity_all(Iwe)=errinfocapacity;
    errinfocapacitytot_all(Iwe)=errinfocapacitytot;
    errinfocapacitytotrnd_all(Iwe)=errinfocapacitytotrnd;
    errinfocapacitySC_all(Iwe)=errinfocapacitySC;
    errsusceptibility_all(Iwe)=errsusceptibility;
    errsusceptibilitytot_all(Iwe)=errsusceptibilitytot;
    errsusceptibilitytotrnd_all(Iwe)=errsusceptibilitytotrnd;
    errsusceptibilitySC_all(Iwe)=errsusceptibilitySC;
    Inflam_all(Iwe,:)=Inflam;
    Err_Inflam_all(Iwe,:)=Err_Inflam;
    Err_Rlam_all(Iwe,:)=Err_Rlam;
    Inflamtot_all(Iwe,:)=Inflamtot;
    Err_Inflamtot_all(Iwe,:)=Err_Inflamtot;
    Err_Rlamtot_all(Iwe,:)=Err_Rlamtot;
    Inflamtotrnd_all(Iwe,:)=Inflamtotrnd;
    Err_Inflamtotrnd_all(Iwe,:)=Err_Inflamtotrnd;
    Err_Rlamtotrnd_all(Iwe,:)=Err_Rlamtotrnd;
    InflamSC_all(Iwe,:)=InflamSC;
    Err_InflamSC_all(Iwe,:)=Err_InflamSC;
    Err_RlamSC_all(Iwe,:)=Err_RlamSC;
    InfoCascade_all(Iwe,:,:)=InfoCascade;
    InfoCascadetot_all(Iwe,:,:)=InfoCascadetot;
    InfoCascadetotrnd_all(Iwe,:,:)=InfoCascadetotrnd;
    InfoCascadeSC_all(Iwe,:,:)=InfoCascadeSC;
    CorrRSN_all(Iwe,:,:)=CorrRSN;
    CorrRSNtot_all(Iwe,:,:)=CorrRSNtot;
    CorrRSNtotrnd_all(Iwe,:,:)=CorrRSNtotrnd;
    CorrRSNSC_all(Iwe,:,:)=CorrRSNSC;
    mutinfo_p_all(Iwe,:,:)=mutinfo_p;
    mutinfo_ptot_all(Iwe,:,:)=mutinfo_ptot;
    mutinfo_ptotrnd_all(Iwe,:,:)=mutinfo_ptotrnd;
    mutinfo_pSC_all(Iwe,:,:)=mutinfo_pSC;
end

figure(1)
shadedErrorBar(G_range,nanmean(err_hete_all,2),nanstd(err_hete_all,[],2),'-k',0.7)
hold on;
shadedErrorBar(G_range,nanmean(err_hetetot_all,2),nanstd(err_hetetot_all,[],2),'-r',0.7)
shadedErrorBar(G_range,nanmean(err_heteSC_all,2),nanstd(err_heteSC_all,[],2),'-b',0.7)
shadedErrorBar(G_range,nanmean(err_hetetotrnd_all,2),nanstd(err_hetetotrnd_all,[],2),'-g',0.7)


figure(1)
shadedErrorBar(G_range,nanmean(err_hete_all,2),nanstd(err_hete_all,[],2),'-k',0.7)
hold on;
shadedErrorBar(G_range,nanmean(err_hetetot_all,2),nanstd(err_hetetot_all,[],2),'-r',0.7)
shadedErrorBar(G_range,nanmean(Rmeta_all,2),nanstd(Rmeta_all,[],2),'-g',0.7)
shadedErrorBar(G_range,nanmean(Rmetatot_all,2),nanstd(Rmetatot_all,[],2),'-b',0.7)


fitting=nanmean(err_hete_all,2);
[aux index]=min(fitting);
fitting=nanmean(err_hetetot_all,2);
[aux indextot]=min(fitting);
fitting=nanmean(err_heteSC_all,2);
[aux indexSC]=min(fitting);
fitting=nanmean(err_hetetotrnd_all,2);
[aux indextotrnd]=min(fitting);

% index=indextot;
% indexSC=indextot;

figure(1)
boxplot([err_hetetot_all(indextot,:)' err_hetetotrnd_all(indextotrnd,:)' err_heteSC_all(indexSC,:)' err_hete_all(index,:)']);

ranksum(err_hetetot_all(indextot,:)',err_hetetotrnd_all(indextotrnd,:)')
ranksum(err_hetetot_all(indextot,:)',err_heteSC_all(indexSC,:)')
ranksum(err_hetetot_all(indextot,:)',err_hete_all(index,:)')

figure(2)
boxplot([Err_Rlamtot_all(indextot,:)' Err_Rlamtotrnd_all(indextotrnd,:)' Err_RlamSC_all(indexSC,:)' Err_Rlam_all(index,:)']);

figure(3)
boxplot([Err_Inflamtot_all(indextot,:)' Err_Inflamtotrnd_all(indextotrnd,:)' Err_InflamSC_all(indexSC,:)' Err_Inflam_all(index,:)']);

figure(4)
boxplot([fcfittlongtot_all(indextot,:)' fcfittlongtotrnd_all(indextotrnd,:)' fcfittlongSC_all(indexSC,:)' fcfittlong_all(index,:)']);

ranksum(fcfittlongtot_all(indextot,:)',fcfittlongtotrnd_all(indextotrnd,:)')
ranksum(fcfittlongtot_all(indextot,:)',fcfittlongSC_all(indexSC,:)')
ranksum(fcfittlongtot_all(indextot,:)',fcfittlong_all(index,:)')

figure(5)
boxplot([Rmetatot_all(indextot,:)' Rmetatotrnd_all(indextotrnd,:)' RmetaSC_all(indexSC,:)' Rmeta_all(index,:)']);

figure(6)
boxplot([Inflamtot_all(indextot,:)' Inflamtotrnd_all(indextotrnd,:)' InflamSC_all(indexSC,:)' Inflam_all(index,:)']);

figure(7)
bar([1 2 3 4],[infocapacitytot_all(indextot) infocapacitytotrnd_all(indextotrnd) infocapacitySC_all(indexSC) infocapacity_all(index)]);
hold on;
errorbar([1 2 3 4],[infocapacitytot_all(indextot) infocapacitytotrnd_all(indextotrnd) infocapacitySC_all(indexSC) infocapacity_all(index)],[errinfocapacitytot_all(indextot) errinfocapacitytotrnd_all(indextot) errinfocapacitySC_all(indexSC) errinfocapacity_all(index)],'linestyle','none');

figure(8)
bar([1 2 3 4],[susceptibilitytot_all(indextot) susceptibilitytotrnd_all(indextotrnd) susceptibilitySC_all(indexSC) susceptibility_all(index)]);
hold on;
errorbar([1 2 3 4],[susceptibilitytot_all(indextot) susceptibilitytotrnd_all(indextotrnd) susceptibilitySC_all(indexSC) susceptibility_all(index)],[errsusceptibilitytot_all(indextot) errsusceptibilitytotrnd_all(indextot) errsusceptibilitySC_all(indexSC) errsusceptibility_all(index)],'linestyle','none');

figure(9)
subplot(7,1,1)
boxplot([squeeze(CorrRSNtot_all(indextot,1,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,1,:)) squeeze(CorrRSNSC_all(indexSC,1,:)) squeeze(CorrRSN_all(index,1,:))]);
subplot(7,1,2)
boxplot([squeeze(CorrRSNtot_all(indextot,2,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,2,:)) squeeze(CorrRSNSC_all(indexSC,2,:)) squeeze(CorrRSN_all(index,2,:))]);
subplot(7,1,3)
boxplot([squeeze(CorrRSNtot_all(indextot,3,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,3,:)) squeeze(CorrRSNSC_all(indexSC,3,:)) squeeze(CorrRSN_all(index,3,:))]);
subplot(7,1,4)
boxplot([squeeze(CorrRSNtot_all(indextot,4,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,4,:)) squeeze(CorrRSNSC_all(indexSC,4,:)) squeeze(CorrRSN_all(index,4,:))]);
subplot(7,1,5)
boxplot([squeeze(CorrRSNtot_all(indextot,5,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,5,:)) squeeze(CorrRSNSC_all(indexSC,5,:)) squeeze(CorrRSN_all(index,5,:))]);
subplot(7,1,6)
boxplot([squeeze(CorrRSNtot_all(indextot,6,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,6,:)) squeeze(CorrRSNSC_all(indexSC,6,:)) squeeze(CorrRSN_all(index,6,:))]);
subplot(7,1,7)
boxplot([squeeze(CorrRSNtot_all(indextot,7,:)) squeeze(CorrRSNtotrnd_all(indextotrnd,7,:)) squeeze(CorrRSNSC_all(indexSC,7,:)) squeeze(CorrRSN_all(index,7,:))]);
%% 

figure(10)
NIwe=size(mutinfo_p_all,1);
for i=1:NIwe
    mip(:,i)=nanmean(squeeze(mutinfo_p_all(i,:,:)))';
    miptot(:,i)=nanmean(squeeze(mutinfo_ptot_all(i,:,:)))';
    miptotrnd(:,i)=nanmean(squeeze(mutinfo_ptotrnd_all(i,:,:)))';
    mipSC(:,i)=nanmean(squeeze(mutinfo_pSC_all(i,:,:)))';
end
plot(mip(:,index),'-k');
hold on;
plot(miptot(:,indextot),'-r');
plot(miptotrnd(:,indextotrnd),'-g');
plot(mipSC(:,indexSC),'-b');


figure(11)
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascade_all(index,:,:)),1),nanstd(squeeze(InfoCascade_all(index,:,:)),[],1),'-k',0.7)
hold on;
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascadetot_all(indextot,:,:)),1),nanstd(squeeze(InfoCascadetot_all(indextot,:,:)),[],1),'-r',0.7)
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascadeSC_all(indexSC,:,:)),1),nanstd(squeeze(InfoCascadeSC_all(indexSC,:,:)),[],1),'-b',0.7)
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascadetotrnd_all(indextotrnd,:,:)),1),nanstd(squeeze(InfoCascadetotrnd_all(indextotrnd,:,:)),[],1),'-g',0.7)

save results_hopf_SClong_3_40.mat G_range LAMBDA DIST InfoCascade_all InfoCascadetot_all InfoCascadetotrnd_all InfoCascadeSC_all mip miptot miptotrnd mipSC mutinfo_p_all mutinfo_ptot_all mutinfo_ptotrnd_all mutinfo_pSC_all CorrRSN_all CorrRSNtot_all CorrRSNtotrnd_all CorrRSNSC_all Inflam_all Err_Inflam_all Err_Rlam_all Inflamtot_all Err_Inflamtot_all Err_Rlamtot_all Inflamtotrnd_all Err_Inflamtotrnd_all Err_Rlamtotrnd_all InflamSC_all Err_InflamSC_all Err_RlamSC_all infocapacity_all infocapacitytot_all infocapacitytotrnd_all infocapacitySC_all susceptibility_all susceptibilitytot_all susceptibilitytotrnd_all susceptibilitySC_all errinfocapacity_all errinfocapacitytot_all errinfocapacitytotrnd_all errinfocapacitySC_all errsusceptibility_all errsusceptibilitytot_all errsusceptibilitytotrnd_all errsusceptibilitySC_all Rmeta_all Rmetatot_all Rmetatotrnd_all RmetaSC_all fcfittlong_all fcfittlongtot_all fcfittlongtotrnd_all fcfittlongSC_all err_hete_all err_hetetot_all err_hetetotrnd_all err_heteSC_all;



figure(1)
boxplot([err_hetetot_all(indextot,:)' err_hete_all(index,:)']);
ranksum(err_hetetot_all(indextot,:)',err_hete_all(index,:)')

figure(2)
boxplot([fcfittlongtot_all(indextot,:)' fcfittlong_all(index,:)']);
ranksum(fcfittlongtot_all(indextot,:)',fcfittlong_all(index,:)')

figure(3)
boxplot([Inflamtot_all(indextot,:)' Inflam_all(index,:)']);
ranksum(Inflamtot_all(indextot,:)',Inflam_all(index,:)')

figure(4)
bar([1 2],[infocapacitytot_all(indextot) infocapacity_all(index)]);
hold on;
errorbar([1 2],[infocapacitytot_all(indextot) infocapacity_all(index)],[errinfocapacitytot_all(indextot) errinfocapacity_all(index)],'linestyle','none');

figure(5)
bar([1 2],[susceptibilitytot_all(indextot) susceptibility_all(index)]);
hold on;
errorbar([1 2],[susceptibilitytot_all(indextot) susceptibility_all(index)],[errsusceptibilitytot_all(indextot) errsusceptibility_all(index)],'linestyle','none');

figure(6)
subplot(7,1,1)
boxplot([squeeze(CorrRSNtot_all(indextot,1,:)) squeeze(CorrRSN_all(index,1,:))]);
subplot(7,1,2)
boxplot([squeeze(CorrRSNtot_all(indextot,2,:))  squeeze(CorrRSN_all(index,2,:))]);
subplot(7,1,3)
boxplot([squeeze(CorrRSNtot_all(indextot,3,:)) squeeze(CorrRSN_all(index,3,:))]);
subplot(7,1,4)
boxplot([squeeze(CorrRSNtot_all(indextot,4,:))  squeeze(CorrRSN_all(index,4,:))]);
subplot(7,1,5)
boxplot([squeeze(CorrRSNtot_all(indextot,5,:)) squeeze(CorrRSN_all(index,5,:))]);
subplot(7,1,6)
boxplot([squeeze(CorrRSNtot_all(indextot,6,:)) squeeze(CorrRSN_all(index,6,:))]);
subplot(7,1,7)
boxplot([squeeze(CorrRSNtot_all(indextot,7,:)) squeeze(CorrRSN_all(index,7,:))]);
%% 
ranksum(squeeze(CorrRSNtot_all(indextot,1,:)),squeeze(CorrRSN_all(index,1,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,2,:)),squeeze(CorrRSN_all(index,2,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,3,:)),squeeze(CorrRSN_all(index,3,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,4,:)),squeeze(CorrRSN_all(index,4,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,5,:)),squeeze(CorrRSN_all(index,5,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,6,:)),squeeze(CorrRSN_all(index,6,:)))
ranksum(squeeze(CorrRSNtot_all(indextot,7,:)),squeeze(CorrRSN_all(index,7,:)))


figure(7)
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascade_all(index,:,:)),1)-nanmean(squeeze(InfoCascade_all(1,:,:)),1),nanstd(squeeze(InfoCascade_all(index,:,:)),[],1),'-k',0.7)
hold on;
shadedErrorBar(LAMBDA(2:end),nanmean(squeeze(InfoCascadetot_all(indextot,:,:)),1)-nanmean(squeeze(InfoCascadetot_all(1,:,:)),1),nanstd(squeeze(InfoCascadetot_all(indextot,:,:)),[],1),'-r',0.7)

figure(8)
plot(mip(:,index),'-k');
hold on;
plot(miptot(:,indextot),'-r');

figure(9)
boxplot([fclongtot_all(indextot,:)' fclong_all(index,:)']);
ranksum(fclongtot_all(indextot,:)',fclong_all(index,:)')

figure(10)
boxplot([Rmetatot_all(indextot,:)' Rmeta_all(index,:)']);
