function [Results, ResultsHeader] = stats_TtestAllStats(Sample1,Sample2)
%TTESTALLSTATS calculates all the possible statistics associated with a
    %two-sample Ttest.

    %invidual sample statstics:
    s1=length(Sample1);
    s2=length(Sample2);

    mean1=mean(Sample1);
    mean2=mean(Sample2);

    std1=std(Sample1);
    std2=std(Sample2);

    %comparison statistics:
    CohenD = stats_computeCohen_d(Sample1,Sample2);

    [~,pVal,ci95,stats]=ttest2(Sample1,Sample2); %two sample (unpaired T-test)
    CI95_l=ci95(1);
    CI95_u=ci95(2);
    Tstat=stats.tstat;

    %put into variable
    Results=[Tstat,pVal,CohenD,CI95_l,CI95_u,mean1,mean2,std1,std2, s1,s2];
    ResultsHeader={'Tstat','pVal','CohenD','CI95_l','CI95_u','mean1','mean2','sd1', 'sd2','N1','N2'};
end

