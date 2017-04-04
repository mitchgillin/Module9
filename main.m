

%{
    ColumnAbsorb: OD Readers
(Elute#)  (1) (2) (X) (4) (5)(5E) (6) (6) (7)  (7)  (X) (X)
           1   2   3   4   5   6   7   8   9   10   11  12
        A 1Fold
        B 2
        C 4
        D 8
        E 16
        F 32
        G 64
        H 128


%}

dil=[1,2,4,8,16,32,64,128];
range=[1:5];
iRange=[1:10];
clear ssres sstot P s 
close all;

%i is how many Ellutions
%Range is how Many dilutoins we do

for i=iRange
    figure;
    [P,s]=polyfit(dil(range)',ColumnAbsorb(range',i),1);
    bfl=P(1)*dil(range)+P(2);
    disp(P);
    sstot = sum((ColumnAbsorb(range',i) - mean(ColumnAbsorb(:,i))).^2);
    
    ssres = sum((ColumnAbsorb(range',i) - bfl').^2);
    
    rsquared = 1 - (ssres / sstot);
    r2=sprintf('R^2 is: %f for %i',rsquared,i);
    
    subplot(2,1,1)
    hold on;
    plot(dil(range),bfl);
    xlabel('Dilution (Fold)');
    ylabel('Absorbance');
    scatter(dil(range),ColumnAbsorb(range',i));
    
    text(max(dil(range))/2,mean([max(ColumnAbsorb(range',i)),min(ColumnAbsorb(range',i))]),['R^2 = ',num2str(rsquared)]);
    str = sprintf('Fit of Elute: %f ',i);
    title(str);
    
    subplot(2,1,2)
    bar(ColumnAbsorb(i,:))
    str2=sprintf('Bar graph of Dillution %f ',dil(i));
    title(str2);
    xlabel('Elution Number');
end
%{
figure;
plot(ColumnAbsorb(5,1:12));
title('Column Elution OD Read: 128 dilution');
%}
