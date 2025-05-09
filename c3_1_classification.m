close all;clc;
%%
% Inputs: 
% --------
% MAVClass1: the features of the VF case (stimulus and rest features)
% MAVClass2: the features of the Pinch case (stimulus and rest features)
% TriggerClass1: labels for VF features (stimulus or rest label)
% TriggerClass2: labels for Pinch features (stimulus or rest label)

MAVClass1 = VFMAVfeatures;
VARClass1 = VFVARfeatures;
MAVClass2 = PinchMAVfeatures;
VARClass2 = PinchVARfeatures;
MAVClass3 = FlexMAVfeatures;
VARClass3 = FlexVARfeatures;

TriggerClass1 = VFlabels;
TriggerClass2 = Pinchlabels;
TriggerClass3 = Flexlabels;

% Build the datasets
MAV_class1 = MAVClass1(find(TriggerClass1==1));
MAV_rest1 = MAVClass1(find(TriggerClass1==0));

VAR_class1 = VARClass1(find(TriggerClass1==1));
VAR_rest1 = VARClass1(find(TriggerClass1==0));

MAV_class2 = MAVClass2(find(TriggerClass2==1));
MAV_rest2 = MAVClass2(find(TriggerClass2==0));

VAR_class2 = VARClass2(find(TriggerClass2==1));
VAR_rest2 = VARClass2(find(TriggerClass2==0));

MAV_class3 = MAVClass3(find(TriggerClass3==1));
MAV_rest3 = MAVClass3(find(TriggerClass3==0));

VAR_class3 = VARClass3(find(TriggerClass3==1));
VAR_rest3 = VARClass3(find(TriggerClass3==0));

% Concantenate the rest classes
MAV_rest = [MAV_rest1 MAV_rest2 MAV_rest3];
VAR_rest = [VAR_rest1 VAR_rest2 VAR_rest3];


%%
% Class1 vs Rest dataset
MAV_Data_Class1vsRest = [MAV_class1 MAV_rest];
MAV_Labels_Class1vsRest = [ones(1,length(MAV_class1)) 2*ones(1,length(MAV_rest))];

VAR_Data_Class1vsRest = [VAR_class1 VAR_rest];
VAR_Labels_Class1vsRest = MAV_Labels_Class1vsRest;

% Class2 vs Rest dataset
MAV_Data_Class2vsRest = [MAV_class2 MAV_rest];
MAV_Labels_Class2vsRest = [ones(1,length(MAV_class2)) 2*ones(1,length(MAV_rest))];

VAR_Data_Class2vsRest = [VAR_class2 VAR_rest];
VAR_Labels_Class2vsRest = MAV_Labels_Class2vsRest;

% Class1 vs Class2 dataset
MAV_Data_Class1vsClass2 = [MAV_class1 MAV_class2];
MAV_Labels_Class1vsClass2 = [ones(1,length(MAV_class1)) 2*ones(1,length(MAV_class2))];

VAR_Data_Class1vsClass2 = [VAR_class1 VAR_class2];
VAR_Labels_Class1vsClass2 = MAV_Labels_Class1vsClass2;



% Class3 vs Rest dataset
MAV_Data_Class3vsRest = [MAV_class3 MAV_rest];
MAV_Labels_Class3vsRest = [ones(1,length(MAV_class3)) 2*ones(1,length(MAV_rest))];

VAR_Data_Class3vsRest = [VAR_class3 VAR_rest];
VAR_Labels_Class3vsRest = MAV_Labels_Class3vsRest;

% Class1 vs Class3 dataset
MAV_Data_Class1vsClass3 = [MAV_class1 MAV_class3];
MAV_Labels_Class1vsClass3 = [ones(1,length(MAV_class1)) 2*ones(1,length(MAV_class3))];

VAR_Data_Class1vsClass3 = [VAR_class1 VAR_class3];
VAR_Labels_Class1vsClass3 = MAV_Labels_Class1vsClass3;

% Class2 vs Class3 dataset
MAV_Data_Class2vsClass3 = [MAV_class2 MAV_class3];
MAV_Labels_Class2vsClass3 = [ones(1,length(MAV_class2)) 2*ones(1,length(MAV_class3))];

VAR_Data_Class2vsClass3 = [VAR_class2 VAR_class3];
VAR_Labels_Class2vsClass3 = MAV_Labels_Class2vsClass3;



%%
% Both feature datasets
MAVVAR_Data_Class1vsRest = [MAV_Data_Class1vsRest; VAR_Data_Class1vsRest];
MAVVAR_Labels_Class1vsRest = MAV_Labels_Class1vsRest;

MAVVAR_Data_Class2vsRest = [MAV_Data_Class2vsRest; VAR_Data_Class2vsRest];
MAVVAR_Labels_Class2vsRest = MAV_Labels_Class2vsRest;

MAVVAR_Data_Class1vsClass2 = [MAV_Data_Class1vsClass2; VAR_Data_Class1vsClass2];
MAVVAR_Labels_Class1vsClass2 = MAV_Labels_Class1vsClass2;

MAVVAR_Data_Class3vsRest = [MAV_Data_Class3vsRest; VAR_Data_Class3vsRest];
MAVVAR_Labels_Class3vsRest = MAV_Labels_Class3vsRest;

MAVVAR_Data_Class1vsClass3 = [MAV_Data_Class1vsClass3; VAR_Data_Class1vsClass3];
MAVVAR_Labels_Class1vsClass3 = MAV_Labels_Class1vsClass3;

MAVVAR_Data_Class2vsClass3 = [MAV_Data_Class2vsClass3; VAR_Data_Class2vsClass3];
MAVVAR_Labels_Class2vsClass3 = MAV_Labels_Class2vsClass3;

%%
% Classify all combinations (training set)
k = 10; % for k-fold cross validation
c1 = cvpartition(length(MAV_Labels_Class1vsRest),'KFold',k);
c2 = cvpartition(length(VAR_Labels_Class1vsRest),'KFold',k);
c3 = cvpartition(length(MAVVAR_Labels_Class1vsRest),'KFold',k);

c4 = cvpartition(length(MAV_Labels_Class2vsRest),'KFold',k);
c5 = cvpartition(length(VAR_Labels_Class2vsRest),'KFold',k);
c6 = cvpartition(length(MAVVAR_Labels_Class2vsRest),'KFold',k);

c7 = cvpartition(length(MAV_Labels_Class1vsClass2),'KFold',k);
c8 = cvpartition(length(VAR_Labels_Class1vsClass2),'KFold',k);
c9 = cvpartition(length(MAVVAR_Labels_Class1vsClass2),'KFold',k);


c10 = cvpartition(length(MAV_Labels_Class3vsRest),'KFold',k);
c11 = cvpartition(length(VAR_Labels_Class3vsRest),'KFold',k);
c12 = cvpartition(length(MAVVAR_Labels_Class3vsRest),'KFold',k);

c13 = cvpartition(length(MAV_Labels_Class1vsClass3),'KFold',k);
c14 = cvpartition(length(VAR_Labels_Class1vsClass3),'KFold',k);
c15 = cvpartition(length(MAVVAR_Labels_Class1vsClass3),'KFold',k);

c16 = cvpartition(length(MAV_Labels_Class2vsClass3),'KFold',k);
c17 = cvpartition(length(VAR_Labels_Class2vsClass3),'KFold',k);
c18 = cvpartition(length(MAVVAR_Labels_Class2vsClass3),'KFold',k);

% Intialize total variables for averaging
totalTstAcc_MAV_C1rest = 0;
totalTstAcc_VAR_C1rest = 0;
totalTstAcc_MAVVAR_C1rest = 0;
totalTstAcc_MAV_C2rest = 0;
totalTstAcc_VAR_C2rest = 0;
totalTstAcc_MAVVAR_C2rest = 0;
totalTstAcc_MAV_C1C2 = 0;
totalTstAcc_VAR_C1C2 = 0;
totalTstAcc_MAVVAR_C1C2 = 0;
totalTstAcc_MAV_C3rest = 0;
totalTstAcc_VAR_C3rest = 0;
totalTstAcc_MAVVAR_C3rest = 0;
totalTstAcc_MAV_C1C3 = 0;
totalTstAcc_VAR_C1C3 = 0;
totalTstAcc_MAVVAR_C1C3 = 0;
totalTstAcc_MAV_C2C3 = 0;
totalTstAcc_VAR_C2C3 = 0;
totalTstAcc_MAVVAR_C2C3 = 0;
totalTstCM_MAV_C1rest = 0;
totalTstCM_VAR_C1rest = 0;
totalTstCM_MAVVAR_C1rest = 0;
totalTstCM_MAV_C2rest = 0;
totalTstCM_VAR_C2rest = 0;
totalTstCM_MAVVAR_C2rest = 0;
totalTstCM_MAV_C1C2 = 0;
totalTstCM_VAR_C1C2 = 0;
totalTstCM_MAVVAR_C1C2 = 0;
totalTstCM_MAV_C3rest = 0;
totalTstCM_VAR_C3rest = 0;
totalTstCM_MAVVAR_C3rest = 0;
totalTstCM_MAV_C1C3 = 0;
totalTstCM_VAR_C1C3 = 0;
totalTstCM_MAVVAR_C1C3 = 0;
totalTstCM_MAV_C2C3 = 0;
totalTstCM_VAR_C2C3 = 0;
totalTstCM_MAVVAR_C2C3 = 0;


% Repeat the following for i=1:k, and average performance metrics across all iterations
i=1;
% loop over all k-folds and avergae the performance
for i=1:k

    [TstMAVFC1Rest TstMAVErrC1Rest] = classify(MAV_Data_Class1vsRest(c1.test(i))',MAV_Data_Class1vsRest(c1.training(i))',MAV_Labels_Class1vsRest(c1.training(i)));
    [TstCM_MAV_C1rest dum1 TstAcc_MAV_C1rest dum2] = confusion(MAV_Labels_Class1vsRest(c1.test(i)), TstMAVFC1Rest);
    totalTstAcc_MAV_C1rest = totalTstAcc_MAV_C1rest + TstAcc_MAV_C1rest;
    totalTstCM_MAV_C1rest = totalTstCM_MAV_C1rest + TstCM_MAV_C1rest;

    [TstVARFC1Rest TstVARErrC1Rest] = classify(VAR_Data_Class1vsRest(c2.test(i))',VAR_Data_Class1vsRest(c2.training(i))',VAR_Labels_Class1vsRest(c2.training(i)));
    [TstCM_VAR_C1rest dum1 TstAcc_VAR_C1rest dum2] = confusion(VAR_Labels_Class1vsRest(c2.test(i)), TstVARFC1Rest);
    totalTstAcc_VAR_C1rest = totalTstAcc_VAR_C1rest + TstAcc_VAR_C1rest;
    totalTstCM_VAR_C1rest = totalTstCM_VAR_C1rest + TstCM_VAR_C1rest;


    [TstMAVVARFC1Rest TstMAVVARErrC1Rest] = classify(MAVVAR_Data_Class1vsRest(:,c3.test(i))',MAVVAR_Data_Class1vsRest(:,c3.training(i))',MAVVAR_Labels_Class1vsRest(c3.training(i)));
    [TstCM_MAVVAR_C1rest dum1 TstAcc_MAVVAR_C1rest dum2] = confusion(MAVVAR_Labels_Class1vsRest(c3.test(i)), TstMAVVARFC1Rest);
    totalTstAcc_MAVVAR_C1rest = totalTstAcc_MAVVAR_C1rest + TstAcc_MAVVAR_C1rest;
    totalTstCM_MAVVAR_C1rest = totalTstCM_MAVVAR_C1rest + TstCM_MAVVAR_C1rest;

    % Class2 vs Rest
    [TstMAVFC2Rest TstMAVErrC2Rest] = classify(MAV_Data_Class2vsRest(c4.test(i))',MAV_Data_Class2vsRest(c4.training(i))',MAV_Labels_Class2vsRest(c4.training(i)));
    [TstCM_MAV_C2rest dum1 TstAcc_MAV_C2rest dum2] = confusion(MAV_Labels_Class2vsRest(c4.test(i)), TstMAVFC2Rest);
    totalTstAcc_MAV_C2rest = totalTstAcc_MAV_C2rest + TstAcc_MAV_C2rest;
    totalTstCM_MAV_C2rest = totalTstCM_MAV_C2rest + TstCM_MAV_C2rest;

    [TstVARFC2Rest TstVARErrC2Rest] = classify(VAR_Data_Class2vsRest(c5.test(i))',VAR_Data_Class2vsRest(c5.training(i))',VAR_Labels_Class2vsRest(c5.training(i)));
    [TstCM_VAR_C2rest dum1 TstAcc_VAR_C2rest dum2] = confusion(VAR_Labels_Class2vsRest(c5.test(i)), TstVARFC2Rest);
    totalTstAcc_VAR_C2rest = totalTstAcc_VAR_C2rest + TstAcc_VAR_C2rest;
    totalTstCM_VAR_C2rest = totalTstCM_VAR_C2rest + TstCM_VAR_C2rest;

    [TstMAVVARFC2Rest TstMAVVARErrC2Rest] = classify(MAVVAR_Data_Class2vsRest(:,c6.test(i))',MAVVAR_Data_Class2vsRest(:,c6.training(i))',MAVVAR_Labels_Class2vsRest(c6.training(i)));
    [TstCM_MAVVAR_C2rest dum1 TstAcc_MAVVAR_C2rest dum2] = confusion(MAVVAR_Labels_Class2vsRest(c6.test(i)), TstMAVVARFC2Rest);
    totalTstAcc_MAVVAR_C2rest = totalTstAcc_MAVVAR_C2rest + TstAcc_MAVVAR_C2rest;
    totalTstCM_MAVVAR_C2rest = totalTstCM_MAVVAR_C2rest + TstCM_MAVVAR_C2rest;

    % Class1 vs Class2
    [TstMAVFC1C2 TstMAVErrC1C2] = classify(MAV_Data_Class1vsClass2(c7.test(i))',MAV_Data_Class1vsClass2(c7.training(i))',MAV_Labels_Class1vsClass2(c7.training(i)));
    [TstCM_MAV_C1C2 dum1 TstAcc_MAV_C1C2 dum2] = confusion(MAV_Labels_Class1vsClass2(c7.test(i)), TstMAVFC1C2);
    totalTstAcc_MAV_C1C2 = totalTstAcc_MAV_C1C2 + TstAcc_MAV_C1C2;
    totalTstCM_MAV_C1C2 = totalTstCM_MAV_C1C2 + TstCM_MAV_C1C2;

    [TstVARFC1C2 TstVARErrC1C2] = classify(VAR_Data_Class1vsClass2(c8.test(i))',VAR_Data_Class1vsClass2(c8.training(i))',VAR_Labels_Class1vsClass2(c8.training(i)));
    [TstCM_VAR_C1C2 dum1 TstAcc_VAR_C1C2 dum2] = confusion(VAR_Labels_Class1vsClass2(c8.test(i)), TstVARFC1C2);
    totalTstAcc_VAR_C1C2 = totalTstAcc_VAR_C1C2 + TstAcc_VAR_C1C2;
    totalTstCM_VAR_C1C2 = totalTstCM_VAR_C1C2 + TstCM_VAR_C1C2;

    [TstMAVVARFC1C2 TstMAVVARErrC1C2] = classify(MAVVAR_Data_Class1vsClass2(:,c9.test(i))',MAVVAR_Data_Class1vsClass2(:,c9.training(i))',MAVVAR_Labels_Class1vsClass2(c9.training(i)));
    [TstCM_MAVVAR_C1C2 dum1 TstAcc_MAVVAR_C1C2 dum2] = confusion(MAVVAR_Labels_Class1vsClass2(c9.test(i)), TstMAVVARFC1C2);
    totalTstAcc_MAVVAR_C1C2 = totalTstAcc_MAVVAR_C1C2 + TstAcc_MAVVAR_C1C2;
    totalTstCM_MAVVAR_C1C2 = totalTstCM_MAVVAR_C1C2 + TstCM_MAVVAR_C1C2;

    %Class 3 vs rest
    [TstMAVFC3Rest TstMAVErrC3Rest] = classify(MAV_Data_Class3vsRest(c10.test(i))',MAV_Data_Class3vsRest(c10.training(i))',MAV_Labels_Class3vsRest(c10.training(i)));
    [TstCM_MAV_C3rest dum1 TstAcc_MAV_C3rest dum2] = confusion(MAV_Labels_Class3vsRest(c10.test(i)), TstMAVFC3Rest);
    totalTstAcc_MAV_C3rest = totalTstAcc_MAV_C3rest + TstAcc_MAV_C3rest;
    totalTstCM_MAV_C3rest = totalTstCM_MAV_C3rest + TstCM_MAV_C3rest;

    [TstVARFC3Rest TstVARErrC3Rest] = classify(VAR_Data_Class3vsRest(c11.test(i))',VAR_Data_Class3vsRest(c11.training(i))',VAR_Labels_Class3vsRest(c11.training(i)));
    [TstCM_VAR_C3rest dum1 TstAcc_VAR_C3rest dum2] = confusion(VAR_Labels_Class3vsRest(c11.test(i)), TstVARFC3Rest);
    totalTstAcc_VAR_C3rest = totalTstAcc_VAR_C3rest + TstAcc_VAR_C3rest;
    totalTstCM_VAR_C3rest = totalTstCM_VAR_C3rest + TstCM_VAR_C3rest;

    [TstMAVVARFC3Rest TstMAVVARErrC3Rest] = classify(MAVVAR_Data_Class3vsRest(:,c12.test(i))',MAVVAR_Data_Class3vsRest(:,c12.training(i))',MAVVAR_Labels_Class3vsRest(c12.training(i)));
    [TstCM_MAVVAR_C3rest dum1 TstAcc_MAVVAR_C3rest dum2] = confusion(MAVVAR_Labels_Class3vsRest(c12.test(i)), TstMAVVARFC3Rest);
    totalTstAcc_MAVVAR_C3rest = totalTstAcc_MAVVAR_C3rest + TstAcc_MAVVAR_C3rest;
    totalTstCM_MAVVAR_C3rest = totalTstCM_MAVVAR_C3rest + TstCM_MAVVAR_C3rest;

    %Class 1 vs Class 3

    [TstMAVFC1C3 TstMAVErrC1C3] = classify(MAV_Data_Class1vsClass3(c13.test(i))',MAV_Data_Class1vsClass3(c13.training(i))',MAV_Labels_Class1vsClass3(c13.training(i)));
    [TstCM_MAV_C1C3 dum1 TstAcc_MAV_C1C3 dum2] = confusion(MAV_Labels_Class1vsClass3(c13.test(i)), TstMAVFC1C3);
    totalTstAcc_MAV_C1C3 = totalTstAcc_MAV_C1C3 + TstAcc_MAV_C1C3;
    totalTstCM_MAV_C1C3 = totalTstCM_MAV_C1C3 + TstCM_MAV_C1C3;

    [TstVARFC1C3 TstVARErrC1C3] = classify(VAR_Data_Class1vsClass3(c14.test(i))',VAR_Data_Class1vsClass3(c14.training(i))',VAR_Labels_Class1vsClass3(c14.training(i)));
    [TstCM_VAR_C1C3 dum1 TstAcc_VAR_C1C3 dum2] = confusion(VAR_Labels_Class1vsClass3(c14.test(i)), TstVARFC1C3);
    totalTstAcc_VAR_C1C3 = totalTstAcc_VAR_C1C3 + TstAcc_VAR_C1C3;
    totalTstCM_VAR_C1C3 = totalTstCM_VAR_C1C3 + TstCM_VAR_C1C3;

    [TstMAVVARFC1C3 TstMAVVARErrC1C3] = classify(MAVVAR_Data_Class1vsClass3(:,c15.test(i))',MAVVAR_Data_Class1vsClass3(:,c15.training(i))',MAVVAR_Labels_Class1vsClass3(c15.training(i)));
    [TstCM_MAVVAR_C1C3 dum1 TstAcc_MAVVAR_C1C3 dum2] = confusion(MAVVAR_Labels_Class1vsClass3(c15.test(i)), TstMAVVARFC1C3);
    totalTstAcc_MAVVAR_C1C3 = totalTstAcc_MAVVAR_C1C3 + TstAcc_MAVVAR_C1C3;
    totalTstCM_MAVVAR_C1C3 = totalTstCM_MAVVAR_C1C3 + TstCM_MAVVAR_C1C3;


    %Class 2 vs Class 3

    [TstMAVFC2C3 TstMAVErrC2C3] = classify(MAV_Data_Class2vsClass3(c16.test(i))',MAV_Data_Class2vsClass3(c16.training(i))',MAV_Labels_Class2vsClass3(c16.training(i)));
    [TstCM_MAV_C2C3 dum1 TstAcc_MAV_C2C3 dum2] = confusion(MAV_Labels_Class2vsClass3(c16.test(i)), TstMAVFC2C3);
    totalTstAcc_MAV_C2C3 = totalTstAcc_MAV_C2C3 + TstAcc_MAV_C2C3;
    totalTstCM_MAV_C2C3 = totalTstCM_MAV_C2C3 + TstCM_MAV_C2C3;

    [TstVARFC2C3 TstVARErrC2C3] = classify(VAR_Data_Class2vsClass3(c17.test(i))',VAR_Data_Class2vsClass3(c17.training(i))',VAR_Labels_Class2vsClass3(c17.training(i)));
    [TstCM_VAR_C2C3 dum1 TstAcc_VAR_C2C3 dum2] = confusion(VAR_Labels_Class2vsClass3(c17.test(i)), TstVARFC2C3);
    totalTstAcc_VAR_C2C3 = totalTstAcc_VAR_C2C3 + TstAcc_VAR_C2C3;
    totalTstCM_VAR_C2C3 = totalTstCM_VAR_C2C3 + TstCM_VAR_C2C3;

    [TstMAVVARFC2C3 TstMAVVARErrC2C3] = classify(MAVVAR_Data_Class2vsClass3(:,c18.test(i))',MAVVAR_Data_Class2vsClass3(:,c18.training(i))',MAVVAR_Labels_Class2vsClass3(c18.training(i)));
    [TstCM_MAVVAR_C2C3 dum1 TstAcc_MAVVAR_C2C3 dum2] = confusion(MAVVAR_Labels_Class2vsClass3(c18.test(i)), TstMAVVARFC2C3);
    totalTstAcc_MAVVAR_C2C3 = totalTstAcc_MAVVAR_C2C3 + TstAcc_MAVVAR_C2C3;
    totalTstCM_MAVVAR_C2C3 = totalTstCM_MAVVAR_C2C3 + TstCM_MAVVAR_C2C3;

end
%%


disp('Before division:')
disp(totalTstAcc_MAV_C1rest)
totalTstAcc_MAV_C1rest = totalTstAcc_MAV_C1rest/10;
disp('After division:')
disp(totalTstAcc_MAV_C1rest)


disp('Before division:')
disp(totalTstAcc_MAV_C1rest)
totalTstAcc_MAV_C1rest = totalTstAcc_MAV_C1rest/10;
disp('After division:')
disp(totalTstAcc_MAV_C1rest)

% %divide all the performance metrics by k to get the average
%totalTstAcc_MAV_C1rest = totalTstAcc_MAV_C1rest/10;
totalTstAcc_VAR_C1rest = totalTstAcc_VAR_C1rest/10;
totalTstAcc_MAVVAR_C1rest = totalTstAcc_MAVVAR_C1rest/10;
totalTstAcc_MAV_C2rest = totalTstAcc_MAV_C2rest/10;
totalTstAcc_VAR_C2rest = totalTstAcc_VAR_C2rest/10;
totalTstAcc_MAVVAR_C2rest = totalTstAcc_MAVVAR_C2rest/10;
totalTstAcc_MAV_C1C2 = totalTstAcc_MAV_C1C2/10;
totalTstAcc_VAR_C1C2 = totalTstAcc_VAR_C1C2/10;
totalTstAcc_MAVVAR_C1C2 = totalTstAcc_MAVVAR_C1C2/10;
totalTstAcc_MAV_C3rest = totalTstAcc_MAV_C3rest/10;
totalTstAcc_VAR_C3rest = totalTstAcc_VAR_C3rest/10;
totalTstAcc_MAVVAR_C3rest = totalTstAcc_MAVVAR_C3rest/10;
totalTstAcc_MAV_C1C3 = totalTstAcc_MAV_C1C3/10;
totalTstAcc_VAR_C1C3 = totalTstAcc_VAR_C1C3/10;
totalTstAcc_MAVVAR_C1C3 = totalTstAcc_MAVVAR_C1C3/10;
totalTstAcc_MAV_C2C3 = totalTstAcc_MAV_C2C3/10;
totalTstAcc_VAR_C2C3 = totalTstAcc_VAR_C2C3/10;
totalTstAcc_MAVVAR_C2C3 = totalTstAcc_MAVVAR_C2C3/10;
totalTstCM_MAV_C1rest = totalTstCM_MAV_C1rest/10;
totalTstCM_VAR_C1rest = totalTstCM_VAR_C1rest/10;
totalTstCM_MAVVAR_C1rest = totalTstCM_MAVVAR_C1rest/10;
totalTstCM_MAV_C2rest = totalTstCM_MAV_C2rest/10;
totalTstCM_VAR_C2rest = totalTstCM_VAR_C2rest/10;
totalTstCM_MAVVAR_C2rest = totalTstCM_MAVVAR_C2rest/10;
totalTstCM_MAV_C1C2 = totalTstCM_MAV_C1C2/10;
totalTstCM_VAR_C1C2 = totalTstCM_VAR_C1C2/10;
totalTstCM_MAVVAR_C1C2 = totalTstCM_MAVVAR_C1C2/10;
totalTstCM_MAV_C3rest = totalTstCM_MAV_C3rest/10;
totalTstCM_VAR_C3rest = totalTstCM_VAR_C3rest/10;
totalTstCM_MAVVAR_C3rest = totalTstCM_MAVVAR_C3rest/10;
totalTstCM_MAV_C1C3 = totalTstCM_MAV_C1C3/10;
totalTstCM_VAR_C1C3 = totalTstCM_VAR_C1C3/10;
totalTstCM_MAVVAR_C1C3 = totalTstCM_MAVVAR_C1C3/10;
totalTstCM_MAV_C2C3 = totalTstCM_MAV_C2C3/10;
totalTstCM_VAR_C2C3 = totalTstCM_VAR_C2C3/10;
totalTstCM_MAVVAR_C2C3 = totalTstCM_MAVVAR_C2C3/10;


% Create a confusion chart
%figure;
%confMatrix = round(TstCM_MAVVAR_C1C2);
%confusionchart(confMatrix, {'Class 1', 'Class 2'}, 'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
%title('Confusion Matrix');

% Optional: Adjust the appearance of the chart
%ax = gca;
%ax.FontSize = 12;
