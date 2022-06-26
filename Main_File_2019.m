
close all
clear
clc
warning off
addpath(genpath('functions'))
global aa;
%part 1
aa='Subject_1.txt';
fileID = fopen(aa,'r');
fileIDD=fileID;
this_line=0;
var1={};
while this_line ~=-1
this_line=fgetl(fileID);
if this_line ~=-1
    var1=[var1;this_line];
end 
end
var2=var1(250);    %enter row number you want to analyze 
dlmwrite('myFile.txt',var2,'delimiter','','roffset',1);
fileID = fopen('myFile.txt','r');
commas = char(44);
sizeA = [1 Inf];
[A] = fscanf(fileID,['%d' commas],sizeA);
fs=500;
A=A';
Z=A(129,1);

%part 2
j=fft(A);
oo=length(j);
L=(0:oo-1)*(fs/oo);
subplot(331)
stem(L,20*log((j)))
title('spectrum of origional data')
xlabel('frequency(Hz)-->')
ylabel('amplitude(db)-->')

%%%filters
%part 3
%%delta

Fp=.5;
Fs=3.75;
Rp=0.057501127785;
wn=[Fp, Fs]/(fs/2);
fs=500;
Rs=0.0001;
[Or,F,po,w] = firpmord(wn, [1 0], [Rp, Rs]);
b1 = firpm(Or, F, po, w);
F1 = dfilt.dffir(b1);
x1=filter(F1,A);
subplot(332)
ts=1;
t=0:ts:128;
stem(t,x1,'r')
title('delta time domain')
xlabel('time(S)-->')
ylabel('amplitude(db)-->')

j=fft(x1);
oo=length(j);
L=(0:oo-1)*(fs/oo);
subplot(333)
stem(L,20*log((j)));
title('delta spectrum')
xlabel('frequency(Hz)-->')
ylabel('amplitude(db)-->')

%theta
Fs1=3.75;
fp1=4;
fp2=7;
fs2=7.75;
Rs1=.001;
Rs2=.0001;
Rp=0.057501127785;
wn=[Fs1 fp1 fp2 fs2]/(fs/2);
[Or, F, po, w] = firpmord(wn, [0 1 0], [Rs1, Rp,Rs2]);
b1 = firpm(Or, F, po, w);
F1 = dfilt.dffir(b1);
x1=filter(F1,A);
subplot(334)
ts=1;
t=0:ts:128;
stem(t,x1,'r')
title('theta in time domain')
xlabel('time(S)-->')
ylabel('amplitude(db)-->')

j=fft(x1);
oo=length(j);
L=(0:oo-1)*(fs/oo);
subplot(335)
stem(L,20*log((j)))
title('theta spectrum')
xlabel('frequency(Hz)-->')
ylabel('amplitude(db)-->')


%alpha
Fs1=7.75;
fp1=8;
fp2=13;
fs2=13.5;
Rs1=.001;
Rs2=.0001;
Rp=0.057501127785;
wn=[Fs1 fp1 fp2 fs2]/(fs/2);
[Or, F, po, w] = firpmord(wn, [0 1 0], [Rs1, Rp,Rs2]);
b1 = firpm(Or, F, po, w);
F1 = dfilt.dffir(b1);
x1=filter(F1,A);
subplot(336)
ts=1;
t=0:ts:128;
stem(t,x1,'r')
title('Alpha in time domain')
xlabel('time(S)-->')
ylabel('amplitude(db)-->')

j=fft(x1);
oo=length(j);
L=(0:oo-1)*(fs/oo);
subplot(337)
stem(L,20*log((j)))
title('Alpha spectrum')
xlabel('frequency(Hz)-->')
ylabel('amplitude(db)-->')

%beta

Fs1=13.5;
fp1=14;
fp2=29.5;
fs2=30;
Rs1=.001;
Rs2=.0001;
Rp=0.057501127785;
wn=[Fs1 fp1 fp2 fs2]/(fs/2);
[Or, F, po, w] = firpmord(wn, [0 1 0], [Rs1, Rp,Rs2]);
b1 = firpm(Or, F, po, w);
F1 = dfilt.dffir(b1);
x1=filter(F1,A);
subplot(338)
ts=1;
t=0:ts:128;
stem(t,x1,'r')
title('beta in time domain')
xlabel('time(S)-->')
ylabel('amplitude(db)-->')

j=fft(x1);
oo=length(j);
L=(0:oo-1)*(fs/oo);
subplot(339)
stem(L,20*log((j)))
title('beta spectrum')
xlabel('frequency(Hz)-->')
ylabel('amplitude(db)-->')

%% Define bands for filtering
Subject_No='A01';
order = 4;   % Order of the filter
band=[8 12; 14 30];

do_preprocess=1;
if(do_preprocess==1)
    %==========================================================================
    %%                   Load DATASET Training Data
    disp('########   Loading Data    ##################')
    disp('Subject: A01')
    load('DataSet2a_A01T_Old.mat');
    load('A01T.mat');
    EEG_SIG_Tr=EEG_SIG;
    Labels_T=classlabel;
    Triggers_T=Triggers_A01T;
    
    load('DataSet2a_A01E_Old.mat');
    load('A01E.mat');
    EEG_SIG_Ts=EEG_SIG;
    Labels_E=classlabel;
    Triggers_E=Triggers_A01E;
    disp('#######  Extracting Features for Each band and Applying CSP ##########')
    [Features]=f_Feature_Extraction(EEG_SIG_Tr, EEG_SIG_Ts,Labels_T, Labels_E,Triggers_T, Triggers_E, band, order, Subject_No, SamplingRate_A01T);
    
    save('Features_A01.mat', 'Features')
else
    load('Features_A01.mat', 'Features')
end
Count=0;
[No_Filters,dim]=size(band);
for k=1:No_Filters
    Temp=band(k,:);
    for l=1:No_Filters
        if(find((Features{1,l}.band )==Temp))
            TEMP_TrX=Features{1,l}.Train_X;
            TEMP_Tr_Y=Features{1,l}.Train_Y;
            TEMP_TsX=Features{1,l}.Test_X;
            TEMP_Ts_Y=Features{1,l}.Test_Y;
            Count=Count+1;
            if(Count==1)
                Train_X=TEMP_TrX;
                Train_Y=TEMP_Tr_Y;
                Test_X=TEMP_TsX;
                Test_Y=TEMP_Ts_Y;
            else
                Train_X=cat(1,Train_X,TEMP_TrX);
                Test_X=cat(2,Test_X,TEMP_TsX);
            end
        else
        end
    end
end
%==========================================================================
%#################### Training #################################
% Train the Classifiers on Training Data
size(Train_X)
disp('#######  Training The SVM Classsifier ##########')
model=fitcsvm(Train_X',Train_Y');

%% Perform 10-Fold Cross_validation
CVSVMModel = crossval(model)
kloss=kfoldLoss(CVSVMModel, 'LossFun', 'classiferror', 'Folds', 10);
kloss_mat=(1-kloss)*100;
disp(kloss_mat)
%% Evaluation or Testing
run('SVM_Data.p');
[Label]=f_Adaptive_Learning_A(Test_X,model);

%%  SVM Classification accuracy---------------------------------------------
cp = classperf(Test_Y, Label.SVM);
acc=cp.CorrectRate*148;
disp(acc)

