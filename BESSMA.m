% ======================================================== % 
% This code was developed by Ahmed A.Raouf Mohamed: ------ %
% ------ Ra2ooof@gmail.com / amohamed06@qub.ac.uk -------- %
% https://pure.qub.ac.uk/en/persons/ahmed-mohamed -------- %
% -------------------------------------------------------- %
% ---  This work is part of INTERREG VA SPIRE2 Project --- %
% EPIC Research Cluster ---:--- Queen's University Belfast %
% Copyright @2020 ---------------------------------------- %
% -------------------------------------------------------- %
% Please check Codeguide.pdf for details on how to run the %
% code. -------------------------------------------------- %
% The measurements should be inserted in data.m ---------- %
% -------------------------------------------------------- %
% ======================================================== %    

clc; close all; clear;
tic;
warning('off')
format long g
global  RE T D EV PV EX PRP BESSU tau gf TD BESSRR k BESSDD XX BESSP SOCG BESS SOCMIN SOCMAX
%% Main Inputs 
SaveR=0; %if 1=save results in excel files, other values=don't save : (saving results will reduce excution time) 
Prog=1;  %1 for the conventional rule-based method (CRBA), 2 for the proposed day-ahead scheduling (PDSA), 3 for the proposed rule-based algorithm (PRBA) 
DataRes=10; %Data resolution 10 for 10 minutes reso, 30 for 30 minutes reso, 60 for 60 minutes(1 hour) reso and so on...
%% Call Data
Data % Insert the simulation data as explained in the Data.m file
T=length(Profile(:,1));
ND=round((DataRes/60)*(T/24)); %Number of days
TD=T/ND; %Length of one day
D=Profile(:,1); %Demand
PV=Profile(:,2);
EV=Profile(:,3);
tau=TD/(24); % {Time interval=1/tau}
PE=24; %end time of peak
%% BESS Inputs
BESS=9.8; %Actual BESS Capacity
DOD=0.8; %BESS DOD
SOCMAX=1; %Max SoC
BESSP=5; %BESS Power
RE=0.95*0.95; %round trip efficency = 0.95(BESS) * 0.95(Inverter)
SOCMIN=SOCMAX-DOD; %Min SoC 
BESSU=BESS*DOD; %Usable BESS Capacity
%% ToU Tariff Data e.g. Economy 7  (https://powerni.co.uk/plan-prices/compare-our-plans/economy-7-unit-rates/)
HR=17.19; %Day Rate 8am-1am 
LR=9.59; % Night Rate 1am-8am 
TLS=1;  %Night rate start time
TLE=8; %Night Rate end time
EX=5; %Export Tariff = 5p/kWh 
TPR=[HR LR LR LR LR LR LR LR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR];%Tariff Profile
PRP = repelem(TPR,tau);
%% Network Inputs for the voltage calculations %% Voltage / Load variance / Losses before BESS 
VT=240; %Transformer Voltage
PF=0.95; % Power Factor
R=0.240979; %Resistance from transformer bus to the household
X=0.0030569; %Reactance from the transformer bus to the household
%% Program 1 Inputs: Conventional Rule-based Algorithm (CRBA) 
PTHD=0; % Specify the Upper threshold for BESS Discharge 
PTHC=-0; % Specify the Lower threshold for BESS Charge  
%% Program 3 Inputs: Proposed Rule-based Algorithm (PRBA)
EVS=3.8; % EV charger in kW
AVGD=10.3; %average daily consumption in kWh
Season=1; % 1 for High PV season - Summer, 0 for Low PV season - Winter
PTHDn=0; %Normal upper threshold
EVA=3.5;%Average Cahrging hours of the EV 
EVC=1;%Electric Vehicle charging next day ? yes:1 , No:0
FPV=9.8; %Forecasted PV daily generation in day ahead 
PVS=7; %start time of PV
PVe=18; %end time of PV
PVL=(PVe-PVS); %Number of hours of PV generation period
EVA=ones(1,ND)*EVA;EVC=ones(1,ND)*EVC;FPV=ones(1,ND)*FPV; PVL=ones(1,ND)*PVL;PVS=ones(1,ND)*PVS*tau;PVe=ones(1,ND)*PVe*tau;
%% Extra part to extract the value of EVC, EVA, and FPV automatically From the data
%  Also, it captures the start and end time and the number
%  of hours forthe PV instead of inserting them manually and the season.
% x=0;FPV=zeros(1,ND);EVC=zeros(1,ND);PVS=zeros(1,ND);PVe=zeros(1,ND);PVL=zeros(1,ND);EVA=zeros(1,ND);
% for i=1:ND 
% FPV(i)=sum(PV(1+TD*x:TD*i))/tau;  
% if sum(EV(1+TD*x:TD*i))>0
%     EVC(i)=1;
% elseif sum(EV(1+TD*x:TD*i))==0
%     EVC(i)=2;
% end
% v=find(PV(1+TD*x:TD*i)>0);
% PVS(i)=min(v); %Start time of PV
% PVe(i)=max(v); %End time of PV
% PVL(i)=(PVe(i)-PVS(i));
% vE=((EV(1+TD*x:TD*i)>0));
% EVA(i)=round(sum(vE)/tau,1);% Number of charging hours EV
% if i<4320
% Season=0;
% elseif i>=5760 && i <12960
%     Season=1;
% elseif i>=12960
%     Season=0;
% end
% x=x+1;
% end
%% 
MAINCODE