%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pnet=zeros(1,T); %Net 
Qnet=zeros(1,T); %Net 
VBMax=zeros(1,ND);
VBMin=zeros(1,ND);
LOADVB=zeros(1,ND);
PlosDaB=zeros(1,ND);
VH=zeros(1,T); %House Voltage 
PlosB=zeros(1,T); %Losses
QPV=0; %PV Reactive Power 
for i=1:T
  Pnet(i)=PV(i)-D(i)-EV(i); 
  Qnet(i)=QPV-(D(i)*PF); 
  VH(i)=VT+((R*(Pnet(i)*1000)+X*(Qnet(i)*1000))/VT);
  PlosB(i)=(abs(Pnet(i)*1000/VH(i))^2)*R; %losses at each time-point
end
gk=0;
for b=1:ND
  VBMax(b)=max(VH(1+TD*gk:TD*b));
  VBMin(b)=min(VH(1+TD*gk:TD*b));
  LOADVB(b)=var(Pnet(1+TD*gk:TD*b),1);
  PlosDaB(b)=(sum(PlosB(1+TD*gk:TD*b))/1000)*1/tau; %total daily losses in kWh
  gk=gk+1;
end
if DIA==1
    DIAA=tau;
    DIX=0;
else
    DIAA=1;
    DIX=tau;
end
if PCN==0
ETOC=0;
end
%% Program 1: Conventional Rule-based Algorithm (CRBA)
if Prog==1
disp('Program 1 - Conventional Rule-based Algorithm (CRBA) = Started ')    
Pnet=zeros(1,T); %Net 
BESSC=zeros(1,T); %Charge
BESSD=zeros(1,T); %Discharge
BESSR=zeros(1,T); %Residual
BESSRR=zeros(1,ND);
SOC=ones(1,T)*SOCMIN; %SOC
SOC(1)=SOCI;
FT=zeros;FCAL=zeros;FDOD=zeros;FSOC=zeros;FC=zeros;FCYC=zeros;FD=zeros;L=zeros; %Degradation Model parameters
g=2;f=1;gf=0;k=1;
for i=1:T
   Pnet(i)=D(i)+EV(i)-PV(i); 
   %Charge
   %Charge at night low ToU
   if i==((TLS)*tau+gf) 
   BESSRR(k)=(SOCMAX-SOC(i))*BESS;
   end
   if (i>=((TLS*tau)+gf)&& i<=(ETOC*tau+gf)) && SOC(i)<SOCMAX && PCN(k)>0
         BESSC(i)=(PCN(k)*BESSRR(k)/RE)/(ETOC-TLS);
   end
      %Charge in normal time according to a Threshold
     if Pnet(i)<PTHC && SOC(i)<SOCMAX
      BESSR(i)=(SOCMAX-SOC(i))*BESS;
     if (PTHC-Pnet(i))*(1/tau)<=BESSR(i) && (PTHC-Pnet(i))<=BESSP 
      BESSC(i)=(PTHC-Pnet(i));
      elseif (PTHC-Pnet(i))>=BESSP  && BESSP*(1/tau)<=BESSR(i)
      BESSC(i)=BESSP;
     else 
      BESSC(i)=BESSR(i)*(tau);
     end
     end
   %Discharge
    if Pnet(i)>PTHD && SOC(i)>SOCMIN && (i>=(TLE*DIAA+gf)) && (i>=(ETOC*DIX+gf))
    BESSR(i)=(SOC(i)-SOCMIN)*BESS;
    if (Pnet(i)-PTHD)*(1/tau)<=BESSR(i) && (Pnet(i)-PTHD)<=BESSP
    BESSD(i)=(Pnet(i)-PTHD);
    elseif (Pnet(i)-PTHD)>=BESSP && BESSP*(1/tau)<=BESSR(i)
    BESSD(i)=BESSP; 
    else
    BESSD(i)=BESSR(i)*(tau);    
    end
    end
   SOC(g)=SOC(i)+(((BESSC(i)*(1/tau)*RE))/BESS)-(((BESSD(i)*(1/tau))/RE)/BESS);  
      if SOC(g)<SOCMIN
      SOC(g)=SOCMIN; 
      end
         if SOC(g)>SOCMAX
      SOC(g)=SOCMAX; 
   end
   %Update days
   if i==TD*k
      gf=TD*k;
      k=k+1;
   end 
g=g+1;
end
%Count Cycles (Rainflow Count Algorithm)
c=rainflow(SOC);
cyc=c(:,1);
S=c(:,4);
E=c(:,5);
DODD=c(:,2);
SOCD=c(:,3);
LL1=S(1:end-1);
LL2=E(1:end-1);
S1=SOC(LL1);
S2=SOC(LL2);
NCYC=round(abs(sum(S1-(S2)))/(DOD)); %Number of full cycles at DoD
NCYCPF=sum(cyc); %Number of all cycles partial and full from the rainflow algorithm 

% Degradation Model (Cycle Aging)
for f=1:length(S1)
    FT(f)=(exp((6.93*10^-2)*(20-25)*(25/20)));
    FC(f)=(exp(0.263*((DODD(f)/(abs(LL1(f)-LL2(f)))*1/tau)-1)));
    FDOD(f)=(((8.95*10^4)*(DODD(f)^(-0.486)))-7.28*10^4)^-1;
    FSOC(f)=(exp(1.04*(SOCD(f)-(0.5))));
    FCYC(f)=FDOD(f)*FSOC(f)*FC(f)*FT(f)*(abs(S1(f)-S2(f))/(DOD));
end
FD=sum(FCYC);
L=1-((0.0575*exp(FD*-121))+((1-0.0575)*exp(-FD)));
L2=FCYC(1:length(S1));
L2 = cumsum(L2);
LFL=(L2/max(L2));
LFL2=LFL*(BESS*L);
BESSn=[BESS BESS-LFL2];  %BESS residual Capcity 
BESSPer=100*BESSn/BESS;
end
%% Program 2: Proposed Day-ahead Scheduling Algorithm (PDSA)
if Prog==2
disp('Program 2 - Optimizaion (PDSA) = Started - Please wait, this program takes time according to the data resolution and number of days')    
SOC=ones(1,T)*SOCMIN; %SOC
SOC(1)=SOCI;
FT=zeros;FCAL=zeros;FDOD=zeros;FSOC=zeros;FC=zeros;FCYC=zeros;FD=zeros;L=zeros; %Degradation Model parameters
BESSC=zeros(1,T); %Charge
BESSD=zeros(1,T); %Discharge
BESSRR=zeros(1,ND);
BESSDD=zeros(1,ND);
SOCG=zeros(1,ND);
gf=0;
f=1;
 for k=1:ND
    SOCG(k)=SOC(TD*k+1-TD); 
    funcworhp
    BESSC(gf+1:TD*k)=XX(1:TD)/RE;
    BESSD(gf+1:TD*k)=XX(TD+1:TD*2)*RE;
    BESSD(TD*k)=0;
    BESSC(TD*k+1-TD)=0;
    g=2;
    for i=1:TD
   SOC(g+gf)=SOC(i+gf)+(((BESSC(i+gf)*(1/tau)*RE))/BESS)-(((BESSD(i+gf)*(1/tau))/RE)/BESS);  
   g=g+1;
    end
   
gf=TD*k;
 end
%Count Cycles (Rainflow Count Algorithm)
c=rainflow(SOC);
cyc=c(:,1);
S=c(:,4);
E=c(:,5);
DODD=c(:,2);
SOCD=c(:,3);
LL1=S(1:end-1);
LL2=E(1:end-1);
S1=SOC(LL1);
S2=SOC(LL2);
NCYC=round(abs(sum(S1-(S2)))/(DOD)); %Number of full cycles at DoD
NCYCPF=sum(cyc); %Number of all cycles partial and full from the rainflow algorithm 

% Degradation Model (Cycle Aging)
for f=1:length(S1)
    FT(f)=(exp((6.93*10^-2)*(20-25)*(25/20)));
    FC(f)=(exp(0.263*((DODD(f)/(abs(LL1(f)-LL2(f)))*1/tau)-1)));
    FDOD(f)=(((8.95*10^4)*(DODD(f)^(-0.486)))-7.28*10^4)^-1;
    FSOC(f)=(exp(1.04*(SOCD(f)-(0.5))));
    FCYC(f)=FDOD(f)*FSOC(f)*FC(f)*FT(f)*(abs(S1(f)-S2(f))/(DOD));
end
FD=sum(FCYC);
L=1-((0.0575*exp(FD*-121))+((1-0.0575)*exp(-FD)));
L2=FCYC(1:length(S1));
L2 = cumsum(L2);
LFL=(L2/max(L2));
LFL2=LFL*(BESS*L);
BESSn=[BESS BESS-LFL2];  %BESS residual Capcity 
BESSPer=100*BESSn/BESS;
end

%% Program 3: Proposed Rule-based Algorithm (PRBA)
%%%%%%%%%Inputs%%%%%%%%%%%%  
if Prog==3
disp('Program 3 -  Proposed Rule-based Algorithm (PRBA) = Started ')   
PTHD=zeros;
PTHC=zeros;
ABGDP=zeros(1,ND);
AGDPF=zeros(1,ND);
RESSD=zeros(1,ND);
EVSs=zeros(1,ND);
EVe=zeros(1,ND);
x=0;
for i=1:ND
vEe=find((EV(1+TD*x:TD*i)>0));
if sum(vEe)>0
    if min(vEe(vEe>16*tau)) ~= 0    
EVSs(i)=min(vEe(vEe>16*tau)); %Start time of EV
    else
EVSs(i)=min(vEe);
    end
EVe(i)=max(vEe); %End time of EV
end
if EVSs(i)<PVe(i) && PVe(i)>18*tau
EVSs(i)=18*tau;
elseif EVSs(i)<PVe(i) && PVe(i)<18*tau
    EVSs(i)=PVe(i);
end
ABGDP(i)=0.5*(AVGD/24)*PVL(i)*1/tau;
AGDPF(i)=FPV(i)-(ABGDP(i)); %Excess PV generation
x=x+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 
Pnet=zeros(1,T); %Net 
BESSC=zeros(1,T); %Charge
BESSD=zeros(1,T); %Discharge
BESSR=zeros(1,T); %Residual
BESSRD=zeros(1,T); %Residual
SOC=ones(1,T)*SOCMIN; %SOC
SOC(1)=SOCI;
RESB=zeros(1,ND);
FT=zeros(1,ND);FCAL=zeros(1,ND);FDOD=zeros(1,ND);FSOC=zeros(1,ND);FC=zeros(1,ND);FCYC=zeros(1,ND);FD=zeros(1,ND);L=zeros(1,ND); %Degradation Model parameters
% Main program 
g=2;ff=1;gf=0;k=1;
fj=0;gk=1;
for i=1:T
   Pnet(i)=D(i)+EV(i)-PV(i); 
   %Charge
   BESSR(i)=(SOCMAX-SOC(i))*BESS;
   BESSRD(i)=(SOC(i)-SOCMIN)*BESS;       
   if AGDPF(k)>((SOCMAX-SOC(1+TD*k-TD))*BESS)/RE   && Pnet(i)<0 
       PTHC(k)=((AGDPF(k)-(((SOCMAX-SOC(1+TD*k-TD))*BESS)/RE))/(PVL(k)/tau));
%        PTHC(k)=3.68/5;
    if abs(Pnet(i))>PTHC(k) && (abs(Pnet(i))-PTHC(k))*(1/tau)<BESSR(i) && (abs(Pnet(i))-PTHC(k))<BESSP 
    BESSC(i)=abs(Pnet(i))-PTHC(k);    
    elseif abs(Pnet(i))>PTHC(k) && BESSP*(1/tau)<BESSR(i) && (abs(Pnet(i))-PTHC(k))>BESSP 
    BESSC(i)=BESSP;      
     
    elseif abs(Pnet(i))>PTHC(k) && (abs(Pnet(i))-PTHC(k))*(1/tau)>BESSR(i) && (abs(Pnet(i))-PTHC(k))<BESSP || abs(Pnet(i))>PTHC(k) && BESSP*(1/tau)>BESSR(i) && (abs(Pnet(i))-PTHC(k))>BESSP  || abs(Pnet(i))-PTHC(k)>0 && abs(Pnet(i))*(1/tau)>BESSR(i) && (abs(Pnet(i)))<BESSP
    BESSC(i)=BESSR(i)*(tau);
    end 
   end  
   if AGDPF(k)<((SOCMAX-SOC(1+TD*k-TD))*BESS)/RE  
              PTHC(k)=0;
       if AGDPF(k)>0
          RESB(k)=(((SOCMAX-SOC(1+TD*k-TD))*BESS)/RE)-AGDPF(k); 
   elseif AGDPF(k)<0
       RESB(k)=(((SOCMAX-SOC(1+TD*k-TD))*BESS)/RE); 
       end
   if (i>=(TLE+gf)&& i<=(TLE*tau+gf)) && (RESB(k))/((TLE-TLS))<BESSR(i)*tau  && SOC(i)<1  && Pnet(i)<EVS*0.5
   BESSC(i) = (RESB(k))/((TLE-TLS));  
   end 
    if i>=(PVS(k)+gf)&& i<=(PVe(k)+gf) && Pnet(i)<0      
    if abs(Pnet(i))>PTHC(k) && (abs(Pnet(i))-PTHC(k))*(1/tau)<BESSR(i) && (abs(Pnet(i))-PTHC(k))<BESSP 
    BESSC(i)=abs(Pnet(i))-PTHC(k);
    elseif abs(Pnet(i))>PTHC(k) && BESSP*(1/tau)<BESSR(i) && (abs(Pnet(i))-PTHC(k))>BESSP 
    BESSC(i)=BESSP;
    elseif abs(Pnet(i))-PTHC(k)>0 && abs(Pnet(i))*(1/tau)<BESSR(i) && (abs(Pnet(i)))<BESSP
    BESSC(i)=abs(Pnet(i));  
    elseif abs(Pnet(i))>PTHC(k) && (abs(Pnet(i))-PTHC(k))*(1/tau)>BESSR(i) && (abs(Pnet(i))-PTHC(k))<BESSP || abs(Pnet(i))>PTHC(k) && BESSP*(1/tau)>BESSR(i) && (abs(Pnet(i))-PTHC(k))>BESSP || abs(Pnet(i))-PTHC(k)>0 && abs(Pnet(i))*(1/tau)>BESSR(i) && (abs(Pnet(i)))<BESSP
    BESSC(i)=BESSR(i)*(tau);
    end
   end
   end
   
   %Discharge
    if EVC(k)==1
           BESSRD(i)=(SOC(i)-SOCMIN)*BESS;
    if i>=(EVSs(k)+gf)  && i<=(EVe(k)*tau+gf) && Pnet(i)>EVS*0.5
    PTHD(gk)=EVS-((BESSRD(i)*RE)/((EVA(k)-((fj)/tau)))); %discharge trigger threshold     

    fj=fj+1;
    if Pnet(i)>abs(PTHD(gk)) && SOC(i)>SOCMIN 
    if (Pnet(i)-abs(PTHD(gk)))*(1/tau)<BESSRD(i) && (Pnet(i)-abs(PTHD(gk)))<BESSP
    BESSD(i)=(Pnet(i)-abs(PTHD(gk)));
    elseif (Pnet(i)-abs(PTHD(gk)))*(1/tau)<BESSRD(i) && (Pnet(i)-abs(PTHD(gk)))>BESSP 
    BESSD(i)=BESSP;
    elseif (Pnet(i)-abs(PTHD(gk)))*(1/tau)>BESSRD(i) && (Pnet(i)-abs(PTHD(gk)))>BESSP && BESSP*(1/tau)<BESSRD(i)
    BESSD(i)=BESSP;
    elseif (Pnet(i)-abs(PTHD(gk)))*(1/tau)>BESSRD(i) && BESSP>BESSRD(i)*tau
    BESSD(i)=BESSRD(i)*(tau);
    end
    end

    end
    if (Pnet(i)>0 && SOC(i)>SOCMIN && (i>=(EVe(k)+gf) && i<=(PE*tau+gf))) || ((Pnet(i)>0 && SOC(i)>SOCMIN && i>=(PVS(k)+gf) && i<=(PVe(k)+gf) && Season==1 && i<(EVSs(k)+gf))) || (Pnet(i)>EVS*0.5 && (i>(PVe(k)+gf)) && EVe(k)==TD && EVSs(k)>20*tau)
    if (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))<BESSP
    BESSD(i)=(Pnet(i));
    elseif (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))>BESSP
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && (Pnet(i))>BESSP && BESSP*(1/tau)<BESSRD(i)
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && BESSP>BESSRD(i)*tau
    BESSD(i)=BESSRD(i)*(tau);
    end
    end
    end
    
    if EVC(k)==2
    if Season==1
        BESSRD(i)=(SOC(i)-SOCMIN)*BESS;
    if Pnet(i)>0 && i<(PVe(k)+gf) && i>=(1+gf)
    if (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))<BESSP 
    BESSD(i)=(Pnet(i));
    elseif (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))>BESSP
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && (Pnet(i))>BESSP && BESSP*(1/tau)<BESSRD(i)
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && BESSP>BESSRD(i)*tau
    BESSD(i)=BESSRD(i)*(tau);
    end
    end

    if i>(PVe(k)+gf) && i<=(PE*tau+gf)
    RESSD(k)=((SOC(PVe(k)+gf)-SOCMIN)*BESS*RE)/(PE-(PVe(k)/tau));
    PTHDs=RESSD(k);
    if Pnet(i)>-PTHDs && SOC(i)>SOCMIN 
    if (Pnet(i)+PTHDs)*(1/tau)<BESSRD(i) && (Pnet(i)+PTHDs)<BESSP
    BESSD(i)=(Pnet(i)+PTHDs);
    elseif (Pnet(i)+PTHDs)*(1/tau)<BESSRD(i) && (Pnet(i)+PTHDs)>BESSP 
    BESSD(i)=BESSP;
    elseif (Pnet(i)+PTHDs)*(1/tau)>BESSRD(i) && (Pnet(i)+PTHDs)>BESSP && BESSP*(1/tau)<BESSRD(i)
    BESSD(i)=BESSP;
    elseif (Pnet(i)+PTHDs)*(1/tau)>BESSRD(i) && BESSP>BESSRD(i)*tau
    BESSD(i)=BESSRD(i)*(tau);
    end
    end    
    end
    end
    
    if Season==0
    if Pnet(i)>0 && SOC(i)>SOCMIN
    BESSRD(i)=(SOC(i)-SOCMIN)*BESS;
    if (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))<BESSP
    BESSD(i)=(Pnet(i));
    elseif (Pnet(i))*(1/tau)<BESSRD(i) && (Pnet(i))>BESSP
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && (Pnet(i))>BESSP && BESSP*(1/tau)<BESSRD(i)
    BESSD(i)=BESSP;
    elseif (Pnet(i))*(1/tau)>BESSRD(i) && BESSP>BESSRD(i)*tau
    BESSD(i)=BESSRD(i)*(tau);
    end
    end
    end
   end
      SOC(g)=SOC(i)+((BESSC(i)*(1/tau)*RE)/BESS)-((BESSD(i)*(1/tau)/RE)/BESS);  
    %Update days
   if i==TD*k
      gf=TD*k;
      k=k+1;
      fj=0;
   end 

      gk=gk+1;
g=g+1;
end
%Count Cycles (Rainflow Count Algorithm)
c=rainflow(SOC);
cyc=c(:,1);
S=c(:,4);
E=c(:,5);
DODD=c(:,2);
SOCD=c(:,3);
LL1=S(1:end-1);
LL2=E(1:end-1);
S1=SOC(LL1);
S2=SOC(LL2);
NCYC=round(abs(sum(S1-(S2)))/(DOD)); %Number of full cycles at DoD
NCYCPF=sum(cyc); %Number of all cycles partial and full from the rainflow algorithm 

% Degradation Model (Cycle Aging)
for f=1:length(S1)
    FT(f)=(exp((6.93*10^-2)*(20-25)*(25/20)));
    FC(f)=(exp(0.263*((DODD(f)/(abs(LL1(f)-LL2(f)))*1/tau)-1)));
    FDOD(f)=(((8.95*10^4)*(DODD(f)^(-0.486)))-7.28*10^4)^-1;
    FSOC(f)=(exp(1.04*(SOCD(f)-(0.5))));
    FCYC(f)=FDOD(f)*FSOC(f)*FC(f)*FT(f)*(abs(S1(f)-S2(f))/(DOD));
end
FD=sum(FCYC);
L=1-((0.0575*exp(FD*-121))+((1-0.0575)*exp(-FD)));
L2=FCYC(1:length(S1));
L2 = cumsum(L2);
LFL=(L2/max(L2));
LFL2=LFL*(BESS*L);
BESSn=[BESS BESS-LFL2];  %BESS residual Capcity 
BESSPer=100*BESSn/BESS;
end
%% Metrices Calculations
BESSNP=BESSC-BESSD;
%Curtailed Power
PNETC1=zeros(1,T);PNETC2=zeros(1,T);PCU1=zeros;PCU2=zeros;PSC1=zeros;PSC2=zeros;p=1;PDD1=zeros(1,T);PDD2=zeros(1,T);n=1;m=1;
for ip=1:T
 PDD1(ip)=D(ip)+EV(ip);
 PNETC1(ip)=D(ip)+EV(ip)-PV(ip)+BESSNP(ip);
 PNETC2(ip)=D(ip)+EV(ip)-PV(ip);
 if PNETC1(ip)<0 && PNETC1(ip)*-1>EXP
   PCU1(ip)= (PNETC1(ip)*-1)-EXP;
 end
 if PNETC2(ip)<0 && PNETC2(ip)*-1>EXP
   PCU2(ip)= (PNETC2(ip)*-1)-EXP;
 end
end
PCUA=sum(PCU1)*1/tau;%power curtailed after BESS
PCUB=sum(PCU2)*1/tau;%power curtailed after BESS
%% Electricity Bill Calculations
totold=zeros(1,ND);
totnew=zeros(1,ND);
totoldb=zeros(1,ND);
P=zeros(1,ND);
for j=1:ND
Pnetbb=zeros(1,TD); %Net
Pneto=zeros(1,TD); %Net
PNETNEW=zeros(1,TD); %Net
g=1;
for k=(1+TD*(j-1)):(TD*(j))
        PVVV(g)=PV(k);

    Pnetbb(g)=D(k)+EV(k); 
    Pneto(g)=D(k)+EV(k)-PV(k); 
    PNETNEW(g)=D(k)+EV(k)-PV(k)+BESSNP(k);
        BESSPOW(g)=BESSNP(k);
    
    if Pnetbb(g)>=PV(k)
    DESAPV(g)=PV(k);
    elseif Pnetbb(g)<PV(k)
    DESAPV(g)=Pnetbb(g);
    end 
    
    
    g=g+1;
end
Pneto(Pneto<-EXP)=0;
PNETNEW(PNETNEW<-EXP)=0;
% w/o BESS
PNETDO=max(0,Pneto);
PNETPVO=min(0,Pneto);
PNETPVO(PNETPVO<-EXP)=0;
EBO=PNETDO.*PRP+PNETPVO*EX;
totold(j)=((sum(EBO)*(1/tau))/100)+(SC/100);
EEr(j)=sum(PNETPVO)*1/tau;

% w/ BESS
PNETD=max(0,PNETNEW);
PNETPV=min(0,PNETNEW);
PNETPV(PNETPV<-EXP)=0;
EBN=PNETD.*PRP+PNETPV*EX;
totnew(j)=((sum(EBN)*(1/tau))/100)+(SC/100);
EErr(j)=sum(PNETPV)*1/tau;

EBBB=Pnetbb.*PRP;
totoldb(j)=((sum(EBBB)*(1/tau))/100)+(SC/100); %without PV or BESS
idx1(j)=find(PVVV > 0, 1, 'first');
idx2(j)=find(PVVV > 0, 1, 'last')+1;

DEDUPV=(DESAPV(idx1(j):idx2(j)));
DECOV(j)=sum(DEDUPV)*1/tau; %Demand capacity covered during PV period
BESSDUPV=(BESSPOW(idx1(j):idx2(j)));
BESSCOV(j)=sum(max(0,BESSDUPV))*1/tau; %BESS capacity charged during PV period
DECOVBESS(j)=abs(sum(min(0,BESSPOW)))*1/tau; %Demand covered by the BESS
DEneBESS(j)=abs(sum(max(0,BESSPOW)))*1/tau;%Demand needed by the BESS

DEMON(j)=sum(Pnetbb)*1/tau; %demand only
%save in file
if SaveR==1
range10 = strcat('A',num2str(j+1));
range11 = strcat('B',num2str(j+1));
range9 = strcat('C',num2str(j+1));
range1 = strcat('D',num2str(j+1));

writematrix(j,'Results.xls','Sheet',2,'Range', range10)
writematrix(totoldb(j),'Results.xls','Sheet',2,'Range', range11)
writematrix(totold(j),'Results.xls','Sheet',2,'Range', range9)
writematrix(totnew(j),'Results.xls','Sheet',2,'Range', range1)
end
end
RESBESSC=sum(BESSC)*1/tau+(sum(EEr)-sum(EErr));
%Self-Consumption(SC) and Self-Sufficiency (SS)

% SC1=(((sum(PV)*1/tau)+sum(EEr))/(sum(PV)*1/tau))*100 %SC before BESS
% SS1=(((sum(PV)*1/tau)+sum(EEr))/(sum(PDD1)*1/tau))*100 %SS before BESS
% 
% SC2=((((sum(PV)*1/tau)+sum(EErr)))/(((sum(PV)*1/tau))))*100 %SC AFTER BESS
% SS2=((((sum(PV)*1/tau)+sum(EErr))+(sum(BESSD)*1/tau))/((sum(PDD1)*1/tau)+(sum(BESSC)*1/tau)))*100 %SS after BESS


SC2=((sum(BESSCOV)+sum(DECOV))/(sum(PV)*1/tau))*100;
SS2=((sum(BESSCOV)+sum(DECOV)+sum(DECOVBESS))/(sum(DEneBESS)+sum(DEMON)))*100;

SC1=((sum(DECOV))/(sum(PV)*1/tau))*100;
SS1=((sum(DECOV))/(sum(DEMON)))*100;


BILL1=sum(round(totold,2)); %Elec bill w/o BESS
BILL2=sum(round(totnew,2)); %Elec bill w BESS
BILL0=sum(round(totoldb,2));%Elec bill w/o PV and BESS
disp('=============================================================================================')
disp('********************           SIMULATION RESULTS OF THE RBMT          **********************')
disp('=============================================================================================')
disp(['Electricity Bill w/o PV w/o BESS:       ',num2str(round(BILL0,3)),' £'])
disp(['Electricity Bill w/ PV w/o BESS:        ',num2str(round(BILL1,3)),' £'])
disp(['Electricity Bill w/ PV w/ BESS:         ',num2str(round(BILL2,3)),' £'])
disp(['Curtailed PV power w/o BESS:            ',num2str(round(PCUB,3)),' kWh'])
disp(['Curtailed PV power w/ BESS:             ',num2str(round(PCUA,3)),' kWh'])
disp(['Exported PV power w/o BESS:             ',num2str(round(sum(EEr*-1),3)),' kWh'])
disp(['Exported PV power w/ BESS:              ',num2str(round(sum(EErr*-1),3)),' kWh'])
disp(['PV Self-Consumption w/o BESS:           ',num2str(round(SC1,3)),' %'])
disp(['PV Self-Consumption w/ BESS:            ',num2str(round(SC2,3)),' %'])
disp(['Premises Self-Sufficiency w/o BESS:     ',num2str(round(SS1,3)),' %'])
disp(['Premises Self-Sufficiency w/ BESS:      ',num2str(round(SS2,3)),' %'])
disp('=============================================================================================')

%% Voltage / Load variance / Losses after BESS 
Pnet=zeros(1,T); %Net 
Qnet=zeros(1,T); %Net 
VHA=zeros(1,T); %House Voltage 
PlosA=zeros(1,T); %House Voltage 
LOADVA=zeros(1,ND);
for i=1:T
  Pnet(i)=PV(i)-D(i)-EV(i)-BESSNP(i); 
  Qnet(i)=QPV-D(i)*PF; 
  VHA(i)=VT+((R*(Pnet(i)*1000)+X*(Qnet(i))*1000)/VT);
  PlosA(i)=(abs(Pnet(i)*1000/VHA(i))^2)*R; %losses at each time-point

end
VAMax=zeros(1,ND);
VAMin=zeros(1,ND);
PlosDaA=zeros(1,ND);
gk=0;
for b=1:ND
      LOADVA(b)=var(Pnet(1+TD*gk:TD*b),1);
      VAMax(b)=max(VHA(1+TD*gk:TD*b));
      VAMin(b)=min(VHA(1+TD*gk:TD*b));
      PlosDaA(b)=(sum(PlosA(1+TD*gk:TD*b))/1000)*1/tau; %total daily losses in kWh
      gk=gk+1;
end
%Voltage and load variance
if SaveR==1
range212 = strcat('A',num2str(2));
range20122 = strcat('B',num2str(2));
range2A22 = strcat('C',num2str(2));
range2B22 = strcat('D',num2str(2));
range2BB22 = strcat('E',num2str(2));
range2BAB22 = strcat('F',num2str(2));
range2LB22 = strcat('G',num2str(2));
range2LA22 = strcat('H',num2str(2));
writematrix(VH','Results.xls','Sheet',4,'Range', range212)
writematrix(VHA','Results.xls','Sheet',4,'Range', range20122)
writematrix(VBMax','Results.xls','Sheet',4,'Range', range2A22)
writematrix(VBMin','Results.xls','Sheet',4,'Range', range2B22)
writematrix(VAMax','Results.xls','Sheet',4,'Range', range2BB22)
writematrix(VAMin','Results.xls','Sheet',4,'Range', range2BAB22)
writematrix(LOADVB','Results.xls','Sheet',4,'Range', range2LB22)
writematrix(LOADVA','Results.xls','Sheet',4,'Range', range2LA22)
end
%% Figures and Results
Pnet=D+EV-PV;
PNETW=D+EV-PV+BESSNP';

if SaveR==1
RES=[(1:T)' Pnet PNETW BESSNP' SOC(1:T)'];
range101 = strcat('A',num2str(2));
writematrix(RES,'Results.xls','Sheet',1,'Range', range101)
%Degrad and losses
range1012 = strcat('A',num2str(2));
range10122 = strcat('B',num2str(2));
rangelb22 = strcat('A',num2str(2));
rangela122 = strcat('B',num2str(2));
writematrix(BESSPer(end),'Results.xls','Sheet',3,'Range', range1012)
writematrix(BESSn(end),'Results.xls','Sheet',3,'Range', range10122)
writematrix(PlosDaB','Results.xls','Sheet',5,'Range', rangelb22)
writematrix(PlosDaA','Results.xls','Sheet',5,'Range', rangela122)
end

figure
yyaxis left
ylabel('Power [kW]')
hold on
plot(PNETW,'r-')
hold on
plot(Pnet,'k')
yyaxis right
ylabel('SOC [%]')
stairs(SOC,'b--')
legend({'w/ BESS','w/o BESS','BESS SOC'},'Location','northeast','AutoUpdate','off')
axis([1 T -inf inf])
set(gca,'FontSize',14,'FontName', 'Times New Roman')
xlabel('Time')

figure
newcolors = [0, 0.4470, 0.7410
0.8500, 0.3250, 0.0980];  
colororder(newcolors)
x=linspace(0,ND,(length(BESSn)));
yyaxis left
ylabel('BESS Capacity [kWh]')
hold on
Y1=plot(x,BESSn,'-');
hold on
yyaxis right
ylabel('SoH [%]')
Y2=plot(x,BESSPer,'--');
legend({'BESS Capacity', 'BESS SoH'},'Location','northeast','AutoUpdate','off')
xlim([0 ND])
set(gca,'FontSize',14,'FontName', 'Times New Roman')
xlabel('Days')
toc;

%% Optimization Function for Program 2
function funcworhp(x)
global BESSP TD  XX SOCMIN SOCMAX
wData.XL = [zeros(1,2*TD)] ;
wData.XU = [ones(1,2*TD)*BESSP]; 
% % % Define bounds for G
wData.GL = [ones(1,TD)*SOCMIN];
wData.GU = [ones(1,TD)*SOCMAX];
% Initial estimates for X and Lambda
wData.xInit = [zeros(1,length(wData.XL))];
wData.lambdaInit = zeros(size(wData.XL));
% and for Mu
wData.muInit = zeros(size(wData.GL));
% Initialise sparsity structure of derivates
% % according to the WORHP user manual (Coordinate Storage format)
wData.DFrow = int32([]);
wData.DGrow = int32([]);
wData.DGcol = int32([]);
wData.HMrow = int32([]);
wData.HMcol = int32([]);
wData.param = 'worhp.xml';
% The callback functions.
wCallback.f           = @objective;
wCallback.g           = @constraints;
wCallback.df          = [];
wCallback.dg          = [];
wCallback.hm          = [];
% Call the solver
[xFinal , ~, ~, ~] = worhp(wData, wCallback);

for i=1:2*TD
XX(i)=xFinal(i);
end
end

function f = objective(x)
global  RE TD D EV PV PRP gf EX 
Pnet=zeros(1,TD); %Net
Bill=zeros(1,TD); %Net
Pnet1=zeros(1,TD); %Net
Pnet2=zeros(1,TD); %Net
for i=1:TD
 Pnet(i)=D(i+gf)+EV(i+gf)-PV(i+gf)+((x(i))/RE)-(x(i+TD)*RE);
 Bill(i)=Pnet(i)*PRP(i);
  if Pnet(i)>0
     Pnet1(i)=Pnet(i);
      Bill(i)=Pnet1(i)*PRP(i);
  elseif Pnet(i)< 0
    Pnet2(i)=Pnet(i);
    Bill(i)=(Pnet2(i))*EX;
  end
end 
f=0.7*(sum(Bill)/1000)+0.05*((var(Pnet,1)/8))+0.25*((var(Pnet2,1)/2));
end

function g = constraints(x)
global tau TD  BESS k SOCG
SOC=zeros(1,TD);
SOC(1)=SOCG(k);
tg=2;
 for i=1:TD
   SOC(tg)=SOC(i)+(((x(i)*(1/tau)))/BESS)-((((x(i+TD))*(1/tau)))/BESS)  ;
   tg=tg+1;
 end
 g=zeros;
 for j=1:TD
 g(j) = (SOC(j));
 end
end

