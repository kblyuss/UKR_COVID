function covid_ukraine

clear all
close all

global Tinf Ntot Tinc Tsev Tdeath K1 K2 K3 K4 pa pm ps pf mix_matrix pop_dist R0

Tinc=4.9; % incubation period
Tinf=5; % infectious period

% numbers of stages in Gamma distribution

K1=2;
K2=4;
K3=24;
K4=2;

Tsev=22.6; % duration of severe cases (hospitalisation)
Tdeath=9; % duration until death after hospitalisation

pop_dist=[1905861	2307946	2083524	1795352	2161728	2866088	3469194	3220048	2957318	2832930	2730967	3068098	2740347	2304448	1449678	3139329];
Ntot=sum(pop_dist);

zz=load('AM3.dat'); % load percentages of symptomatic/asymptomatic/hospitalisation rates/mortality rates

psymp(1:16)=zz(1:16,7)/100;
pasymp(1:16)=zz(1:16,8)/100;
hosp_rates(1:16)=zz(1:16,4)/100;
mort_rates1(1:16)=zz(1:16,3)/100;

pa=pasymp; % proportion of asymptomatic
ph=hosp_rates; % proportion of those who require hospitalisation
pm=psymp-ph; % proportion of those who are symptomatic but will recover without requiring hospitalisation

%%% Load mixing matrices %%%%

school=load('ukr_school.txt');
work=load('ukr_work.txt');
other=load('ukr_other.txt');
hhome=load('ukr_home.txt');

%total cases/deaths in Ukraine by 6 July 2020
cases=[1	1	1	1	1	1	1	1	1	3	3	3	3	7	14	16	26	41	47	73	73	97	145	156	237	311	418	480	645	794	897	1072	1225	1251	1319	1462	1668	1892	2227	2511	2801	3104	3420	3764	4161	4662	5106	5449	5710	6125	6592	7170	7647	8125	8617	9009	9410	9866	10406	10861	11411	11913	12331	12697	13184	13691	14195	14710	15232	15648	16023	16425	16847	17330	17858	18291	18616	18876	19230	19706	20148	20580	20986	21245	21584	21905	22382	22811	23204	23672	24012	24340	24823	25411	25964	26514	26999	27462	27856	28381	29070	29753	30506	31154	31810	32476	33234	34063	34984	35825	36560 37241	38074	39014 40008	41117	42065	42982	43628	44334	44998	45887 46763 47677 48500 49043];
deaths=[0	0	0	0	0	0	0	0	0	0	1	1	1	1	2	2	3	3	3	3	3	3	5	5	5	8	8	11	17	20	22	27	32	37	38	45	52	57	70	74	84	93	99	108	116	125	133	141	151	161	174	187	193	201	209	220	239	250	261	272	279	288	303	316	327	340	361	376	391	408	425	439	456	476	497	514	535	548	564	579	588	605	617	623	644	658	669	679	696	708	718	727	735	747	762	777	788	797	810	833	854	870	880	889	901	912	943	966	985	994	1002 1012	1035	1051 1067	1086	1110	1129 1147	1159	1173	1185 1212 1227 1249 1262];

%total cases/deaths in Ukraine by 13 July 2020
cases_updated=[	1	1	1	1	1	1	1	1	1	3	3	3	3	7	14	16	26	41	47	73	73	97	145	156	237	311	418	480	645	794	897	1072	1225	1251	1319	1462	1668	1892	2227	2511	2801	3104	3420	3764	4161	4662	5106	5449	5710	6125	6592	7170	7647	8125	8617	9009	9410	9866	10406	10861	11411	11913	12331	12697	13184	13691	14195	14710	15232	15648	16023	16425	16847	17330	17858	18291	18616	18876	19230	19706	20148	20580	20986	21245	21584	21905	22382	22811	23204	23672	24012	24340	24823	25411	25964	26514	26999	27462	27856	28381	29070	29753	30506	31154	31810	32476	33234	34063	34984	35825	36560	37241	38074	39014	40008	41117	42065	42982	43628	44334	44998	45887	46763 47677 48500 49043 49607 50414 51224 52043 52843 53521 54133];
deaths_updated=[	0	0	0	0	0	0	0	0	0	0	1	1	1	1	2	2	3	3	3	3	3	3	5	5	5	8	8	11	17	20	22	27	32	37	38	45	52	57	70	74	84	93	99	108	116	125	133	141	151	161	174	187	193	201	209	220	239	250	261	272	279	288	303	316	327	340	361	376	391	408	425	439	456	476	497	514	535	548	564	579	588	605	617	623	644	658	669	679	696	708	718	727	735	747	762	777	788	797	810	833	854	870	880	889	901	912	943	966	985	994	1002	1012	1035	1051	1067	1086	1110	1129	1147	1159	1173	1185	1212 1227 1249 1262 1283 1306 1327 1345 1372 1383 1398];

N=length(cases);
N1=N-6;

avs=zeros(N1,1);
avd=zeros(N1,1);

for i=1:N1
    for j=-3:1:3
     avs(i)=avs(i)+cases(i+j+3);
     avd(i)=avd(i)+deaths(i+j+3);
    end;
end;
avs=avs/7;

for j=1:N1
    if avs(j)>100
        break;
    end;
end;

days2=linspace(1,N1-j+1,N1-j+1);
avs2=avs(j:N1);

NNL=length(cases_updated)-length(cases);

diffs=cases(j+4:N1+4)-cases(j+3:N1+3);
cases2=cases(j+3:N1+3);
cases2upd=cases_updated(N1+3+1:N1+3+NNL);
days2upd=days2(end)+linspace(1,NNL,NNL);
deaths_upd=deaths_updated(N1+3+1:N1+3+NNL);
LL=length(avs2);
diffs2=avs2(2:LL)-avs2(1:LL-1);
diffs_upd=cases_updated(N1+4:N1+3+NNL)-cases_updated(N1+3:N1+2+NNL);

deaths2=deaths(j+3:N1+3);
dediffs=deaths(j+4:N1+4)-deaths(j+3:N1+3);
dediffs_upd=deaths_updated(N1+3+1:N1+3+NNL)-deaths_updated(N1+3:N1+2+NNL);

% compute with SEIR stages

IC0=zeros((K1+K2+K3+K4+5)*16,1);

IC0(1:16)=pop_dist;
IC0(9)=IC0(9)-1;
IC0(25)=1;

R0=2.5;

% coefficients for scaling mixing matrices
a1=1; % school
a2=1; % work
a3=1; % other
a4=1; % home

mix_matrix=a1*school+a2*work+a3*other+a4*hhome;

pfc=0.16;
mort_rates=(pfc/0.0269)*mort_rates1;

pf=mort_rates; % proportion of those who will die
ps=ph-pf; % proportion of those who will go to hospital but will recover

opts=odeset('AbsTol',1e-8,'RelTol',1e-8);
[t,x]=ode45(@seir_age,[0 50],IC0,opts);

% break once the total cases reach 100
for ii=1:length(t)
    if sum(x(ii,17:end))>100
        break;
    end;
end;

IC1d=x(ii,:);

totE=sum(IC1d(16*(K1+1)+1:16*(K1+K2+1)));

IC1d(16*(K1+K2+K3+K4+4)+7)=1;
IC1d(16*(K1+K2+K3+K4+4)+11)=1;
IC1d(16*(K1+K2+K3+K4+4)+15)=1;

for i=1:16
    IC1d(16*(K1+K2+K3+4)+i)=2*totE/(3*K4*16); %IF1
end;

for j=1:(K4-1)
    for i=1:16
        IC1d(16*(K1+K2+K3+4+j)+i)=2*totE/(3*K4*16); % IF2 to IF_K3
    end;
end;

for i=16*(K1+1)+1:16*(K1+K2+1)
    IC1d(i)=1*IC1d(i)/3;
end;

xpld=1;
ypld=avs2(1);
Ind=sum(IC1d(16*(K1+1)+1:16*(K1+K2+1)));
ypld_d=deaths2(1);

for k=1:(length(diffs2)-1)    
    R0d=diffs2(k)*Tinf/Ind;
    % update CFR along the timeline, fix 30 days after 100 cases
    if k<10
        pfc=0.3;
        mort_rates=(pfc/0.0269)*mort_rates1;
        pf=mort_rates; % proportion of those who will die
        ps=ph-pf; % proportion of those who will go to hospital but will recover
    elseif k<30
        pfc=0.3-0.28*k/32;
        mort_rates=(pfc/0.0269)*mort_rates1;
        pf=mort_rates; % proportion of those who will die
        ps=ph-pf; % proportion of those who will go to hospital but will recover
    else
        mort_rates=1.65*mort_rates1;
        pf=mort_rates; % proportion of those who will die
        ps=ph-pf; % proportion of those who will go to hospital but will recover
    end;
    if k<12 % from 25 Mar until 6 April
        a1=0;a2=0.4;a3=0.4;a4=1;
        other_kg=other;
        other_kg(1,1)=0;
        mix_matrix=a1*school+a2*work+a3*other_kg+a4*hhome;
    elseif k<58 % shielding of over 60+ introduced on 6 April + other restrictions until 25 May 
        % 60+ at 0.2 for work/other compared to the rest, who are at 0.4
        work_sh=work;
        other_sh=other;
        for ii=1:12
            for jj=13:16
                work_sh(ii,jj)=0.5*work_sh(ii,jj);
                other_sh(ii,jj)=0.5*other_sh(ii,jj);
            end;
        end;
        for ii=1:12
            for jj=13:16
                work_sh(jj,ii)=0.5*work_sh(jj,ii);
                other_sh(jj,ii)=0.5*other_sh(jj,ii);
            end;
        end;
        for ii=13:16
            for jj=13:16
                work_sh(ii,jj)=0.5*work_sh(ii,jj);
                other_sh(ii,jj)=0.5*other_sh(ii,jj);
            end;
        end;
        other_sh(1,1)=0;
        mix_matrix=a1*school+a2*work_sh+a3*other_sh+a4*hhome;
    elseif k<72 % some restrictions are lifted on 25 May, but 60+ are still shielded
        other_sh(1,1)=0.5*other(1,1); % opening of kindergardens
        a1=0.1;a2=0.6;a3=0.6;a4=1;
        mix_matrix=a1*school+a2*work_sh+a3*other_sh+a4*hhome;
    else % 60+ not shielded any longer
        a1=0.1;a2=0.6;a3=0.6;a4=1;
        other_kg(1,1)=0.5*other(1,1);
        mix_matrix=a1*school+a2*work+a3*other_kg+a4*hhome;
    end;

    R0=R0d;
    [td,xd]=ode45(@seir_age,[0 1],IC1d,[]);
    R0x(k)=k;
    R0y(k)=R0d;
    xpld=[xpld (k+1)];
    ypld=[ypld sum(xd(end,17:end))];
    ypld_d=[ypld_d sum(xd(end,16*(K1+K2+K3+K4+4)+1:16*(K1+K2+K3+K4+5)))];
    IC1d=xd(end,:);
    Ind=sum(xd(end,16*(K1+1)+1:16*(K1+K2+1)));
end;

xpld_end=xpld(end);

ICfull_d=IC1d;

R0=R0d;

% compute short-term forecast

[t2,x2]=ode45(@seir_age,[0:1:10],ICfull_d,opts);

Ytot_d=zeros(length(t2),1);
for i=1:length(t2)
    Ytot_d(i)=sum(x2(i,17:end));
end;

xpld2=t2'+xpld_end;
ypld2=Ytot_d';
ypld_d2=zeros(length(t2),1);
for i=1:length(ypld_d2)
    ypld_d2(i)=sum(x2(i,16*(K1+K2+K3+K4+4)+1:16*(K1+K2+K3+K4+5)));
end;

%%%%%%%%%% Plot total casess %%%%%%
figure
subplot(2,2,1)
plot(days2-1,cases2,'ko','MarkerSize',8); % plot raw data
hold on
ylabel('Total cases','FontSize',24);
axis([0 108 0 5.6e4]);
set(gca,'LineWidth',2,'FontSize',24,'xtick',[0 20 40 60 80 100],'xticklabel',{'25 Mar','14 Apr','4 May','24 May','13 Jun','3 Jul'},'ytick',[0 2e4 4e4]);

plot(xpld,ypld,'r','LineWidth',2); % plot numerical solution

plot(xpld2,ypld2,'g','LineWidth',2); % plot forecast

plot(days2upd-1,cases2upd,'b*','MarkerSize',8);
% 
xpld=[xpld xpld2(2:end)]; % combined time points
ypld=[ypld Ytot_d(2:end)']; % combined total cases
ypld_d=[ypld_d ypld_d2(2:end)']; % combined total deaths

%%%%%%%%%% Plot new cases %%%%%%
subplot(2,2,3)

b1=bar(days2,diffs,'FaceColor',[0.5 0.5 0.5]);
b1.FaceAlpha=0.2;
hold on
ylabel('New cases','FontSize',24);
axis([0 108 0 1200])
set(gca,'LineWidth',2,'FontSize',24,'xtick',[0 20 40 60 80 100],'xticklabel',{'25 Mar','14 Apr','4 May','24 May','13 Jun','3 Jul'},'ytick',[0 400 800 1200]);

NNd=length(ypld);
ypld2=ypld(2:NNd)-ypld(1:(NNd-1));

plot(xpld(1:length(ypld2)),ypld2,'r-','LineWidth',2);%'MarkerSize',8) days3d

plot(xpld2-1,ypld2(end-length(xpld2)+1:end),'g','LineWidth',2); % plot forecast

b2=bar(days2upd,diffs_upd,'FaceColor','b');
b2.FaceAlpha=0.2;

%%%%%%%%%% Plot new deaths %%%%%%
subplot(2,2,4)

NNd_d=length(ypld_d);
ypld2_d=ypld_d(2:NNd_d)-ypld_d(1:(NNd_d-1));

LLd_d=length(ypld2_d);

hold on

zpl=linspace(1,LLd_d,LLd_d);
plot(zpl(1:end),ypld2_d(1:end),'r-','LineWidth',2); %,'MarkerSize',8,

b1=bar(days2,dediffs,'FaceColor',[0.5 0.5 0.5]);
b1.FaceAlpha=0.2;

axis([0 108 0 35]);
set(gca,'LineWidth',2,'FontSize',24,'Box','on','xtick',[0 20 40 60 80 100],'xticklabel',{'25 Mar','14 Apr','4 May','24 May','13 Jun','3 Jul'},'ytick',[0 10 20 30]);
ylabel('New deaths','FontSize',24);

b2=bar(days2upd,dediffs_upd,'FaceColor','b');
b2.FaceAlpha=0.2;

plot(xpld2-1,ypld2_d(end-length(xpld2)+1:end),'g','LineWidth',2); % plot forecast

%%%%%%%%% Plot total deaths %%%%%%%%%
subplot(2,2,2)
plot(xpld,ypld_d,'r-','LineWidth',2);
hold on
axis([0 108 0 1450])
set(gca,'LineWidth',2,'FontSize',24,'Box','on','xtick',[0 20 40 60 80 100],'xticklabel',{'25 Mar','14 Apr','4 May','24 May','13 Jun','3 Jul'},'ytick',[0 400 800 1200]);
plot(days2-1,deaths2,'ok','MarkerSize',8);
plot(days2upd-1,deaths_upd,'b*','MarkerSize',8);

plot(xpld2,ypld_d(end-length(xpld2)+1:end),'g','LineWidth',2); % plot forecast

ylabel('Total deaths','FontSize',24);

%%% Plot age distribution of cases %%%%%%%%%

figure

subplot(1,2,1);

cases_age=zeros(16,1);
for i=1:16
    for j=1:(K1+K2+K3+K4+4)
        cases_age(i)=cases_age(i)+x2(end,16*j+i);
    end;
end;

%%% plot age distribution of cases

z=zeros(8,1);
z(1)=cases_age(1)+cases_age(2);
z(2)=cases_age(3)+cases_age(4);
z(3)=cases_age(5)+cases_age(6);
z(4)=cases_age(7)+cases_age(8);
z(5)=cases_age(9)+cases_age(10);
z(6)=cases_age(11)+cases_age(12)
z(7)=cases_age(13)+cases_age(14);
z(8)=cases_age(15)+cases_age(16);
ages_for_cases=[1;2;3;4;5;6;7;8];

actual_cases=[1982;2418;5452;9011;9575;11265;8764;4926];

zsum=sum(z);
z=z/zsum;
acsum=sum(actual_cases);
actual_cases=actual_cases/acsum;

yy=[z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8);actual_cases(1) actual_cases(2) actual_cases(3) actual_cases(4) actual_cases(5) actual_cases(6) actual_cases(7) actual_cases(8)];
b=bar(ages_for_cases,100*yy');
b(1).BarWidth=1;
b(1).LineWidth=2;
b(2).BarWidth=1;
b(2).LineWidth=2;
% b.FaceAlpha=0.2;
axis([0.5 8.5 0 25])
set(gca,'xticklabel',{'0-10','10-20','20-30','30-40','40-50','50-60','60-70','70+'},'FontSize',20,'LineWidth',2,'ytick',[0 10 20],'yticklabel',{'0', '10%', '20%'});
title('Total cases','FontSize',24,'FontWeight','Normal');
legend('Simulation','Data');
xlabel('Age','FontSize',24);

%%% Plot age distribution of deaths %%%%%%%%%

deaths_age=zeros(16,1);
for i=1:16
    deaths_age(i)=x2(end,16*(K1+K2+K3+K4+4)+i);
end;

subplot(1,2,2)
z=zeros(8,1);
z(1)=deaths_age(1)+deaths_age(2)+deaths_age(3)+deaths_age(4);
z(2)=deaths_age(5)+deaths_age(6);
z(3)=deaths_age(7)+deaths_age(8);
z(4)=deaths_age(9)+deaths_age(10);
z(5)=deaths_age(11)+deaths_age(12)
z(6)=deaths_age(13)+deaths_age(14);
z(7)=deaths_age(15)+deaths_age(16);

actual_deaths=[4;6;33;104;257;415;571];

zsum=sum(z);
z=z/zsum;
acsum=sum(actual_deaths);
actual_deaths=actual_deaths/acsum;

yy=[z(1) z(2) z(3) z(4) z(5) z(6) z(7);actual_deaths(1) actual_deaths(2) actual_deaths(3) actual_deaths(4) actual_deaths(5) actual_deaths(6) actual_deaths(7)];
ages_for_deaths=[1;2;3;4;5;6;7];
b=bar(ages_for_deaths,100*yy');%,'FaceColor','b');
b(1).BarWidth=1;
b(1).LineWidth=2;
b(2).BarWidth=1;
b(2).LineWidth=2;
axis([0.5 7.5 0 50])
set(gca,'xticklabel',{'0-20','20-30','30-40','40-50','50-60','60-70','70+'},'FontSize',20,'LineWidth',2,'ytick',[0 20 40],'yticklabel',{'0', '20%', '40%'});
title('Total deaths','FontSize',24,'FontWeight','Normal');
legend('Simulation','Data');
xlabel('Age','FontSize',24);

%%% Next function defines the model

function tata=seir_age(t,x)

global Tinf Tinc Tsev Tdeath K1 K2 K3 K4 R0 pa pf ps pm mix_matrix

tata=zeros((K1+K2+K3+K4+5)*16,1);

NN=zeros(16,1);

for i=1:16
    for j=0:(K1+K2+K3+K4+4)
        NN(i)=NN(i)+x(16*j+i);
    end;
end;

r1=max(real(eig(mix_matrix)));

Icur=zeros(16,1);

for i=1:16
    for j=1:K2
        Icur(i)=Icur(i)+x(16*(K1+1)+16*(j-1)+i);
    end;
end;

for i=1:16
    for j=1:16
        tata(i)=tata(i)-R0*x(i)*mix_matrix(i,j)*Icur(j)/(Tinf*NN(j)*r1); % S
    end;
end;

for i=1:16
    tata(16+i)=-tata(i)-K1*x(16+i)/Tinc; % E1, the first age group of the incubating
end;

for i=1:16
    for j=1:(K1-1)
        tata(16*(j+1)+i)=K1*x(16*j+i)/Tinc-K1*x(16*(j+1)+i)/Tinc; % E2 to E_K1
    end;
end;

for i=1:16
    tata(16*(K1+1)+i)=K1*x(16*K1+i)/Tinc-K2*x(16*(K1+1)+i)/Tinf; % I1, the first age group of infected
end;

for i=1:16
    for j=1:(K2-1)
        tata(16*(K1+j+1)+i)=K2*x(16*(K1+j)+i)/Tinf-K2*x(16*(K1+j+1)+i)/Tinf; % I2 to I_K2
    end;
end;

for i=1:16
    tata(16*(K1+K2+1)+i)=K2*pa(i)*x(16*(K1+K2)+i)/Tinf; %Ra
end;

for i=1:16
    tata(16*(K1+K2+2)+i)=K2*pm(i)*x(16*(K1+K2)+i)/Tinf; %Ra
end;

for i=1:16
    tata(16*(K1+K2+3)+i)=K2*ps(i)*x(16*(K1+K2)+i)/Tinf-K3*x(16*(K1+K2+3)+i)/Tsev; %IS1
end;

for j=1:(K3-1)
    for i=1:16
        tata(16*(K1+K2+3+j)+i)=K3*x(16*(K1+K2+2+j)+i)/Tsev-K3*x(16*(K1+K2+3+j)+i)/Tsev; % IS2 to IS_K3
    end;
end;

for i=1:16
    tata(16*(K1+K2+K3+3)+i)=K3*x(16*(K1+K2+K3+2)+i)/Tsev; %RS
end;

for i=1:16
    tata(16*(K1+K2+K3+4)+i)=K2*pf(i)*x(16*(K1+K2)+i)/Tinf-K4*x(16*(K1+K2+K3+4)+i)/Tdeath; %IF1
end;

for j=1:(K4-1)
    for i=1:16
        tata(16*(K1+K2+K3+4+j)+i)=K4*x(16*(K1+K2+K3+3+j)+i)/Tdeath-K4*x(16*(K1+K2+K3+4+j)+i)/Tdeath; % IF2 to IF_K3
    end;
end;

for i=1:16
    tata(16*(K1+K2+K3+K4+4)+i)=K4*x(16*(K1+K2+K3+K4+3)+i)/Tdeath; %DD
end;



