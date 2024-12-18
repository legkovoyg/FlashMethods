clear all
clc
%%   Задание параметров модели
%mixture 1
%z = [0.228	0.605	83.578	7.4	3.345	0.755	0.962	0.338	0.316	0.356	0.483	0.593	0.275	0.297	0.163	0.112	0.09	0.035	0.039	0.021	0.008	0.001];

%mixture 2
%z = [0.105	0.574	52.685	6.788	3.872	1.008	1.406	0.575	0.56	0.73	1.569	2.193	1.153	1.484	1.08	1.036	1.338	0.834	1.55	2.2	2.96	14.3];

%mixture 3
%z = [0.258	0.457	55.06	6.964	5.087	0.808	2.442	0.614	0.829	0.675	1.839	1.933	1.074	1.536	1.192	1.108	1.446	0.927	1.681	2.21	3.08	8.78];

%mixture 4 h=3575.0
%z = [0.032	0.506	88.114	6.483	1.93	0.532	0.484	0.207	0.147	0.208	0.232	0.338	0.189	0.227	0.1	0.073	0.081	0.028	0.039	0.026	0.011	0.013];

%mixture 5 h=3652.8
%z = [0.092	0.644	84.412	7.031	2.729	0.675	0.791	0.307	0.268	0.347	0.458	0.61	0.33	0.411	0.226	0.174	0.181	0.074	0.102	0.072	0.037	0.029];

%mixture 4 h=3666.9
%z = [0.11	0.672	83.51	7.13	2.906	0.705	0.866	0.33	0.299	0.382	0.52	0.682	0.368	0.461	0.263	0.205	0.211	0.089	0.122	0.088	0.047	0.034];

%mixture 27 h=3751.7
%z = [0.302	0.833	75.697	7.624	4.187	0.903	1.477	0.509	0.577	0.672	1.155	1.381	0.718	0.932	0.676	0.558	0.541	0.28	0.381	0.298	0.201	0.098];

%mixture 28 h=3758.8
%z = [0.328	0.846	74.922	7.646	4.3	0.918	1.537	0.525	0.606	0.7	1.227	1.454	0.753	0.979	0.723	0.599	0.577	0.302	0.41	0.322	0.222	0.104];

%mixture 29 h=3765.9
%z = [0.355	0.858	74.172	7.663	4.412	0.931	1.597	0.54	0.636	0.727	1.3	1.526	0.786	1.024	0.769	0.639	0.61	0.324	0.438	0.345	0.241	0.107];

%mixture 30 h=3772.9
z = [0.386	0.871	73.444	7.677	4.522	0.944	1.656	0.554	0.665	0.754	1.372	1.596	0.818	1.066	0.814	0.677	0.641	0.345	0.464	0.366	0.259	0.109];

%z = [0.386	0.871	73.444	7.677	4.522	0.944	1.656	0.554	0.665	0.754	1.372	1.596	0.818	1.066	0.814	0.677	0.641	0.345	0.464	0.366	0.259	0.109];

%mixture Methane for test H Enthalpy
%z = [0.00001	0.00001	  0.99999  0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001	0.00001];


sum(z);
%z=z./100;
z=z./sum(z);

if (sum(z)~=1.0) 
    %fprintf('ОШИБКА В ЗАДАНИИ СОСТАВА СМЕСИ!!! \n'); 
end
N = size(z,2);                                      % количество компонент N

mass = [28.014	44.01	16.043	30.07	44.097	58.124	58.124	72.151	72.151	86.178	91.925	105.718	120.343	140.56	169.086	199.098	236.892	269.752	310.637	372.292	448.232	701.123];
Pkr = [3.394388	7.376459	4.600155	4.883864	4.245517	3.647701	3.799688	3.384255	3.374123	2.968823	3.574564	3.122868	2.771944	2.50968	2.19141	1.972617	1.792949	1.694122	1.573101	1.488414	1.424685	1.351062];
Tkr = [126.2	304.19995	190.6	305.39996	369.79996	408.1001	425.2001	460.3999	469.6	507.4001	546.6401	569.9407	591.9464	600.8699	633.8047	664.5171	699.4725	727.1604	767.4636	813.2253	865.1088	1039.5869];
Vkr = [89.8	94	99.00002	148	203	263	255	306	304.0001	370.0001	425.4445	471.2443	523.4991	600.5144	712.9235	838.2478	1003.597	1147.239	1337.256	1629.407	1997.867	3509.115];
w = [0.04	0.225	0.008	0.098	0.152	0.176	0.193	0.227	0.251	0.296	0.436319	0.472542	0.510256	0.616047	0.692018	0.768374	0.859284	0.931732	0.972036	1.076324	1.174966	1.109354];
cpen=[0.92	3.03	0.63	2.63	5.06	7.29	7.86	10.93	12.18	17.98	6.72	13.03	19.41	15.12	19.66	20.81	17.47	10.18	9.11	-14.21	-49.83	-186.03];

%Константы для энтальпии идеального-газового состояния
Cp1_id = zeros(1,N);
Cp2_id = zeros(1,N);
Cp3_id = zeros(1,N);
Cp4_id = zeros(1,N);

% Константы по Reid, 1977
 Cp1_id = [31.1000000000	19.8000000000	19.2500000000	5.4100000000	-4.2200000000	-1.3900000000	9.4900000000	-9.5200000000	-3.6300000000	-4.4100000000	-5.1500000000	-6.1000000000	3.1400000000	-7.9100000000	-9.3300000000	-11.0000000000	-13.0000000000	-15.5000000000	-17.2000000000	-20.7000000000	-25.0000000000	-30.1400000000];
 Cp2_id = [-0.0135600000	0.0734300000	0.0521200000	0.1781000000	0.3063000000	0.3847000000	0.3313000000	0.5066000000	0.4873000000	0.5819000000	0.6761000000	0.7712000000	0.6774000000	0.9608000000	1.1490000000	1.3380000000	1.5290000000	1.7740000000	2.0020000000	2.3810000000	2.8540000000	3.4220000000];
 Cp3_id = [0.0000267900	-0.0000560200	0.0000119700	-0.0000693700	-0.0001586000	-0.0001846000	-0.0001108000	-0.0002729000	-0.0002580000	-0.0003119000	-0.0003651000	-0.0004195000	-0.0001928000	-0.0005288000	-0.0006347000	-0.0007423000	-0.0008537000	-0.0010150000	-0.0011230000	-0.0013390000	-0.0016100000	-0.0018810000];
 Cp4_id = [-0.0000000117	0.0000000172	-0.0000000113	0.0000000087	0.0000000322	0.0000000290	-0.0000000028	0.0000000572	0.0000000531	0.0000000649	0.0000000766	0.0000000886	-0.0000000298	0.0000001131	0.0000001359	0.0000001598	0.0000001850	0.0000002205	0.0000002445	0.0000002925	0.0000003526	0.0000004246];

% Константы по формуле 11, 2019
%Cp1_id = [31.1000000000	19.8000000000	19.2500000000	5.4100000000	-4.2200000000	-1.3900000000	9.4900000000	-9.5200000000	-3.6300000000	2.070882024	2.680427097	3.061998984	3.448838617	3.963983093	4.661621733	5.369993367	6.246422599	7.027781774	7.977955941	9.402230224	11.19448541	17.11742253];
%Cp2_id = [-0.0135600000	0.0734300000	0.0521200000	0.1781000000	0.3063000000	0.3847000000	0.3313000000	0.5066000000	0.4873000000	0.537931391	0.57139549	0.657236586	0.748346359	0.874392967	1.052392918	1.239796259	1.475872416	1.68103	1.9364032	2.321552548	2.795743977	4.375101534];
%Cp3_id = [0.0000267900	-0.0000560200	0.0000119700	-0.0000693700	-0.0001586000	-0.0001846000	-0.0001108000	-0.0002729000	-0.0002580000	-0.000224192	-0.000234127	-0.000269476	-0.000307147	-0.000359429	-0.000433511	-0.000511723	-0.000610382	-0.000695956	-0.000802662	-0.000963668	-0.001161573	-0.001821109];
%Cp4_id = [-0.0000000117	0.0000000172	-0.0000000113	0.0000000087	0.0000000322	0.0000000290	-0.0000000028	0.0000000572	0.0000000531	0	0	0	0	0	0	0	0	0	0	0	0	0];

cpen=cpen/1000;
c = zeros(N,N);
T0 = 273.15; % Tempreture for calculate Cp, H in K
R = 0.00831675;                % универсальная газовая МПа*м3/(кмоль*К)
V_pkr=z*Vkr';
T_pkr=z*Tkr';
Control_phase = V_pkr*T_pkr^2;

%% Задание температуры и давления
T = 373.15;    % Температура, К
P = 3;    % Давление, МПа

Ttest=382.15:50:382.15 %100:10:200;
Ptest=0.1:0.6:0.1;

Ptest = [5]

%Ptest(1)=0.1;
fid2 = fopen('W(P).txt','w');
fid5 = fopen('Mesh_Mix2.txt','w');
fid10 = fopen('All.txt','w');
fidH = fopen('Entalphy.txt','w');

step=0;
Mesh=zeros(size(Ttest,2)*size(Ptest,2),3);

tic
for hT=1:size(Ttest,2)
    T=Ttest(hT);
    for hP=1:size(Ptest,2)
        P=Ptest(hP);
        fprintf(fid10,'Температура =  %3.2f K \n', T); 
        fprintf(fid10,'Давление =  %3.2f МПа \n', P); 
        step=step+1;

        %% Расчетный модуль Соаве-Редлиха-Квонга
        ac_i = 0.42747*R^2*Tkr.^2./Pkr;                       
        psi_i = 0.48 + 1.574*w - 0.176*w.^2;
        alpha_i = (1+psi_i.*(1-(T./Tkr).^.5)).^2;
        a_i = ac_i.*alpha_i;
        b_i = 0.08664*R*Tkr./Pkr;
        c_i = cpen;

%% Шаг 0 Начало алгоритма

%  начальное приближение К0
K_i = (exp(5.373*(1+w).*(1-Tkr/T)).*Pkr/P).^(1.0);


aw=0;
bw=0;
ave = 0;
for i=1:N
    for j=1:N
        aw = aw + z(i)*z(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end

bw=z*b_i';
Aw = aw*P/(R^2*T^2);
Bw = bw*P/(R*T);
cw=cpen*z';
Cw = cw*P/(R*T);

Biw = b_i*P/(R*T);
Ciw = c_i*P/(R*T);

% SRK-Peneloux
ww = roots([1 (3*Cw-1) (3*Cw^2-Bw^2-2*Cw-Bw+Aw) (Cw^3-Bw^2*Cw-Cw^2-Bw*Cw+Aw*Cw-Aw*Bw)]);
Z_v = ww(imag(ww)==0);
Z_v=max(Z_v);
Z_init=Z_v; 

%% Шаг 4 Расчет летучестей состава z
avvv=zeros(1,N);
for i=1:N
    avv = 0;
    for j=1:N
        avv = avv + z(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    avvv(i)=avv;
end
% SRK-Peneloux
fz_i = exp(log(z*P) - log(Z_v + Cw - Bw) + (Biw-Ciw)/(Z_v + Cw - Bw) - (Aw/Bw)*((2*avvv/aw) - (b_i/bw)) * log( (Z_v + Bw + Cw) / (Z_v + Cw )) -...
    (Aw/Bw)*(Biw+Ciw)/(Z_v + Bw + Cw ) + (Aw/Bw)*Ciw/(Z_v + Cw ));


m=0;
indx=0;
eps_f=1;
Ri_v=1;
TS_v_flag=0;
TS_l_flag=0;

%% Часть 1 Проверка газовой фазы
%%  Паровая фаза
while (m<30)

Yi_v=z.*K_i;
Sv=sum(Yi_v);
y_i=Yi_v/Sv;
Yi_v
Sv
y_i
aw=0;
bw=0;
% for i=1:N
%     for j=1:N
%         aw = aw + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
%     end
% end

% Оптимизация циклов для случая c(i,j)=0
aw = (y_i*(a_i.^0.5)')^2;


bw=y_i*b_i';
Aw = aw*P/(R^2*T^2);
Bw = bw*P/(R*T);

cw=cpen*y_i';
Cw = cw*P/(R*T);

Biw = b_i*P/(R*T);
Ciw = c_i*P/(R*T);

% SRK-Peneloux
ww_y = roots([1 (3*Cw-1) (3*Cw^2-Bw^2-2*Cw-Bw+Aw) (Cw^3-Bw^2*Cw-Cw^2-Bw*Cw+Aw*Cw-Aw*Bw)]);

Z_v = ww_y(imag(ww_y)==0);
Z_v=Z_v(Z_v>0);
Z_v=max(Z_v);

%% Шаг 4 Расчет летучестей в паровой фазе

avvv=zeros(1,N);
% for i=1:N
%     avv = 0;
%     for j=1:N
%         avv = avv + y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
%     end
%     avvv(i)=avv;
% end

% Оптимизация циклов для случая c(i,j)=0
avvv=(a_i.^0.5)*(aw)^0.5;



% SRK-Peneloux
fw_i = exp(log(y_i*P) - log(Z_v + Cw - Bw) + (Biw-Ciw)/(Z_v + Cw - Bw) - (Aw/Bw)*((2*avvv/aw) - (b_i/bw)) * log( (Z_v + Bw + Cw) / (Z_v + Cw )) -...
    (Aw/Bw)*(Biw+Ciw)/(Z_v + Bw + Cw ) + (Aw/Bw)*Ciw/(Z_v + Cw ));


Ri=fz_i./(Sv*fw_i);
Ri_v=sum((Ri-1).^2);

if (Ri_v<10^(-12))
   %fprintf('Сходимость достигнута по газу!!! \n');
   Sv;
   m=30;
end

K_i=K_i.*Ri;
TS_v=sum(log(K_i).^2);

if (TS_v<10^(-4))
   %fprintf('TS по газу найдено!!! \n') ;
   TS_v_flag=1;
   m=30;
end

m=m+1;
end

K_iv=K_i;

%%% Проверка на жидкость (генерация начального приближения)

%  начальное приближение К0
K_i = (exp(5.373*(1+w).*(1-Tkr/T)).*Pkr/P).^1.0;

aw=0;
bw=0;
for i=1:N
    for j=1:N
        aw = aw + z(i)*z(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end

bw=z*b_i';
Aw = aw*P/(R^2*T^2);
Bw = bw*P/(R*T);
cw=cpen*z';
Cw = cw*P/(R*T);

Biw = b_i*P/(R*T);
Ciw = c_i*P/(R*T);

% SRK-Peneloux
ww = roots([1 (3*Cw-1) (3*Cw^2-Bw^2-2*Cw-Bw+Aw) (Cw^3-Bw^2*Cw-Cw^2-Bw*Cw+Aw*Cw-Aw*Bw)]);
Z_v = ww(imag(ww)==0);
Z_v=max(Z_v);


%% Шаг 4 Расчет летучестей состава z
avvv=zeros(1,N);
for i=1:N
    avv = 0;
    for j=1:N
        avv = avv + z(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    avvv(i)=avv;
end

% SRK-Peneloux
fz_i = exp(log(z*P) - log(Z_v + Cw - Bw) + (Biw-Ciw)/(Z_v + Cw - Bw) - (Aw/Bw)*((2*avvv/aw) - (b_i/bw)) * log( (Z_v + Bw + Cw) / (Z_v + Cw )) -...
    (Aw/Bw)*(Biw+Ciw)/(Z_v + Bw + Cw ) + (Aw/Bw)*Ciw/(Z_v + Cw ));


%% Часть 2 Проверка жидкой фазы

ml=0;
Ri_l=1;

while (ml<30)

    Yi_l=z./K_i;
    Sl=sum(Yi_l);
    x_i=Yi_l/Sl;
    Yi_l
    Sl
    x_i
    al=0;
    bl=0;

    % for i=1:N
    %     for j=1:N
    %         al = al + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    % end
    
    % Оптимизация циклов для случая c(i,j)=0
    al = (x_i*(a_i.^0.5)')^2;

    bl=x_i*b_i';
    Al = al*P/(R^2*T^2);
    Bl = bl*P/(R*T);
    cl = cpen*x_i';
    Cl = cl*P/(R*T);
    
    Bil = b_i*P/(R*T);
    Cil = c_i*P/(R*T);
    
    ww_x = roots([1 (3*Cl-1) (3*Cl^2-Bl^2-2*Cl-Bl+Al) (Cl^3-Bl^2*Cl-Cl^2-Bl*Cl+Al*Cl-Al*Bl)]);
    
    Z_l = ww_x(imag(ww_x)==0);
    Z_l=Z_l(Z_l>0);
    Z_l=min(Z_l);
    
    %% Шаг 6 Расчет летучестей в жидкой фазе
    
    alll=zeros(1,N);
    % for i=1:N
    %     all = 0;
    %     for j=1:N
    %         all = all + x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    %     alll(i)=all;
    % end
    
    % Оптимизация циклов для случая c(i,j)=0
    alll = (a_i.^0.5)*(al)^0.5;


    % SRK-Peneloux
    fl_i = exp(log(x_i*P) - log(Z_l + Cl - Bl) + (Bil-Cil)/(Z_l + Cl - Bl) - (Al/Bl)*((2*alll/al) - (b_i/bl)) * log( (Z_l + Bl + Cl) / (Z_l + Cl )) -...
        (Al/Bl)*(Bil+Cil)/(Z_l + Bl + Cl ) + (Al/Bl)*Cil/(Z_l + Cl ));
    
    Ri=Sl*fl_i./fz_i;
    Ri_l=sum((Ri-1).^2);

    xoy = sum(fz_i);
    xoy

    if (Ri_l<10^(-12))
        %fprintf('Сходимость достигнута по жидкости!!! \n');
        Sl;
        ml=30;
    end
    
    K_i=K_i.*Ri;
    TS=sum(log(K_i).^2);
    
    if (TS<10^(-4))
        %fprintf('TS по жидкости найдено!!! \n');
        TS_l_flag=1;
        ml=30;
    end
    
    ml=ml+1;
    K_i
    TS


end

K_il=K_i;

        if ((TS_l_flag==1 && TS_v_flag==1) ||  (Sv<=1 && TS_l_flag==1) || (Sl<=1 && TS_v_flag==1) || (Sv<1 && Sl<=1))
            %fprintf('Стабильное состояние!!! \n');
            Stable=1;
            TestPTF = 1; 
        else
            %fprintf('Нестабильное состояние!!! \n');
            Stable = 0;
            TestPTF = 0; 
        end

%Mesh(step,1)=T;
%Mesh(step,2)=P;
%Mesh(step,3)=Stable;
%fprintf(fid5,'%3.5f %3.5f %3.5f \t \n', Mesh(step,1), Mesh(step,2), Mesh(step,3));
Stable
TestPTF
K_iv
K_il
Sv
Sl
if (Stable==0 || TestPTF==1)

    %% Расчетный модуль Соаве-Редлиха-Квонга
    ac_i = 0.42747*R^2*Tkr.^2./Pkr;                       % Коэффициенты
    psi_i = 0.48 + 1.574*w - 0.176*w.^2;
    alpha_i = (1+psi_i.*(1-(T./Tkr).^.5)).^2;
    a_i = ac_i.*alpha_i;
    b_i = 0.08664*R*Tkr./Pkr;
    c_i=cpen;

    Kst_v=sum((K_iv-1).^2);
    Kst_l=sum((K_il-1).^2);
    if (Kst_l>Kst_v)
        K_i=K_il;
    else
        K_i=K_iv;
    end

m=0;
eps_f=1;

while (eps_f>0.000001 & m<400)

    %% Шаг 1 Нахождение общей доли пара

    [W] = findroot(z,K_i);

    %% Шаг 2 Нахождение мольных долей xi, yi
    x_i = z./(1 + W*(K_i-1)); 
    y_i = K_i.*x_i;

    %% Шаг 3 Нахождение z-фактора
    %  Паровая фаза

    aw=0;
    bw=0;
    % for i=1:N
    %     for j=1:N
    %         aw = aw + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    % end

    aw = (y_i*(a_i.^0.5)')^2;
    bw=y_i*b_i';
    Aw = aw*P/(R^2*T^2);
    Bw = bw*P/(R*T);
    cw=cpen*y_i';
    Cw = cw*P/(R*T);
    
    Biw = b_i*P/(R*T);
    Ciw = c_i*P/(R*T);
    
    ww_y = roots([1 (3*Cw-1) (3*Cw^2-Bw^2-2*Cw-Bw+Aw) (Cw^3-Bw^2*Cw-Cw^2-Bw*Cw+Aw*Cw-Aw*Bw)]);
    
    Z_v = ww_y(imag(ww_y)==0);
    Z_v=Z_v(Z_v>0);
    Z_v=max(Z_v);


    %% Шаг 4 Расчет летучестей в паровой фазе
    
    avvv=zeros(1,N);
    % for i=1:N
    %     avv = 0;
    %     for j=1:N
    %         avv = avv + y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    %     avvv(i)=avv;
    % end
    
    % Оптимизация циклов для случая c(i,j)=0
    avvv=(a_i.^0.5)*(aw)^0.5;

    % SRK-Peneloux
    fw_i = exp(log(y_i*P) - log(Z_v + Cw - Bw) + (Biw-Ciw)/(Z_v + Cw - Bw) - (Aw/Bw)*((2*avvv/aw) - (b_i/bw)) * log( (Z_v + Bw + Cw) / (Z_v + Cw )) -...
        (Aw/Bw)*(Biw+Ciw)/(Z_v + Bw + Cw ) + (Aw/Bw)*Ciw/(Z_v + Cw ));
    
    %% Шаг 5 Нахождение z-фактора
    % В жидкой фазе
    
    al=0;
    bl=0;
    % for i=1:N
    %     for j=1:N
    %         al = al + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    % end
    
    al = (x_i*(a_i.^0.5)')^2;

    bl=x_i*b_i';
    Al = al*P/(R^2*T^2);
    Bl = bl*P/(R*T);
    cl = cpen*x_i';
    Cl = cl*P/(R*T);
    
    Bil = b_i*P/(R*T);
    Cil = c_i*P/(R*T);
    
    ww_x = roots([1 (3*Cl-1) (3*Cl^2-Bl^2-2*Cl-Bl+Al) (Cl^3-Bl^2*Cl-Cl^2-Bl*Cl+Al*Cl-Al*Bl)]);
    
    Z_l = (ww_x(imag(ww_x)==0));
    Z_l=Z_l(Z_l>0);
    Z_l=min(Z_l);
    
    
    %% Шаг 6 Расчет летучестей в жидкой фазе
    
    alll=zeros(1,N);
    % for i=1:N
    %     all = 0;
    %     for j=1:N
    %         all = all + x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    %     end
    %     alll(i)=all;
    % end
    
    % Оптимизация циклов для случая c(i,j)=0
    alll = (a_i.^0.5)*(al)^0.5;

    % SRK-Peneloux
    fl_i = exp(log(x_i*P) - log(Z_l + Cl - Bl) + (Bil-Cil)/(Z_l + Cl - Bl) - (Al/Bl)*((2*alll/al) - (b_i/bl)) * log( (Z_l + Bl + Cl) / (Z_l + Cl )) -...
        (Al/Bl)*(Bil+Cil)/(Z_l + Bl + Cl ) + (Al/Bl)*Cil/(Z_l + Cl ));
    
    %% Корректировка распределения Ki
    
    df_lv = zeros(1,N);
    if (m<=N)
    for t=1:N
        if (fl_i(t)~=0)
            K_i(t) = K_i(t)*(fl_i(t)/fw_i(t))^1;
            df_lv(t) = (fl_i(t)/fw_i(t))-1;
        end
    end
    end
    
    Rr = fl_i./fw_i;
     
    if (m>1)
         Crit3 = sum((Rr-1).^2);
         Crit1 = Crit3/sum((Rr_old-1).^2);
         Crit2 = abs(W-W_old);
    end

    if (m>N & (Crit1 > 0.8) & (Crit2<0.1) &  (Crit3<0.001))
        
         for t=1:N
            if (fl_i(t)~=0)
            K_i(t) = K_i(t)*(fl_i(t)/fw_i(t))^(6.0); %^2
            df_lv(t) = (fl_i(t)/fw_i(t))-1;
            end
         end

    elseif (m>N)
        for t=1:N
            if (fl_i(t)~=0)
            K_i(t) = K_i(t)*(fl_i(t)/fw_i(t))^1;
            df_lv(t) = (fl_i(t)/fw_i(t))-1;
            end
         end
     end
     W_old=W;
     Rr_old = Rr;
     eps_f = max(abs(df_lv));
     m=m+1;

     if ((m>5) & (abs(W)>2))
        m=400;
        Stable=1;
     else
        Stable=0;
     end

     if ((m>80) & (abs(W-0.5)>0.501))
        m=400;
        Stable=1;
     else
        Stable=0;
     end

end

% Запись в файл результатов: давление Р и доля пара W
%fprintf(fid2,'%3.5f %3.9f \t \n', P, m);

end

if (Stable==1 || (abs(W-0.5)>0.5))

%% Старая проверка фазы    
%     K_i = (exp(5.373*(1+w).*(1-Tkr/T)).*Pkr/P).^(1.0);
%     ZK_Mult = sum(z.*K_i)-1;
%     ZK_Div = 1-sum(z./K_i);
% 
% %% Проверка на однофазное состояние жидкости
%  if (ZK_Mult<0)    % & ZK_Div<0
%       %fprintf('Состояние жидкости!!! \n') 
%       W=0;
%       x_i=z;
%       y_i=0;
%  end
% 
% %% Проверка на однофазное состояние газ
% if (ZK_Mult>0)   % & ZK_Div>0
%     %fprintf('Состояние газа!!! \n') 
%     W=1;
%     y_i=z;
%     x_i=0;    
% end

%% Новая проверка фазы

% if (W>1)
%     fprintf('Состояние газа!!! \n')
%      W=1;
%      x_i=0;
%      y_i=z;
%      Z_v = Z_init;
%      Z_l = 0;
% elseif (W<0)
%    fprintf('Состояние жидкость!!! \n') 
%        W=0;
%        x_i=z;
%        y_i=0;
%        Z_v = 0;
%        Z_l = Z_init;
% end

Volume = 1000*(Z_init*R*T/P);
if (T_pkr<260)
    %fprintf('Состояние газа!!! \n')
    W=1;
    x_i=zeros(1,N);
    y_i=z;
    Z_v = Z_init;
    Z_l = 0;
else
    if (Volume*T^2 > Control_phase)
          %fprintf('Состояние газа!!! \n')
          W=1;
          x_i=zeros(1,N);
          y_i=z;
          Z_v = Z_init;
          Z_l = 0;
    else
          %fprintf('Состояние жидкость!!! \n') 
          W=0;
          x_i=z;
          y_i=zeros(1,N);
          Z_v = 0;
          Z_l = Z_init;
    end
end

end

%% Свойства фаз !!!

%% Молекулярный вес, г/моль
molweight_v = y_i*mass';
molweight_l = x_i*mass';

%% Объем, м3/кмоль
V_v = Z_v*R*T/P;
V_l = Z_l*R*T/P;

%% Плотность, кг/м3
rho_v = molweight_v/V_v;
rho_l = molweight_l/V_l;

%% Расчет Н энтальпии газовой фазы, размерность Дж/моль
Hid = Cp1_id*(T-T0) + 0.5*Cp2_id*(T^2-T0^2) + Cp3_id*(T^3-T0^3)/3 + 0.25*Cp4_id*(T^4-T0^4);
dYi_dT = - psi_i.*(Tkr*T).^(-0.5).*(1+psi_i.*(1-(T./Tkr).^0.5));
d2Yi_dT2 = 0.5*psi_i.*(Tkr).^(-0.5).*(1+psi_i)*T^(-1.5);   % Для расчета Ср
dai_dT=ac_i.*dYi_dT;
d2ai_dT2=ac_i.*d2Yi_dT2; % Для расчета Ср

am=0;
am_l=0;
for i=1:N
    for j=1:N
        am = am + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
        am_l = am_l + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end

dam_dT=0;
dam_dT_l=0;
    for i=1:N
        for j=1:N
            dam_dT = dam_dT + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^(-0.5)*(dai_dT(i)*a_i(j) + dai_dT(j)*a_i(i));
            dam_dT_l = dam_dT_l + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^(-0.5)*(dai_dT(i)*a_i(j) + dai_dT(j)*a_i(i));
        end
    end
dam_dT = 0.5*dam_dT; 
dam_dT_l = 0.5*dam_dT_l;

dAm_dT = dam_dT*P/(R^2*T^2);
dAm_dT_l = dam_dT_l*P/(R^2*T^2);
bm=y_i*b_i';
bm_l=x_i*b_i';
cm = cpen*y_i';
cm_l = cpen*x_i';
Am=am*P/(R^2*T^2);
Aml=am_l*P/(R^2*T^2);
Bm = bm*P/(R*T);
Bml = bm_l*P/(R*T);
Cm = cm*P/(R*T);
Cml = cm_l*P/(R*T);

Hres = 1000*R*T*( Z_v - 1 - ((Am - T*dAm_dT)/Bm)*log((Z_v+Cm+Bm)/(Z_v+Cm)));
Hres_l = 1000*R*T*( Z_l - 1 - ((Aml - T*dAm_dT_l)/Bml)*log((Z_l+Cml+Bml)/(Z_l+Cml)));

H = y_i*Hid' + Hres;
Hl = x_i*Hid' + Hres_l;

%% Расчет S энтропии, размерность Дж/(моль*К)

Sid = Cp1_id*(log(T)-log(T0)) + Cp2_id*(T-T0) + 0.5*Cp3_id*(T^2-T0^2) + Cp4_id*(T^3-T0^3)/3;   %- 1000*R*log(P/0.1) - 1000*R*log(x_i);
Sres = 1000*R*(log(Z_v+Cm-Bm) + T*dAm_dT/Bm*log((Z_v+Cm+Bm)/(Z_v+Cm)));
Sres_l = 1000*R*(log(Z_l+Cml-Bml) + T*dAm_dT_l/Bml*log((Z_l+Cml+Bml)/(Z_l+Cml)));

S = y_i*Sid' - 1000*R*(y_i*log(y_i)') - 1000*R*log(P/0.101325) + Sres; 
Sl = x_i*Sid' - 1000*R*(x_i*log(x_i)') - 1000*R*log(P/0.101325) + Sres_l;

%% Расчет Cp и Cv теплоемкости газовой фазы, размерность Дж/(моль*К)

Cid = Cp1_id + Cp2_id*T + Cp3_id*T^2 + Cp4_id*T^3;

d2am_dT2=0;
d2am_dT2_l=0;
    for i=1:N
        for j=1:N
            d2am_dT2 = d2am_dT2 + y_i(i)*y_i(j)*(1-c(i,j))*( (a_i(i)/a_i(j))^(0.5)*d2ai_dT2(j) + (a_i(j)/a_i(i))^(0.5)*d2ai_dT2(i) + ...
                       0.5*(a_i(i)/a_i(j))^(0.5)*dai_dT(i)*(dai_dT(j)*a_i(i)-a_i(j)*dai_dT(i))/a_i(i)^2 + ...
                       0.5*(a_i(j)/a_i(i))^(0.5)*dai_dT(j)*(dai_dT(i)*a_i(j)-a_i(i)*dai_dT(j))/a_i(j)^2);
            d2am_dT2_l = d2am_dT2_l + x_i(i)*x_i(j)*(1-c(i,j))*( (a_i(i)/a_i(j))^(0.5)*d2ai_dT2(j) + (a_i(j)/a_i(i))^(0.5)*d2ai_dT2(i) + ...
                       0.5*(a_i(i)/a_i(j))^(0.5)*dai_dT(i)*(dai_dT(j)*a_i(i)-a_i(j)*dai_dT(i))/a_i(i)^2 + ...
                       0.5*(a_i(j)/a_i(i))^(0.5)*dai_dT(j)*(dai_dT(i)*a_i(j)-a_i(i)*dai_dT(j))/a_i(j)^2);

        end
    end
d2am_dT2 = 0.5*d2am_dT2;
d2am_dT2_l = 0.5*d2am_dT2_l;

d2Am_dT2 = d2am_dT2*P/(R^2*T^2);
d2Am_dT2_l = d2am_dT2_l*P/(R^2*T^2);

dCp1 = R*T^2*d2Am_dT2*log((Z_v+Cm+Bm)/(Z_v+Cm))/Bm;
dCp1_l = R*T^2*d2Am_dT2_l*log((Z_l+Cml+Bml)/(Z_l+Cml))/Bml;

dCp2 = ((P/T)*(1/(Z_v+Cm-Bm)-T*dAm_dT/((Z_v+Cm)*(Z_v+Cm+Bm))))^2;
dCp2_l = ((P/T)*(1/(Z_l+Cml-Bml)-T*dAm_dT_l/((Z_l+Cml)*(Z_l+Cml+Bml))))^2;

dCp3 = (P^2/(R*T))* (1/(Z_v+Cm-Bm)^2 + (Am/Bm)*(1/(Z_v+Cm+Bm)^2 - 1/(Z_v+Cm)^2));
dCp3_l = (P^2/(R*T))* (1/(Z_l+Cml-Bml)^2 + (Aml/Bml)*(1/(Z_l+Cml+Bml)^2 - 1/(Z_l+Cml)^2));

dCp = 1000*(dCp1 + T*dCp2/dCp3 - R); 
dCpl = 1000*(dCp1_l + T*dCp2_l/dCp3_l - R); 

dCv = 1000*(dCp1 - R);
dCvl = 1000*(dCp1_l - R);

Cp = y_i*Cid' + dCp;
Cpl = x_i*Cid' + dCpl;

Cv = y_i*Cid' + dCv;
Cvl = x_i*Cid' + dCvl;

%% Коэффициент Джоуля-Томсона, размерность К/МПа

JTh = -(1/Cp)*1000*((Z_v*R*T/P)-T*(dCp2)^0.5/dCp3);
JThl = -(1/Cpl)*1000*((Z_l*R*T/P)-T*(dCp2_l)^0.5/dCp3_l);

%% Скорость звука, размерность м/с
VelSound_v = 1000*V_v/(molweight_v)^0.5*(dCp3*Cp/Cv)^0.5;
VelSound_l = 1000*V_l/(molweight_l)^0.5*(dCp3_l*Cpl/Cvl)^0.5;

%% Расчёт CVD
Nmol_mix = 1;
Nmol_l = Nmol_mix*(1-W);
Vol_l = Nmol_l*molweight_l/rho_l
1000*Vol_l / 82.49
Vol_v = 82.49 - 1000*Vol_l

%% Вывод информации

L=1-W;
W
Z_v
Z_l
Mesh(step,1)=T;
Mesh(step,2)=P;
Mesh(step,3)=W;
fprintf(fid5,'%3.5f %3.5f %3.5f \t \n', Mesh(step,1), Mesh(step,2), Mesh(step,3));

fprintf(fid10,'---------------------- \n'); 
fprintf(fid10,'W =  %3.4f \n', W);
fprintf(fid10,'\n'); 
fprintf(fid10,'z Пар      z Жидкость \n'); 
fprintf(fid10,'%3.5f %3.5f \t \n', Z_v, Z_l); 
fprintf(fid10,'\n'); 
fprintf(fid10,'y Пар   x Жидкость \n'); 
for i=1:22
fprintf(fid10,'%3.5f %3.5f \t \n', 100*y_i(i), 100*x_i(i)); %, x_i); 
end
fprintf(fid10,'---------------------- \n'); 

fprintf(fid2,'%3.5f %3.9f \t \n', P, m);

fprintf(fidH,'%3.5f %5.5f %5.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f\t \n', T, H, Hl, Cp, Cpl, Cv, Cvl, S, Sl);


end   % Цикл по Р (давлению)
end   % Цикл по Т (температуре)
toc
