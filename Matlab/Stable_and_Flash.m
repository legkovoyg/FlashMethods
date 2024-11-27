clear all
clc
%%   Задание параметров модели
%mixture 1
%z = [0.228	0.605	83.578	7.4	3.345	0.755	0.962	0.338	0.316	0.356	0.483	0.593	0.275	0.297	0.163	0.112	0.09	0.035	0.039	0.021	0.008	0.001];
%mixture 2
%z = [0.105	0.574	52.685	6.788	3.872	1.008	1.406	0.575	0.56	0.73	1.569	2.193	1.153	1.484	1.08	1.036	1.338	0.834	1.55	2.2	2.96	14.3];
%mixture 30 h=3772.9
z = [0.386	0.871	73.444	7.677	4.522	0.944	1.656	0.554	0.665	0.754	1.372	1.596	0.818	1.066	0.814	0.677	0.641	0.345	0.464	0.366	0.259	0.109];


z=z./100;

if (sum(z)~=1.0) 
    fprintf('ОШИБКА В ЗАДАНИИ СОСТАВА СМЕСИ!!! \n') 
end
N = size(z,2);                                      % количество компонент N

Pkr = [3.394388	7.376459	4.600155	4.883864	4.245517	3.647701	3.799688	3.384255	3.374123	2.968823	3.574564	3.122868	2.771944	2.50968	2.19141	1.972617	1.792949	1.694122	1.573101	1.488414	1.424685	1.351062];
Tkr = [126.2	304.19995	190.6	305.39996	369.79996	408.1001	425.2001	460.3999	469.6	507.4001	546.6401	569.9407	591.9464	600.8699	633.8047	664.5171	699.4725	727.1604	767.4636	813.2253	865.1088	1039.5869];
w = [0.04	0.225	0.008	0.098	0.152	0.176	0.193	0.227	0.251	0.296	0.436319	0.472542	0.510256	0.616047	0.692018	0.768374	0.859284	0.931732	0.972036	1.076324	1.174966	1.109354];
c = zeros(N,N);
R = 0.0083675;                % универсальная газовая МПа*м3/(кмоль*К)

%% Задание температуры и давления
T = 273.15;    % Температура, К
P = 6;    % Давление, МПа


Ttest=200:10:700 %100:10:200;
Ptest=1:1:60;
fid2 = fopen('W(P).txt','w');
fid3 = fopen('Stable.txt','w');

step = 0;
Mesh=zeros(size(Ttest,2)*size(Ptest,2),3);

tic
for hT=1:size(Ttest,2)
    T=Ttest(hT);
    %fprintf(fid2,'Температура =  %3.2f K \n', T); 
    for hP=1:size(Ptest,2)
        P=Ptest(hP);
        step=step+1;

%% Расчетный модуль Пенга-Робинсона
ac_i = 0.457235*R^2*Tkr.^2./Pkr;                       % Коэффициенты
psi_i = 0.37464 + 1.54226*w - 0.26992*w.^2;
alpha_i = (1+psi_i.*(1-(T./Tkr).^.5)).^2;
a_i = ac_i.*alpha_i;
b_i = 0.077796*R*Tkr./Pkr;

%% Шаг 0 Начало алгоритма

%  начальное приближение К0
K_i = (exp(5.373*(1+w).*(1-Tkr/T)).*Pkr/P).^(1.0);

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

ww = roots([1 -(1-Bw) (Aw-3*Bw^2-2*Bw) -Aw*Bw+Bw^2+Bw^3]);
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

fz_i = exp(log(z*P) - log(Z_v - Bw) + b_i*(Z_v-1)/bw - (Aw/(2^(1.5)*Bw)) * ((2*avvv/aw) - (b_i/bw)) * log( (Z_v+(1+sqrt(2))*Bw) / (Z_v+(1-sqrt(2))*Bw) ));

m=0;
indx=0;
eps_f=1;
Ri_v=1;
TS_v_flag=0;
TS_l_flag=0;

%% Шаг 1 Проверка газовой фазы
%  Паровая фаза
while (m<30)
    m;

Yi_v=z.*K_i;
Sv=sum(Yi_v);
y_i=Yi_v/Sv;

aw=0;
bw=0;
for i=1:N
    for j=1:N
        aw = aw + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end

bw=y_i*b_i';
Aw = aw*P/(R^2*T^2);
Bw = bw*P/(R*T);

ww_y = roots([1 -(1-Bw) (Aw-3*Bw^2-2*Bw) -Aw*Bw+Bw^2+Bw^3]);
Z_v = ww_y(imag(ww_y)==0);
Z_v=Z_v(Z_v>0);
Z_v=max(Z_v);

%% Шаг 4 Расчет летучестей в паровой фазе

avvv=zeros(1,N);
for i=1:N
    avv = 0;
    for j=1:N
        avv = avv + y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    avvv(i)=avv;
end

fw_i = exp(log(y_i*P) - log(Z_v - Bw) + b_i*(Z_v-1)/bw - (Aw/(2^(1.5)*Bw)) * ((2*avvv/aw) - (b_i/bw)) * log( (Z_v+(1+sqrt(2))*Bw) / (Z_v+(1-sqrt(2))*Bw) ));

Ri=fz_i./(Sv*fw_i);
Ri_v=sum((Ri-1).^2);

if (Ri_v<10^(-12))
   %fprintf('Сходимость достигнута по газу!!! \n');
   Sv;
end

K_i=K_i.*Ri;
TS_v=sum(log(K_i).^2);

if (TS_v<10^(-4))
   %fprintf('TS по газу найдено!!! \n') ;
   TS_v_flag=1;
end

m=m+1;
end

K_iv=K_i;

%%% Проверка на жидкость

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

ww = roots([1 -(1-Bw) (Aw-3*Bw^2-2*Bw) -Aw*Bw+Bw^2+Bw^3]);
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

fz_i = exp(log(z*P) - log(Z_v - Bw) + b_i*(Z_v-1)/bw - (Aw/(2^(1.5)*Bw)) * ((2*avvv/aw) - (b_i/bw)) * log( (Z_v+(1+sqrt(2))*Bw) / (Z_v+(1-sqrt(2))*Bw) ));

% %% Шаг 2 Проверка жидкой фазы

ml=0;
Ri_l=1;

while (ml<30)
    ml;
    
Yi_l=z./K_i;
Sl=sum(Yi_l);
x_i=Yi_l/Sl;

al=0;
bl=0;
for i=1:N
    for j=1:N
        al = al + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end
bl=x_i*b_i';
Al = al*P/(R^2*T^2);
Bl = bl*P/(R*T);

ww_x = roots([1 -(1-Bl) (Al-3*Bl^2-2*Bl) -Al*Bl+Bl^2+Bl^3]);
Z_l = ww_x(imag(ww_x)==0);
Z_l=Z_l(Z_l>0);
Z_l=max(Z_l);
%% Шаг 6 Расчет летучестей в жидкой фазе

alll=zeros(1,N);
for i=1:N
    all = 0;
    for j=1:N
        all = all + x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    alll(i)=all;
end

fl_i = exp(log(x_i*P) - log(Z_l - Bl) + b_i*(Z_l-1)/bl - (Al/(2^(1.5)*Bl)) * ((2*alll/al) - (b_i/bl)) * log( (Z_l+(1+sqrt(2))*Bl) / (Z_l+(1-sqrt(2))*Bl) ));

Ri=Sl*fl_i./fz_i;
Ri_l=sum((Ri-1).^2);
if (Ri_l<10^(-12))
    %fprintf('Сходимость достигнута по жидкости!!! \n');
    Sl;
end
K_i=K_i.*Ri;
TS=sum(log(K_i).^2);
if (TS<10^(-4))
    %fprintf('TS по жидкости найдено!!! \n');
    TS_l_flag=1;
end

ml=ml+1;
end

K_il=K_i;

        if ((TS_l_flag==1 && TS_v_flag==1) ||  (Sv<=1 && TS_l_flag==1) || (Sl<=1 && TS_v_flag==1) || (Sv<1 && Sl<=1))
            %fprintf('Стабильное состояние!!! \n');
            Stable=1;
        else
            %fprintf('Нестабильное состояние!!! \n');
            Stable=0;
        end

if (Stable==2)

%% Расчетный модуль Пенга-Робинсона
    ac_i = 0.457235*R^2*Tkr.^2./Pkr;                       % Коэффициенты
    psi_i = 0.37464 + 1.54226*w - 0.26992*w.^2;
    alpha_i = (1+psi_i.*(1-(T./Tkr).^.5)).^2;
    a_i = ac_i.*alpha_i;
    b_i = 0.077796*R*Tkr./Pkr;

    Kst_v=sum((K_iv-1).^2);
    Kst_l=sum((K_il-1).^2);
    if (Kst_l>Kst_v)
        K_i=K_il;
    else
        K_i=K_iv;
    end
    

m=0;
eps_f=1;
while (eps_f>0.000001 & m<50)
    m;
eps_f;

%% Коэффициенты для проверки состояния фазы

ZK_Mult = sum(z.*K_i);
ZK_Div = sum(z./K_i);

%% Шаг 1 Нахождение общей доли пара

[W] = findroot(z,K_i);
W;

   if (((W<0) || (W>1)) & (m>3) & (m<5))
       W=0.5;
       fprintf('Изменил значение W!!! \n') 
   end



if (abs(W)>1000 || ((W>1) & m>5))
    fprintf('Состояние однофазное!!! \n') 
    W=1;
    break;
end

   
%% Шаг 2 Нахождение мольных долей xi, yi
x_i = z./(1 + W*(K_i-1)); 
y_i = K_i.*x_i;


%% Шаг 3 Нахождение z-фактора
%  Паровая фаза

aw=0;
bw=0;
for i=1:N
    for j=1:N
        aw = aw + y_i(i)*y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end

bw=y_i*b_i';
Aw = aw*P/(R^2*T^2);
Bw = bw*P/(R*T);

%[Z_v]=zfactor(Aw,Bw,1);
ww_y = roots([1 -(1-Bw) (Aw-3*Bw^2-2*Bw) -Aw*Bw+Bw^2+Bw^3]);
Z_v = ww_y(imag(ww_y)==0);
Z_v=Z_v(Z_v>0);
Z_v=max(Z_v);


%% Шаг 4 Расчет летучестей в паровой фазе

avvv=zeros(1,N);
for i=1:N
    avv = 0;
    for j=1:N
        avv = avv + y_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    avvv(i)=avv;
end

fw_i = exp(log(y_i*P) - log(Z_v - Bw) + b_i*(Z_v-1)/bw - (Aw/(2^(1.5)*Bw)) * ((2*avvv/aw) - (b_i/bw)) * log( (Z_v+(1+sqrt(2))*Bw) / (Z_v+(1-sqrt(2))*Bw) ));

%% Шаг 5 Нахождение z-фактора
% В жидкой фазе

al=0;
bl=0;
for i=1:N
    for j=1:N
        al = al + x_i(i)*x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
end
bl=x_i*b_i';
Al = al*P/(R^2*T^2);
Bl = bl*P/(R*T);

%[Z_l]=zfactor(Al,Bl,2);
ww_x = roots([1 -(1-Bl) (Al-3*Bl^2-2*Bl) -Al*Bl+Bl^2+Bl^3]);
Z_l = (ww_x(imag(ww_x)==0));
Z_l=Z_l(Z_l>0);
Z_l=max(Z_l);


%% Шаг 6 Расчет летучестей в жидкой фазе

alll=zeros(1,N);
for i=1:N
    all = 0;
    for j=1:N
        all = all + x_i(j)*(1-c(i,j))*(a_i(i)*a_i(j))^0.5;
    end
    alll(i)=all;
end

fl_i = exp(log(x_i*P) - log(Z_l - Bl) + b_i*(Z_l-1)/bl - (Al/(2^(1.5)*Bl)) * ((2*alll/al) - (b_i/bl)) * log( (Z_l+(1+sqrt(2))*Bl) / (Z_l+(1-sqrt(2))*Bl) ));

%% Корректировка распределения Ki

df_lv = zeros(1,N);

for t=1:N
    if (fl_i(t)~=0)
        K_i(t) = K_i(t)*(fl_i(t)/fw_i(t));
        df_lv(t) = (fl_i(t)/fw_i(t))-1;
    end
end

eps_f = max(abs(df_lv));


m=m+1;
end


%% Проверка на однофазное состояние жидкости
if (ZK_Mult<1 & ZK_Div>1)
     fprintf('Состояние жидкости!!! \n') 
     W=0;
     x_i=z;
     y_i=0;
end

%% Проверка на однофазное состояние газ
if (ZK_Mult>1 & ZK_Div<1)
    fprintf('Состояние газа!!! \n') 
    W=1;
    y_i=z;
    x_i=0;    
end


TS=sum(log(K_i).^2);
if (TS<10^(-4))
    fprintf('TS найдено!!! \n');
    TS_l_flag=1;
end

fprintf(fid2,'%3.5f %3.5f \t \n', P, W);

end

if (Stable==1)
    Kst_v=sum((K_iv-1).^2);
    Kst_l=sum((K_il-1).^2);
    if (Kst_l>Kst_v)
        K_i=K_il;
        ZK_Mult = sum(z.*K_i);
        ZK_Div = sum(z./K_i);
        %fprintf('Состояние жидкость!!! \n') 

%% Проверка на однофазное состояние жидкости
if (ZK_Mult<1 & ZK_Div>1)
     fprintf('Состояние жидкости!!! \n') 
     W=0;
     x_i=z;
     y_i=0;
end

%% Проверка на однофазное состояние газ
if (ZK_Mult>1 & ZK_Div<1)
    fprintf('Состояние газа!!! \n') 
    W=1;
    y_i=z;
    x_i=0;    
end

    else
        K_i=K_iv;
        ZK_Mult = sum(z.*K_i);
        ZK_Div = sum(z./K_i);
   %     fprintf('Состояние газ!!! \n') 

   %% Проверка на однофазное состояние жидкости
if (ZK_Mult<1 & ZK_Div>1)
     fprintf('Состояние жидкости!!! \n') 
     W=0;
     x_i=z;
     y_i=0;
end

%% Проверка на однофазное состояние газ
if (ZK_Mult>1 & ZK_Div<1)
    fprintf('Состояние газа!!! \n') 
    W=1;
    y_i=z;
    x_i=0;    
end
    end
end

Mesh(step,1)=T;
Mesh(step,2)=P;
Mesh(step,3)=Stable;
fprintf(fid3,'%3.5f %3.5f %3.5f \t \n', Mesh(step,1), Mesh(step,2), Mesh(step,3));



    end   % Цикл по Р (давлению)
end   % Цикл по Т (температуре)
toc