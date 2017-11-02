function lab03_main
%=== Дисципліна:Основи обробки біомедичної інформації ===
%--- Лабораторна робота #3 ФІЛЬТРАЦІЯ БІОСИГНАЛІВ ФІЛЬТРАМИ З НІХ
%
% Використовуйте файли даних: 
%   ecg105.txt - сигнал ЕКГ
%   ecg2x60.dat - сигнал ЕКГ з мережевою перешкодою частотою 60 Гц
%
%----------------------------------------------------------
 
clear, clc, close all
disp('Лабораторна робота #3')
disp('ФІЛЬТРАЦІЯ БІОСИГНАЛІВ ФІЛЬТРАМИ З НІХ')
disp('Виконав: Шелемба П.В., група БМ-462 ННІІДС')
 
%=== Завдання #1.1 ===
% Синтезування смугового фільтру 2-го порядку 
fs = 200;
r = 0.6;
fc = 20; % центральна частота
phi = 2*pi.*fc/fs;
b = [1 0 -1]; % чисельник CФ 
a = 1-2*r*cos(phi)*r^2; % знаменник

%=== Завдання #1.2 ===
% АЧХ та ФЧХ смугового фільтру 2-го порядку
n = 512; % кількість точок, що розраховуються
figure(1)
[h, f] = freqz(b, a, n);
mag = abs(h);
phase = angle(h)*180/pi;
subplot(2, 1, 1); plot(f/(2*pi)*fs, mag), grid on;
title('АЧХ смугового фільтру при r = 0.6'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f/(2*pi)*fs, unwrap(phase)), grid on;
title('ФЧХ смугового фільтру при r = 0.6'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
disp('Нулі смугового фільтру при r = 0.6:');
x = roots(b);
disp(x);
disp('Полюси смугового фільтру при r = 0.6:');
y = poly(a);
disp(y);

% Карта нулів та полюсів фільтру
figure(2)
z = zplane(b);
title('Карта нулів та полюсів смугового фільтру при r = 0.6');

%=== Завдання #1.3 ===
% r = 0.7
r1 = 0.7;
b1 = [1 0 -1];          
a1 = 1-2*r1*cos(phi)*r1^2;     
figure(3)
[h1, f1] = freqz(b1, a1, n);
mag2 = abs(h1);
phase2 = angle(h1)*180/pi;
subplot(2, 1, 1); plot(f1/(2*pi)*fs, mag2), grid on;
title('АЧХ смугового фільтру при r = 0.7'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f1/(2*pi)*fs, unwrap(phase2)), grid on;
title('ФЧХ смугового фільтру при r = 0.7'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
disp('Нулі смугового фільтру при r = 0.7:');
x1 = roots(b1);
disp(x1);
disp('Полюси смугового фільтру при r = 0.7:');
y1 = poly(a1);
disp(y1);

% Карта нулів та полюсів фільтру
figure(4)
z1 = zplane(b1);
title('Карта нулів та полюсів смугового фільтру при r = 0.7');

% r = 0.9
r2 = 0.9;
b2 = [1 0 -1];          
a2 = 1-2*r2*cos(phi)*r2^2;     
figure(5)
[h2, f2] = freqz(b2, a2, n);
mag3 = abs(h2);
phase3 = angle(h2)*180/pi;
subplot(2, 1, 1); plot(f2/(2*pi)*fs, mag3), grid on;
title('АЧХ смугового фільтру при r = 0.9'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f2/(2*pi)*fs, unwrap(phase3)), grid on;
title('ФЧХ смугового фільтру при r = 0.9'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
disp('Нулі смугового фільтру при r = 0.9:');
x2 = roots (b2);
disp(x2);
disp('Полюси смугового фільтру при r = 0.9:');
y2 = poly (a2);
disp(y2);

% Карта нулів та полюсів фільтру
figure(6)
z2 = zplane(b2);
title('Карта нулів та полюсів смугового фільтру при r = 0.9:');

%=== Завдання #1.5 ===
%Побудува графіків перехідних процесів
n2 = 5;
[h3, t3] = stepz(b1, a1, n2, fs);
figure(7)
stepz(t3, h3);
[h3, t3] = stepz(b2, a2, n2, fs);
figure(8)
stepz(t3, h3);

%=== Завдання #2.1 ===
% Фільтрація ЕКГ при r = 0.6 
r3 = 0.6;
b3 = [1 0 -1];          
a3 = 1-2*r3*cos(phi)*r3^2; 
ecg = load('ecg105.txt'); % сигнал ЕКГ
ecgd1 = detrend(ecg);
ecgf1 = filter(b3, a3, ecgd1);
t1 = (0:length(ecgf1)-1)/fs;
figure(9)
subplot(2, 1, 1); plot(t1, ecg), grid on;
title('Сигнал з шумом');
xlim([0 2]);
ylabel('Амплітуда');
subplot(2, 1, 2); plot(t1, ecgd1); grid on;
title('Відфільтрований сигнал при r = 0.6');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%=== Завдання #2.2 ===
% Фільтрація ЕКГ при r = 0.7 та r = 0.8
% r = 0.7
r4 = 0.7;
b4 = [1 0 -1];          
a4 = 1-2*r4*cos(phi)*r4^2; 
ecgf2 = filter(b4, a4, ecgd1);
t2 = (0:length(ecgf2)-1)/fs;
figure(10)
subplot(2, 1, 1); plot(t2, ecgd1), grid on;
title('Сигнал з шумом');
xlim([0 2]);
ylabel('Амплітуда');
subplot(2, 1, 2); plot(t2, ecgf2), grid on;
title('Відфільтрований сигнал при r = 0.7');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

% r = 0.8
r5 = 0.8;
b5 = [1 0 -1];          
a5 = 1-2*r5*cos(phi)*r5^2; 
ecgf3 = filter(b5, a5, ecgd1);
t3 = (0:length(ecgf3)-1)/fs;
figure(11)
subplot(2, 1, 1); plot(t3, ecgd1), grid on;
title('Сигнал з шумом');
xlim([0 2]);
ylabel('Амплітуда');
subplot(2, 1, 2); plot(t3, ecgf3), grid on;
title('Відфільтрований сигнал при r = 0.8');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%Графіки результатів фільтрації
figure(12)
subplot(2, 1, 1); plot(t2, ecgf2); grid on;
title('Відфільтрований сигнал при r = 0.7');
xlim([0 2]);
ylabel('Амплітуда');
subplot(2, 1, 2); plot(t3, ecgf3); grid on;
title('Відфільтрований сигнал при r = 0.8');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%=== Завдання #3.1 ===
% Дослідження властивостей режекторного НІХ-фільтру
r6 = 0.8;
phi1 = 110*pi/180;
phi2 = 130*pi/180; 
ax = [1-2*r6*cos(phi1)*r6^2] ;
ax2 = [1-2*r6*cos(phi2)*r6^2];
a6 = conv(ax, ax2);
b6 = [1 1 1]; 

% Передавальна функція
disp('Передавальна функція режекторного НІХ-фільтру:');
H1 = filt(b6, a6) % передавальна функція

%=== Завдання #3.2 ===
% АЧХ та ФЧХ смугового режекторного НІХ-фільтру
figure(12)
[h4, f4] = freqz(a6, b6, n);
mag4 = abs(h4);
phase4 = angle(h4)*180/pi;
subplot(2, 1, 1); plot(f4/(2*pi)*fs, mag4), grid on;
title('АЧХ смугового режекторного НІХ-фільтру'); ylabel ('Підсилення');
subplot(2, 1, 2); plot(f4/(2*pi)*fs, unwrap(phase4)), grid on;
title('ФЧХ смугового режекторного НІХ-фільтру'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
disp('Нулі смугового режекторного НІХ-фільтру:');
x3 = roots(b);
disp(x3);
disp('Полюси смугового режекторного НІХ-фільтру:');
y3 = poly(a);
disp(y3);
 
% Карта нулів та полюсів фільтру
figure(13)
z3 = zplane(b,a);
title('Карта нулів та полюсів смугового режекторного НІХ-фільтру:');

%=== Завдання #3.3 ===
% Порівняння АЧХ і ФЧХ режекторних НІХ і СІХ-фільтрів
figure(14)
subplot(4, 1, 1); plot(f4/(2*pi)*fs, mag4), grid on;
title('АЧХ режекторного НІХ-фільтру'); ylabel('Підсилення');
subplot(4, 1, 2); plot(f4/(2*pi)*fs, unwrap(phase4)), grid on;
title('ФЧХ режекторного НІХ-фільтру'); xlabel('Частота'); ylabel('Фаза');
b7 = [1,0.618, 1]; % коефіцієнти різницевого рівняння 
a7 = 1;
[h5, f5] = freqz(b7, a7, n, fs);
mag5 = abs(h5);
phase5 = angle(h5)*180/pi;
subplot(4, 1, 3); plot(f5/(2*pi)*fs, mag5), grid on;
title('АЧХ режекторного СІХ-фільтру'); ylabel('Підсилення');
subplot(4, 1, 4); plot(f5/(2*pi)*fs, unwrap(phase5)), grid on;
title('ФЧХ режекторного СІХ-фільтру'); xlabel('Частота'); ylabel('Фаза');

%=== Завдання #3.4 ===
% Фільтрацію сигналу ЕКГ(файл ecg2x60.dat) режекторним фільтром
ecg2 = load('ecg2x60.dat'); % сигнал ЕКГ
ecgd2 = detrend(ecg2);
ecgf4 = filter(b6, a6, ecgd2);
t4 = (0:length(ecgf4)-1)/fs;
figure(15)   
subplot(2, 1, 1); plot(t4, ecgd2), grid on;
title('Сигнал з мережевою перешкодою 60 Гц');
xlim([1 2]);
ylabel('Амплітуда');
subplot(2, 1, 2); plot (t4, ecgf4), grid on;
title('Відфільтрований сигнал');
xlim([1 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%=== Завдання #4.1 ===
% АЧХ і ФЧХ цифрових інтеграторів
% Інтегрування методом прямокутників 
fs = 300; % частота дискретизації
T = 1/fs; 
bi1 = T; % чисельник ПФ
ai1 = [1 -1]; % знаменник ПФ

% Інтегрування методом трапецій.
bi2 = [1 1]*T/2;
ai2 = [1 -1];

% Інтегрування методом парабол (Сімпсона).
bi3 = [1 4 1]*T/3;
ai3 = [1 0 -1];

[hi1, wi1] = freqz(bi1, ai1, n);
magi1 = abs(hi1);
phasei1 = angle(hi1)*180/pi;

[hi2, wi2] = freqz(bi2, ai2, n);
magi2 = abs(hi2);
phasei2 = angle(hi2)*180/pi;

[hi3, wi3] = freqz(bi3, ai3, n);
magi3 = abs(hi3);
phasei3 = angle(hi3)*180/pi;

figure(16)
subplot(3, 1, 1); plot(wi1/(2*pi)*fs, magi1), grid on;
title('АЧХ інтегратора (метод прямокутників)'); ylabel('Підсилення');
subplot(3, 1, 2); plot(wi2/(2*pi)*fs, magi2), grid on;
title ('АЧХ інтегратора (метод трапецій)'); ylabel ('Підсилення');
subplot(3, 1, 3); plot(wi3/(2*pi)*fs, magi3), grid on;
title ('АЧХ інтегратора (метод парабол)'); xlabel ('Частота'); ylabel ('Підсилення');

figure(17)
subplot(3, 1, 1); plot(wi1/(2*pi)*fs, unwrap(phasei1)), grid on;
title('ФЧХ інтегратора (метод прямокутників)'); ylabel('Фаза');
subplot(3, 1, 2); plot(wi2/(2*pi)*fs, unwrap(phasei2)), grid on;
title('ФЧХ інтегратора (метод трапецій)'); ylabel('Фаза');
subplot(3, 1, 3); plot(wi3/(2*pi)*fs, unwrap(phasei3)), grid on;
title('ФЧХ інтегратора (метод парабол)'); xlabel('Частота'); ylabel('Фаза');
