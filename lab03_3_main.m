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
disp('Виконав: Шелемба П.В., група БМ-462 ІІДС')
 
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
[h,f] = freqz(b, a, n);
mag = abs(h);
phase = angle(h)*180/pi;
subplot(2, 1, 1); plot(f/(2*pi)*fs, mag), grid on;
title('АЧХ смугового фільтру 2-го порядку'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f/(2*pi)*fs, unwrap(phase)), grid on;
title('ФЧХ смугового фільтру 2-го порядку'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
x = roots (b)
y = poly (a)

% Карта нулів та полюсів фільтру
figure(2)
z = zplane(b);

%=== Завдання #1.3 ===
% r = 0.7
r1 = 0.7;
b1 = [1 0 -1];          
a1 = 1-2*r1*cos(phi)*r1^2;     
figure(3)
[h,f] = freqz(b1, a1, n);
mag2 = abs(h);
phase = angle(h)*180/pi;
subplot(2, 1, 1); plot(f/(2*pi)*fs, mag2), grid on;
title('АЧХ фільтру'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f/(2*pi)*fs, unwrap(phase)), grid on;
title('ФЧХ фільтру'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
x1 = roots (b1)
y1 = poly (a1)

% Карта нулів та полюсів фільтру
figure(4)
z1 = zplane(b1);

% r = 0.9
r2 = 0.9;
b2 = [1 0 -1];          
a2 = 1-2*r2*cos(phi)*r2^2;     
figure(5)
[h,f] = freqz(b2, a2, n);
mag3 = abs(h);
phase = angle(h)*180/pi;
subplot(2, 1, 1); plot(f/(2*pi)*fs, mag3), grid on;
title('АЧХ фільтру'); ylabel('Підсилення');
subplot(2, 1, 2); plot(f/(2*pi)*fs, unwrap(phase)), grid on;
title('ФЧХ фільтру'); xlabel('Частота'); ylabel('Фаза');

% Обчислення нулів та полюсів фільтру
x2 = roots (b2)
y2 = poly (a2)

% Карта нулів та полюсів фільтру
figure(6)
z2 = zplane(b2);

%=== Завдання #1.4 ===
%Побудува графіків перехідних процесів
n = 5;
[h,t] = stepz(b1, a1, n, fs)
figure(7)
plot (t,h);
[h,t] = stepz(b2, a2, n, fs)
figure(8)
plot (t,h);

%=== Завдання #2.1 ===
% Фільтрація ЕКГ при r = 0.6
fs = 200; 
r3 = 0.6;
b3 = [1 0 -1];          
a3 = 1-2*r3*cos(phi)*r3^2; 
ecg = load('ecg105.txt'); % сигнал ЕКГ
ecg1 = detrend(ecg);
ecgf1 = filter(b3, a3, ecg1);
t1 = (0:length(ecgf1)-1)/fs;
figure(9)
[h,f] = freqz(b3, a3, n);
mag3 = abs(h);
phase = angle(h)*180/pi;
subplot (2, 1, 1); plot(t1, ecg1), grid on;
title('Нефільтрований сигнал');
xlim([0 2]);
ylabel('Амплітуда');
subplot (2, 1, 2); plot (t1, ecgf1); grid on;
title('Відфільтрований сигнал');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%=== Завдання #2.2 ===
% Фільтрація ЕКГ при r = 0.7 та r = 0.8
% r = 0.7
fs = 200; 
r4 = 0.7;
b4 = [1 0 -1];          
a4 = 1-2*r4*cos(phi)*r4^2; 
ecg = load('ecg105.txt'); % сигнал ЕКГ
ecg2 = detrend(ecg);
ecgf2 = filter(b4, a4, ecg2);
t2 = (0:length(ecgf2)-1)/fs;
figure(10)
[h,f] = freqz(b4, a4, n);
mag4 = abs(h);
phase = angle(h)*180/pi;
subplot (2, 1, 1); plot(t2, ecg2), grid on;
title('Нефільтрований сигнал');
xlim([0 2]);
ylabel('Амплітуда');
subplot (2, 1, 2); plot(t2, ecgf2); grid on;
title('Відфільтрований сигнал');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

% r = 0.8
fs = 200; 
r5 = 0.8;
b3 = [1 0 -1];          
a3 = 1-2*r5*cos(phi)*r5^2; 
ecg = load('ecg105.txt'); % сигнал ЕКГ
ecg3 = detrend(ecg);
ecgf3 = filter(b3, a3, ecg3);
t3 = (0:length(ecgf3)-1)/fs;
figure(11)
[h,f] = freqz(b3, a3, n);
mag3 = abs(h);
phase = angle(h)*180/pi;
subplot (2, 1, 1); plot(t3, ecg3), grid on;
title('Нефільтрований сигнал');
xlim([0 2]);
ylabel('Амплітуда');
subplot (2, 1, 2); plot(t3, ecgf3); grid on;
title('Відфільтрований сигнал');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');

%Графіки результатів фільтрації
figure(12)
subplot(2, 1, 1); plot(t2, ecgf2); grid on;
title('Відфільтрований сигнал при r = 0.7');
xlim([0 2]);
ylabel('Амплітуда');
subplot (2, 1, 2); plot (t3, ecgf3); grid on;
title('Відфільтрований сигнал при r = 0.8');
xlim([0 2]);
xlabel('Відліки'); ylabel('Амплітуда');
