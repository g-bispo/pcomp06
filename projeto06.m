function projeto06_final()
% projeto06_final - Cálculo de trabalho usando malha fina de N pontos

clc; clear; close all;

%% --- Dados ---
m = 70;       % kg
g = 9.8;      % m/s^2

x = [0.000, 0.226, 0.331, 0.437, 0.539, 0.647, 0.748, 0.845, 0.954, 1.055, ...
     1.153, 1.253, 1.360, 1.462, 1.565, 1.665, 1.768, 1.867, 1.971, 2.064, ...
     2.179, 2.270, 2.371, 2.497, 2.606, 2.694, 2.801, 2.924, 3.025, 3.148, 3.242]' * 1000;

y = [832.9, 822.8, 821.0, 817.2, 815.2, 807.1, 811.7, 822.1, 834.3, 838.7, ...
     845.5, 849.2, 856.1, 859.1, 864.7, 867.0, 865.1, 865.3, 857.1, 850.6, ...
     846.8, 841.0, 833.4, 829.0, 835.4, 839.4, 846.0, 851.0, 849.0, 847.6, 842.2]';

%% --- Interpolação (Lagrange)

interp_lagrange(x,y)


%% --- Interpolação (splines) ---
hlin = interp1(x, y, 'linear', 'pp'); % spline linear em forma pp
hcub = spline(x, y);                   % spline cúbica (pp)


% Plotagem do perfil
xq_plot = linspace(min(x), max(x), 1000);
y_lin_plot = ppval(hlin, xq_plot);
y_cub_plot = ppval(hcub, xq_plot);

figure('Name','Perfil - Splines');
plot(x, y, 'ok', 'MarkerFaceColor','g'); hold on;
plot(xq_plot, y_lin_plot, '--r', 'LineWidth', 1.2);
plot(xq_plot, y_cub_plot, '-b', 'LineWidth', 1.2);
xlabel('Distância (m)'); ylabel('Altitude (m)');
legend('Pontos','Spline Linear','Spline Cúbica','Location','best'); grid on;

%% --- Malha fina ---
N_pontos_integracao = 1001;
xq_fine = linspace(x(1), x(end), N_pontos_integracao);

%% --- Derivadas (FD + Complex-step) ---
h_fd = 1e-5;   % passo
modo_fd = 'central';          % tipo
fprintf('Usando h = %.3g m para diferenças finitas (%s)\n', h_fd, modo_fd);

% Derivadas FD e Complex-step
dh_lin_fd  = derivada_fd(hlin, xq_fine, h_fd, modo_fd);
dh_cub_fd  = derivada_fd(hcub, xq_fine, h_fd, modo_fd);
dh_lin_cs  = derivada_complexa(hlin, xq_fine);
dh_cub_cs  = derivada_complexa(hcub, xq_fine);

% Comparação e escolha
rms_lin = sqrt(mean((dh_lin_fd - dh_lin_cs).^2));
rms_cub = sqrt(mean((dh_cub_fd - dh_cub_cs).^2));
fprintf('RMS(erro FD vs Complex-step): Linear = %.3e | Cúbica = %.3e\n', rms_lin, rms_cub);

dh_lin_fine = dh_lin_fd;
dh_cub_fine = dh_cub_fd;

% Plotagem das derivadas
figure('Name','Derivadas - Comparação FD x Complex-step');
subplot(2,1,1);
plot(xq_fine, dh_lin_cs, '-k');hold on;
plot(xq_fine, dh_lin_fd, '--r'); 
title('Spline Linear - Derivadas');
legend('FD','Complex-step'); grid on;

subplot(2,1,2);
plot(xq_fine, dh_cub_cs, '-k');hold on;
plot(xq_fine, dh_cub_fd, '--b'); 
title('Spline Cúbica - Derivadas');
legend('FD','Complex-step'); grid on;

%% --- Integrandos ---
sinth = @(dh) dh ./ sqrt(1 + dh.^2);

f_lin_full_fine = sinth(dh_lin_fine) .* sqrt(1 + dh_lin_fine.^2);
f_cub_full_fine = sinth(dh_cub_fine) .* sqrt(1 + dh_cub_fine.^2);

f_lin_pos_fine = max(0, f_lin_full_fine);
f_lin_neg_fine = min(0, f_lin_full_fine);
f_cub_pos_fine = max(0, f_cub_full_fine);
f_cub_neg_fine = min(0, f_cub_full_fine);

f_lin_pos_d2 = segunda_derivada(xq_fine, f_lin_pos_fine);
f_cub_pos_d2 = segunda_derivada(xq_fine, f_cub_pos_fine);

%% --- Integrações ---
ngauss = 4;

% Trapézio
[area_lin_pos_T, Et_lin_pos] = trapezio(f_lin_pos_fine, xq_fine, f_lin_pos_d2);
[area_lin_neg_T, ~] = trapezio(f_lin_neg_fine, xq_fine, f_lin_pos_d2);
[area_cub_pos_T, Et_cub_pos] = trapezio(f_cub_pos_fine, xq_fine, f_cub_pos_d2);
[area_cub_neg_T, ~] = trapezio(f_cub_neg_fine, xq_fine, f_cub_pos_d2);

% Simpson
[area_lin_pos_S, ok_lin_S] = simpson13(xq_fine, f_lin_pos_fine);
[area_lin_neg_S, ~] = simpson13(xq_fine, f_lin_neg_fine);
[area_cub_pos_S, ok_cub_S] = simpson13(xq_fine, f_cub_pos_fine);
[area_cub_neg_S, ~] = simpson13(xq_fine, f_cub_neg_fine);

% Gauss global (G1)
f_lin_pos_spline = spline(xq_fine, f_lin_pos_fine);
f_lin_neg_spline = spline(xq_fine, f_lin_neg_fine);
f_cub_pos_spline = spline(xq_fine, f_cub_pos_fine);
f_cub_neg_spline = spline(xq_fine, f_cub_neg_fine);

area_lin_pos_G1 = gauss_quad(@(xi) ppval(f_lin_pos_spline, xi), x(1), x(end), ngauss);
area_lin_neg_G1 = gauss_quad(@(xi) ppval(f_lin_neg_spline, xi), x(1), x(end), ngauss);
area_cub_pos_G1 = gauss_quad(@(xi) ppval(f_cub_pos_spline, xi), x(1), x(end), ngauss);
area_cub_neg_G1 = gauss_quad(@(xi) ppval(f_cub_neg_spline, xi), x(1), x(end), ngauss);

% Gauss por segmento (G2)
area_lin_pos_G2 = gauss_quad_spline(f_lin_pos_spline);
area_lin_neg_G2 = gauss_quad_spline(f_lin_neg_spline);
area_cub_pos_G2 = gauss_quad_spline(f_cub_pos_spline);
area_cub_neg_G2 = gauss_quad_spline(f_cub_neg_spline);

%% --- Integrações para o método Complex-Step ---
f_lin_full_cs = sinth(dh_lin_cs) .* sqrt(1 + dh_lin_cs.^2);
f_cub_full_cs = sinth(dh_cub_cs) .* sqrt(1 + dh_cub_cs.^2);

f_lin_pos_cs = max(0, f_lin_full_cs);
f_lin_neg_cs = min(0, f_lin_full_cs);
f_cub_pos_cs = max(0, f_cub_full_cs);
f_cub_neg_cs = min(0, f_cub_full_cs);

f_lin_pos_d2_c = segunda_derivada(xq_fine, f_lin_pos_cs);
f_cub_pos_d2_c = segunda_derivada(xq_fine, f_cub_pos_cs);

% --- Trapézio ---
[area_lin_pos_T_c, Et_lin_pos_c] = trapezio(f_lin_pos_cs, xq_fine, f_lin_pos_d2_c);
[area_lin_neg_T_c, ~] = trapezio(f_lin_neg_cs, xq_fine, f_lin_pos_d2_c);
[area_cub_pos_T_c, Et_cub_pos_c] = trapezio(f_cub_pos_cs, xq_fine, f_cub_pos_d2_c);
[area_cub_neg_T_c, ~] = trapezio(f_cub_neg_cs, xq_fine, f_cub_pos_d2_c);

% --- Simpson ---
[area_lin_pos_S_c, ok_lin_S_c] = simpson13(xq_fine, f_lin_pos_cs);
[area_lin_neg_S_c, ~] = simpson13(xq_fine, f_lin_neg_cs);
[area_cub_pos_S_c, ok_cub_S_c] = simpson13(xq_fine, f_cub_pos_cs);
[area_cub_neg_S_c, ~] = simpson13(xq_fine, f_cub_neg_cs);

% --- Gauss global (G1) ---
f_lin_pos_spline_c = spline(xq_fine, f_lin_pos_cs);
f_lin_neg_spline_c = spline(xq_fine, f_lin_neg_cs);
f_cub_pos_spline_c = spline(xq_fine, f_cub_pos_cs);
f_cub_neg_spline_c = spline(xq_fine, f_cub_neg_cs);

area_lin_pos_G1_c = gauss_quad(@(xi) ppval(f_lin_pos_spline_c, xi), x(1), x(end), ngauss);
area_lin_neg_G1_c = gauss_quad(@(xi) ppval(f_lin_neg_spline_c, xi), x(1), x(end), ngauss);
area_cub_pos_G1_c = gauss_quad(@(xi) ppval(f_cub_pos_spline_c, xi), x(1), x(end), ngauss);
area_cub_neg_G1_c = gauss_quad(@(xi) ppval(f_cub_neg_spline_c, xi), x(1), x(end), ngauss);

% --- Gauss por segmento (G2) ---
area_lin_pos_G2_c = gauss_quad_spline(f_lin_pos_spline_c);
area_lin_neg_G2_c = gauss_quad_spline(f_lin_neg_spline_c);
area_cub_pos_G2_c = gauss_quad_spline(f_cub_pos_spline_c);
area_cub_neg_G2_c = gauss_quad_spline(f_cub_neg_spline_c);



%% --- Trabalhos ---
W_mgh = m*g*(y(end) - y(1));

% --- Diferenças Finitas ---
% Trapézio
Wfd_lin_pos_T = m*g*area_lin_pos_T;
Wfd_lin_neg_T = m*g*area_lin_neg_T;
Wfd_lin_liq_T = Wfd_lin_pos_T + Wfd_lin_neg_T;
Wfd_cub_pos_T = m*g*area_cub_pos_T;
Wfd_cub_neg_T = m*g*area_cub_neg_T;
Wfd_cub_liq_T = Wfd_cub_pos_T + Wfd_cub_neg_T;

% Simpson
Wfd_lin_pos_S = m*g*area_lin_pos_S;
Wfd_lin_neg_S = m*g*area_lin_neg_S;
Wfd_lin_liq_S = Wfd_lin_pos_S + Wfd_lin_neg_S;
Wfd_cub_pos_S = m*g*area_cub_pos_S;
Wfd_cub_neg_S = m*g*area_cub_neg_S;
Wfd_cub_liq_S = Wfd_cub_pos_S + Wfd_cub_neg_S;

% Gauss G1
Wfd_lin_pos_G1 = m*g*area_lin_pos_G1;
Wfd_lin_neg_G1 = m*g*area_lin_neg_G1;
Wfd_lin_liq_G1 = Wfd_lin_pos_G1 + Wfd_lin_neg_G1;
Wfd_cub_pos_G1 = m*g*area_cub_pos_G1;
Wfd_cub_neg_G1 = m*g*area_cub_neg_G1;
Wfd_cub_liq_G1 = Wfd_cub_pos_G1 + Wfd_cub_neg_G1;

% Gauss G2
Wfd_lin_pos_G2 = m*g*area_lin_pos_G2;
Wfd_lin_neg_G2 = m*g*area_lin_neg_G2;
Wfd_lin_liq_G2 = Wfd_lin_pos_G2 + Wfd_lin_neg_G2;
Wfd_cub_pos_G2 = m*g*area_cub_pos_G2;
Wfd_cub_neg_G2 = m*g*area_cub_neg_G2;
Wfd_cub_liq_G2 = Wfd_cub_pos_G2 + Wfd_cub_neg_G2;


% ---Diferenças complexas ---
Wc_lin_pos_T = m*g*area_lin_pos_T_c;
Wc_lin_neg_T = m*g*area_lin_neg_T_c;
Wc_lin_liq_T = Wc_lin_pos_T + Wc_lin_neg_T;
Wc_cub_pos_T = m*g*area_cub_pos_T_c;
Wc_cub_neg_T = m*g*area_cub_neg_T_c;
Wc_cub_liq_T = Wc_cub_pos_T + Wc_cub_neg_T;

Wc_lin_pos_S = m*g*area_lin_pos_S_c;
Wc_lin_neg_S = m*g*area_lin_neg_S_c;
Wc_lin_liq_S = Wc_lin_pos_S + Wc_lin_neg_S;
Wc_cub_pos_S = m*g*area_cub_pos_S_c;
Wc_cub_neg_S = m*g*area_cub_neg_S_c;
Wc_cub_liq_S = Wc_cub_pos_S + Wc_cub_neg_S;

Wc_lin_pos_G1 = m*g*area_lin_pos_G1_c;
Wc_lin_neg_G1 = m*g*area_lin_neg_G1_c;
Wc_lin_liq_G1 = Wc_lin_pos_G1 + Wc_lin_neg_G1;
Wc_cub_pos_G1 = m*g*area_cub_pos_G1_c;
Wc_cub_neg_G1 = m*g*area_cub_neg_G1_c;
Wc_cub_liq_G1 = Wc_cub_pos_G1 + Wc_cub_neg_G1;

Wc_lin_pos_G2 = m*g*area_lin_pos_G2_c;
Wc_lin_neg_G2 = m*g*area_lin_neg_G2_c;
Wc_lin_liq_G2 = Wc_lin_pos_G2 + Wc_lin_neg_G2;
Wc_cub_pos_G2 = m*g*area_cub_pos_G2_c;
Wc_cub_neg_G2 = m*g*area_cub_neg_G2_c;
Wc_cub_liq_G2 = Wc_cub_pos_G2 + Wc_cub_neg_G2;


%% --- Impressão dos resultados ---
fprintf('\n=== RESULTADOS: TRABALHO (J) [N=%d] ===\n', N_pontos_integracao);
fprintf('Analítico (m*g*Δh): %.6f J\n\n', W_mgh);

fprintf('==================== Diferenças finitas ====================\n');
fprintf('--- TRAPÉZIO ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f | Erro_T = %.3e\n\n', Wfd_lin_pos_T, Wfd_lin_neg_T, Wfd_lin_liq_T,Et_lin_pos);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f | Erro_T = %.3e\n\n', Wfd_cub_pos_T, Wfd_cub_neg_T, Wfd_cub_liq_T, Et_cub_pos);

fprintf('--- SIMPSON 1/3 ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wfd_lin_pos_S, Wfd_lin_neg_S, Wfd_lin_liq_S);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wfd_cub_pos_S, Wfd_cub_neg_S, Wfd_cub_liq_S);

fprintf('--- GAUSS-LEGENDRE (G1: %d pts global) ---\n', ngauss);
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wfd_lin_pos_G1, Wfd_lin_neg_G1, Wfd_lin_liq_G1);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wfd_cub_pos_G1, Wfd_cub_neg_G1, Wfd_cub_liq_G1);

fprintf('--- GAUSS-LEGENDRE (G2: 2-pt/segmento) ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wfd_lin_pos_G2, Wfd_lin_neg_G2, Wfd_lin_liq_G2);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wfd_cub_pos_G2, Wfd_cub_neg_G2, Wfd_cub_liq_G2);


fprintf('==================== Diferenças complexas ====================\n');
fprintf('--- TRAPÉZIO ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f | Erro_T = %.3e\n\n', Wc_lin_pos_T, Wc_lin_neg_T, Wc_lin_liq_T,Et_lin_pos_c);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f | Erro_T = %.3e\n\n', Wc_cub_pos_T, Wc_cub_neg_T, Wc_cub_liq_T, Et_cub_pos_c);

fprintf('--- SIMPSON 1/3 ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wc_lin_pos_S, Wc_lin_neg_S, Wc_lin_liq_S);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wc_cub_pos_S, Wc_cub_neg_S, Wc_cub_liq_S);

fprintf('--- GAUSS-LEGENDRE (G1: %d pts global) ---\n', ngauss);
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wc_lin_pos_G1, Wc_lin_neg_G1, Wc_lin_liq_G1);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wc_cub_pos_G1, Wc_cub_neg_G1, Wc_cub_liq_G1);

fprintf('--- GAUSS-LEGENDRE (G2: 2-pt/segmento) ---\n');
fprintf('Linear:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n', Wc_lin_pos_G2, Wc_lin_neg_G2, Wc_lin_liq_G2);
fprintf('Cúbica:  W+ = %10.4f | W- = %10.4f | Wliq = %10.4f\n\n', Wc_cub_pos_G2, Wc_cub_neg_G2, Wc_cub_liq_G2);

%% --- Erros relativos (%) ---
fprintf('==================== Erros Relativos ====================\n');
err = @(W) 100*abs((W - W_mgh)/W_mgh);
fprintf('Linear (FD): Trap = %.4f%% | Simp = %.4f%% | G1 = %.4f%% | G2 = %.4f%%\n', ...
    err(Wfd_lin_liq_T), err(Wfd_lin_liq_S), err(Wfd_lin_liq_G1), err(Wfd_lin_liq_G2));
fprintf('Cúbica (FD): Trap = %.4f%% | Simp = %.4f%% | G1 = %.4f%% | G2 = %.4f%%\n', ...
    err(Wfd_cub_liq_T), err(Wfd_cub_liq_S), err(Wfd_cub_liq_G1), err(Wfd_cub_liq_G2));
fprintf('Linear (CS): Trap = %.4f%% | Simp = %.4f%% | G1 = %.4f%% | G2 = %.4f%%\n', ...
    err(Wc_lin_liq_T), err(Wc_lin_liq_S), err(Wc_lin_liq_G1), err(Wc_lin_liq_G2));
fprintf('Cúbica (CS): Trap = %.4f%% | Simp = %.4f%% | G1 = %.4f%% | G2 = %.4f%%\n', ...
    err(Wc_cub_liq_T), err(Wc_cub_liq_S), err(Wc_cub_liq_G1), err(Wc_cub_liq_G2));


%% --- Gráficos Adicionais (Análise - Itens 2d e 5d) ---

% 1. Plot da Derivada h'(x) (Perfil de Subidas/Descidas)
figure('Name', 'Perfil - Derivada h''(x) (Malha Fina)');
plot(xq_fine, dh_lin_fine, '--r', 'LineWidth', 1.2); hold on;
plot(xq_fine, dh_cub_fine, '-b', 'LineWidth', 1.2);
plot(x, zeros(size(x)), 'k--'); % Linha zero
xlabel('Distância (m)'); ylabel('Inclinação (m/m) - h''(x)');
title('Perfil de Inclinação (h''(x) vs. x)');
legend('h''(x) Linear','h''(x) Cúbica', 'Nível (h''=0)', 'Location','best');
grid on;
% Discussão (Item 2d): A derivada da spline linear é constante por trechos.
% A da spline cúbica é suave e mais realista.

% 2. Plot do Integrando f(x) (Contribuição ao Trabalho W+)

figure('Name', 'Integrando - Trabalho nas Subidas (Malha Fina)');
area(xq_fine, f_lin_pos_fine, 'FaceColor', [0.5 0.9 0.5], 'EdgeColor', 'none');
hold on;
plot(xq_fine, f_lin_pos_fine, '-b', 'LineWidth', 1.2);
xlabel('Distância (m)'); ylabel('Integrando f(x)');
title('Contribuição ao Trabalho W+ (f(x) - Spline Linear)');
legend('Área = f(x) > 0', 'f(x) (Cúbica)', 'Location','best');
grid on;


figure('Name', 'Integrando - Trabalho nas Subidas (Malha Fina)');
area(xq_fine, f_cub_pos_fine, 'FaceColor', [0.5 0.9 0.5], 'EdgeColor', 'none');
hold on;
plot(xq_fine, f_cub_pos_fine, '-b', 'LineWidth', 1.2);
xlabel('Distância (m)'); ylabel('Integrando f(x)');
title('Contribuição ao Trabalho W+ (f(x) - Spline Cúbica)');
legend('Área = f(x) > 0', 'f(x) (Cúbica)', 'Location','best');
grid on;





end


%% ------------------------------------------------------------------------
%% --- Funções auxiliares
%% ------------------------------------------------------------------------

function d = derivada_fd(spline_fun, x, h, tipo)

n = length(x);
d = zeros(size(x));

% função genérica para avaliar spline
try
    f = @(xx) fnval(spline_fun, xx);
catch
    f = @(xx) ppval(spline_fun, xx);
end

switch lower(tipo)
    case 'forward'
        for i = 1:n
            d(i) = (f(x(i) + h) - f(x(i))) / h;
        end
    case 'backward'
        for i = 1:n
            d(i) = (f(x(i)) - f(x(i) - h)) / h;
        end
    otherwise % 'central'
        for i = 1:n
            d(i) = (f(x(i) + h/2) - f(x(i) - h/2)) / h;
        end
end
end


function d = derivada_complexa(spline_fun, x)
% DERIVADA_COMPLEXA - Derivada via complex-step
h = 1e-20;
try
    vals = fnval(spline_fun, x + 1i*h);
catch
    vals = ppval(spline_fun, x + 1i*h);
end
d = imag(vals) / h;
end


function d2 = segunda_derivada(x, y)
% SEGUNDA_DERIVADA - 2ª derivada numérica (diferenças centrais)
n = length(x);
d2 = zeros(size(y));
for i = 2:n-1
    h1 = x(i) - x(i-1);
    h2 = x(i+1) - x(i);
    y1 = (y(i+1) - y(i)) / h2;
    y0 = (y(i) - y(i-1)) / h1;
    d2(i) = 2 * (y1 - y0) / (h1 + h2);
end
d2(1) = d2(2);
d2(n) = d2(n-1);
end


function [I, Et] = trapezio(f, x, f2)
% TRAPÉZIO COM ESTIMATIVA DE ERRO
I = trapz(x, f);
h = x(2) - x(1);
M2 = max(abs(f2));
Et = abs(((x(end) - x(1))/12) * h^2 * M2);
end


function [I, ok] = simpson13(x, f)
% SIMPSON 1/3 COMPOSTO (para malha uniforme)
N = length(x);
if mod(N,2) == 0
    warning('Simpson requer número ímpar de pontos.');
    ok = false; I = NaN; return;
end
ok = true;
h = x(2) - x(1);
I = (h/3)*(f(1) + 4*sum(f(2:2:end-1)) + 2*sum(f(3:2:end-2)) + f(end));
end


function I = gauss_quad(f, a, b, n)
% QUADRATURA DE GAUSS GLOBAL (N pontos)
[xi, wi] = gauss_legendre_nodes_weights(n);
s = 0;
for k = 1:n
    xk = ((b - a)/2)*xi(k) + (a + b)/2;
    s = s + wi(k) * f(xk);
end
I = s * (b - a)/2;
end


function [x, w] = gauss_legendre_nodes_weights(n)
switch n
    case 2
        x = [-1/sqrt(3), 1/sqrt(3)];
        w = [1, 1];
    case 3
        x = [-sqrt(3/5), 0, sqrt(3/5)];
        w = [5/9, 8/9, 5/9];
    case 4
        x = [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7), ...
              sqrt((3-2*sqrt(6/5))/7),  sqrt((3+2*sqrt(6/5))/7)];
        w = [(18 - sqrt(30))/36, (18 + sqrt(30))/36, ...
             (18 + sqrt(30))/36, (18 - sqrt(30))/36];
    otherwise
        error('Número de pontos de Gauss não suportado.');
end
end


function I = gauss_quad_spline(pp)
% GAUSS_QUAD_SPLINE - 2-pt exato por segmento para spline cúbica
breaks = pp.breaks;
coefs  = pp.coefs;
nPieces = size(coefs, 1);
xi = [-1/sqrt(3), 1/sqrt(3)];
wi = [1, 1];
I = 0;
for k = 1:nPieces
    a = breaks(k); b = breaks(k+1);
    h = (b - a)/2; mid = (a + b)/2;
    c = coefs(k,:);
    for j = 1:2
        xg = h*xi(j) + mid;
        dx = xg - a;
        f = ((c(1)*dx + c(2))*dx + c(3))*dx + c(4);
        I = I + wi(j)*f*h;
    end
end
end

function interp_lagrange(x, f)
tic
    % entrada
    % x: nós da interpolação
    % f: valores de f(x)
    
    nnos = length(x); % número de nós
    
    % verificar se há valores repetidos em x
    if nnos ~= length(unique(x))
        error('x tem entradas repetidas')
    end
    
    % coeficientes do polinômio interpolante (são os próprios f)
    c = f;
    
    toc
    
    % domínio de avaliação dentro do intervalo dos dados
    refino = 100;
    xg = linspace(x(1), x(end), refino);
    fg = zeros(1, refino);
    
    % cálculo do polinômio interpolante de Lagrange
    for i = 1:refino
        for j = 1:nnos
            Lj = 1;
            for k = 1:nnos
                if k ~= j
                    Lj = Lj * (xg(i) - x(k)) / (x(j) - x(k));
                end
            end
            fg(i) = fg(i) + c(j) * Lj;
        end
    end
    
    % gráfico
    figure;
    plot(x, f, '*r', 'MarkerSize', 8); % pontos originais em preto
    hold on
    plot(xg, fg, 'b', 'LineWidth', 1.5); % interpolação
    grid on
    xlabel('x')
    ylabel('h(x)')
    legend('Pontos originais', 'Interpolação de Lagrange', 'Location', 'best')
    grid on
    hold off
end

