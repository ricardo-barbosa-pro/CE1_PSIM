%% =======================================================================
%  RLC Series Circuit - Step Response Analysis
%  Superamortecido / Criticamente Amortecido / Subamortecido
%
%  Circuito: fonte Vs em degrau (t=0), R, L, C em série
%  Condicoes iniciais: i(0) = 0 A ,  vC(0) = 0 V
%
%  Equacao caracteristica:  s^2 + (R/L) s + 1/(LC) = 0
%  Raizes:  s1,2 = -R/(2L) +- sqrt( (R/(2L))^2 - 1/(LC) )
% =======================================================================
clear; clc; close all;

% Forca tema claro (fundo branco, texto preto) em todas as figuras,
% mesmo que o MATLAB Online esteja com tema escuro ativado.
% (aplicado diretamente em cada eixo/legenda logo apos serem criados)

%% ---------------- Parametros do circuito (fixos) -----------------------
L  = 100e-3;    % Indutor  = 100 mH
C  = 100e-6;    % Capacitor = 100 uF
Vs = 100;       % Degrau de tensao da fonte = 100 V

% Resistencia critica teorica: Rc = 2*sqrt(L/C)
Rc_teorico = 2*sqrt(L/C);
fprintf('Resistencia critica teorica: Rc = %.4f ohm\n', Rc_teorico);

% Os tres casos simulados no PSIM
R_casos   = [20, 62, 100];                 % Ohms
nomes     = {'Subamortecido','Criticamente Amortecido','Superamortecido'};
cores     = {[0.85 0.10 0.10],[0.10 0.55 0.10],[0.10 0.30 0.85]};

% Vetor de tempo (mesma janela das medicoes do PSIM: 0 a 100 ms)
t = linspace(0, 0.1, 5000);   % 5000 pontos, resolucao fina p/ curvas suaves

%% ---------------- Loop pelos tres casos ---------------------------------
results = struct();

for k = 1:3
    R = R_casos(k);

    alpha = R/(2*L);              % coeficiente de amortecimento
    wn    = 1/sqrt(L*C);          % frequencia natural nao amortecida
    disc  = alpha^2 - wn^2;       % discriminante (em termos de alpha, wn)

    fprintf('\n--- Caso: %s (R = %d ohm) ---\n', nomes{k}, R);
    fprintf('alpha = %.4f , wn = %.4f , alpha^2 - wn^2 = %.4f\n', alpha, wn, disc);

    % ------------- Solucao analitica (resposta forcada ao degrau) -------
    if disc > 1e-6
        % ---------- SUPERAMORTECIDO: raizes reais distintas ----------
        s1 = -alpha + sqrt(disc);
        s2 = -alpha - sqrt(disc);

        vC = Vs*( 1 + (s1.*exp(s2*t) - s2.*exp(s1*t))/(s2 - s1) );
        i  = (Vs/L) * (exp(s2*t) - exp(s1*t)) / (s2 - s1);

        tipo = 'Superamortecido (raizes reais distintas)';

    elseif abs(disc) <= 1e-6
        % ---------- CRITICAMENTE AMORTECIDO: raiz dupla ----------
        s = -alpha;

        vC = Vs*( 1 - (1 + alpha*t).*exp(-alpha*t) );
        i  = (Vs/L) * t .* exp(-alpha*t);

        tipo = 'Criticamente amortecido (raiz dupla)';

    else
        % ---------- SUBAMORTECIDO: raizes complexas conjugadas ----------
        wd = sqrt(-disc);   % frequencia natural amortecida

        vC = Vs*( 1 - exp(-alpha*t).*( cos(wd*t) + (alpha/wd)*sin(wd*t) ) );
        i  = (Vs/(wd*L)) * exp(-alpha*t) .* sin(wd*t);

        tipo = 'Subamortecido (raizes complexas conjugadas)';
    end

    vR = R * i;                 % Lei de Ohm
    vL = Vs - vR - vC;          % KVL: Vs = vR + vL + vC

    fprintf('Tipo identificado: %s\n', tipo);

    % Armazena para plot e para comparacao
    results(k).nome = nomes{k};
    results(k).R    = R;
    results(k).t    = t;
    results(k).i    = i;
    results(k).vC   = vC;
    results(k).vR   = vR;
    results(k).vL   = vL;

    % ------------- Solucao numerica via ODE45 (equacoes de estado) -------
    % Variaveis de estado: x1 = vC ,  x2 = i
    % dx1/dt = i/C
    % dx2/dt = (Vs - R*i - vC)/L
    odefun = @(t,x) [ x(2)/C ;
                      (Vs - R*x(2) - x(1))/L ];
    x0 = [0; 0];   % vC(0)=0, i(0)=0
    [t_ode, x_ode] = ode45(odefun, t, x0);

    results(k).t_ode  = t_ode;
    results(k).vC_ode = x_ode(:,1);
    results(k).i_ode  = x_ode(:,2);
end

%% ---------------- Figura 1: Corrente i(t) - analitico x ODE45 ----------
fig1 = figure('Color','w');
for k = 1:3
    ax = subplot(3,1,k);
    plot(results(k).t*1e3, results(k).i, '-', 'Color', cores{k}, 'LineWidth', 1.6); hold on;
    plot(results(k).t_ode*1e3, results(k).i_ode, 'k--', 'LineWidth', 1.0);
    grid on; box on;
    xlabel('Tempo (ms)'); ylabel('i(t)  [A]');
    title(sprintf('%s  (R = %d \\Omega)', results(k).nome, results(k).R));
    lg = legend('Analitico', 'ODE45', 'Location', 'best');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    set(get(ax,'Title'),'Color','k'); set(get(ax,'XLabel'),'Color','k'); set(get(ax,'YLabel'),'Color','k');
    set(lg, 'TextColor', 'k', 'Color', 'w');
end
sgtitle('Corrente no circuito RLC serie - solucao analitica vs numerica (ODE45)', 'Color', 'k');

% Exporta em formato vetorial (PDF)
exportgraphics(fig1, 'corrente_i_t.pdf', 'ContentType', 'vector');

%% ---------------- Figura 2: Tensao no capacitor vC(t) ------------------
fig2 = figure('Color','w');
ax2 = gca;
for k = 1:3
    plot(results(k).t*1e3, results(k).vC, 'LineWidth', 1.8, 'Color', cores{k}); hold on;
end
yline(Vs, 'k:', 'V_s = 100 V');
grid on; box on;
xlabel('Tempo (ms)'); ylabel('v_C(t)  [V]');
title('Tensao no capacitor - comparacao dos tres regimes de amortecimento');
lg2 = legend(nomes{:}, 'Location', 'southeast');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
set(get(ax2,'Title'),'Color','k'); set(get(ax2,'XLabel'),'Color','k'); set(get(ax2,'YLabel'),'Color','k');
set(lg2, 'TextColor', 'k', 'Color', 'w');

exportgraphics(fig2, 'tensao_vC_t.pdf', 'ContentType', 'vector');

%% ---------------- Figura 3: vR, vL, vC juntos (por caso) ---------------
for k = 1:3
    figk = figure('Color','w');
    axk = gca;
    plot(results(k).t*1e3, results(k).vR, 'LineWidth', 1.6); hold on;
    plot(results(k).t*1e3, results(k).vL, 'LineWidth', 1.6);
    plot(results(k).t*1e3, results(k).vC, 'LineWidth', 1.6);
    grid on; box on;
    xlabel('Tempo (ms)'); ylabel('Tensao (V)');
    title(sprintf('%s (R = %d \\Omega) - v_R, v_L, v_C', results(k).nome, results(k).R));
    lgk = legend('v_R(t)','v_L(t)','v_C(t)', 'Location', 'best');
    set(axk, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    set(get(axk,'Title'),'Color','k'); set(get(axk,'XLabel'),'Color','k'); set(get(axk,'YLabel'),'Color','k');
    set(lgk, 'TextColor', 'k', 'Color', 'w');

    fname = sprintf('formas_onda_%s.pdf', strrep(lower(results(k).nome),' ','_'));
    exportgraphics(figk, fname, 'ContentType', 'vector');
end

%% ---------------- (Opcional) Comparacao com dados do PSIM --------------
% Se voce exportar o MEDICOES.txt de cada pasta do PSIM para o mesmo
% diretorio do MATLAB Online, descomente o bloco abaixo para sobrepor
% os pontos medidos no PSIM com a curva analitica em MATLAB.
%
% psim = readmatrix('SUPERAMORTECIDO_MEDICOES.txt'); % ajuste o nome/colunas
% figure;
% plot(results(3).t*1e3, results(3).i, 'b-', 'LineWidth', 1.5); hold on;
% plot(psim(:,1)*1e3, psim(:,3), 'ro', 'MarkerSize', 4); % coluna IT do PSIM
% legend('MATLAB (analitico)', 'PSIM'); grid on;
% xlabel('Tempo (ms)'); ylabel('i(t) [A]');
% title('Comparacao MATLAB x PSIM - caso superamortecido');

fprintf('\nConcluido. PDFs vetoriais exportados no diretorio atual.\n');
