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

%% ---------------- Exportacao de tabelas para comparacao com o PSIM -----
% Gera um arquivo .txt por caso, com os valores calculados (analitico)
% EXATAMENTE nos mesmos instantes de tempo do MEDICOES.txt do PSIM
% (1 ms a 100 ms, passo de 1 ms) -> facil de colocar lado a lado.
%
% Se o arquivo MEDICOES.txt do PSIM correspondente estiver na mesma pasta
% do MATLAB Online (renomeie/ajuste o nome no vetor 'arquivos_psim' abaixo),
% o script tambem calcula o erro absoluto e o erro percentual automaticamente.

t_psim = (0.001:0.001:0.1)';   % mesma janela e passo do MEDICOES.txt

% Nomes dos arquivos MEDICOES.txt do PSIM (ajuste se necessario).
% Deixe '' (vazio) se voce nao for subir o arquivo do PSIM no MATLAB Online.
arquivos_psim = {'MEDICOES_subamortecido.txt', ...
                 'MEDICOES_criticamente_amortecido.txt', ...
                 'MEDICOES_superamortecido.txt'};

% Coluna (dentro do MEDICOES.txt) onde estao IT e VC em cada pasta do
% repositorio original (a ordem das colunas NAO e igual em todas!):
%   SUBAMORTECIDO / SUPERAMORTECIDO : Time, VI, IT, VC, VR  -> IT=col3, VC=col4
%   CRITICAMENTE_AMORTECIDO         : Time, IT, VI, VR, VC  -> IT=col2, VC=col5
col_IT = [3, 2, 3];   % k = 1 (Sub), 2 (Critico), 3 (Super)
col_VC = [4, 5, 4];

for k = 1:3
    R = results(k).R;
    alpha = R/(2*L);
    wn    = 1/sqrt(L*C);
    disc  = alpha^2 - wn^2;

    if disc > 1e-6
        s1 = -alpha + sqrt(disc); s2 = -alpha - sqrt(disc);
        i_calc  = (Vs/L) * (exp(s2*t_psim) - exp(s1*t_psim)) / (s2 - s1);
        vC_calc = Vs*( 1 + (s1.*exp(s2*t_psim) - s2.*exp(s1*t_psim))/(s2 - s1) );
    elseif abs(disc) <= 1e-6
        i_calc  = (Vs/L) * t_psim .* exp(-alpha*t_psim);
        vC_calc = Vs*( 1 - (1 + alpha*t_psim).*exp(-alpha*t_psim) );
    else
        wd = sqrt(-disc);
        i_calc  = (Vs/(wd*L)) * exp(-alpha*t_psim) .* sin(wd*t_psim);
        vC_calc = Vs*( 1 - exp(-alpha*t_psim).*( cos(wd*t_psim) + (alpha/wd)*sin(wd*t_psim) ) );
    end
    vR_calc = R * i_calc;
    vL_calc = Vs - vR_calc - vC_calc;

    fname_out = sprintf('calculado_%s.txt', strrep(lower(results(k).nome),' ','_'));
    fid = fopen(fname_out, 'w');

    % Tenta carregar o MEDICOES.txt do PSIM correspondente (opcional)
    tem_psim = false;
    if ~isempty(arquivos_psim{k}) && isfile(arquivos_psim{k})
        try
            psim = readmatrix(arquivos_psim{k}, 'NumHeaderLines', 1);
            IT_psim = psim(:, col_IT(k));
            VC_psim = psim(:, col_VC(k));
            tem_psim = true;
        catch
            tem_psim = false;
        end
    end

    if tem_psim
        erro_IT = i_calc - IT_psim;
        erro_VC = vC_calc - VC_psim;
        fprintf(fid, '%-14s %-14s %-14s %-14s %-14s %-14s %-14s\n', ...
            'Time(s)','IT_calc(A)','IT_psim(A)','Erro_IT(A)', ...
            'VC_calc(V)','VC_psim(V)','Erro_VC(V)');
        for n = 1:length(t_psim)
            fprintf(fid, '%-14.6e %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', ...
                t_psim(n), i_calc(n), IT_psim(n), erro_IT(n), ...
                vC_calc(n), VC_psim(n), erro_VC(n));
        end
        fprintf('  -> %s comparado ao PSIM: erro max IT = %.4g A | erro max VC = %.4g V\n', ...
            results(k).nome, max(abs(erro_IT)), max(abs(erro_VC)));
    else
        fprintf(fid, '%-14s %-14s %-14s %-14s %-14s\n', ...
            'Time(s)','IT_calc(A)','VC_calc(V)','VR_calc(V)','VL_calc(V)');
        for n = 1:length(t_psim)
            fprintf(fid, '%-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', ...
                t_psim(n), i_calc(n), vC_calc(n), vR_calc(n), vL_calc(n));
        end
    end
    fclose(fid);
    fprintf('Tabela exportada: %s\n', fname_out);
end

fprintf('\nConcluido. PDFs vetoriais exportados no diretorio atual.\n');
