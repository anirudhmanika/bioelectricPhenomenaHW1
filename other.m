function other()
    % Time span
    t_span = 0:0.01:20; 
    
    % Resting potential and initial gating variables
    V_rest = -65; 
    [~, ninf, ~, ~] = n_gate(V_rest);
    n0 = ninf;
    [~, minf, ~, ~] = m_gate(V_rest);
    m0 = minf;
    [~, hinf, ~, ~] = h_gate(V_rest);
    h0 = hinf;

    % Injected current functions
    Iinj_big = @(t) (t >= 5 & t <= 5.5) * 25; 
    Iinj_small = @(t) (t >= 5 & t <= 5.5) * 6;

    % Solve ODEs
    [t_big, y_big] = ode45(@(t, y) hh(t, y, Iinj_big(t)), t_span, [n0, m0, h0, V_rest]);
    [t_small, y_small] = ode45(@(t, y) hh(t, y, Iinj_small(t)), t_span, [n0, m0, h0, V_rest]);

    % Calculate currents
    INa_big = arrayfun(@(i) calculateCurrents(y_big(i, :)), 1:length(t_big));
    IK_big = arrayfun(@(i) calculateCurrents(y_big(i, :), 'K'), 1:length(t_big));
    INa_small = arrayfun(@(i) calculateCurrents(y_small(i, :)), 1:length(t_small));
    IK_small = arrayfun(@(i) calculateCurrents(y_small(i, :), 'K'), 1:length(t_small));

    % Plot results
    figure(1);
    plot(t_big, INa_big, 'b', t_big, IK_big, 'r', ...
         t_small, INa_small, '--b', t_small, IK_small, '--r');
    xlabel('Time (ms)');
    ylabel('Current (nA)');
    title('Sodium and Potassium Currents');
    legend('I_{Na, big}', 'I_{K, big}', 'I_{Na, small}', 'I_{K, small}');

    figure(2);
    plot(t_big, y_big(:, 4), 'b', 'DisplayName', 'V_m (big pulse)');
    hold on;
    plot(t_small, y_small(:, 4), 'r', 'DisplayName', 'V_m (small pulse)');
    xlabel('Time (ms)');
    ylabel('Membrane Potential (mV)');
    title('Membrane Potential over Time');
    legend;
    hold off;

    figure(3);
    plot(t_big, y_big(:, 1), 'b', t_big, y_big(:, 2), 'g', t_big, y_big(:, 3), 'r');
    hold on;
    plot(t_small, y_small(:, 1), '--b', t_small, y_small(:, 2), '--g', t_small, y_small(:, 3), '--r');
    xlabel('Time (ms)');
    ylabel('Gating Variables');
    title('Gating Variables (n, m, h)');
    legend('n (big)', 'm (big)', 'h (big)', 'n (small)', 'm (small)', 'h (small)');
    hold off;

    % Hodgkin-Huxley Equations
    function dydt = hh(t, y, I)
        gNa = 120; gK = 36; gL = 0.3;
        ENa = 50; EK = -77; EL = -54.4;
        Cm = 1;

        n = y(1); m = y(2); h = y(3); Vm = y(4);

        [~, ~, an, bn] = n_gate(Vm);
        [~, ~, am, bm] = m_gate(Vm);
        [~, ~, ah, bh] = h_gate(Vm);

        INa = gNa * m^3 * h * (Vm - ENa);
        IK = gK * n^4 * (Vm - EK);
        IL = gL * (Vm - EL);
        Iion = INa + IK + IL;

        dydt = zeros(4, 1);
        dydt(1) = an * (1 - n) - bn * n;
        dydt(2) = am * (1 - m) - bm * m;
        dydt(3) = ah * (1 - h) - bh * h;
        dydt(4) = (I - Iion) / Cm;
    end

    % Calculate currents
    function [INa, IK] = calculateCurrents(y, type)
        gNa = 120; gK = 36;
        ENa = 50; EK = -77;
        Vm = y(4);
        m = y(2); h = y(3); n = y(1);
        switch type
            case 'Na'
                INa = gNa * m^3 * h * (Vm - ENa);
                IK = 0;
            case 'K'
                INa = 0;
                IK = gK * n^4 * (Vm - EK);
        end
    end

    % Gating Variable Functions
    function [tn, ninf, an, bn] = n_gate(Vm)
        an = 0.01 * (Vm + 55) / (1 - exp(-(Vm + 55) / 10));
        bn = 0.125 * exp(-(Vm + 65) / 80);
        ninf = an / (an + bn);
        tn = 1 / (an + bn);
    end

    function [tm, minf, am, bm] = m_gate(Vm)
        am = 0.1 * (Vm + 40) / (1 - exp(-(Vm + 40) / 10));
        bm = 4 * exp(-Vm / 18);
        minf = am / (am + bm);
        tm = 1 / (am + bm);
    end

    function [th, hinf, ah, bh] = h_gate(Vm)
        ah = 0.07 * exp(-Vm / 20);
        bh = 1 / (exp((30 - Vm) / 10) + 1);
        hinf = ah / (ah + bh);
        th = 1 / (ah + bh);
    end
end
