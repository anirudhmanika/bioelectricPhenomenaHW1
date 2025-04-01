% Solve the Hodgkin-Huxley equations using the modified current step
[t,Vm,n,m,h,INa,IK] = solve_hh(100, @currentStep);

% Generate the current vector for plotting
current_inj = arrayfun(@currentStep, t);

% Plot membrane potential (Vm) and injected current (inj) as subplots
figure;

% Subplot 1: Membrane Potential (Vm)
subplot(2,1,1);  % 2 rows, 1 column, 1st plot
plot(t, Vm, 'b', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Over Time with Hyperpolarizing Stimulus');
grid on;

% Subplot 2: Injected Current (inj)
subplot(2,1,2);  % 2 rows, 1 column, 2nd plot
plot(t, current_inj, 'r', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Injected Current (nA)');
title('Injected Hyperpolarizing Current');
grid on;



%% Function for hyperpolarizing stimulus
function inj = currentStep(t)
    % Hyperpolarizing square pulse
    if t >= 11 && t < (11 + 5)  % You can adjust the 5 to increase the stimulus duration
        inj = -20;  % Hyperpolarizing current (negative stimulus)
    else
        inj = 0;
    end
end

%% Hodgkin-Huxley solver
function [t,Vf,nf,mf,hf,INaf,IKf] = solve_hh(dt, currentStep)
    % Constants
    gNa_max = 120; gK_max = 36; gLeak = 0.3; 
    ELeak = -50; ENa = 57; EK = -75;
    Cm = 1; Vrest = -61.855;

    % Initial n, m, h
    [an0,bn0,am0,bm0,ah0,bh0] = ab(0);
    [~,n0] = ti(an0,bn0); [~,m0] = ti(am0,bm0); [~,h0] = ti(ah0,bh0);
    y0 = [Vrest, n0, m0, h0];

    % ODE solver options
    options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',0.01);
    [t,y] = ode45(@(t,y) hh(t,y,currentStep(t)), [0 dt], y0, options);

    % Extract results
    Vf = y(:,1); nf = y(:,2); mf = y(:,3); hf = y(:,4);
    INaf = gNa_max * (mf.^3).*hf.*(Vf-ENa);
    IKf = gK_max * (nf.^4).*(Vf-EK);
end

%% Hodgkin-Huxley ODE system
function dy = hh(t,y,inj)
    % Constants
    gNa_max = 120; gK_max = 36; gLeak = 0.3; 
    ENa = 57; EK = -75; ELeak = -50; 
    Cm = 1; Vrest = -61.855;

    % Unpack y
    V = y(1); n = y(2); m = y(3); h = y(4);

    % Gating variable rate constants
    [an,bn,am,bm,ah,bh] = ab(V-Vrest);
    dn = an*(1-n)-bn*n;
    dm = am*(1-m)-bm*m;
    dh = ah*(1-h)-bh*h;

    % Currents
    INa = gNa_max*(m^3)*h*(V-ENa);
    IK = gK_max*(n^4)*(V-EK);
    ILeak = gLeak*(V-ELeak);

    % Membrane potential derivative
    dV = (inj-INa-IK-ILeak)/Cm;
    dy = [dV; dn; dm; dh];
end

%% Gating variables
function [an,bn,am,bm,ah,bh] = ab(Vm)
    an = 0.01 * (10-Vm)/(exp((10-Vm)/10)-1);
    bn = 0.125 * exp(-Vm/80);
    if Vm == 10 an = 0.1; end

    am = 0.1 * (25-Vm)/(exp((25-Vm)/10)-1);
    bm = 4 * exp(-Vm/18);
    if Vm == 25 am = 1; end

    ah = 0.07*exp(-Vm/20);
    bh = 1 / (exp((30-Vm)/10) + 1);
end

%% Gating variable time constants and steady-state values
function [tau, inf] = ti(a,b)
    tau = 1/(a+b);
    inf = a/(a+b);
end
