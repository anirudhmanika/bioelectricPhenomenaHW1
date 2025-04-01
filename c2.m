% HW1 #1
% .m function for solving the hodgkin huxley model
% param: dt -> the total simulation duration, in ms
% Param: currentStep -> function for input current
% returns: membrane voltage, ionic currents, and gating variable
% probabilities throughout the simulation

function [t,Vf,nf,mf,hf,INaf,IKf] = solve_hh(dt, currentStep)
    % constants
    gNa_max = 120; gK_max = 36; gLeak = 0.3; 
    ELeak = -50; ENa = 57; EK = -75;
   
    Cm = 1; Vrest = -63;

    % initial n,m,h
    [an0,bn0,am0,bm0,ah0,bh0] = ab(0);
    [~,n0] = ti(an0,bn0); [~,m0] = ti(am0,bm0); [~,h0] = ti(ah0,bh0);
    y0 = [Vrest, n0, m0, h0];

    options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',0.01);
    [t,y] = ode45(@(t,y) hh(t,y,currentStep(t)), [0 dt], y0, options);

    Vf = y(:,1); nf = y(:,2); mf = y(:,3); hf = y(:,4);
    INaf = gNa_max * (mf.^3).*hf.*(Vf-ENa);
    IKf = gK_max * (nf.^4).*(Vf-EK);

    % derivative function
    function dy = hh(t,y,inj)
        % unpack y
        V = y(1); n = y(2); m = y(3); h = y(4);

        % gate rate constants
        [an,bn,am,bm,ah,bh] = ab(V-Vrest);
        % gate derivatives
        dn = an*(1-n)-bn*n;
        dm = am*(1-m)-bm*m;
        dh = ah*(1-h)-bh*h;

        % currents
        INa = gNa_max*(m^3)*h*(V-ENa);
        IK = gK_max*(n^4)*(V-EK);
        ILeak = gLeak*(V-ELeak);

        dV = (inj-INa-IK-ILeak)/Cm;
        dy = [dV; dn; dm; dh];
    end
    
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
    
    function [tau, inf] = ti(a,b)
        tau = 1/(a+b);
        inf = a/(a+b);
    end
end