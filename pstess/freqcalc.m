function freqcalc(k,t,vmod)
% Author: Felipe Wilches-Bernal (fwilche@sandia.gov)
% Sandia National Laboratories
% Feb. 2017
% Function developed to calculate bus frequencies

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

Nbf_hpf = 10;
K1bf_hpf = 1;

% b0 = 0.00000000000000;
b1 = 277.237187626600;
a0 = 277.237187626600;
a1 = 28.8394098913241;

if (vmod == 0)

    g.freq.kpx = ones(size(t));
    g.freq.bus_freq = zeros(size(g.bus.theta));
    g.freq.bus_freq(:,1) = 1;

    %=====================================================================%

    g.freq.bf_hpf = zeros(size(g.bus.theta));
    g.freq.dbf_hpf = zeros(size(g.bus.theta));
    g.freq.bus_freqf = zeros(size(g.bus.theta));
    g.freq.bus_freqf(:,1) = 1;
    g.freq.o2p = zeros(size(g.bus.theta,1),1);
    g.freq.bf_hpf(:,1) = K1bf_hpf.*g.bus.theta(:,1);

    %** Sandia filter
    g.freq.bus_freqsnl = zeros(size(g.bus.theta));
    g.freq.bus_freqsnl(:,1) = 1;
    g.freq.x1_snlf = zeros(size(g.bus.theta));
    g.freq.dx1_snlf = zeros(size(g.bus.theta));
    g.freq.x2_snlf = zeros(size(g.bus.theta));
    g.freq.dx2_snlf = zeros(size(g.bus.theta));

elseif (vmod == 1)

    if (k ~= 1)

        dVangle = g.bus.theta(:,k) - g.bus.theta(:,k-1);
        dVangle(dVangle >= pi) = dVangle(dVangle >= pi) - 2*pi;
        dVangle(dVangle <= -pi) = dVangle(dVangle <= -pi) + 2*pi;

        dbus_freq = dVangle/(t(k)-t(k-1))*(1/g.sys.basrad);

        g.freq.bus_freq(:,k) = g.freq.bus_freq(:,1) + dbus_freq;

        %=====================================================================%

        ixp2p = g.bus.theta(:,k)-g.bus.theta(:,k-1) >= pi;
        ixm2p = g.bus.theta(:,k)-g.bus.theta(:,k-1) <= -pi;
        g.freq.o2p(ixp2p) = g.freq.o2p(ixp2p) + 1*g.freq.kpx(k);
        g.freq.o2p(ixm2p) = g.freq.o2p(ixm2p) - 1*g.freq.kpx(k);
        theta_un_k = g.bus.theta(:,k) - 2*pi.*g.freq.o2p;

        g.freq.dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta_un_k ...
                              - Nbf_hpf.*g.freq.bf_hpf(:,k);
        g.freq.bus_freqf(:,k) = g.freq.bus_freqf(:,1) ...
                                + g.freq.dbf_hpf(:,k)*(1/g.sys.basrad);
        g.freq.kpx(k) = 0;

        %** Sandia filter
        g.freq.dx1_snlf(:,k) = -a0*g.freq.x2_snlf(:,k);
        g.freq.dx2_snlf(:,k) = b1*theta_un_k + g.freq.x1_snlf(:,k) ...
                               - a1*g.freq.x2_snlf(:,k);
        g.freq.bus_freqsnl(:,k) = g.freq.bus_freqsnl(:,1) ...
                                  + g.freq.x2_snlf(:,k)*(1/g.sys.basrad);

    else

        g.freq.dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*g.bus.theta(:,k) ...
                              - Nbf_hpf.*g.freq.bf_hpf(:,k);
        g.freq.bus_freqf(:,k) = g.freq.bus_freqf(:,1) ...
                                + g.freq.dbf_hpf(:,k)*(1/g.sys.basrad);

        % when k==1, the vector is already initialized

        %** Sandia filter
        g.freq.x1_snlf(:,1) = -b1*g.bus.theta(:,k);
        g.freq.x2_snlf(:,1) = 0;

        g.freq.dx1_snlf(:,k) = -a0*g.freq.x2_snlf(:,k);
        g.freq.dx2_snlf(:,k) = b1*g.bus.theta(:,k) ...
                               + g.freq.x1_snlf(:,k) - a1*g.freq.x2_snlf(:,k);
        g.freq.bus_freqsnl(:,k) = g.freq.bus_freqsnl(:,1) ...
                                  + g.freq.x2_snlf(:,k)*(1/g.sys.basrad);

    end

end

end  % function end

% eof
