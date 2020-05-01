clear all;
close all;
%rho_s = [1;10^0.5;1e1;10^1.5;1e2];
rho_s = [1;1e1;1e2;1e3;1e4];
sinr = zeros(length(rho_s), 3);
sinr_onet = zeros(length(rho_s), 2);
bound = zeros(length(rho_s), 4);
Nh = 5;
Nv = 20;
M = Nh*Nv;
K = 7;
Lp = 10;%3
NZC = 3;
Nite = 1e2;
Xgp = zeros(K, Lp);
for m = 1 : K
    for l = 1 : Lp
        cml = m + l - 2;
        Xgp(m, l) = exp(1i * pi * NZC * cml * (cml + mod(Lp, 2)) / Lp);
    end
end
for k = 1 : K
    Xgp(k, :) = Xgp(k, :) / norm(Xgp(k, :));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
BS_height = 10;
ISD = 200;
Min_d = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kc = 0;
Wjj = zeros(M, M);
Kx = M;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y<=0 && (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%%%%%%cell 1
        Kc = Kc + 1;
        a = arrster(x, y, M, Nh, Nv, BS_height);
        Wjj = Wjj + (a * a');
    else
    end
end
Wjj = Wjj / Kx;
Wjl1 = zeros(M, M);
Kc = 0;
posiclx1 = ISD * cos(pi/6);
posicly1 = -ISD * sin(pi/6);
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && x<=0 && atan(abs(y/x))<=pi/3%cell 2
        Kc = Kc + 1;
        x = x + posiclx1;
        y = y + posicly1;
        a = arrster(x, y, M, Nh, Nv, BS_height);
        Wjl1 = Wjl1 + (a * a');
    else
    end
end
Wjl1 = Wjl1 / Kx;
Wjl2 = zeros(M, M);
Kc = 0;
posiclx2 = 0;
posicly2 = -ISD;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        x = x + posiclx2;
        y = y + posicly2;
        a = arrster(x, y, M, Nh, Nv, BS_height);
        Wjl2 = Wjl2 + (a * a');
    else
    end
end
Wjl2 = Wjl2 / Kx;
[Ujj, Djj] = eig(Wjj);
Ujj = Ujj(:, 1:K);
Djj = Djj(1:K, 1:K);
for n = 1 : Nite
    Hqjj = (randn(M, K) + 1i* randn(M, K)) / sqrt(2);
    Hjj = Wjj^0.5 * Hqjj;
    Hbarjj = Djj^0.5 * Ujj' * Hqjj;
    Hqjl1 = (randn(M, K) + 1i* randn(M, K)) / sqrt(2);
    Hjl1 = Wjl1^0.5 * Hqjl1;
    Hbarjl1 = Ujj' * Hjl1;
    Hqjl2 = (randn(M, K) + 1i* randn(M, K)) / sqrt(2);
    Hjl2 = Wjl2^0.5 * Hqjl2;
    Hbarjl2 = Ujj' * Hjl2;
    Zp = (randn(M, Lp) + 1i* randn(M, Lp)) / sqrt(2);
    for snr_it = 1 : length(rho_s)
        rho = rho_s(snr_it, 1);
        Sigmajk = Djj * inv(Djj + Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj + 1/rho * eye(K));
        Xgpx = Xgp * sqrt(rho);
        Ygp = Hbarjj + Hbarjl1 + Hbarjl2 + 1 / rho * Ujj' * Zp * Xgpx';
        Hbathatg = zeros(K, K);
        for m = 1 : K
            Hbathatg(:, m) = Sigmajk * Ygp(:, m);
        end
        Xij = rho * Djj * (K * eye(K) - K * Sigmajk') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
        Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
        Vxxj = inv(eye(K)+Vj);
        sinr(snr_it, 1) = sinr(snr_it, 1) + 1/real(Vxxj(1,1))-1;
        Xqjk = Xij  / rho ;
        for nxx = 2 : K
            Xqjk = Xqjk + Hbathatg(:, nxx)' * Hbathatg(:, nxx);
        end
        ajk = Xqjk \  Hbathatg(:, 1);
        Rxjq = zeros(K, K);
        for kkk = 1 : K
            xijk = Hbarjj(:, 1) - Hbathatg(:, 1);
            Rxjq = Rxjq + rho * xijk * xijk';
        end
        Rxjq = Rxjq + rho * (Hbarjl1*Hbarjl1' + Hbarjl2*Hbarjl2') + eye(K);
        interno = real(ajk' * Rxjq * ajk);
        for nxx = 2 : K
            interno = interno + rho * (abs(ajk' * Hbathatg(:, nxx)))^2;
        end
        sinr(snr_it, 2) = sinr(snr_it, 2) + rho * (abs(ajk' * Hbathatg(:, 1)))^2 / interno;
        Vj = rho * Hbathatg' * inv(Rxjq) * Hbathatg;
        Vxxj = inv(eye(K)+Vj);
        sinr(snr_it, 3) = sinr(snr_it, 3) + 1/real(Vxxj(1,1))-1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         Ygpx = Hjj + Hjl1 + Hjl2 + 1 / rho * Zp * Xgpx';
        %         Hhatg = Ygpx;
        %         Rxijk = 2 * rho * K * (Wjl1 + Wjl2) + (K+1) * eye(M);
        %         noiseandint = abs(Hhatg(:, 1)' * Rxijk * Hhatg(:, 1));
        %         for nxx = 2 : K
        %             noiseandint = noiseandint + rho * (abs(Hhatg(:, 1)' * Hhatg(:, nxx)))^2;
        %         end
        %         sinr_onet(snr_it, 1) = sinr_onet(snr_it, 1) + rho * (abs(Hhatg(:, 1)' * Hhatg(:, 1)))^2 / noiseandint;
        %%%%%%%%%%%%%%%%%%%%%
        Sigmajkx = Wjj * inv(Wjj + Wjl1 + Wjl2 + 1/rho * eye(M));
        Ygpx = Hjj + Hjl1 + Hjl2 + 1 / rho * Zp * Xgpx';
        Hhatg = zeros(M, K);
        for m = 1 : K
            Hhatg(:, m) = Sigmajkx * Ygpx(:, m);
        end
        Xijx = rho * Wjj * (K * eye(M) - K * Sigmajkx') + rho * K * (Wjl1 + Wjl2) + eye(M);
        Vjx = rho * Hhatg' * inv(Xijx) * Hhatg;
        Vxxjx = inv(eye(K)+Vjx);
        sinr_onet(snr_it, 1) = sinr_onet(snr_it, 1) + 1/real(Vxxjx(1,1))-1;
        %%%%%%%%%%%%%%
        Sigmajkx = Wjj * inv(Wjj + Wjl1 + Wjl2 + 1/rho * eye(M));
        Ygpx = Hjj + Hjl1 + Hjl2 + 1 / rho * Zp * Xgpx';
        Hhatg = zeros(M, K);
        for m = 1 : K
            Hhatg(:, m) = Sigmajkx * Ygpx(:, m);
        end
        Xijx = rho * Wjj * (K * eye(M) - K * Sigmajkx') + rho * K * (Wjl1 + Wjl2) + eye(M);
        Vjx = rho * Hhatg' * inv(Xijx) * Hhatg;
        Vxxjx = inv(eye(K)+Vjx);
        Xqjk = Xijx  / rho ;
        for nxx = 2 : K
            Xqjk = Xqjk + Hhatg(:, nxx)' * Hhatg(:, nxx);
        end
        ajk = Xqjk \  Hhatg(:, 1);
        Rxjq = zeros(M, M);
        for kkk = 1 : K
            xijk = Hjj(:, 1) - Hhatg(:, 1);
            Rxjq = Rxjq + rho * xijk * xijk';
        end
        Rxjq = Rxjq + rho * (Hjl1*Hjl1' + Hjl2*Hjl2') + eye(M);
        interno = real(ajk' * Rxjq * ajk);
        for nxx = 2 : K
            interno = interno + rho * (abs(ajk' * Hhatg(:, nxx)))^2;
        end
        sinr_onet(snr_it, 2) = sinr_onet(snr_it, 2) + rho * (abs(ajk' * Hhatg(:, 1)))^2 / interno;
        %%%%%%%%%%%%%%%%%%5
    end
    disp([n])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
sinr = sinr / Nite;
sinr_onet = sinr_onet / Nite;
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_ss = 10 * log10(rho_s);
sinr_ss = 10 * log10(sinr);
sinr_onet_ss = 10 * log10(sinr_onet);
figure
h1=plot(rho_ss, sinr_ss(:, 1), 'k-^','LineWidth',1,'MarkerSize',10);
hold on
plot(rho_ss, sinr_onet_ss(:, 1), 'k-*','LineWidth',1,'MarkerSize',10);
plot(rho_ss, sinr_ss(:, 2), 'g-^','LineWidth',1,'MarkerSize',10);
plot(rho_ss, sinr_onet_ss(:, 2), 'g-*','LineWidth',1,'MarkerSize',10);
plot(rho_ss, sinr_ss(:, 3), 'r-o','LineWidth',1,'MarkerSize',10);
xlim([min(rho_ss), max(rho_ss)])
le = legend('Approx-two-tier','Approx-full-dim', 'True-two-tier','True-full-dim','halftrue-two-tier','Location', 'northwest');
set(gca,'XTick',rho_ss)
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('SNR (dB)','Fontname','Times','Fontsize',12)
ylabel('SINR (dB)','Fontname','Times','Fontsize',12)
title('Uniform')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
