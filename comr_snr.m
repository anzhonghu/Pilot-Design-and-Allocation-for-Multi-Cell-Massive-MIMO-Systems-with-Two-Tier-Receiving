clear all;
close all;
rho_s = [1;1e1;1e2;1e3;1e4];
sumrate = zeros(length(rho_s), 8);
bound = zeros(length(rho_s), 2);
rho = rho_s(5, 1);
Nh = 5;
Nv = 20;
M = Nh*Nv;
K = 7;
Lp = 10;%3
L = 3;
NZC = 3;
Nite = 1e2;
Npd = 100;
Xgp = zeros(K, Lp);
for m = 1 : K
    for l = 1 : Lp
        cml = m + l - 2;
        Xgp(m, l) = exp(1i * pi * NZC * cml * (cml + mod(Lp, 2)) / Lp);
    end
end
Fn = zeros(K, Lp*Npd);
for k = 1 : K
    Xgp(k, :) = Xgp(k, :) / norm(Xgp(k, :));
end
for n = 1 : Npd
    an = randn;
    for k = 1 : K
        for m = 1 : Lp
            Fn(k, (n-1)*Lp+m) = exp(1i*2*pi/Lp*(k+an)*m);
        end
        Fn(k, (n-1)*Lp+1:n*Lp) = Fn(k, (n-1)*Lp+1:n*Lp) / norm(Fn(k, (n-1)*Lp+1:n*Lp));
    end
end
XgpF = zeros(K, Lp*L);
Xgpmin = minphase(K, Lp, L);
BS_height = 10;
ISD = 200;
Min_d = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
posi_sto = zeros(M*L, 2);
bs_posi(1,:) = [0,0];
bs_posi(2,:) = [ISD * cos(pi/6),-ISD * sin(pi/6)];
bs_posi(3,:) = [0,-ISD];
Kc = 0;
a = zeros(M, 1);
Kx = M;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y<=0 && (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%%%%%%cell 1
        Kc = Kc + 1;
        posi_sto(Kc,:) = [x,y];
    else
    end
end
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && x<=0 && atan(abs(y/x))<=pi/3%cell 2
        Kc = Kc + 1;
        posi_sto(M+Kc,:) = [x,y];
    else
    end
end
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        posi_sto(2*M+Kc,:) = [x,y];
    else
    end
end
W11 = zeros(M, M);
W12 = zeros(M, M);
W13 = zeros(M, M);
W22 = zeros(M, M);
W21 = zeros(M, M);
W23 = zeros(M, M);
W33 = zeros(M, M);
W31 = zeros(M, M);
W32 = zeros(M, M);
for l = 1 : L
    for k = 1 : M
        x = posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(1,1);
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(1,2);
        a = arrster(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W11 = W11 + (a * a');
            case 2
                W12 = W12 + (a * a');
            case 3
                W13 = W13 + (a * a');
        end
    end
end
for l = 1 : L
    for k = 1 : M
        x = abs(posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(2,1));
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(2,2);
        a = arrsterc2(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W21 = W21 + (a * a');
            case 2
                W22 = W22 + (a * a');
            case 3
                W23 = W23 + (a * a');
        end
    end
end
for l = 1 : L
    for k = 1 : M
        x = abs(posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(3,1));
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(3,2);
        a = arrsterc3(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W31 = W31 + (a * a');
            case 2
                W32 = W32 + (a * a');
            case 3
                W33 = W33 + (a * a');
        end
    end
end
W11 = W11 / Kx;
W12 = W12 / Kx;
W13 = W13 / Kx;
W21 = W21 / Kx;
W22 = W22 / Kx;
W23 = W23 / Kx;
W31 = W31 / Kx;
W32 = W32 / Kx;
W33 = W33 / Kx;
%%%%%%%%%%%%%%%%%%%%%%
UUcor = zeros(L*2,1);
for l = 1 : L
    switch l
        case 1
            Wjj = W11;
            Wjl1 = W12;
            Wjl2 = W13;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 2
            Wjj = W22;
            Wjl1 = W21;
            Wjl2 = W23;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 3
            Wjj = W33;
            Wjl1 = W31;
            Wjl2 = W32;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
    end
    UUcor((l-1)*2+1:l*2,:) = [norm(Ujj' * Ujl1); norm(Ujj' * Ujl2)];
end
[~, cor_in] = max(UUcor);
j0 = floor((cor_in-1)/2)+1;
j1 = rem(cor_in,2);
if j1==0
    j1=2;
end
switch j0
    case 1
        j1 = j1 + 1;
    case 2
        if j1 < 2
        else
            j1 = j1 + 1;
        end
end
j2 = 6-j0-j1;
XgpF(:,(j0-1)*Lp+1:j0*Lp) = Fn(:, 1:Lp);
for cellr = 1 : 2
    upbur = zeros(Npd, 1);
    switch cellr
        case 1
            jxx = j1;
        case 2
            jxx = j2;
    end
    for n = 1 : Npd
        XgpF(:,(jxx-1)*Lp+1:jxx*Lp) = Fn(:, (n-1)*Lp+1:n*Lp);
        for l = 1 : L
            switch l
                case 1
                    Wjj = W11;
                    Wjl1 = W12;
                    Wjl2 = W13;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
                case 2
                    Wjj = W22;
                    Wjl1 = W21;
                    Wjl2 = W23;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
                case 3
                    Wjj = W33;
                    Wjl1 = W31;
                    Wjl2 = W32;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
            end
            Xgpx = XgpF(:,(l-1)*Lp+1:l*Lp) * sqrt(rho);
            Sigtem = zeros(K,K*K);
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpF(:,(ll-1)*Lp+1:ll*Lp)' * sqrt(rho);
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Sigsum = zeros(K, K);
            for m = 1 : K
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Xjl1 = ones(K, K);
            Xjl2 = ones(K, K);
            for k = 1 : K
                for kk = 1 : K
                    Xjl1(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl1);
                    Xjl2(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl2);
                end
            end
            Yj = eye(K);
            for k = 1 : K
                Yj(k,k) = (rho * trace(Djj * Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K))...
                    + trace(Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K)));
            end
            Rj = Yj;
            llc = 0;
            for ll = 1 : L
                if ll==l
                else
                    Cjl1 = Xgpx * XgpF(:,(ll-1)*Lp+1:ll*Lp)' * sqrt(rho);
                    llc = llc + 1;
                    if llc==1
                        Rj = Rj + Xjl1 .* (Cjl1'* Cjl1);
                    else
                        Rj = Rj + Xjl2 .* (Cjl1'* Cjl1);
                    end
                end
            end
            Rj = real(Rj);
            inR = inv(Rj);
            for k = 1 : K
                Rxx = [Rj(1:k-1,:);Rj(k+1:K,:)];
                Ryy = [Rxx(:,1:k-1),Rxx(:,k+1:K)];
                upbur(n,1) = upbur(n,1) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/K/inR(k,k));
            end
            disp('first-pilot')
            disp([cellr,n,l])
        end
    end
    [~,nq] = max(upbur);
    XgpF(:,(jxx-1)*Lp+1:jxx*Lp) = Fn(:, (nq-1)*Lp+1:nq*Lp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snr_it = 1 : length(rho_s)
    rho = rho_s(snr_it, 1);
    XgpFt = XgpF * sqrt(rho);
    for l = 1 : L
        switch l
            case 1
                Wjj = W11;
                Wjl1 = W12;
                Wjl2 = W13;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
            case 2
                Wjj = W22;
                Wjl1 = W21;
                Wjl2 = W23;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
            case 3
                Wjj = W33;
                Wjl1 = W31;
                Wjl2 = W32;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
        end
        Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
        Sigtem = zeros(K,K*K);
        for k = 1 : K
            clco = 0;
            for ll = 1 : L
                if l==ll
                    Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                else
                    clco = clco + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                end
                if 1==clco
                    Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                else
                    if 2==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    end
                end
            end
            Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
        end
        Sigsum = zeros(K, K);
        for m = 1 : K
            Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
        end
        Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
        Xjl1 = ones(K, K);
        Xjl2 = ones(K, K);
        for k = 1 : K
            for kk = 1 : K
                Xjl1(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl1);
                Xjl2(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl2);
            end
        end
        Yj = eye(K);
        for k = 1 : K
            Yj(k,k) = (rho * trace(Djj * Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K))...
                + trace(Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K)));
        end
        Rj = Yj;
        llc = 0;
        for ll = 1 : L
            if ll==l
            else
                Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                llc = llc + 1;
                if llc==1
                    Rj = Rj + Xjl1 .* (Cjl1'* Cjl1);
                else
                    Rj = Rj + Xjl2 .* (Cjl1'* Cjl1);
                end
            end
        end
        Rj = real(Rj);
        inR = inv(Rj);
        for k = 1 : K
            Rxx = [Rj(1:k-1,:);Rj(k+1:K,:)];
            Ryy = [Rxx(:,1:k-1),Rxx(:,k+1:K)];
            bound(snr_it,1) = bound(snr_it,1) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/K/inR(k,k));
        end
        disp('first-bound')
        disp([snr_it,l])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1 : L
    switch l
        case 1
            Wjj = W11;
            Wjl1 = W12;
            Wjl2 = W13;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 2
            Wjj = W22;
            Wjl1 = W21;
            Wjl2 = W23;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 3
            Wjj = W33;
            Wjl1 = W31;
            Wjl2 = W32;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
    end
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Sigmajk = Djj * inv(Djj + Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj +1/rho * eye(K));
            Xgpx = Xgp * sqrt(rho);
            Ygp = Hbarjj + Hbarjl1 + Hbarjl2 + 1 / rho * Ujj' * Zp * Xgpx';
            Hbathatg = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigmajk * Ygp(:, m);
            end
            Xij = rho * Djj * (K * eye(K) - K * Sigmajk') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 1) = sumrate(snr_it, 1) + log2(1/real(Vxxj(k,k)));%%reused
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XgpFt = Xgpmin * sqrt(rho);
            Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Sigtem = zeros(K,K*K);
            cll = 0;
            for ll = 1 : L
                if ll==l
                else
                    cll = cll + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    if cll==1
                        Ygp = Ygp + 1 / rho * Hbarjl1 * Cjl1;
                    else
                        if cll==2
                            Ygp = Ygp + 1 / rho * Hbarjl2 * Cjl1;
                        end
                    end
                end
            end
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Hbathatg = zeros(K, K);
            Sigsum = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigtem(:,(k-1)*K+1:k*K) * Ygp(:, m);
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 2) = sumrate(snr_it, 2) + log2(1/real(Vxxj(k,k)));%%proposed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XgpFt = XgpF * sqrt(rho);
            Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Sigtem = zeros(K,K*K);
            cll = 0;
            for ll = 1 : L
                if ll==l
                else
                    cll = cll + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    if cll==1
                        Ygp = Ygp + 1 / rho * Hbarjl1 * Cjl1;
                    else
                        if cll==2
                            Ygp = Ygp + 1 / rho * Hbarjl2 * Cjl1;
                        end
                    end
                end
            end
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Hbathatg = zeros(K, K);
            Sigsum = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigtem(:,(k-1)*K+1:k*K) * Ygp(:, m);
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 3) = sumrate(snr_it, 3) + log2(1/real(Vxxj(k,k)));%%proposed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Sigmajk = Djj * inv(Djj + 1/rho * eye(K));
            Xgpx = Xgp * sqrt(rho);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Hbathatg = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigmajk * Ygp(:, m);
            end
            Xij = rho * Djj * (K * eye(K) - K * Sigmajk') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 4) = sumrate(snr_it, 4) + log2(1/real(Vxxj(k,k)));%%reused
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('first-cal')
            disp([l,n,snr_it])
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
posi_sto = zeros(K*L, 2);
bs_posi(1,:) = [0,0];
bs_posi(2,:) = [ISD * cos(pi/6),-ISD * sin(pi/6)];
bs_posi(3,:) = [0,-ISD];
Kc = 0;
a = zeros(M, 1);
Kx = M;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (50)^2 && x^2+y^2>=Min_d^2 && y<=0 && (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%%%%%%cell 1
        Kc = Kc + 1;
        posi_sto(Kc,:) = [x,y];
    else
    end
end
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (50)^2 && x^2+y^2>=Min_d^2 && x<=0 && atan(abs(y/x))<=pi/3%cell 2
        Kc = Kc + 1;
        posi_sto(M+Kc,:) = [x,y];
    else
    end
end
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (50)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        posi_sto(2*M+Kc,:) = [x,y];
    else
    end
end
W11 = zeros(M, M);
W12 = zeros(M, M);
W13 = zeros(M, M);
W22 = zeros(M, M);
W21 = zeros(M, M);
W23 = zeros(M, M);
W33 = zeros(M, M);
W31 = zeros(M, M);
W32 = zeros(M, M);
for l = 1 : L
    for k = 1 : M
        x = posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(1,1);
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(1,2);
        a = arrster(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W11 = W11 + (a * a');
            case 2
                W12 = W12 + (a * a');
            case 3
                W13 = W13 + (a * a');
        end
    end
end
for l = 1 : L
    for k = 1 : M
        x = abs(posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(2,1));
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(2,2);
        a = arrsterc2(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W21 = W21 + (a * a');
            case 2
                W22 = W22 + (a * a');
            case 3
                W23 = W23 + (a * a');
        end
    end
end
for l = 1 : L
    for k = 1 : M
        x = abs(posi_sto((l-1)*M+k,1)+bs_posi(l,1)-bs_posi(3,1));
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2)-bs_posi(3,2);
        a = arrsterc3(x, y, M, Nh, Nv, BS_height);
        switch l
            case 1
                W31 = W31 + (a * a');
            case 2
                W32 = W32 + (a * a');
            case 3
                W33 = W33 + (a * a');
        end
    end
end
W11 = W11 / Kx;
W12 = W12 / Kx;
W13 = W13 / Kx;
W21 = W21 / Kx;
W22 = W22 / Kx;
W23 = W23 / Kx;
W31 = W31 / Kx;
W32 = W32 / Kx;
W33 = W33 / Kx;
%%%%%%%%%%%%%%%%%%%%%%
UUcor = zeros(L*2,1);
for l = 1 : L
    switch l
        case 1
            Wjj = W11;
            Wjl1 = W12;
            Wjl2 = W13;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 2
            Wjj = W22;
            Wjl1 = W21;
            Wjl2 = W23;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 3
            Wjj = W33;
            Wjl1 = W31;
            Wjl2 = W32;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
    end
    UUcor((l-1)*2+1:l*2,:) = [norm(Ujj' * Ujl1); norm(Ujj' * Ujl2)];
end
[~, cor_in] = max(UUcor);
j0 = floor((cor_in-1)/2)+1;
j1 = rem(cor_in,2);
if j1==0
    j1=2;
end
switch j0
    case 1
        j1 = j1 + 1;
    case 2
        if j1 < 2
        else
            j1 = j1 + 1;
        end
end
j2 = 6-j0-j1;
XgpF(:,(j0-1)*Lp+1:j0*Lp) = Fn(:, 1:Lp);
for cellr = 1 : 2
    upbur = zeros(Npd, 1);
    switch cellr
        case 1
            jxx = j1;
        case 2
            jxx = j2;
    end
    for n = 1 : Npd
        XgpF(:,(jxx-1)*Lp+1:jxx*Lp) = Fn(:, (n-1)*Lp+1:n*Lp);
        for l = 1 : L
            switch l
                case 1
                    Wjj = W11;
                    Wjl1 = W12;
                    Wjl2 = W13;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
                case 2
                    Wjj = W22;
                    Wjl1 = W21;
                    Wjl2 = W23;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
                case 3
                    Wjj = W33;
                    Wjl1 = W31;
                    Wjl2 = W32;
                    [Ujj, Djj] = eig(Wjj);
                    Ujj = Ujj(:, 1:K);
                    Djj = Djj(1:K, 1:K);
                    [Ujl1, ~] = eig(Wjl1);
                    [Ujl2, ~] = eig(Wjl2);
            end
            Xgpx = XgpF(:,(l-1)*Lp+1:l*Lp) * sqrt(rho);
            Sigtem = zeros(K,K*K);
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpF(:,(ll-1)*Lp+1:ll*Lp)' * sqrt(rho);
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Sigsum = zeros(K, K);
            for m = 1 : K
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Xjl1 = ones(K, K);
            Xjl2 = ones(K, K);
            for k = 1 : K
                for kk = 1 : K
                    Xjl1(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl1);
                    Xjl2(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl2);
                end
            end
            Yj = eye(K);
            for k = 1 : K
                Yj(k,k) = (rho * trace(Djj * Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K))...
                    + trace(Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K)));
            end
            Rj = Yj;
            llc = 0;
            for ll = 1 : L
                if ll==l
                else
                    Cjl1 = Xgpx * XgpF(:,(ll-1)*Lp+1:ll*Lp)' * sqrt(rho);
                    llc = llc + 1;
                    if llc==1
                        Rj = Rj + Xjl1 .* (Cjl1'* Cjl1);
                    else
                        Rj = Rj + Xjl2 .* (Cjl1'* Cjl1);
                    end
                end
            end
            Rj = real(Rj);
            inR = inv(Rj);
            for k = 1 : K
                Rxx = [Rj(1:k-1,:);Rj(k+1:K,:)];
                Ryy = [Rxx(:,1:k-1),Rxx(:,k+1:K)];
                upbur(n,1) = upbur(n,1) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/K/inR(k,k));
            end
            disp('second-pilot')
            disp([cellr,n,l])
        end
    end
    [~,nq] = max(upbur);
    XgpF(:,(jxx-1)*Lp+1:jxx*Lp) = Fn(:, (nq-1)*Lp+1:nq*Lp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snr_it = 1 : length(rho_s)
    rho = rho_s(snr_it, 1);
    XgpFt = XgpF * sqrt(rho);
    for l = 1 : L
        switch l
            case 1
                Wjj = W11;
                Wjl1 = W12;
                Wjl2 = W13;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
            case 2
                Wjj = W22;
                Wjl1 = W21;
                Wjl2 = W23;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
            case 3
                Wjj = W33;
                Wjl1 = W31;
                Wjl2 = W32;
                [Ujj, Djj] = eig(Wjj);
                Ujj = Ujj(:, 1:K);
                Djj = Djj(1:K, 1:K);
                [Ujl1, ~] = eig(Wjl1);
                [Ujl2, ~] = eig(Wjl2);
        end
        Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
        Sigtem = zeros(K,K*K);
        for k = 1 : K
            clco = 0;
            for ll = 1 : L
                if l==ll
                    Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                else
                    clco = clco + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                end
                if 1==clco
                    Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                else
                    if 2==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    end
                end
            end
            Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
        end
        Sigsum = zeros(K, K);
        for m = 1 : K
            Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
        end
        Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
        Xjl1 = ones(K, K);
        Xjl2 = ones(K, K);
        for k = 1 : K
            for kk = 1 : K
                Xjl1(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl1);
                Xjl2(k, kk) = 1 / rho * trace(Ujj *Sigtem(:,(k-1)*K+1:k*K)' * inv(Xij)*Sigtem(:,(kk-1)*K+1:kk*K) * Ujj'*Wjl2);
            end
        end
        Yj = eye(K);
        for k = 1 : K
            Yj(k,k) = (rho * trace(Djj * Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K))...
                + trace(Sigtem(:,(k-1)*K+1:k*K)'* inv(Xij)*Sigtem(:,(k-1)*K+1:k*K)));
        end
        Rj = Yj;
        llc = 0;
        for ll = 1 : L
            if ll==l
            else
                Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                llc = llc + 1;
                if llc==1
                    Rj = Rj + Xjl1 .* (Cjl1'* Cjl1);
                else
                    Rj = Rj + Xjl2 .* (Cjl1'* Cjl1);
                end
            end
        end
        Rj = real(Rj);
        inR = inv(Rj);
        for k = 1 : K
            Rxx = [Rj(1:k-1,:);Rj(k+1:K,:)];
            Ryy = [Rxx(:,1:k-1),Rxx(:,k+1:K)];
            bound(snr_it,2) = bound(snr_it,2) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/K/inR(k,k));
        end
        disp('second-bound')
        disp([snr_it,l])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1 : L
    switch l
        case 1
            Wjj = W11;
            Wjl1 = W12;
            Wjl2 = W13;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 2
            Wjj = W22;
            Wjl1 = W21;
            Wjl2 = W23;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
        case 3
            Wjj = W33;
            Wjl1 = W31;
            Wjl2 = W32;
            [Ujj, Djj] = eig(Wjj);
            Ujj = Ujj(:, 1:K);
            Djj = Djj(1:K, 1:K);
            [Ujl1, ~] = eig(Wjl1);
            [Ujl2, ~] = eig(Wjl2);
    end
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            for k = 1 : K
                sumrate(snr_it, 5) = sumrate(snr_it, 5) + log2(1/real(Vxxj(k,k)));%%reused
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XgpFt = Xgpmin * sqrt(rho);
            Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Sigtem = zeros(K,K*K);
            cll = 0;
            for ll = 1 : L
                if ll==l
                else
                    cll = cll + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    if cll==1
                        Ygp = Ygp + 1 / rho * Hbarjl1 * Cjl1;
                    else
                        if cll==2
                            Ygp = Ygp + 1 / rho * Hbarjl2 * Cjl1;
                        end
                    end
                end
            end
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Hbathatg = zeros(K, K);
            Sigsum = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigtem(:,(k-1)*K+1:k*K) * Ygp(:, m);
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 6) = sumrate(snr_it, 6) + log2(1/real(Vxxj(k,k)));%%proposed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XgpFt = XgpF * sqrt(rho);
            Xgpx = XgpFt(:,(l-1)*Lp+1:l*Lp);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Sigtem = zeros(K,K*K);
            cll = 0;
            for ll = 1 : L
                if ll==l
                else
                    cll = cll + 1;
                    Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    if cll==1
                        Ygp = Ygp + 1 / rho * Hbarjl1 * Cjl1;
                    else
                        if cll==2
                            Ygp = Ygp + 1 / rho * Hbarjl2 * Cjl1;
                        end
                    end
                end
            end
            for k = 1 : K
                clco = 0;
                for ll = 1 : L
                    if l==ll
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Djj;
                    else
                        clco = clco + 1;
                        Cjl1 = Xgpx * XgpFt(:,(ll-1)*Lp+1:ll*Lp)';
                    end
                    if 1==clco
                        Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl1 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                    else
                        if 2==clco
                            Sigtem(:,(k-1)*K+1:k*K) = Sigtem(:,(k-1)*K+1:k*K) + Ujj' * Wjl2 * Ujj * (Cjl1(:,k)'*Cjl1(:,k))/rho^2;
                        end
                    end
                end
                Sigtem(:,(k-1)*K+1:k*K)  = Djj * inv(Sigtem(:,(k-1)*K+1:k*K) + 1/rho * eye(K));
            end
            Hbathatg = zeros(K, K);
            Sigsum = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigtem(:,(k-1)*K+1:k*K) * Ygp(:, m);
                Sigsum = Sigsum + Sigtem(:,(k-1)*K+1:k*K);
            end
            Xij = rho * Djj * (K * eye(K) - Sigsum') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 7) = sumrate(snr_it, 7) + log2(1/real(Vxxj(k,k)));%%proposed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Sigmajk = Djj * inv(Djj + 1/rho * eye(K));
            Xgpx = Xgp * sqrt(rho);
            Ygp = Hbarjj + 1 / rho * Ujj' * Zp * Xgpx';
            Hbathatg = zeros(K, K);
            for m = 1 : K
                Hbathatg(:, m) = Sigmajk * Ygp(:, m);
            end
            Xij = rho * Djj * (K * eye(K) - K * Sigmajk') + rho * K * (Ujj' * Wjl1 * Ujj + Ujj' * Wjl2 * Ujj) + eye(K);
            Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
            Vxxj = inv(eye(K)+Vj);
            for k = 1 : K
                sumrate(snr_it, 8) = sumrate(snr_it, 8) + log2(1/real(Vxxj(k,k)));%%reused
            end
            %%%%%%%%%%%%%%%%%%%%
            disp('second-cal')
            disp([l,n,snr_it])
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumrate = sumrate / Nite;
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_ss = 10 * log10(rho_s);
figure;
subplot(1,2,1)
h1=semilogy(rho_ss, sumrate(:, 1), 'k--s','LineWidth',1,'MarkerSize',10);
hold on
semilogy(rho_ss, sumrate(:, 2), 'k--x','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, sumrate(:, 3), 'k-p','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, sumrate(:, 4), 'k-^','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, bound(:, 1), 'k-o','LineWidth',1,'MarkerSize',10)%%upper bound
xlim([min(rho_ss), max(rho_ss)])
ylim([1, 1000])
le = legend('Reuse','Min-phase [6][7]','Proposed','Ideal','Upper bound (39)','Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',rho_ss)
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('SNR (dB)','Fontname','Times','Fontsize',12)
ylabel('Sum rate (bps/Hz)','Fontname','Times','Fontsize',12)
title('Uniform')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
semilogy(rho_ss, sumrate(:, 5), 'k--s','LineWidth',1,'MarkerSize',10);
hold on
semilogy(rho_ss, sumrate(:, 6), 'k--x','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, sumrate(:, 7), 'k-p','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, sumrate(:, 8), 'k-^','LineWidth',1,'MarkerSize',10);
semilogy(rho_ss, bound(:, 2), 'k-o','LineWidth',1,'MarkerSize',10)%%upper bound
xlim([min(rho_ss), max(rho_ss)])
ylim([1, 1000])
le = legend('Reuse','Min-phase [6][7]','Proposed','Ideal','Upper bound (39)','Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',rho_ss)
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('SNR (dB)','Fontname','Times','Fontsize',12)
ylabel('Sum rate (bps/Hz)','Fontname','Times','Fontsize',12)
title('Nonuniform')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
