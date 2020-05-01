clear all;
close all;
ant_s = [10;30;50;70;90;110;130];
rho = 1e2;
sumrate = zeros(length(ant_s), 5);
bound = zeros(length(ant_s), 2);
Nh = 1;
Nv = 20;
M = Nh*Nv;
K = 6;
G = 2;
Kq = 3;
Lp = 10;%3
L = 3;
NZC = 3;
Nite = 1e3;
Npd = 100;
disp('pilot-b')
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
Xgpp = zeros(K, Lp*L);
for l = 1 : L
    Xgpp(:, (l-1)*Lp+1:l*Lp) = Xgp;
end
Fn1 = zeros(K, Lp*Npd);
for n = 1 : Npd
    an = randn;
    for k = 1 : K
        for m = 1 : Lp
            Fn1(k, (n-1)*Lp+m) = exp(1i*2*pi/Lp*(k+an)*m);
        end
        Fn1(k, (n-1)*Lp+1:n*Lp) = Fn1(k, (n-1)*Lp+1:n*Lp) / norm(Fn1(k, (n-1)*Lp+1:n*Lp));
    end
end
disp('pilot-o')
Xgpmin = minphase(K, Lp, L);
XgpF1 = zeros(K, Lp*L);
BS_height = 25;
ISD = 500;
Min_d = 35;
range = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
posi_sto = zeros(2*K*L, 2);
bs_posi(1,:) = [0,0];
bs_posi(2,:) = [ISD * cos(pi/6),-ISD * sin(pi/6)];
bs_posi(3,:) = [0,-ISD];
midd = [ISD * cos(pi/6) * 0.5, -ISD * 0.5];
Kc = 0;
pcg11 = midd(1, 1) / 3;
pcg12 = midd(1, 2) / 3;
pcg21 = pcg11 * 2;
pcg22 = pcg12 * 2;
pcg31  = (midd(1,1)  - bs_posi(2,1)) / 3;
pcg32  = (midd(1,2)  - bs_posi(2,2)) / 3;
pcg41 = pcg31 * 2;
pcg42 = pcg32 * 2;
pcg51 = (midd(1,1)  - bs_posi(3,1)) / 3;
pcg52 = (midd(1,2)  - bs_posi(3,2)) / 3;
pcg61 = pcg51 * 2;
pcg62 = pcg52 * 2;
Kx = 100;
disp('posi-b')
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg11;
    y = (rand(1,1)-0.5) * range+pcg12;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y<=0 && (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%%%%%%cell 1
        Kc = Kc + 1;
        posi_sto(Kc,:) = [x,y];
    else
    end
end
disp('posi-o1')
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg21;
    y = (rand(1,1)-0.5) * range+pcg22;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y<=0 && (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%%%%%%cell 1
        Kc = Kc + 1;
        posi_sto(Kx+Kc,:) = [x,y];
    else
    end
end
disp('posi-o2')
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg31;
    y = (rand(1,1)-0.5) * range+pcg32;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && x<=0 && atan(abs(y/x))<=pi/3%cell 2
        Kc = Kc + 1;
        posi_sto(2*Kx+Kc,:) = [x,y];
    else
    end
end
disp('posi-o3')
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg41;
    y = (rand(1,1)-0.5) * range+pcg42;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && x<=0 && atan(abs(y/x))<=pi/3%cell 2
        Kc = Kc + 1;
        posi_sto(3*Kx+Kc,:) = [x,y];
    else
    end
end
disp('posi-o4')
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg51;
    y = (rand(1,1)-0.5) * range+pcg52;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        posi_sto(4*Kx+Kc,:) = [x,y];
    else
    end
end
disp('posi-o5')
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * range+pcg61;
    y = (rand(1,1)-0.5) * range+pcg62;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        posi_sto(5*Kx+Kc,:) = [x,y];
    else
    end
end
disp('posi-o6')
pilot_seg;
for ant_t = 1 : length(ant_s)
    Nv = ant_s(ant_t,1);
    M = Nv*Nh;
    a = zeros(M, 1);
    disp('wss-b1')
    Wss = zeros(M, M*G*L*L);
    for l = 1 : L
        for g = 1 : G
            for k = 1 : Kx
                x = posi_sto((l-1)*G*Kx+(g-1)*Kx+k,1)+bs_posi(l,1);
                y = posi_sto((l-1)*G*Kx+(g-1)*Kx+k,2)+bs_posi(l,2);
                a = arrster(x, y, M, Nh, Nv, BS_height);
                Wss(:,M*G*(l-1)+(g-1)*M+1:M*G*(l-1)+g*M) = Wss(:,M*G*(l-1)+(g-1)*M+1:M*G*(l-1)+g*M) + (a * a');
            end
        end
    end
    disp('wss-o1')
    for l = 1 : L
        for g = 1 : G
            for k = 1 : Kx
                x = abs(posi_sto((l-1)*G*Kx+(g-1)*Kx+k,1)+bs_posi(l,1)-bs_posi(2,1));
                y = posi_sto((l-1)*G*Kx+(g-1)*Kx+k,2)+bs_posi(l,2)-bs_posi(2,2);
                a = arrsterc2(x, y, M, Nh, Nv, BS_height);
                Wss(:,M*G*L+M*G*(l-1)+(g-1)*M+1:M*G*L+M*G*(l-1)+g*M) = Wss(:,M*G*L+M*G*(l-1)+(g-1)*M+1:M*G*L+M*G*(l-1)+g*M) + (a * a');
            end
        end
    end
    disp('wss-o2')
    for l = 1 : L
        for g = 1 : G
            for k = 1 : Kx
                x = abs(posi_sto((l-1)*G*Kx+(g-1)*Kx+k,1)+bs_posi(l,1)-bs_posi(3,1));
                y = posi_sto((l-1)*G*Kx+(g-1)*Kx+k,2)+bs_posi(l,2)-bs_posi(3,2);
                a = arrsterc3(x, y, M, Nh, Nv, BS_height);
                Wss(:,2*M*G*L+M*G*(l-1)+(g-1)*M+1:2*M*G*L+M*G*(l-1)+g*M) = Wss(:,2*M*G*L+M*G*(l-1)+(g-1)*M+1:2*M*G*L+M*G*(l-1)+g*M) + (a * a');
            end
        end
    end
    disp('wss-o3')
    Wss = Wss / Kx;
    Cjl = zeros(Kq*(G*L)^2, Kq);
    Xgpx = XgpF1 * sqrt(rho);
    for lc = 1 : L
        for gc = 1 : G
            for lcc = 1 : L
                for gcc = 1 : G
                    if lc==lcc && gc==gcc
                        Cjl((lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+(gcc-1)*Kq+1:(lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+gcc*Kq,:)...
                            =eye(Kq)*rho;
                    else
                        Cjl((lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+(gcc-1)*Kq+1:(lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+gcc*Kq,:)...
                            = Xgpx((gcc-1)*Kq+1:gcc*Kq,(lcc-1)*Lp+1:lcc*Lp)*Xgpx((gc-1)*Kq+1:gc*Kq,(lc-1)*Lp+1:lc*Lp)';
                    end
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%
    for l = 1 : L
        for g = 1 : G
            Wjjg = Wss(:,(l-1)*M*G*L+M*G*(l-1)+(g-1)*M+1:(l-1)*M*G*L+M*G*(l-1)+g*M);
            switch g
                case 1
                    gp = 2;
                case 2
                    gp = 1;
            end
            Wjjgp = Wss(:,(l-1)*M*G*L+M*G*(l-1)+(gp-1)*M+1:(l-1)*M*G*L+M*G*(l-1)+gp*M);
            [Ujj, Djj] = eig(Wjjg);
            Ujj = Ujj(:, 1:Kq);
            Djj = Djj(1:Kq, 1:Kq);
            Hbathatg = zeros(Kq, Kq);
            Sigmajksto = zeros(Kq, Kq*Kq);
            Sigmajksum = zeros(Kq, Kq);
            for k = 1 : Kq
                Sigsi = 1/rho * eye(Kq);
                for ll = 1 : L
                    for gg = 1 : G
                        Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                        Sigsi = Sigsi + Ujj' * Wjlt * Ujj * ...
                            norm(Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+(gg-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+gg*Kq,k))^2/rho^2;
                    end
                end
                Sigmajk = Djj * inv(Sigsi);
                Sigmajksto(:, (k-1)*Kq+1:k*Kq) = Sigmajk;
                Sigmajksum = Sigmajksum + Sigmajk;
            end
            Xij = rho * Djj * (Kq * eye(Kq) - Sigmajksum') + eye(Kq);
            Xij = Xij + rho * Kq * Ujj' * Wjjgp * Ujj;
            for ll = 1 : L
                if ll == l
                else
                    for gg = 1 : G
                        Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                        Xij = Xij + rho * Kq * Ujj' * Wjlt * Ujj;
                    end
                end
            end
            Yj = zeros(Kq, Kq);
            for k = 1 : Kq
                Yj(k,k) = rho * trace(Djj * Sigmajksto(:, (k-1)*Kq+1:k*Kq)'* inv(Xij)*Sigmajksto(:, (k-1)*Kq+1:k*Kq))...
                    + trace(Sigmajksto(:, (k-1)*Kq+1:k*Kq)'* inv(Xij)*Sigmajksto(:, (k-1)*Kq+1:k*Kq));
            end
            Rj = Yj;
            Xjl = ones(Kq, Kq);
            for k = 1 : Kq
                for kk = 1 : Kq
                    Xjl(k, kk) = 1 / rho * trace(Ujj *Sigmajksto(:, (k-1)*Kq+1:k*Kq)' * inv(Xij)*Sigmajksto(:, (kk-1)*Kq+1:kk*Kq)' * Ujj'*Wjjgp);
                end
            end
            Rj = Rj + Xjl .* (Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+(gp-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+gp*Kq,:)'...
                *Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+(gp-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+gp*Kq,:));
            for ll = 1 : L
                if ll==l
                else
                    for gg = 1 : G
                        Cjl1 = Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+(gg-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+gg*Kq,:);
                        Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                        Xjl = ones(Kq, Kq);
                        for k = 1 : Kq
                            for kk = 1 : Kq
                                Xjl(k, kk) = 1 / rho * trace(Ujj *Sigmajksto(:, (k-1)*Kq+1:k*Kq)' ...
                                    * inv(Xij)*Sigmajksto(:, (kk-1)*Kq+1:kk*Kq)' * Ujj'*Wjlt);
                            end
                        end
                        Rj = Rj + Xjl .* (Cjl1'* Cjl1);
                    end
                end
            end
            Rj = real(Rj);
            inR = inv(Rj);
            for k = 1 : Kq
                Rxx = [Rj(1:k-1,:);Rj(k+1:Kq,:)];
                Ryy = [Rxx(:,1:k-1),Rxx(:,k+1:Kq)];
                bound(ant_t, 2) = bound(ant_t, 2) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/Kq/inR(k,k));
            end
            disp('first')
            disp([ant_t,n,l,g])
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    for l = 1 : L
        for g = 1 : G
            Wjjg = Wss(:,(l-1)*M*G*L+M*G*(l-1)+(g-1)*M+1:(l-1)*M*G*L+M*G*(l-1)+g*M);
            switch g
                case 1
                    gp = 2;
                case 2
                    gp = 1;
            end
            Wjjgp = Wss(:,(l-1)*M*G*L+M*G*(l-1)+(gp-1)*M+1:(l-1)*M*G*L+M*G*(l-1)+gp*M);
            [Ujj, Djj] = eig(Wjjg);
            Ujj = Ujj(:, 1:Kq);
            Djj = Djj(1:Kq, 1:Kq);
            for n = 1 : Nite
                Hqjjg = (randn(M, Kq) + 1i* randn(M, Kq)) / sqrt(2);
                Hjjg = Wjjg^0.5 * Hqjjg;
                Hbarjjg = Djj^0.5 * Ujj' * Hqjjg;
                Hqjjgp = (randn(M, Kq) + 1i* randn(M, Kq)) / sqrt(2);
                Hjjgp = Wjjgp^0.5 * Hqjjgp;
                Hbarjgpq = Ujj' * Hjjgp;
                Hbarjlt = (randn(M, K*2) + 1i* randn(M, K*2)) / sqrt(2);
                Hbarjl = zeros(Kq, K*2);
                lc = 0;
                for ll = 1 : L
                    if l == ll
                    else
                        lc = lc + 1;
                        for gg = 1 : G
                            Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                            Hbarjl(:, (lc-1)*K+(gg-1)*Kq+1:(lc-1)*K+gg*Kq) = Ujj' * Wjlt^0.5 * Hbarjlt(:, (lc-1)*K+(gg-1)*Kq+1:(lc-1)*K+gg*Kq);
                        end
                    end
                end
                Zp = (randn(M, Lp) + 1i* randn(M, Lp)) / sqrt(2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for pic = 1 : 3
                    switch pic
                        case 1
                            Xgpx = Xgpp * sqrt(rho);
                        case 2
                            Xgpx = Xgpmin * sqrt(rho);
                        case 3
                            Xgpx = XgpF1 * sqrt(rho);
                    end
                    for lc = 1 : L
                        for gc = 1 : G
                            for lcc = 1 : L
                                for gcc = 1 : G
                                    Cjl((lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+(gcc-1)*Kq+1:(lc-1)*G*G*L*Kq+(gc-1)*G*L*Kq+(lcc-1)*G*Kq+gcc*Kq,:)...
                                        = Xgpx((gcc-1)*Kq+1:gcc*Kq,(lcc-1)*Lp+1:lcc*Lp)*Xgpx((gc-1)*Kq+1:gc*Kq,(lc-1)*Lp+1:lc*Lp)';
                                end
                            end
                        end
                    end
                    Hbathatg = zeros(Kq, Kq);
                    Sigmajksto = zeros(Kq, Kq*Kq);
                    Sigmajksum = zeros(Kq, Kq);
                    for k = 1 : Kq
                        Sigsi = 1/rho * eye(Kq);
                        for ll = 1 : L
                            for gg = 1 : G
                                Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                                Sigsi = Sigsi + Ujj' * Wjlt * Ujj * ...
                                    norm(Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+(gg-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+gg*Kq,k))^2/rho^2;
                            end
                        end
                        Sigmajk = Djj * inv(Sigsi);
                        Sigmajksto(:, (k-1)*Kq+1:k*Kq) = Sigmajk;
                        Sigmajksum = Sigmajksum + Sigmajk;
                        Ygp = rho*Hbarjjg + Ujj' * Zp * Xgpx((g-1)*Kq+1:g*Kq,(l-1)*Lp+1:l*Lp)';
                        Ygp = Ygp + Hbarjgpq*...
                            Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+(gp-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(l-1)*G*Kq+gp*Kq,:);
                        lc = 0;
                        for ll = 1 : L
                            if l == ll
                            else
                                lc = lc + 1;
                                for gg = 1 : G
                                    Ygp = Ygp + Hbarjl(:, (lc-1)*K+(gg-1)*Kq+1:(lc-1)*K+gg*Kq)*...
                                        Cjl((l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+(gg-1)*Kq+1:(l-1)*G*G*L*Kq+(g-1)*G*L*Kq+(ll-1)*G*Kq+gg*Kq,:);
                                end
                            end
                        end
                        Hbathatg(:, k) = Sigmajk * Ygp(:, k)/rho;
                    end
                    Xij = rho * Djj * (Kq * eye(Kq) - Sigmajksum') + eye(Kq);
                    Xij = Xij + rho * Kq * Ujj' * Wjjgp * Ujj;
                    for ll = 1 : L
                        if ll == l
                        else
                            for gg = 1 : G
                                Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                                Xij = Xij + rho * Kq * Ujj' * Wjlt * Ujj;
                            end
                        end
                    end
                    Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
                    Vxxj = inv(eye(Kq)+Vj);
                    for m = 1 : Kq
                        sumrate(ant_t, pic) = sumrate(ant_t, pic) + log2(1/real(Vxxj(m,m)));%%reused
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    if 1 == pic
                        Hbathatg = zeros(Kq, Kq);
                        Sigmajksto = zeros(Kq, Kq*Kq);
                        Sigmajksum = zeros(Kq, Kq);
                        for k = 1 : Kq
                            Sigsi = 1/rho * eye(Kq) + Djj;
                            Sigmajk = Djj * inv(Sigsi);
                            Sigmajksto(:, (k-1)*Kq+1:k*Kq) = Sigmajk;
                            Sigmajksum = Sigmajksum + Sigmajk;
                            Ygp = rho*Hbarjjg + Ujj' * Zp * Xgpx((g-1)*Kq+1:g*Kq,(l-1)*Lp+1:l*Lp)';
                            Hbathatg(:, k) = Sigmajk * Ygp(:, k)/rho;
                        end
                        Xij = rho * Djj * (Kq * eye(Kq) - Sigmajksum') + eye(Kq);
                        Xij = Xij + rho * Kq * Ujj' * Wjjgp * Ujj;
                        for ll = 1 : L
                            if ll == l
                            else
                                for gg = 1 : G
                                    Wjlt = Wss(:,(l-1)*M*G*L+M*G*(ll-1)+(gg-1)*M+1:(l-1)*M*G*L+M*G*(ll-1)+gg*M);
                                    Xij = Xij + rho * Kq * Ujj' * Wjlt * Ujj;
                                end
                            end
                        end
                        Vj = rho * Hbathatg' * inv(Xij) * Hbathatg;
                        Vxxj = inv(eye(Kq)+Vj);
                        for m = 1 : Kq
                            sumrate(ant_t, 4) = sumrate(ant_t, 4) + log2(1/real(Vxxj(m,m)));%%reused
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    disp('cal')
                    disp([ant_t,l,g,n,pic])
                end
                if length(ant_s) == ant_t
                    Hbathatg = Djj / (Djj+1/rho*eye(Kq)) * (Hbarjjg+1/rho*Ujj' * Zp * Xgpx((g-1)*Kq+1:g*Kq,(l-1)*Lp+1:l*Lp)');
                    Xij = rho * Kq * Djj / (rho * Djj+ eye(Kq)) + eye(Kq);
                    Vxxj =  inv(eye(Kq)+rho *Hbathatg'*inv(Xij)*Hbathatg);
                    for m = 1 : Kq
                        bound(:, 1) = bound(:, 1) + log2(1/real(Vxxj(m,m)));
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumrate = sumrate / Nite;
bound(:,1) = bound(:,1) / Nite;
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1=plot(ant_s, sumrate(:, 1), 'k--s','LineWidth',1,'MarkerSize',10);
hold on
plot(ant_s, sumrate(:, 2), 'k--x','LineWidth',1,'MarkerSize',10);
plot(ant_s, sumrate(:, 3), 'k-p','LineWidth',1,'MarkerSize',10);
plot(ant_s, sumrate(:, 4), 'k-^','LineWidth',1,'MarkerSize',10);
plot(ant_s, bound(:, 2), 'k-o','LineWidth',1,'MarkerSize',10);
plot(ant_s, bound(:, 1), 'k-+','LineWidth',1,'MarkerSize',10);
xlim([min(ant_s), max(ant_s)])
ylim([0,90])
le = legend('Reuse','Min-phase [6][7]','Proposed','Ideal','Upper bound (52)','Asymptotic (56)','Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',ant_s)
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('Number of vertical antennas','Fontname','Times','Fontsize',12)
ylabel('Sum rate (bps/Hz)','Fontname','Times','Fontsize',12)
grid on%%%%%%%%%%%%%%%%%%%%%%%%
