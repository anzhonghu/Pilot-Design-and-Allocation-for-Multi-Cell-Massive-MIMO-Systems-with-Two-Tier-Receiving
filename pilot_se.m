ant_it =  length(ant_s);
Nv = ant_s(ant_it, 1);
M = Nh * Nv;
Kx = M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kc = 0;
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
        posi_sto(Kx+Kc,:) = [x,y];
    else
    end
end
Kc = 0;
while Kc<Kx
    x = (rand(1,1)-0.5) * ISD;
    y = (rand(1,1)-0.5) * ISD;
    if x^2+y^2 <= (0.5*ISD)^2 && x^2+y^2>=Min_d^2 && y>=0 &&  (x>=0 || (x<=0 && atan(abs(y/x))>=pi/3))%cell 3
        Kc = Kc + 1;
        posi_sto(2*Kx+Kc,:) = [x,y];
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W11 = zeros(M, M);
W12 = zeros(M, M);
W13 = zeros(M, M);
W22 = zeros(M, M);
W21 = zeros(M, M);
W23 = zeros(M, M);
W33 = zeros(M, M);
W31 = zeros(M, M);
W32 = zeros(M, M);
a = zeros(M, 1);
for l = 1 : L
    for k = 1 : Kx
        x = posi_sto((l-1)*M+k,1)+bs_posi(l,1);
        y = posi_sto((l-1)*M+k,2)+bs_posi(l,2);
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
    for k = 1 : Kx
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
    for k = 1 : Kx
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
W11 = W11  / Kx;
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
        end
        disp('first-pilot')
        disp([ant_it,cellr,n])
    end
    [maxupbur,nq] = max(upbur);
    XgpF(:,(jxx-1)*Lp+1:jxx*Lp) = Fn(:, (nq-1)*Lp+1:nq*Lp);
end