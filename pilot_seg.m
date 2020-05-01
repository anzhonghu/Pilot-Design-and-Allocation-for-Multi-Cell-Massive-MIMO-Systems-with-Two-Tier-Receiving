ant_t = length(ant_s);
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
%%%%%%%%%%%%%%%%%%%%%%
for cellr = 1 : L
    upbur = zeros(Npd, 1);
    for n = 1 : Npd
        XgpF1(:,(cellr-1)*Lp+1:cellr*Lp) = Fn1(:, (n-1)*Lp+1:n*Lp);
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
                    upbur(n,1) = upbur(n,1) + log2(1+Rj(k,k)-trace(Ryy)/ (1 + trace(Ryy))/Kq/inR(k,k));
                end
                disp('first')
                disp([ant_t, cellr,n,l,g])
            end
        end
    end
    [upburmaxv,nq] = max(upbur);
    XgpF1(:,(cellr-1)*Lp+1:cellr*Lp) = Fn1(:, (nq-1)*Lp+1:nq*Lp);
end