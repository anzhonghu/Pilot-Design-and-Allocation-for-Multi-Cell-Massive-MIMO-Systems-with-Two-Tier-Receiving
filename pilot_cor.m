Lp = 10;
N = 100;
F = zeros(Lp, Lp*N);
for n = 1 : N
    an = randn;
    for lp = 1 : Lp
        for lpp = 1 : Lp
            F(lp, (n-1)*Lp+lpp) = 1/sqrt(Lp)*exp(1i*2*pi/Lp*(lpp+an)*lp);
        end
    end
end
cor = abs(F(:,1)'*F(:,Lp+1:N*Lp));
h1=histogram(cor,'FaceColor','k');
set(gca,'FontName','Times New Roman','FontSize',12)
h1.Normalization = 'probability';
xlabel('Correlation value','Fontname','Times','Fontsize',15)
ylabel('Proportion','Fontname','Times','Fontsize',15)