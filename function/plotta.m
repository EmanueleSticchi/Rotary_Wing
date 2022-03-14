function plotta(x,y,label,formatspec,nome)
% label = {labelx ; labely};
set(0,'DefaultAxesFontName', 'Times New Roman');
if nargin ==3
    plot(x,y,'k')
elseif nargin==5
    plot(x,y,formatspec,'DisplayName',nome)
else
    plot(x,y,formatspec)
end
xlabel(label(1),'Interpreter','latex','FontSize',18)
ylabel(label(2),'Interpreter','latex','FontSize',18)
hold on
grid on
end