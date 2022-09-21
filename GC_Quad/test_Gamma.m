function test_Gamma

Gamma_list = linspace(0.4,0.6,100);
total_count_list = zeros(size(Gamma_list));
for i = 1:length(Gamma_list)
    total_count_list(i) = GC_Quad(Gamma_list(i));
end

plot(Gamma_list,total_count_list)
print('TotalNumofFun_vs_Gamma','-depsc2')
end