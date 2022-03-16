Answer_Vs_Expected = readtable("Answer_Vs_Expected");
Variance = readtable("Variance");
RR_Variance = readtable("RR_Variance");


x = linspace(10, 1000, 100);
y = 0.0075./x;
plot_variance(Variance.No_Samples, Variance.Variance, x, y);

y = 0.03./x;
plot_variance(RR_Variance.No_Samples, RR_Variance.Variance, x, y);

figure
plot(Answer_Vs_Expected.No_Bounces, Answer_Vs_Expected.Expected_Answer, Answer_Vs_Expected.No_Bounces, Answer_Vs_Expected.Actual_Answer);