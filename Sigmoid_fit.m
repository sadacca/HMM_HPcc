function sigmoid_output = Sigmoid_fit(alpha, xdata)

sigmoid_output = alpha(2) + 1 ./ (1 + exp(xdata * alpha(1))) ;

