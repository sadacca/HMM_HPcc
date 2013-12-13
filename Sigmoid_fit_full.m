function sigmoid_output = Sigmoid_fit_full(alpha, xdata)

sigmoid_output = alpha(1) + alpha(2) ./ ( 1 + exp( -( xdata - alpha(3) ) / alpha(4) ) ) ;