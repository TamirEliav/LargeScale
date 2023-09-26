function rounded = roundToNearestValue(x, targetValue, rounding_func)
arguments
    x
    targetValue
    rounding_func = @round
end
    if targetValue == 0
        error('Target value cannot be zero.');
    end
    
    rounded = rounding_func(x / targetValue) * targetValue;
end
