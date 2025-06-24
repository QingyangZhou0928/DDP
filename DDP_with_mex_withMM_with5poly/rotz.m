function T = rotz(alpha)
    T = [cos(alpha), -sin(alpha), 0, 0;
        sin(alpha), cos(alpha), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];
end