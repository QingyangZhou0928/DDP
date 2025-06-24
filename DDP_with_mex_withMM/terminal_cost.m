function terminalcost = terminal_cost(x, xgoal, QN)
    terminalcost = 0.5*(x-xgoal)'*QN*(x-xgoal);
end