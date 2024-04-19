function memory_check(c1, c2)
    % memory_check.m
    %
    % Checks memory usage, displaying a warning below a certain level, and
    % exitting if the memory usage is too high.
    %
    % Inputs:
    %   c1: The percentage of total memory used after which to display a warning.
    %   c2: The percentage of total memory used after which to throw an error.
    %
    % Author: Aidan Gerkis
    % Date: 15-04-2024
    
    switch nargin
        case 0 % No arguments needed
            c1 = 0.8; % Percent of memory after which to display a warning
            c2 = 0.95; % Percent of memory after which to give an error
        case 1 % Only c1 passed
            c2 = 0.95;
        case 2 % Both arguments passed, does nothing
        otherwise
            error("Error in function call, too many input arguments");
    end

    m = memory; % Call MATLABs default memory evaluation code
    max_bytes = m.MaxPossibleArrayBytes; % Total available memory
    mem_used = m.MemUsedMATLAB; % Memory used

    if mem_used > max_bytes*c1
        warning("IMOPT memory usage high, performance may be poor.");
    elseif mem_used >max_bytes*c2
        error("IMOPT memory usage exceeded" + num2str(max_bytes*c2) + "Bytes, exitting program.");
    end
end