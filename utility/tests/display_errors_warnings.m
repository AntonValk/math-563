function display_errors_warnings(err, warn)
    % Extract field names
    err_fields = fieldnames(err);
    warn_fields = fieldnames(warn);

    % Display all errors
    for i=1:length(err_fields)
        disp("  " + err.(err_fields{i}));
    end

    % Display all warnings
    % Display all errors
    for i=1:length(warn_fields)
        disp("  " + warn.(warn_fields{i}));
    end
end