function SOI_Mjday = calculate_SOI(year_CA,month_CA_num,date_CA,PHA_bsp)
    %% Calculate SOI date
    for PHA = 33 : 33
        % CA Date
        [years1,months1, days1, fds1] = iauJd2cal( 2400000.5, Mjday(year_CA, month_CA_num, date_CA));
        [hours1, minutes1, secs1] = fd_to_hms(fds1);
        [years1, months1, days1, hours1, minutes1, secs1] = fix_seconds(years1, months1, days1, hours1, minutes1, secs1);
        secs1 = round(secs1,3);
        months1_string = monthToString(months1);
        str = sprintf(' %s %g , %g %g:%g:%g', months1_string, days1, years1, hours1, minutes1, secs1);
        et2 = cspice_str2et(str);
        
        % one year previous
        [years1,months1, days1, fds1] = iauJd2cal( 2400000.5, Mjday(year_CA-1, month_CA_num, date_CA));
        [hours1, minutes1, secs1] = fd_to_hms(fds1);
        [years1, months1, days1, hours1, minutes1, secs1] = fix_seconds(years1, months1, days1, hours1, minutes1, secs1);
        secs1 = round(secs1,3);
        months1_string = monthToString(months1);
        str = sprintf(' %s %g , %g %g:%g:%g', months1_string, days1, years1, hours1, minutes1, secs1);
        et1 = cspice_str2et(str);
        
        step_size = 3600; % s
        
        [Y_PHA, ~] = cspice_spkezr(char(PHA_bsp),  et2 : -step_size : et1 , 'J2000', 'NONE', 'SUN');
        
        
        [Y_Earth, ~] = cspice_spkezr('Earth', et2 : -step_size : et1 , 'J2000', 'NONE', 'SUN');
        
        relative_distance = zeros(length(Y_Earth),1);
        for i = 1 : length(Y_Earth)
            relative_distance(i) = norm( Y_Earth(1:3,i) - Y_PHA(1:3,i) ) ;  % km
        end
        
        SOI_Mjday = Mjday(year_CA, month_CA_num, date_CA) - ( ( find(relative_distance > 2.6e6,1) - 1 ) * step_size ) / 86400 - 1; % minus one more extra day 
    
    end
    
end