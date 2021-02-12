for space = 1:LDT
    A(space) = ((R * T * unDelP) / (P_atm * (1 - (A2/A1)^2) * sqrt ((2 * DelP(space) * R * T) / (P_atm * (1 - (A2/A1)^2)))))

    B(space) = (((-(1 - (A2/A1)^2)) * DelP(space) * R * T * unP_atm) / ((P_atm * (1 - (A2/A1)^2))^2 * sqrt ((2 * DelP(space) * R * T) / (P_atm * (1 - (A2/A1)^2)))))

    C(space) = ((DelP(space) * R * unT) / (P_atm * (1 - (A2/A1)^2) * sqrt ((2 * DelP(space) * R * T) / (P_atm * (1 - (A2/A1)^2)))))

    Venturi_Error(space) = sqrt((A^2)+(B^2)+(C^2))
end

