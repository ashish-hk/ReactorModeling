module Molecular

export Species

mutable struct Species
    Name::String
    Molecular_weight::Float64
    Phase::Char
    High_temperature_coeff::Array{Float64,1}
    Low_temperature_coeff::Array{Float64,1}
    Species() = new()
    Species(Name::String) = new(Name)
end

f1 = open("therm.dat")
data=readlines(f1)

R=8.314

# function1
function mol_sp_heat(T,array)
    sp_heat=Dict()
    for ele in array
        for i in range(2, stop=length(data)-1, step=4)
            if strip(data[i][1:18]) == ele
                if T >= parse(Float64,strip(data[i][66:73])) # Upper temp Range
                    a1 = parse(Float64, strip(data[i+1][1:15]))
                    a2 = parse(Float64, strip(data[i+1][16:30]))
                    a3 = parse(Float64, strip(data[i+1][31:45]))
                    a4 = parse(Float64, strip(data[i+1][46:60]))
                    a5 = parse(Float64, strip(data[i+1][61:75]))
                    a6 = parse(Float64, strip(data[i+2][1:15]))
                    sp_heat[ele] = R*(a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4))
                                        
                else  # Lower temp Range
                    a1 = parse(Float64, strip(data[i+2][31:45]))
                    a2 = parse(Float64, strip(data[i+2][46:60]))
                    a3 = parse(Float64, strip(data[i+2][61:75]))
                    a4 = parse(Float64, strip(data[i+3][1:15]))
                    a5 = parse(Float64, strip(data[i+3][16:30]))
                    sp_heat[ele] = R*(a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4))
                end
            end
        end
    end
    return sp_heat 
end


# function2
function mol_enthaply(T,array)
    enthalpy=Dict()
    R=8.314
    for ele in array
        for i in range(2, stop=length(data)-1, step=4)
            if strip(data[i][1:18]) == ele
                if T >= parse(Float64,strip(data[i][66:73])) # Upper temp Range
                    a1 = parse(Float64, strip(data[i+1][1:15]))
                    a2 = parse(Float64, strip(data[i+1][16:30]))
                    a3 = parse(Float64, strip(data[i+1][31:45]))
                    a4 = parse(Float64, strip(data[i+1][46:60]))
                    a5 = parse(Float64, strip(data[i+1][61:75]))
                    a6 = parse(Float64, strip(data[i+2][1:15]))
                    enthalpy[ele] = R*T*(a1 + a2*T/2 + a3*(T^2/3) + a4*(T^3/4) + a5*(T^4/5) + a6/T)/1000  # In KJ
                                        
                else  # Lower temp Range
                    a1 = parse(Float64, strip(data[i+2][31:45]))
                    a2 = parse(Float64, strip(data[i+2][46:60]))
                    a3 = parse(Float64, strip(data[i+2][61:75]))
                    a4 = parse(Float64, strip(data[i+3][1:15]))
                    a5 = parse(Float64, strip(data[i+3][16:30]))
                    a6 = parse(Float64, strip(data[i+3][31:45]))
                    enthalpy[ele] = R*T*(a1 + a2*T/2 + a3*(T^2/3) + a4*(T^3/4) + a5*(T^4/5) + a6/T)/1000  # In KJ
                end

            end
        end
    end
    return enthalpy   # in KJ/mol
end

# function3
function mol_entropy(T,array)
    entropy=Dict()
    R=8.314
    for ele in array
        for i in range(2, stop=length(data)-1, step=4)
            if strip(data[i][1:18]) == ele
                    
                if T >= parse(Float64,strip(data[i][66:73]))
                    a1 = parse(Float64, strip(data[i+1][1:15]))
                    a2 = parse(Float64, strip(data[i+1][16:30]))
                    a3 = parse(Float64, strip(data[i+1][31:45]))
                    a4 = parse(Float64, strip(data[i+1][46:60]))
                    a5 = parse(Float64, strip(data[i+1][61:75]))
                    a7 = parse(Float64, strip(data[i+2][16:30]))
                    entropy[ele] = R*(a1*log(T) +a2*T+ a3*(T^2/2) + a4*(T^3/3) + a5*(T^4/4) + a7 )                  
                else
                    a1 = parse(Float64, strip(data[i+2][31:45]))
                    a2 = parse(Float64, strip(data[i+2][46:60]))
                    a3 = parse(Float64, strip(data[i+2][61:75]))
                    a4 = parse(Float64, strip(data[i+3][1:15]))
                    a5 = parse(Float64, strip(data[i+3][16:30]))
                    a7 = parse(Float64, strip(data[i+3][46:60]))
                    entropy[ele] = R*(a1*log(T) +a2*T+ a3*(T^2/2) + a4*(T^3/3) + a5*(T^4/4) + a7 )
                end    
            end
        end
    end
    return entropy   # in J/mol.K
end



# function4
function mixture_enthalpy(T, mol_fracs,species_objects)
    mol_frac_dict=Dict()
    
    enthaplies = mol_enthaply(T,species_objects)
    mixture_enthalpy=0
    
    for i =1:length(mol_fracs)
            mol_frac_dict[species_objects[i]] = mol_fracs[i]
    end

    
    for ele in species_objects
        mixture_enthalpy += enthaplies[ele].*mol_frac_dict[ele]
    end
    return mixture_enthalpy

end


# function5
function mixture_specifc_heat(T, mol_fracs,species_objects)
    mol_frac_dict=Dict()
    mol_sp_heats= mol_sp_heat(T,species_objects)
    mixture_specifc_heat=0
    for i =1:length(species_objects)
            mol_frac_dict[species_objects[i]] = mol_fracs[i]
    end
        
    for ele in species_objects
        mixture_specifc_heat += mol_sp_heats[ele].*mol_frac_dict[ele]
    end
    return mixture_specifc_heat
end

end