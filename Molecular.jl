module Molecular

export Species

mutable struct Species
    Name::String
    Molecular_weight::Float64
    Phase::Char
    High_temperature_coeff::Array{Float64,1}
    Low_temperature_coeff::Array{Float64,1}
    Common_Temp ::Float64
    Species() = new()
    Species(Name::String) = new(Name)
end

# f1 = open("therm.dat")
# data=readlines(f1)

R=8.314




# function1
function mol_sp_heat(T,array)
    sp_heat=Dict()
    R=8.314
    for ele in array
        if T >= ele.Common_Temp        # Upper temp Range
            a1 = ele.High_temperature_coeff[1]
            a2 = ele.High_temperature_coeff[2]
            a3 = ele.High_temperature_coeff[3]
            a4 = ele.High_temperature_coeff[4]
            a5 = ele.High_temperature_coeff[5]
            sp_heat[ele.Name] = R*(a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4))

        else  # Lower temp Range
            a1 = ele.Low_temperature_coeff[1]
            a2 = ele.Low_temperature_coeff[2]
            a3 = ele.Low_temperature_coeff[3]
            a4 = ele.Low_temperature_coeff[4]
            a5 = ele.Low_temperature_coeff[5]
            sp_heat[ele.Name] = R*(a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4))
        end
end
    return sp_heat 
end


# function2
function mol_enthaply(T,array)
    enthalpy=Dict()
    R=8.314
    for ele in array
        if T >= ele.Common_Temp         # Upper temp Range
            a1 = ele.High_temperature_coeff[1]
            a2 = ele.High_temperature_coeff[2]
            a3 = ele.High_temperature_coeff[3]
            a4 = ele.High_temperature_coeff[4]
            a5 = ele.High_temperature_coeff[5]
            a6 = ele.High_temperature_coeff[6]

            enthalpy[ele.Name] = R*T*(a1 + a2*T/2 + a3*(T^2/3) + a4*(T^3/4) + a5*(T^4/5) + a6/T)/1000  # In KJ

        else  # Lower temp Range
            
            a1 = ele.Low_temperature_coeff[1]
            a2 = ele.Low_temperature_coeff[2]
            a3 = ele.Low_temperature_coeff[3]
            a4 = ele.Low_temperature_coeff[4]
            a5 = ele.Low_temperature_coeff[5]
            a6 = ele.Low_temperature_coeff[6]

            enthalpy[ele.Name] = R*T*(a1 + a2*T/2 + a3*(T^2/3) + a4*(T^3/4) + a5*(T^4/5) + a6/T)/1000  # In KJ
        end
    end
    return enthalpy   # in KJ/mol
end





# function3
function mol_entropy(T,array)
    entropy=Dict()
    R=8.314
    for ele in array
        if T >= ele.Common_Temp         # Upper temp Range
            a1 = ele.High_temperature_coeff[1]
            a2 = ele.High_temperature_coeff[2]
            a3 = ele.High_temperature_coeff[3]
            a4 = ele.High_temperature_coeff[4]
            a5 = ele.High_temperature_coeff[5]
            a7 = ele.High_temperature_coeff[7]
            entropy[ele.Name] = R*(a1*log(T) +a2*T+ a3*(T^2/2) + a4*(T^3/3) + a5*(T^4/4) + a7 )                  

        else  # Lower temp Range
            
            a1 = ele.Low_temperature_coeff[1]
            a2 = ele.Low_temperature_coeff[2]
            a3 = ele.Low_temperature_coeff[3]
            a4 = ele.Low_temperature_coeff[4]
            a5 = ele.Low_temperature_coeff[5]
            a7 = ele.Low_temperature_coeff[7]
            
            entropy[ele.Name] = R*(a1*log(T) +a2*T+ a3*(T^2/2) + a4*(T^3/3) + a5*(T^4/4) + a7 )                  
        end
    end
    return entropy  
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
        mixture_enthalpy += enthaplies[ele.Name].*mol_frac_dict[ele]
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
        mixture_specifc_heat += mol_sp_heats[ele.Name].*mol_frac_dict[ele]
    end
    return mixture_specifc_heat
end

end