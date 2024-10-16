"""
Main module for PetroBase, provides basic types and functions to solve petrological problems.

# Exports
$(EXPORTS)
"""
module PetroBase

export 
    Chemical, 
    Component, 
    TraceElement, 
    Phase, 
    PetroSystem,
    ≃, 
    name, 
    concentration, 
    sum_mass, 
    isunique,
    findchemical,
    gibbs,
    mol,
    sum_mols,
    majorcation

using
    DocStringExtensions
# Write your package code here.

const O_MASS = 15.999
const H_MASS = 1.00784
"""
Abstract supertype for chemical components and trace elements
"""
abstract type Chemical end

"""
$(SIGNATURES)
Basic unit in describing composition of a phase or system. Critical to most calculations. 
Typically the only values that will ever changed are mol and chemical potential.

$(TYPEDFIELDS)
"""
struct Component <: Chemical
    "Name of the component"
    name::String
    "Molar mass (g/mol)"
    molarmass::Float64 
    "Moles of the component"
    mol::Float64
    "Number of cations"
    cat::Real
    "Number of oxygen atoms"
    ox::Real
    "Charge on cations"
    catcharge::Real
    "Chemical potential (J/mol)"
    μ::Float64 #Chemical potential in J/mol
end

Component(name,molarmass,mol,cat,ox,catcharge) = Component(name,molarmass,mol,cat,ox,catcharge,0.0)
Component(name,molarmass,mol) = Component(name,molarmass,mol,0,0,0,0.0)

Component(name,molarmass,cat,ox,catcharge;mass=0,μ=0) = Component(name,molarmass,mass/molarmass,cat,ox,catcharge,μ)
"""
$(SIGNATURES)
Clones the parameters of a 'Component' but with change of 'mol' and/or 'μ'
"""
# Component(clone, mol = clone.mol;μ=clone.μ) = Component(clone.name, clone.molarmass, mol,clone.cat,clone.ox,clone.catcharge, μ)
# Component(clone;mass = 0, μ=clone.μ) = Component(clone.name, clone.molarmass, mass/clone.molarmass,clone.cat,clone.ox,clone.catcharge, μ)

function Component(clone; mol = Nothing,  μ = clone.μ, mass = Nothing)
    if mass isa Number && mol isa Number
        throw(ArgumentError("Cannot define component with both mol and mass, please choose one"))
    end

    if mass isa Number
        return Component(clone.name, clone.molarmass, mass/clone.molarmass,clone.cat,clone.ox,clone.catcharge, μ)
    elseif mol isa Number
        return Component(clone.name, clone.molarmass, mol,clone.cat,clone.ox,clone.catcharge, μ)
    else
        return Component(clone.name,clone.molarmass,clone.mol,clone.cat,clone.ox,clone.catcharge,μ)
    end
end

"""
$(SIGNATURES)
Describes the trace element concentration of a phase or system. Typically the only value that will change is concentration

$(TYPEDFIELDS)
"""
struct TraceElement <: Chemical
    "Name of the element"
    name::String
    "Molar mass (g/mol)"
    molarmass::Float64
    "Concentration in the system or phase (µg/g)"
    concentration::Float64
    
end

"""
$(SIGNATURES)
Clones the parameters of a 'TraceElement' but with change of 'concentration'
"""
TraceElement(clone,concentration) = TraceElement(clone.name,clone.molarmass,concentration)

"""
$(TYPEDSIGNATURES)

Returns 'name' of any PetroBase struct with a 'name' parameter, including 'Component' and 'Phase' useful for broadcasting
"""
function name(item) 
    return item.name
end

"""
$(TYPEDSIGNATURES)

Returns the 'mol' of 'comp', useful for broadcasting
"""
function concentration(comp::Component)
    return comp.mol
end

"""
$(TYPEDSIGNATURES)

Returns the 'concentration' of a 'TraceElement', useful for broadcasting
"""
function concentration(te::TraceElement)
    return te.concentration
end
#An operator defined to compare if two components are the same except for the amount of moles
"""
$(TYPEDSIGNATURES)

This is a simple boolean operator that returns true if two 'Chemical' variables have the same name and molar mass
"""
function ≃(comp1::Chemical, comp2::Chemical)
    if comp1.name == comp2.name && comp1.molarmass ≈ comp2.molarmass
        return true
    end
    
    return false
end

"""
$(TYPEDSIGNATURES)

This is a simple boolean operator that returns true if two 'Component' variables have the same name. molar mass, and amount of moles
"""
function Base.:≈(comp1::Component, comp2::Component)

    if comp1 ≃ comp2 && comp1.mol ≈ comp2.mol
        return true
    end
    
    return false
end

"""
$(TYPEDSIGNATURES)

This is a simple boolean operator that returns true if two 'TraceElement' variables have the same name. molar mass, and amount of moles
"""
function Base.:≈(te1::TraceElement, te2::TraceElement)

    if te1 ≃ te2 && te1.mol ≈ te2.mol
        return true
    end
    
    return false
end
#Basic operations for QoL addition of components,
#currently have it return an array if components arent compatible but it might make more sense
#to just throw an error
"""
$(TYPEDSIGNATURES)
    
Adds together the 'mol' parameter of two 'Component' variables that have the same name and molar mass.
This will return a 'Component' with the chemical potential of the first argument.
"""
function Base.:+(comp1::Component, comp2::Component)
    if comp1 ≃ comp2
        return Component(comp1,mol = comp1.mol + comp2.mol)
    else
        throw(ArgumentError("comp1 and comp2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Adds 'num' to the 'mol' of 'comp1'
"""
function Base.:+(comp1::Component, num::Real)
    return Component(comp1, mol = comp1.mol+num)
end

Base.:+(num::Real, comp1::Component) = Base.:+(comp1,num)

"""
$(TYPEDSIGNATURES)
    
Adds together the 'concentration' parameter of two 'TraceElement' variables that have the same name and molar mass
"""
function Base.:+(te1::TraceElement,te2::TraceElement)
    if te1 ≃ te2
        return TraceElement(te1,te1.concentration + te2.concentration)
    else
        throw(ArgumentError("te1 and te2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Adds 'num' to the 'concentration' paremeter of 'te1'
"""
function Base.:+(te1::TraceElement,num::Real)
    return TraceElement(te1,te1.concentration+num)
end

Base.:+(num::Real,te1::TraceElement) = Base.:+(te1,num)

"""
$(TYPEDSIGNATURES)

Subtracts the 'mol' parameter of 'comp2' from 'comp1' as long as they have the same name and molar mass.
This will return a 'Component' with the chemical potential of the first argument.
"""
function Base.:-(comp1::Component, comp2::Component)
    if comp1 ≃ comp2
        return Component(comp1,mol = comp1.mol - comp2.mol)
    else
        throw(ArgumentError("comp1 and comp2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Subtracts 'num' to the 'mol' parameter of 'comp1'
"""
function Base.:-(comp1::Component, num::Real)
    return Component(comp1, mol = comp1.mol-num)
end

"""
$(TYPEDSIGNATURES)

Subtracts the 'concentration' parameter of 'te2' from 'te1' as long as they have the same name and molar mass
"""
function Base.:-(te1::TraceElement,te2::TraceElement)
    if te1 ≃ te2
        return TraceElement(te1,te1.concentration - te2.concentration)
    else
        throw(ArgumentError("te1 and te2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Subtracts 'num' to the 'concentration' parameter of 'te1'
"""
function Base.:-(te1::TraceElement, num::Real)
    return TraceElement(te1, te1.concentration-num)
end

"""
$(TYPEDSIGNATURES)

Multiplies the 'mol' parameter of 'comp1' by 'num'
"""
function Base.:*(comp1::Component, num::Real)
    return Component(comp1, mol = comp1.mol*num)
end

Base.:*(num::Real, comp1::Component) = Base.:*(comp1,num)
   
"""
$(TYPEDSIGNATURES)

Multiplies the 'concentration' parameter of 'te1' by 'num'
"""
function Base.:*(te1::TraceElement, num::Real)
    return TraceElement(te1, te1.concentration*num)
end

Base.:*(num::Real, te1::TraceElement) = Base.:*(te1,num)

"""
$(TYPEDSIGNATURES)

Divides the 'mol' parameter of 'comp1' by 'num'
"""
function Base.:/(comp1::Component, num::Real)
    return Component(comp1, mol = comp1.mol/num)
end

"""
$(TYPEDSIGNATURES)

Divides the 'concentration' parameter of 'te1' by 'num'
"""
function Base.:/(te1::TraceElement, num::Real)
    return TraceElement(te1, te1.concentration/num)
end

"""
$(TYPEDSIGNATURES)
Calculates the molar mass of an array of 'Component' variables by adding up the product of the molar mass and mol of each 'Component' in the array.
Intent of use is calculating molar mass of a phase that is described using a 'Component' array
"""
function sum_mass(components)
    molarmass = 0
    for comp in components
        molarmass += comp.molarmass*comp.mol
    end
    return molarmass
end

function sum_mols(components)
    mol = 0

    for comp in components
        mol += concentration(comp)
    end

    return mol
end
"""
$(TYPEDSIGNATURES)
Checks if each cell in a 'Chemical' array is not repeated elsewhere in the array.
"""
function isunique(chemicals)
    for i in 1:(lastindex(chemicals)-1), j in (i+1):lastindex(chemicals) #Does not check against element that already checked the whole array
        if chemicals[i] ≃ chemicals[j]
            return false
        end
    end
    return true
end


"""
$(TYPEDSIGNATURES)
Finds the first index of 'fChem' in the 'chem' array. Returns 0 if 'fChem' isnt present. Best used with arrays of unique 'Chemical' variables.
"""
function findchemical(chemicals,fchem::Chemical)
    for i in 1:lastindex(chemicals)
        if chemicals[i] ≃ fchem
            return i
        end
    end
    return 0
end


"""
$(TYPEDSIGNATURES)
Finds the first index of an element in the 'chem' array with the name of 'fChem'. Returns 0 if 'fChem' isnt present. 
Best used with arrays of unique 'Chemical' variables.
"""
function findchemical(chemicals, fchem)
    for i in 1:lastindex(chemicals)
        if chemicals[i].name == fchem
            return i
        end
    end
    return 0
end

"""
$(SIGNATURES)
This type is defined to contain all the most relevant properties of a phase, any number of variables can be initialized and defaults will be 0
or empty arrays/strings. Units are selected based on convenience for petrological modelling.
$(TYPEDFIELDS)
"""
@kwdef struct Phase
    "Name of the phase"
    name::String
    "Composition of the phase as defined by a 'Component' array"
    composition::Array{Component}
    "Concentration of trace elements in the phase"
    traceelements::Array{TraceElement} = Array{TraceElement}([])#In general I anticipate it is more useful to have TEs and Components seperate
    "Number of moles of the phase present in the system"
    mol::Float64 = 0
    "Volume of phase present in the system (J/bar [m^3/10^5])"
    vol::Float64 = 0#In J/bar
    "Total mass of the phase present in the system (g)"
    mass::Float64 = 0 #In g
    "Density (kg/m^3)"
    ρ::Float64 = 0 #Density in kg/m^3
    "Molar mass (g/mol)"
    molarmass::Float64 = 0
    #Thermo properties
    "Gibbs free energy (J/mol)"
    G::Float64 = 0 #J/mol
    "Enthalpy (J/mol)"
    H::Float64 = 0 #J/mol
    "Entropy (J/mol)"
    S::Float64 = 0#In J/Kmol
    "Heat capacity (J/K·mol)"
    Cp::Float64 =0 #In J/Kmol
    "Molar volume (J/bar·mol)"
    Vmol::Float64 = 0 #InJ/barmol
    "Cp/Cv (heat capacity ratio)"
    Cp_Cv::Float64 = 0 #Cp/Cv = heat capacity ratio
    "Thermal expansion coefficient(1/K)"
    α::Float64 = 0#Thermal expansion coeficiient in 1/K
    "Compressibility (1/bar)"
    β::Float64 = 0#Compressibility in 1/bar
end


"""
$(TYPEDSIGNATURES)

Returns 'G' of 'phase', useful for broadcasting
"""
function gibbs(phase)
    return phase.G
end


"""
$(TYPEDSIGNATURES)

Returns the 'mol' of 'phase', useful for broadcasting
"""
function mol(phase)
    return phase.mol
end

"""
$(SIGNATURES)
This type is defined to contain all the most relevant properties of a petrological system, any number of variables can be initialized 
and defaults will be 0 or empty arrays/strings. Units are selected based on convenience for petrological modelling.
$(TYPEDFIELDS)
"""
@kwdef struct PetroSystem #A lot of these I can probably remove
    "Composition of the system as defined by a 'Component' array"    
    composition::Array{Component} = Array{Component}([])
    "Phases within the system"
    phases::Array{Phase} = Array{Phase}([])
    "Concentration of trace elements in the system"
    traceelements::Array{TraceElement} = Array{TraceElement}([])
    "Total moles of all components in the system"#Is this actually useful?
    mol::Float64 = 0
    "Total volume of all phases in the system(J/bar [m^3/10^5])"
    vol::Float64 = 0#In J/bar
    "Total mass of the system(g)"
    mass::Float64 = 0#In g
    "Density of the system (kg/m^3)"
    ρ::Float64 = 0#Density in kg/m^3
    "Molar mass of the system (kg/m^3)"
    molarmass::Float64 = 0
    #Thermo properties
    "Gibbs free energy (J/mol)"
    G::Float64 = 0#J/mol
    "Enthalpy (J/mol)"
    H::Float64 = 0#J/mol
    "Entropy (J/K·mol)"
    S::Float64 = 0#In J/Kmol
    "Heat capacity (J/K·mol)"
    Cp::Float64 = 0#In J/Kmol
    "Molar volume (J/bar·mol)"
    Vmol::Float64 = 0#InJ/barmol
    "Cp/Cv (heat capacity ratio)"
    Cp_Cv::Float64 = 0#Cp/Cv = heat capacity ratio
    "Thermal expansion coefficient(1/K)"
    α::Float64 = 0 #Thermal expansion coeficiient in 1/K
    "Compressibility (1/bar)"
    β::Float64 = 0#Compressibility in 1/bar

end

"""
$(TYPEDSIGNATURES)

Calculates the major cations of a phase given the expected number of cations, oxygens and hydroxides
"""
function majorcation(phase, cat, ox, hydrox)
    totalcharge = 2*ox-hydrox
    catcomponents = Array{Component}([])

    for component in phase.composition
        if component.ox > 0
            cat_name = match(r"[^O|\d]*",component.name).match
            if cat_name == "Fe"
                cat_name = "Fe2+"
            end
            cat_mol = concentration(component) * component.cat
            cat_molarmass = (component.molarmass - component.ox*O_MASS)/component.cat
            cat_component = Component(cat_name,cat_molarmass,cat_mol,component.cat,0,component.catcharge)
            push!(catcomponents,cat_component)
        else
            push!(catcomponents,component)
        end
    end

    cat_mol_total = 0
    for component in catcomponents
        if component.name != "H" || component.name != "F" || component.name != "Cl"
            cat_mol_total += concentration(component)
        end
    end

    # cat_mol_total = sum_mols(catcomponents)
    cat_charge_total = 0
    for i in 1:lastindex(catcomponents)
        # mol_norm = concentration(component)*cat/cat_mol_total
        
        catcomponents[i] = catcomponents[i] *(cat/cat_mol_total)
        if catcomponents[i].name != "H" || catcomponents[i].name != "F" || catcomponents[i].name != "Cl"
            cat_charge_total += concentration(catcomponents[i])*catcomponents[i].catcharge
        end
    end

    charge_def = totalcharge-cat_charge_total
    
    fe_index = findchemical(catcomponents, "Fe2+")
    if fe_index > 0
        if charge_def > concentration(catcomponents[fe_index])
            fe3_mol = concentration(catcomponents[fe_index])
            push!(catcomponents,Component("Fe3+",catcomponents[fe_index].molarmass,fe3_mol,1,0,3))
            popat!(catcomponents,fe_index)
            
        elseif charge_def >= 0 
            fe3_mol = charge_def
            push!(catcomponents,Component("Fe3+",catcomponents[fe_index].molarmass,fe3_mol,1,0,3))
            catcomponents[fe_index] -= fe3_mol
        end
    end

    if hydrox > 0
        h_index = findchemical(catcomponents, "H")
        if h_index == 0
            push!(catcomponents, Component("H",H_MASS,0,1,0,1))
            h_index = lastindex(catcomponents)
        end
        f_index = findchemical(catcomponents, "F")
        cl_index = findchemical(catcomponents, "Cl")
        h_total = hydrox

        if f_index > 0
            h_total -= concentration(catcomponents[f_index])
        end
        if cl_index > 0
            h_total -= concentration(catcomponents[cl_index])
        end

        catcomponents[h_index] = Component(catcomponents[h_index],mol = h_total)
   
        
            
    end

    return catcomponents
end

end
