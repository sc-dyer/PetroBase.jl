"""
Main module for PetroBase, provides basic types and functions to solve petrological problems.

# Exports
$(EXPORTS)
"""
module PetroBase

export 
    Chemical, 
    Component, 
    TraceElem, 
    Phase, 
    PetroSystem, 
    name, 
    mol, 
    sumMass, 
    sumMol, 
    checkUnique,
    findComp,
    gibbs

using
    DocStringExtensions
# Write your package code here.

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
    mMass::Float64 #short for Molar Mass
    "Moles of the component"
    mol::Float64
    "Chemical potential (J/mol)"
    μ::Float64 #Chemical potential in J/mol
end
#Simple method for being able to copy all values and change the moles
"""
$(SIGNATURES)
Clones the parameters of a 'Component' but with change of 'mol' and/or 'μ'
"""
Component(clone::Component, mol::Real = clone.mol;μ::Real=clone.μ) = Component(clone.name, clone.mMass, mol, μ)


"""
$(SIGNATURES)
Describes the trace element concentration of a phase or system. Typically the only value that will change is concentration

$(TYPEDFIELDS)
"""
struct TraceElem <: Chemical
    "Name of the element"
    name::String
    "Molar mass (g/mol)"
    mMass::Float64
    "Concentration in the system or phase (µg/g)"
    concentration::Float64
    
end

"""
$(SIGNATURES)
Clones the parameters of a 'TraceElem' but with change of 'concentration'
"""
TraceElem(clone::TraceElem,conc::Real) = TraceElem(clone.name,clone.mMass,conc)

"""
$(TYPEDSIGNATURES)

Returns 'name' of a 'chem', useful for broadcasting
"""
function name(chem::Chemical) 
    return chem.name
end

"""
$(TYPEDSIGNATURES)

Returns the 'mol' of 'comp', useful for broadcasting
"""
function conc(comp::Component)
    return comp.mol
end

"""
$(TYPEDSIGNATURES)

Returns the 'concentration' of a 'TraceElem', useful for broadcasting
"""
function conc(te::TraceElem)
    return te.concentration
end
#An operator defined to compare if two components are the same except for the amount of moles
"""
$(TYPEDSIGNATURES)

This is a simple boolean operator that returns true if two 'Chemical' variables have the same name and molar mass
"""
function ≃(comp1::Chemical, comp2::Chemical)
    if comp1.name == comp2.name && comp1.mMass ≈ comp2.mMass
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
        return Component(comp1,comp1.mol + comp2.mol)
    else
        throw(ArgumentError("comp1 and comp2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Adds 'num' to the 'mol' of 'comp1'
"""
function Base.:+(comp1::Component, num::Real)
    return Component(comp1, comp1.mol+num)
end

Base.:+(num::Real, comp1::Component) = Base.:+(comp1,num)

"""
$(TYPEDSIGNATURES)
    
Adds together the 'concentration' parameter of two 'TraceElem' variables that have the same name and molar mass
"""
function Base.:+(te1::TraceElem,te2::TraceElem)
    if te1 ≃ te2
        return TraceElem(te1,te1.concentration + ste2.concentration)
    else
        throw(ArgumentError("te1 and te2 must have the same name and molar mass"))
end

"""
$(TYPEDSIGNATURES)
    
Adds 'num' to the 'concentration' paremeter of 'te1'
"""
function Base.:+(te1::TraceElem,num::Real)
    return TraceElem(te1,te1.mol+num)
end

Base.:+(num::Real,te1::TraceElem) = Base.:+(te1,num)

"""
$(TYPEDSIGNATURES)

Subtracts the 'mol' parameter of 'comp2' from 'comp1' as long as they have the same name and molar mass.
This will return a 'Component' with the chemical potential of the first argument.
"""
function Base.:-(comp1::Component, comp2::Component)
    if comp1 ≃ comp2
        return Component(comp1,comp1.mol - comp2.mol)
    else
        throw(ArgumentError("comp1 and comp2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Subtracts 'num' to the 'mol' parameter of 'comp1'
"""
function Base.:-(comp1::Component, num::Real)
    return Component(comp1, comp1.mol-num)
end

"""
$(TYPEDSIGNATURES)

Subtracts the 'concentration' parameter of 'te2' from 'te1' as long as they have the same name and molar mass
"""
function Base.:-(te1::TraceElem,te2::TraceElem)
    if te1 ≃ te2
        return TraceElem(te1,te1.mol - te2.mol)
    else
        throw(ArgumentError("te1 and te2 must have the same name and molar mass"))
    end
end

"""
$(TYPEDSIGNATURES)
    
Subtracts 'num' to the 'concentration' parameter of 'te1'
"""
function Base.:-(te1::TraceElem, num::Real)
    return TraceElem(te1, te1.mol-num)
end

"""
$(TYPEDSIGNATURES)

Multiplies the 'mol' parameter of 'comp1' by 'num'
"""
function Base.:*(comp1::Component, num::Real)
    return Component(comp1, comp1.mol*num)
end

Base.:*(num::Real, comp1::Component) = Base.:*(comp1,num)
   
"""
$(TYPEDSIGNATURES)

Multiplies the 'concentration' parameter of 'te1' by 'num'
"""
function Base.:*(te1::TraceElem, num::Real)
    return Component(te1, te1.mol*num)
end

Base.:*(num::Real, te1::TraceElem) = Base.:*(te1,num)

"""
$(TYPEDSIGNATURES)

Divides the 'mol' parameter of 'comp1' by 'num'
"""
function Base.:/(comp1::Component, num::Real)
    return Component(comp1, comp1.mol/num)
end

"""
$(TYPEDSIGNATURES)

Divides the 'concentration' parameter of 'te1' by 'num'
"""
function Base.:/(te1::TraceElem, num::Real)
    return Component(te1, te1.mol/num)
end

"""
$(TYPEDSIGNATURES)
Calculates the molar mass of an array of 'Component' variables by adding up the product of the molar mass and mol of each 'Component' in the array.
Intent of use is calculating molar mass of a phase that is described using a 'Component' array
"""
function sumMass(comps::Array{Component})
    mMass = 0
    for comp in comps
        mMass += comp.mMass*comp.mol
    end
    return mMass
end


"""
$(TYPEDSIGNATURES)
Checks if each cell in a 'Chemical' array is not repeated elsewhere in the array.
"""
function checkUnique(chem::Array{<:Chemical})
    for i in 1:(lastindex(chem)-1), j in (i+1):lastindex(chem) #Does not check against element that already checked the whole array
        if chem[i] ≃ chem[j]
            return false
        end
    end
    return true
end


"""
$(TYPEDSIGNATURES)
Finds the first index of 'fChem' in the 'chem' array. Returns 0 if 'fChem' isnt present. Best used with arrays of unique 'Chemical' variables.
"""
function findChem(chem::Array{<:Chemical},fChem::Chemical)
    for i in 1:lastindex(chem)
        if chem[i] ≃ fChem
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
function findChem(chem::Array{Chemical}, fChem::String)
    for i in 1:lastindex(chem)
        if chem[i].name == fChem
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
struct Phase
    "Name of the phase"
    name::String
    "Composition of the phase as defined by a 'Component' array"
    compo::Array{Component}
    "Concentration of trace elements in the phase"
    traceElems::Array{TraceElem}#In general I anticipate it is more useful to have TEs and Components seperate
    "Number of moles of the phase present in the system"
    mol::Float64
    "Volume of phase present in the system (J/bar [m^3/10^5])"
    vol::Float64 #In J/bar
    "Total mass of the phase present in the system (g)"
    mass::Float64 #In g
    "Density (kg/m^3)"
    ρ::Float64 #Density in kg/m^3
    "Molar mass (g/mol)"
    mMass::Float64
    #Thermo properties
    "Gibbs free energy (J/mol)"
    G::Float64 #J/mol
    "Enthalpy (J/mol)"
    H::Float64 #J/mol
    "Entropy (J/mol)"
    S::Float64 #In J/Kmol
    "Heat capacity (J/K·mol)"
    Cp::Float64 #In J/Kmol
    "Molar volume (J/bar·mol)"
    Vmol::Float64 #InJ/barmol
    "Cp/Cv (heat capacity ratio)"
    Cp_Cv::Float64 #Cp/Cv = heat capacity ratio
    "Thermal expansion coefficient(1/K)"
    α::Float64 #Thermal expansion coeficiient in 1/K
    "Compressibility (1/bar)"
    β::Float64 #Compressibility in 1/bar
    Phase(;name = "none", compo = Array{Component}([]),traceElems = Array{TraceElem}([]),mol = 0, vol = 0, 
            mass = 0, ρ = 0, mMass = 0, G = 0, H = 0, S = 0, Cp = 0, Vmol = 0, 
            Cp_Cv = 0, α = 0, β = 0) = new(name,compo,traceElems,mol,vol,mass,ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
end


"""
$(TYPEDSIGNATURES)

Returns 'G' of 'phase', useful for broadcasting
"""
function gibbs(phase::Phase)
    return phase.G
end


"""
$(TYPEDSIGNATURES)

Returns the 'mol' of 'phase', useful for broadcasting
"""
function mol(phase::Phase)
    return phase.mol
end

"""
$(SIGNATURES)
This type is defined to contain all the most relevant properties of a petrological system, any number of variables can be initialized 
and defaults will be 0 or empty arrays/strings. Units are selected based on convenience for petrological modelling.
$(TYPEDFIELDS)
"""
struct PetroSystem #A lot of these I can probably remove
    "Composition of the system as defined by a 'Component' array"    
    compo::Array{Component}
    "Phases within the system"
    phases::Array{Phase}
    "Concentration of trace elements in the system"
    traceElems::Array{TraceElem}
    "Total moles of all components in the system"#Is this actually useful?
    mol::Float64
    "Total volume of all phases in the system(J/bar [m^3/10^5])"
    vol::Float64 #In J/bar
    "Total mass of the system(g)"
    mass::Float64 #In g
    "Density of the system (kg/m^3)"
    ρ::Float64 #Density in kg/m^3
    "Molar mass of the system (kg/m^3)"
    mMass::Float64
    #Thermo properties
    "Gibbs free energy (J/mol)"
    G::Float64 #J/mol
    "Enthalpy (J/mol)"
    H::Float64 #J/mol
    "Entropy (J/K·mol)"
    S::Float64 #In J/Kmol
    "Heat capacity (J/K·mol)"
    Cp::Float64 #In J/Kmol
    "Molar volume (J/bar·mol)"
    Vmol::Float64 #InJ/barmol
    "Cp/Cv (heat capacity ratio)"
    Cp_Cv::Float64 #Cp/Cv = heat capacity ratio
    "Thermal expansion coefficient(1/K)"
    α::Float64 #Thermal expansion coeficiient in 1/K
    "Compressibility (1/bar)"
    β::Float64 #Compressibility in 1/bar
    PetroSystem(;compo = Array{Component}([]),phases = Array{Phase}([]), traceElems = Array{TraceElem}([]),mol = 0, vol = 0, 
            mass = 0, ρ = 0, mMass = 0, G = 0, H = 0, S = 0, Cp = 0, Vmol = 0, 
            Cp_Cv = 0, α = 0, β = 0) = new(compo,phases,traceElems,mol,vol,mass,ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
end


end
