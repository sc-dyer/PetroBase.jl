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
    Component(name, mMass, mol, μ)
    

Basic unit in describing composition in all phases. Critical to most calculations. 
A component has a name, molar mass, an amount (in mol), and a chemical potential
typically the only values that will ever changed are mol and chemical potential.

"""
struct Component <: Chemical
    name::String
    mMass::Float64 #short for Molar Mass
    mol::Float64
    μ::Float64 #Chemical potential in J/mol
end

struct TraceElem <: Chemical
    name::String
    mMass::Float64
    concentration::Float64
end



#Simple method for being able to copy all values and change the moles
"""
    Component(clone, mol)
    Component(clone;mol,μ)

Additional constructor for Component that allows copying of all parameters while changing the moles or chemical potential
"""
Component(clone::Component, mol::Real) = Component(clone.name, clone.mMass, mol, clone.μ) 
Component(clone::Component;mol::Real=clone.mol,μ::Real=clone.μ) = Component(clone.name, clone.mMass,mol,μ)

"""
    name(val)
    mol(val)
Simple functions for broadcasting on Component arrays
"""
function name(val::Chemical) 
    return val.name
end

function mol(val::Component)
    return val.mol
end
#An operator defined to compare if two components are the same except for the amount of moles
"""
    ≃(comp1, comp2)

This is a simple boolean operator that returns true if the two components are identical except for the moles.
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
    +(comp1, comp2)
    
When two components are added together, they will add their moles to each other if Component1 ≃ Component2.
Otherwise, it will throw an exception.
"""
function Base.:+(comp1::Component, comp2::Component)
    if comp1 ≃ comp2
        return Component(comp1,comp1.mol + comp2.mol)
    else
        throw(ArgumentError("Component objects must match in all parameters except moles"))
    end
end

"""
    +(comp1, num)
    +(num, comp1)

Adds a real number to the moles in the Component
"""
function Base.:+(comp1::Component, num::Real)
    return Component(comp1, comp1.mol+num)
end

Base.:+(num::Real, comp1::Component) = Base.:+(comp1,num)
    
"""
    -(comp1, comp2)
    
When one component is subtracted from another, it will subtract the moles of Component2 from Component1 if Component1 ≃ Component2.
Otherwise, it will throw an exception.
"""
function Base.:-(comp1::Component, comp2::Component)
    if comp1 ≃ comp2
        return Component(comp1,comp1.mol - comp2.mol)
    else
        throw(ArgumentError("Component objects must match in all parameters except mol"))
    end
end

"""
    -(comp1, num)
    
Subtracts a real number from the moles in the Component
"""
function Base.:-(comp1::Component, num::Real)
    return Component(comp1, comp1.mol-num)
end

"""
    *(comp1, num)
    *(num, comp1)
Multiplies the moles of the Component by a Real number
"""
function Base.:*(comp1::Component, num::Real)
    return Component(comp1, comp1.mol*num)
end

Base.:*(num::Real, comp1::Component) = Base.:*(comp1,num)
   

"""
    /(comp1, num)
    
Divides the moles of the Component by a Real number
"""
function Base.:/(comp1::Component, num::Real)
    return Component(comp1, comp1.mol/num)
end
"""
    sumMass(comps)
Calculates the molar mass of an array of Component objects
"""
function sumMass(comps::Array{Component})
    mMass = 0
    for comp in comps
        mMass += comp.mMass*comp.mol
    end
    return mMass
end

"""
    sumMol(comps)
Adds together the moles of all Component objects in comps 
"""
function sumMol(comps::Array{Component})
    mol = 0
    for comp in comps
        mol += comp.mol
    end
    return mol
end
"""
    checkUnique(comps)
Checks if each cell in an array of Component objects is a unique component and not elsewhere in the array
"""
function checkUnique(comps::Array{Component})
    for i in 1:(length(comps)-1), j in (i+1):length(comps) #Does not check against element that already checked the whole array
        if comps[i] ≃ comps[j]
            return false
        end
    end
    return true
end


"""
    findComp(comps, fComp)
Finds the index of fComp in the comps array, where fComp can be either a Component object or the name of a component.
"""
function findComp(comps::Array{Component},fComp::Component)
    for i in 1:length(comps)
        if comps[i] ≃ fComp
            return i
        end
    end
    return 0
end

function findComp(comps::Array{Component}, fComp::String)
    for i in 1:length(comps)
        if comps[i].name == fComp
            return i
        end
    end
    return 0
end

"""
    Phase(name, compo, mol, vol, mass, ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
Phase type contains most relevant properties of a phase, any number of variables can be initialized and defaults will be 0 or empty.
The variable descriptions are as follows:
name::String = name of the phase
comp::Array{Component} = Composition of the phase as defined by array of components
mol::Float64 = Number of moles of the phase present
vol::Float64 = volume of phase present in system (J/bar [m^3/10^5])
mass::Float64 = mass of phase preset in systen (g)
ρ::Float64 = Density (kg/m^3)
mMass::Float64 = molar mass (g/mol)
G::Float64 = Gibbs free energy (J/mol)
H::Float64 = Enthalpy (J/mol)
S::Float64 = Entropy (J/K·mol)
Cp::Float64 = Heat capacity (J/K·mol)
Vmol::Float64 = molar volume (J/bar·mol)
Cp_Cv::Float64 = Cp/Cv or heat capacity ratio
α::Float64 = Thermal expansion coeficiient (1/K)
β::Float64 = Compressibility (1/bar)
"""
struct Phase
    name::String
    compo::Array{Component}
    mol::Float64
    vol::Float64 #In J/bar
    mass::Float64 #In g
    ρ::Float64 #Density in kg/m^3
    mMass::Float64
    #Thermo properties
    G::Float64 #J/mol
    H::Float64 #J/mol
    S::Float64 #In J/Kmol
    Cp::Float64 #In J/Kmol
    Vmol::Float64 #InJ/barmol
    Cp_Cv::Float64 #Cp/Cv = heat capacity ratio
    α::Float64 #Thermal expansion coeficiient in 1/K
    β::Float64 #Compressibility in 1/bar
    Phase(;name = "none", compo = Array{Component}([]),mol = 0, vol = 0, 
            mass = 0, ρ = 0, mMass = 0, G = 0, H = 0, S = 0, Cp = 0, Vmol = 0, 
            Cp_Cv = 0, α = 0, β = 0) = new(name,compo,mol,vol,mass,ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
end

function gibbs(phase::Phase)
    return phase.G
end

function mol(phase::Phase)
    return phase.mol
end
"""
    PetroSystem(name, compo, mol, vol, mass, ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
PetroSystem type contains most relevant properties of a metamorphic system, any number of variables can be initialized and defaults will be 0 or empty.
The variable descriptions are as follows:
comp::Array{Component} = Composition of the system as defined by array of components
phases::Array{Phase} = Phases in the system
mol::Float64 = Number of moles of the phase present
vol::Float64 = volume of phase present in system (J/bar [m^3/10^5])
mass::Float64 = mass of phase preset in systen (g)
ρ::Float64 = Density (kg/m^3)
mMass::Float64 = molar mass (g/mol)
G::Float64 = Gibbs free energy (J/mol)
H::Float64 = Enthalpy (J/mol)
S::Float64 = Entropy (J/K·mol)
Cp::Float64 = Heat capacity (J/K·mol)
Vmol::Float64 = molar volume (J/bar·mol)
Cp_Cv::Float64 = Cp/Cv or heat capacity ratio
α::Float64 = Thermal expansion coeficiient (1/K)
β::Float64 = Compressibility (1/bar)
"""
struct PetroSystem
    compo::Array{Component}
    phases::Array{Phase}
    mol::Float64
    vol::Float64 #In J/bar
    mass::Float64 #In g
    ρ::Float64 #Density in kg/m^3
    mMass::Float64
    #Thermo properties
    G::Float64 #J/mol
    H::Float64 #J/mol
    S::Float64 #In J/Kmol
    Cp::Float64 #In J/Kmol
    Vmol::Float64 #InJ/barmol
    Cp_Cv::Float64 #Cp/Cv = heat capacity ratio
    α::Float64 #Thermal expansion coeficiient in 1/K
    β::Float64 #Compressibility in 1/bar
    PetroSystem(;compo = Array{Component}([]),phases = Array{Phase}([]), mol = 0, vol = 0, 
            mass = 0, ρ = 0, mMass = 0, G = 0, H = 0, S = 0, Cp = 0, Vmol = 0, 
            Cp_Cv = 0, α = 0, β = 0) = new(compo,phases,mol,vol,mass,ρ,mMass,G,H,S,Cp,Vmol,Cp_Cv,α,β)
end


end
