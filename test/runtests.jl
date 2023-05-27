using PetroBase
using Test

@testset "PetroBase.jl" begin
    # Write your tests here.

    #Make components with different constructors to make sure they compile properly
    comp1 = Component("SiO2",60.08,1.0,0.0)
    comp2 = Component(comp1,3.0)
    comp3 = Component(comp1,μ=300.0)
    comp4 = Component("MgO",40.304,2.0,0.0)
    #Same with trace elements
    te1 = TraceElem("Y",88.9059,400.0)
    te2 = TraceElem(te1,900.0)
    te3 = TraceElem("Lu", 174.967, 50.0)

    #Testing the operators and simple method calls
    @test name(comp2) == "SiO2"
    @test name(te2)  == "Y"
    @test comp1 ≃ comp2
    @test comp2 ≃ comp3
    @test !(comp1 ≃ comp4)
    @test te1 ≃ te2
    @test conc(comp3) ≈ 1
    @test conc(te1) ≈ 400
    @test conc(comp1 + comp2) ≈ 4
    @test conc(1 + comp2) ≈ 4
    @test conc(te1+te2) ≈ 1300
    @test conc(comp2 - comp1) ≈ 2
    @test conc(comp2 - 1) ≈ 2
    @test conc(te2 - te1) ≈ 500
    @test conc(te2 - 400) ≈ 500
    @test conc(10*comp1) ≈ 10
    @test conc(comp1/10) ≈ 0.1
    @test conc(10*te1) ≈ 4000
    @test conc(te1/10) ≈ 40
    @test_throws ArgumentError comp1 + comp4
    @test_throws ArgumentError te1 + te3
    @test_throws ArgumentError comp4 - comp1
    @test_throws ArgumentError te3-te1

    
    compoList1 = [comp1,comp4]
    compoList2 = [comp1,comp2,comp3]
   
    teList1 = [te1,te3]
    teList2 = [te1,te2,te3]
    #Now testing some more complicated functions
    @test sumMass(compoList1) ≈ 140.688
    @test checkUnique(compoList1)
    @test !checkUnique(compoList2)
    @test checkUnique(teList1)
    @test !checkUnique(teList2)
    @test findChem(compoList1,comp4) == 2
    @test findChem(compoList1,"MgO") == 2
    @test findChem(compoList1,Component("K2O",94.2,3.0,0.0)) == 0
    @test findChem(compoList1, "K2O") == 0
    @test findChem(teList1,te3) == 2
    @test findChem(teList1,"Lu") ==2

    phase1 = Phase(name = "Forsterite", compo = compoList1, traceElems = teList1, mol = 3.0, G = -2403.2)
    @test gibbs(phase1) ≈ -2403.2
    @test mol(phase1) ≈ 3

    system = PetroSystem(compo = compoList1, phases = [phase1], traceElems = teList1, mol = 3.0, G = -2403.2)

end
