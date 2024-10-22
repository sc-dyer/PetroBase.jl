using PetroBase
using Test

@testset "PetroBase.jl" begin
    # Write your tests here.

    #Make components with different constructors to make sure they compile properly
    comp1 = Component("SiO2",60.08,1.0,1,2,4,0.0)
    comp1a = Component("SiO2",60.08,1)
    comp2 = Component(comp1,mol=3.0)
    comp3 = Component(comp1,μ=300.0)
    comp4 = Component("MgO",40.304,2.0,1,1,2,0.0)
    #Same with trace elements
    te1 = TraceElement("Y",88.9059,400.0)
    te2 = TraceElement(te1,900.0)
    te3 = TraceElement("Lu", 174.967, 50.0)

    #Testing the operators and simple method calls
    @test name(comp2) == "SiO2"
    @test name(te2)  == "Y"
    @test comp1 ≃ comp2
    @test comp2 ≃ comp3
    @test comp1 ≈ comp1a
    @test !(comp1 ≈ comp2)
    @test !(comp1 ≃ comp4)
    @test te1 ≃ te2
    @test concentration(comp3) ≈ 1
    @test concentration(te1) ≈ 400
    @test concentration(comp1 + comp2) ≈ 4
    @test concentration(1 + comp2) ≈ 4
    @test concentration(te1+te2) ≈ 1300
    @test concentration(400 + te2) ≈ 1300
    @test concentration(comp2 - comp1) ≈ 2
    @test concentration(comp2 - 1) ≈ 2
    @test concentration(te2 - te1) ≈ 500
    @test concentration(te2 - 400) ≈ 500
    @test concentration(10*comp1) ≈ 10
    @test concentration(comp1/10) ≈ 0.1
    @test concentration(10*te1) ≈ 4000
    @test concentration(te1/10) ≈ 40
    @test_throws ArgumentError comp1 + comp4
    @test_throws ArgumentError te1 + te3
    @test_throws ArgumentError comp4 - comp1
    @test_throws ArgumentError te3-te1

    
    compoList1 = [comp1,comp4]
    compoList2 = [comp1,comp2,comp3]
   
    teList1 = [te1,te3]
    teList2 = [te1,te2,te3]
    #Now testing some more complicated functions
    @test sum_mass(compoList1) ≈ 140.688
    @test sum_mols(compoList1) ≈ 3
    @test isunique(compoList1)
    @test !isunique(compoList2)
    @test isunique(teList1)
    @test !isunique(teList2)
    @test findchemical(compoList1,comp4) == 2
    @test findchemical(compoList1,"MgO") == 2
    @test findchemical(compoList1,Component("K2O",94.2,3.0)) == 0
    @test findchemical(compoList1, "K2O") == 0
    @test findchemical(teList1,te3) == 2
    @test findchemical(teList1,"Lu") ==2

    compstart = Component("SiO2",60.08,1)
    compend = Component("SiO2",60.08,10)
    comprange = range(compstart,compend,length=10)
    @test concentration(comprange[4]) ≈ 4.0

    phase1 = Phase(name = "Forsterite", composition = compoList1, traceelements = teList1, mol = 3.0, G = -2403.2, vol = 3.0)
    @test gibbs(phase1) ≈ -2403.2
    @test mol(phase1) ≈ 3

    system = PetroSystem(composition = compoList1, phases = [phase1], traceelements = teList1, mol = 3.0, G = -2403.2)

    @test mol(getphase(system,"Forsterite")[1]) ≈ 3.0
    @test mol(getphase(system, "Ilmenite")[1]) ≈ 0.0

    c1 = Component("FeO",71.850,1,1,2,mass=45.333)
    c2 = Component("TiO2",79.90,1,2,4,mass=52.7124)
    c3 = Component("MnO",70.937,1,1,2,mass=2.64039)

    ilm = Phase(name = "Ilmenite",composition = [c1,c2,c3], mol = 1.0,vol = 7.0)
    ilmcat = majorcation(ilm,2,3,0)
    @test length(ilmcat) == 4
    @test round(concentration(ilmcat[findchemical(ilmcat,"Ti")]),sigdigits=3) ≈ 0.994
    @test round(concentration(ilmcat[findchemical(ilmcat,"Mn")]),sigdigits=2) ≈ 0.056
    @test round(concentration(ilmcat[findchemical(ilmcat,"Fe2+")]),sigdigits=3) ≈ 0.938
    @test round(concentration(ilmcat[findchemical(ilmcat,"Fe3+")]),sigdigits=2) ≈ 0.013

    system = PetroSystem(composition = compoList1, phases = [phase1,ilm], traceelements = teList1, mol = 3.0, G = -2403.2)

    @test get_volprop(system,"forsterite") ≈ 0.3
    @test get_volprop(system,"ilmenite") ≈ 0.7
end
