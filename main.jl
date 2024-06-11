module SPAT

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")
    include("modules/psrsalsa.jl")


    function J0820()

        #Data.psrfit2ascii("input/1.lowf", "input/1_low.txt")
        #Data.psrfit2ascii("input/1.highf", "input/1_high.txt")

        data1 = Data.load_ascii("input/1_low.txt")
        data2 = Data.load_ascii("input/1_high.txt")
        #data3 = Data.load_ascii("input/1.txt")



        #Plot.single(data1, "output"; bin_st=430, bin_end=570, name_mod="low")
        #Plot.single(data2, "output"; bin_st=430, bin_end=570, name_mod="high")
        
        #Plot.average(data2, "output"; bin_st=430, bin_end=570, name_mod="high")
        Plot.averageX([data1, data2], "output"; bin_st=475, bin_end=550, name_mod="all")
        #=
        Plot.single(data3, "output"; bin_st=430, bin_end=570, name_mod="1")
        =#
    end

    function J0820_tests()

        data1 = Data.load_ascii("input/1.txt")
        data2 = Data.load_ascii("input/2.txt")
        data3 = Data.load_ascii("input/1.debase.p3fold")

        # convert template to ASCII
        #Data.psrfit2ascii("input/J0820-1350.std", "input/J0820-1350.std.ascii") # TODO why it does not work?

        #Plot.test_fft(data1)


        # RESAMPLE test
        #data2a = Tools.resample(data2, [470, 550], verbose=1)

        #Plot.single0(data1, "output"; number=150, bin_st=470, bin_end=550, name_mod="J0820")

        #Plot.single0(data1, "output"; bin_st=470, bin_end=550, name_mod="1_0")
        #Plot.single(data1, "output"; start=100, number=50, bin_st=470, bin_end=550, name_mod="1")
        #Plot.single0(data2, "output"; bin_st=470, bin_end=550, name_mod="2_0")
        
        #Plot.single(data2, "output"; bin_st=470, bin_end=550, name_mod="2")
        #Plot.single(data2a, "output"; bin_st=470, bin_end=550, name_mod="2a")

        #Plot.single(data2b, "output"; bin_st=470, bin_end=550, name_mod="2b")

        #Plot.single(data1, "output"; start=20, bin_st=470, bin_end=550, name_mod="1", number=150)

        # P3 Folds 
        #Plot.lrfs(data2, "output"; bin_st=470, bin_end=550, name_mod="2", fix_fftphase=true)
        #Plot.lrfs(data2a, "output"; bin_st=470, bin_end=550, name_mod="2a")
        #Plot.p3folded(data1, "output", 4.8; ybins=18, times=10,bin_st=470, bin_end=550, name_mod="1") # simple folding
        #Plot.single(data3, "output"; times=10, bin_st=470, bin_end=550, name_mod="1_folded", number=18) # PSRSALSA check /home/szary/work3/MeerKAT/data/J0820/single/salsa1.sh # not the pulse number on the left (ybin number)!

        # moving to scripts...
        #PSRSALSA.debase("input/1.spCF", [473, 556])
        #PSRSALSA.debase("input/1.spCF", nothing)

    end


    function main()

        J0820()
        #J0820_tests()
        println("Bye")

    end

end # module

SPAT.main()
