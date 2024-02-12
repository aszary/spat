module SPAT

    include("modules/data.jl")
    include("modules/plot.jl")



    function J0820()

        data1 = Data.load_ascii("input/1.txt")
        data2 = Data.load_ascii("input/2.txt")
        data3 = Data.load_ascii("input/1.debase.p3fold")

        #Plot.test_fft(data1)


        #Plot.single0(data1, "output"; number=150, bin_st=400, bin_end=600, name_mod="J0820")

        #Plot.single0(data1, "output"; bin_st=400, bin_end=600, name_mod="1_0")
        #Plot.single(data1, "output"; start=100, number=50, bin_st=470, bin_end=570, name_mod="1")
        #Plot.single0(data2, "output"; bin_st=400, bin_end=600, name_mod="2_0")
        #Plot.single(data2, "output"; bin_st=400, bin_end=600, name_mod="2")


        #Plot.single(data1, "output"; start=20, bin_st=470, bin_end=550, name_mod="1", number=50)

        # P3 Folds 
        #Plot.lrfs(data1, "output"; bin_st=470, bin_end=550, name_mod="1")
        Plot.p3folded(data1, "output", 4.8; ybins=18, times=10,bin_st=470, bin_end=550, name_mod="1") # simple folding
        Plot.single(data3, "output"; times=10, bin_st=470, bin_end=550, name_mod="1_folded", number=18) # PSRSALSA check /home/szary/work3/MeerKAT/data/J0820/single/salsa1.sh


    end


    function main()

        J0820()
        println("Bye")

    end

end # module

SPAT.main()
