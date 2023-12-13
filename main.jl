module SPAT

    include("modules/data.jl")
    include("modules/plot.jl")



    function J0820()

        data = Data.load_ascii("input/1.txt")
        data2 = Data.load_ascii("input/2.txt")

        #Plot.test_fft(data)


        #Plot.single0(data, "output"; number=150, bin_st=400, bin_end=600, name_mod="J0820")

        #Plot.single0(data, "output"; bin_st=400, bin_end=600, name_mod="1_0")
        Plot.single(data, "output"; start=100, number=50, bin_st=470, bin_end=570, name_mod="1")
        #Plot.single0(data2, "output"; bin_st=400, bin_end=600, name_mod="2_0")
        Plot.single(data2, "output"; bin_st=400, bin_end=600, name_mod="2")

    end


    function main()

        J0820()
        println("Bye")

    end

end # module

SPAT.main()
