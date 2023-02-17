module SPAT

    include("modules/data.jl")
    include("modules/plot.jl")


    function J0820()

        data = Data.load_ascii("input/low.txt")
        data2 = Data.load_ascii("input/high.txt")

        #Plot.single0(data, "output"; number=150, bin_st=400, bin_end=600, name_mod="J0820")

        Plot.single(data, "output"; bin_st=400, bin_end=600, name_mod="low")
    end


    function main()

        J0820()
        println("Bye")

    end

end # module

SPAT.main()
