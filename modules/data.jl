module Data


    """
    Loads PSRCHIVE ASCII file
    """
    function load_ascii(infile)
        f = open(infile)
        lines = readlines(f)
        res = split(lines[1])
        pulses = parse(Int, res[6])
        bins = parse(Int, res[12])
        data = Array{Float64}(undef, pulses, bins)
        for i in 2: length(lines)
            try
                res = split(lines[i])
                bin = parse(Int, res[3]) + 1
                pulse = parse(Int, res[1]) + 1
                inte = parse(Float64, res[4])
                data[pulse, bin] = inte
            catch
            end
        end
        close(f)
        return data
    end


    """
    Converts PSRFIT file to ASCII (using PSRCHIVE tools)
    """
    function psrfit2ascii(infile, outfile)
        run(pipeline(`pdv -t -F -p $infile`, stdout="$outfile", stderr="errs.txt"))
        # change -t to -A to get frequancy information
    end


end # module
