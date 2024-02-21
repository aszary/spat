module PSRSALSA
    using Expect

    """
    Dabase signal using PSRSALSA 
    
    This approach may be usefull for multiple sources, but for one source just write .sh scripts for PSRSALSA
    
        https://docs.juliahub.com/General/Expect/stable/#
    """
    function debase(infile, on_bins; verbose=1)
        if on_bins === nothing # enter interactive mode
            cmd = `pmod -debase $infile`
            proc = ExpectProc(cmd, 100)
            lines = expect!(proc, "/NULL):")
            if verbose > 0
                println(lines,  "/NULL):")
            end
            println(proc, "/xw")
            lines2 = expect!(proc, "PGPLOT window.")
            println(lines2, " PGPLOT window.")
        else # just run it when on_bins are specified
            cmd = `pmod -debase -onpulse "$(on_bins[1]) $(on_bins[2])" $infile`
            proc = ExpectProc(cmd, 100)
            println(proc, "/NULL")
            expect!(proc, "done")
            println(proc.before, " done")
        end
    end

end # module