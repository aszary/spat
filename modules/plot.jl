module Plot
    using CairoMakie # 2D vector plots
    # using GLMakie # 3D plots (pngs) # https://docs.makie.org/stable/documentation/figure_size/

    function single0(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, norm=2.0, name_mod="PSR_NAME")

        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        bin_numbers = bin_st:1:bin_end

        size_inches = (8/2.54, 11/2.54) # 8cm x 11cm
        size_pt = 72 .* size_inches

        fig = Figure(resolution = size_pt, fontsize = 8)
        ax = Axis(fig[1, 1], xlabel="bin number", ylabel="Pulse number", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(ax, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = false, minorticks = false)
        hideydecorations!(ax, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = false, minorticks = false)

        for i in start:1:start +number-1
            da = data[i,:] .* norm .+ i
            da = da[bin_st:bin_end]
            lines!(ax, bin_numbers, da, color=:grey, linewidth=0.3)
            #band!(ax, bin_numbers,  ones(length(da)) * i,  da, color=:white) # work on this one day
        end

        filename = "$outdir/$(name_mod)_single0.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)

    end


end # module
