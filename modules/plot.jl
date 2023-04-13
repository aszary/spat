module Plot
    using CairoMakie # 2D vector plots
    using FFTW
    using GLMakie # 3D plots (pngs) # https://docs.makie.org/stable/documentation/figure_size/
    GLMakie.activate!()
    #CairoMakie.activate!()

    include("tools.jl")


    function single0(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, norm=2.0, name_mod="PSR_NAME")

        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        bin_numbers = bin_st:1:bin_end

        # Figure size
        size_inches = (8/2.54, 11/2.54) # 8cm x 11cm
        size_pt = 72 .* size_inches

        fig = Figure(resolution = size_pt, fontsize = 8)
        ax = Axis(fig[1, 1], xlabel="bin number", ylabel="Pulse number", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(ax, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = false, minorticks = false)
        hideydecorations!(ax, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = false, minorticks = false)

        for i in start:1:start+number-1
            da = data[i,:] .* norm .+ i
            da = da[bin_st:bin_end]
            lines!(ax, bin_numbers, da, color=:grey, linewidth=0.3)
            #band!(ax, bin_numbers,  ones(length(da)) * i,  da, color=:white) # work on this one day
        end

        filename = "$outdir/$(name_mod)_single0.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)

    end


    function single(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)
        num, bins = size(data)
        if number === nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        intensity, pulses = Tools.pulses_intensity(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start - 1  # julia

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # Figure size
        size_inches = (8/2.54, 11/2.54) # 8cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution = size_pt, fontsize = 8)

        left = Axis(fig[1:3, 1], ylabel=L"Pulse number $$", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(left,label = true, ticklabels = true, ticks = true, grid = true,    minorgrid = true, minorticks = true)
        hideydecorations!(left, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = true, minorticks = false)
        lines!(left, intensity, pulses, color=:grey, linewidth=0.5)
        ylims!(left, [pulses[1], pulses[end]])

        right = Axis(fig[1:3, 2:3], xlabel="bin number", ylabel="Pulse number", xminorticksvisible=true, yminorticksvisible=true)
        hidedecorations!(right, grid=true, label=true, ticks=true, ticklabels=true)
        #image!(right, transpose(da))
        heatmap!(right, transpose(da))

        bottom = Axis(fig[4, 2:3], xlabel=L"longitude ($^\circ$)", ylabel=L"$$ intensity", xminorticksvisible=true, yminorticksvisible=true) #, alignmode=Mixed()) #, alignmode=Inside()) #Mixed()) #, alignmode = Outside()) # TODO check alignmode
        hidexdecorations!(bottom, label = false, ticklabels = false, ticks = false, grid = true,    minorgrid = true, minorticks = false)
        hideydecorations!(bottom, grid=true, label=false, ticks=false, ticklabels=false,  minorticks=false)
        lines!(bottom, longitude, average, color=:grey, linewidth=0.5)
        xlims!(bottom, [longitude[1], longitude[end]])
        bottom.alignmode = Mixed(left=-34)

        #println(dump(bottom, maxdepth=1))
        #bottom.height=30

        #linkyaxes!(left, right) # nope
        #linkxaxes!(right, bottom) # nope

        rowgap!(fig.layout, 0)
        colgap!(fig.layout, 0)

        filename = "$outdir/$(name_mod)_single.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)

    end

    function test(data)

        #da = transpose(data)
        da = data
        #da = vcat(da[120, :], da[121, :], da[122, :])
        da = da[1:256, :]
        res = []
        for i in 25:55
            res = vcat(res, da[i, :])
        #println(size(res))
        end
        da = convert(Array{Float64}, res)

        println(size(da))
        println(typeof(da))
        #return

        sz = size(da, 1)
        half = floor(Int, sz/2) 
        
        ff = abs.(fft(da))
        freq = fftfreq(sz)

        fig = Figure()
        
        ax = Axis(fig[1, 1])
        lines!(ax, da, color=:black, linewidth=1)

        ax2 = Axis(fig[2, 1])
        #lines!(ax2, ff, color=:red, linewidth=1)
        #lines!(ax2, freq, ff, color=:red, linewidth=1)
        lines!(ax2, freq[1:half], ff[1:half], color=:red, linewidth=1)
        #lines!(fig, 2*da, color=:grey, linewidth=0.5)

        #save("output/test.pdf", fig)
        display(fig)
    end


    #=
    function single(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start - 1  # julia

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))
        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071), frameon=true)  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.09, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da))
        #axvline(x=563, lw=2)
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_single.pdf")
        savefig("$outdir/$(name_mod)_single.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end
    =#



end # module
