module Plot
using CairoMakie # 2D vector plots
using GLMakie # 3D plots (pngs) # https://docs.makie.org/stable/documentation/figure_size/
#GLMakie.activate!()
CairoMakie.activate!()
using FFTW
using StatsBase
using DSP

include("tools.jl")


struct Panels
    left
    right
    top
    bottom
    center
end


function triple_panels()
    
    # Figure size
    size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
    size_pt = 72 .* size_inches
    #println(size_pt)

    fig = Figure(resolution=size_pt, fontsize=8)

    left = Axis(fig[1:6, 1], xminorticksvisible=true, yminorticksvisible=true, xticks=[0.5]) # 2:6, 1
    #top = Axis(fig[1, 2:3], xaxisposition=:top, yaxisposition = :right)
    center = Axis(fig[1:6, 2:3]) # 2:6, 2:3
    bottom = Axis(fig[7, 2:3], yaxisposition = :left, xminorticksvisible=true, yminorticksvisible=true, yticklabelsvisible=false)

    left.xreversed=true

    hidedecorations!.(center)
    hidedecorations!.(left, grid=true, ticks=false, ticklabels=false, label=false, minorticks=false)
    #hidedecorations!.(top, grid=false, ticks=false, ticklabels=false)
    hidedecorations!.(bottom, grid=true, ticks=false, ticklabels=false,label=false, minorticks=false)
    left.alignmode = Mixed(bottom = MakieLayout.Protrusion(0))
    #top.alignmode = Mixed(left = MakieLayout.Protrusion(0))
    bottom.alignmode = Mixed(left = MakieLayout.Protrusion(0))

    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)

    return fig, Panels(left, nothing, nothing, bottom, center)

end


function quad_panels()

    # Figure size
    size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
    size_pt = 72 .* size_inches
    #println(size_pt)

    fig = Figure(resolution=size_pt, fontsize=8)

    left = Axis(fig[2:6, 1], xminorticksvisible=true, yminorticksvisible=true, xticks=[0.5])
    top = Axis(fig[1, 2:3], xaxisposition=:top, yaxisposition = :left, xminorticksvisible=true, yminorticksvisible=true)
    center = Axis(fig[2:6, 2:3])
    bottom = Axis(fig[7, 2:3], yaxisposition = :left, xminorticksvisible=true, yminorticksvisible=true, yticklabelsvisible=false)

    left.xreversed=true

    hidedecorations!.(top, grid=true, ticks=false, ticklabels=false, label=false, minorticks=false)
    hidedecorations!.(center)
    hidedecorations!.(left, grid=true, ticks=false, ticklabels=false, label=false, minorticks=false)
    hidedecorations!.(bottom, grid=true, ticks=false, ticklabels=false,label=false, minorticks=false)
    left.alignmode = Mixed(bottom = MakieLayout.Protrusion(0))
    top.alignmode = Mixed(left = MakieLayout.Protrusion(0))
    bottom.alignmode = Mixed(left = MakieLayout.Protrusion(0))

    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)

    return fig, Panels(left, nothing, top, bottom, center)

end


function single0(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, norm=2.0, name_mod="PSR_NAME")

    num, bins = size(data)
    if number == nothing
        number = num - start  # missing one?
    end
    if bin_st == nothing
        bin_st = 1
    end
    if bin_end == nothing
        bin_end = bins
    end

    bin_numbers = bin_st:1:bin_end

    # Figure size
    size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=8)
    ax = Axis(fig[1, 1], xlabel=L"bin number $$", ylabel=L"Pulse number $$", xminorticksvisible=true, yminorticksvisible=true)
    hidexdecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)
    hideydecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

    for i in start:1:start+number-1
        da = data[i, :] .* norm .+ i
        da = da[bin_st:bin_end]
        lines!(ax, bin_numbers, da, color=:grey, linewidth=0.3)
        #band!(ax, bin_numbers,  ones(length(da)) * i,  da, color=:white) # work on this one day
    end

    filename = "$outdir/$(name_mod)_single0.pdf"
    println(filename)
    save(filename, fig, pt_per_unit=1)

end


function single_old(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)
    num, bins = size(data)
    if number === nothing
        number = num - start  # missing one?
    end
    if bin_st == nothing
        bin_st = 1
    end
    if bin_end == nothing
        bin_end = bins
    end
    da = data[start:start+number-1, bin_st:bin_end]
    average = Tools.average_profile(da)
    intensity, pulses = Tools.pulses_intensity(da)
    intensity .-= minimum(intensity)
    intensity ./= maximum(intensity)

    pulses .+= start - 1  # julia

    # Pulse longitude
    db = (bin_end + 1) - bin_st  # yes +1
    dl = 360.0 * db / bins
    longitude = collect(range(-dl / 2.0, dl / 2.0, length=db))

    # Figure size
    size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
    size_pt = 72 .* size_inches
    #println(size_pt)
    fig = Figure(resolution=size_pt, fontsize=8)

    left = Axis(fig[1:3, 1], ylabel=L"Pulse number $$", xminorticksvisible=true, yminorticksvisible=true)
    hidexdecorations!(left, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
    hideydecorations!(left, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
    lines!(left, intensity, pulses, color=:grey, linewidth=0.5)
    ylims!(left, [pulses[1], pulses[end]])

    right = Axis(fig[1:3, 2:3], xlabel="bin number", ylabel="Pulse number", xminorticksvisible=true, yminorticksvisible=true)
    hidedecorations!(right, grid=true, label=true, ticks=true, ticklabels=true)
    #image!(right, transpose(da))
    heatmap!(right, transpose(da))

    bottom = Axis(fig[4, 2:3], xlabel=L"longitude ($^\circ$)", ylabel=L"$$ intensity", xminorticksvisible=true, yminorticksvisible=true) #, alignmode=Mixed()) #, alignmode=Inside()) #Mixed()) #, alignmode = Outside()) # TODO check alignmode
    hidexdecorations!(bottom, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
    hideydecorations!(bottom, grid=true, label=false, ticks=false, ticklabels=false, minorticks=false)
    lines!(bottom, longitude, average, color=:grey, linewidth=0.5)
    xlims!(bottom, [longitude[1], longitude[end]])
    bottom.alignmode = Mixed(left=-34) # hack 

    #println(dump(bottom, maxdepth=1))
    #bottom.height=30

    #linkyaxes!(left, right) # nope
    #linkxaxes!(right, bottom) # nope

    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 0)

    filename = "$outdir/$(name_mod)_single.pdf"
    println(filename)
    save(filename, fig, pt_per_unit=1)

end


function single(data, outdir; start=1, number=100, times=1, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)

    # PREPARE DATA
    num, bins = size(data)
    if number === nothing
        number = num - start  # missing one?
    end
    if bin_st === nothing
        bin_st = 1
    end
    if bin_end === nothing
        bin_end = bins
    end
    da = data[start:start+number-1, bin_st:bin_end]
    da = repeat(da, times) # repeat data X times
    average = Tools.average_profile(da)
    intensity, pulses = Tools.pulses_intensity(da)
    intensity .-= minimum(intensity)
    intensity ./= maximum(intensity)

    pulses .+= start - 1  # julia

    # Pulse longitude
    db = (bin_end + 1) - bin_st  # yes +1
    dl = 360.0 * db / bins
    longitude = collect(range(-dl / 2.0, dl / 2.0, length=db))

    # CREATE FIGURE
    fig, p = triple_panels()
    p.left.ylabel = L"Pulse number $$"
    p.left.xlabel = L"intensity $$"
    p.bottom.xlabel = L"longitude ($^\circ$)"

    # PLOTTING DATA
    lines!(p.left, intensity, pulses, color=:grey, linewidth=0.5)
    #xlims!(left, [0.01, 1.01])
    ylims!(p.left, [pulses[1]-0.5, pulses[end]+0.5])

    heatmap!(p.center, transpose(da))

    lines!(p.bottom, longitude, average, color=:grey, linewidth=0.5)
    xlims!(p.bottom, [longitude[1], longitude[end]])

    #screen = display(fig)
    #resize!(screen, 500, 800)
    filename = "$outdir/$(name_mod)_single.pdf"
    println(filename)
    save(filename, fig, pt_per_unit=1)

end


function lrfs(data, outdir; start=1, number=nothing, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", skip_firstfreq=true, fix_fftphase=true)

    # PREPARE DATA
    num, bins = size(data)
    println("$num $bins")
    if number === nothing
        number = num - start + 1   # missing one?
    end
    if bin_st === nothing bin_st = 1 end
    if bin_end === nothing bin_end = bins end

    da = data[start:start+number-1, bin_st:bin_end]

    average = Tools.average_profile(da)
    lrfs, intensity, freq, peak = Tools.lrfs(da)
    println("\tpeak freq $(freq[peak]) ")
    println("\tpeak P3 $(1/freq[peak])")
   
    # skip freq = 0 and normalize intensity to 1
    if skip_firstfreq
        inten = intensity[2:end]
        fre = freq[2:end]
    else
        inten = intensity[1:end]
        fre = freq[1:end]
    end
    inten .-= minimum(inten)
    inten ./= maximum(inten)
    
    # Pulse longitude
    db = (bin_end + 1) - bin_st  # yes +1
    dl = 360. * db / bins
    longitude = collect(range(-dl/2., dl/2., length=db))

    # Finding P3
    pars, errs = Tools.fit_gaussian(fre, inten; μ=freq[peak-1])  # skip zero freq
    fr = pars[2]
    frer = pars[3] # errs[2] # greater taken?

    println("\tFrequency (gaussian fit): $fr, P3: $(1/fr)")
    f_min = fr - frer
    f_max = fr + frer
    # P3 error is it OK?
    dP1 = 1 / fr - 1 / f_max
    dP2 = 1 / f_min  - 1 / fr
    dP = maximum([dP1, dP2])
    println("\tFrequency error (gaussian fit): $frer, P3 error: $dP")

    # LRFS phase  
    phase_ = rad2deg.(angle.(view(lrfs, peak, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data

    if fix_fftphase == true
        Tools.fix_fftphase!(phase_)
    end
    
    # bootstrap scheme to evaluate uncertinies
    phases, phase_, ephase_minus, ephase_plus = Tools.phase_errors(data[start:start+number-1,:], [473, 556]; num=100, bin_range=[bin_st, bin_end], fix_fftphase=fix_fftphase)

    # CREATE FIGURE
    fig, p = quad_panels()
    p.left.ylabel = L"frequency $(1/P)$"
    p.left.xlabel = L"intensity $$"
    p.bottom.xlabel = L"longitude ($^\circ$)"
    p.top.xlabel = L"longitude ($^\circ$)"
    p.top.ylabel = L"FFT phase ($^\circ$)"

    # PLOTTING DATA
    # all bootstrap data
    for phase in phases
        scatter!(p.top, longitude, phase, color=:grey, markersize=0.7)
    end
    # errors bars based on bootstarping scheme # very small?
    #scatter!(p.top, longitude, phase_, color=:grey, markersize=1)
    #errorbars!(p.top, longitude, phase_, ephase_minus, ephase_plus, color=:red, whiskerwidth = 0, linewidth=0.2)


    xlims!(p.top, [longitude[1], longitude[end]])
    lines!(p.left, inten, fre, color=:grey, linewidth=0.5)
    lines!(p.left, Tools.gauss(fre, pars), fre, color=(:red, 0.3), linewidth=0.3)
    #xlims!(left, [0.01, 1.01])
    ylims!(p.left, [fre[1], fre[end]])
    if skip_firstfreq
        heatmap!(p.center, transpose(abs.(lrfs[2:end,:])))
    else
        heatmap!(p.center, transpose(abs.(lrfs)))
    end
    lines!(p.bottom, longitude, average, color=:grey, linewidth=0.5)
    xlims!(p.bottom, [longitude[1], longitude[end]])

    filename = "$outdir/$(name_mod)_lrfs.pdf"
    println(filename)
    save(filename, fig, pt_per_unit=1)
end


function p3folded(data, outdir, p3; ybins=18, start=1, number=nothing, times=10, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)

    # PREPARE DATA
    num, bins = size(data)
    if number === nothing
        number = num - start  # missing one?
    end
    if bin_st == nothing
        bin_st = 1
    end
    if bin_end == nothing
        bin_end = bins
    end
    da = data[start:start+number-1, bin_st:bin_end]
    p3fold = Tools.p3fold(da, p3; ybins=ybins)
    da = repeat(p3fold, times) # p3 folds data from now on..

    average = Tools.average_profile(da)
    intensity, pulses = Tools.pulses_intensity(da)
    intensity .-= minimum(intensity)
    intensity ./= maximum(intensity)
    pulses .+= start - 1  # julia

    # Pulse longitude
    db = (bin_end + 1) - bin_st  # yes +1
    dl = 360.0 * db / bins
    longitude = collect(range(-dl / 2.0, dl / 2.0, length=db))

    # CREATE FIGURE
    fig, p = triple_panels()
    p.left.ylabel = L"ybin number $$"
    p.left.xlabel = L"intensity $$"
    p.bottom.xlabel = L"longitude ($^\circ$)"

    # PLOTTING DATA
    lines!(p.left, intensity, pulses, color=:grey, linewidth=0.5)
    #xlims!(p.left, [0.01, 1.01])
    ylims!(p.left, [pulses[1]-0.5, pulses[end]+0.5])

    heatmap!(p.center, transpose(da))

    lines!(p.bottom, longitude, average, color=:grey, linewidth=0.5)
    xlims!(p.bottom, [longitude[1], longitude[end]])

    #screen = display(fig)
    #resize!(screen, 500, 800)
    filename = "$outdir/$(name_mod)_folded.pdf"
    println(filename)
    save(filename, fig, pt_per_unit=1)
end


function test_fft(data)

    #da = transpose(data)
    da = data
    #da = vcat(da[120, :], da[121, :], da[122, :])
    da = da[1:256, :]
    res = []
    for i in 25:100
        res = vcat(res, da[i, :])
        #println(size(res))
    end
    da = convert(Array{Float64}, res)

    # Fourier Transform 
    sz = size(da, 1)
    half = floor(Int, sz / 2)
    ff0 = fft(da)
    ff = abs.(fft(da))
    #ff = real.(fft(da))
    freq = fftfreq(sz)
    #iff = real.(ifft(ff0))

    #=
    ff_new = []
    st = 3
    en = 10
    for i in 1:sz
        #push!(ff_new, 0)
        if (freq[i] >= 0.5 * st/1024) && (freq[i] <= 0.5 * en/1024)
            push!(ff_new, 0)
        elseif freq[i] < 0
            push!(ff_new, 0)
        else
            push!(ff_new, ff0[i])
        end
    end
    ff_new = convert(Array{ComplexF64}, ff_new)
    #println(typeof(ff0))
    iff = real.(ifft(ff_new))
    =#
    iff = real.(ifft(ff0))

    # Autocorrelation
    lags = collect(1:2048)
    r = autocor(da, lags) #; demean=false)

    # periodogram
    p = periodogram(da)
    #p = welch_pgram(da)


    fig = Figure()

    ax = Axis(fig[1, 1])
    lines!(ax, da, color=:black, linewidth=2)
    lines!(ax, iff, color=(:red, 0.4), linewidth=5)

    ax2 = Axis(fig[2, 1])
    #lines!(ax2, ff, color=:red, linewidth=1)
    #lines!(ax2, freq, ff, color=:red, linewidth=1)
    lines!(ax2, freq[1:half], ff[1:half], color=:red, linewidth=1)
    #lines!(ax2, freq[1:half], abs.(ff_new)[1:half], color=:green, linewidth=1)
    #lines!(ax2, freq, abs.(ff_new), color=:green, linewidth=1)
    #lines!(ax2, freq, ff, color=:red, linewidth=1)
    #lines!(ax2, freq, imag.(ff0), color=:blue, linewidth=1)

    vlines!(ax2, 1/1024, color=:blue, linewidth=1)
    #vlines!(ax2, 2/1024, color=:blue, linewidth=1)
    #vlines!(ax2, 3/1024, color=:blue, linewidth=1)
    #lines!(fig, 2*da, color=:grey, linewidth=0.5)

    ax3 = Axis(fig[3, 1])
    v, ind = findmax(r[3:end])
    lines!(ax3, r , color=:blue, linewidth=1)
    vlines!(ax3, ind+2, color=:red, linewidth=1)

    ax4 = Axis(fig[4, 1])
    #lines!(ax4, p.freq, DSP.pow2db.(p.power), color=:black, linewidth=1)
    lines!(ax4, p.freq, p.power, color=:black, linewidth=1)
    vlines!(ax4, 1/1024, color=(:red, 0.4), linewidth=4)

    #save("output/test.pdf", fig)
    display(fig)
end


end # module
