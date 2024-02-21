module Tools
    using CairoMakie # 2D vector plots
    CairoMakie.activate!()
    using Statistics
    using LsqFit
    using FFTW
    using Peaks


    """
    Calculates RMS
    """
    function rms(data)
        s = 0.0
        for a in data
            s += a * a
        end
        return sqrt(s / length(data))
    end


    """
    Calculates average profile
    """
    function average_profile(data; norm=true)
        pulses, bins = size(data)
        ave = zeros(bins)
        for i in 1:pulses
            for j in 1:bins
                ave[j] += data[i,j]
            end
        end
        if norm == true
            ma = maximum(ave)
            ave = ave / ma
        elseif length(norm) == 2
            me = mean(ave[norm[1]:norm[2]])
            ave .-= me
            ma = maximum(ave)
            ave = ave ./ ma
        end
        return ave
    end


    """
    Calculates intensity of single pulses
    """
    function pulses_intensity(data)
        pulse_num, bins = size(data)
        intensity = zeros(pulse_num)
        for i in 1:pulse_num
            intensity[i] = sum(data[i, :])
        end
        pulses = collect(1:pulse_num)
        mi = minimum(intensity)
        if mi < 0
            intensity .+= -mi
        end
        intensity /= maximum(intensity)
        return (intensity, pulses)
    end


    """
    Finds peaks in 1D data
    """
    function peaks(intensity)
        y = intensity
        pks, vals = findmaxima(y)
        pks, proms = peakproms(pks, y)
        maxprom, ind = findmax(proms)
        return pks[ind] # index of peak with maximum prominence
    end


    """
    Calculates LRFS
    """
    function lrfs(data)
        da = transpose(data)
        bins, pulse_num = size(da)
        half = floor(Int, pulse_num / 2) # one side frequency range?
        lrfs = fft(da, 2)[:, 1:half] # second dim! important!
        lrfs = transpose(lrfs)
        freq = fftfreq(pulse_num)[1:half]  # the [1:half] added. fftfreq changed?
        intensity = zeros(half)
        ab = abs.(lrfs)
        for i in 1:half
            intensity[i] = sum(ab[i,:]) # this is important!
        end
        pk = peaks(intensity)
        return lrfs, intensity, freq, pk
    end

    """
    Gaussian function
    """
    @. gauss(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]


    """
    Fits Gaussian
    """
    function fit_gaussian(xdata, ydata; a=nothing, μ=nothing, σ=nothing, baselevel=nothing)
        if baselevel === nothing
            baselevel = mean(ydata)
        end
        if a === nothing
            a = maximum(ydata) - baselevel
        end
        if μ === nothing
            ind = trunc(Int, length(xdata) / 2)
            μ = xdata[ind]
        end
        if σ === nothing
            σ = xdata[20] - xdata[1] # TODO fix this
        end
        p0 = [a, μ, σ, baselevel]  # look here
        fit = curve_fit(gauss, xdata, ydata, p0)
        p = coef(fit)
        err = stderror(fit)
        return p, err
    end


    """
    Folds single pulse data using P3 (simplest)
    """
    function p3fold(single_pulses, p3; ybins=9)
        pulses, bins = size(single_pulses)
        folded = zeros(ybins, bins)
        for i in 1:pulses
            frac = i / p3 - trunc(Int, i /p3)
            j = trunc(Int, frac * ybins) + 1 # start index = 1
            for k in 1:bins
                folded[j,k] += single_pulses[i,k]
            end
        end
        return folded
    end


    """
    Resample single pulse data for bootstarping 
    on_bins - on pulse bins
    """
    function resample(single_pulses, on_bins; off_bins=nothing, verbose=0)
        pulses, bins = size(single_pulses)

        if off_bins === nothing
            # choose random off regions
            dbins = on_bins[2] - on_bins[1]
            half = trunc(Int, dbins / 2 + 0.5) # +1 for odd numbers
            ranges = vcat(2*half : on_bins[1]-2*half, on_bins[2]+2*half : bins-2*half) # away from signal and borders
            mid = rand(ranges)
            off_bins = [mid-half, mid+half]
        end

        # pick random pulse
        pulse = rand(1:pulses)

        on_rms = Tools.rms(single_pulses[pulse, on_bins[1]:on_bins[2]])
        off_rms = Tools.rms(single_pulses[pulse, off_bins[1]:off_bins[2]])
        if verbose > 0
            println("on pulse $pulse $on_bins RMS ", on_rms)
            println("off pulse $pulse $off_bins RMS ", off_rms)
        end

        resampled_data =  Array{Float64}(undef, size(single_pulses))

        for i in 1:pulses
            for j in 1:bins
                noise =  (rand() - 0.5) * off_rms
                resampled_data[i, j] = single_pulses[i, j] + noise
            end
        end
        return resampled_data
    end

    
    """
    Fixes fft phase (strange approach but works for now)
    """
    function fix_fftphase!(phase)
        smo = length(phase) - 1
        # calculate phase differences
        dps = Array{Float64}(undef, smo)
        for i in 1:smo
            dp = abs(phase[i+1] - phase[i])
            dps[i] = dp
        end
        # finding breakpoints
        bps = []
        med = median(dps)
        for i in 1:smo
            if (dps[i] > 360 - 2*med) && (dps[i-1] < 2 * med) && (dps[i+1] < 2 * med)
                push!(bps, i)
            end
        end
        #println(bps)
        # changing phase
        for i in eachindex(phase)
            # find segment
            seg = nothing
            for j in eachindex(bps)
                if i <= bps[j]
                    seg = (j-1)
                    break
                end
            end
            if seg === nothing
                seg = length(bps)
            end
            phase[i] += seg * 360.
            #println("$i $seg")
        end
        #=
        me = mean(dps)
        med = median(dps)
        fig = Figure()
        ax = Axis(fig[1,1])
        scatter!(ax, dps)
        hlines!(ax, me, color=:blue)
        hlines!(ax, med, color=:red)
        hlines!(ax, 360 - 2*med, color=:red)
        vlines!(ax, bps, color=:red)
        display(fig)
        =#
    end


    """
    Gets FFT phase errors using bootstarping method
    num - number of lrfses to calculate
    bin_range - changes bin range of the resultant data
    """
    function phase_errors(single_pulses, on_bins; num=10, bin_range=nothing, fix_fftphase=false)

        phases = []
        raw_phases = []
        for i in 1:num
            data = resample(single_pulses, [on_bins[1], on_bins[2]], verbose=0)
            lrfs, intensity, freq, peak = Tools.lrfs(data)
            phase_ = rad2deg.(angle.(lrfs[peak, :]))  # fft phase variation 
            if bin_range !== nothing
                phase_ = phase_[bin_range[1]:bin_range[2]]
            end
            push!(raw_phases, copy(phase_))
            if fix_fftphase == true
                fix_fftphase!(phase_)
            end
            push!(phases, phase_)
        end

        phase = Array{Float64}(undef, length(phases[1]))
        ephase_plus = Array{Float64}(undef, length(phases[1]))        
        ephase_minus = Array{Float64}(undef, length(phases[1]))        

        # manual calculation # raw phases (for proper errors)
        for i in 1:length(phase)
            mean = 0
            mi, ma = 1e50, -1e50
            for j in 1:num
                mean += raw_phases[j][i]
                if raw_phases[j][i] < mi 
                    mi = raw_phases[j][i]
                end
                if phases[j][i] > ma 
                    ma = raw_phases[j][i]
                end
            end
            mean /= num
            #phase[i] = mean
            ephase_minus[i] = mean - mi
            ephase_plus[i] = ma - mean
        end
        # manual calculation # for phases
        for i in 1:length(phase)
            mean = 0
            for j in 1:num
                mean += phases[j][i]
            end
            mean /= num
            phase[i] = mean
        end

        #=
        fig = Figure()
        ax = Axis(fig[1,1])
        errorbars!(ax, collect(1:length(phase)), phase, ephase_minus, ephase_plus)
        #for ph in phases
        #    scatter!(ax, ph)
        #end
        display(fig)
        =#

        return phases, phase, ephase_minus, ephase_plus

    end

end # module
