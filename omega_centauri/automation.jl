# --- Automated snapshot -> diffusion-coefficient pipeline --- #
#
# Requires the functions/globals set up by main.jl (read_file_h5, find_psi_arrays,
# psi_exact, orbit_sampled_kde_rv_df, activate_orbit_psi!, save_orbit_psi,
# generate_coeffs, save_coeffs, snapshot_keys_by_time). Include main.jl first.

using Printf

## Filename time tag "1dot00" from a Gyr value (2 decimals, "." -> "dot") ##
_time_tag(t_gyr) = replace(@sprintf("%.2f", t_gyr), "." => "dot")

## Parse the Gyr time out of a snapshot key "<idx>(t=<time>Gyr)" ##
function _snapshot_time(key::AbstractString)
    m = match(r"t=([-+0-9.eE]+)\s*Gyr", key)
    m === nothing && error("automate_dc: could not parse a time from key \"$key\"")
    return parse(Float64, m.captures[1])
end

## IMBH mass (code units) at a given time from a CMC centmass.dat file ##
#
# The file has a header line "#1:t #2:cenma.m #3:Dt ..." then whitespace-separated rows,
# with the MBH mass in column `cenma.m` and time in column `t`, both in CODE units.
# `target_time_myr` is in Myr; `myr_per_code` is the Myr per code time unit (e.g. t_conv),
# used to convert it to code time before interpolating cenma.m linearly (nearest-edge
# outside the recorded range).
function get_MBH_mass(centmass_path::AbstractString, target_time_myr::Real, myr_per_code::Real)

    centmass_path = expanduser(centmass_path)

    t_idx = 0; m_idx = 0
    times = Float64[]; masses = Float64[]

    for line in eachline(centmass_path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#")
            if occursin(r"#\d+:", s)
                names = [String(mm.captures[1]) for mm in eachmatch(r"#\d+:(\S+)", s)]
                ti = findfirst(==("t"), names);        t_idx = ti === nothing ? 0 : ti
                mi = findfirst(==("cenma.m"), names);  m_idx = mi === nothing ? 0 : mi
            end
            continue
        end
        (t_idx == 0 || m_idx == 0) &&
            error("get_MBH_mass: 't'/'cenma.m' columns not found before data in $centmass_path")
        f = split(s)
        push!(times,  parse(Float64, f[t_idx]))
        push!(masses, parse(Float64, f[m_idx]))
    end

    isempty(times) && error("get_MBH_mass: no numeric rows in $centmass_path")

    target_code = target_time_myr / myr_per_code   # Myr -> code time
    return _interp1(target_code, times, masses)     # linear interp, nearest at the ends
end

## Process every snapshot spaced ~`sep` Gyr apart: build the orbit-sampled potential
## and diffusion coefficients for each, saving into `out_folder`:
##   <prefix>_<time>Gyr.h5          -- diffusion coefficients (grid_size x grid_size)
##   <prefix>_<time>Gyr_rv_psi.h5   -- the matching orbit-sampled potential
## and one shared record file
##   <prefix>_lowest101_potentials.h5  -- ETab[1:101] (lowest-101 particle energies,
##                                        snapshot potential) and the time, per snapshot.
## `time` is formatted with 2 decimals and a "dot" (e.g. 1.00 Gyr -> "IMBH01_1dot00Gyr").
## Reassigns the same globals main.jl uses; run it fresh (e.g. via run_automated_dc.jl).
function automate_dc(filename::AbstractString, sep::Real, out_folder::AbstractString,
                     name_prefix::AbstractString; centmass_path::AbstractString,
                     myr_per_code::Real, n_rv::Integer = 1000, grid_size::Integer = 200)

    global M_bh, dat, psi_tab, psi_tot_tab, M_tab, psi_rtab, psi_Mtot, psi_calc, ETab, DF

    filename = expanduser(filename)
    isdir(out_folder) || mkpath(out_folder)

    snap_keys = snapshot_keys_by_time(filename; sep = sep)
    nsnap = length(snap_keys)
    println("automate_dc: $nsnap snapshots (sep = $sep Gyr) from $filename")
    println("  output -> $out_folder,  prefix = $name_prefix,  n_rv = $n_rv,  grid = $grid_size")
    flush(stdout)

    potential_file = joinpath(out_folder, "$(name_prefix)_lowest101_potentials.h5")

    h5open(potential_file, "w") do pf
        for (i, key) in enumerate(snap_keys)
            t_gyr = _snapshot_time(key)
            tag   = _time_tag(t_gyr)

            println("[$i/$nsnap] $key  (t = $t_gyr Gyr)  ->  $(name_prefix)_$(tag)Gyr")
            flush(stdout)

            try
                # Update the IMBH mass for this snapshot time (potential + coeffs use M_bh)
                M_bh = get_MBH_mass(centmass_path, t_gyr * 1000.0, myr_per_code)
                println("    M_bh = $M_bh (code units)")
                flush(stdout)

                # Snapshot + its (Henon) potential
                dat = read_file_h5(filename, key)
                psi_tab, psi_tot_tab, M_tab = find_psi_arrays(dat.r, dat.m)
                psi_rtab = dat.r
                psi_Mtot = sum(dat.m)
                psi_calc = psi_exact(psi_rtab, psi_tab, M_tab)
                ETab = sort(0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r))

                # Record the lowest-101 particle energies (snapshot potential) + time + M_bh
                nlow = min(101, length(ETab))
                g = create_group(pf, "$(tag)Gyr")
                g["ETab_low101"] = ETab[1:nlow]
                attrs(g)["time"] = t_gyr
                attrs(g)["key"]  = String(key)
                attrs(g)["M_bh"] = M_bh
                flush(pf)

                # Orbit-sampled DF + potential, then switch the active potential to it
                DF = orbit_sampled_kde_rv_df(dat, n_rv)
                activate_orbit_psi!()

                # Save the orbit-sampled potential and the diffusion coefficients
                coef_file = joinpath(out_folder, "$(name_prefix)_$(tag)Gyr.h5")
                psi_file  = joinpath(out_folder, "$(name_prefix)_$(tag)Gyr_rv_psi.h5")
                save_orbit_psi(psi_file)
                coef_new = generate_coeffs(grid_size)
                save_coeffs(coef_file, coef_new)

                println("    saved $(basename(coef_file)) + $(basename(psi_file))")
                flush(stdout)
            catch err
                err isa InterruptException && rethrow()
                println(stderr, "automate_dc: FAILED on $key (t = $t_gyr Gyr): $err")
                flush(stderr)
            end
        end
    end

    println("automate_dc: complete ($nsnap snapshots); potentials in $(basename(potential_file))")
    return potential_file
end


# --- Power-law extrapolation of the diffusion coefficients below the data --- #
#
# The coefficient grid only reaches emin (the most-bound particle); below that our
# particles are discrete so there is no data. These helpers extend the grid down to a
# fixed deep energy by fitting a power law in |E| per j column and extrapolating.

## Linear interpolation with nearest-value extrapolation at the ends (like np.interp) ##
function _interp1(x::Real, xs::Vector{Float64}, ys::Vector{Float64})
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs, x)
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + t * (ys[k+1] - ys[k])
end

## Replace NaN/Inf entries by linear interpolation across the valid ones ##
function _fill_invalid(v::Vector{Float64})
    valid = isfinite.(v)
    s = count(valid)
    s == 0 && error("_fill_invalid: no valid values to interpolate from")
    s == 1 && return fill(v[findfirst(valid)], length(v))
    idx = collect(1:length(v))
    xs  = Float64.(idx[valid]); ys = v[valid]
    out = copy(v)
    for i in idx[.!valid]
        out[i] = _interp1(Float64(i), xs, ys)
    end
    return out
end

## 1D Gaussian smoothing with "nearest" boundary (like scipy gaussian_filter1d) ##
function _gaussian_smooth(v::Vector{Float64}, sigma::Real)
    sigma <= 0 && return copy(v)
    n = length(v)
    radius = max(1, round(Int, 4 * sigma))
    kernel = [exp(-0.5 * (i / sigma)^2) for i in -radius:radius]
    kernel ./= sum(kernel)
    out = similar(v)
    @inbounds for i in 1:n
        acc = 0.0
        for (ki, k) in enumerate(-radius:radius)
            acc += kernel[ki] * v[clamp(i + k, 1, n)]
        end
        out[i] = acc
    end
    return out
end

## Extend a DiffusionCoeffs grid in E via a per-j power-law fit + extrapolation ##
#
# Keeps rows with E >= E_fit_deep; below that (down to E_end) every coefficient matrix
# is refilled from a power law fit over E_fit_deep <= E <= E_fit_shallow (independently
# per j column, then the slope/amplitude smoothed across j). E_fit_deep/E_fit_shallow are
# energies (e.g. the 2nd and 11th lowest particle energies). The new grid has n_E_total
# rows in E (j unchanged), giving finer resolution in the extrapolated region. Returns a
# new DiffusionCoeffs; nothing is written to disk.
function extrapolate_coeffs(coef::DiffusionCoeffs, E_fit_deep::Real, E_fit_shallow::Real;
                            E_end::Real = -1.5e6, n_E_total::Integer = 300,
                            smooth_sigma::Real = 2.0, value_floor::Real = 1e-15,
                            min_fit_points::Integer = 3)

    E_tab = collect(float.(coef.E_tab))     # ascending: emin ... ~-0.1
    j_tab = collect(float.(coef.j_tab))
    n_j   = length(j_tab)

    E_end < E_fit_deep      || error("extrapolate_coeffs: E_end ($E_end) must be below E_fit_deep ($E_fit_deep)")
    E_fit_deep < E_fit_shallow || error("extrapolate_coeffs: E_fit_deep must be below E_fit_shallow")

    # Rows to keep (reliable data) and the fit interval.
    keep_mask = E_tab .>= E_fit_deep
    E_keep    = E_tab[keep_mask]
    n_keep    = length(E_keep)
    n_keep > 0 || error("extrapolate_coeffs: no grid rows at/above E_fit_deep = $E_fit_deep")

    n_tail = n_E_total - n_keep
    n_tail >= 1 || error("extrapolate_coeffs: n_E_total ($n_E_total) <= kept rows ($n_keep)")

    fit_mask = (E_tab .>= E_fit_deep) .& (E_tab .<= E_fit_shallow)
    n_fit    = count(fit_mask)
    n_fit >= min_fit_points ||
        error("extrapolate_coeffs: only $n_fit grid rows in fit interval [$E_fit_deep, $E_fit_shallow]")

    # New deep tail: geomspace(E_end, E_keep[1]) excluding the E_keep[1] endpoint.
    E_tail = -exp10.(range(log10(abs(E_end)), log10(abs(E_keep[1])), length = n_tail + 1))[1:end-1]
    E_new  = vcat(E_tail, E_keep)           # ascending, length n_E_total

    logE_fit  = log10.(abs.(E_tab[fit_mask]))
    logE_tail = log10.(abs.(E_tail))
    fit_cols  = findall(fit_mask)

    mats = (coef.DE1_tab, coef.DE2_tab, coef.Dj1_tab, coef.Dj2_tab,
            coef.DEE_tab, coef.Djj_tab, coef.DEj_tab)
    out = Matrix{Float64}[]

    for M in mats                            # M is [n_j x n_E]
        Mn = Matrix{Float64}(undef, n_j, n_E_total)
        Mn[:, (n_tail + 1):end] = M[:, keep_mask]

        slopes = fill(NaN, n_j); inters = fill(NaN, n_j); signs = fill(NaN, n_j)

        for r in 1:n_j                       # fit each j row over the fit columns
            vals  = Float64[M[r, c] for c in fit_cols]
            valid = isfinite.(vals) .& (abs.(vals) .> value_floor)
            count(valid) >= min_fit_points || continue
            vv = vals[valid]; le = logE_fit[valid]
            dom  = count(>(0), vv) >= count(<(0), vv) ? 1.0 : -1.0
            same = sign.(vv) .== dom
            count(same) >= min_fit_points || continue
            A   = hcat(le[same], ones(count(same)))
            sol = A \ log10.(abs.(vv[same]))
            slopes[r] = sol[1]; inters[r] = sol[2]; signs[r] = dom
        end

        count(isfinite, slopes) > 0 || error("extrapolate_coeffs: power-law fit failed for every j")

        slopes = _fill_invalid(slopes); inters = _fill_invalid(inters)
        signs  = ifelse.(_fill_invalid(signs) .>= 0, 1.0, -1.0)
        slopes = _gaussian_smooth(slopes, smooth_sigma)
        inters = _gaussian_smooth(inters, smooth_sigma)

        for r in 1:n_j                       # evaluate the smoothed law on the tail
            Mn[r, 1:n_tail] = signs[r] .* exp10.(inters[r] .+ slopes[r] .* logE_tail)
        end
        all(isfinite, Mn) || error("extrapolate_coeffs: non-finite tail values")
        push!(out, Mn)
    end

    DE1, DE2, Dj1, Dj2, DEE, Djj, DEj = out
    interp(Z) = (j, E) -> bilinear_interp_clamped(j, E, j_tab, E_new, Z)
    return DiffusionCoeffs(DE1, DE2, Dj1, Dj2, DEE, Djj, DEj,
        E_new, j_tab, E_new[1], coef.emax,
        interp(DE1), interp(DE2), interp(Dj1), interp(Dj2),
        interp(DEE), interp(Djj), interp(DEj))
end


## Run the Monte Carlo for every snapshot, using per-snapshot coeffs/psi from automate_dc ##
#
# For each snapshot (spaced ~`sep` Gyr): load <prefix>_<time>Gyr.h5 (coeffs) and
# <prefix>_<time>Gyr_rv_psi.h5 (potential) from `dc_folder`, optionally extrapolate the
# coefficients below the data (method = 2, default), seed one walker per bound BH/NS from
# `snapshot_file`, integrate to (max_t, max_steps), and store the FINAL state per walker.
# Output: <dc_folder>/<prefix>_mc.h5, one group per snapshot.
#   method = 1 : no coefficient extrapolation (assume no diffusion below the data)
#   method = 2 : power-law extrapolate over [ETab[level_low], ETab[level_high]] down to E_end
# max_t is in Myr (run_mc_rp_ra now converts to code units internally).
function automate_mc(snapshot_file::AbstractString, dc_folder::AbstractString, sep::Real,
                     name_prefix::AbstractString, max_t::Real, max_steps::Integer;
                     centmass_path::AbstractString, myr_per_code::Real,
                     method::Integer = 2, level_low::Integer = 2, level_high::Integer = 11,
                     E_end::Real = -1.5e6, n_E_total::Integer = 300, smooth_sigma::Real = 2.0)

    global M_bh

    snapshot_file = expanduser(snapshot_file)
    dc_folder     = expanduser(dc_folder)

    snap_keys = snapshot_keys_by_time(snapshot_file; sep = sep)
    nsnap = length(snap_keys)
    out_file = joinpath(dc_folder, "$(name_prefix)_mc.h5")

    println("automate_mc: $nsnap snapshots (sep = $sep Gyr), method = $method")
    println("  snapshots: $snapshot_file")
    println("  coeffs/psi: $dc_folder,  output: $out_file")
    flush(stdout)

    h5open(out_file, "w") do of
        attrs(of)["snapshot_file"] = snapshot_file
        attrs(of)["sep"]        = float(sep)
        attrs(of)["max_t"]      = float(max_t)
        attrs(of)["max_steps"]  = max_steps
        attrs(of)["method"]     = method
        attrs(of)["level_low"]  = level_low
        attrs(of)["level_high"] = level_high
        attrs(of)["E_end"]      = float(E_end)
        attrs(of)["n_E_total"]  = n_E_total
        attrs(of)["nthreads"]   = Threads.nthreads()  # M_bh is per-snapshot (see each group)

        for (i, key) in enumerate(snap_keys)
            t_gyr = _snapshot_time(key)
            tag   = _time_tag(t_gyr)
            coef_file = joinpath(dc_folder, "$(name_prefix)_$(tag)Gyr.h5")
            psi_file  = joinpath(dc_folder, "$(name_prefix)_$(tag)Gyr_rv_psi.h5")

            println("[$i/$nsnap] $key  (t = $t_gyr Gyr)")
            flush(stdout)

            if !isfile(coef_file) || !isfile(psi_file)
                println(stderr, "automate_mc: missing coeff/psi for $tag; skipping.")
                flush(stderr)
                continue
            end

            try
                # IMBH mass for this snapshot time -- must match what automate_dc used, and
                # feeds loss_cone / gw_rates / find_Lc and the psi IMBH term below.
                M_bh = get_MBH_mass(centmass_path, t_gyr * 1000.0, myr_per_code)
                println("    M_bh = $M_bh (code units)")
                flush(stdout)

                # Match the potential + coefficients this snapshot was built with.
                coef = load_coeffs(coef_file)
                load_orbit_psi!(psi_file)
                activate_orbit_psi!()
                dat_s = read_file_h5(snapshot_file, key)

                # Optional power-law extrapolation of the coefficients below the data.
                if method == 2
                    ETab_s = sort(0.5 .* (dat_s.vr .^ 2 .+ dat_s.vt .^ 2) .+ psi_calc.(dat_s.r))
                    coef = extrapolate_coeffs(coef, ETab_s[level_low], ETab_s[level_high];
                        E_end = E_end, n_E_total = n_E_total, smooth_sigma = smooth_sigma)
                end

                # One walker per bound BH/NS; integrate and keep the final state.
                ics = black_hole_ics(dat_s)
                N   = length(ics.indices)
                Ef  = Vector{Float64}(undef, N)
                jf  = Vector{Float64}(undef, N)
                rpf = Vector{Float64}(undef, N)
                raf = Vector{Float64}(undef, N)
                reasons = Vector{Int}(undef, N)

                Threads.@threads :dynamic for k in 1:N
                    ts, es, js, rps, ras, reason =
                        run_mc_rp_ra(ics.E0[k], ics.J0[k], coef, ics.m0[k];
                                     max_step = max_steps, max_t = max_t)
                    Ef[k] = es[end]; jf[k] = js[end]
                    rpf[k] = rps[end]; raf[k] = ras[end]
                    reasons[k] = reason
                end

                g = create_group(of, "$(tag)Gyr")
                attrs(g)["time"] = t_gyr
                attrs(g)["M_bh"] = M_bh
                attrs(g)["N"]    = N
                g["indices"] = ics.indices
                g["E0s"] = ics.E0;  g["J0s"] = ics.J0;  g["m0s"] = ics.m0
                g["Ef"]  = Ef;      g["jf"]  = jf
                g["rpf"] = rpf;     g["raf"] = raf
                g["reason_to_end"] = reasons
                flush(of)

                println("    $N walkers integrated")
                flush(stdout)
            catch err
                err isa InterruptException && rethrow()
                println(stderr, "automate_mc: FAILED on $key (t = $t_gyr Gyr): $err")
                flush(stderr)
            end
        end
    end

    println("automate_mc: complete -> $(basename(out_file))")
    return out_file
end
