# --- Snapshot Read in --- #

## Struct to hold snapshot data ##
struct SnapDat
    id::Vector{Float64}
    startype::Vector{Float64}
    binflag::Vector{Float64}
    r::Vector{Float64}
    m::Vector{Float64}
    vr::Vector{Float64}
    vt::Vector{Float64}
    v::Vector{Float64}
    N::Int
    rmin::Float64
    rmax::Float64
    vmin::Float64
    vmax::Float64
    vrmin::Float64
    vrmax::Float64
    vtmin::Float64
    vtmax::Float64
end

## Read in snapshot data as struct ##
function read_file(filename::String)

    # Read data into DF
    filename = expanduser(filename)
    df = CSV.read(filename, DataFrame)

    # Extract columns
    id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag_tab = eachcol(df[:, 1:7])

    # Compute total velocity, N
    v_tab = sqrt.(vr_tab .^ 2 .+ vt_tab .^ 2)
    N = length(r_tab)

    # Compute extrema
    rmin, rmax = extrema(r_tab)
    vmin, vmax = extrema(v_tab)
    vrmin, vrmax = extrema(vr_tab)
    vtmin, vtmax = extrema(vt_tab)
    
    # Return struct
    return SnapDat(
        id_tab, startype_tab, binflag_tab,
        r_tab, m_tab, vr_tab, vt_tab, v_tab,
        N, rmin, rmax, vmin, vmax, vrmin, vrmax, vtmin, vtmax)
end

function read_file_h5(filename::String, dataset::Union{Nothing, String} = nothing,
                      conv_path::Union{Nothing, AbstractString} = nothing)

    filename = expanduser(filename)

    # If a .conv.sh path is given, set the model constants from it first (c, t_conv,
    # r_conv, timeunitsmyr, massunitmsun). Each model has its own massunitmsun -- the
    # M_sun-per-code-mass factor that converts particle masses -- so hard-coding it
    # (was 5.52429e6, the Elena value) silently corrupts masses for any other model.
    if conv_path !== nothing
        set_model_constants!(conv_path)
    end
    conversion = massunitmsun   # M_sun per code mass unit (from the active model constants)

    dat = h5open(filename, "r") do f

        # Select the snapshot dataset.
        dsname = dataset
        rows = read(f[dsname])

        have = fieldnames(eltype(rows))
        for fld in (:m_MSUN, :r, :vr, :vt, :startype, :binflag)
            fld in have || error(
                "read_file_h5: field \"$fld\" missing in \"$dsname\". Available: $have")
        end

        # Extract columns by name (order-independent, unlike the CSV path).
        m_tab        = Float64[row.m_MSUN  for row in rows] ./ conversion
        r_tab        = Float64[row.r       for row in rows]
        vr_tab       = Float64[row.vr      for row in rows]
        vt_tab       = Float64[row.vt      for row in rows]
        startype_tab = Float64[row.startype for row in rows]
        binflag_tab  = Float64[row.binflag for row in rows]

        (m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag_tab)
    end

    m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag_tab = dat

    N = length(r_tab)
    id_tab = collect(0.0:(N - 1))  # 0-based id, matching the pandas index

    # Compute total velocity, extrema (identical to read_file)
    v_tab = sqrt.(vr_tab .^ 2 .+ vt_tab .^ 2)

    rmin, rmax = extrema(r_tab)
    vmin, vmax = extrema(v_tab)
    vrmin, vrmax = extrema(vr_tab)
    vtmin, vtmax = extrema(vt_tab)

    return SnapDat(
        id_tab, startype_tab, binflag_tab,
        r_tab, m_tab, vr_tab, vt_tab, v_tab,
        N, rmin, rmax, vmin, vmax, vrmin, vrmax, vtmin, vtmax)
end

## Set the model-dependent constants from a CMC .conv.sh file ##
#
# The .conv.sh file (e.g. output.conv.sh) is a list of `key=value` bash assignments.
# We read three of them and derive the code speed of light:
#
#   nbtimeunitsmyr   -> t_conv        N-body time unit [Myr]  (the MC / orbit clock)
#   lengthunitparsec -> r_conv        code length unit [pc]
#   timeunitsmyr     -> timeunitsmyr  code time unit [Myr]    (converts centmass.dat times)
#
#   c = c_phys / v_unit,   v_unit = (r_conv pc) / (t_conv Myr)   [the code velocity unit]
#
# G and M are always 1 and are NOT touched. M_bh is per-snapshot (from centmass.dat) and
# is likewise not set here. Reassigns the typed globals c, t_conv, r_conv, timeunitsmyr
# (see constants.jl) and prints ONE summary line. Call once per model import -- every
# snapshot of a given model shares one .conv.sh. Returns the constants as a NamedTuple.
function set_model_constants!(conv_path::AbstractString)

    conv_path = expanduser(conv_path)

    # Parse `key=value` lines (skip comments, blanks, and non-numeric right-hand sides
    # such as `outprefix=output`). Exact-key match avoids `timeunitsmyr` colliding with
    # the `nbtimeunitsmyr` / `timeunitcgs` substrings.
    vals = Dict{String,Float64}()
    for line in eachline(conv_path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        eq = findfirst('=', s)
        eq === nothing && continue
        key = strip(s[1:prevind(s, eq)])
        rhs = strip(s[nextind(s, eq):end])
        v = tryparse(Float64, rhs)
        v === nothing && continue
        vals[key] = v
    end

    for k in ("nbtimeunitsmyr", "lengthunitparsec", "timeunitsmyr",
              "lengthunitcgs", "nbtimeunitcgs", "massunitmsun")
        haskey(vals, k) || error("set_model_constants!: '$k' not found in $conv_path")
    end

    global t_conv       = vals["nbtimeunitsmyr"]     # N-body time unit [Myr]
    global r_conv       = vals["lengthunitparsec"]   # code length unit [pc]
    global timeunitsmyr = vals["timeunitsmyr"]       # code time unit [Myr] (centmass.dat)
    global massunitmsun = vals["massunitmsun"]       # code mass unit [M_sun] (particle m_MSUN -> code)

    # Code speed of light from the code velocity unit, computed straight from the file's
    # own cgs quantities so we inherit exactly CMC's pc/year definitions (no locally
    # assumed conversion constants). The velocity unit is the DYNAMICAL one:
    # length / N-body time = lengthunitcgs / nbtimeunitcgs -- the same time unit as
    # t_conv = nbtimeunitsmyr. NOT timeunitcgs (the code/FP time unit, ~1e6x larger,
    # which would give a nonsensical c ~ 1e9). c_phys is therefore in cgs (cm/s).
    c_phys_cgs = 2.99792458e10                                   # cm/s
    v_unit_cgs = vals["lengthunitcgs"] / vals["nbtimeunitcgs"]   # code velocity unit [cm/s]
    global c   = c_phys_cgs / v_unit_cgs

    println(">>> Set model constants from conv.sh: $conv_path")
    println(">>>   c = $c,  t_conv = $t_conv Myr/nbody,  r_conv = $r_conv pc/code,  " *
            "timeunitsmyr = $timeunitsmyr Myr/code,  massunitmsun = $massunitmsun Msun/code,  " *
            "v_unit = $(round(v_unit_cgs/1e5, digits=3)) km/s")
    flush(stdout)

    return (c = c, t_conv = t_conv, r_conv = r_conv,
            timeunitsmyr = timeunitsmyr, massunitmsun = massunitmsun)
end

## Snapshot keys spaced ~`sep` Gyr apart, in time order ##
#
# Keys look like "<index>(t=<time>Gyr)"; the index prefix is NOT time-ordered, so the
# time is parsed out and used for sorting. Starting from `sep`, the first snapshot at or
# after each target (sep, 2*sep, ...) is selected, giving snapshots spaced >= `sep` Gyr.
# Returns the selected keys as Strings in increasing time order. `sep` defaults to
# 0.1 Gyr (100 Myr). Feed each returned key to read_file_h5(filename, key).
function snapshot_keys_by_time(filename::String; sep::Real = 0.1)

    filename = expanduser(filename)

    all_keys = h5open(filename, "r") do f
        collect(keys(f))
    end

    # Parse the time (Gyr) out of each key; skip non-matching keys and the t = 0
    # snapshot (always excluded).
    times     = Float64[]
    keys_kept = String[]
    for key in all_keys
        m = match(r"t=([-+0-9.eE]+)\s*Gyr", key)
        m === nothing && continue
        t = parse(Float64, m.captures[1])
        t > 0 || continue
        push!(times, t)
        push!(keys_kept, key)
    end

    isempty(times) && return String[]

    # Sort by time (the numeric index prefix is unordered).
    ord       = sortperm(times)
    times     = times[ord]
    keys_kept = keys_kept[ord]

    # Fixed thresholds k*sep (k = 1, 2, ...): select the first snapshot at/after each,
    # deduplicating when snapshots are sparser than `sep`. Fixed thresholds (rather than
    # chaining off the last picked time) keep every ~sep-spaced snapshot even when the
    # times jitter around round values, e.g. 0.10007, 0.20002, 0.40001.
    tol      = 1e-6 * sep
    selected = String[]
    tmax     = times[end]
    last_idx = 0
    k        = 1
    while k * sep <= tmax + tol
        idx = searchsortedfirst(times, k * sep - tol)   # first time >= k*sep
        if idx <= length(times) && idx != last_idx
            push!(selected, keys_kept[idx])
            last_idx = idx
        end
        k += 1
    end

    return selected
end
