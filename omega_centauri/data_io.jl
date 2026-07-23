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

function read_file_h5(filename::String, dataset::Union{Nothing, String} = nothing)

    filename = expanduser(filename)
    conversion = 5.52429e6

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
