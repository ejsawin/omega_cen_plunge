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
