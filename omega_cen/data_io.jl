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
        N, rmin, rmax, vmin, vmax, vrmin, vrmax, vtmin, vtmax
    )
end
