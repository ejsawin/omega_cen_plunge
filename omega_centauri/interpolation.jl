# --- Necessary functions for interpolation, integration --- #

## Bilinear interpolation (weighted means) ##
function bilinear_interp(x, y, x_grid, y_grid, f_grid)
    nx = length(x_grid)
    ny = length(y_grid)

    # Explicit bounds check -> Only sample within [rmin, rmax], [vmin, vmax]
    if x < x_grid[1] || x > x_grid[end] || y < y_grid[1] || y > y_grid[end]
        return 0.0
    end

    # Find index of left/bottom center
    ix = searchsortedlast(x_grid, x)
    iy = searchsortedlast(y_grid, y)

    # Check index bounds
    if ix < 1 || iy < 1 || ix >= nx || iy >= ny
        return 0.0
    end

    @inbounds begin # safe
        x1, x2 = x_grid[ix], x_grid[ix + 1]
        y1, y2 = y_grid[iy], y_grid[iy + 1]

        f11 = f_grid[ix, iy]
        f12 = f_grid[ix, iy + 1]
        f21 = f_grid[ix + 1, iy]
        f22 = f_grid[ix + 1, iy + 1]

        denom = (x2 - x1) * (y2 - y1)

        return ((x2 - x) * (y2 - y) * f11 +
                (x2 - x) * (y - y1) * f12 +
                (x - x1) * (y2 - y) * f21 +
                (x - x1) * (y - y1) * f22) / denom
    end
end

## Bilinear interpolation with nearest-edge extrapolation ##
# Out-of-bounds queries are clamped to the grid, so they take the closest edge value
# instead of 0. Used for the (E,j) diffusion-coefficient samplers in mc_step; the DF
# samplers keep the zero-outside bilinear_interp. Assumes ascending grids.
function bilinear_interp_clamped(x, y, x_grid, y_grid, f_grid)
    nx = length(x_grid)
    ny = length(y_grid)

    # Clamp the query into the grid (nearest-edge extrapolation)
    xc = clamp(x, x_grid[1], x_grid[end])
    yc = clamp(y, y_grid[1], y_grid[end])

    # Bracket indices, clamped to a valid interior cell so the exact upper edge works
    ix = clamp(searchsortedlast(x_grid, xc), 1, nx - 1)
    iy = clamp(searchsortedlast(y_grid, yc), 1, ny - 1)

    @inbounds begin
        x1, x2 = x_grid[ix], x_grid[ix + 1]
        y1, y2 = y_grid[iy], y_grid[iy + 1]

        f11 = f_grid[ix, iy]
        f12 = f_grid[ix, iy + 1]
        f21 = f_grid[ix + 1, iy]
        f22 = f_grid[ix + 1, iy + 1]

        denom = (x2 - x1) * (y2 - y1)

        return ((x2 - xc) * (y2 - yc) * f11 +
                (x2 - xc) * (yc - y1) * f12 +
                (xc - x1) * (y2 - yc) * f21 +
                (xc - x1) * (yc - y1) * f22) / denom
    end
end

## Bilinear interp: clamp x to the grid (nearest-edge), but return 0 if y is off its grid. ##
# For the (E,j) diffusion-coefficient samplers the closure is itp(j, E, j_tab, E_tab, Z), so
# x = j and y = E: clamp in j (a very small / near-radial j takes the nearest j-edge value, so
# a walker there keeps diffusing toward the loss cone instead of freezing) but return 0 when E
# leaves the data range (method 1's "no diffusion beyond the E grid"). Assumes ascending grids.
function bilinear_interp_clampj_zeroE(x, y, x_grid, y_grid, f_grid)
    ny = length(y_grid)
    if y < y_grid[1] || y > y_grid[end]        # E off the grid -> no diffusion
        return 0.0
    end
    nx = length(x_grid)
    xc = clamp(x, x_grid[1], x_grid[end])      # clamp j (nearest-edge)
    ix = clamp(searchsortedlast(x_grid, xc), 1, nx - 1)
    iy = clamp(searchsortedlast(y_grid, y),  1, ny - 1)

    @inbounds begin
        x1, x2 = x_grid[ix], x_grid[ix + 1]
        y1, y2 = y_grid[iy], y_grid[iy + 1]

        f11 = f_grid[ix, iy]
        f12 = f_grid[ix, iy + 1]
        f21 = f_grid[ix + 1, iy]
        f22 = f_grid[ix + 1, iy + 1]

        denom = (x2 - x1) * (y2 - y1)

        return ((x2 - xc) * (y2 - y) * f11 +
                (x2 - xc) * (y - y1) * f12 +
                (xc - x1) * (y2 - y) * f21 +
                (xc - x1) * (y - y1) * f22) / denom
    end
end

## 1D midpoint integration ##
function midpoint(func,x_start,x_end,steps)

    delx=(x_end-x_start)/steps
    integral=0.0

    for i in 0:steps-1

        x0=x_start+i*delx # Define subinterval
        x1=x0+delx

        mid=(x0+x1)/2 # Find midpoint

        integral+=func(mid)*delx
    end

    return integral
end


## 1D trapezoidal method ##
function trapezoid(func, x_start, x_end, steps)

    delx = (x_end - x_start) / steps
    integral = 0.0

    for i in 0:steps-1

        # Define subinterval
        x0 = x_start + i * delx
        x1 = x0 + delx

        # Avg func vals at endpoints
        integral += (func(x0) + func(x1)) / 2 * delx
    end
    return integral
end

## 2D midpoint integration ##
function integrate_2d(func, xmin, xmax, ymin, ymax; nx=1000, ny=1000)
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    integral = 0.0

    for i in 0:(nx - 1)
        x0 = xmin + i * dx
        x1 = x0 + dx
        xmid = (x0 + x1) / 2

        for j in 0:(ny - 1)
            y0 = ymin + j * dy
            y1 = y0 + dy
            ymid = (y0 + y1) / 2

            integral += func(xmid, ymid) * dx * dy
        end
    end

    return integral
end
