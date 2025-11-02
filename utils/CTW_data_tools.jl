# version 003: Tools for the "diagonal" version.
# version 002: Just pick one barotropic mode
# version 001: Adds more diagnosis tools
#
module CTW_data_tools

export read_fort_data, sort_modes, calc_grid, save_modes, orthogonality,
  slope_vector, bottom_vector, map_to_xz, normfactor,
  remap_diag
#  remap_to_vertices

using Match
using NCDatasets
using OffsetArrays
using StrFormat


# == Read bin data ===================
"""
Read one record from a Fortran sequential binary file.
InOut:
  ist: input stream
Input (optional):
  mes: message to print out.
Output:
  buffer
"""
function read_seq_rec!(ist, buffer::AbstractArray; mes = nothing)
  header = read(ist, Int32)
  read!(ist, buffer)
  trailer = read(ist, Int32)
  #@assert header == trailer
  if header != trailer
    a = isnothing(mes) ? "" : "$(mes): "
    error("$(a)header = $(header), trailer = $(trailer)")
  end
  esz = sizeof(eltype(buffer))
  sz = header ÷ esz
  #(! isnothing(mes)) && println("$(mes): $(header) bytes")
  (! isnothing(mes)) && println("$(mes): $(sz) elements")
  return buffer
end

"Same as read_seq_rec! except that this one is for stream binary."
function read_stream_rec!(ist, buffer::AbstractArray; mes = nothing)
  read!(ist, buffer)
  sz = length(buffer)
  (! isnothing(mes)) && println("$(mes): $(sz) elements")
  return buffer
end


"""
Input:
  types: `DataType` or Tuple of `DataType`s. For example,
    types = Int32
    types = (Int32, Float64)
"""
function read_seq_rec(ist, types; mes = nothing)
  header = read(ist, Int32)
  ts = (isa(types, AbstractArray) || isa(types, Tuple)) ? types : (types,)
  buf = Any[]
  for t in ts # collection of DataTypes
    v = read(ist, t)
    push!(buf, v)
  end
  trailer = read(ist, Int32)
  (header != trailer) && error("header==$(header), trailer==$(trailer)")
  (! isnothing(mes)) && println("$(mes): $(header) bytes")
  if length(buf) == 1
    return buf[1]
  else
    return tuple(buf...)
  end
end

"Same as read_seq_rec except that this one is for stream binary."
function read_stream_rec(ist, types; mes = nothing)
  ts = (isa(types, AbstractArray) || isa(types, Tuple)) ? types : (types,)
  buf = Any[]
  sz = 0
  for t in ts # collection of DataTypes
    v = read(ist, t)
    push!(buf, v)
    sz += sizeof(t)
  end
  (! isnothing(mes)) && println("$(mes): $(sz) bytes")
  if length(buf) == 1
    return buf[1]
  else
    return tuple(buf...)
  end
end


"""
Read original binary data
 (N, im, km, num, alphaR, alphaI, pr, dx, dz, f0, bvf2e)
  and calculate kocn, iocn.

num[0:im+1, 0:km+1]
dx[0:im+1]
dz[0:km+1]

See ctwmodes-006.f90 for the format of the binary file.
The default floating-point type in the file is Float32
but there is a switch in ctwmodes_pars.f90 to change that.
Correspondingly, you need to use `intype = Float64` to read
a double precision fortran binary file.
"""
#function read_fort_data(fnam = "Eigenmodes-z.fort.bin") # intype = Float32
function read_fort_data(fnam; is_streambin = true)
  s = (is_streambin) ? "stream binary" : "Fortran sequential binary"
  println("Reading from $(fnam), $(s)")

  read_rec  = is_streambin ? read_stream_rec  : read_seq_rec
  read_rec! = is_streambin ? read_stream_rec! : read_seq_rec!

  ist = open(fnam, "r")

  (N, im, km) = read_rec(ist, (Int32,Int32,Int32))
  @show N, im, km

  num = OffsetArray(Array{Int32,2}(undef, im+2, km+2), 0:im+1, 0:km+1)
  kocn = OffsetArray(Array{Int32,1}(undef, im+2), 0:im+1)
  iocn = OffsetArray(Array{Int32,1}(undef, km+2), 0:km+1)

  read_rec!(ist, num;  mes="num")
  read_rec!(ist, kocn; mes="kocn")
  read_rec!(ist, iocn; mes="iocn")

  outtype = read_rec(ist, Int32)
  intype = @match outtype begin
    4 => Float32
    8 => Float64
    _ => error("unknown outtype = $(outtype)")
  end

  alphaR = Vector{intype}(undef, N)
  alphaI = Vector{intype}(undef, N)
  pr = Array{intype,2}(undef, N, N)
  dx = OffsetArray{intype,1}(undef, 0:im+1)
  dz = OffsetArray{intype,1}(undef, 0:km+1)
  # eitch = OffsetArray{Float32,1}(undef, 0:im+1)
  # dHdx  = Vector{Float32}(undef, im)
  bvf2e = OffsetArray{intype,1}(undef, 0:km)

  read_rec!(ist, alphaR; mes="alphaR")
  read_rec!(ist, alphaI; mes="alphaI")
  read_rec!(ist, pr;     mes="pr")

  read_rec!(ist, dx; mes="dx vector")
  read_rec!(ist, dz; mes="dz vector")
  # read_rec!(ist, eitch; mes="eitch(0:im+1)")
  # read_rec!(ist, dHdx;  mes="dHdx")
  f0 = read_rec(ist, intype; mes="f0")
  read_rec!(ist, bvf2e; mes="bvf2e(0:km)") # cell edges
  close(ist)

#  # kocn[i] is the last (deepmost) ocean point, cell center, for each i.
#  # iocn[k] is the first (westernmost) ocean point for each k.
#  kocn = Vector{Int32}(undef,im) # 1:im only
#  iocn = Vector{Int32}(undef,km) # 1:km only
#  for i = axes(kocn,1)
#    kocn[i] = findlast(ik -> ik > 0, num[i, 1:km])
#  end
#  for k = axes(iocn,1)
#    iocn[k] = findfirst(ik -> ik > 0, num[1:im, k])
#  end
#  #@show kocn
#  #@show iocn
  return (num=num, kocn=kocn, iocn=iocn
          ,alphaR=alphaR, alphaI=alphaI
          ,pr=pr, dx=dx, dz=dz
          ,f0=f0, bvf2e=bvf2e)
  #kocn, iocn)
end


# == eigenvalues ================

"""
Examine the imaginary parts (alphaI) of the eigenvalues.
We don't want imaginary eigenvalues.
"""
function check_all_real(; alphaI)
  allzero = true
  for n in axes(alphaI,1)
    if alphaI[n] != 0
      println(stderr, "alphaI[$(n)] = $(alphaI[n])")
      allzero = false
    end
  end
  if allzero
    println("all eigenvalues are real.")
  else
    println(stderr, "There are non-zero imaginary components as above.")
  end
end

"""
c_all[:]: unsorted c's in the original order.
Returns the index of barotropic mode.

First search for huge values.
If there are, they are barotropic modes.
If there aren't, then the largest value is the barotropic mode.
"""
function get_bt_mode(;c_all)
  c_abs = abs.(c_all)
  nsorted = sortperm(c_abs, rev=true) # larger to smaller |c|
  nns = findall(n -> c_abs[n] > 1.0e10, nsorted) # returns indices of indices
  ns = nsorted[nns]
  m = length(ns)
  if m == 0
    ns = [nsorted[1]]
  elseif m == 1
    n = ns[1]
    println("Barotropic: c[$(n)] = $(c_all[n])")
  else
    map(n->println(stderr,"Barotropic?: c[$(n)] = $(c_all[n])"), ns)
    if m == 2
      c1 = c_all[ns[1]]
      c2 = c_all[ns[2]]
      (c1*c2 > 0) && error("Same sign.")
    else
      println(stderr, "Are these barotropic modes?")
    end
  end
  return ns[1]
end


"""
Sort the modes into western (alphaR < 0)
 and eastern (alphaR > 0) modes  separately
 in the descending order of abs(c).

Input:
  `alphaR`, `alphaI`: real and imaginary parts of the eigenvalues.
  `inv` indicates whether `alphaR` is 1/c or is c.

Output:
  Extracts only western-boundary (alphaR < 0) and eastern-boundasry (alphaR > 0) modes
  separately. But we define cee > 0.
  Returns (cee, nidx), in such a way that
```
  cee[0], pr[:, nidx[0]] . . . mode 0
  cee[1], pr[:, nidx[1]] . . . mode 1
  cee[2], pr[:, nidx[2]] . . . mode 2
  . . .
```
That is,
`cee` is already ordered in the mode-number order.
`nidx` is the vector of indices into pr to extract the modes
   in the mode-number order.
"""
function sort_modes(alphaR, alphaI; inv = false)
  c_thresh = -1.0e-12 # unphysical modes
  check_all_real(; alphaI=alphaI)
  c_all = (inv) ? (1.0 ./ alphaR) : alphaR
  #nsbt = get_bt_modes(;c_all=c_all) # barotropic mode
  #@assert length(nsbt) == 1
  #nbt = nsbt[1]
  nbt = get_bt_mode(;c_all=c_all) # barotropic mode
  print("barotropic c = $(c_all[nbt])")
  if c_all[nbt] > 0 # sign is unstable for bt
    c_all[nbt] = - c_all[nbt] # We want western boundary
    print(" ... flipping sign ...")
  end
  println()
  idx = sortperm(c_all)
  i = findfirst(c -> c >= c_thresh, c_all[idx])
  idx_west = idx[1:(i-1)] # larger |c| to smaller |c|
  idx_east = reverse(idx[i:end]) # larger c to smaller c
  nax_west = 0:(length(idx_west)-1)
  nax_east = 0:(length(idx_east)-1)
  cee_west = OffsetArray(-c_all[idx_west], nax_west)
  cee_east = OffsetArray(c_all[idx_east], nax_east)
  nidx_west = OffsetArray(idx_west,nax_west)
  nidx_east = OffsetArray(idx_east, nax_east)
  (cee_west, nidx_west, cee_east, nidx_east)
end



# == map the eigenvectors to the gridded F(x,z) field. ======
"""
Map (expand) pr to the xz grid:
```
   n = nidx[m]
   ik = num[i,k]
   prij[i,j,m] = pr[ik,n]
```
where and m = 0, 1, 2, . . . is mode number.

prij[i,k,m] = NaN where it's not defined.

Note that nidx[m] == nothing and prij[:,:,m] = NaN if mode m is missing.
"""
function map_to_xz(; pr, nidx, num, outtype=Float64)
  #@show axes(nidx)
  #nmodes = length(nidx)
  #nax = 0:(nmodes-1)
  nax = axes(nidx,1)
  #iax = 0:(im+1)
  #kax = 0:(km+1)
  (iax, kax) = axes(num)
  #@show nax
  #@show iax
  #@show kax

  prij = fill(outtype(NaN), iax, kax, nax) # OffsetArray
  # @show typeof(prij)
  #@show axes(pr)
  #@show axes(prij)

  for (m, n) in pairs(nidx) # (mode number, original array index)
    isnothing(n) && continue # mode m is missing
    for k in kax # 0:(km+1)
      for i in iax # 0:(im+1)
        ik = num[i,k] # positive means defined.
        (ik > 0) && (prij[i,k,m] = pr[ik,n])
      end
    end
  end

  #>>>>
  #@show prij[end,end,1]
  #<<<<<
  return prij
end

"""
## This mapping is desinged for the old steppy slope configuration. ###

Map p[i,j] on to the vertices of gridboxes for plotting.
Along the slope, the vertex values are taken from the wall
 value (p[0,k] + p[iocn[k],k])/2, ignoring neiboring values
 in the x direction, because the vertical distances
 are much smaller.

Input:
  p2d[0:im+1,0:km+1] == prij[:,:,n] for an n.
  dx[0:im+1]
  dz[0:km+1]
Output:
  pv[0:im, 0:km]
  (iv=0 is to the east of i=0; iv=im to the east of i=im)
  (kv=0 is below k=0; kv=km below k=km)
"""
function remap_to_vertices(p2d; dx, dz, iocn, kocn)
  error("obsolete. missing isn't handled. ctwmodes-007 needs different remapping.")
  (iax,kax) = axes(p2d)
  im = iax[end-1]
  km = kax[end-1]
  @assert iax == 0:(im+1)
  @assert kax == 0:(km+1)
  @assert axes(dx,1) == iax
  @assert axes(dz,1) == kax
  ivax = 0:im
  kvax = 0:km

  #-- Intermediate gridded values:
  #  u[xxax,zax] ... value at the vertical edges
  #  v[xax,zzax] ... value at the horizontal edges
  u = fill(NaN, ivax, kax)
  v = fill(NaN, iax,  kvax)
  for k in 1:km # omit k=0,km+1
    i = iocn[k]
    u[i-1,k] = (p2d[0,k] + p2d[i,k])/2 # slope-wall
    for i in (iocn[k]+1):(im+1)
      xl = dx[i-1]; xr = dx[i]
      u[i-1,k] = (xr*p2d[i-1,k] + xl*p2d[i,k])/(xl+xr)
    end
  end
  for i in 1:im # omit i=0,im+1
    for k in 1:kocn[i]
      za = dz[k-1]; zb = dz[k]
      v[i,k-1] = (zb*p2d[i,k-1] + za*p2d[i,k])/(za+zb)
    end
    k = kocn[i] + 1
    v[i,k-1] = (p2d[i,k-1] + p2d[i,km+1])/2 # bottom
  end

  val_wei(a,w) = isnan(a) ? (0.0,0.0) : (a,w)

  pv = fill(NaN, ivax, kvax)
  for kv in kvax, iv in ivax
    (a,wa) = val_wei(u[iv,kv  ], 1/dz[kv]  )
    (b,wb) = val_wei(u[iv,kv+1], 1/dz[kv+1])
    (l,wl) = val_wei(v[iv,  kv], 1/dx[iv]  )
    (r,wr) = val_wei(v[iv+1,kv], 1/dx[iv+1])
    (wa==0 && wb==0 && wl==0 && wr==0) && continue
    pv[iv,kv] = (a*wa + b*wb + l*wl + r*wr)/(wa+wb+wl+wr)
  end
  xxax = OffsetArray(vcat(0.0,cumsum(dx[1:im])), ivax)
  zzax = OffsetArray(vcat(0.0, - cumsum(dz[1:km])), kvax)
  return xxax, zzax, pv
end
#=
function remap_for_plotting(p2d; dx, dz, xax, zax, iocn)
  (iax,kax) = axes(p2d)
  im = iax[end-1]
  km = kax[end-1]
  @assert iax == 0:(im+1)
  @assert kax == 0:(km+1)
  @assert axes(dx,1) == iax
  @assert axes(dz,1) == kax
  ivax = 0:im
  kvax = 0:km
  xxax = vcat(0.0, cumsum(dx[1:im]) )
  zzax = vcat(0.0, - cumsum(dz[1:km]) )
  pv = fill(NaN, ivax, kvax)
  #-- surface ---
  pv[0, 0] = (p2d[0, 1] + p2d[1, 1])/2 # wall value below
  pv[im,0] = (p2d[im,0] + p2d[im,1])/2 # surface value to the left
  kv = 0
  for iv in 1:im-1 # surface-interior
    il = iv
    ir = iv + 1
    l = (p2d[il,0] + p2d[il,1])/2 # surface value to the left
    r = (p2d[ir,0] + p2d[ir,1])/2 # surface value to the right
    xl = xxax[iv] - xax[il]
    xr = xax[ir] - xxax[iv]
    pv[iv,0] = (l*xr + r*xl)/(xr + xl)
  end
  #-- right edge ---
  pv[im,km] = (p2d[im,km] + p2d[im,km+1])/2 # bottom value to the left
  kv = 0
  for kv in 1:km-1 # right-edge--interior
    ka = kv
    kb = kv + 1
    a = (p2d[im,ka] + p2d[im+1,ka])/2 # x=Lx value above
    b = (p2d[im,kb] + p2d[im+1,kb])/2 # x=Lx value below
    za = zax[ka] - zzax[kv]
    zb = zzax[kv] - zax[kb]
    pv[im,kv] = (a*zb + b*za)/(za + zb)
  end
  #-- slope -> bottom ---
  for k in 1:km # cell center
    i = iocn[k] # cell center
    w = (p2d[0,k] + p2d[i,k])/2
    ??????????
  end
  for kv in kvax, iv in ivax
  end
end
=#


"""
This function maps the irregular grid of ctwmodes-007.f90
  onto the 1/2-spacing regular grid.
Input:
  p2d[0:im+1,0:km+1]: p of a mode already mapped to x-z coordinates.
  dx[0:im+1], dz[0:km+1]
The new gridpoints are at vertices.
      ii=0          ii=2*im
kk=0    +--+-- . . . +
        |
        :
kk=2*km +

This function assumes that p2d[i,k] == NaN where it's undefined;
 it exploits the fact that "v `op` NaN == NaN".

Note:
 At a vertical wall, the wall value is stored to the left; it's not average.
 At a flat bottom, the bottom value is stored below bottom; it's not average.
"""
function remap_diag(p2d; dx, dz, iocn, kocn)
  iax = axes(p2d,1)
  kax = axes(p2d,2)
  im = length(iax) - 2
  km = length(kax) - 2
  @assert iax == 0:im+1
  @assert kax == 0:km+1
  @assert axes(dx,1) == iax
  @assert axes(dz,1) == kax
  iiax = 0:(2*im)
  kkax = 0:(2*km)
  iiext = -1:(2*im+1)
  kkext = -1:(2*km+1)

  xx = OffsetArray{eltype(dx),1}(undef, iiax)
  zz = OffsetArray{eltype(dz),1}(undef, kkax)
  xx[0] = 0
  for i in 1:im
    ii = i*2
    xx[ii-1] = xx[ii-2] + dx[i]/2
    xx[ii]   = xx[ii-2] + dx[i]
  end
  zz[0] = 0
  for k in 1:km
    kk = k*2
    zz[kk-1] = zz[kk-2] - dz[k]/2
    zz[kk]   = zz[kk-2] - dz[k]
  end

  pp = fill(eltype(p2d)(NaN), iiext, kkext)

  # Scan original grid cells and fill values
  #   at the centers of the edges.
  for k in 1:km # scan the cell centers
    kk = 2*k - 1 # cell center in the finer grid
    for i in 1:im # scan the cell centers
      ii = 2*i - 1 # cell center in the finer grid

      if i >= iocn[k] && k <= kocn[i]
        # just copy cell center
        pp[ii,kk] = p2d[i,k]
        # center of right edge of the cell
        pp[ii+1,kk] = (dx[i]*p2d[i+1,k] + dx[i+1]*p2d[i,k])/(dx[i]+dx[i+1])
        # center of top edge of the cell
        pp[ii,kk-1] = (dz[k-1]*p2d[i,k] + dz[k]*p2d[i,k-1])/(dz[k-1]+dz[k])
      end

      #-- wall-slope-bottom and x=Lx boundary
      if iocn[k] == i && kocn[i] != k # we're to the right of vertical wall
        pp[ii-1,kk] = p2d[i-1,k] # copy from left lobe
      end
      if kocn[i] == k && iocn[k] != i # we're at flat bottom
        pp[ii,kk+1] = p2d[i,k+1] # copy from bottom lobe
      end

    end
  end

  # Special case: top edge of cell (i,k)=(im+1,1)
  let k = 1, i = im+1
    kk = 2*k - 1 # cell center in the finer grid
    ii = 2*i - 1 # cell center in the finer grid
    pp[ii,kk-1] = pp[ii-2,kk-1] # copy from the left defined above
  end

  # Special case: bottom edge of cell (i,k)=(im+1,km)
  let k = km, i = im+1
    kk = 2*k - 1 # cell center in the finer grid
    ii = 2*i - 1 # cell center in the finer grid
    pp[ii,kk+1] = pp[ii-2,kk+1] # copy from the left defined above
  end

  # Fill values at the vertices of the original grid cells.

  # Upper-left corner is a special case.
  # If the gridpoint below is defined, pp[0,0] = pp[0,1],
  #   assuming that ∂(pp)/∂z = 0 even with the free-surface condition;
  #   if it's not defined, take the value from the right, pp[0,0] = pp[1,0].
  #
  pp[0,0] = isnan(pp[0,1]) ? pp[1,0] : pp[0,1]

  # Go over the vertices
  for kk in 0:2:(2*km)
    for ii in 0:2:(2*im)
      (ii == 0 && kk == 0) && continue
      dza = dz[kk÷2]
      dzb = dz[kk÷2 + 1]
      dxl = dx[ii÷2]
      dxr = dx[ii÷2 + 1]
      a  = pp[ii,kk-1] # above
      b  = pp[ii,kk+1] # below
      r  = pp[ii+1,kk] # right
      l  = pp[ii-1,kk] # left
      rb = pp[ii+1,kk+1] # right below
      la = pp[ii-1,kk-1] # left above
      il = ii÷2     # cell center to the left
      ir = ii÷2 + 1 # cell center to the right
      ka = kk÷2 # cell center above
      if isnan(a) && isnan(b) && isnan(l) && isnan(r) # below the slope-bottom.
        continue
      elseif kk == 0 # surface
        pp[ii,kk] = (dxl*r + dxr*l)/(dxl + dxr) # horizontal interpolation
      elseif ir == iocn[ka] && ka != kocn[ir] # at the lower edge of a vertical
        if isnan(b) # we're at the upper edge of a diagonal slope
          @assert !isnan(rb)
          ddiag = sqrt(dzb^2 + dxr^2)
          pp[ii,kk] = (dza*rb + ddiag*a)/(dza + ddiag) # along the slope
        else # below is also vertical
          pp[ii,kk] = (dza*b + dzb*a)/(dza + dzb) # vertical interpolation
        end
      elseif il == iocn[ka] && ka == kocn[il] # lower edge of diagonal
        ddiag = sqrt(dza^2 + dxl^2)
        if isnan(b) # we're at the left edge a flat bottom
          #isnan(r) && (@show ii, kk, il, ka)
          @assert !isnan(r)
          #pp[ii,kk] = (dxl*la + ddiag*r)/(dxl + ddiag) # along the slope
          @assert !isnan(a)
          pp[ii,kk] = a # interior value: vertical gradient must be small.
        else # below is a vertical
          pp[ii,kk] = (dzb*la + ddiag*b)/(dzb + ddiag) # along the slope
        end
      elseif il != iocn[ka] && ka == kocn[il] # right edge of flat bot
        pp[ii,kk] = (dxl*r + dxr*l)/(dxl + dxr)
      else # interior -> vertical interpolation
        pp[ii,kk] = (dza*b + dzb*a)/(dza + dzb)
      end
    end
  end

  # take the [iiax, kkax] section of pp[iiext, kkext]
  pp_shrunk = OffsetArray(pp[iiax, kkax], iiax, kkax)

  return (xx=xx, zz=zz, pp=pp_shrunk)
end


# == Inner product & orthogonality ==================

"""
Calculate the wall values for pre ver.007 grid,
  where the wall-boundary values were the average
  between the left lobe (i=0) and the next "interior" (i=iocn[k]) value.
input: F[i,k] ... x-z field
output: Fwall[1:kmax]
"""
function wall_vector(F2d; iocn, kmax)
  [ (F2d[0,k] + F2d[iocn[k],k])/2 for k in 1:kmax]
end


"""
Pick up values on the slope for the "diagonal" case (ctwmodes-007.f90).
input: F[i,k] ... x-z field
  iocn[k]: i index of the leftmost ocean gridpoint at depth k.
  kocn[i]: k index of the bottom-most ocean gridpoint at position i.
output: Fslope[1:kmax]
"""
function slope_vector(F2d; iocn, kocn, kmax)
#  @show axes(F2d)
  Fslope = Vector{eltype(F2d)}(undef, kmax)
  for k in 1:kmax
    i = iocn[k]
    # @show i, k
    if kocn[i] == k # diagonal boundary cell
      Fslope[k] = F2d[i,k]
    else
      Fslope[k] = F2d[i-1,k] # left lobe because vertical wall
    end
  end
  #@show Fslope
  return Fslope
end

"""Counterpart of slope_vector() to pick up "bottom" values."""
function bottom_vector(F2d; iocn, kocn, imax)
  Fbot = Vector{eltype(F2d)}(undef, imax)
  for i in 1:imax
    k = kocn[i]
    if iocn[k] == i # diagonal boundary cell
      Fbot[i] = F2d[i,k]
    else
      Fbot[i] = F2d[i,k+1] # below bottom because flat bottom
    end
  end
  return Fbot
end


"""
Boundary inner product.
u[i,k], v[i,k]: x-z field
iocn[k]: i index of the leftmost ocean gridpoint at depth k.
"""
#function binner(u,v; iocn, dz, kmax)
#  um = wall_vector(u; iocn, kmax)
#  vm = wall_vector(v; iocn, kmax)
function binner(u,v; get_slope_vec, dz, kmax)
  um = get_slope_vec(u)
  vm = get_slope_vec(v)
  s = 0.0
  for k in 1:kmax
    #um = (u[0,k] + u[iocn[k],k])/2 # wall value
    #vm = (v[0,k] + v[iocn[k],k])/2 # wall value
    s += dz[k] * um[k] * vm[k]
  end
  return s
end

"""
Calculate or show the correlation matrix between the modes.
Input: pij[i,j,n]
If info == false, calculate the matrix and return it.
If info == true, print out the matrix components on the screen.
"""
function orthogonality(pij; iocn, kocn, dz, eps = 1.0e-14, info = true)
  nax = axes(pij,3)
  zax = axes(pij,2)
  kmax = length(zax) - 2
  @assert zax == 0:(kmax+1)
  #get_slope_vec(u) = wall_vector(u; iocn, kmax)
  get_slope_vec(u) = slope_vector(u; iocn, kocn, kmax)
  inner(l,n) = binner(view(pij,:,:,l), view(pij,:,:,n)
                      ; get_slope_vec, dz, kmax)
  sz = [sqrt(inner(n,n)) for n in nax] # norms of the vectors
  corr = info ? nothing : OffsetArray{Float64}(undef,nax,nax)
  for l in nax
    info && print("$(l): ")
    for n in l:nax[end] # upper triangle
      s = inner(l,n) / (sz[l] * sz[n])
      if ! info
        corr[l,n] = s
        continue
      end
      if (n == l)
        print("1,")
      elseif (abs(s) > eps)
        print(f"\%8.1e(s),")
      else
        print(f"-,")
      end
    end # for n
    info && println()
  end # for l
  # fill the lower triangle
  if ! info
    for l in nax
      for n in nax[begin]:(l-1)
        corr[l,n] = corr[n,l]
      end
    end
  end
  return corr
end

"""
Returns normalization factor such that F(x,z)/factor be normalized.
This function makes binner(u,u)/totdep == 1.
"""
function normfactor(p2d; iocn, kocn, dz, kmax)
  totdep = sum(dz[1:kmax])
  get_slope_vec(u) = slope_vector(u; iocn, kocn, kmax)
  a2 = binner(p2d, p2d; get_slope_vec, dz, kmax) / totdep
  a = sqrt(a2)
  sign_ul = sign(p2d[1,1]) # upper left corner
  return a * sign_ul
end

"""
Normalize p2d[0:im+1, 0:km+1] so that p2d[1,1] = 1.
(Another common choice is that binner(u,u) == 1
 but we don't need it.)
"""
function normalize!(p2d; mes = nothing)
  a = p2d[1,1]
  # (a == 0) && error("p2d[1,1] == 0")
  if a == 0
    (! isnothing(mes)) && print(mes)
    println(stderr, "p2d[1,1] == 0")
  else
    @. p2d = p2d / a
  end
end


# == Verify the solutions ===
"""
Measure the accuracy of p2d[:,:].
  p2d[:,:] == prij[:,:,n]
  ginv = 1/g; ginv = 0 means rigid-lid
  Assume the given c > 0, whereas the actual c < 0.
If the |values| < eps, accuracy is not assessed.
Otherwise, use relative error.
"""
function accuracy(p2d; dx, dz, bvf2e, c, f0, ginv, iocn, kocn, eps = 1.0e-14)
  function err(a,b,eps)
    m = max(abs(a),abs(b))
    (m < eps) ? m : abs(a - b)/m
  end
  fsq = f0*f0
  foN2 = fsq ./ bvf2e # (f^2/N^2)
  im = size(dx,1) - 2
  km = size(dz,1) - 2
  # interior
  for k in 1:km
    i1 = iocn[k]
    for i in i1:im
      k1 = kocn[i]
      il = (i==i1) ? 0 : (i-1)
      ir = i+1
      ka = k-1
      kb = (k==k1) ? (km+1) : (k+1)
      #m = maximum(map(abs, (p2d[i,k], p2d[il,k]
      #                      ,p2d[ir,k], p2d[i,ka], p2d[i,kb])))
      #(m < eps) && continue
      pxl = (p2d[i, k] - p2d[il,k])/((dx[i]   + dx[i-1])*0.5)
      pxr = (p2d[ir,k] - p2d[i, k])/((dx[i+1] + dx[i]  )*0.5)
      pxx = (pxr - pxl)/dx[i]
      mx = maximum(map(abs, (p2d[i,k], p2d[il,k], p2d[ir,k])))
      mx /= (dx[i]*dx[i])
      pza = (p2d[i,ka] - p2d[i,k ])/((dz[k-1] + dz[k  ])*0.5)*foN2[k-1]
      pzb = (p2d[i,k ] - p2d[i,kb])/((dz[k  ] + dz[k+1])*0.5)*foN2[k]
      pzz = (pza - pzb)/dz[k]
      mz = maximum(map(abs, (p2d[i,k], p2d[i,ka], p2d[i,kb])))
      mz = mz / (dz[k]*dz[k]) * foN2[k-1]
      s = pxx + pzz
      e = abs(s/max(mx,mz)) / max(abs(pxx/mx), abs(pzz/mz))
      (e > eps) && println("i,k,err,pxx/mx,pzz/mz = $(i), $(k), $(e), $(pxx/mx), $(pzz/mz)")
    end
  end
  # surface
  for i in 1:im
    lhs = (p2d[i,0] - p2d[i,1])/((dz[0] + dz[1])*0.5)*foN2[0]
    rhs = - (p2d[i,0] + p2d[i,1])/2*(fsq*ginv)
    e = err(lhs, rhs, eps)
    # err = abs(lhs-rhs) / max(abs(lhs),abs(rhs))
    (e > eps) &&
      println("surface: i,err,lhs,rhs = $(i), $(e), $(lhs), $(rhs)")
  end
  # right edge
  for k in 1:km
    a = p2d[im+1,k]
    b = p2d[im,k]
    e = err(a,b,eps)
    (e > eps) && println("right: i,err = $(k), $(e)")
  end
  # bottom
  for i in 1:im
    k1 = kocn[i]
    a = p2d[i,k1]
    b = p2d[i,km+1]
    e = err(a,b,eps)
    (e > eps) && println("bottom: i,err = $(i), $(e)")
  end
  # slope-wall
  for k in 1:km
    i1 = iocn[k]
    lhs = (p2d[i1,k] - p2d[0,k])/((dx[i1-1] + dx[i1])*0.5)
    rhs = ((p2d[i1,k] + p2d[0,k])/2) * (-f0/c) # flip the sign of c
    e = err(lhs, rhs, eps)
    (e > eps) &&
      println("surface: k,e,lhs,rhs = $(k), $(e), $(lhs), $(rhs)")
  end
end


# == grid data and netCDF files ==================

"""
Calculate cell-center points, xax[0:im+1] and zax[0:km+1]
 from cell-widths dx[0:im+1] and dz[0:km+1]
"""
function calc_grid(; dx, dz)
  iax = axes(dx,1)
  kax = axes(dz,1)
  im = length(iax) - 2
  km = length(kax) - 2
  @assert iax == 0:(im+1)
  @assert kax == 0:(km+1)
  xax = OffsetArray(Vector{Float64}(undef, length(iax)), iax)
  zax = OffsetArray(Vector{Float64}(undef, length(kax)), kax)
  xax[0] = -0.5 * dx[0]
  for i = 1:im+1
    xax[i] = xax[i-1] + (dx[i-1] + dx[i])/2
  end
  zax[0] = 0.5 * dz[0]
  for k = 1:km+1
    zax[k] = zax[k-1] - (dz[k-1] + dz[k])/2
  end
  # @show xax
  # @show zax
  @show xax[1:3] xax[end-2:end]
  @show zax[1:3] zax[end-2:end]

  return (xax, zax)
end

"calc lower edges zzax[0:km] of grid cells."
function calc_edges(; dz)
  km = size(dz,1) - 2
  zzax = OffsetArray(Vector{Float64}(undef, km+1), 0:km)
  zzax[0] = 0.0
  for k in 1:km
    zzax[k] = zzax[k-1] - dz[k]
  end
  zzax
end


"""
Save the eigenvalues and eigenvectors as a netCDF file.
#pr[ik,n]: the original array of eigenvectors.
prij[i,k,nn]: array of eigenvectors in x-z grid
#num[i,k]: map from (i,k) to the ik index of pr.
#nidx[nn]: extracts modes in the mode-number order.
cee[nn]: c's in the mode-number order.
xax[i], zax[k]: cell-center coordinates of prij[i,k,:]
dx[i], dz[k]: cell widths and heights
f0: Coriolis parameter (scalar)
bvf2e[kk]: N^2(z) defined at edges xxax[kk]
normfactor[nn]: normalization factor such that
  integral (p/normfactor)^2 dz/totdep == 1.
"""
function save_modes(fnam; prij, cee, xax, zax
                    ,pslope = nothing, pbot = nothing
                    ,dx = nothing, dz = nothing, f0 = nothing, bvf2e = nothing
                    ,kocn = nothing, iocn = nothing
                    ,outtype = Float32
                    ,title = nothing
                    ,normfactor = nothing
                    )
  nooffs = OffsetArrays.no_offset_view
  fillvalue = outtype(NaN)
  #=
  prij = map_to_xz(; pr=pr, nidx=nidx, num=num, outtype=outtype)
  #foreach(normalize!, eachslice(prij,dims=3))
  for m in axes(prij,3)
    #println("m=$(m): normalizing")
    normalize!(view(prij,:,:,m); mes = "mode $(m): ")
    #>>>>
    #@show prij[end,end,m]
    #<<<<<
  end
  =#
  @assert size(prij,1) == size(xax,1)
  @assert size(prij,2) == size(zax,1)
  @assert axes(prij,3) == axes(cee,1)
  mode_numbers = collect(axes(prij,3))
  NCDataset(fnam, "c") do ds
    (! isnothing(title)) && (ds.attrib["title"] = title)
    defVar(ds, "xax", nooffs(xax), ("xax",),
           attrib=(axis="X", units="meters",
                   long_name = "x cell center") )
    defVar(ds, "zax", nooffs(zax), ("zax",),
           attrib=(axis="Z", units="meters",
                   long_name = "z cell center") )
    defVar(ds, "mode", Int32.(mode_numbers), ("mode",),
           attrib=(axis="E", long_name="mode number"))
    if ! isnothing(bvf2e)
      zzax = calc_edges(dz = dz)
      defVar(ds, "zzax", nooffs(zzax), ("zzax",),
             attrib=(axis="Z", units="meters",
                     long_name="z lower edges") )
      defVar(ds, "bvf2e", nooffs(bvf2e), ("zzax",)
             ,attrib=(units="1/s^2"
                      ,long_name="N^2(z)") )
    end
    (! isnothing(dx)) &&
      defVar(ds, "dx", nooffs(dx), ("xax",),
             attrib=(units="meters",long_name="x gridcell width") )
    (! isnothing(dz)) &&
      defVar(ds, "dz", nooffs(dz), ("zax",),
             attrib=(units="meters", long_name="z gridcell width") )
    (! isnothing(f0)) &&
      defVar(ds, "f0", f0, (), # scalar
             attrib=(units="1/s", long_name="Coriolis parameter") )
    if ! isnothing(kocn)
      #u = typemin(eltype(kocn))
      #kocn1 = vcat(u, kocn, u)
      #@assert length(kocn1) == length(dx)
      @assert length(kocn) == length(dx)
      u = kocn[end]
      @assert u < 0
      defVar(ds, "kocn", nooffs(kocn), ("xax",),
             attrib = (_FillValue = u, long_name="max k of ocean interior") )
    end
    if ! isnothing(iocn)
      #u = typemin(eltype(iocn))
      #iocn1 = vcat(u, iocn, u)
      #@assert length(iocn1) == length(dz)
      @assert length(iocn) == length(dz)
      u = iocn[end]
      @assert u < 0
      defVar(ds, "iocn", nooffs(iocn), ("zax",),
             attrib = (_FillValue = u, long_name="min i of ocean interior") )
    end
    (! isnothing(normfactor)) &&
      defVar(ds, "normfactor", nooffs(normfactor), ("mode",),
             attrib=(long_name="Normalization factor"
                     ,comment="integral (p/normfactor)^2 dz/totdep == 1") )
    F = defVar(ds, "pmode", outtype, ("xax", "zax", "mode"),
               attrib=Dict("_FillValue"=>fillvalue,
                       "long_name" => "CTW p mode (x,z)"))
    c = defVar(ds, "c", Float64, ("mode",),
               attrib=Dict("long_name"=>"gravity-wave speed",
                       "units"=>"m/s"))
    for (i, n) in pairs(mode_numbers) # (index, value)
      F[:,:,i] = prij[:,:,n]
      c[i] = cee[n]
    end

    (!isnothing(pslope)) &&
      defVar(ds, "pslope", nooffs(pslope), ("zax","mode"),
             attrib=Dict("_FillValue" => fillvalue
                     , "long_name"=>"slope boundary values of mode") )
    (!isnothing(pbot)) &&
      defVar(ds, "pbot", nooffs(pbot), ("xax","mode"),
             attrib=Dict("_FillValue" => fillvalue
                     , "long_name"=>"'bottom' boundary values of mode") )
  end
  println("save_modes: $(fnam)")
end

"""
A function to test, and demonstrate the use of, this module,
  not intended for outside use.
"""
function test_save_modes()
  (; num, kocn, iocn, alphaR, alphaI, pr, dx, dz, f0, bvf2e) = read_fort_data() # ; intype=Float64
  (cee_west, idx_west, cee_east, idx_east) = sort_modes(alphaR, alphaI)
  (xax, zax) = calc_grid(; dx=dx, dz=dz)
  save_modes("tmp-west.nc"; prij=pr, cee=cee_west
             ,xax=xax, zax=zax, dx=dx, dz=dz, f0=f0, bvf2e=bvf2e
             ,kocn=kocn, iocn=iocn)
  save_modes("tmp-east.nc"; prij=pr, cee=cee_east,
             xax=xax, zax=zax, dx=dx, dz=dz, title="test_save_modes()")
end


"""
Open the netCDF file and wrap the variables in OffsetArrays.
(The dataset is closed because the data are copied into
 arrays. We could return ds in the NamedTuple to save it.)
"""
function read_modes(ncfnam)
  gravit_default = 9.8 # m/s^2

  ds = NCDataset(ncfnam, "r")

  dx0 = ds["dx"]
  dz0 = ds["dz"]
  imax = size(dx0,1) - 2
  kmax = size(dz0,1) - 2
  iax = 0:(imax+1)
  kax = 0:(kmax+1)
  kkax = 0:kmax
  nax = let mode0 = ds["mode"]
    nsize = size(mode0,1)
    nax = 0:(nsize-1)
    @assert nax == mode0[:]
    nax
  end

  dx = OffsetArray(dx0[:], iax)
  dz = OffsetArray(dz0[:], kax)
  xax = OffsetArray(ds["xax"][:], iax)
  zax = OffsetArray(ds["zax"][:], kax)
  kocn = OffsetArray(ds["kocn"][:], iax)
  iocn = OffsetArray(ds["iocn"][:], kax)
  pij = OffsetArray(ds["pmode"][:,:,:], iax, kax, nax)
  #@show iocn
  #@show pij[1,1,1]
  cee = OffsetArray(ds["c"][:], nax)
  bvf2e = OffsetArray(ds["bvf2e"][:], kkax)
#  return (ds=ds

  f0 = ds["f0"][]
  gravit = if haskey(ds,"gravit")
    ds["gravit"][:]
  else
    println("$(ncfnam) doesn't contain gravit. Set gravit=$(gravit_default)")
    gravit_default
  end

  return (pij=pij
          ,cee=cee
          ,bvf2e=bvf2e
          ,kocn=kocn
          ,iocn=iocn
          ,dx=dx
          ,dz=dz
          ,xax=xax
          ,zax=zax
          ,nax=nax
          ,f0=f0
          ,gravit=gravit
          )
end


end
