#!/usr/bin/env julia
# =============================================================================
# Tool Name:    CTW configure script
# Description:  Script to produce configuration files
#                 (Grid.txt, Ne.txt, and ctwmodes_pars.f90).
# Author:       Ryo Furue <ryofurue@gmail.com>
# License:      MIT License
# Version:      v"0.1"
# Created:      2025-??-??
# Updated:      2025-07-01
# =============================================================================
"""
"""

const SCRIPT_VERSION = v"0.1" # `Base.VERSION` is used by the interpreter.

# Optional: Print version info at runtime
println("CTW configure script, version ", SCRIPT_VERSION)


# conf.jl v0.1:
#   Adding functions that determine dx, dz, etc. from h(x).
#   Adding functions that determine Ne(z) from user-provided gridded Ne.
#
# Instead of adding sophisticated input and output functionalities
#   to ctwmodes.f90, we provide scripts to prepare the input files
#   (this script) and to process the output (fortbin ??????????????????)
#
# This script generates
#
# - ctwmodes_pars.f90 . . . configuration of the Fortran program
#     . . . for grid sizes im, km, to be included by ctwmodes.f90.
#
# - Grid.txt . . . topography and geometry
#     . . . determines dx[i], dz[k], and kocn[i], to be read by ctwmodes.f90.
#
# - Ne.txt . . . Brunt-Väisälä frequency
#     . . . determines Ne[k], to be read by ctwmodes.f90.
#
# Topography and geometry can be determined either
#   1) by a function defined in here, or
#   2) from the text file "Conf_topo.txt".
#
# Ne can be determined
#   1) by a simple function Nfunc(z), or
#   2) from the text file "Conf_Ne.txt".
#
# This script isn't necessary: The user can generate these files
#    for themself. This script just serves as a summary of the procedure.
#
# This script is a rework of gen-h-N-???.jl (gen-h-N-001.jl is the latest),
#   which has over time become too complex.
#

# Experimental Julia language extension. For the @np macro.
using NamedPositionals

# To interpolate h(x) and N(z) onto gridpoints.
using Interpolations

# To use the consistent indexing convention with the CTW code:
#  xx[0]--dx[1]--xx[1]--dx[2]-- . . . --dx[nn-1]--xx[nn-1]
# we need zero-based arrays.
using OffsetArrays

# Internal use, for testing
using Random  # We generate a ramdom slope.
using GLMakie # Plot the slope for testing.

#=
"""
For internal use. Split the text line if it's not a comment line.
"""
function split_if_not_comment(str)
  isnothing(match(r"\s*#", str)) ? split(str) : nothing
end
=#


"""
For internal use.
Prints "prefix: message"
prefix is nameof(var"#self#") when inside a function, otherwise basename(@__FILE__)
"""
macro say(msg)
  return esc(
    :(let __msg__ = $msg
        prefix = try
          nameof(var"#self#") # works only inside real function bodies
        catch
          basename(@__FILE__) # fallback: filename where macro is used
        end
        println(prefix, ": ", __msg__)
      end))
end


#=
"""
Like `@show` but returns the string instead of printing it out.
"""
macro sshow(ex)
    # turn the expression into a string, same as @show does
    ex_str = sprint(show, ex)
    return esc(quote
        string($ex_str, " = ", $ex)
    end)
end
=#

"""
For internal use.
Like `@show` but returns the string instead of printing it out.
"""
macro sshow(args...)
  parts = []
  for ex in args
    tmp = gensym(:_sshow_tmp)
    # Use sprint with :source, which matches `@show` style better
    ex_str = sprint(io -> Base.show_unquoted(io, ex))
    push!(parts, :(let $tmp = $ex
                     string($(QuoteNode(ex_str)), " = ", $tmp)
                   end))
  end
  return esc(:(join([$(parts...)], ", ")))
end


"""
For internal use.
Read a line and parse it.
"""
function read_line_parse(is, typelist)
  ll = readline(is)
  strlist = split(ll)
  if length(strlist) == 0
    error("empty line")
  elseif length(strlist) > length(typelist)
    @say("Types are $(typelist) but there are more items on the line: $(ll)")
  elseif length(strlist) < length(typelist)
    error("Types are $(typelist) but there are fewer items on the line: $(ll)")
  end
  res = map(x -> parse(x[1], x[2]), zip(typelist, strlist))
  if length(res) == 0
    error("parse failed: $(ll)")
  elseif length(res) == 1
    return res[1]
  else
    return Tuple(x for x in res) # Vector -> Tuple
  end
end


"""
For internal use.
Read the contents of the next line if it exists and is non-empty;
 returns nothing if there is none or it's empty.
"""
function read_line_if_exists(is)
  eof(is) && return nothing
  s = strip(readline(is)) # remove white space characters
  return (s != "") ? s : nothing
end


#--- helper functions
"""
Helper function to interpolate the "function" (xs0[:],ys0[:])
   on to the xs[i] points.
Outside the xs0[:] range, extrapolate assuming constant.
  If x is to the left  of the xs0 range, y = ys0[begin];
  If x is to the right of the xs0 range, y = ys0[end].

The interpolator allows only for strictly increasing
 x coordinates. If it's decreasing, reverse all vectors.
"""
function interpolate_helper(;xs0,ys0,xs)
  inc = is_increasing(xs0; strict=true, errorstrength=0)
  if ! inc
    xs0 = reverse(xs0)
    ys0 = reverse(ys0)
    xs = reverse(xs)
  end

  ys = similar(xs)

  # This special-casing is necessary only because
  #   linear_interpolation() doesn't support single-knot cases.
  if length(ys0) == 1
    ys  .= ys0[begin] # same value everywhere
  else
    itp = linear_interpolation(xs0,ys0; extrapolation_bc = NaN)
    for i in axes(xs,1)
      ys[i] = itp(xs[i])
    end
    let i = findfirst(!isnan, ys)
      ys[begin:(i-1)] .= ys0[begin]
    end
    let i = findlast(!isnan, ys)
      ys[(i+1):end]   .= ys0[end]
    end
  end

  return inc ? ys : reverse(ys)
end

"""
Helper funciton: Are the vector elements in the increasing order?
errorstrength == 0: no message.
errorstrength == 1: print message.
errorstrength == 2: stop on error.
"""
function is_increasing(v; strict = true, errorstrength = 2)
  length(v) == 0 && error("zero length vector")
  length(v) == 1 && error("one length vector")
  lte = strict ? ((x,y) -> x < y) : ((x,y) -> x <= y)
  grt = ∘(!,lte) # not smaller
  a = v[begin]
  isgood = true
  for i in axes(v,1)[begin+1:end]
    if grt(a,v[i])
      (errorstrength > 0) &&
        println(stderr, "vec must increase monotonically:" *
        " v[$(i-1)],v[$(i)]=$(v[i-1]),$(v[i])" )
      (errorstrength > 1) && error("")
      isgood = false
    end
    a = v[i]
  end
  return isgood
end


#--- Topography -----------
# - Some of the functions use Nfunc(z) to determine dz.
#    Nfunc(z) is the Brunt-Väisälä freuency as function of depth, N(z),
#    where z = 0 is the sea surface and z < 0.
# - Some of the functions read h(x) and dx from a file and determine dz.

"""
Read h(x) from a file.
Output: dx, dz, kocn
  (extra outputs xx and hh are provided just for convenience in testing.)

Input file format 1:
```
  nn
  xx[0]    h(xx[0])
  xx[1]    h(xx[1])
  . . .
  xx[nn-1] h(xx[nn-1])
  xmax
```
Input file format 2:
```
  nn
  xx[0]    h(xx[0])
  xx[1]    h(xx[1])
  . . .
  xx[nn-1] h(xx[nn-1])
  xmax
  im
  dx[1]
  dx[2]
  . . .
  dx[im]
```

For both formats, the given x axis is shifted (including xmax)
  so that xx[0] = 0.

For format 1, `dx[i] = xx[i] - xx[i-1]`.
The region beyond `xx[nn]` is the flat-bottom region, where `dx` is gradually stretched until the gridpoint runs beyond `xamx`.

For format 2, h(x) is determined by linear interpolation at
  the gridpoints
```
     xx'[1] = 0
     xx'[2] = xx'[1] + dx[1]
     xx'[3] = xx'[2] + dx[2]
    ...
```
or where the gridpoint xx'[i] is beyond the slope (xx'[i] > xx[nn]),
`h(xx'[i]) = h(xx[nn])`, which is the flat-bottom region.

dz_samp is used only when hh[0] > 0 (depth is finite at the left edge).
"""
function slope_from_file(fnam; dz_samp = 20) # sample dz in meters
  (;xx,hh,dx) = get_hh_from_file(fnam)
  im = length(xx) - 1
  @assert axes(xx,1) == 0:im
  @assert axes(hh,1) == 0:im
  @assert axes(dx,1) == 1:im
  @say @sshow hh
  (hh[0] < 0) && error("hh[0]=$(hh[0]) < 0.")
  dz = # initial step(s)
    if hh[0] > 0
      n = round(Int, hh[0]/dz_samp, RoundUp) # ceiling
      dz0 = hh[0] / n
      fill(dz0, n) # repeat same value
    else
      (eltype(hh))[] # empty vector
    end
  (length(dz) != 0) && @say ("Initial step(s): " * @sshow dz)
  @say "Generating dz and kocn from hh ..."
  k = length(dz)
  kocn = fill(typemin(Int), im)
  for i in axes(dx,1)
    if hh[i] == hh[i-1] # flat
      kocn[i] = k
    else # diagonal
      k += 1
      kocn[i] = k
      push!(dz, hh[i] - hh[i-1])
    end
  end
  #-- check internal consistency --
  cnt = 0
  for i in axes(kocn,1)
    if kocn[i] <= 0
      println(stderr,"i,kocn[i] = $(i), $(kocn[i])")
      cnt += 1
    end
  end
  (cnt > 0) && error("kocn includes non positive integer(s).")
  @say @sshow dx
  @say @sshow kocn
  @say @sshow dz
  (xx[begin] != 0) &&
    @say "Initial position x = $(xx[begin]) will be shifted to x = 0"
  return (dx = dx, dz = dz, kocn = kocn) #, hh=hh, xx=xx)
end

"""
Helper function to sope_from_file(): get and calculate (xx,hh,dx).
Regardless of the format of the original data file,
in the output, the wall-slope-bottom is defined to be
  [(xx[i], hh[i]) for i in 0:im ]
The x axis is shifted so that xx[i] = 0.
The cell widths
  dx[i] = xx[i] - xx[i-1] for i = 1:im
is returned just for convenience although it's redundant.
"""
function get_hh_from_file(fnam)
  @say "Reading h(x) from $(fnam)"
  is = open(fnam; write=false)

  nn = @np read_line_parse(is, typelist = (Int,);)

  xx0 = OffsetArray{Float64,1}(undef, 0:nn-1)
  hh0 = OffsetArray{Float64,1}(undef, 0:nn-1)
  #@show typeof(nn), nn
  xxax = axes(xx0,1)

  for i in xxax
    (x, h) = read_line_parse(is, (Float64,Float64))
    xx0[i] = x
    hh0[i] = h
  end

  is_increasing(xx0; strict=true)
  is_increasing(hh0; strict=false)
  (hh0[0] < 0) && (xx0[0] >= 0) &&
    error("hh0[0] = $(hh0[0]) < 0, xx0[0] = $(xx0[0]) >= 0")

  xmax = read_line_parse(is, (Float64,))

  if (xx0[begin] != 0)
    xleft = xx0[begin]
    @say "left edge xx[0] == $(xleft) is shifted to x = 0 ."
    xx0 = xx0 .- xleft
    xmax = xmax - xleft
  end

  line = read_line_if_exists(is)
  ret =
    if isnothing(line) # Format 1: no dx provided
      @say "Format 1 detected: Generating dx from xx[:] ..."
      #calc dx from cell edges
      dx = [(xx0[i] - xx0[i-1]) for i in xxax[begin+1:end]]
      #@show typeof(dx)
      (;hh,xx,dx) = extend_domain(;hh=hh0,xx=xx0,dx,xmax)
    else # Format 2
      @say "Format 2 detected: Reading dx ..."
      im = parse(Int,line)
      @say @sshow im
      dx = [read_line_parse(is, (Float64,)) for i in 1:im]
      xx = OffsetArray(vcat(0.0, cumsum(dx)), -1)
      hh = interpolate_helper(;xs0=xx0, ys0=hh0, xs=xx)
      (;hh=hh,xx=xx,dx=dx)
    end
  close(is)
  return ret
end


"""
Helper function to extend dx in the flat-bottom region.
"""
function extend_domain(;hh,xx,dx,xmax)
  (xx[end] >= xmax) && error("Already beyond xmax: xx=$(xx[end]),xmax=$(xmax)")
  hhn = copy(hh)
  xxn = copy(xx)
  dxn = copy(dx)
  while xxn[end] < xmax
    d1 = 1.1 * dxn[end]
    x1 = xxn[end] + d1
    push!(dxn,d1)
    push!(xxn,x1)
    push!(hhn,hh[end]) # flat bottom
  end
  (hh = hhn, xx = xxn, dx = dxn)
end

#=
"""
Helper function to interpolate (hh0[:],xx0[:]) on to xx[i] points.
Outside dht xx0 range, extrapolate assuming constant.
  If xx is left  of xx0[:], hh = hh0[begin] .
  If xx is right of xx0[:], hh = hh0[end]
"""
function interpolate_hh(;hh0,xx0,xx)
  hhint = linear_interpolation(xx0,hh0; extrapolation_bc = NaN)
  hh = similar(xx)
  for i in axes(xx,1)
    hh[i] = hhint(xx[i])
  end
  i = findfirst(!isnan, hh)
  hh[begin:(i-1)] .= hh0[begin]
  i = findlast(!isnan, hh)
  hh[(i+1):end] .= hh0[end]
  return hh
end
=#


# =============
"""
Helper:
Determine dz between z1 and z2 (z1 > z2) such that dz ∝ 1/Nfunc(zz)
  The resultant dz[:]'s exactly split (z1, z2):
    z2 + sum(dz) == z1
"""
function dz_from_N(; Nfunc, z1 = 0.0, z2, dz_sample = 20.0)
  N_samp = Nfunc(z1)
  dz = Vector{Float64}(undef, 0)

  # initial zz[:] so that dz ∝ 1/N(zz)
  zz = z1
  while zz >= z2
    dz0 = dz_sample * N_samp / Nfunc(zz)
    push!(dz, dz0)
    zz = zz - dz0
  end

  # shrink dz so that sum(dz) == z1 - z2
  shrink = (z1 - z2) / sum(dz)
  # @show shrink
  dz .*= shrink
  return dz
end


"""
Helper:
Determine dx and dx for a linear slope spanning from z1 to z2.
We first determine dz[:] from Nfunc(z) and then
  simply use this relation h_x = dz / dx for all dz[k].
Slope part only.
"""
function dz_dx_linslope(; z1, z2, dz_sample, h_x, Nfunc)
  dz = dz_from_N(; Nfunc, z1, z2, dz_sample = dz_sample)
  dx = dz ./ h_x #
  (dx=dx, dz=dz)
end


"""
Example:
Determine a grid for a simple h(x):
  - a little vertical wall at x = 0,
  - followed by a linear slope,
  - followed by a flat bottom.
The vertical wall is necessary to compare our results
  with those from the Brink-Champman code.

Output: dx, dz, kocn
"""
function grid_linslope(; z1, z2, dz_sample, h_x, Nfunc, xmax)
  @assert z1 <= 0.0
  (; dx, dz) = dz_dx_linslope(; z1, z2, dz_sample = dz_sample, h_x, Nfunc)
  ksz = length(dz)
  kocn = collect(1:ksz) # the slope starts at k = 1
  if z1 < 0
    dz_init = -z1 # little vertical wall
    dz = vcat(dz_init, dz) # prepent dz_init
    ksz += 1
    kocn = collect(2:ksz) # the slope starts at k = 2
    println("vertical wall at x = 0; dz_init = $(dz_init)")
  end
  xx = sum(dx) # bottom of the slope
  d = dx[end]
  while xx < xmax
    d = 1.1*d
    push!(dx, d)
    xx += d
    push!(kocn, ksz)
  end
  (dx=dx, dz=dz, kocn=kocn)
end


"""
Example:
"""
function grid_linslope_example()
  x1   = 100.0e3
  xmax = 300.0e3
  hbot = 5000.0
  h_x = hbot/x1
  nf(z) = 0.003
  grid_linslope(; z1=-60.0, z2=-5000.0
                ,dz_sample=30.0, h_x=h_x, Nfunc=nf, xmax=xmax)
end


"""
Testing:
Create an artificial topography to test the Fortran code.
"""
function slope_for_testing()
  Random.seed!(1234) # just for repreducibility

  xmax = 300.0e3
  h_x = totdep / 150.0e3 # a uniform slope
  dx = Float64[]
  dz = Float64[]
  kocn = Int[]
  push!(dz,40.0) # down
  zz = -dz[1]
  k = 1
  xx = 0.0
  while zz > -totdep
    if rand() < 0.5
      d = dz[end] * 1.1
      push!(dz, d)
      zz -= d
      k += 1
    else
      push!(kocn, k)
      d = dz[end] / h_x
      push!(dx, d)
      xx += d
    end
  end
  (xx > xmax) && error("there is no flat bottom region.")
  kmax = length(dz)

  while xx < xmax
    push!(kocn, kmax)
    d = 1.1*dx[end]
    push!(dx, d)
    xx += d
  end

  (dx=dx, dz=dz, kocn=kocn)
end


"""
Produces the Grid file that ctwmode.f90 reads.

Same as the corresponding funciton in gen-h-N-001.jl
  *except* that this one doesn't produce lobe regions.
"""
function write_grid(fnam; dx, dz, kocn)
  im = size(dx,1)
  km = size(dz,1)
  @say "Writing $(fnam) from dx, dz, kocn ."
  @assert size(kocn) == size(dx)
  open(fnam, "w") do os
    println(os, "$(im)  $(km)")
    for i in axes(dx,1)
      println(os, "$(i)  $(dx[i])  $(kocn[i])")
    end
    for k in axes(dz,1)
      println(os, "$(k)  $(dz[k])")
    end
  end
end


#--- Brunt-Väisälä frequency -----------

"""
Helper function: calculate the edge coordinates zzax[0:km] from dz[1:km],
  assuming that zzax[0] = 0 is the surface.
"""
function calc_zzax(dz)
  km = length(dz)
  @assert axes(dz,1) == 1:km
  zzax = OffsetArray{Float64,1}(undef, 0:km)
  zzax[0] = 0.0
  for (k,d) in pairs(dz) # (index, value)
    zzax[k] = zzax[k-1] - d
  end
  zzax
end

"""
Calculate the gridded values of Ne[0:km] at zzax[0:km] from Nfunc(z).
"""
function grid_Nfunc(Nfunc; zzax)
  @say "Generating Ne[:] from Nfunc() and the z axis zzax[:]."
  Nfunc.(zzax)
end

"""
Helper function: Regrid the gridded Ne profile (zs0, Ne0) onto
 the new gridpoints zzax[:].
See `interpolate_helper()` for what happens if zzax points
 fall outside the range of zs0[:].
"""
function regrid_Ne(; zs0, Ne0, zzax)
  interpolate_helper(;xs0=zs0,ys0=Ne0,xs=zzax)
end


"""
Read a N(z) profile from the file and regrid it onto zzax[:].
The file format is
  nn
  dep[1] Ne(-dep[1])
  dep[2] Ne(-dep[2])
  . . .
  dep[nn] Ne(-dep[nn])
Note that depth >= 0 and z = -depth <= 0.
"""
function get_Ne_from_file(fnam; zzax)
  println("Reading Ne(z) from $(fnam) .")
  (zz, bvf) = open(fnam; write=false) do is
    nn = read_line_parse(is, (Int,))
    dep = fill(NaN,nn)
    bvf = fill(NaN,nn)
    for i in 1:nn
      (dep[i],bvf[i]) = read_line_parse(is, (Float64,Float64))
    end
    zz = -dep
    (zz, bvf)
  end
  println("Mapping Ne(z) onto the model grid.")
  regrid_Ne(; zs0=zz, Ne0=bvf, zzax=zzax)
end


"""
Write Ne[:] values defined at zzax[:] points.
The format is
  num_of_points
  z1  Ne1
  z2  Ne2
  . . .
"""
function write_Ne(fnam; Ne, zzax)
  @assert axes(Ne,1) == axes(zzax,1)
  @say "Writing $(fnam) from Ne, zzax ."
  open(fnam; truncate=true) do os
    println(os, size(zzax,1)) # write num of points
    for k in axes(zzax,1)
      println(os, zzax[k], " ", Ne[k])
    end
  end
end

#----------------------------------

"""
Produces the Fortran parameter file.
   outtype = :single | :double
"""
function write_fort_pars(fnam; im, km, f0
                         ,grid_file, Ne_file
                         ,freesurface=true, outtype = :single
                         ,verbose=true)
  @say "Writing $(fnam) ..."
  fs = freesurface ? ".true." : ".false."
  vb = verbose     ? ".true." : ".false."
  write(fnam,
  """
  !! generated with conf.jl
  module ctwmodes_pars
  implicit NONE
  integer, parameter:: single = kind(1.0)
  integer, parameter:: double = selected_real_kind(2*precision(1.0_single))
  integer, parameter:: im = $(im)
  integer, parameter:: km = $(km)
  real(double), parameter:: f0 = $(f0)_double
  logical, parameter:: freesurface = $(fs)
  integer, parameter:: outtype = $(outtype)
  logical, parameter:: verbose = $(vb)
  character(*), parameter:: grid_file = "$(grid_file)"
  character(*), parameter:: Ne_file   = "$(Ne_file)"
  end module ctwmodes_pars
  """ ) # write()
end

# === Testing ===
# This section is for internal use, to test the code.

using Distributions

function testing()
  #topo_fmt = :format1 # (xx, hh) only
  topo_fmt = :format2 # dx is included

  x1 = 20.e3
  x2 = 200.e3
  h1 = 300
  h2 = 2000

  h(x) = (x < x1) ? h1*(x/x1) :
    (x < x2) ? h1 + (h2-h1)*(x - x1)/(x2 - x1) : h2
  xmax = x2 + 500.e3

  #path,io = mktemp()
  test_topo = "tmp.txt"
  open(test_topo; truncate=true) do io
    @show test_topo
    dxc = 10.e3 # just a sample dx
    xx = [20e3 + dxc*i for i in 0:20]
    println(io, length(xx))
    map(x->println(io, "$(x) $(h(x))"), xx)
    println(io, xmax)
    if topo_fmt == :format2
      dx0 = let
        xx0 = [-10e3 + (dxc*1.2)*i for i in 0:150] # format 2
        [(xx0[i+1] - xx0[i]) for i in axes(xx0,1)[begin:end-1]]
      end
      println(io, length(dx0)) # save dx0
      map(x->println(io, "$(x)"), dx0)
    end
  end
  println("testing: a test topo file $(test_topo) generated.")

  (; hh,xx,dx) = slope_from_file(test_topo)
  @show dx
  f = let
    n2 = 10
    n3 = 10
    xxs = [i*x2/n2 for i in -2:n2] # slope
    xxr = [xxs[end] + i*(xmax - x2)/n3 for i in 1:n3] # flat bottom
    xxq = vcat(xxs,xxr)
    hhq = h.(xxq) # original f(x)
    f = lines(xxq,hhq) # original f(x)
    #lines!(xx,hh)
    scatter!(xx,hh; color = Makie.wong_colors()[2])
    f
  end
  display(f)
end

"""
Produce a sample Conf_topo.txt
  See slope_from_file() for the format of the file
"""
function testing_topo()
  #format = :format1
  format = :format2
  dep1 = 150.0
  dep2 = 5500.0
  x1 = 100.0e3
  x2 = 200.0e3
  h_x1 = dep1 / x1
  h_x2 = (dep2 - dep1) / (x2 - x1)
  h(x) =
    (x < x1) ? (h_x1 * x) :
    (x < x2) ? (h_x1 * x1 + h_x2 * (x-x1)) : dep2
    xxs = let dx = sort(rand(Uniform(7.0e3, 20.0e3), 11))
    cumsum(dx)
  end
  xmax = xxs[end] + 100.0e3
  @say "Generating $(conf_topo_file)"
  open(conf_topo_file; truncate=true) do os
    println(os, length(xxs))
    for xx in xxs
      println(os, "$(xx)  $(h(xx))")
    end
    println(os, xmax)
    if (format == :format2)
      println("testing_topo: Format 2: Adding dx")
      dx = sort(rand(Uniform(2.0e3, 21.0e3), 30))
      println(os, length(dx))
      for d in dx
        println(os, "$(d)")
      end
    else
      println("testing_topo: Format 1")
    end
  end # open
end

# === Main: Build config files ===


# I use N_southAus2 and f0
#   from SeaParams.jl but that's not necessary.
# You can just write down your own N(z) funciton here.
# You can just set your own f0 here.
push!(LOAD_PATH, pwd())
#using SeaParams: N_southAus2, totdep, f0
using SeaParams: SeaParams

Nfunc = SeaParams.N_southAus2 # function N(z)

#=
"""
An example of function h(x)
   - A little vertical wall followed by a linear continental slope
     followe by a flat bottom.
   - Undefined for x < 0
"""
function hfunc(x)
  x1 = 100.0e3 # end of continental slope
  h1   = 60.0 # initial vertical wall
  hbot = 5000.0 # dep of flat bottom region.
  h_x = (hbot - h1)/x1
  (x < 0) ? NaN :
    (x < x1) ? (h1 + h_x * x) : hbot
end
=#


# Input files
const conf_topo_file = "Conf_topo.txt"
const conf_Ne_file   = "Conf_Ne.txt"

# Output files
const grid_file = "Grid.txt"
const Ne_file = "Ne.txt"
const par_file = "ctwmodes_pars.f90"

#=
"""
Depth of at x = 0.
  We can set z1 = 0 to remove this initial step.
"""
#z1 = 0.0 # [m]
z1 = -30.0 # [m], initial depth at x = 0.

const xmax = 300e3 # meters
const wid_slope = 150e3 # [m], width of the linear slope.

"""
Sample dz.
  dz[:] is determined in the function on the basis of dz_sample.
"""
dz_sample = 40.0 # [m]
(; dx, dz, kocn) = let
  totdep1 = totdep + z1
  h_x = totdep1 / wid_slope
  grid_ linslope(; z1 = z1, z2=-totdep
                ,dz_sample = dz_sample
                ,h_x, Nfunc = Nfunc, xmax = xmax)
end
# Because there is one step dz[1] at x = 0,
#   there are only km-1 horizontal steps along the slope.
@assert wid_slope ≈ sum(dx[1:(km-1)])
=#

function main()
#  testing_topo()
#  (; dx, dz, kocn) = slope_for_testing()

  # Determine the grid:
  #   If the conf file exists, read it;
  #   if not, use an example.
  (; dx, dz, kocn) = isfile(conf_topo_file) ?
    slope_from_file(conf_topo_file) :
    grid_linslope_example()

  # Write grid_file:
  write_grid(grid_file; dx, dz, kocn)

  # Generate N(z) and write the Ne_file:
  #   If the conf file exists, read it;
  #   if not, use the Nfunc(z) function.
  let zzax = calc_zzax(dz)
    Ne = isfile(conf_Ne_file) ?
      get_Ne_from_file(conf_Ne_file; zzax) :
      grid_Nfunc(Nfunc; zzax)
    write_Ne(Ne_file; Ne, zzax)
  end

  # Write par_file.
  let im = length(dx), km = length(dz)
    write_fort_pars(par_file; im, km, SeaParams.f0, grid_file, Ne_file)
  end
end

main()
