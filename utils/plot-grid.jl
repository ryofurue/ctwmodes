#!/usr/bin/env julia
# Just to visualize the grid configuration,
# read Grid.txt (see conf.jl) and plot it.

const fnam_txt_default = "Grid.txt"
const fnam_nc_default  = "Eigenmodes-z.nc"

#using CairoMakie
using GLMakie
using Match
using NamedPositionals # experimental
using NCDatasets
using OffsetArrays


"""
Copied from conf.jl
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


"""
Copied from conf.jl .
"""
#function read_line_parse(is, typelist)
#  strlist = split(readline(is))
#  res = map(x -> parse(x[1], x[2]), zip(typelist, strlist))
#  Tuple(x for x in res) # Vector -> Tuple
#end
function read_line_parse(is, typelist)
  strlist = split(readline(is))
  res = map(x -> parse(x[1], x[2]), zip(typelist, strlist))
  if length(res) == 0
    error("parse failed")
  elseif length(res) == 1
    return res[1]
  else
    return Tuple(x for x in res) # Vector -> Tuple
  end
end

"Copied from conf.jl"
function read_line_if_exists(is)
  eof(is) && return nothing
  s = strip(readline(is)) # remove white space characters
  return (s != "") ? s : nothing
end


# ---
"""
Copied from `conf.jl`.
If the file format is updated there,
this one has to be updated too.
"""
function get_hh_from_file(fnam)
  @say "Reading (xx[:], hh[:]) from $(fnam)"
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
  #is_increasing(xx0; strict=true)
  #is_increasing(hh0; strict=false)
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
      @say "Format 1 detected."
      #calc dx from cell edges
      #dx = [(xx0[i] - xx0[i-1]) for i in xxax[begin+1:end]]
      @say "x axis isn't extended."
      #(;hh,xx,dx) = extend_domain(;hh=hh0,xx=xx0,dx,xmax)
      (;hh0=hh0,xx0=xx0)
    else # Format 2
      @say "Format 2 detected, but dx is ignored ..."
      #im = parse(Int,line)
      #@say @sshow im
      #dx = [read_line_parse(is, (Float64,)) for i in 1:im]
      #xx = OffsetArray(vcat(0.0, cumsum(dx)), -1)
      #hh = interpolate_helper(;xs0=xx0, ys0=hh0, xs=xx)
      #(;hh=hh,xx=xx,dx=dx)
      (;hh0=hh0,xx0=xx0)
    end
  close(is)
  return ret
end


# ---

"""
Read "Grid.txt"
"""
function read_grid_txt(fnam)
  println("Reading dx, dz, kocn from text file $(fnam) .")
  is = open(fnam,"r")

  im, km = read_line_parse(is, (Int,Int))
  kocn = OffsetArray{Int}(undef, 0:(im+1))
  dx = OffsetArray{Float64}(undef, 0:(im+1))
  dz = OffsetArray{Float64}(undef, 0:(km+1))
  for i in 1:im
    (ii, dx[i], kocn[i]) = read_line_parse(is, (Int,Float64,Int))
    @assert ii == i
  end
  for k in 1:km
    (kk, dz[k]) = read_line_parse(is, (Int,Float64))
    @assert kk == k
  end
  dx[0] = dx[1]
  dx[im+1] = dx[im]
  dz[0] = dz[1]
  dz[km+1] = dz[km]

  close(is)

  dx ./= 1000 # meters -> kilometers
  (dx=dx, dz=dz, kocn=kocn,
   im=im, km=km)
end

function read_grid_nc(fnam)
  println("Reading dx, dz, kocn from netCDF file $(fnam) .")
  (dx, dz, kocn, im, km) = NCDataset(fnam, "r") do ds
    dx0 = copy(ds["dx"][:]) ./ 1000 # meters -> kilometers
    dz0 = copy(ds["dz"][:])
    kocn0 = copy(ds["kocn"][:])
    im = length(dx0) - 2
    km = length(dz0) - 2
    dx   = OffsetArray(dx0, 0:im+1)
    dz   = OffsetArray(dz0, 0:km+1)
    kocn = OffsetArray(kocn0, 0:im+1)
    (dx=dx, dz=dz, kocn=kocn, im=im, km=km)
  end
end

function read_grid(fnam)
  m = match(r".*\.nc$", fnam)
  isnothing(m) ? read_grid_txt(fnam) : read_grid_nc(fnam)
end

"NOT tested"
function calc_iocn(; kocn, im, km)
  iocn = OffsetArray{Int}(undef, 0:(km+1))
  kk = 1
  for i in 1:im
    kbot = kocn[i]
    for k in kk:kbot # no iteration if kk > kbot
      iocn[k] = i
    end
    kk = kbot + 1
  end
  iocn
end

"""
Calculatex vertex coordinates xxs[0:im], zzs[0:km]
"""
function calc_vertices(; dx, dz, im, km)
  @assert axes(dx,1) == 0:im+1
  @assert axes(dz,1) == 0:km+1
  xxs = OffsetArray(vcat(-dx[0], cumsum(dx[0:im+1]).-dx[0]), -1:im+1)
  zzs = OffsetArray(vcat(dz[0], cumsum(-dz[0:km+1]).+dz[0]), -1:km+1)
  (xxs=xxs, zzs=zzs)
end

function celltype(i, k; kocn, iocn, im, km)
  #(k == 7) && (i == 1) && @show i, k, iocn[k], kocn[i+1]
  (k == 0) && return :undefined
  (i == im+1) && return :undefined
  if (k <= km && i == iocn[k] - 1) && (i >= 1 && k == kocn[i] + 1) &&
    (k > 1 && i != iocn[k-1]) && (i < im && k != kocn[i+1])
    return :protruding # ctwmodes-007.f90 will eliminate this
  end
  if (k == km+1)
    if i >= 1  && k == kocn[i] + 1 && i != iocn[k-1]
      return :below_flatbottom
    else
      return :undefined
    end
  elseif (i == 0)
    if k <= km && i == iocn[k] - 1 && k != kocn[i+1]
      return :left_of_vertical
    else
      return :undefined
    end
  elseif i == iocn[k] && k == kocn[i] # corner point
    return :corner
  elseif i == iocn[k] # vertical wall
    return :vertical
  elseif k == kocn[i] # flat bottom
    return :flatbottom
  elseif i == iocn[k] - 1 && k != kocn[i+1]
    return :left_of_vertical
  elseif k == kocn[i] + 1 && i != iocn[k-1]
    return :below_flatbottom
  elseif i > iocn[k]
    return :interior
  else
    return :undefined
  end
end


function slope_bottom(; dx, dz, kocn, iocn, im, km)
  pnts = [ (xx = 0.0, zz =0.0) ] # initial point
  # Scan all cell centers. Inefficient but easy to program.
  for k in 1:km
    for i in 1:im
      if i != iocn[k] && k != kocn[i]; continue; end
      @match celltype(i,k; kocn, iocn, im, km) begin
        #if i == iocn[k] && k == kocn[i] # corner point
        :corner => begin
          #@show i, k
          xx = pnts[end].xx + dx[i]
          zz = pnts[end].zz - dz[k]
          push!(pnts, (xx=xx,zz=zz))
        end
      #elseif i == iocn[k] # vertical wall
        :vertical => begin
          #@show i, k, kocn[i]
          xx = pnts[end].xx
          zz = pnts[end].zz - dz[k]
          push!(pnts, (xx=xx,zz=zz))
          break # i loop
        end
      #elseif k == kocn[i] # flat bottom
        :flatbottom => begin
          #@show i, k, iocn[k]
          xx = pnts[end].xx + dx[i]
          zz = pnts[end].zz
          push!(pnts, (xx=xx,zz=zz))
        end
        # otherwise don't advance.
      end#match
    end
  end
  Tuple.(pnts) # drop keys
end

# ==========================================
function plot_grid!(ax; dx, dz)
  lineatts = (linewidth=0.5, color=:red, linestyle=:dash)
  xx0 = -dx[0]
  zz0 = dz[0]
  xx1 = xx0 + sum(dx)
  zz1 = zz0 - sum(dz)

  #--- hor. grid lines ---
  lines!(ax, [(xx0,zz0), (xx1,zz0)]; lineatts...)
  zz = zz0
  for k in axes(dz,1)
    zz -= dz[k]
    lines!(ax, [(xx0,zz), (xx1,zz)]; lineatts...)
  end

  #--- ver. grid lines ---
  lines!(ax, [(xx0,zz0), (xx0,zz1)]; lineatts...)
  xx = xx0
  for i in axes(dx,1)
    xx += dx[i]
    lines!(ax, [(xx,zz0), (xx,zz1)]; lineatts...)
  end

  #--- ocean frame
  let
    im = length(dx) - 2
    xx0 = 0.0; xx1 = sum(dx[1:im])
    km = length(dz) - 2
    zz0 = 0.0; zz1 = - sum(dz[1:km])
    lines!(ax, [(xx0,zz0), (xx0,zz1), (xx1,zz1), (xx1,zz0), (xx0,zz0)];
           linewidth=0.7, color=:blue, linestyle=:dash)
  end
end

function plot_cell!(ax,i,k; xxs, zzs, atts)
  ps = Point2f[(xxs[i-1], zzs[k-1]), (xxs[i-1], zzs[k])
               , (xxs[i], zzs[k]), (xxs[i], zzs[k-1])]
  poly!(ax, ps; atts...)
end

function plot_bottom_cells!(ax; kocn, im, xxs, zzs)
  atts = (color = :red, strokecolor = :gray, strokewidth = 1)
  for i in 1:im
    k = kocn[i]
#    ps = Point2f[(xxs[i-1], zzs[k-1]), (xxs[i-1], zzs[k]), (xxs[i], zzs[k]), (xxs[i], zzs[k-1])]
#    poly!(ax, ps; atts...)
    plot_cell!(ax,i,k; xxs, zzs, atts)
  end
end

function plot_slope_cells!(ax; kocn, iocn, xxs, zzs, im, km)
  color4 = 0.5 * Makie.to_color(:red) + 0.5 * Makie.to_color(:blue)
  atts1 = (color = (:blue,  0.5),)# strokecolor = :gray, strokewidth = 1)
  atts2 = (color = (:red,   0.5),)# strokecolor = :gray, strokewidth = 1)
  atts3 = (color = (:green, 0.5),)# strokecolor = :gray, strokewidth = 1)
  atts4 = (color = color4, )# strokecolor = :gray, strokewidth = 1)
  for k in 0:km+1
    for i in 0:im+1
      @match celltype(i,k; kocn, iocn, im, km) begin
        :left_of_vertical => plot_cell!(ax,i,k; xxs, zzs, atts=atts1)
        :below_flatbottom => plot_cell!(ax,i,k; xxs, zzs, atts=atts2)
        :corner           => plot_cell!(ax,i,k; xxs, zzs, atts=atts3)
        :protruding       => plot_cell!(ax,i,k; xxs, zzs, atts=atts4)
        _ => nothing
      end#match
      #(k == 7) && (i == 1) && @show celltype(i,k; kocn, iocn, im, km)
    end
  end
end

function plot_indices!(ax; dx, dz)
  zz = dz[0]
  for k in axes(dz,1)[begin:end-1]
    if (k%5 == 1 || k == lastindex(dz)-1)
      z = zz - dz[k]/2
      text!(ax, -1.3*dx[0], z; text = string(k), align=(:right,:center))
    end
    zz -= dz[k]
  end
  xx = -dx[0]
  for i in axes(dx,1)[begin:end-1]
    if (i%5 == 1 || i == lastindex(dx)-1)
      x = xx + dx[i]/2
      text!(ax, x, 1.3*dz[0]; text = string(i), align=(:center,:baseline))
    end
    xx += dx[i]
  end
end

#---
function determine_fnam()
  isfile(fnam_nc_default) && return fnam_nc_default
  isfile(fnam_txt_default) && return fnam_txt_default
  error("Neither $(fnam_nc_default) or $(fnam_txt_default) exits.")
end

fnam = @match length(ARGS) begin
  0 => determine_fnam()
  1 => ARGS[1]
  _ => error("Give one or zero input file name.")
end

(; dx, dz, kocn, im, km) = read_grid(fnam)
iocn = calc_iocn(; kocn, im, km)
(xxs,zzs) = calc_vertices(; dx, dz, im, km)
pnts = slope_bottom(; dx, dz, kocn, iocn, im, km)

const conf_topo = "Conf_topo.txt"
(; xx0, hh0) = isfile(conf_topo) ?
  get_hh_from_file(conf_topo) : (nothing,nothing)

fig = Figure()
ax = Axis(fig[1,1]
          ;xlabel="x [km]", ylabel="z [m]"
          ,xgridvisible = true, xgridcolor = :gray
          ,xgridwidth = 0.3, xgridstyle = :solid
          ,ygridvisible = true, ygridcolor = :gray
          ,ygridwidth = 0.3, ygridstyle = :solid
          ,xticks=-100:50:1000
          ,yticks=500:-500:-7000
           )
plot_grid!(ax; dx, dz)
#plot_bottom_cells!(ax; kocn, im, xxs, zzs)
plot_slope_cells!(ax; kocn, iocn, xxs, zzs, im, km)
lines!(ax, pnts; linewidth=2, color=:gray20)
#scatter!(ax, pnts)
if !isnothing(hh0)
#  @show xx0
#  @show hh0
  col = Makie.wong_colors()[1]
  lines!(  ax, xx0 ./ 1000, -hh0; color=(col,0.5), linewidth=4)
  scatter!(ax, xx0 ./ 1000, -hh0; color=col)
end
plot_indices!(ax; dx, dz)

display(fig)
#save("tmp-nc.pdf", fig)
