#!/usr/bin/env julia
# =============================================================================
# Tool Name:    CTW data conversion script
# Description:  Script to convert the output of ctwmodes.f90
#                 to a netCDF file.
# Author:       Ryo Furue <ryofurue@gmail.com>
# License:      MIT License
# Version:      v"0.7"
# Created:      2025-??-??
# Updated:      2025-09-04
# =============================================================================
const SCRIPT_VERSION = v"0.7" # `Base.VERSION` is used by the interpreter.

println("CTW data conversion script, version ", SCRIPT_VERSION)

"""
*.bin      is the "stream" binary.
*.fort.bin is the Fortran sequential binary.
"""
const binext     = "bin"
const fortbinext = "fort.bin"
function filetype(fnam)
  if     ! isnothing(match(r".*\.fort\.bin$", fnam))
    :sequential
  elseif ! isnothing(match(r".*\.bin$"      , fnam))
    :stream
  else
    :unkown
  end
end

"""
Default names for the input and output files.
You can specify the names by command line options. See the usage.
"""
const filestem = "Eigenmodes-z"
const inf_def1 = "$(filestem).$(binext)"
const inf_def2 = "$(filestem).$(fortbinext)"
const ouf_def = "$(filestem).nc"

# --- Determine input and output files ---
#  This section is a bit complicated to provide
#  multiple ways of specifying input or output filenames.
using ArgMacros
using Match

function myerror(exitcode,mes)
  println(stderr,mes)
  exit(exitcode)
end

const usage =
  "fortbin2nc [-i infile] [-o outfile] [infile [outfile]]"

#We used to have '--double' or '-d' ... both input and output are double precition.
args = @tuplearguments begin
  @helpusage """
    \n\tfortbin2nc [-i infile] [-o outfile]
    \tfortbin2nc [infile [outfile]]
    \tfortbin2nc -i infile [outfile]
    \tfortbin2nc -o outfile [infile]"""
  @helpdescription """
  Also supported is
    '--title=<title>' or '-t <title>'

  If the extension of the input file is ".bin",
    a raw (stream) binary file is assumed; if it's "fort.bin",
    a Fortran sequential binary file is assumed.
    """
#  @argumentflag     dble  "-d" "--double"
  @argumentoptional String inf1  "-i" "--input"
  @argumentoptional String ouf1  "-o" "--output"
  @argumentoptional String title "-t" "--title"
#  @argumentflag verbose "-v"
  @positionaloptional String par1
  @positionaloptional String par2
end


"""
Determine input and output files from the command line parameters.
"""
function pick_files(inf1,ouf1,par1,par2;inf_def1, inf_def2, ouf_def)
  bad(x)  = isnothing(x)
  good(x) = ! bad(x)
  inf = good(inf1) ? inf1 : nothing
  ouf = good(ouf1) ? ouf1 : nothing

  if (good(inf) && good(ouf)) # if both are determined,
    if good(par1) || good(par2)
      myerror(1,"extra parameter(s): $(par1), $(par2)")
    else
      return (inf,ouf) # -> RETRUN
    end
  end

  if good(inf) # means bad(ouf)
    ouf = par1
  else
    inf = par1
  end

  if bad(inf)
    if isfile(inf_def1)
      inf = inf_def1
    elseif isfile(inf_def2)
      inf = inf_def2
    else
      myerror(2, "Neither $(inf_def1) or $(inf_def2) exists.")
    end
    isfile(inf_def1) && isfile(inf_def2) &&
      println("## Both $(inf_def1) and $(inf_def2) exist. I've taken the former.")
  end

  if bad(ouf)
    ouf = par2
  else
    good(par2) && myerror(3,"infile=$(inf), outfile=$(ouf), extra parameter $(par2)")
  end

  if bad(ouf)
    ouf = ouf_def
  end

  return (inf,ouf)
end

# --- determine inf and ouf ---
(inf,ouf) = pick_files(args.inf1, args.ouf1, args.par1, args.par2;
                       inf_def1, inf_def2, ouf_def)

println("infile=$(inf); outfile=$(ouf)")

# --- determine whether inf is stream binary or not.
is_streambin = @match filetype(inf) begin
  :stream => true
  :sequential => false
  _ => myerror(5,"$(inf): unknown file type.")
end

# --- check the files ---
"""
Check if the owner has write permission.
"""
function has_write_permission(filepath) # ChatGPT gave me this:
  file_stat = stat(filepath)
  return (file_stat.mode & 0o200) != 0
end

(! isfile(inf)) && myerror(5, "$(inf) doesn't exist.")

if isfile(ouf) && (! has_write_permission(ouf))
  myerror(7, "$(ouf) exists but doesn't have write permission.")
end

# === read, sort, and save ===
using OffsetArrays
push!(LOAD_PATH, pwd())
using CTW_data_tools

#intype  = args.dble ? (intype  = Float64, ) : ()
#outtype = args.dble ? (outtype = Float64, ) : ()

(; num, alphaR, alphaI, pr, dx, dz, f0, bvf2e, kocn, iocn) =
  read_fort_data(inf; is_streambin) # ; intype...

@show kocn
@show iocn

(cee_west, idx_west, cee_east, idx_east) = sort_modes(alphaR, alphaI)
let nc = length(cee_west), km  = length(dz) - 2
  if nc != km
    println(stderr, "nc = $(nc) != km = $(km)")
  else
    println("nc = km = $(nc)")
  end
end

(xax, zax) = calc_grid(; dx=dx, dz=dz)

outtype = eltype(pr)

nidx = idx_west
cee = cee_west

prij = map_to_xz(; pr=pr, nidx=nidx, num=num, outtype=outtype)

# normalization is optional
for m in axes(prij,3)
  normalize!(view(prij,:,:,m); mes = "mode $(m): ")
end


(pslope, pbot) = let km = length(dz) - 2, im = length(dx) - 2
  nan = eltype(prij)(NaN)
  pslope = OffsetArray{eltype(prij),2}(undef, axes(prij,2), axes(prij,3))
  pbot   = OffsetArray{eltype(prij),2}(undef, axes(prij,1), axes(prij,3))
  for m in axes(prij,3)
    v = slope_vector(view(prij,:,:,m); iocn=iocn, kocn=kocn, kmax=km)
    v2 = vcat(nan, v, nan)
    for (k1, k2) in zip(axes(pslope,1), axes(v2,1))
      pslope[k1,m] = v2[k2]
    end
    v = bottom_vector(view(prij,:,:,m); iocn=iocn, kocn=kocn, imax=im)
    v2 = vcat(nan, v, nan)
    for (i1, i2) in zip(axes(pbot,1), axes(v2,1))
      pbot[i1,m] = v2[i2]
    end
    OffsetArrays.no_offset_view(pbot[:,m]) .= vcat(nan, v, nan)
  end
  (pslope,pbot)
end

#println(stderr, "### Normal output temporarily omitted ###")
save_modes(ouf; prij=prij, cee=cee
           ,pslope=pslope, pbot=pbot
           ,xax=xax, zax=zax, dx=dx, dz=dz, f0=f0, bvf2e=bvf2e
           ,kocn=kocn, iocn=iocn
           ,outtype=outtype
           ,title=args.title)

modeax = axes(prij,3)
#@show typeof(modeax)

#---- map to the (dx/2, dz/2) vertex grid.
(; xx, zz, pp) = remap_diag(view(prij,:,:,1); dx, dz, iocn, kocn)
prij2 = OffsetArray{eltype(prij)}(undef
  , axes(pp,1), axes(pp,2), modeax)
prij2[:,:,modeax[1]] .= pp

for m in modeax[2:end]
  local xx, zz, pp
  (; xx, zz, pp) = remap_diag(view(prij,:,:,m); dx, dz, iocn, kocn)
  prij2[:,:,m] .= pp
end

(stem,_) = splitext(ouf)
ouf_finer = stem * "-finer" * ".nc"
save_modes(ouf_finer; prij=prij2, cee=cee
           ,xax=xx, zax=zz, f0=f0
           ,outtype=outtype
           ,title=args.title)

#-- This section is optional: Examine the unphysical modes --
found = false
for m in axes(cee_east,1)
  global found
  if abs(cee_east[m]) > 1.0e-14
    println(stderr, "cee_east[$(m)] = $(cee_east[m])")
    found = true
  end
end
(! found) && println("All extra eigenvalues are essentially 0.")
