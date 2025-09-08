"""
This module is not an essential part of ctwmodes.
It just provides ad-hoc parameters and functions
such as f0, N(z), . . .
"""
module SeaParams

"Standard depth in the flat-bottom region"
const totdep = 5500.0

"Coriolis parameter"
const f0 = 0.0001 # [1/s]

"Determine N so that c1 = 2.65 m/s"
const Nconst =
  let c1_samp = 2.65 # m/s
    N = c1_samp * pi / totdep
  end

"""
A constant N(z) function.
Uses a different, simpler, standard value.
"""
function N_const2( _ ) # N(z) but z is not used
  return 3.0e-3
end

"""
Simple exponential profile that fits the south Australian
profile Lei gave me.
"""
function N_southAus(z)
  bvf0 = 4.5e-3 # 1/s, N(0)
  bvf1 = 0.3e-3 # 1/s, N(-maxdep)
  scale = 1700.0 # exponential scale
  maxdep = 5500
  func0(z) = exp((maxdep + z)/scale) - 1 # func0(z=-maxdep) = 0
  return (bvf0 - bvf1) * func0(z) / func0(0) + bvf1
end

"""
Simple exponential profile that fits observed N(z) profiles
  near (142E, 40S) (off southeastern Australia).
  Similar to N_southAus(z) in SeaParams.jl, but adjusted
  to fit WOA18 1deg x 1deg annual climatology over 2005-2017.
  It fits well at (142E, 40S) and the fitting is good
  in the surrounding area, too. I looped over 140E to 142E and 42S to 38S
  and made plots to see that. The fitting was purely *subjective* (by eye).
"""
function N_southAus2(z)
  bvf0 = 4.5e-3 # 1/s, N(0)
  bvf1 = 0.5e-3 # 1/s, N(-maxdep)
  scale = 1700.0 # exponential scale
  maxdep = 5500
  func0(z) = exp((maxdep + z)/scale) - 1 # func0(z=-maxdep) = 0
  return (bvf0 - bvf1) * func0(z) / func0(0) + bvf1
end


end
