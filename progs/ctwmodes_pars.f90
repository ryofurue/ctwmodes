!! generated with conf.jl
module ctwmodes_pars
implicit NONE
integer, parameter:: single = kind(1.0)
integer, parameter:: double = selected_real_kind(2*precision(1.0_single))
integer, parameter:: im = 30
integer, parameter:: km = 17
real(double), parameter:: f0 = 0.0001_double
logical, parameter:: freesurface = .true.
integer, parameter:: outtype = single
logical, parameter:: verbose = .true.
character(*), parameter:: grid_file = "Grid.txt"
character(*), parameter:: Ne_file   = "Ne.txt"
end module ctwmodes_pars
