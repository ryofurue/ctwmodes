# Numerical calculation of coastal trapped wave modes

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17507173-blue)](https://doi.org/10.5281/zenodo.17507173)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ryofurue/ctwmodes/blob/main/LICENSE)
<!--
[![Build Status](https://github.com/ryofurue/ctwmodes/actions/workflows/ci.yml/badge.svg)](https://github.com/ryofurue/ctwmodes/actions)
-->


Provides a program to calculate baroclinic coastal trapped (CTW) wave modes
in the low-frequency limit (ω ≪ f).
Currently, we are writing a paper to fully describe our numerical method.
We will update this document when we have made the paper public.
<!--
See [Furue & Tanaka (2025)?????](???? No URL yet ?????)
for a full description and discussion
of our numerical method.
-->


<!-- These links are invisible in the rendered version. -->
[progs]: progs
[utils]: utils
[project]: project

[utils/conf.jl]:       utils/conf.jl
[utils/fortbin2nc.jl]: utils/fortbin2nc.jl
[utils/plot-grid.jl]:  utils/plot-grid.jl

[utils/Conf_topo.txt]: utils/Conf_topo.txt
[utils/Conf_Ne.txt]: utils/Conf_Ne.txt

[progs/ctwmodes.f90]: progs/ctwmodes.f90
[progs/ctwmodes_pars-sample.f90]: progs/ctwmodes_pars-sample.f90
[progs/Grid.txt]: progs/Grid.txt
[progs/Ne.txt]: progs/Ne.txt

## 1. Description

[Brink & Chapman (1987)](https://www.whoi.edu/cms/files/Fortran_30425.htm)
published their Fortran program,
`bigload4.for`,
to calculate baroclinic CTW modes.

Comparison of the two programs:

- **BC87** is much more versatile:
High frequencies are fine, bottom friction can be included,
h(x) doesn't have to be monotonic.

- **BC87**: The user gives an initial guess of an eigenvalue and other parameters. The program searches a solution around the initial guess by iteration. If the iteration converges, the user obtains a CTW mode.

- **Our program** gives all the CTW modes at once and is fast.

- **BC87** uses the sigma coordinates; solutions are likely more accurate when the slope is gentle.

- **Our program** uses the z coordinates and therefore it works
well for very steep slopes. When the "slope" is purely vertical,
it gives the Kelvin-wave modes.

- **Our program**: When there are `km` gridpoints in the vertical,
  there are `km` modes, which are guaranteed to satisfy an
  extremely simple orthogonality relation.

Yuki Tanaka wrote the original Fortran program
using the sigma coordinates
([Tanaka 2023](https://doi.org/10.1175/JPO-D-22-0201.1))
and modified it to use the z coordinates.
The latter version has evolved into the program we provide here.


## 2. Getting Started

### Dependencies

* A Fortran compiler: The main solver is written in Fortran.
* `LAPACK` the linear algebra library:
   Used in the Fortran program. (`LAPACK` comes with macOS.)
* (optional) Julia interpreter: Various supporting tools
  are written in the Julia language.

Notes

* On Mac, `LAPACK` comes with the operating system and
  can be used by `gfortran` without any special configuration.
  (See the `progs/Makefile`.)
* On Mac or Linux, `juliaup` installs the Julia environment
  under the user's home directory. No superuser privilege is needed.

### Quick Start

1. [Download the repository as a ZIP](https://github.com/ryofurue/ctwmodes/archive/refs/heads/main.zip)
   (or click the green **Code** button above → **Download ZIP**).
2. Unzip it.
3. Enter the project directory:

    ```shell
    $ cd ctwmodes-main
    $ cd project
	```

4. (Linux) Modify `Makefile` in such a way that `LAPACK` is linked.
5. Build and run:

    ```shell
    $ make       #(1) Build the solver.
    $ ./ctwmodes #(2) Run the solver.
    ```
will produce a solution using the sample topography h(x) and the sample
stratification N(z) in the `project` directory.

4. Optionally, you can convert the binary output to a netCDF file:

    ```shell
    $ bash ../utils/install-julia.bash #(3) Install Julia and some packages.

    . . . Restart the terminal to activate new PATH . . .

    $ julia fortbin2nc.jl              #(4) Convert binary -> netCDF.
    ```
    You need to do step (3) only once; next time, you just skip it.

    > ℹ️ The script `utils/install-julia.bash` is provided to install
    > a Julia environment and then install the packages our scripts use.
	> The julia installer will modify your PATH, to activate which
	> you can just restart you terminal.
    > This script isn't essential at all.
    > We provide
    > [a brief tutorial](#appendix-c-installing-and-using-the-julia-interpreter)
	> to install Julia and the packages.

5. To change the parameters or configuration, edit
[`ctwmodes_pars.f90`](#specification-of-ctwmodes_pars.f90),
[`Grid.txt`](#specification-of-grid.txt), and
[`Ne.txt`](#specification-of-ne.txt) under the `project` directory.
You can do this
by [using the script `conf.jl`](#optional-generating-the-configuration-files).

### Notes on the units
Among the input parameters,
only `dx[:]`, `dz[:]`, `N[:]`, and `f0` are dimensioned.
Use meters for `dx` and `dz`
and radian/seconds for `N` and `f0`.

This restriction to MKS units is solely due to the hard-coded
gravitational constant g in [progs/ctwmodes.f90].
(In the rigid-lid mode, this restriction is lifted.)


## 3. Details

The above [Quick Start](#quick-start) shows a minimal procedure,
where all files are placed in the [project] directory.

The subsequent explanation uses our original directory structure
and does not assume that the user did the "Quick Start".

<!--
But, we do not necessarily recommend putting everything in one directory.
This github repository reflects our original directory structure:
the Fortran solver resides in the [progs] directory,
the other supporting tools in the [utils] directory,
and the [project] directory just houses symlinks to the files.
The following explanation uses [progs], [utils], and a third directory where you run the solver.
-->


### Installing & preparing

Solver

* The [progs] directory houses the solver.
* (Linux/Unix) Install `LAPACK` or make sure it is available to the Fortran compiler you use.
* The current `Makefile` assumes `gfortran`. If you use a different Fortran compiler, edit `Makefile` accordingly, where you may needed to alter options related to `LAPACK`.
* If you are on Unix/Linux, you may need to know how to link `LAPACK` and modify `Makefile` accordingly.


(Optional) Tools

* The [utils] directory houses several tools.
*  Install [`juliaup`](https://github.com/JuliaLang/juliaup) and run it. It will create a Julia environment under your home directory (Mac and Linux). You may need to add `julia` (the Julia interpreter) to your command search path `PATH`. See [the tutorial](#appendix-c-installing-and-using-the-julia-interpreter).
* Our Julia programs use various packages. See below for how to install them.
* Some of the tools ([utils/fortbin2nc.jl] and [utils/conf.jl]) are designed to run as a command
and so you may want to give them an `x` (executable) permission.

<!--
### Running the program
## 3. Running
-->



### (Optional) Generating the configuration files
You can easily write or edit
the configuration files, `ctwmodes_pars.f90`, `Grid.txt`, and `Ne.txt`,
as their formats are simple
([Appendix A](#appendix-a-parameter-file-formats)).


But, if you want to create these files from observational data of N(z) or real topography h(x) or from functional forms of h(x) and N(z), you might want to use [utils/conf.jl].

`conf.jl` reads a text file `Conf_topo.txt`
describing a discrete h(x) profile
and another text file `Conf_Ne.txt`
describing a discrete N(z) profile.
The formats of these files are
described [below](#hx-input-file-for-conf.jl)
and [below](#nz-input-file-for-conf.jl).
Examples are provided as [utils/Conf_topo.txt]
and [utils/Conf_Ne.txt].

If the corresponding data file does not exist,
`conf.jl` uses a continuous function for h(x) or N(z).
Currently, `conf.jl` uses functions defined in
`SeaParams.jl`.
So, you can add your function to this module
edit `conf.jl` to use it.
But `SeaParams.jl` isn't an essential part of
the tools at all.
You can directly define your function in `conf.jl`.



### (Optional) Visualizing the topography and grid

[utils/plot-grid.jl] reads `Grid.txt` or the netCDF file
(generated by [utils/fortbin2nc.jl])
and graphically shows the topography and the grid configuration. This may be useful to verify your configuration:
```shell
$ julia
julia> include("plot-grid.jl")
```
Either `Grid.txt` or the netCDF file must be in the current working directory.
Like `fortbin2nc.jl`, you'll have missing-package errors first time you run this.
See [the explanation of `fortbin2nc.jl`](#optional-converting-the-binary-output-to-a-netcdf-file).


### Running the solver
Prepare the three parameter files (all text): `ctwmodes_pars.f90`, `Grid.txt`, and `Ne.txt`.

* Copy `ctwmodes_pars.f90` into the `progs` directory, where the Fortran source code resides.
* Copy `Grid.txt` and `Ne.txt` into the directory where you run the solver.

The formats of parameter files are described in
[Appendix A](#appendix-a-parameter-file-formats).
Examples are provided as
[progs/ctwmodes_pars-sample.f90], [progs/Grid.txt], and [progs/Ne.txt].

Make the executable
```shell
$ cd progs
$ make
```

Optionally, copy the executable `ctwmodes` to the directory you want to run the solve in.
By default, `Grid.txt` and `Ne.txt` need to be in the same directory.
(Alternatively, the paths of them can be specified in `ctwmodes_pars.f90`.
See [the explanation of the file](#specification-of-ctwmodes_pars.f90).)


Run the Fortran program
```shell
$ ./ctwmodes
```

The format of the output binary file is explained below
[below](#binary-output-format).

### (Optional) Converting the binary output to a netCDF file
To convert the binary output from the solver, run the script [utils/fortbin2nc.jl].
If it's in the same directory as the binary file (and if it's made executable):
```shell
$ ./fortbin2nc.jl
```
or on the interpreter:
```shell
$ julia
. . . Julia's logo . . .

julia> include("fortbin2nc.jl")
```

First time you run the program, you'll get errors for missing packages.
To resolve them,

* You first run the script within the interpreter (the second method above).
* On each missing package, you press the `]` key to enter the `pkg` mode and `add` the missing package there.
* Press `^C` (control-C) to return to the `julia>` prompt.

Repeat this until you stop seeing missing package errors.

There are some optional features, which you can learn by
```shell
$ ./fortbin2nc.jl --help
```


<!--
## Help

To run a Julia script
```
$ julia thescript.jl
```
-->

## 4. Project information

### Authors
* Ryo Furue
[<ryofurue@gmail.com>](mailto:ryofurue@gmail.com)
[![ORCiD](https://img.shields.io/badge/ORCID-0000--0003--0563--1189-A6CE39)](https://orcid.org/0000-0003-0563-1189)
* Yuki Tanaka

<!--
[![email](https://img.shields.io/badge/email-ryofurue%40gmail.com-blue)](mailto:ryofurue@gmail.com)
-->


### Version History

* 0.7.1: Cleaned up for public release.
* 0.7.0: Complete for internal use.

### License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

<!--
### Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
-->

### References
* Tanaka Y. 2023. Energy conversion rate from subinertial surface tides to internal tides. *J. of Phys. Oceanogr.* 53 (5): 1355–1374.
[https://doi.org/10.1175/JPO-D-22-0201.1](https://doi.org/10.1175/JPO-D-22-0201.1) .
* Furue R. & Tanaka Y. 2025. ???????????????????

---
## Appendix A: Parameter file formats

### Binary output format
(If you don't use Fortran,
you can still read the output from `ctwmodes.f90`
because it's a "stream" ("plain", "raw") binary.)


The best way to specify the format is to copy
the write statements from [progs/ctwmodes.f90]:
```fortran
  integer:: M, im, km
  integer:: num(0:im+1,0:km+1)
  integer:: kocn(0:im+1)
  integer:: icon(0:km+1)
  integer:: outtype   ! usually 4 or 8
  real(?):: alphaR(M) ! c(M)
  real(?):: alphaI(M) ! == 0
  real(?):: VR(M,M)   ! CTW modes
  real(?):: dx(0:im+1)
  real(?):: dz(0:km+1)
  real(?):: f0
  real(?):: N2e(0:km)
  open(newunit=uni, file=ofile &
	,access="STREAM", form="UNFORMATTED" &
    ,status="REPLACE", position="REWIND")
  write(uni) M, im, km
  write(uni) num
  write(uni) kocn
  write(uni) iocn
  write(uni) outtype
  write(uni) real(alphaR, kind=outtype)
  write(uni) real(alphaI, kind=outtype)
  write(uni) real(VR,     kind=outtype)
  write(uni) real(dx,     kind=outtype)
  write(uni) real(dz,     kind=outtype)
  write(uni) real(f0,     kind=outtype)
  write(uni) real(N2e,    kind=outtype)
```
* `M`: the number of grid cells within the "ocean"
[ M = im × km − (# of landcells) ]. Equal to the
total number of modes, physical or spurious.
* `im`, `km`: grid sizes.
See [the explanation of `ctwmodes_pars.f90`](#specification-of-ctwmodes_pars.f90).
* `num`: num[i,k] is the "serial number" of the grid cell (i,k).
See below.
* `kocn`: the k index of the deepest ocean grid cell at location i.
* `iocn`: the i index of the leftmost ocean grid cell at depth k.
* `outtype`: (usually) 4 or 8, specifying the precision of output.
* `alphaR`, `alphaI`: Real and imaginary parts of the eigenvalue c(n).
  Currently, `alphaI` == 0 for all the modes.
* `VR`: `VR[ik,n]` is the value of F(x,z) of mode n
   at the grid cell `ik`, where ik = 1,...,M and n = 0,...,M-1.
* `dx`: `dx[i]` is the horizontal width of the i-th cell
     with i = 1,...,im.
* `dz`: `dz[k]` is the height of the k-th cell.
     with k = 1,...,km.
* `f0`: the constant Coriolis parameter.
* `N2e`: `N2e[k]` is the N^2 value at the cell edge below cell k,
  with k = 1,...,km, except `N2e[0]` is the value at z = 0.


`num` is the map from (i,k) indices to the 1d (serialized) indices;
that is, F[num[i,k]] == F(i,k). Therefore, you can recover
the 2d array of the mode function by
```fortran
  real(double):: F2d(im,km,mmax) ! 2d fields of modes
  F2d = ieee_value(1.d0, IEEE_QUIET_NAN) ! Set to NaN
  do n=1,mmax
  do k=1,km
  do i=1,im
    ik = num(i,k)
	if (ik < 0) cycle ! cell is not ocean.
    F2d(i,k,n) = VR(num(i,k), n) ! VR(ik,n) is F(i,k) of mode n.
  enddo; enddo; enddo
```
Note that `num(i,k) < 0` if the cell (i,k) is not ocean
(not within the ocean).
The `F2d(i,k,n)` array element
does not get assigned for non-ocean cells.

It's proven that there are exactly `km` physical modes
and the remaining `M - km` modes are all unphysical with zero
eigenvalues (`alphaR(n)` == 0).

**Notes** on the binary format. By default, the output is "stream" binary:
it's just a linear sequence of binary values.
To switch to Fortran's "sequential" binary,
set
```fortran
   logical, parameter:: streambin = .false. ! output is sequential binary.
```
in [progs/ctwmodes.f90].
The filename of a stream binary ends in `.bin` whereas
that of a sequential binary ends in `.fort.bin`.



### Specification of `ctwmodes_pars.f90`
This small Fortran module provides constants to the main program.
The most important constants are `im` and `km`,
the grid sizes in the x and z directions.

An example is provided as [progs/ctwmodes_pars-sample.f90].
Here it is reproduced for explanation:
```fortran
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
```
* `single` and `double`:
Throughout the program, single and double precision real variables
are declared as `real(single):: r` and `real(double):: d`.
* `im` and `km`: the numbers of grid cells in the x and z directions.
* `f0`: constant Coriolis parameter
* `freesurface`:
If set to `.false.`, a rigid-lid approximation will be used,
but this isn't recommended because the program runs much more slowly.
See Furue & Tanaka (2025).
* `outtype`: Determines whether the binary output is in single or double
precision. If you further analyze the result,
you might prefer double precision.
* `verbose`: If set to `.true.`, a lot of messages will be printed.
* `grid_file` and `Ne_file`: the paths of the grid and Ne files.
Their usual names are `Grid.txt` and `Ne.txt`.
See the specifications of [`Grid.txt`](#specification-of-grid.txt)
and [`Ne.txt`](#specification-of-ne.txt)


### Specification of `Grid.txt`
An example is provided as [progs/Grid.txt].
```
im km
1   dx[1]   kocn[1]
2   dx[2]   kocn[2]
. . .
im  dx[im]  kocn[im]
1   dz[1]
2   dz[2]
. . .
km  dz[km]
```
* `im` and `km`: the numbers of grid cells in the x and z directions.
These numbers have to agree with those in
[`ctwmodes_pars.f90`](#specification-of-ctwmodes_pars.f90).
* `dx[i]`: the width of the i-th grid cell.
* `dz[k]`: the height of the k-th grid cell.
* `kocn[i]`: the "k" index of the deepest "ocean" grid cell
at horizontal position i.


### Specification of `Ne.txt`
An example is provided as [progs/Ne.txt].

This file lists N(z) values at cell edges:
```
nn
zz[0]     N(zz[0])
zz[1]     N(zz[1])
. . .
zz[nn-1]  N(zz[nn-1])
```
This data file must use the grid scheme of the `ctwmodes` solver:

* `nn`: the number of the points. Must satisfy `nn == km + 1`
  because N is defined at cell edges.
* `zz[k]`: the vertical coordinate (zz < 0). Must satisfy
    - `zz[0] == 0`: sea surface
    - `zz[k] == zz[k-1] - dz[k]`: cell edges

In other words, the `zz` coordinates here must agree
with `dz` specified in [`Grid.txt`](#specification-of-grid.txt).


## Appendix B: Input files to `conf.jl`

### h(x) input file for `conf.jl`
`conf.jl` reads a topography file `Conf_topo.txt`.
It allows for two formats.
Read the description within the script, but here is minimal information.

**Input file format 1:**
```
nn
xx[0]    h(xx[0])
xx[1]    h(xx[1])
. . .
xx[nn-1] h(xx[nn-1])
xmax
```
where `nn` is the number of the points where the pair (x, h(x)) is defined
and `xmax` is the position of the right edge beyond the last topography point
`xx[nn-1]`. Between `xx[nn-1]` and `xmax` a flat bottom of `h(xx[nn-1])` is assumed.

Before processing,
the `xx` coordinates (including `xmax`) are first shifted so that
the 0th position coincides with x = 0.

The grid spacings for the CTW calculation will be determined by the intervals
between `xx`'s.

**Input file format 2** extends format 1:
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
where `im` is the number of grid cells in the x direction
and `dx` are the cell widths.
Before processing,
the `xx` coordinates (including `xmax`) are first shifted so that
the 0th position coincides with x = 0.

In this case `h(x)` values are interpolated onto the gridpoints
defined by `dx` starting from `x = 0`.
If the gridpoint falls outside the range of (the shifted) `xx`'s,
the h value
at the closest defined point is used (that is a "constant extrapolation").
This extrapolation can happen only to the right of `xx[n-1]`
because `xx[0]` is shifted to `x = 0`,
and it will create the flat-bottom region, assuming
that the grid created by `dx` extends beyond `x = xx[nn-1]`.

Currently, `xmax` is ignored for Format 2 because `dx` determines the grid.
We could extend `conf.jl` so that the grid determined
by `dx` be automatically extended to reach `xmax`.


### N(z) input file for `conf.jl`

`conf.jl` reads a N(z) data file `Conf_Ne.txt`.
The file format is
```
nn
dep[1]  Ne(-dep[1])
dep[2]  Ne(-dep[2])
. . .
dep[nn] Ne(-dep[nn])
```
Note the sign convention: dep ≥ 0 and z = −dep ≤ 0.


In the solver, N(z) values are defined at cell boundaries:
z = 0, -dz[1], -(dz[1]+dz[2]), . . .
But, `conf.jl` maps the given data onto the gridpoints of the solver.
If the gridpoint falls outside the range of `dep[:]`,
the N value
at the closed defined point is used (that is a "constant extrapolation").

Therefore, if you want to use a constant N, just give something like
```
1
100.0  0.003
```
In the above, the dep doesn't matter at all.



## Appendix C: Installing and using the Julia interpreter
Tutorials are available

* [https://book.jinguo-group.science/stable/chap2/julia-setup/](https://book.jinguo-group.science/stable/chap2/julia-setup/)
* [https://github.com/JuliaLang/juliaup](https://github.com/JuliaLang/juliaup)

The below is an outline. If problems strike,
please refer to the above tutorials.

### Installing

Install [`juliaup`](https://github.com/JuliaLang/juliaup):
```shell
curl -fsSL https://install.julialang.org | sh
```
Perhaps you need to modify your search `PATH` to find the interpreter,
which is named `julia`.

To update the interpreter
```shell
juliaup update
```

### Running

Launch the REPL (interactive interpreter):
```shell
$ julia
. . . Julia's logo etc. . . .
julia>
```

To run your script within the REPL:
```julia
julia> include("yourscript.jl")
```
To run it from the shell command line:
```shell
$ julia "yourscript.jl"
```
As usual for scripting languages, you can turn your scripts into a command
by starting the script with
```shell
#!/usr/bin/env julia
```
and give the script file the exec permission by
```shell
$ chmod +x yourscript.jl
```

### Installing and updating packages
The julia interpreter includes its own package manager.
After launching the REPL, press the `]` key:
```julia
julia> ]
(@v1.11) pkg>
```
where `@v1.11` indicates the current version of the interpreter.

To exit the package manager, type `^C` (Ctrl+C):
```julia
(@v1.11) pkg> [Ctrl+C]
julia>
```

To install a package
```julia
(@v1.11) pkg> add PackageYouWant
```
To see the installed packages
```julia
(@v1.11) pkg> status
```
To update all the packages
```julia
(@v1.11) pkg> update
```

## Appendix D: Notes on the structure of `ctwmodes.f90`
The core of `ctwmodes.f90` is where the matrices
A and B are defined.
This part resides in a module within the source file.

After version 006b, we tested a new discretization scheme.
We anticipated that we would perhaps provide the two schemes
in such a way the user can switch between them.
To make this switching easy,
we extracted the matrix-definition part from 006b
into an external module and modified it for the new scheme.

That process resulted in `ctwmodes-006c.f90`, which handled only logistics,
plus the external module for the new scheme.

After comparing this with results from the Brink & Champan code,
we decided the new scheme is much better and there is no point in providing
the original scheme.

When we created version 007 to be published as 0.7.0, we included the
module in the source file of `ctwmodes.f90`
instead of pasting the contents of the module into the main program.
If we ever need to extend our program, we will again separate
the source file into the main program and the module
and implement new schemes in the latter.
