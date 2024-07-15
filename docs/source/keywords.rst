ExoCross keywords
=================

.. _keywords:

.. note:: The position of the keywords is not important. The input is not key-sensitive.

* GNS: Nuclear statistical weight `gns 2.0`


* Temperature:  Temperature (K):

::

    Temperature  3000.0


`Aliases`: Temp, Temperature

+ Temperature can be split into rotational temperature and vibrational temperature, which in input is given as

::

    Temperature  (Rotational) 700 Vib (vibrational) 2000 (K)


In this case it is also important to provide the addresses of the vibrational columns using the ``QN`` section.


+ A reference temperature (used e.g. with `SPECTRA`) can be specifed using the `REF` keyword as part of the `Temperature` line, e.g:

::

    Temperature 1000 Ref 296



* QN: The designation of the quantum numbers  (QN) to columns can be done using the QN section. In this section `vib` is used for the range of the vibrational QNs; `K` is for the column with the rotational "K"; `Sym` is for symmetry; `Nsym` is number of the symmetries; `Nmodes` stands for the number of vibrational modes. The column numbering starts after the J column (4th). The vibrational addressing is required for when vibrational temperature is used.
`density` can be used to give the non-LTE vibrational populations. ``Lifetime`` is to define the column with lifetimes, needed for the pre-dissociation line broadening calculations.

::

      QN
        Vib 2 7
        Nmodes 9
        K  7
      end


::

      QN
        Vib 2 7
        density 9
      end


* PF: input value of the partition function.

::

    PF  1000.0


+ A reference partition function (used e.g. with `SPECTRA`) can be specified using the `REF` keyword as part of the `PFe` line, e.g:

::

    PF 4000 Ref 100


* Range:  Wavenumber (cm-1) or wavelength (um or micron) range, i.e. from nu1 to nu2:

::


    Range 0.0 2000.000


Here is example for microns (wavelength), which however does nor work will all profiles:

::

     Range 100.0 200.000 um


`Aliases`: Window, Frequencies, Wavenumbers, Range

* `um` or `micron`: units of wavelength used in combination with Range

* `cm-1`: units of wavenumber used in combination with Range, can be omitted as default.

* `abundance`: is the molecular abundance, which is the factor used in intensities or cross-sections.

::

     abundance 0.97


* `Cut-off`: to use the HITRAN-cut-off scheme (see HITRAN 2012 paper). `nu_crit` is to give the corresponding parameter in the HITRAN cut-off scheme.

::

     cutoff HITRAN (S_crit) 1e-29  nu_crit 2000

* `Cutoff exp`: is used together with `trans` for partitioning into strong and weak parts for a reference temperature using the
  following dynamic threshhold: :math:`I_{\rm cut} = I_0 \exp(-\nu/\alpha)`.

Example
::

    cutoff Exp 1e-25  alpha  2000



* Npoints: Number of grid points (usually an odd number)

::

    Npoints 100001


> Aliases: Npoints, Number-of-points

* `Absorption`: type of spectra

::

    Absorption

* `Emission`: type of spectra`

::

    Emission


* `lifetime`: type of spectra, to compute life time of Ej as 1/sum_i A_ji; saved into a new .life file.

Only one of these three types is needed, it can appear outside any sections anywhere in the program. 'Absorption' is default.

* Partfunc or Partition-function: To compute the partition function and its contributions for a set of temperature; it also starts  a section defining the corresponding parameters

- Ntemps: Number of temperature steps

::

    Ntemps 500


- TempMax, maxtemps: Maximal temperature

::

    TempMax 5000.0

- Moment: Moment to compute (0=none,1=Q1,2=Q2,3=CP)

::

    moment 1

- CP: the same as Moment 3

::

    CP

* output: File name for crosssections, stick-spectra etc.

::

    output CH4_stick_T300K


* Cooling is a section name to compute the cooling function on a grid of temperatures. The parameters used in the Cooliing section are TempMax (maxtemps) and Ntemps, see Partfunc.

* `States`, `StatesFile`, `StateFile`, `States_file`, `States-File`: The name of the .states file `States CH4_linelist.states`

* Transitions, TransitionFiles: The name of the .trans file or the section name for the list of the .trans files.

::

    Transitions CH4_linelist.trans


or:
::

    Transitions
     CH4_100.trans
     CH4_200.trans
     CH4_300.trans
     CH4_400.trans
    end


In combination with the keyword `histogram-J` the transition filename can be followed by the J-value this file is associated with.


* Threshold: Intensity threshold. A line is skipped from line profile evaluation, or simply from the output if the corresponding absorption coefficient/emissivity is smaller than `Threshold`.

::

    Threshold 1e-28


* Enermax: Energy Threshold used to select states below some energy value

::

    enermax  20000.0


* Gaussian, Gauss, Gauss0, Doppler, Doppl, Doppl0, Box, Bin, Rect,  Sticks, Stick, Voigt, pseudo, pse-Rocco, pse-Liu, Voig-Quad, Lorentz, elorentz:  types of the line profile.


::

    gauss

    Doppler

    stick

    bin

    trans

* Sampling: used together with the line profile to indicate that a sampling (not averaging) version will be used. For example

::

    gaussian sampling


Currently this makes sense only in combination with Gaussian and Doppler.


* HWHM, HalfWidth: Half width at the half maximum, used for Gauss, Lorenz.


::

    `HWHM  0.1`


* Mass, Masses: Effective molecular mass (u0), used for Doppl.

::

    mass 16.0


* Verbose: verbose level, defines the amount of data to be printed out in the output `Verbose 3`

* `histogram` computing crosssections using intensities saved on (usually) an equidistant grid.


* `histogram-J` computing cross-sections using intensities saved on a grid with J-dependent histogram files.


* `NRAM` or `Ncash`: number of transitions to put into RAM; Alias is `LINES-TO-CASH`


* `Nprocs`: Number of OMP processes. Ideally should be the same as the number of omp-processors allocated, but would work with any number. Aliases are `OMP_NUM_PROCS` and `OMP_PROCS`.


* `mem`: maximal memory allocated for the job. It is used to estimate how many transitions can be put into RAM. Should be less than or equal too the memory of the system.




* ``cutoff`` (alias `line-cutoff`) is to define line-width cut-off in wavenumbers (cm-1). The default value is 25 cm-1.

::

    cutoff 25 (cm-1)


* `line-cutoff` (alias `cutoff`) is to define line-width cut-off in wavenumbers (cm-1). The default value is 25 cm-1.

::

    cutoff 25 (cm-1)



`cutoff` or `line-cutoff` (line-wings cut-off)  can be defined in terms of HWHM:

::

    cutoff 50 HWHM


`cutoff`s can allow using different values for different regions, see the multi-grid section
`grid`.



* `pressure` value in bar must be specified (otherwise P=1 bar is assumed) used for Voigt.


* `Species` or  `Broadener` starts a section to define the Voigt-type broadening parameters

:math:`\gamma = \sum_i \gamma_i (T^0_i/T)^n P/P^0_i r_i`


`gamma` or `gamma0` is the reference HWHM (cm-1) for Voigt used in `Species`

`n` is the exponent n_i for Voigt used in `Species`

`delta` is the pressure shift.

`T0` is the reference T (K), usually 298 for Voigt used in `Species`.

`P0` is the reference pressure in bar, usually 1 for Voigt used in ``Species``.

`ratio` is the mixing ratio of the species (unitless) for example the solar mixing ratio of H2 and He is 0.9 and 0.1 (`species`).

`file` is the name of the file with broadening parameters

`model` is the broadening model, e.g.  ``a0``, ``a1``, ``box`` etc. 

``nquad`` is the number of quadrature points used for `Voigt-Quad`.

``box`` is a ``model`` of the line broadening where the HWHM depends on the state number as in a particle in a box,  linearly. 


The name of the species should be the first thing on the line.

Example:
::


     mass 16.0
     pressure 0.5
     Temperature 1300.0
     Species
       H2  gamma 0.05 n 0.4 t0 298.0 ratio 0.9
       He  gamma 0.04 n 1.0 t0 298.0 ratio 0.1
     end



Example 2:

::

     Species
       H2  gamma 0.0207 n 0.44 t0 298.0 file 1H2-16O__H2__a1.broad model JJ ratio 0.84
       He  gamma 0.043 n 0.02 t0 298.0 file  1H2-16O__He__a1.broad model JJ ratio 0.16
     end


* `unc` is the keyword for specifying the uncertainty threshold. it should appear in the filer section.


* `filter` is the section name for specifying the filters

* `Phoenix` is the keyword for converting ExoMol line list to the Phoenix format. The `species` are expected to specify the Voigt parameters of the broadeners. `Phoenix` should appear anywhere in the main body of the input in the way as a line profile keyword.

* `HITRAN` is to use the HITRAN-format of the transition file or output.  Reading from hitran (.par) requires also the definition of the partition function `pf` and the isotopologue number `iso`. No .states is needed. To read from HITRAN use `HITRAN READ`

* `OXFORD` is to convert to the Oxford-format. The input structure is the same as used for writing in the HITRAN format.

* `TRANS` is used to partitioning the line list (.trans) files into weak and strong part defined for a reference temperature. Currently used in combination with
   `cutoff exp`.

* `WATT` or `WATT/STR/MOLECULE` can be used after the `emission` keyword to switch to the watt/str/molecule untis. E.g.


    emission watt


* `PHOTONS` or `PHOTONS/S` can be used after the `emission` keyword to switch to the photons/s untis. E.g.


    emission photons/s



* `SPECTRA` is to use the SPECTRA-format (IAO.ru, Tomsk) of the transition file.  This will also require the definition of the (i) reference temperature, (ii)  partition function for the target temperature, (iii)
partition function for the reference temperature, and (iv) the molecule/isotope pair  (`iso`).

Example:
::

    Temperature  500.0  ref 296.0
    Range 0.0  10000.0

    Npoints 10001

    absorption
    Voigt

    spectra
    iso 26 1
    pf 1000.0 ref 500

    output CH4_voigt_T500K

    Transitions  SpectraMol_CH4_296K.txt



* `iso` is to define the isotopologue for HITRAN or SPECTRA molecule, e.g. 26 1 for 12C2H2.

Example
::

    hitran
    iso 26 1
    pf 1000.0
    output C2H2_ab_g0.5
    Transitions  26_hit12.par


 `HITRAN` can also form a section where there HITRAN-error codes are specified, e.g.
::

     hitran
      write
      error-E  qn 4 ierr 4  vmax  10 ierr 3  vmax  20 ierr 2  vmax  30  ierr 1  vmax  40 ierr 0  vmax  100
      error-S  qn 4 ierr 5  vmax  10 ierr 4  vmax  20 ierr 3  vmax  100
      error-Air   ierr 4
      error-Self  ierr 4
      error-N      ierr 4
      error-delta  ierr 4
      error-unc
     end

Here ``error-unc`` is to define the energies/frequency error using the ExoMol uncertanties from column 5 in States file.

``error-E`` and ``error-S`` define uncertanties of the frequencies and intensities based on the quantum numbers.

``error-Air``, ``error-self``, ``error-N``, ``error-delta`` define  ucertanties of the corresponding line shape parameters.


* `grid` is to define a multi-grid with different resolutions in different sub-grids as follows

Example
::

    grid
      Range   0    100   Npoints 10000 cutoff 25.0
      Range 100   1000  Npoints 1000
      Range 1000 10000  Npoints 100
    end


Here `cutoff` can be substituted by `LINE-CUTOFF`.


+ `error-E` and `error-S` are used to specify the ranges for the quantum numbers for different Energy and intensity error-codes, respectively.

+ `qn` is the quantum number (the number of the QN-column after J) used for the error-specification.

+ `ierr` is followed by the error code, followed by

+ `vmax` (keyword) followed by the maximum value of qn this error code applies, which is followed by another error code.

+ `error-Air ierr` to give the error code for the HITRAN air-broadening (one constant value).

+ `error-self ierr` to give the error code for the HITRAN self-broadening (one constant value).

+ `error-n ierr` to give the error code for the HITRAN n-broadening (one constant value).

+ `error-delta ierr` to give the error code for the HITRAN line-shift (one constant value).



* `non-LTE`: A name of the section for non-LTE calculations containing the infomation on the columns with the vibrational quantum numbers, vibrational densites as well as the reference vaue of :math:`J`.

* `Jref`: a keyword specifying the reference :math:`J` value used to define reference vibrational energies taken as rovibrational (rovibronic) term values from the States file at this value of :math:`J` as part of the non-LTE calculations.

* `density`: a keyword specifying the column that contains the custom vibrational densities used in the non-LTE calculations.


* `Vald`: a single keyword trigering the conversion to the VALD format with log10(gf) factors using the astronomical convention for the statistical weight.
   To convert to the Astrophysical convention of the degeneracy factors, use `gf_factor`. For example, for water the convertion factor is 0.25.
   This will modify both gf and the g-factor.


::

    gf_factor 0.25



* `GF`: computer gf-factors. See the note about the gf_factor in `Vald` to convert to the astrophyical convention.

* `gf_factor`: a keyword to convert gf and g-factors to the Astrophysical convention of the degeneracy. For example for water it should 0.25 (see note for Vald).


* `predissociation`: a free floating single-valued keyword to invoke the pre-dissociative line broadening, where the lifetimes are used to estimate the pre-dissociative HWHM. For this option, the States file must contain the lifetimes and its column number must be specified in the ``QN`` (``nono-LTE``) section using the keyword ``lifetime``, unless it is the default column 6:

::

     predissociation

     QN
      lifetime 5
     END


.. note:: The lifetime column specification can be combined with the ``non-LTE`` section, since ``QN`` and ``non-LTE`` are essentially aliases of each other.



* Voigt-unc: Line profile used for to calculate cross sections broadening by the uncertainty of the line positions.

* elorents: Line profile to propagate the uncertainties of the line positions to the cross-sections used in combination with ``error`` 

* error: type of the cross sections describing the uncertainty of the cross sections as propagation of the uncertainty of the line positions. 




