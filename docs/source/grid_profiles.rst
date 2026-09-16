.. _grid_profiles:

User-supplied GRID profiles
==========================

A ``PROFILES`` block selects a line profile supplied as a table of wavenumber
offsets and normalized densities. Each entry specifies a temperature in K,
a band-symmetry label, and a filename. Band symmetry is the direct product of
the upper- and lower-state irreps::

    Temperature 1000
    Range 0 10000
    Npoints 10001
    absorption
    Symmetry C2v
    Nirreps = 4

    QN
      IRREP 9
    END

    PROFILES
      1000 A1 profile_H2S_A1_T1000K.prof
      1000 B2 profile_H2S_B2_T1000K.prof
    END

    cutoff 1000
    output H2S_grid_T1000
    States i-H2S_p24_linearised_0_0.states
    Transitions i-H2S_p24_linearised_0_0.trans

The top-level ``SYMMETRY`` selects the group. Inside QN, ``IRREP`` and
``SYMMETRY`` are aliases for the state-file column. That column is read for
both states. The column number counts from the state ID in column 1 and must
be at least 5. The former UPPER/LOWER selection is not used. Labels can contain
up to 20 characters, matching the existing QN storage.

Supported groups are C2v (alias C2v(M), irreps A1, A2, B1, B2), C3v
(alias C3v(M), irreps A1, A2, E), and Cs (alias Cs(M), irreps A', A").
Cs also accepts two apostrophes, A'', for double prime. Group and irrep names
are case-insensitive. ``Nirreps`` is inferred as 4, 3 or 2 respectively; the
optional ``Nirreps = 4`` or ``Nirreps 4`` checks the count. The group and count
may appear before or after PROFILES and QN.

The internal profile name is ``GRID``, but no standalone GRID keyword should
be added: that keyword is already an alias for the existing multiple-output-grid
block. A PROFILES block selects this profile automatically. HWHM is unused.

Temperature arrays
------------------

For ``TEMPERATURE-LIST`` or ``TEMPERATURE-ARRAY``, provide one profile for every
requested temperature and every profile label, for example::

    temperature-list
      1000
      2000
    end

    PROFILES
      1000 A1 profile_H2S_A1_T1000K.prof
      2000 A1 profile_H2S_A1_T2000K.prof

      1000 B2 profile_H2S_B2_T1000K.prof
      2000 B2 profile_H2S_B2_T2000K.prof
    END

Profile entries may occur in any order. Output columns follow the temperature
list. Duplicate or missing temperature/label pairs and unrequested temperatures
are errors. Matching uses a relative temperature tolerance of 1e-10; no
temperature interpolation or extrapolation is performed. An unknown state irrep
or a missing product-component profile is an error if used by a retained
transition. Only the band irreps actually needed require profiles; a state
irrep need not have a profile under its own label.

Irrep products and degenerate bands
-----------------------------------

The product table for C2v is:

==== ==== ==== ==== ====
x    A1   A2   B1   B2
==== ==== ==== ==== ====
A1   A1   A2   B1   B2
A2   A2   A1   B2   B1
B1   B1   B2   A1   A2
B2   B2   B1   A2   A1
==== ==== ==== ==== ====

For C3v, A1 is the identity, A2 x A2 = A1, A2 x E = E, and
E x E = A1 + A2 + E. For Cs, A' is the identity and A" x A" = A'.

For the reducible E x E product, the initial profile prescription is the equal
average ``(f_A1 + f_A2 + f_E)/3``. Each component is interpolated on its own
grid. All three component profiles must be supplied at every requested
temperature. Their weights sum to one, preserving the original line intensity;
the E profile is not given an additional degeneracy factor.

This average is a provisional shape model. Group theory supplies the product
decomposition, not the relative transition strengths. No additional dipole
selection rules or physical branch weights are inferred.

The new ``symmetry.f90`` module follows TROVE's SymmetryT structure. It stores
``sym%Nirreps``, ``sym%label``, ``sym%degen`` and integer product multiplicities
``sym%product(gamma,upper,lower)``. ``gamma_lookup(upper,lower)`` returns the
multiplicity vector; for C3v E x E it returns [1,1,1] in A1,A2,E order.
The profile averaging policy is separate, in ``configure_band_weights`` in
``grid_profiles.f90``, so a future physical weighting can replace it without
changing the group algebra.

Profile files and sampling
--------------------------

Files have two whitespace-separated columns: offset in cm-1 and profile density
per unit wavenumber. Offsets must strictly increase and span zero; values must
be finite and nonnegative. At least two points are required. All files must
have the same point count, but their grids may differ and may be nonuniform.
Asymmetric profiles are allowed; they are not automatically recentered.
Blank lines and parenthesised comments are allowed. Paths are relative to the
working directory and should be quoted when they contain spaces.

The trapezoidal integral must be positive and within 1e-3 of one. The original
integral is printed, and accepted profiles are divided by that integral.
Profiles are then evaluated by linear interpolation at output wavenumber minus
line centre. No extrapolation is made beyond the supplied support.

Uniform grids use direct arithmetic indexing with spacing determined once
when each profile is loaded. Nonuniform profiles use binary search when a new
interpolation interval cannot be reached by advancing one neighbouring interval.
Both paths interpolate using the actual tabulated wavenumbers and values.

An explicit CUTOFF is in cm-1. If omitted, the tables' full support is used.
Lines centred outside RANGE are included when their wings can overlap it.
Truncation by CUTOFF or RANGE does not cause renormalization. This is point
sampling: an output grid too coarse to resolve a profile need not preserve its
integral exactly.

Supported calculations
----------------------

The initial implementation supports LTE absorption using ExoMol states and
transitions, one uniform output grid in cm-1, scalar or array temperatures,
and constant intensity thresholds. The existing restriction on combining
temperature arrays and filters remains. Pressure lists, super-lines, additional
analytic broadening, non-LTE, emission, HITRAN/SPECTRA input, and multiple or
resolving-power output grids are not supported.
