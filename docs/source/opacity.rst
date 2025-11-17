Opacities
=========

Here we list features and constructs designed specifically to facilitate opacity produciton: computation of cross sections on a large grid of temperatures and pressures. 


The method is described in:
K. L. Chubb, M. Rocchetto, S. N. Yurchenko, M. Min, I. Waldmann, J. K. Barstow, P. Molliére, A. F. Al-Refaie, M. Phillips, and J. Tennyson. "The ExoMolOP Database: Cross-sections and k-tables for Molecules of Interest in High-Temperature Exoplanet Atmospheres". A&A, 646:A21, 2020. 

`See the paper here`_.

.. _See the paper here: http://dx.doi.org/10.1051/0004-6361/202038350


Computing cross sections on a grid of temperatures (array job)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For opacity productions when cross sections are computed for a large set of temperatures, it can be more efficient to process cross sections for all temperatures simultaneously, similarly how the partition and cooling functions are computed. To this end, the following ``ARRAY`` block structure can be used. The saving is due to a reduced IO: the line list files need to be read only ones and used for all the temperatures.

::

     Range 0 5900

     Npoints 5901

     array
      n 10
      tmin 100
      tmax 1000
     end

     absorption
     voigt

     mass 26

     species
        air   gamma 0.007  n 0.65 t0 296.0
     end

     pressure 1

     output 12C-14N__5000K_test
     NCache 500
     NPROCS 36

     States 12C-14N__MoLLIST.states
     Transitions 12C-14N__MoLLIST_10lines.trans

Here


* ``N``  is the number of temperature steps (e.g. 10)
* ``Tmin``  is the minima temperature in K
* ``Tmax``  is the maximal temperature in K

All cross sections are printed out into a single ``output`` file as different columns following the frequency column, e.g.

::

     ....
     5.88500000E+003   7.97795033E-036  1.72172105E-030  1.50934418E-028  1.84056995E-027  8.75964032E-027  2.47807332E-026  5.14418985E-026
     5.88700000E+003   1.17182482E-035  1.75296505E-030  1.32946282E-028  1.57424065E-027  7.50135150E-027  2.14298647E-026  4.51073365E-026
     5.88800000E+003   2.12842862E-035  2.98897895E-030  1.56465959E-028  1.28410634E-027  4.90675516E-027  1.23246563E-026  2.40312414E-026
     5.88900000E+003   6.41335689E-035  9.31393311E-030  4.36473210E-028  2.97204673E-027  9.52430574E-027  2.08438541E-026  3.65880214E-026
     ....

The main (standard) output will also contain values of the total absorption for each temperature.

The partition function will be recomputed based on the .states file provided.

.. warning:: This feature works currently only for the basic case of ``Absorption`` and ``Voigt`` with no ``Filters``, ``non-LTE`` etc.


Two-step production of cross sections  using super-lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For very large linel lists, it is more efficient (but less accurate) to (i) compute temperature-dependent super-lines and then use them to produce  temperature/pressure-dependent cross-sections. This can be done using the line profile ``Voigt-super``. Everything else is as for the normal case of the ``Voigt`` profile. This procedure can only use single line-broadening parameters. By defaults, the cross sections are computed on the same grid as the super-lines, i.e. the grid spacings are the same. Here is a basic example:

::

    Temperature 2000
    Range 2000 2100

    Npoints 1001

    pffile 12C-14N__KTPSYT_500.pf

    absorption
    voigt-super

    mass 26

    species
         air   gamma 0.007  n 0.65 t0 296.0 ratio 1.00 delta 0.000000
    end

    pressure 1

    output 12C-14N__2000K_Voigt


    States 12C-14N__MoLLIST.states

    Transitions 12C-14N__MoLLIST.trans



It is possible to use the resolving-power grid for super-lines and then mapped it onto an equidistant map of grid of ``Npoints``, for example

    Temperature 2000
    Range 2000 2100

    R  100000
    Npoints 10001

    pffile 12C-14N__KTPSYT_500.pf

    absorption
    voigt-super
    .....


Here, the super-lines are computed on the :math:`R=100000` grid and then transformed to :math:`N_{\rm points} = 10001` grid.


.. warning:: This feature works only currently only for the basic case of ``Absorption`` and ``Voigt`` with no ``Filters``, ``non-LTE`` etc.


Super-lines for an array of temperatures.
-----------------------------------------

The ``Voigt-super`` can be combined with the ``Array`` structure to be processed on a grid of temperatures by adding an ``Array`` section, e.g.


::

    Temperature 2000
    Range 2000 2100

    Npoints 1001

    pffile 12C-14N__KTPSYT_500.pf

    absorption
    voigt-super

    array
     N  10
     tmin 100
     tmax 1000
    end

    mass 26

    species
         air   gamma 0.007  n 0.65 t0 296.0 ratio 1.00 delta 0.000000
    end

    pressure 1

    output 12C-14N__2000K_Voigt

    States 12C-14N__MoLLIST.states
    Transitions 12C-14N__MoLLIST.trans



Computing cross sections on a list of pressures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: This feature is currently only available for ``Absorption`` and ``Voigt-super``, with or without ``Array`` and with no ``Filters``, ``non-LTE`` etc.

If the goal is to compute opacities on a list of temperatures and pressures and in order to take the full advantage of the two-step procedure (super-lines+Voigt)  ``Voigt-super``, it can be processed for a list of pressures as part of a single cross sections job. To this end, a block structure ``pressure-list`` is introduced as given by

::

   pressure-list
    1.00000000e-05
    2.15443469e-05
    4.64158883e-05
    1.00000000e-04
    2.15443469e-04
    4.64158883e-04
    1.00000000e-03
   end

The result is a set of output files with cross sections written for each of the pressure values given. The filename is appended with a suffix containing the  associated value of pressure, e.g. ``12C-14N__2000K_Voigt__1.00000000e-05.xsec``.

The ``pressure-list`` structure can be combined with the ``array`` of temperatures block, e.g.

::

   mem 16.0 gb

   Temperature 2000
   Range 2010 2070

   Npoints 90001

   R 100000
   pffile 12C-14N__KTPSYT_500.pf

   absorption
   voigt-super

   array
   n  2
   tmin 1000
   tmax 2000
   end

   pressure-list
   1.00000000e-05
   2.15443469e-05
   4.64158883e-05
   1.00000000e-04
   2.15443469e-04
   4.64158883e-04
   1.00000000e-03
   end

   mass 26

   species
       air   gamma 0.007  n 0.65 t0 296.0 ratio 1.00 delta 0.000000
   end

   pressure 1

   output 12C-14N__2000K_Voigt

   Ncache 500

   NPROCS 36

   verbose 4

   States 12C-14N__MoLLIST.states
   Transitions 12C-14N__MoLLIST_line.trans


