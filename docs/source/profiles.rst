Line profiles
=============

Doppler
^^^^^^^

This is a binned Doppler line profile which conserves the area at each grid point for any resolution used.
It requires a mass of the molecule (Dalton) defined.


Example
::

     temperature 2000.0
     absorption
     Doppler
     mass 20.0



Doppler Sampling or Doppl0
^^^^^^^^^^^^^^^^^^^^^^^^^^
::

     temperature 2000.0
     absorption
     Doppler Sampling
     mass 20.0



`Doppler Sampling` is used for the simple sampling method.



Gaussian
^^^^^^^^
This is a binned Gaussian line profile which conserves the area at each grid point for any resolution used.
It requires an `HWHM` value to be defined.

Example
::

     temperature 2000.0
     absorption
     gaussian
     hwhm .321


Gaussian Sampling or Gauss0
^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Gaussian Sampling` is used for the simple sampling method.


Example
::

     temperature 2000.0
     absorption
     gaussian Sampling
     hwhm .321



Voigt (Sampling)
^^^^^^^^^^^^^^^^

Example
::

     temperature 2000.0
     absorption
     voigt
     hwhm .321
     mass 16.0313
     cutoff 25 (cm-1)



Lorentzian
^^^^^^^^^^

This profile uses a sampling method, where the cross sections at a given wavenumber points represents an aveerage over the wavenumber bin.

Example
::


     temperature 2000.0
     absorption
     lorentzian
     hwhm .321
     cutoff 25 (cm-1)



A binned version can be invoked using the ``BINNING`` keyword or with the ``Lorentz0`` (``LORENTZIAN0'``)  line profile, e.g.
::

     temperature 2000.0
     absorption
     lorentz0
     hwhm .321
     cutoff 25 (cm-1)


or
::

     temperature 2000.0
     absorption
     lorentzian binning
     hwhm .321
     cutoff 25 (cm-1)



For the binned profile method, the cross sections at a given wavenumber point represents an aveerage over the wavenumber bin.


Error cross sections
--------------------


To invoke the error cross section calculations, a free floating ``ERROR`` keyword is used.


Here, ExoCross uses the energies uncertainties to define uncertainties of the cross-sections in a form of absorption (emission) error cross-sections for different line profiles as given by

.. math::

   \left(\Delta \sigma(\tilde{\nu})\right)^2 = \sum_{i,j} \left(\frac{\partial \sigma(\tilde{\nu})}{\partial \tilde{\nu}_{ij}}\right)^2 \left[(\Delta \tilde{E}_i)^2 + (\Delta \tilde{E}_j)^2  \right],


where :math:`\Delta \tilde{E}_i` and :math:`\Delta \tilde{E}_j` are the uncertainties of the upper and lower states. Assuming a given line-profile :math:`f(\tilde\nu)`, the derivative wrt the energy is given by

.. math::

    \frac{\partial \sigma(\tilde{\nu})}{\partial \tilde{\nu}_{ij}} = I_{if} \frac{\partial f(\tilde{\nu})}{\partial \tilde{\nu}_{ij}}.


Here is an input example:
::

     elorentz
     error
     hwhm 0.2


It is important that the energy uncertainties are provided in the column 5 of the States file. The error cross sections can be combined with the uncertainty filters: ::

   filter
     unc 0.01
   end




Error cross sections for the Lorentzian profile (Elorentz)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


For the Lorentzian line profile centred at :math:`\tilde{\nu}_{ij}` with HWHM :math:`\gamma` given by

.. math::

  f(\tilde\nu,\tilde\nu_{ij},\gamma)_{\rm Lo} = \frac{\gamma}{\pi} \frac{1}{(\tilde{\nu}-\tilde{\nu}_{ij})^2+\gamma^2}


the corresponding derivative wrt :math:`\tilde{\nu}_{ij}` is given by

.. math::

    \frac{\partial f(\tilde\nu)_{\rm Lo}}{\partial \tilde\nu_{ij}} = \frac{\gamma}{\pi} \frac{2(\tilde{\nu}_{ij}-\tilde{\nu})}{\left[(\tilde{\nu}-\tilde{\nu}_{ij})^2+\gamma^2\right]^2}.



Cross sections with energy uncertainties as line broadening (Voigt-unc)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we assume that the energy uncertainties distributed according with the normal distribution with :math:`\sigma` as the uncertainty of the corresponding line position :math:`\tilde\nu_{ij}`:

.. math::

    \sigma_{ij}^{\rm unc} = \sqrt{\Delta \tilde{E}_{i}^2+ \Delta \tilde{E}_{j}^2},

we can obtain the cross sections with the energy uncertainties included by convolving a Gaussian with :math:`\alpha_{ij}^{\rm G} = \sqrt{2\ln(2)} \sigma_{ij}`  with the Lorentzian (Voigt) of $\gamma_{\rm V}$. The resulted line profile us Voigt. The Doppler broadening can be added to the uncertainty broadening as follows

.. math::

    \alpha_{ij}^{\rm G} =  \sqrt{2\ln(2) (\sigma_{ij}^{\rm unc})^2 + (\alpha_{ij}^{\rm Doppl})^2}.

This feature is implemented as the ``Voigt-unc`` type and can be used as in the following example:
::

     temperature 1000 Range  0 20000

     Npoints 2000000

     mass 64

     absorption
     Voigt-unc

     pressure  1

     species
      H2  gamma 0.0468  n 0.500 t0 296.0  ratio 0.85
      He  gamma 0.0468  n 0.5   t0 296.0  ratio 0.15
     end


     Output TiO_Voigt-unc_T1000K_P1atm

     States 48Ti-16O__Toto.states

     Transitions 48Ti-16O__Toto.trans







Pre-dissociative line profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In case of pre-dissociative effects, the lines van be broadened beyond the collisional or Doppler broadening. For these cases, ExoCross can use the lifetimes :math:`\tau`  from the States file (usually column 6 after the uncertainty column) to estimate the pre-dissociative line broadening (HWHM) via


.. math::

    \gamma_{\rm prediss}  = \frac{1}{2\pi c \tau} \frac{1}{2}

where :math:`c` is the speed of light in cm/s, :math:`\tau` is the lifetime in s and the factor :math:`1/2` is to convert to the HWHM :math:`\gamma`. ExoCross will apply the largest of the two line broadening values, :math:`\gamma_{\rm prediss}` and :math:`\gamma_{\rm collis}`.  This option is activated via a free floating keyword ``PREDISSOCIATION``. By default, the lifetimes column is assumed to be column 6. Otherwise it is important to specify the number of the lifetime column as part of the ``QN`` section using the keyword ``lifetime``
::

     predissociation

     QN
      lifetime 5
     END


.. note:: The lifetime column specification can be combined with the ``non-LTE`` section, since ``QN`` and ``non-LTE`` are essentially aliases of each other.




Pseudo-Voigt profiles (less accurate)
-------------------------------------


Pseudo-Voigt
^^^^^^^^^^^^

See  wiki_

.. _wiki: https://en.wikipedia.org/wiki/Voigt_profile](https://en.wikipedia.org/wiki/Voigt_profile

The Pseudo-Voigt Profile (or Pseudo-Voigt Function) is an approximation of the Voigt Profile V(x), using a linear combination of a Gaussian curve G(x) and a Lorentzian curve L(x) instead of their convolution.

The Pseudo-Voigt Function is often used for calculations of experimental Spectral line shape profiles.

The mathematical definition of the normalized Pseudo-Voigt profile is given by:


:math:`V_p(x)= \eta \cdot L(x) + (1-\eta) \cdot G(x)`


with  :math:`0 < \eta < 1`

There is several possible choices for the  \eta  parameter. A simple formula, accurate to 1%, is:
:math:`\eta = 1.36603` (:math:`f_L/f`) - 0.47719 :math:`(f_L/f)^2 + 0.11116(f_L/f)^3`
where

:math:`f = [f_G^5 + 2.69269 f_G^4 f_L + 2.42843 f_G^3 f_L^2 + 4.47163 f_G^2 f_L^3 + 0.07842 f_G f_L^4 + f_L^5]^{1/5}`

Example
::

     absorption
     pseudo
     hwhm .321
     mass 16.0313
     cutoff 25 (cm-1)


Pseudo-Liu
^^^^^^^^^^

(Liu_Lin_JOptSocAmB_2001)


Example
::

     absorption
     pseudo-Liu
     hwhm .321
     mass 16.0313
     cutoff 25 (cm-1)

Pseudo-Rocco
^^^^^^^^^^^^

(Rocco_Cruzado_ActaPhysPol_2012)

Example
::

     absorption
     pseudo-Rocco
     hwhm .321
     mass 16.0313
     cutoff 25 (cm-1)


Voigt-parameters
^^^^^^^^^^^^^^^^

`Species` or  `Broadener` starts a section to define the Voigt-type broadening parameters

     :math:`\gamma(Voigt) = \sum_i \gamma_i (T^0_i/T)^n P/P^0_i {\rm ratio}_i`

The keywords are:

`gamma` or `gamma0` is the reference HWHM (cm-1), `n` is the exponent n_i, `T0` is the reference T (K),usually 298, `P0` is the reference pressure in bar, usually 1, `ratio` is the mixing ratio of the species (unitless), for example the solar mixing ratio of H2 and He is 0.9 and 0.1.

The name of the species should be the first thing on the line.

The `pressure` value in bar must be specified (otherwise P=1 bar is assumed).

The effective molar `mass` of the molecule/atom mass  be specified (1.0 is the default).

Example
::

     mass 16.0
     pressure 0.5
     Temperature 1300.0
     Species
       H2  gamma 0.05 n 0.4 t0 298.0 ratio 0.9
       He  gamma 0.04 n 1.0 t0 298.0 ratio 0.1
     end




Voigt-Quad
^^^^^^^^^^

`Voigt-Quad` is the Voigt obtained using the Guass-Hermite quadrature integrations. An analytical integration of the Lorentzian is used for the average contribution for each bin. The effect of the line truncation with cutoff parameter is folded back into the main part using the analytical expression. The line guarantees the area to conserve.


Example
::


     Temperature   500  (K)
     pressure 10. (bar)
     absorption
     Voigt-Quad
     mass 16.0313
     cutoff 25 (cm-1)
     nquad   20   (N quadrature points)

     Species
       H2  gamma 0.05 n 0.4 t0 298.0 ratio 0.9
       He  gamma 0.04 n 1.0 t0 298.0 ratio 0.1
     end


ExoMol diet and broadening recipes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ExoCross is equipped to work with the ExoMol diet, which provides line broadening parameters in .broad files. The ExoMol diet allows for the line broadening parameters :math:`\gamma` and :math:`n` to depend, at least in principle, on any quantum numbers used in .states files. In practice, only simplest quantum number cases are currently in use, mainly for the rotational quantum number :math:`J`. Different quantum number cases ('recipes') are labeled with two character strings including ``a0`` (:math:`\gamma` and :math:`n` depend on :math:`J''` only), ``a1`` (:math:`J''` and :math:`J'`), ``m0`` and ``m1`` (rotational index :math:`|m|` and :math:`m`, respectively), ``v0`` (:math:`v'`), ``v1`` (:math:`J''` and :math:`v'`), ``k1`` (:math:`J''` and :math:`k'`). Currently, the vibrational :math:`v` or rotational quantum :math:`k` numbers can only be specified for the upper states.


A :math:`J''`-dependent set of broadening parameters can be thus provided in an external file using the ExoMol Diet structure, recipe ``a0``, e.g.
::

     mass 16.0
     pressure 0.5
     Temperature 1300.0

     species
       H2  gamma 0.0207 n 0.44 t0 298.0 file 1H2-16O__H2.broad  ratio 0.84
       He  gamma 0.043  n 0.02 t0 298.0 file  1H2-16O__He.broad ratio 0.16
     end


where `file` is the filename with the parameters :math:`\gamma` and :math:`n`. An example of the ``a0`` recipe is as follows
::

    a0   0.0145 0.500       0
    a0   0.0156 0.417       1
    a0   0.0164 0.350       2


where the 1st column describes the Diet model, the following two columns are the Voigt's :math:`\gamma` and :math:`n`, and the last one is  :math:`J''` . The values ``gamma`` and ``n`` in the ``species`` section are the default values in case of missing :math:`J` in the broadening file. For :math:`J>J_{\rm max}`, the values :math:`\gamma` and :math:`n` the values of :math:`J=J_{\rm max}` are assumed.


The  ``a1`` recipe allows providing both the upper and lower state :math:`J`, expecting :math`|J'-J''|\le 2`.
An example of an ``a1`` recipe broadening file has is as follows:
::

    a1   0.0145 0.500       0       1
    a1   0.0156 0.417       1       2
    a1   0.0164 0.350       2       3


Here the last two columns list :math:`J''` and :math`J'` (i.e. opposite to the conventional order).

Another  example  of .broad for ``m1``:
::

    m1   0.0156 0.417      -2
    m1   0.0164 0.350      -1
    m1   0.0145 0.500       0
    m1   0.0156 0.417       1
    m1   0.0164 0.350       2



'Particle in the box'
---------------------

This 'recipe' is for continuum states spectra calculations, i.e. for transitions between bound and unbound (usually upper) states. The corresponding line lists are represented by discretized transitions with their intensities to be re-distributed across the simulation wavenumber regions. This is typically done with a Gaussian line profile with a larger width (50-200 cm-1) to cover gaps between these discrete lines. In the case of a very large spectroscopic range, the distance between these discrete continuum lines rapidly increases. For example, for the case of the particle-in-a-box solution, the energy separation increases with the excitation number linearly: 

.. math::


    \Delta \tilde{E}_n^{\rm box} = \tilde{E}_{n+1}^{\rm box} - \tilde{E}_n^{\rm box} = \frac{h (2n+1)}{8 \mu L^2 c}.

where :math:`\mu` is the reduced mass, :math:`L` is the box size) and :math:`n\ge 1` is the state counting number. 
Due to the linear dependence on :math:`n`,  there is no a single optimal value of :math:`\alpha_{\rm G}` (Gaussian HWHM)  for the entire region. In the ``BOX`` recipe, we define  :math:`alpha_{\rm G}` to be :math:`\Delta \tilde{E}_n^{\rm box}`. Here is an example of the ExoCross input:
:: 

   Temperature  5000
   Range    45000 85000 
   
   npoints  40000
   
   QN
     K 8
   end
   
   absorption
   gaussian
   
   offset 5000
   
   species
        box  gamma 1.0 mass 0.937  Lbox 20  model box
   end
   
   output NH_box_5000K_L20
   
   States NH.states
   
   Transitions NH.trans
   

Here  ``model box`` defines the recipe type ``box``, ``Lbox`` is the size of the box and ``mass`` specifying the reduced mass. 




