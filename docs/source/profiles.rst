Line profiles 
=============

Doppler 
^^^^^^^

This is the averaged Doppler line profile which conserves the area at each grid point for any resolution used. 
It requires a mass of the molecule (Dalton) defined.  


Example 
::

     temperature 2000.0 
     absorption
     Doppler
     mass 20.0 



Doppler Sampling
^^^^^^^^^^^^^^^^
::

     temperature 2000.0 
     absorption
     Doppler Sampling
     mass 20.0 




`Doppler Sampling` is used for the simple sampling method. 



Gaussian
^^^^^^^^
This is the averaged Gaussian line profile which conserves the area at each grid point for any resolution used. 
It requires an `HWHM` value to be defined.  

Example 
::

     temperature 2000.0 
     absorption
     gaussian
     hwhm .321


Gaussian Sampling
^^^^^^^^^^^^^^^^^

`Gaussian Sampling` is used for the simple sampling method. 


Example 
::

     temperature 2000.0 
     absorption
     gaussian Sampling
     hwhm .321



Voigt
^^^^^

Example 
::

     temperature 2000.0 
     absorption
     voigt
     hwhm .321
     mass 16.0313
     offset 25 (cm-1)



Lorentzian 
^^^^^^^^^^

Example 
::


     temperature 2000.0 
     absorption
     lorentzian
     hwhm .321
     offset 25 (cm-1)


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
     offset 25 (cm-1)


Pseudo-Liu 
^^^^^^^^^^

(Liu_Lin_JOptSocAmB_2001)


Example 
::

     absorption
     pseudo-Liu
     hwhm .321
     mass 16.0313
     offset 25 (cm-1)

Pseudo-Rocco
^^^^^^^^^^^^

(Rocco_Cruzado_ActaPhysPol_2012)

Example 
::

     absorption
     pseudo-Rocco
     hwhm .321
     mass 16.0313
     offset 25 (cm-1)


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


A :math:`J`-dependent set of broadening parameters can be provided in an external file, e.g.
:: 
     
     mass 16.0
     pressure 0.5 
     Temperature 1300.0 

     species 
       H2  gamma 0.0207 n 0.44 t0 298.0 file 1H2-16O__H2.broad model JJ ratio 0.84
       He  gamma 0.043  n 0.02 t0 298.0 file  1H2-16O__He.broad model JJ ratio 0.16
     end


where `file` is the filename with parameters and JJ (alias a1) is the name of  the model. Two models are available: 
J (or a0) and JJ (or a1), which stand for the broadening dependent on the lower only and the lower/upper Js. 

The broadening file has the following structure 
::

     0.0145 0.500       0       1
     0.0156 0.417       1       2
     0.0164 0.350       2       3


where the first two columns are Voigt's gamma and n, and the last two are J" and J' (i.e. in the opposite to the conventional order). The values `gamma` and `n` in the `species` section are the default values in case of missing Js in the broadening file. 



Voigt-Quad 
^^^^^^^^^^

`Voigt-Quad` is the Voigt obtained using the Guass-Hermite quadrature integrations. An analytical integration of the Lorentzian is used for the average contribution for each bin. The effect of the line truncation with offset parameter is folded back into the main part using the analytical expression. The line guarantees the area to conserve. 


Example 
::

     
     Temperature   500  (K)
     pressure 10. (bar)
     absorption
     Voigt-Quad
     mass 16.0313
     offset 25 (cm-1)
     nquad   20   (N quadrature points)
     
     Species
       H2  gamma 0.05 n 0.4 t0 298.0 ratio 0.9
       He  gamma 0.04 n 1.0 t0 298.0 ratio 0.1
     end
     
