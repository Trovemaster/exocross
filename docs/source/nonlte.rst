Non-LTE
=======

There currentoy two main methods to model non-LTE spectra with ExoCross: (i) 
using the two-temperature Treanor distribution and (ii) custom vibrational 
densities.  

Treanor distribution 
^^^^^^^^^^^^^^^^^^^^

The non-LTE spectra are modelled using the Treanor, non-Boltzmann distribution, which that the rotational and vibrational modes themselves are in LTE and the non-LTE population of a given state is taken as the product of the two Boltzmann distributions \citep{ExoCross,19PaLaxx}

:math:`F_{J,v,k}(T_{\rm vib},T_{\rm rot}) = e^{-c_2 \tilde{E}_{v}^{\rm vib}/T_{\rm vib}} e^{-c_2 \tilde{E}_{J,k}^{v,\rm rot}/T_{\rm rot}},`

where  :math:`c_2= hc / k_B` is the second radiation constant (cm K), :math:`\tilde{E}_i = E_i/h c` is the energy term value, and :math:`T` is the temperature in K.


An absorption line intensity :math:`I_{\rm fi}` (cm/molecule) is then given by

:math:`I({\rm f} \gets {\rm i}) = \frac{g_f^{\rm tot} A_{\rm fi}}{8 \pi c \tilde{\nu}_{\rm fi}^2}  \frac{F_{J,v,k}(T_{\rm vib},T_{\rm rot}) \left( 1-e^{-c_2\tilde{\nu}_{\rm fi}/T} \right)}{Q(T)},`

where :math:`A_{\rm fi}` is the Einstein-A coefficient (:math:`s^{-1}`), :math:`\tilde{\nu}_{\rm fi}` is the transition wavenumber, :math:`Q(T)` is the non-LTE partition function defined as a sum over states

:math:` Q(T) =\sum_{n}  g_n^{\rm tot} F_{J,v,k}(T_{\rm vib},T_{\rm rot}),`

:math:`g_n^{\rm tot}` is the total degeneracy given by 

:math:`g_n^{\rm tot} = g^{\rm ns}_n (2 J_n+1),`

:math:`J_n` is the corresponding total angular momentum, :math:`g^{\rm ns}_n` is the nuclear-spin statistical weight factor,  :math:`c_2= hc / k_B` is the second radiation constant (cm K), :math:`\tilde{E}_i = E_i/h c` is the energy term value, and :math:`T`  is the temperature in K.


For this approach ExoCross equires the vibrational and rotational temperatures as well as the rotational and vibational energies. The separate vibrational and rotational temperatures for the molecule and approximating the total energy as the sum of the respective rotational and vibrational (or vibronic) energies:

:math:`\tilde{E}_{v,J,k} = \tilde{E}_{v}^{\rm vib} + \tilde{E}_{J,k}^{v,\rm rot},`

where $v$ and $k$ are generic vibrational (vibronic) and rotational quantum numbers, respectively. The rotational contribution is  given by

:math:`\tilde{E}_{J,k}^{v,\rm rot} = \tilde{E}_{v,J,k} - \tilde{E}_{v}^{\rm vib}.`


An example of a non-LTE ExoCross input for a rotational tempeature of 296 K and vibrational temperature of 2000 K in this case is given by

::

    Temperature  296.0 vib 2000 
    Range 0.0  10000.0
    
    Npoints 10001

    NON-LTE
      JREF 1.5
      vib 6 7
    end

    absorption
    gauss
    hwhm 0.5 (cm-1)
    
    output OH_T296_2000K_non-LTE

    States      16O-1H__MoLLIST.states
    Transitions 16O-1H__MoLLIST.trans
    
    

Here `NON-LTE` has the same stucture and meaning as a section `QN` used to specify the columns with the (vibrational) quantum numbers. 
The keyword `Jref` specifies  the reference values of :math:`J` used to select the vibrational energies with the default value of 0. 
The keyword `vib` defines the range of columns in the states file containing the vibrational quantum numbers, which are also used 
fo building the reference set of vibrational energies. 



Custom vibrational densities 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vibrational densities :math:`N_{\rm vib}` can be also inputed directly into the calculations using the states file. In this case the non-LTE population density is 
given by 

:math:`F_{J,v,k}(T_{\rm vib},T_{\rm rot}) = N_{\rm v} e^{-c_2 \tilde{E}_{J,k}^{v,\rm rot}/T_{\rm rot}}.`

The non-LTE partition function :math:`Q(T)` is defined as 

:math:`Q(T) =\sum_{n}  g_n^{\rm tot} F_{J,v,k}(T_{\rm vib},T_{\rm rot}).`


Here is an example of a non-LTE ExoCross input for a rotational tempeature of 296 K and the normalised vibrational denisty 
listed in column 8:


::

    Temperature  296.0
    Range 0.0  10000.0
    
    Npoints 10001

    NON-LTE
      JREF 1.5
      denstity 8
      vib 6 7
    end

    absorption
    gauss
    hwhm 0.5 (cm-1)
    
    output OH_T296_2000K_non-LTE

    States      16O-1H__MoLLIST.states
    Transitions 16O-1H__MoLLIST.trans
    


     


