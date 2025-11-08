Cooling function
================

The cooling function is computed as emissivity (erg/molecule/sterradian) on a grid of temperatures (K). Here is an example of the input: 

::

    cooling
     Tmax     5000.0
     N  5000
    end

    output      ScH_adj_16665     

    Transitions   ScH_adj_16665.trans
    States        ScH_adj_16665.states
 


Keywords: 
^^^^^^^^^

* ``N`` or ``Ntemps`` is the number of temperature steps 
* ``Tmax``, ``tempmax``, ``Maxtemp`` or ``Max-Temperature``  is the maximal temperature in K (minimal T = 0K ) 

