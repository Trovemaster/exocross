Cooling function
================

The cooling function is computed as emissivity (erg/molecule/sterradian) on a grid of temperatures (K). Here is an example of the input: 

::

    cooling
     Ntemps     500
     tempmax  5000.0
    end

    output      ScH_adj_16665     

    Transitions   ScH_adj_16665.trans
    States        ScH_adj_16665.states
 


Keywords: 
^^^^^^^^^

* `Ntemp` or `Npoints` is the number of temperature steps 
* `tempmax`, `Maxtemp` or `Max-Temperature`  is the maximal temperature in K (minimal T = 1K ) 

