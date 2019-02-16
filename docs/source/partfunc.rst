Partition functions
===================

Example:
::

    partfunc
     Ntemps     500
     tempmax  5000.0
     Moment 0
    end
    
    output    ScH_adj_16665 
    States    ScH_adj_16665.states
 

Specific heat
^^^^^^^^^^^^^ 
::


    partfunc
     Ntemps     500
     tempmax  5000.0
     Cp
    end



Keywords: 
^^^^^^^^^

* `Ntemp` or `Npoints` is the number of temperature steps 
* `tempmax`, `Maxtemp` or `Max-Temperature`  is the maximal temperature in K (minimal T = 0K ) 
* `Moment` (0,1,2, or 3) is the moment to compute, 0 is the partition function (default and can be omitted), 1 is the 1st moment, 2 is the 2nd moment and 3 is CP.
* `CP` is equivalent to Moment 3

