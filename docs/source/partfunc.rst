Partition functions
===================

Generation of partition funcitons
---------------------------------

A partition functions can be computed on a grid of temperatures using the  block ``partfunc``, which specifies the temperature grid via ``Tmax`` (maximal temperature) and ``N`` (the number of the temperatures).


Example:
::

    partfunc
     N     500
     TMax  5000.0
    end

    output    ScH_adj_16665
    States    ScH_adj_16665.states


Different moments (1,2,3) can be computed as follows

Example:
::

    partfunc
     N     500
     TMax  5000.0
     moment 1
    end

    output    ScH_adj_16665
    States    ScH_adj_16665.states


Alternatively, the temperature grid can be  given as  ``TEMPERATURE-LIST``, see the corresponding documentation. The structure of the input in this case is as follows: 


    partfunc
      list
    end
    
    temperature-list
     200
     296
     300
     1000
    end
    
where the keyword ``list`` is used to indicate the temperature will be read from  a ``temperature-list``.


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

* ``N``, ``Ntemp`` or ``Npoints`` is the number of temperature steps
* ``Tmax``, ``tempmax``, ``Maxtemp`` or ``Max-Temperature``  is the maximal temperature in K (minimal T = 0K )
* `Moment` (0,1,2, or 3) is the moment to compute, 0 is the partition function (default and can be omitted), 1 is the 1st moment, 2 is the 2nd moment and 3 is CP.
* `CP` is equivalent to Moment 3


Providing partition external function functions
-----------------------------------------------

An single partition function value can be defined using the keyword ``pf``:
::

    pf 3817.0


Alternatively, a file containing a partition function on a grid of temperatures (assuming the standard format, :math:`T\,\,\,\,\,Q(T)`) can be provided using the keyword ``pffile``, e.g.
::

    pffile  28Si-16O__SiOUVenIR.pf

.. note: If the temperature(s) requested is outside (larger than) the temperature range in the partition function file, the value will be computed using the energies provided in the .states file.



