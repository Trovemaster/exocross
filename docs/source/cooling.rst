Cooling function
================

The cooling function is computed as emissivity (erg/molecule/sterradian or Watt/Str/Molecule) on a grid of temperatures (K). Here is an example of the input: 

::

    cooling
     Tmax     5000.0
     N  5000
    end

    output      ScH_adj_16665     

    Transitions   ScH_adj_16665.trans
    States        ScH_adj_16665.states
 

The ``colling`` functionality can be combined with the ``temperature-list`` functionality similar to the computation of the parition funciton with ``partfunc``. 
::
    
    cooling
     list
    end
    
   temperature-list
   200
   300
   1000
   end
    
The default units are erg/molecule/sterradian as the units of emissiom, which can be changed to Watt/Str/Molecule as follows:
::

    cooling
     Tmax     5000.0
     N  5000
     unit Watts 
    end

    
Keywords: 
^^^^^^^^^

* ``N`` or ``Ntemps`` is the number of temperature steps 
* ``Tmax``, ``tempmax``, ``Maxtemp`` or ``Max-Temperature``  is the maximal temperature in K (minimal T = 0K ) 
* ``list``: switch use the ``temperature-list`` construct. 
* ``units``: swtich to different unins Watts: Watt/Str/Molecule.
