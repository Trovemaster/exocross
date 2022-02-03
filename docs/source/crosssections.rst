Cross-Sections
==============

The method is described in: 
C. Hill, S. N. Yurchenko, J. Tennyson, "Temperature-dependent molecular absorption cross sections for exoplanets and other atmospheres",  Icarus, 226, 1673-1677 (2013). 
`See the paper here`_.


.. _See the paper here: http://www.sciencedirect.com/science/article/pii/S0019103512003041



Gaussian profile
^^^^^^^^^^^^^^^^

Example (HWHM is the half-width at half-maximum):

::
    
    Temperature  2000.0
    Range 0.0  12000.0
    
    Npoints 200001
    absorption
    Gaussian  
    HWHM 0.5 (cm-1)
    threshold 1e-40 (to skip weak lines)
    
    output abs_gauss_0.5_T2000.0
    
    States  "../all/ch4-50.states-all"
    
    Transitions
     a-03000.dat
     a-03100.dat
     a-03200.dat
     a-03300.dat
     a-03400.dat
     a-03500.dat
     a-03600.dat
     a-03700.dat
     a-03800.dat
     a-03900.dat
    end
    


Doppler profile 
^^^^^^^^^^^^^^^

Doppler_ is the effective mass of the molecule in amu. 

.. _Doppler: https://github.com/Trovemaster/exocross/blob/master/img/alpha.png

Example:

::
    
    Temperature  1500.0
    Range 0.0  12000.0
    Npoints 200001
    
    emission
    doppl
    mass 16.0313

    output dop_emiss_1500.0
    States ch4-50.states-all

    Transitions a-02100.dat
    

Lorentzian profile
^^^^^^^^^^^^^^^^^^

Here b is the normalization factor. 

Example:

::

    Temperature  300
    Range 0.0  10000.0
    Npoints 10001
 
    absorption
    Loren  
    HWHM 0.1 (cm-1)
    
    output abs_lor_0.1_T300.0
    
    States  NaH.states
    
    Transitions NiH.trans


Stick spectrum
^^^^^^^^^^^^^^

Example:
::
    
    (ScH stick spectrum)
    Temperature 1500.0
    Range 0.  16000.0
    
    Npoints 16001
    
    absorption
    stick 
    threshold 1e-29
    
    output ScH_1500K_stick
    States       ScH.states
    Transitions  ScH.trans
        

bin 
^^^

is to produce average intensity per the wavenumber or wavelength interval as defined by 
Range/(Npoints-1). The wavelength is invoked by adding um to the range values. 

Example:
::

    
    (ScH bin spectrum)
    Temperature 1500.0
    Range 0.  16000.0
    
    Npoints 16001
    
    absorption
    bin  
    
    output ScH_1500K_bin_stick
    States       ScH.states
    Transitions  ScH.trans
    
 
or:
::
    
    (ScH bin spectrum)
    Temperature 1500.0
    Range 1.  100.0 um (or micron)
    

Box
^^^
Is to plot the maximal transition intensity per wavenumber interval, which is a cheaper alternative for the stick spectrum

Example:
::
    
    (ScH box spectrum)
    Temperature 1500.0
    Range 0.  16000.0
    
    Npoints 16001

    abundance 0.97
    
    absorption
    box
    threshold 1e-29
    
    output ScH_1500K_box_stick
    States       ScH.states
    Transitions  ScH.trans
    

Line-width cut-offs
^^^^^^^^^^^^^^^^^^^


A line width cut-off can be defined using `offset` or `line-cutoff`


::

    offset 25 (cm-1)


::

    line-cutoff 25 (cm-1)


where the offset value is in wavenumbers (cm-1). The default value is 25 cm-1. Alternatively, 
one can define the cut-off in terms of the HWHM as follows: 

::

    offset 50 HWHM



multi-grid
^^^^^^^^^^

    
A multi-grid with regions of different resolutions can be defined using the following `grid` section:

::     
    
    grid
      Range   0    100   Npoints 10000 offset 10 
      Range 100   1000  Npoints 1000  offset 25
      Range 1000 10000  Npoints 100
    end
     

The maximal number of sub-grids is 100. Currently this option only works with 
simple sampling-type profiles, such as `Voigt`, `Doppler Sampling`,  `Gaussian Sampling` or `Bin`. 
The latter is commonly used to generate super-lines.  
`offset` or `line-cutoff` is an optional keyword to allow region-dependent cutoffs for line profiles. If undefined, the value of the 
global keyword `offset`  the corresponding default value (25 cm-1) is used.



gf line list
^^^^^^^^^^^^

A stick spectrum is produced with the gf-factors in place of the Einstein coefficients. Here is the example for the VALD format, where 
the columns are the wavelength in Angstrom, the lower state energy in eV, log10(gf), 0.0, the statistical weight 2J'+1 (upper state J') and zero. 

::

    temperature 5000     
    Range 100.  16000.0
    
    gf
    vald
    threshold 1e-29
    
    output ScH_gf
    States       ScH.states
    Transitions  ScH.trans
