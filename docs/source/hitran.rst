HITRAN
======

ExoCross can be used with the standard HITRAN line list by adding the `HITRAN` keyword anywhere in the input file (outside any sections). This will also require the partition function (`pf`) and isotopologue number (`iso`) defined. 

HITRAN broadening parameters will be used unless the species-section is given, which specifies the broadening.  The .states file is not required and ignored if given. The HITRAN total statistical weights are used directly.

For example: 
::

    Temperature  296.0
    Range 0.0  10000.0
    
    Npoints 10001
    
    absorption
    gauss
    hwhm 0.5 (cm-1)
    
    hitran
    iso 26 1
    pf 1000.0
    
    output C2H2_ab_g0.5
    Transitions  26_hit12.par
    
    

Here `iso`  can appear in a HITRAN form as, e.g. 261. 

HITRAN keyword can also also used for writing in the HITRAN output. In this case the broadening parameters are expected to be specified for air and self as part of the SPECIES section, including delta. To invoke HITRAN output use WRITE next to the HITRAN keyword, which starts a section, which should be ended either by END or an emty line. The ID_ISO should be also specified using ISO keyword.
::

    
    Temperature  1900.0
    Range 0.0  10000.0
    
    absorption
    stick
    threshold -1e-25
    
    hitran write

    iso 28
    
    pressure  1.0
    species 
         air   gamma 0.0155 n 0.41 t0 298.0  ratio 0.82 delta 0.000
         self  gamma 0.1070 n 0.77 t0 298.0  ratio 0.18 delta 0.000
    end
    
    mass 40.0
    
    output NO_1900K_HITRAN
    
    States 15N16O.states
    
    Transitions  NO_128_N15.trans.head
    

    
The following example allows one to specify the error codes (6 values) for the HITRAN format outputs. It is still a part of the HITRAN WRITE section, the six keywords (`error-E`, `error-S`, `error-Air`, `error-self`, `error-delta`) start lines with the corresponding specifications. The Energy and Intensity (S) lines are used to give the ranges of he quantum numbers for different error codes. The other four can hold only one error code.

::
          
     hitran write
      error-E  qn 4 ierr 4  vmax  10 ierr 3  vmax  20 ierr 2  vmax  30  ierr 1  vmax  40 ierr 0  vmax  100 
      error-S  qn 4 ierr 5  vmax  10 ierr 4  vmax  20 ierr 3  vmax  100 
      error-Air   ierr 4
      error-self  ierr 4
      error-n     ierr 4
      error-delta ierr 0
     end
     
     
Here `write` indicates that the HITRAN-stick-like list will be printed. 

`qn 4` indicates the quantum number (4th column in States after J's columns);

`ierror 4 vmax 10` means that the error code 4 is applied for all values of `qn` less or equal to 10 etc. 

`error-Air   ierr 4` indicates that the error code for the Air-broadening is 4. 


The absorption spectrum using HITRAN can be computed using the HITRAN reference intensity:

:math:`I(T)=I(T_{\rm ref}) \frac{ \exp(-c_0 \tilde{E}''/T) \left(1-\exp(-c_0\tilde{\nu}/T)\right) Q(T_{\rm ref})}{\exp(-c_0 \tilde{E}''/T_{\rm ref}) \left(1-\exp(-c_0\tilde{\nu}/T_{\rm ref})\right) Q(T)}`

where :math:`T_{\rm ref}` and :math:`Q(T_{\rm ref})` are the reference temeprature (296 K) and the corresponding reference partition function. 
This expression is always used if :math:`Q(T_{\rm ref})` is provided:
The reference temperature is assumed 296 K (HITRAN standard). It can be changed via 
::

    
    pf 100000. ref 30000.

Otherwsie the standard absorption expression is amployed. The reference temperature is assumed 296 K (HITRAN standard). It can be changed via 
::

    
    Temperature  1000 ref 173 



The intensity cut-off (stick) can be done using the HITRAN method: 
:math:`S=S_{0} \tanh(c_2 \nu/2T)` for :math:`\nu\le 2000` cm :math:`^{-1}` and :math:`10^{-29}` cm/molecule above. 

::     
    
    absorption
    stick
    cutoff HITRAN (S_crit) 1e-29  nu_crit 2000 
    
    output ScH_1500K_box_stick
    States       ScH.states
    Transitions  ScH.trans
     


VALD
====

The Vald format has the following structure: 

lambda   E_low_eV   loggf  0.0  gtot_f  0.0

where lambda is the wavelength in nm, E_low_eV is the energy of the lower state in eV, loggf is  log10(gf), gf is the gf factor 
in the Physics (ExoMol) convention, gtot_f is the total degeneracy in the ExoMol convention. To convert to the Astrophysics convention, 
the keyword `gf_factor` canm be used. 

Here is an example to generate a Vald-formatted line list for 13C__8states using `gf_factor 0.25` to convert to he Astrophysics convention. Here, the  
nuclear statistical weights are 1 (A1,A2) and 3 (B1,B3), which gives the factor 1/(1+3) = 1/4. 


::
    
    Temperature  8000
    Range 1000 30000
    
    vald
    gf
    
    gf_factor 0.25
    threshold 1e-99
    
    output 13C2__8states_vald_T8000K_1e-99_10000-30000

    States 13C2__8states.states
    Transitions 13C2__8states.trans
    





