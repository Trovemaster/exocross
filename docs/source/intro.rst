Introduction 
============

ExoCross is a Fortran code for generating spectra (emission, absorption) and 
thermodynamic properties (partition function, specific heat etc.) from 
molecular line lists. Input is taken in several formats, including ExoMol and 
HITRAN formats. ExoCross is efficiently parallelized showing also a high 
degree of vectorization. It can work with several line profiles such as 
Doppler, Lorentzian and Voigt and support several broadening schemes. Voigt 
profiles are handled by several methods allowing fast and accurate 
simulations. Two of these methods are new. ExoCross is capable of working 
with the recently proposed method of super-lines. It supports calculations of 
lifetimes, cooling functions, specific heats and other properties. ExoCross 
can be used to convert between different formats, such as HITRAN, ExoMol and 
Phoenix. It is capable of simulating non-LTE spectra using a simple 
two-temperature approach. Different electronic, vibronic or vibrational bands 
can be simulated separately using an efficient filtering scheme based on the 
quantum numbers.  
