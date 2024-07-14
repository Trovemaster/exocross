Modelling continuum
===================

For continuum states spectra calculations, i.e. for transitions between bound and unbound (usually upper) states,  the corresponding line lists are represented by discretized transitions with their intensities, which need to be re-distributed across the simulation wavenumber regions. This is typically done with a Gaussian line profile with a larger width (50-200 cm-1) to cover gaps between these discrete lines. For example, 
::

   Temperature  5000
   Range    45000 85000

   npoints  40000

   absorption
   gaussian
   
   hwhm 200 

   cutoff 5000

   output NH_box_5000K_L20

   States NH.states

   Transitions NH.trans
   

It is also important to set ``cutoff`` (wing cut-off) to a large number.

'Particle in the box'
---------------------

In the case of a very large spectroscopic range, the distance between discrete continuum lines rapidly increases. For example, for the case of the particle-in-a-box solution, the energy separation increases with the excitation number linearly:

.. math::


    \Delta \tilde{E}_n^{\rm box} = \tilde{E}_{n+1}^{\rm box} - \tilde{E}_n^{\rm box} = \frac{h (2n+1)}{8 \mu L^2 c}.

where :math:`\mu` is the reduced mass, :math:`L` is the box size) and :math:`n\ge 1` is the state counting number, :math:`c` is the speed of light in cgs.
Due to the linear dependence on :math:`n`,  there is no a single optimal value of :math:`\alpha_{\rm G}` (Gaussian HWHM)  for the entire region. 
In order to account for such an increase, the 'recipe' ``Box'' for continuum states spectra calculations can be used as part of the ``Species`` construct. The corresponding line lists are represented by discretized transitions with their intensities to be re-distributed across the simulation wavenumber regions. In the ``Box`` recipe, we define  :math:`\alpha_{\rm G}` to be :math:`\gamma_n^{\rm box}  = \gamma_0 \Delta \tilde{E}_n^{\rm box}`, i.e.

.. math::


    \gamma_n^{\rm box} = \gamma_0 \frac{h (2n+1)}{8 \mu L^2 c},


where :math:`n` is the counting number of the continuum states (from the lowest) and :math:`\gamma_0` is the initial value of HWHM for :math:`n=1`.


Here is an example of the ExoCross input:
::

   Temperature  5000
   Range    45000 85000

   npoints  40000

   QN
     K 8
   end

   absorption
   gaussian

   cutoff 5000

   species
        particle  gamma 1.0 mass 0.937  Lbox 20  model box
   end

   output NH_box_5000K_L20

   States NH.states

   Transitions NH.trans


Here  ``model`` ``box`` defines the recipe type ``box``, ``Lbox`` is the size of the box and ``mass`` specifying the reduced mass, ``particle`` is an arbitrary dummy name for the broadening and ``gamma`` defines the initial HWHM value :math:`\gamma_0`. 

