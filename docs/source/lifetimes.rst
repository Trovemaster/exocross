Computing Lifetimes
===================

The lifetimes are computed by summing Einstein A coefficients as

:math:`\tau_{i'}= \sum_{i''< i'} \frac{1}{A_{i',i''}}`

For states with no channels to lower states, the lifetimes are considered undefined and assign NaN. This feature is based on the external
library IEEE. In the older version of ExoCross, such states were assigned -1.0000.



Example  of the input file for lifetimes::


    lifetime

    output 28SiO-16

    States 28Si-16O__EBJT.states
    Transitions 28Si-16O__EBJT.trans


By default, the lifetimes are written into the column 6 assuming that the uncertainty column is present which should be always stays as column 5. If however the uncertainties are not given and the lifetimes should be written into column, this option can be specified using the QN section:
:: 
       
     QN
      lifetime 5
     END
     

Lifetimes as input for line profiles 
------------------------------------

If pre-dissociation is considered, the lifetime can be used to estimate the effective line width. This option can be activated using the following input
::
     
     predissociation

     QN
      lifetime 5
     END
     

For more details, see the "Line profiles" section. 



