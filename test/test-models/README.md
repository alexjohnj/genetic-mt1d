Several different models for testing the forward modelling code. The expected
responses have been calculated using a working version of `calc_response`. I've
checked that the phase and apparent resistivity curves are self-consistent, that
the curves actually make sense, and visually compared the output to the results
from a similar forward modelling implementation by Martyn Unsworth (University
of Alberta).

## Model 1

```plain
-------------------------------------------------------------------------------
                                     100 Ωm
```

100 Ωm half-space model. Expected resistivity response is a constant 100 Ωm
across all frequencies. Expected phase response is a constant 45° across all
frequencies.

## Model 2

```plain
-------------------------------------------------------------------------------
                                     100 Ωm
--------------------------------------------------------------------------z=1km
                                    1000 Ωm
```

2 layer model with resistivity increasing with depth. Resistivity graph should
start at 100 Ωm and end on 1000 Ωm. Change should occur around 25 Hz based on
skin-depth calculations. Phase angle should dip below 45° at 25 Hz and return to
45° at the lowest frequencies.

## Model 3

```plain
--------------------------------------------------------------------------------
                                     100 Ωm
---------------------------------------------------------------------------z=1km
                                     10 Ωm
---------------------------------------------------------------------------z=2km
                                     100 Ωm
```

3 layer model with (relatively) conductive layer sandwiched between two
resistive layers. Resistivity response should be constant up to 25 Hz, decrease
towards 10 Ωm and then increase towards 100 Ωm around 10 Hz. Phase should be
constant up to 25 Hz, then increase above 45°, decrease below 45° and finally
flatten out at 45°.

## Model 4

```plain
--------------------------------------------------------------------------------
                                     100 Ωm
---------------------------------------------------------------------------z=1km
                                      1 Ωm
-------------------------------------------------------------------------z=1.1km
                                     100 Ωm
--------------------------------------------------------------------------z=40km
                                     0.1 Ωm
------------------------------------------------------------------------z=40.1km
                                     100 Ωm
```

5 layer model featuring two highly conductive layers sandwiched between
resistive layers. At high frequencies, resistivity response should be 100
Ωm. From skin-depth calculations, resistivity should start to decrease at 25 Hz,
and as the 2nd layer has a high conductance, should get really low. Resistivity
should start to return to 100 Ωm before decreasing again at around 0.01
Hz. Finally, it should start to increase towards 100 Ωm again.

The phase curve should be constant up to the skin-depth. It should then show two
peaks > 45° at around 10 Hz and 0.01 Hz. The higher frequency peak should be
larger than the lower frequency peak.

## Model 5


```plain
--------------------------------------------------------------------------------
                                     100 Ωm
---------------------------------------------------------------------------z=1km
                                      1 Ωm
-------------------------------------------------------------------------z=1.1km
                                     100 Ωm
--------------------------------------------------------------------------z=40km
                                      1 Ωm
------------------------------------------------------------------------z=40.1km
                                     100 Ωm
```

This model is similar to model 4 except the 4th layer's resistivity is increased
by a factor of ten. The conductance of this layer is now the same as the 2nd
layer, so it shouldn't produce any response in the apparent resistivity or phase
curves.
