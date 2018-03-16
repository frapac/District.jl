# Models

District.jl implements various devices, and allows to combine them
together to build complex networks.


## Time period

First, we specify the time period of interest by:
```@docs
TimeSpan(day, ndays)
```


## Weather data

Load given weather data:
```@docs
loadweather
```


## Uncertainties modeling


```@docs
loadnoise
```

```@docs
District.nnoise
```

```@docs
District.optscenarios
```

```@docs
District.genscenarios
```

```@docs
District.genforecast
```

```@docs
District.fit
```

```@docs
District.getname
```


### Quantization functions

```@docs
District.normalquantization(n)
```


## Devices

```@docs
District.elecload
```
```@docs
District.gasload
```
```@docs
District.thermalload
```
```@docs
District.parsedevice
```
```@docs
District.nstates
```
```@docs
District.ncontrols
```
```@docs
District.xbounds
```
```@docs
District.ubounds
```
```@docs
District.isstock
```
```@docs
District.getname
```
```@docs
District.reset!
```
