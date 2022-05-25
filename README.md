
# Deuteron Solution

The final write-up of the project can be found [here in the repository](https://github.com/eli-temanson/Nuc1-Deuteron-Project/blob/main/Nuclear_Physcs_1_Project.pdf).

## installation
### get julia

https://julialang.org/
https://medium.com/coffee-in-a-klein-bottle/install-julia-1-5-on-ubuntu-bb8be4b2571d

### install some useful packages

start julia and install depencies (this may take time)

```
julia dep.jl
```
## prepare fortran potentials to be called by julia

```
gfortran -shared -fPIC -o av18pot.so av18pot.f
```

## running the program

Since julia uses just intime compilation, loading packages and reloading functions actually can take some time. It's best to start julia and then inside the terminal run...

```
julia
```

```
inclue("deuteron.jl")
``` 
