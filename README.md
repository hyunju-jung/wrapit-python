# Wrap UCLA-LES subroutines

## how to wrap it

```
f2py -c ucla_srfc.f90 -m ucla_srfc.so
```

## call the module in Python

```
import ucla_srfc
```
or 

```
from ucla_srfc import srfcscls

ustart, tstart, rstart = srfcscls(n2, n3, z, u, dth, drt)
```

[The original UCLA-LES code](https://github.com/uclales/uclales)
