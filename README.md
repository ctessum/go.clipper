go.clipper
==========

BUILD
-----

*int32 version*
```
$ go build -tags 'use_int32'
```

*int64 version*
```
$ go build -tags '!use_int32'
# or combine build tags
$ go build -tags '!use_int32,!use_xyz'
```
