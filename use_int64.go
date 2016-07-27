// +build !use_int32

package clipper

type CInt int64

const loRange CInt = 0x3FFFFFFF
const hiRange CInt = 0x3FFFFFFFFFFFFFFF //L;  // it won't compile with the L there
