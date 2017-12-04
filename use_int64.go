// +build !use_int32

package clipper

type CInt int64

const loRange CInt = 0x3FFFFFFF
const hiRange CInt = 0x3FFFFFFFFFFFFFFF //L;  // it won't compile with the L there
const CIntMin = CInt(-0x8000000000000000)
const CIntMax = CInt(0x7FFFFFFFFFFFFFFF)
