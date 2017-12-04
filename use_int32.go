// +build use_int32

package clipper

type CInt int32

const loRange CInt = 46340
const hiRange CInt = 46340
const CIntMin = CInt(-0x80000000)
const CIntMax = CInt(0x7FFFFFFF)
