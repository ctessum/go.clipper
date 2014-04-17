// +build !use_int32

package clipper

type cInt int64

const loRange cInt = 0x3FFFFFFF
const hiRange cInt = 0x3FFFFFFFFFFFFFFF //L;  // it won't compile with the L there
