/*
	COMPRESS.H
	----------
*/
#ifndef COMPRESS_H_
#define COMPRESS_H_

#include <stdint.h>
#define ANT_compressable_integer uint32_t
/*
	class ANT_COMPRESS
	------------------
*/
class ANT_compress
{
public:
	ANT_compress() {}
	virtual ~ANT_compress() {}

	virtual void encodeArray(const uint32_t *in, uint64_t len, uint32_t *out, uint64_t *nvalue) = 0;
	virtual void decodeArray(const uint32_t *in, uint64_t len, uint32_t *out, uint64_t nvalue) = 0;

} ;

#endif  /* COMPRESS_H_ */

