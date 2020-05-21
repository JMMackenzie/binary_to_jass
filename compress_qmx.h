/*
	COMPRESS_QMX.H
	--------------
*/
#ifndef COMPRESS_QMX_H_
#define COMPRESS_QMX_H_

#include "compress.h"

/*
	class ANT_COMPRESS_QMX
	----------------------
*/
class ANT_compress_qmx : public ANT_compress
{
private:
	uint8_t *length_buffer;
	uint64_t length_buffer_length;

public:
	ANT_compress_qmx();
	virtual ~ANT_compress_qmx();

	virtual void encodeArray(const uint32_t *in, uint64_t len, uint32_t *out, uint64_t *nvalue);
	virtual void decodeArray(const uint32_t *in, uint64_t len, uint32_t *out, uint64_t nvalue);
} ;

#endif

