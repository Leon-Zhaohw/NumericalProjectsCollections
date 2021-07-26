#ifndef __INCLUDED_RENDERLINES_H__
#define __INCLUDED_RENDERLINES_H__

namespace planning {

void renderLines( 
	const std::string& linesFilename,
	const vl::Mat4f& transformMatrix,
	const std::string outFilePrefix,
	size_t width,
	size_t height );

} // namespace planning

#endif
