/*
 *	exporter3.cpp
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "exporter3.h"
#include "particle3.h"
#include "octree3.h"
#include "object3.h"
#include "camera3.h"
#include "bcc3.h"
#include "bccmesher3.h"
#include "util3.h"
#include "extsurf3.h"
#include <sys/stat.h>
#include <string>
#include <zlib.h>
#include <stdint.h>
using namespace std;

void exporter3::writeMesh( uint frame, const surf3 &surf, const std::vector<particle3 *> &particles, FLOAT64 dpx, const mesher3 *mesher, const octree3 *octree,
						   const camera3 *camera, const std::vector<object3 *> *objects, bool writeParticleData ) {
	// Make sure 'src' directory exist
	if( ! is_exist("src" )) {
		printf( "src directory is missing !\n" );
		exit(0);
	}
	
	if( frame == 0 ) {
		// Make Directory
		run("cp src/mkmovie.sh %s/mkmovie.sh", root_path );
	}
	
	// Write Mitsuba file
	tick(); dump("Exporting %s Mitsuba files...", nth(frame));
	write_mitsuba(frame,surf,particles,dpx,writeParticleData,camera,objects,mesher);
	dump("Done. Took %s.\n",stock());
	
	// Write particles
	if( writeParticleData ) {
		tick(); dump("Exporting %s particle file...", nth(frame));
		if( ! is_exist("%s/particles",root_path) ) {
			run("mkdir %s/particles",root_path);
		}
		extsurf3 extsurf;
		extsurf.writeFile(format_str("%s/particles/%d_particles_tmp.dat.gz",root_path,frame),particles,dpx,octree);
		run("mv %s/particles/%d_particles_tmp.dat.gz %s/particles/%d_particles.dat.gz", root_path, frame, root_path, frame );
		dump("Done. Took %s.\n",stock());
	}
}

static void append( const void *data, uint size, std::vector<uint8_t> &buffer ) {
	buffer.insert(buffer.end(),(uint8_t *)data,(uint8_t *)data+size);
}

static void compress_memory(void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data) {
	std::vector<uint8_t> buffer;
	
	const size_t BUFSIZE = 128 * 1024;
	uint8_t temp_buffer[BUFSIZE];
	
	z_stream strm;
	strm.zalloc = 0;
	strm.zfree = 0;
	strm.next_in = reinterpret_cast<uint8_t *>(in_data);
	strm.avail_in = in_data_size;
	strm.next_out = temp_buffer;
	strm.avail_out = BUFSIZE;
	
	deflateInit(&strm, Z_BEST_COMPRESSION);
	
	while (strm.avail_in != 0)
	{
		int res = deflate(&strm, Z_NO_FLUSH);
		assert(res == Z_OK);
		if (strm.avail_out == 0)
		{
			buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUFSIZE);
			strm.next_out = temp_buffer;
			strm.avail_out = BUFSIZE;
		}
	}
	
	int deflate_res = Z_OK;
	while (deflate_res == Z_OK)
	{
		if (strm.avail_out == 0)
		{
			buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUFSIZE);
			strm.next_out = temp_buffer;
			strm.avail_out = BUFSIZE;
		}
		deflate_res = deflate(&strm, Z_FINISH);
	}
	
	assert(deflate_res == Z_STREAM_END);
	buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUFSIZE - strm.avail_out);
	deflateEnd(&strm);	
	out_data.swap(buffer);
}

void exporter3::write_mitsuba(uint frame, const surf3 &surf, const std::vector<particle3 *> &particles, FLOAT64 dpx, bool writeParticleData,
							  const camera3 *camera, const std::vector<object3 *> *objects, const mesher3 *tet ) {
	if( ! is_exist("%s/mitsuba", root_path) ) {
		// Make Directory
		run("mkdir %s/mitsuba/", root_path );
		run("mkdir %s/mitsuba/img/", root_path );
		run("cp src/mitsuba/wireframe.xml %s/mitsuba/", root_path);
		if( writeParticleData ) run("cp src/mitsuba/particles.xml %s/mitsuba/", root_path);
		run("cp src/mitsuba/opaque.xml %s/mitsuba/", root_path);
		run("cp src/mitsuba/fancy.xml %s/mitsuba/", root_path);
		run("cp src/mitsuba/envmap_wire.exr %s/mitsuba/", root_path);
		if( writeParticleData ) run("cp src/mitsuba/envmap_particle.exr %s/mitsuba/", root_path);
		run("cp src/mitsuba/render_o.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/render_f.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/render_w.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/render_view.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/render_top.sh %s/mitsuba/", root_path);
		if( writeParticleData ) run("cp src/mitsuba/render_p.sh %s/mitsuba/", root_path);
		run("chmod 755 %s/mitsuba/render_w.sh", root_path);
		run("chmod 755 %s/mitsuba/render_o.sh", root_path);
		run("chmod 755 %s/mitsuba/render_f.sh", root_path);
		run("chmod 755 %s/mitsuba/render_view.sh", root_path);
		run("chmod 755 %s/mitsuba/render_top.sh", root_path);
		if( writeParticleData ) run("chmod 755 %s/mitsuba/render_p.sh", root_path);
		run("cp src/mitsuba/fix_o.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/fix_view.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/fix_f.sh %s/mitsuba/", root_path);
		run("cp src/mitsuba/fix_w.sh %s/mitsuba/", root_path);
		if( writeParticleData ) run("cp src/mitsuba/fix_p.sh %s/mitsuba/", root_path);
		run("chmod 755 %s/mitsuba/fix_w.sh", root_path);
		run("chmod 755 %s/mitsuba/fix_o.sh", root_path);
		run("chmod 755 %s/mitsuba/fix_view.sh", root_path);
		run("chmod 755 %s/mitsuba/fix_f.sh", root_path);
		if( writeParticleData ) run("chmod 755 %s/mitsuba/fix_p.sh", root_path);
		
		FILE *empty_fp = fopen( format_str("%s/mitsuba/empty.xml", root_path), "w" );
		fprintf(empty_fp,"<?xml version='1.0' encoding='utf-8'?>\n" );
		fprintf(empty_fp,"<scene version=\"0.4.0\">\n" );
		fprintf(empty_fp,"</scene>\n" );
		fclose(empty_fp);
	}
	
	// Write camera info like this: <lookat target="0.5, 0.5, 0.3" origin="0.5, -1.3, 2.5" up="0, 0, 1"/>
	if( camera && camera->isEnabled() ) {
		FILE *camera_fp = fopen( format_str("%s/mitsuba/camera_%d.txt", root_path, frame), "w" );
		vec3d origin = camera->origin;
		vec3d target = camera->target;
		vec3d up = camera->up;
		fprintf(camera_fp, "%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n",
				target[0], target[2], target[1], origin[0], origin[2], origin[1], 0.0, 0.0, 1.0 );
		fclose(camera_fp);
	}
		
	// If objects dir is not found
	if( ! is_exist("%s/mitsuba/obj", root_path) ) {
		run("mkdir %s/mitsuba/obj", root_path );
		
		FILE *objects_fp = fopen( format_str("%s/mitsuba/objects.xml", root_path), "w" );
		fprintf(objects_fp,"<?xml version='1.0' encoding='utf-8'?>\n" );
		fprintf(objects_fp,"<scene version=\"0.4.0\">\n\n" );		
		// Export objects available
		if( objects ) for( uint n=0; n<objects->size(); n++ ) {
			const object3 &object = *(*objects)[n];
			std::vector<vec3d> vertices;
			std::vector<vec3i> faces;
			object.getPolygon(vertices,faces);
			if( object.type == object3::SOLID && vertices.size() && faces.size() ) {
				// Write an entry
				fprintf(objects_fp,"<shape type=\"ply\">\n");
				fprintf(objects_fp,"<!-- <transform name=\"toWorld\"><translate x=\"0\" y=\"0\" z=\"-0.75\"/></transform> -->\n");
				fprintf(objects_fp,"<boolean name=\"faceNormals\" value=\"true\"/>\n");
				fprintf(objects_fp,"<string name=\"filename\" value=\"obj/objects_%d.ply\"/>\n",n);
				fprintf(objects_fp,"<bsdf type=\"twosided\">\n");
				fprintf(objects_fp,"<bsdf type=\"diffuse\">\n");
				fprintf(objects_fp,"<srgb name=\"reflectance\" value=\"#FFFF77\"/>\n");
				fprintf(objects_fp,"</bsdf>\n");
				fprintf(objects_fp,"</bsdf>\n");
				fprintf(objects_fp,"</shape>\n\n");
				
				// Now write PLY file
				FILE *ply_fp = fopen( format_str("%s/mitsuba/obj/objects_%d.ply", root_path, n), "w" );
				
				// Write header
				fprintf(ply_fp,"ply\n");
				fprintf(ply_fp,"format binary_little_endian 1.0\n");
				fprintf(ply_fp,"element vertex %d\n", (int)vertices.size());
				fprintf(ply_fp,"property float x\n");
				fprintf(ply_fp,"property float y\n");
				fprintf(ply_fp,"property float z\n");
				fprintf(ply_fp,"element face %d\n", (int)faces.size());
				fprintf(ply_fp,"property list uchar int vertex_indices\n");
				fprintf(ply_fp,"end_header\n");
				fclose(ply_fp);
				
				// Write vertices
				ply_fp = fopen( format_str("%s/mitsuba/obj/objects_%d.ply", root_path, n), "ab" );
				for( uint n=0; n<vertices.size(); n++ ) {
					FLOAT32 v[3] = { vertices[n][0], vertices[n][2], vertices[n][1] };
					fwrite( &v[0], 1, sizeof(FLOAT32), ply_fp );
					fwrite( &v[1], 1, sizeof(FLOAT32), ply_fp );
					fwrite( &v[2], 1, sizeof(FLOAT32), ply_fp );
				}
				
				// Write faces
				for( uint n=0; n<faces.size(); n++ ) {
					unsigned char num = 3;
					fwrite( &num, 1, sizeof(unsigned char), ply_fp );
					fwrite( &faces[n][2], 1, sizeof(int), ply_fp );
					fwrite( &faces[n][1], 1, sizeof(int), ply_fp );
					fwrite( &faces[n][0], 1, sizeof(int), ply_fp );
				}
				fclose(ply_fp);
			}
		}
		
		fprintf(objects_fp,"</scene>\n" );
		fclose(objects_fp);
	}
	
	// Get surface mesh
	const std::vector<vec3d> &vertices = surf.vertices;
	const std::vector<vec3i> &faces = surf.faces;
	
#if 1
	// Now write PLY file
	FILE *ply_fp = fopen( format_str("%s/mitsuba/%d_scene_tmp.ply", root_path, frame), "w" );
	
	// Write header
	fprintf(ply_fp,"ply\n");
	fprintf(ply_fp,"format binary_little_endian 1.0\n");
	fprintf(ply_fp,"element vertex %d\n", (int)vertices.size());
	fprintf(ply_fp,"property float x\n");
	fprintf(ply_fp,"property float y\n");
	fprintf(ply_fp,"property float z\n");
	fprintf(ply_fp,"element face %d\n", (int)faces.size());
	fprintf(ply_fp,"property list uchar int vertex_indices\n");
	fprintf(ply_fp,"end_header\n");
	fclose(ply_fp);
	
	// Write vertices
	ply_fp = fopen( format_str("%s/mitsuba/%d_scene_tmp.ply", root_path, frame), "ab" );
	for( uint n=0; n<vertices.size(); n++ ) {
		FLOAT32 v[3] = { vertices[n][0], vertices[n][2], vertices[n][1] };
		fwrite( &v[0], 1, sizeof(FLOAT32), ply_fp );
		fwrite( &v[1], 1, sizeof(FLOAT32), ply_fp );
		fwrite( &v[2], 1, sizeof(FLOAT32), ply_fp );
	}
	
	// Write faces
	for( uint n=0; n<faces.size(); n++ ) {
		unsigned char num = 3;
		fwrite( &num, 1, sizeof(unsigned char), ply_fp );
		fwrite( &faces[n][2], 1, sizeof(int), ply_fp );
		fwrite( &faces[n][1], 1, sizeof(int), ply_fp );
		fwrite( &faces[n][0], 1, sizeof(int), ply_fp );
	}
	fclose(ply_fp);
	run("mv %s/mitsuba/%d_scene_tmp.ply %s/mitsuba/%d_scene.ply", root_path, frame, root_path, frame );
#endif
	
#if 1
	// Now write serialized file
	FILE *serialized_fp = fopen( format_str("%s/mitsuba/%d_scene_tmp.serialized", root_path, frame), "wb" );
	uint16_t format = 0x041C;
	uint16_t version = 0x0004;
	fwrite(&format,1,sizeof(uint16_t),serialized_fp);
	fwrite(&version,1,sizeof(uint16_t),serialized_fp);
	
	std::vector<uint8_t> buffer;
	uint32_t bit = 0x2000;
	append(&bit,sizeof(uint32_t),buffer);
	const char *name = "liquid";
	//gzwrite(serialized_zfp, name, strlen(name)+1);
	append(name, strlen(name)+1, buffer);
	uint64_t v_number = vertices.size();
	uint64_t f_number = faces.size();
	append(&v_number,sizeof(uint64_t),buffer);
	append(&f_number,sizeof(uint64_t),buffer);
	
	for( uint n=0; n<vertices.size(); n++ ) {
		FLOAT64 v[3] = { vertices[n][0], vertices[n][2], vertices[n][1] };
		append(&v[0],sizeof(FLOAT64),buffer);
		append(&v[1],sizeof(FLOAT64),buffer);
		append(&v[2],sizeof(FLOAT64),buffer);
	}
	for( uint n=0; n<faces.size(); n++ ) {
		uint32_t f[3] = { faces[n][0], faces[n][2], faces[n][1] };
		append(&f[0],sizeof(uint32_t),buffer);
		append(&f[1],sizeof(uint32_t),buffer);
		append(&f[2],sizeof(uint32_t),buffer);
	}
	std::vector<uint8_t> out_data;
	compress_memory(&buffer[0], buffer.size(), out_data);
	fwrite(&out_data[0],1,out_data.size(),serialized_fp);
	fclose(serialized_fp);
	
	run("mv %s/mitsuba/%d_scene_tmp.serialized %s/mitsuba/%d_scene.serialized", root_path, frame, root_path, frame );
#endif
	
	// Write BCC mesh if requested
	if( tet ) {
		// Get surface mesh
		const std::vector<vec3d> &vertices = tet->nodes;
		const std::vector<std::vector<uint> > &elements = tet->elements;
		
		// Now write tet PLY file
		FILE *tet_fp = fopen( format_str("%s/mitsuba/%d_tet_tmp.ply", root_path, frame), "w" );
		
		std::vector<int> vremap(vertices.size());
		for( uint n=0; n<vertices.size(); n++ ) {
			vremap[n] = -1;
		}
		
		// Count the number of faces
		int cut_dim = 2;
		FLOAT64 cut_pos = 0.35;
		int num_elements = 0;
		for( uint n=0; n<elements.size(); n++ ) {
			const std::vector<uint> &element = elements[n];
			FLOAT64 max_pos = 0.0;
			FLOAT64 min_pos = 1.0;
			for( uint m=0; m<4; m++ ) {
				max_pos = fmax(max_pos,vertices[element[m]][cut_dim]);
				min_pos = fmin(min_pos,vertices[element[m]][cut_dim]);
			}
			if( max_pos >= cut_pos && min_pos <= cut_pos ) {
				for( uint m=0; m<4; m++ ) {
					for( uint v0=0; v0<4; v0++ ) for( uint v1=v0+1; v1<4; v1++ ) for( uint v2=v1+1; v2<4; v2++ ) {
						num_elements ++;
						vremap[element[v0]] = 1;
						vremap[element[v1]] = 1;
						vremap[element[v2]] = 1;
					}
				}
			}
		}
		
		uint v_index = 0;
		for( uint n=0; n<vertices.size(); n++ ) {
			if(vremap[n]==1) {
				vremap[n] = v_index;
				v_index ++;
			}
		}
		
		// Write header
		fprintf(tet_fp,"ply\n");
		fprintf(tet_fp,"format binary_little_endian 1.0\n");
		fprintf(tet_fp,"element vertex %d\n", v_index);
		fprintf(tet_fp,"property float x\n");
		fprintf(tet_fp,"property float y\n");
		fprintf(tet_fp,"property float z\n");
		fprintf(tet_fp,"element face %d\n", num_elements);
		fprintf(tet_fp,"property list uchar int vertex_indices\n");
		fprintf(tet_fp,"end_header\n");
		fclose(tet_fp);
		
		// Write vertices
		tet_fp = fopen( format_str("%s/mitsuba/%d_tet_tmp.ply", root_path, frame), "ab" );
		for( uint n=0; n<vertices.size(); n++ ) {
			if( vremap[n] >= 0 ) {
				FLOAT32 v[3] = { vertices[n][0], vertices[n][2], vertices[n][1] };
				fwrite( &v[0], 1, sizeof(FLOAT32), tet_fp );
				fwrite( &v[1], 1, sizeof(FLOAT32), tet_fp );
				fwrite( &v[2], 1, sizeof(FLOAT32), tet_fp );
			}
		}
		
		// Write Faces
		for( uint n=0; n<elements.size(); n++ ) {
			const std::vector<uint> &element = elements[n];
			vec3d avg_pos;
			FLOAT64 max_pos = 0.0;
			FLOAT64 min_pos = 1.0;
			for( uint m=0; m<4; m++ ) {
				max_pos = fmax(max_pos,vertices[element[m]][cut_dim]);
				min_pos = fmin(min_pos,vertices[element[m]][cut_dim]);
				avg_pos += vertices[element[m]] / 4.0;
			}
			if( max_pos >= cut_pos && min_pos <= cut_pos ) {
				for( uint m=0; m<4; m++ ) {
					for( uint v0=0; v0<4; v0++ ) for( uint v1=v0+1; v1<4; v1++ ) for( uint v2=v1+1; v2<4; v2++ ) {
						unsigned char num = 3;
						int face[3] = { element[v0], element[v1], element[v2] };
						if( vremap[face[2]] < 0 || vremap[face[1]] < 0 || vremap[face[0]] < 0 ) {
							dump( "Vertex remap faild...\n");
							exit(0);
						}
						
						vec3d center = (vertices[face[0]]+vertices[face[1]]+vertices[face[2]])/3.0;
						bool flip = ((vertices[face[2]]-vertices[face[0]]) ^ (vertices[face[1]]-vertices[face[0]])) * (avg_pos-center) > 0.0;
						fwrite( &num, 1, sizeof(unsigned char), tet_fp );
						if( flip ) {
							fwrite( &vremap[face[2]], 1, sizeof(int), tet_fp );
							fwrite( &vremap[face[1]], 1, sizeof(int), tet_fp );
							fwrite( &vremap[face[0]], 1, sizeof(int), tet_fp );
						} else {
							fwrite( &vremap[face[0]], 1, sizeof(int), tet_fp );
							fwrite( &vremap[face[1]], 1, sizeof(int), tet_fp );
							fwrite( &vremap[face[2]], 1, sizeof(int), tet_fp );
						}
					}
				}
			}
		}
		
		fclose(tet_fp);
		run("mv %s/mitsuba/%d_tet_tmp.ply %s/mitsuba/%d_tet.ply", root_path, frame, root_path, frame );
	}
	
	// Write particles outside the surface
	if( writeParticleData ) {
		FILE *particle_fp = fopen( format_str("%s/mitsuba/%d_particles_tmp.xml", root_path, frame), "w" );
		fprintf(particle_fp,"<scene version=\"0.4.0\">\n");
		FOR_EACH_PARTICLE(particles) {
			const char *colors[] = { "blue", "green", "red", "yellow", "orange", "purple", "violet"  };
			uint color_index = 0;
			if( p.n <= 1 ) color_index = 0;
			else if( p.n <= 8 ) color_index = 1;
			else if( p.n <= 64 ) color_index = 2;
			else if( p.n <= 512 ) color_index = 3;
			else if( p.n <= 4096 ) color_index = 4;
			else if( p.n <= 32768 ) color_index = 5;
			else color_index = 6;
			
			fprintf(particle_fp,"<shape type=\"sphere\">\n" );
			fprintf(particle_fp,"<point name=\"center\" x=\"%f\" y=\"%f\" z=\"%f\"/> <float name=\"radius\" value=\"%f\"/>\n", p.p[0], p.p[2], p.p[1], 0.4*dpx*p.r );
			fprintf(particle_fp,"<ref name=\"bsdf\" id=\"%s\"/>\n", colors[color_index] );
			fprintf(particle_fp,"</shape>\n");
		} END_FOR
		fprintf(particle_fp,"</scene>\n");
		fclose(particle_fp);
		run("mv %s/mitsuba/%d_particles_tmp.xml %s/mitsuba/%d_particles.xml", root_path, frame, root_path, frame );
	}
	
#if 0
	FILE *isolated_fp = fopen( format_str("%s/mitsuba/%d_isolated_tmp.xml", root_path, frame), "w" );
	fprintf(isolated_fp,"<scene version=\"0.4.0\">\n");
	FOR_EACH_PARTICLE(particles) {
		if( p.isolated ) {
			fprintf(particle_fp,"<shape type=\"sphere\">\n" );
			fprintf(particle_fp,"<point name=\"center\" x=\"%f\" y=\"%f\" z=\"%f\"/> <float name=\"radius\" value=\"%f\"/>\n", p.p[0], p.p[2], p.p[1], 0.4*dpx*p.r );
			fprintf(particle_fp,"<ref name=\"bsdf\" id=\"blue\"/>\n" );
			fprintf(particle_fp,"</shape>\n");
		}
	} END_FOR
	fprintf(isolated_fp,"</scene>\n");
	fclose(isolated_fp);
	
	run("mv %s/mitsuba/%d_isolated_tmp.xml %s/mitsuba/%d_isolated.xml", root_path, frame, root_path, frame );
#endif
}
