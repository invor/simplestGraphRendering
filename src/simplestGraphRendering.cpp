#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <array>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include<iomanip>
#include <memory>

typedef unsigned int uint;

namespace Math
{
	#define PI 3.141592653589793238462643383279502884197169399375105820f

	struct Vec3
	{
		Vec3() : x(0.0), y(0.0), z(0.0) {}
		Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

		float x, y, z;

		Vec3 operator* (const float a)
		{
			return Vec3(a*x,a*y,a*z);
		}

		Vec3 operator+ (const Vec3& v)
		{
			return Vec3(x+v.x,y+v.y,z+v.z);
		}

		float length() const
		{
			return std::sqrt(x*x+y*y+z*z);
		}
	};
	
	struct Mat4x4
	{
		Mat4x4() : data() {}
		Mat4x4(const std::array<float,16> data) : data(data) {}
	
		std::array<float,16> data;
	
		float& operator[] (int index)
	{
		return data[index];
	}
	
		Mat4x4 operator* (const Mat4x4& m)
	{
		Mat4x4 b(m);

		return Mat4x4(
			std::array<float,16>({	// first row
				data[0]*b[0] + data[4]*b[1] + data[8]*b[2] + data[12]*b[3],
				data[1]*b[0] + data[5]*b[1] + data[9]*b[2] + data[13]*b[3],
				data[2]*b[0] + data[6]*b[1] + data[10]*b[2] + data[14]*b[3],
				data[3]*b[0] + data[7]*b[1] + data[11]*b[2] + data[15]*b[3],
				// second row
				data[0]*b[4] + data[4]*b[5] + data[8]*b[6] + data[12]*b[7],
				data[1]*b[4] + data[5]*b[5] + data[9]*b[6] + data[13]*b[7],
				data[2]*b[4] + data[6]*b[5] + data[10]*b[6] + data[14]*b[7],
				data[3]*b[4] + data[7]*b[5] + data[11]*b[6] + data[15]*b[7],
				// third row
				data[0]*b[8] + data[4]*b[9] + data[8]*b[10] + data[12]*b[11],
				data[1]*b[8] + data[5]*b[9] + data[9]*b[10] + data[13]*b[11],
				data[2]*b[8] + data[6]*b[9] + data[10]*b[10] + data[14]*b[11],
				data[3]*b[8] + data[7]*b[9] + data[11]*b[10] + data[15]*b[11],
				// fourth row
				data[0]*b[12] + data[4]*b[13] + data[8]*b[14] + data[12]*b[15],
				data[1]*b[12] + data[5]*b[13] + data[9]*b[14] + data[13]*b[15],
				data[2]*b[12] + data[6]*b[13] + data[10]*b[14] + data[14]*b[15],
				data[3]*b[12] + data[7]*b[13] + data[11]*b[14] + data[15]*b[15],
			})
		);
	}
	
		Mat4x4 inverse()
		{
			float a00 = data[0], a01 = data[1], a02 = data[2], a03 = data[3],
		    a10 = data[4], a11 = data[5], a12 = data[6], a13 = data[7],
		    a20 = data[8], a21 = data[9], a22 = data[10], a23 = data[11],
		    a30 = data[12], a31 = data[13], a32 = data[14], a33 = data[15],

		    b00 = a00 * a11 - a01 * a10,
		    b01 = a00 * a12 - a02 * a10,
		    b02 = a00 * a13 - a03 * a10,
		    b03 = a01 * a12 - a02 * a11,
		    b04 = a01 * a13 - a03 * a11,
		    b05 = a02 * a13 - a03 * a12,
		    b06 = a20 * a31 - a21 * a30,
		    b07 = a20 * a32 - a22 * a30,
		    b08 = a20 * a33 - a23 * a30,
		    b09 = a21 * a32 - a22 * a31,
		    b10 = a21 * a33 - a23 * a31,
		    b11 = a22 * a33 - a23 * a32,

		    // Calculate the determinant
		    det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

		    //if (!det) {
		    //    return null;
		    //}
			// The inverse is used as multiplication factor
		    det = 1.0f / det;

		    return Mat4x4(std::array<float,16>({	(a11 * b11 - a12 * b10 + a13 * b09) * det,
						(a02 * b10 - a01 * b11 - a03 * b09) * det,
						(a31 * b05 - a32 * b04 + a33 * b03) * det,
						(a22 * b04 - a21 * b05 - a23 * b03) * det,
						(a12 * b08 - a10 * b11 - a13 * b07) * det,
						(a00 * b11 - a02 * b08 + a03 * b07) * det,
						(a32 * b02 - a30 * b05 - a33 * b01) * det,
						(a20 * b05 - a22 * b02 + a23 * b01) * det,
						(a10 * b10 - a11 * b08 + a13 * b06) * det,
						(a01 * b08 - a00 * b10 - a03 * b06) * det,
						(a30 * b04 - a31 * b02 + a33 * b00) * det,
						(a21 * b02 - a20 * b04 - a23 * b00) * det,
						(a11 * b07 - a10 * b09 - a12 * b06) * det,
						(a00 * b09 - a01 * b07 + a02 * b06) * det,
						(a31 * b01 - a30 * b03 - a32 * b00) * det,
						(a20 * b03 - a21 * b01 + a22 * b00) * det}));
		}
	};

	Vec3 operator* (const float a, const Vec3& v)
	{
		return Vec3(a*v.x,a*v.y,a*v.z);
	}

	Vec3 operator/ (const Vec3& v,const float a)
	{
		return Vec3(v.x/a,v.y/a,v.z/a);
	}

	float dot(const Vec3& u, const Vec3& v)
	{
		return (u.x*v.x + u.y*v.y + u.z*v.z);
	}

	double dot64(const Vec3& u, const Vec3& v)
	{
		double ux = (double)u.x;
		double uy = (double)u.y;
		double uz = (double)u.z;
		double vx = (double)v.x;
		double vy = (double)v.y;
		double vz = (double)v.z;

		return (ux*vx+uy*vy+uz*vz);
	}

	Vec3 cross(const Vec3& u, const Vec3& v)
	{
		return Vec3(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
	}

	Vec3 normalize(const Vec3& v)
	{
		float l = v.length();

		return ( v/l );
	}
}

/**
 * Node struct for parsing preprocessed osm graph from file
 */
struct Node
{
	Node() : lat(0), lon(0) {}
	Node(double la, double lo) :
	lat(la), lon(lo) {}

	double lat;
	double lon;
};

/**
 * Edge struct for parsing preprocessed osm graph from file
 */
struct Edge
{
	Edge() : source(0), target(0), width(0), color(0) {}
	Edge(uint s, uint t, uint w, uint c) :
	source(s), target(t), width(w), color(c) {}

	uint source;
	uint target;
	uint width;
	int color;
};

/*
 * Each vertex contains the geo coordinates of a node and the color properties of adjacent edges.
 * Thus, for each node with adjacent egdes of different color multiple vertices are constructed.
 */
struct Vertex
{
	Vertex() : longitude(0.0), latitude(0.0), color(-1) {}
	Vertex(float lon, float lat) : longitude(lon), latitude(lat), color(-1) {}

	float longitude;
	float latitude;
	float color;
};

/*
 * 
 */
struct GeoBoundingBox
{
	float min_longitude;
	float min_latitude;

	float max_longitude;
	float max_latitude;
};

/**
 * Function to simply read the string of a shader source file from disk
 */
const std::string readShaderFile(const char* const path)
{
	std::ifstream inFile( path, std::ios::in );

	std::ostringstream source;
	while( inFile.good() ) {
		int c = inFile.get();
		if( ! inFile.eof() ) source << (char) c;
	}
	inFile.close();

	return source.str();
}

/**
 * Function for compiling shader source code. Returns the handle of the compiled shader
 */
GLuint compileShader(const std::string * const source, GLenum shaderType)
{
	/* Check if the source is empty */
	if (source->empty())
	{
		//TODO exit program?
	}

	/* Create shader object */
	const GLchar* c_source = source->c_str();
	GLuint shader = glCreateShader(shaderType);
	glShaderSource(shader, 1, &c_source, NULL);

	/* Compile shader */
	glCompileShader(shader);

	/* Check for errors */
	std::string shaderlog;
	GLint compile_ok = GL_FALSE;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_ok);

	GLint logLen = 0;
	shaderlog = "";
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
	if(logLen > 0)
	{
		char *log = new char[logLen];
		GLsizei written;
		glGetShaderInfoLog(shader, logLen, &written, log);
		shaderlog = log;
		delete [] log;

		std::cout<<shaderlog;
	}

	if(compile_ok == GL_FALSE)
	{
		glDeleteShader(shader);
		return -1;
	}

	return shader;
}

/**
 * Load a shader program
 * \attribute vs_path Path to vertex shader source file
 * \attribute fs_path Path to fragement shader source file
 * \attribute attributes Vertex shader input attributes (i.e. vertex layout)
 * \return Returns the handle of the created GLSL program
 */
GLuint createShaderProgram(const char* const vs_path, const char* const fs_path, std::vector<const char* const> attributes)
{
	/* Create a shader program object */
	GLuint handle;
	handle = glCreateProgram();

	/* Set the location (i.e. index) of the attribute (basically the input variable) in the vertex shader.
	 * The vertices intended to be used with this program will have to match that index in their
	 * attribute decription, so that a connection between the vertex data and the shader input can be made.
	 */
	for(size_t i=0; i<attributes.size(); i++)
		glBindAttribLocation(handle, i, attributes[i]);

	/* Read the shader source files */
	std::string vs_source = readShaderFile(vs_path);

	GLuint vertex_shader = compileShader(&vs_source, GL_VERTEX_SHADER);

	/* Attach shader to program */
	glAttachShader(handle, vertex_shader);

	/* Flag shader program for deletion.
	 * It will only be actually deleted after the program is deleted. (See destructor for program deletion.
	 */
	glDeleteShader(vertex_shader);


	/* Load, compile and attach fragment shader */
	std::string fs_source = readShaderFile(fs_path);

	GLuint fragment_shader = compileShader(&fs_source,GL_FRAGMENT_SHADER);

	/* Attach shader to program */
	glAttachShader(handle, fragment_shader);

	/* Flag shader program for deletion.
	 * It will only be actually deleted after the program is deleted. (See destructor for program deletion.
	 */
	glDeleteShader(fragment_shader);


	/* Link program */
	glLinkProgram(handle);

	/* Check if linking was successful */
	std::string shader_log;
	GLint status = GL_FALSE;
	glGetProgramiv(handle, GL_LINK_STATUS, &status);

	GLint logLen = 0;
	shader_log = "";
	glGetProgramiv(handle, GL_INFO_LOG_LENGTH, &logLen);
	if(logLen > 0)
	{
		char *log = new char[logLen];
		GLsizei written;
		glGetProgramInfoLog(handle, logLen, &written, log);
		shader_log = log;
		delete [] log;
	}

	std::cout<<shader_log<<std::endl;

	if(status == GL_FALSE)
		return -1;

	return handle;

}

/**
 * Collection of functions for loading graphic resources
 */
namespace ResourceLoader
{

	/**
	 * \brief Read a the header of a ppm image file. Courtesy to the computer vision lecture I attended.
	 * \param filename Location of the image file
	 * \param headerEndPos Out parameter, marks the point where the header of the ppm file ends
	 * \param imgDimX Out parameter, containing the dimension of the image in X direction in pixels
	 * \param imgDimY Out parameter, containing the dimension of the image in Y direction in pixels
	 * \return Returns true if the ppm header was succesfully read, false otherwise
	 */
	bool readPpmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY)
	{
		int currentComponent = 0;
		bool firstline = false;
		std::string::iterator itr1;
		std::string::iterator itr2;
		std::string buffer;
		std::string compBuffer;
		std::ifstream file (filename,std::ios::in | std::ios::binary);
	
		/*
		/	Check if the file could be opened.
		*/
		if(!( file.is_open() ))return false;
	
		/*
		/	Go to the beginning of the file and read the first line.
		*/
		file.seekg(0, file.beg);
		std::getline(file,buffer,'\n');
		itr1 = buffer.begin();
		for(itr2 = buffer.begin(); itr2 != buffer.end(); itr2++)
		{
			/*
			/	Check if the first line contains more than just ppm's magic number.
			/	If it does, it should look like this:
			/	"magic_number image_dimension_x image_dimension_y maximum_value"
			/	Therefore we scan the string for a space character and start parsing it.
			*/
			if(*itr2 == ' ')
			{
				if(currentComponent == 0)
				{
					/*	The first component is the magic number. We don't need it.	*/
					currentComponent++;
					firstline = true;
					itr1 = (itr2 + 1);
				}
				else if(currentComponent == 1)
				{
					/*	Get the image dimension in x.	*/
					compBuffer.assign(itr1, itr2);
					imgDimX = atoi(compBuffer.c_str());
					currentComponent++;
					itr1 = (itr2 + 1);
				}
				else if(currentComponent == 2)
				{
					/*	Get the image dimension in y.	*/
					compBuffer.assign(itr1, itr2);
					imgDimY = atoi(compBuffer.c_str());
					currentComponent++;
					itr1 = (itr2 + 1);
				}
			}
		}
	
		/*
		/	If the information we were looking for was inside the first line, we are done here.
		/	Note the position where we left off and exit with return true after closing the file.
		*/
		if(firstline)
		{
			headerEndPos = static_cast<long>(file.tellg());
			file.close();
			return true;
		}
	
		/*
		/	If the information wasn't inside the first line we have to keep reading lines.
		/	Skip all comment lines (first character = '#').
		*/
		std::getline(file,buffer,'\n');
		while( buffer[0]=='#' || (buffer.size() < 1) )
		{
			std::getline(file,buffer,'\n');
		}
	
		/*
		/	Now we should have a string containing the image dimensions and can extract them.
		*/
		itr1 = buffer.begin();
		for(itr2 = buffer.begin(); itr2 != buffer.end(); itr2++)
		{
			/*	Get the image dimension in x.	*/
			if(*itr2 == ' ')
			{
				compBuffer.assign(itr1, itr2);
				imgDimX = atoi(compBuffer.c_str());
				currentComponent++;
				itr1 = (itr2 + 1);
			}
		}
	
		/*
		/	The last component of a line can't be parsed within the loop since it isn't followed by
		/	a space character, but an end-of-line.
		/
		/	Get the image dimension in x.
		*/
		compBuffer.assign(itr1, itr2);
		imgDimY = atoi(compBuffer.c_str());
	
		/*
		/	Read one more line. This should contain the maximum value of the image, but we don't need
		/	that.
		/	Note down the position after this line and exit with return true after closing the file.
		*/
		std::getline(file,buffer,'\n');
		headerEndPos = static_cast<unsigned long>(file.tellg());
		file.close();
		return true;
	}
	
	/**
	 * \brief Read a the data of a ppm image file. Courtesy to the computer vision lecture I attended.
	 * \param filename Location of the image file
	 * \param imageData Pointer to the data buffer, that the image data will be written to
	 * \param dataBegin Marks the location within the ppm file, where the data block begins
	 * \param imgDimX Dimension of the image in X direction in pixels
	 * \param imgDimY Dimension of the image in Y direction in pixels
	 * \return Returns true if the ppm header was succesfully read, false otherwise
	 */
	bool readPpmData(const char* filename, char* imageData, unsigned long dataBegin, int imgDimX, int imgDimY)
	{
		std::ifstream file (filename,std::ios::in | std::ios::binary);
	
		/*
		/	Check if the file could be opened.
		*/
		if(!( file.is_open() ))return false;
	
		/*
		/	Determine the length from the beginning of the image data to the end of the file.
		*/
		file.seekg(0, file.end);
		unsigned long length = static_cast<unsigned long>(file.tellg());
		length = length - dataBegin;
		char* buffer = new char[length];
	
		file.seekg(dataBegin,std::ios::beg);
		file.read(buffer,length);
	
		/*
		/	Rearrange the image information so that the data begins with the lower left corner.
		*/
		int k = 0;
		for(int i=0; i < imgDimY; i++)
		{
			int dataLoc = (imgDimY-1-i)*imgDimX*3;
			for(int j=0; j < imgDimX; j++)
			{
				imageData[k]=buffer[dataLoc+(j*3)];
				k++;
				imageData[k]=buffer[dataLoc+(j*3)+1];
				k++;
				imageData[k]=buffer[dataLoc+(j*3)+2];
				k++;
			}
		}
	
		file.close();
		delete[] buffer;
		return true;
	}

}

/*
 * Camera (for OpenGL) orbiting a sphere that is centered on the origin
 */
struct OrbitalCamera
{
	OrbitalCamera() : longitude(0.0f), latitude(0.0f), orbit(5.0f), near(0.1f), far(10.0f), fovy(1.0), aspect_ratio(1.0) {}
	~OrbitalCamera() {}

	float longitude;
	float latitude;
	float orbit;

	float near;
	float far;
	float fovy;
	float aspect_ratio;

	void moveInOrbit(float delta_lat, float delta_lon, float delta_height)
	{
		latitude += delta_lat;
		longitude += delta_lon;

		if(longitude < -180.0f)
			longitude = 360.0f + longitude;

		if(longitude > 180.0f)
			longitude = -180.0f + (longitude-180.0f);

		//m_orbit += delta_height;

		orbit =  std::max(1.0f+0.000025f,orbit + delta_height);

		updateViewMatrix();
	}

	Math::Mat4x4 view_matrix;
	Math::Mat4x4 projection_matrix;

	void updateProjectionMatrix()
	{
		float f = 1.0f / std::tan(fovy / 2.0f);
        float nf = 1.0f / (near - far);
		projection_matrix[0] = f / aspect_ratio;
		projection_matrix[1] = 0.0f;
		projection_matrix[2] = 0.0f;
		projection_matrix[3] = 0.0f;
		projection_matrix[4] = 0.0f;
		projection_matrix[5] = f;
		projection_matrix[6] = 0.0f;
		projection_matrix[7] = 0.0f;
		projection_matrix[8] = 0.0f;
		projection_matrix[9] = 0.0f;
		projection_matrix[10] = (far + near) * nf;
		projection_matrix[11] = -1.0f;
		projection_matrix[12] = 0.0f;
		projection_matrix[13] = 0.0f;
		projection_matrix[14] = (2.0f * far * near) * nf;
		projection_matrix[15] = 0.0f;
	}

	void updateViewMatrix()
	{
		float lat_sin = sin( (PI/180.0f) * latitude);
		float lon_sin = sin( (PI/180.0f) * longitude);

		float lat_cos = cos( (PI/180.0f) * latitude);
		float lon_cos = cos( (PI/180.0f) * longitude);

		float camera_position[3];
		camera_position[0] = (lon_sin * lat_cos * orbit);
		camera_position[1] = (lat_sin * orbit);
		camera_position[2] = (lat_cos * lon_cos * orbit);

		Math::Mat4x4 lat_rotation(std::array<float,16>({1.0, 0.0, 0.0, 0.0,
								0.0, lat_cos, -lat_sin, 0.0,
								0.0, lat_sin, lat_cos, 0.0,
								0.0, 0.0, 0.0, 1.0}));

		Math::Mat4x4 lon_rotation(std::array<float,16>({lon_cos, 0.0, -lon_sin, 0.0,
								0.0, 1.0, 0.0, 0.0,
								lon_sin, 0.0, lon_cos, 0.0,
								0.0, 0.0, 0.0 , 1.0}));

		Math::Mat4x4 rotation_matrix = lon_rotation * lat_rotation;

		rotation_matrix = rotation_matrix.inverse();

		Math::Mat4x4 translation_matrix(std::array<float,16>({1.0, 0.0, 0.0, 0.0,
										0.0, 1.0, 0.0, 0.0,
										0.0, 0.0, 1.0, 0.0,
										-camera_position[0], -camera_position[1], -camera_position[2], 1.0}));

		view_matrix = rotation_matrix * translation_matrix;
	}

	GeoBoundingBox computeVisibleArea()
	{
		GeoBoundingBox bbox;

		// computer euklid. position of camera
		float lat_sin = sin( (PI/180.0f) * latitude);
		float lon_sin = sin( (PI/180.0f) * longitude);
	
		float lat_cos = cos( (PI/180.0f) * latitude);
		float lon_cos = cos( (PI/180.0f) * longitude);
	
		float r = orbit;
	
		Math::Vec3 world_position( lon_sin * lat_cos * r,
									lat_sin * r,
									lat_cos * lon_cos * r );

		// Build camera vector space
		Math::Vec3 front_vec =  normalize(-1.0 * world_position);
		Math::Vec3 right_vec = normalize(cross(front_vec, Math::Vec3(0.0f,1.0f,0.0f) ));
		Math::Vec3 up_vec = normalize(cross(right_vec,front_vec));

		// generate rays from the camera position along the right and upper frustrum plane
		Math::Vec3 upper_ray = normalize(front_vec + (std::tan(fovy*0.5f) * up_vec));
		Math::Vec3 right_ray = normalize(front_vec + (std::tan(fovy*0.5f)*aspect_ratio * right_vec));

		// intersect with unit sphere

		Math::Vec3 upper_intersection;
		float a = dot(upper_ray, upper_ray);
		float b = 2.0f * dot(upper_ray, world_position);
		float c = dot(world_position, world_position) - 1.0f * 1.0f;
		float discr = b * b - 4.0f * a * c;
		if (discr < 0.0f) {
			// no intersection, use tangent point from camera to sphere in up direction
			upper_intersection = up_vec;
		}
		else
		{
			float lambda = (-b - sqrt(discr)) / (2.0f * a);
			upper_intersection = world_position + (lambda * upper_ray); 
		}

		Math::Vec3 right_intersection;
		a = dot(right_ray, right_ray);
		b = 2.0f * dot(right_ray, world_position);
		c = dot(world_position, world_position) - 1.0f * 1.0f;
		discr = b * b - 4.0f * a * c;
		if (discr < 0.0f) {
			// no intersection, use tangent point from camera to sphere in right direction
			right_intersection = right_vec;
		}
		else
		{
			float lambda = (-b - sqrt(discr)) / (2.0f * a);
			right_intersection = world_position + (lambda * right_ray);
		}

		// using dot product suffers from nasty precision problems...
		// compute angle between Vec(SphereCenter,Camera) and Vec(SphereCenter,RightIntersection)
		//double longitude_angle = acos(dot64(normalize(world_position),normalize(right_intersection)));

		// compute angle between Vec(SphereCenter,Camera) and Vec(SphereCenter,UpperIntersection)
		//double latitude_angle = acos(dot64(normalize(world_position),normalize(upper_intersection)));


		// convert from euklidean to geo coordinates
		float upper_intersection_lat = (180.0f * std::asin(upper_intersection.y))/PI;

		float right_intersection_lat_cos = std::cos(std::asin(right_intersection.y));
		float right_intersection_lon = (180.0f * std::asin(right_intersection.x/right_intersection_lat_cos))/PI;

		if(dot(right_intersection,Math::Vec3(0.0,0.0,-1.0)) > 0.0)
		{
			if(dot(right_intersection,Math::Vec3(1.0,0.0,0.0)) > 0.0)
				right_intersection_lon = 180.0f - right_intersection_lon;
			else
				right_intersection_lon = -180.0f - right_intersection_lon;
		}

		bbox.min_longitude = longitude - abs(right_intersection_lon - longitude);
		bbox.min_longitude = (bbox.min_longitude < -180.0f) ? 360.0f + bbox.min_longitude : bbox.min_longitude;

		bbox.max_longitude = right_intersection_lon;

		bbox.min_latitude = latitude - abs(upper_intersection_lat -latitude);
		bbox.max_latitude = upper_intersection_lat;

		return bbox;
	}
};

/**
 * This struct essentially holds a renderable representation of a subgraph as a mesh, which is made up from
 * a set of vertices and a set of indices (the latter describing the mesh connectivity).
 * In this case the mesh uses line primitives, i.e. two succesive indices describe a single line segment.
 *
 * From a programming point of view, it makes sense to keep the three OpenGL handles required for a mesh obejct
 * organised together in a struct as most high level operations like "send the mesh data to the GPU" require
 * several OpenGL function calls and all of these handles.
 */
struct Subgraph
{
	Subgraph() : va_handle(0), vbo_handle(0), ibo_handle(0), isVisible(true), index_offsets(), line_widths() {}
	Subgraph(const Subgraph&) = delete;
	~Subgraph()
	{
		if( va_handle != 0 )
		{
			// delete mesh resources
			glBindVertexArray(va_handle);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			glDeleteBuffers(1, &ibo_handle);
			glBindBuffer(GL_ARRAY_BUFFER,0);
			glDeleteBuffers(1, &vbo_handle);
			glBindVertexArray(0);
			glDeleteVertexArrays(1, &va_handle);
		}
	}

	/* Handle for the vertex array object */
	GLuint va_handle;

	/* Handle for the vertex buffer object (allows access to vertex data in GPU memory) */
	GLuint vbo_handle;

	/* Handle for the index buffer objects (allows access to index data in GPU memory) */
	GLuint ibo_handle;

	bool isVisible;

	/* To draw lines of different type, i.e of different width seperatly but still store them
	 * in the same index buffer object, offsets into the buffer are used to only draw a subset of the index buffer
	 * in each draw call.
	 */
	std::vector<uint> index_offsets;
	/* Stores the width of each subset of lines */
	std::vector<float> line_widths;


	void loadGraphData(std::vector<Node>& nodes, std::vector<Edge>& edges)
	{
		std::vector<Vertex> vertices;
		std::vector<uint> indices;

		// At least as many vertices as there are nodes are required
		vertices.reserve(nodes.size());
	
		// Each edge contributes two indices
		indices.reserve(edges.size()*2);
	
		// Copy geo coordinates from input nodes to vertices
		for(auto& node : nodes)
		{
			vertices.push_back(Vertex(node.lon,node.lat));
		}
	
		std::sort(edges.begin(),edges.end(), [](Edge u, Edge v) { return u.width < v.width; } );
	
		// Copy indices from edge array to index array
		std::vector<bool> has_next(nodes.size(), false);
		std::vector<uint> next(nodes.size(),0);
		uint width = 0;
		uint counter = 0;
		for(auto& edge : edges)
		{
			uint src_id = edge.source;
			uint tgt_id = edge.target;
	
			while(has_next[src_id] && (vertices[src_id].color != edge.color))
			{
				src_id = next[src_id];
			}
	
			if(vertices[src_id].color == -1)
			{
				vertices[src_id].color = edge.color;
			}
	
			if(vertices[src_id].color != edge.color)
			{
				uint next_id = vertices.size();
				vertices.push_back(Vertex(vertices[src_id].longitude,vertices[src_id].latitude));
				vertices[next_id].color = edge.color;
				has_next.push_back(false);
				next.push_back(0);
	
				next[src_id] = next_id;
				has_next[src_id] = true;
				src_id = next_id;
			}
	
			while(has_next[tgt_id] && (vertices[tgt_id].color != edge.color))
			{
				tgt_id = next[tgt_id];
			}
	
			if(vertices[tgt_id].color == -1)
			{
				vertices[tgt_id].color = edge.color;
			}
	
			if(vertices[tgt_id].color != edge.color)
			{
				uint next_id = (uint)vertices.size();
				vertices.push_back(Vertex(vertices[tgt_id].longitude,vertices[tgt_id].latitude));
				vertices[next_id].color = edge.color;
				has_next.push_back(false);
				next.push_back(0);
	
				next[tgt_id] = next_id;
				has_next[tgt_id] = true;
				tgt_id = next_id;
			}
	
			//std::cout<<"Edge color: "<<edge.color<<std::endl;
			//std::cout<<"Source color: "<<vertices[src_id].color<<std::endl;
			//std::cout<<"Target color: "<<vertices[tgt_id].color<<std::endl;
	
			if(width != edge.width)
			{
				index_offsets.push_back(counter);
				line_widths.push_back((float)edge.width);
				width = edge.width;
			}
	
			indices.push_back(src_id);
			indices.push_back(tgt_id);
	
			counter += 2;
		}
		index_offsets.push_back(indices.size());


		// Allocate GPU memory and send data
		if(vertices.size() < 1 || indices.size() < 1)
			return;

		auto va_size = sizeof(Vertex) * vertices.size();
		auto vi_size = sizeof(uint) * indices.size();

		if(va_handle == 0 || vbo_handle == 0 || ibo_handle == 0)
		{
			glGenVertexArrays(1, &va_handle);
			glGenBuffers(1, &vbo_handle);
			glGenBuffers(1, &ibo_handle);
		}

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glBufferData(GL_ARRAY_BUFFER, va_size, vertices.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_handle);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, vi_size, indices.data(), GL_STATIC_DRAW);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(Vertex), 0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 1, GL_FLOAT, false, sizeof(Vertex), (GLvoid*) (sizeof(GL_FLOAT)*2));
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	
		std::cout<<"GfxGraph consisting of "<<vertices.size()<<" vertices and "<<indices.size()<<" indices"<<std::endl;
	}

	void draw(float scale)
	{
		//glBindVertexArray(va_handle);
		//glDrawElements(GL_LINES, indices.size(),  GL_UNSIGNED_INT,  0 );

		for(size_t i=0; i< index_offsets.size()-1; i++)
		{
			glLineWidth(std::max(1.0f,line_widths[i] * scale));
			//glLineWidth(line_widths[i]);

			glBindVertexArray(va_handle);
			glDrawElements(GL_LINES,  index_offsets[i+1]-index_offsets[i],  GL_UNSIGNED_INT,  (void*)(index_offsets[i] * sizeof(GLuint)) );
		}
	}
};

/**
 * A graph made up from subgraphs. Limited to rendering edges.
 * This struct primarily holds a set of subgraphs and offers the neccessary functionality to add and change subgraphs.
 */
struct Graph
{
	std::vector<std::unique_ptr<Subgraph>> subgraphs;

	void addSubgraph(std::vector<Node>& nodes, std::vector<Edge>& edges)
	{
		std::unique_ptr<Subgraph> subgraph(new Subgraph);
		subgraphs.push_back(std::move(subgraph));

		subgraphs.back()->loadGraphData(nodes,edges);
	}

	void setVisibilty(uint index, bool visibility)
	{
		if(index < subgraphs.size())
			subgraphs[index]->isVisible = visibility;
	}

	void draw(float scale)
	{
		for(auto& subgraph : subgraphs)
		{
			if(subgraph->isVisible)
				subgraph->draw(scale);
		}
	}
};

struct TextLabels
{
	std::array<float,255> u;
	std::array<float,255> v;

	TextLabels()
	{
		std::vector<std::string> atlas_rows;

		atlas_rows.push_back("ABCDEFGHIJKLMN");
		atlas_rows.push_back("OPQRSTUVWXYZab");
		atlas_rows.push_back("cdefghijklmnop");
		atlas_rows.push_back("qrstuvwxyz1234");
		atlas_rows.push_back("567890&@.,?!'\"");
		atlas_rows.push_back("\"()*ßöäü-_");


		float u_value = 1.0/16.0;
		float v_value = 5.0/6.0;

		for(std::string& s : atlas_rows)
		{
			for(char c : s)
			{
				u[c] = u_value;
				u_value += 1.0/16.0;
				v[c] = v_value;
			}

			u_value = 1.0f/16.0f;
			v_value -= 1.0f/6.0f;
		}

		// Create basic quad mesh for rendering characters
		std::array< float, 16 > vertex_array = {{ -0.03,-0.1,0.0,0.0,
												-0.03,0.1,0.0,1.0,
												0.03,0.1,1.0,1.0,
												0.03,-0.1,1.0,0.0 }};

		std::array< GLuint, 6 > index_array = {{ 0,1,2,2,0,3 }};

		glGenVertexArrays(1, &va_handle);
		glGenBuffers(1, &vbo_handle);
		glGenBuffers(1, &ibo_handle);

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_array), vertex_array.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_handle);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(index_array), index_array.data(), GL_STATIC_DRAW);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(GL_FLOAT)*4, 0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, false, sizeof(GL_FLOAT)*4, (GLvoid*) (sizeof(GL_FLOAT)*2));
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// Load text label shader program
		prgm_handle = createShaderProgram("../src/textLabel_v.glsl","../src/textLabel_f.glsl",{"v_position","v_uv"});

		// Load font atlas
		unsigned long begin_pos;
		int x_dim, y_dim;
		char* img_data;
		ResourceLoader::readPpmHeader("../resources/font_atlas.ppm", begin_pos, x_dim, y_dim);
		img_data = new char[x_dim * y_dim * 3];
		ResourceLoader::readPpmData("../resources/font_atlas.ppm", img_data, begin_pos, x_dim, y_dim);

		glGenTextures(1, &font_atlas_handle);
		//assert(font_atlas_handle > 0);
		glBindTexture(GL_TEXTURE_2D, font_atlas_handle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, x_dim, y_dim, 0, GL_RGB, GL_UNSIGNED_BYTE, img_data);
		glBindTexture(GL_TEXTURE_2D,0);

		delete[] img_data;

	}
	/* Structs and classes that free OpenGL resources within their destructor
	 * should be non-copyable to prevent bad stuff from happening
	 */
	TextLabels(const TextLabels&) = delete;
	~TextLabels()
	{
		// delete mesh resources
		glBindVertexArray(va_handle);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glDeleteBuffers(1, &ibo_handle);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glDeleteBuffers(1, &vbo_handle);
		glBindVertexArray(0);
		glDeleteVertexArrays(1, &va_handle);

		// delete GLSL program
		glDeleteProgram(prgm_handle);

		// delete font atlas
		glDeleteTextures(1,&font_atlas_handle);

		for(GLuint tx_handle : text_texture_handles)
			glDeleteTextures(1,&tx_handle);
	}

	GLuint va_handle;

	GLuint vbo_handle;

	GLuint ibo_handle;

	GLuint prgm_handle;

	GLuint font_atlas_handle;

	/** Geo-Cooridnates of each label */
	std::vector<float> geoCoordinates;
	/** String of of each label encoded as look-up-texture for font altas */
	std::vector<GLuint> text_texture_handles;
	/** Number of characters of each label */
	std::vector<unsigned int> lengths;
	/** Relative scale of each label */
	std::vector<float> scales;
	/** Visibility of each label i.e. rendered or not */
	std::vector<bool> visibility;

	void addLabel(std::string label_text, float latitude, float longitude, float scale)
	{
		geoCoordinates.push_back(longitude);
		geoCoordinates.push_back(latitude);

		// convert string to text texture (i.e. character to uv position in texture atlas)
		float* data = new float[label_text.length()*2];
		for(size_t i=0; i<label_text.length(); i++)
		{
			data[i*2] = u[(int)label_text[i]];
			data[i*2 + 1] = v[(int)label_text[i]];
		}

		GLuint texture_handle = 0;
		glGenTextures(1, &texture_handle);
		//assert(font_atlas_handle > 0);
		glBindTexture(GL_TEXTURE_2D, texture_handle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, label_text.length(), 1, 0, GL_RG, GL_FLOAT, data);
		glBindTexture(GL_TEXTURE_2D,0);

		text_texture_handles.push_back(texture_handle);

		lengths.push_back(label_text.length());
		scales.push_back(scale);
		visibility.push_back(true);
	}

	void draw(OrbitalCamera& camera)
	{
		glUseProgram(prgm_handle);

		// bind font atlas texture
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D,font_atlas_handle);
		int zero = 0;
		glUniform1iv(glGetUniformLocation(prgm_handle,"fontAtlas_tx2D"),1,&zero);

		for(size_t i=0; i<visibility.size(); i++)
		{
			if(visibility[i])
			{
				// set uniforms
				glUniformMatrix4fv(glGetUniformLocation(prgm_handle, "view_matrix"), 1, GL_FALSE, camera.view_matrix.data.data());
				glUniformMatrix4fv(glGetUniformLocation(prgm_handle, "projection_matrix"), 1, GL_FALSE, camera.projection_matrix.data.data());

				glUniform2fv(glGetUniformLocation(prgm_handle,"label_geoCoords"),1,&geoCoordinates[i*2]);
				float charCount = (float)lengths[i];
				glUniform1fv(glGetUniformLocation(prgm_handle,"label_charCount"),1,&charCount);
				glUniform1fv(glGetUniformLocation(prgm_handle,"label_scale"),1,&scales[i]);

				// bind content texture
				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D,text_texture_handles[i]);
				int one = 1;
				glUniform1iv(glGetUniformLocation(prgm_handle,"label_text_tx2D"),1,&one);

				glBindVertexArray(va_handle);
				glDrawElementsInstanced(GL_TRIANGLES,  6,  GL_UNSIGNED_INT,  (void*)NULL, lengths[i]);
			}
		}
	}


};

struct DebugSphere
{
	DebugSphere() : va_handle(0), vbo_handle(0), ibo_handle(0) {}
	DebugSphere(const DebugSphere&) = delete;
	~DebugSphere()
	{
		// delete mesh resources
		glBindVertexArray(va_handle);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glDeleteBuffers(1, &ibo_handle);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glDeleteBuffers(1, &vbo_handle);
		glBindVertexArray(0);
		glDeleteVertexArrays(1, &va_handle);
	}

	/* Handle for the vertex array object */
	GLuint va_handle;

	/* Handle for the vertex buffer object (allows access to vertex data in GPU memory) */
	GLuint vbo_handle;

	/* Handle for the index buffer objects (allows access to index data in GPU memory) */
	GLuint ibo_handle;

	void create()
	{
		std::vector<float> vertices;
		std::vector<uint> indices;

		float latitude = -90.0;
		float longitude = -180.0;
		for(int i=0; i<50; i=i+1)
		{
			for(int j=0; j<100; j=j+1)
			{
				float lat_sin = sin( (PI/180.0f) * latitude);
				float lon_sin = sin( (PI/180.0f) * longitude);

				float lat_cos = cos( (PI/180.0f) * latitude);
				float lon_cos = cos( (PI/180.0f) * longitude);

				float r = 1.0f; //6378137.0;

				vertices.push_back(lon_sin * lat_cos * r);
				vertices.push_back(lat_sin * r);
				vertices.push_back(lat_cos * lon_cos * r);

				longitude += (1.0f/100.0f) * 360.0f;
			}
			latitude += (1/50.0f) * 180.0f;
		}

		for(int i=0; i<50*100; i=i+1)
		{
			indices.push_back((uint)i);
		}

		if(vertices.size() < 1 || indices.size() < 1)
			return;

		auto va_size = sizeof(float) * vertices.size();
		auto vi_size = sizeof(uint) * indices.size();

		if(va_handle == 0 || vbo_handle == 0 || ibo_handle == 0)
		{
			glGenVertexArrays(1, &va_handle);
			glGenBuffers(1, &vbo_handle);
			glGenBuffers(1, &ibo_handle);
		}

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glBufferData(GL_ARRAY_BUFFER, va_size, vertices.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_handle);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, vi_size, indices.data(), GL_STATIC_DRAW);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glBindVertexArray(va_handle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(GL_FLOAT)*3, 0);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void draw()
	{
		glBindVertexArray(va_handle);
		glDrawElements(GL_POINTS,  5000,  GL_UNSIGNED_INT,  NULL );
	}
};

void windowSizeCallback(GLFWwindow *window, int width, int height)
{
	OrbitalCamera* active_camera = reinterpret_cast<OrbitalCamera*>(glfwGetWindowUserPointer(window));

	active_camera->aspect_ratio =  ((float)width/(float)height);
	active_camera->updateProjectionMatrix();
}

namespace Controls {

	namespace
	{
		std::array<float,2> latest_cursor_position = {{0.0,0.0}};
	}

	void mouseScrollFeedback(GLFWwindow *window, double x_offset, double y_offset)
	{
		OrbitalCamera* active_camera = reinterpret_cast<OrbitalCamera*>(glfwGetWindowUserPointer(window));

		float camera_height_inertia = std::pow( (active_camera->orbit - 1.0f )*0.1f, 1.0);

		active_camera->moveInOrbit(0.0,0.0,-camera_height_inertia * (float)y_offset);
	}

	void updateOrbitalCamera(GLFWwindow *window)
	{
		OrbitalCamera* active_camera = reinterpret_cast<OrbitalCamera*>(glfwGetWindowUserPointer(window));

		/*	Get current cursor position on the screen */
		double pos_x, pos_y;
		glfwGetCursorPos(window, &pos_x, &pos_y);
		std::array<float,2> currentCursorPosition = {{(float) pos_x, (float) pos_y}};
		std::array<float,2> cursor_movement = {{latest_cursor_position[0] - currentCursorPosition[0],latest_cursor_position[1] - currentCursorPosition[1]}};

		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS)
		{
			//camera_lon += cursor_movement.x;
			//camera_lat -= cursor_movement.y;

			GLfloat camera_vertical_inertia = std::pow( (active_camera->orbit - 1.0)*0.1, 1.2f);

			active_camera->moveInOrbit(-cursor_movement[1]*camera_vertical_inertia,cursor_movement[0]*camera_vertical_inertia,0.0);
		}

		latest_cursor_position = currentCursorPosition;
	}

}

/*
 * Collection of functions for parsing graph files. Taken from lasagne and modified to meet expected graph input format
 */
namespace Parser
{
	/**
	 * Creates a new graph node
	 * @param input_string Input string containing node data
	 * @param r_node Node created from input string
	 */
	void createNode(std::string input_string, Node& r_node)
	{
		std::stringstream ss(input_string);
		std::string buffer;

		ss >> buffer;
		r_node.lat = atof(buffer.c_str());

		ss >> buffer;
		r_node.lon = atof(buffer.c_str());
	}

	/**
	 * createEdge - erstellt eine neue Kante
	 * @param input_string Input string containing edge data
	 * @param r_edge Edge created from input string
	 */
	void createEdge(std::string input_string, Edge& r_edge)
	{
		std::stringstream ss(input_string);
		std::string buffer;

		ss >> buffer;
		r_edge.source = (uint)atoi(buffer.c_str());

		ss >> buffer;
		r_edge.target = (uint)atoi(buffer.c_str());

		ss >> buffer;
		r_edge.width = (uint)atoi(buffer.c_str());

		ss >> buffer;
		r_edge.color = (int)atoi(buffer.c_str());
	}

	/**
	 * Parse node and egde from input file
	 * @param graphfile Path to the graphfile
	 * @param n Vector for storing the parsed nodes
	 * @param e Vector for storing the parsed edges
	 */
	bool parseTxtGraphFile(std::string graphfile, std::vector<Node>& n, std::vector<Edge>& e)
	{
		std::string buffer;
		std::ifstream file;

		file.open(graphfile.c_str(), std::ios::in);

		if( file.is_open())
		{
			file.seekg(0, std::ios::beg);
			getline(file,buffer,'\n');
			uint node_count = (uint)atoi(buffer.c_str());
			getline(file,buffer,'\n');
			uint edge_count = (uint)atoi(buffer.c_str());

			n.reserve(node_count);
			for(uint i=0; i<node_count; i++)
			{
				getline(file,buffer,'\n');
				n.push_back(Node());
				createNode(buffer, n.back() );
			}

			e.reserve(edge_count);
			for(uint j=0; j<edge_count; j++)
			{
				getline(file,buffer,'\n');
				e.push_back(Edge());
				createEdge(buffer, e.back() );
			}
			file.close();

			return true;
		}

		return false;
	}
}

int main(int argc, char*argv[])
{
	/////////////////////////////////////
	// overly simple command line parsing
	/////////////////////////////////////

	std::string filepath;

	int i=1;
    if (argc < 3) {
		std::cout<<"Supply a graph with -gf <graph.gl>"<<std::endl; return 0;
	}
	while(i<argc)
	{
		if(argv[i] == (std::string) "-gf")
		{
			i++;
			if(i<argc) { filepath = argv[i]; i++; }
			else { std::cout<<"Missing parameter for -gf"<<std::endl; return 0; }
		}
		else
		{
			i++;
		}
	}


	/////////////////////////////////////
	// Window and OpenGL Context creation
	/////////////////////////////////////

    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    window = glfwCreateWindow(1600, 900, "Simple Graph Renderer", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

	/*	Initialize glew */
	//glewExperimental = GL_TRUE;
	GLenum error = glewInit();
	if( GLEW_OK != error)
	{
		std::cout<<"-----\n"
				<<"The time is out of joint - O cursed spite,\n"
				<<"That ever I was born to set it right!\n"
				<<"-----\n"
				<<"Error: "<<glewGetErrorString(error);
		return false;
	}
	/* Apparently glewInit() causes a GL ERROR 1280, so let's just catch that... */
	glGetError();

	/* Intialize controls */
	glfwSetWindowSizeCallback(window,windowSizeCallback);
	glfwSetScrollCallback(window, Controls::mouseScrollFeedback);
	/* Hide cursor */
	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);


	////////////////////////////
	// Load graph data from file
	////////////////////////////

	std::vector<Node> nodes;
	std::vector<Edge> edges;

	Parser::parseTxtGraphFile(filepath,nodes,edges);


	/////////////////////////////////////////////////////////////////////
	// Creation of graphics resources, i.e. shader programs, meshes, etc.
	/////////////////////////////////////////////////////////////////////

	/*
	 * OpenGL Objects are usaully represented by handles.
	 * From my understanding, these are basically pointers/adresses to something
	 * the GPU (or GPU driver) owns.
	 *
	 * Instances of structs holding OpenGL handles are created within an additional
	 * scope, so that they are destroyed while the OpenGL context is still alive
	 */	
	{
		/* Create GLSL programs */
		GLuint shader_prgm_handle = createShaderProgram("../src/edge_v.glsl","../src/edge_f.glsl",{"v_geoCoords","v_color"});
		GLuint debug_prgm_handle = createShaderProgram("../src/debug_v.glsl","../src/debug_f.glsl",{"v_position"});

		//GLenum glerror = glGetError();
		//std::cout<<glerror<<std::endl;

		/* Create renderable graph (mesh) */
		Graph lineGraph;
		lineGraph.addSubgraph(nodes,edges);

		/* Create a orbital camera */
		OrbitalCamera camera;
		camera.longitude = 0.0f;
		camera.latitude = 0.0f;
		camera.orbit = 5.0f;
		camera.near = 0.0001f;
		camera.far = 10.0f;
		camera.fovy = 30.0f * PI/180.0f;
		camera.aspect_ratio = 16.0f/9.0f;


		/* Make camera accessable in window callbacks */
		glfwSetWindowUserPointer(window,&camera);

		/* Create text labels */
		TextLabels labels;

		//TODO FIND BUG: camera matrices have to be updated after label object creation....
		camera.updateViewMatrix();
		camera.updateProjectionMatrix();

		for(int lon=-180; lon<=180 ; lon++)
			labels.addLabel(std::to_string(lon),0.0,lon,0.25);

		for(int lat=-90; lat<=90 ; lat++)
			labels.addLabel(std::to_string(lat),lat,0.0,0.25);

		/* Create the debug sphere */
		DebugSphere db_sphere;
		db_sphere.create();


		//////////////////////////////////////////////////////
		// Render Loop - Each loop iteration renders one frame
		//////////////////////////////////////////////////////

		glEnable(GL_BLEND);
		//glDisable(GL_CULL_FACE);
		//glDisable(GL_DEPTH_TEST);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		/* Loop until the user closes the window */
		while (!glfwWindowShouldClose(window))
		{
			Controls::updateOrbitalCamera(window);

		    /* Render here */
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			int width, height;
			glfwGetFramebufferSize(window, &width, &height);
			glViewport(0, 0, width, height);

			glUseProgram(debug_prgm_handle);

			glUniformMatrix4fv(glGetUniformLocation(debug_prgm_handle, "view_matrix"), 1, GL_FALSE, camera.view_matrix.data.data());
			glUniformMatrix4fv(glGetUniformLocation(debug_prgm_handle, "projection_matrix"), 1, GL_FALSE, camera.projection_matrix.data.data());

			glPointSize(2.0);
			db_sphere.draw();


			glUseProgram(shader_prgm_handle);

			glUniformMatrix4fv(glGetUniformLocation(shader_prgm_handle, "view_matrix"), 1, GL_FALSE, camera.view_matrix.data.data());
			glUniformMatrix4fv(glGetUniformLocation(shader_prgm_handle, "projection_matrix"), 1, GL_FALSE, camera.projection_matrix.data.data());

			float scale = std::min((0.0025/(camera.orbit - 1.0f)),2.0);

			lineGraph.draw( scale );

			labels.draw(camera);

			GeoBoundingBox bbox = camera.computeVisibleArea();
			//std::cout << std::fixed;
			//std::cout << std::setprecision(20);
			//std::cout<<"min latitude: "<<bbox.min_latitude<<" max latitude: "<<bbox.max_latitude<<std::endl;
			//std::cout<<"min longitude: "<<bbox.min_longitude<<" max longitude: "<<bbox.max_longitude<<std::endl;
			//std::cout<<"camera lon: "<<camera.longitude<<" "<<camera.latitude<<std::endl;

		    /* Swap front and back buffers */
		    glfwSwapBuffers(window);

		    /* Poll for and process events */
		    glfwPollEvents();
		}


		// delete/free graphics resources
		glDeleteProgram(shader_prgm_handle);
		glDeleteProgram(debug_prgm_handle);
	}

    glfwTerminate();
    return 0;
}
