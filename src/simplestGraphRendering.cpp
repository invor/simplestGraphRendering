#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <array>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>

typedef unsigned int uint;

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
			{{	// first row
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
			}}
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
	    det = 1.0 / det;

	    return Mat4x4({{	(a11 * b11 - a12 * b10 + a13 * b09) * det,
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
					(a20 * b03 - a21 * b01 + a22 * b00) * det	}});
	}
};

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
	Edge() : source(0), target(0), width(0.0), color(0) {}
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
	Vertex(float lon, float lat) : longitude(lon), latitude(lat), color(-1) {}

	float longitude;
	float latitude;
	float color;
};

/*
 * Camera (for OpenGL) orbiting a sphere that is centered on the origin
 */
struct OrbitalCamera
{
	OrbitalCamera() {}
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
		//m_orbit += delta_height;

		orbit =  std::max(1.0f+0.000025f,orbit + delta_height);

		updateViewMatrix();
	}

	Mat4x4 view_matrix;
	Mat4x4 projection_matrix;

	void updateProjectionMatrix()
	{
		float f = 1.0f / std::tan(fovy / 2.0f);
        float nf = 1.0f / (near - far);
		projection_matrix[0] = f / aspect_ratio;
		projection_matrix[1] = 0;
		projection_matrix[2] = 0;
		projection_matrix[3] = 0;
		projection_matrix[4] = 0;
		projection_matrix[5] = f;
		projection_matrix[6] = 0;
		projection_matrix[7] = 0;
		projection_matrix[8] = 0;
		projection_matrix[9] = 0;
		projection_matrix[10] = (far + near) * nf;
		projection_matrix[11] = -1;
		projection_matrix[12] = 0;
		projection_matrix[13] = 0;
		projection_matrix[14] = (2 * far * near) * nf;
		projection_matrix[15] = 0;
	}

	void updateViewMatrix()
	{
		float PI = 3.141592653589793238462643383279502884197169399375105820;

		float lat_sin = sin( (PI/180.0f) * latitude);
		float lon_sin = sin( (PI/180.0f) * longitude);

		float lat_cos = cos( (PI/180.0f) * latitude);
		float lon_cos = cos( (PI/180.0f) * longitude);

		float camera_position[3];
		camera_position[0] = (lon_sin * lat_cos * orbit);
		camera_position[1] = (lat_sin * orbit);
		camera_position[2] = (lat_cos * lon_cos * orbit);

		Mat4x4 lat_rotation({{1.0, 0.0, 0.0, 0.0,
								0.0, lat_cos, -lat_sin, 0.0,
								0.0, lat_sin, lat_cos, 0.0,
								0.0, 0.0, 0.0, 1.0}});

		Mat4x4 lon_rotation({{lon_cos, 0.0, -lon_sin, 0.0,
								0.0, 1.0, 0.0, 0.0,
								lon_sin, 0.0, lon_cos, 0.0,
								0.0, 0.0, 0.0 , 1.0}});

		Mat4x4 rotation_matrix = lon_rotation * lat_rotation;

		rotation_matrix = rotation_matrix.inverse();

		Mat4x4 translation_matrix({{1.0, 0.0, 0.0, 0.0,
										0.0, 1.0, 0.0, 0.0,
										0.0, 0.0, 1.0, 0.0,
										-camera_position[0], -camera_position[1], -camera_position[2], 1.0}});

		view_matrix = rotation_matrix * translation_matrix;
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
	Subgraph() : va_handle(0), vbo_handle(0), ibo_handle(0), index_offsets(), line_widths() {}

	/* Handle for the vertex array object */
	GLuint va_handle;

	/* Handle for the vertex buffer object (allows access to vertex data in GPU memory) */
	GLuint vbo_handle;

	/* Handle for the index buffer objects (allows access to index data in GPU memory) */
	GLuint ibo_handle;

	/* To draw lines of different type, i.e of different width seperatly but still store them
	 * in the same index buffer object, offsets into the buffer are used to only draw a subset of the index buffer
	 * in each draw call.
	 */
	std::vector<uint> index_offsets;
	/* Stores the width of each subset of lines */
	std::vector<float> line_widths;

	/**
	 * Allocate GPU memory and send data
	 */
	void bufferGraphData(std::vector<Vertex>& vertices, std::vector<uint>& indices)
	{
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
struct GfxGraph
{
	GfxGraph() : vertices(), indices() {}

	/* CPU-side storage for vertex data. Used when converting graph data. */
	std::vector<Vertex> vertices;
	/* CPU-side storage for index data. Used when converting graph data. */
	std::vector<uint> indices;

	/* Visibility information for each subgraph */
	std::vector<bool> isVisible;

	/* The subgraphs that make up the graph itself */
	std::vector<Subgraph> subgraphs;

	/**
	 * Add a new subgraph to the graph
	 */
	void addSubgraph(std::vector<Node>& nodes, std::vector<Edge>& edges)
	{
		uint new_subgraph_index = subgraphs.size();

		isVisible.push_back(true);
		subgraphs.push_back(Subgraph());

		convertGraphData(new_subgraph_index,nodes,edges);

		subgraphs[new_subgraph_index].bufferGraphData(vertices,indices);
	}

	/**
	 * Update a specific subgraph with new data
	 */
	void reloadSubgraph(uint index, std::vector<Node>& nodes, std::vector<Edge>& edges)
	{
		if( index >= subgraphs.size() )
			return;

		subgraphs[index].index_offsets.clear();
		subgraphs[index].line_widths.clear();

		convertGraphData(index,nodes,edges);

		subgraphs[index].bufferGraphData(vertices,indices);
	}

	/**
	 * Set visibility for a specific subgraph
	 */
	void setVisibility(uint index, bool visibility)
	{
		if( index < isVisible.size() )
			isVisible[index] = visibility;
	}

	/**
	 * Draw all subgraphs set to visible
	 */
	void draw(float scale)
	{
		for(int i=0; i<subgraphs.size(); i++)
		{
			if(isVisible[i])
				subgraphs[i].draw(scale);
		}
	}

	private:

	/**
	* Converts a graph from the input format into the format that is send to the GPU for rendering, i.e. a mesh build from vertex and index data.
	* A purely CPU based preprocessing function that converts the data.
	* Note: If the updated graph sections are not in the same format that is used for rendering, this function has to be called each time the graph
	* or a section of the graph is updated, before the updated graph data can be send to the GPU for rendering!
	*/
	void convertGraphData(uint subgraph_index, std::vector<Node>& nodes, std::vector<Edge>& edges)
	{
		vertices.clear();
		indices.clear();

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
		uint width = 0.0;
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
				uint next_id = vertices.size();
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
				subgraphs[subgraph_index].index_offsets.push_back(counter);
				subgraphs[subgraph_index].line_widths.push_back((float)edge.width);
				width = edge.width;
			}

			indices.push_back(src_id);
			indices.push_back(tgt_id);

			counter += 2;
		}
		subgraphs[subgraph_index].index_offsets.push_back(indices.size());

		std::cout<<"GfxGraph consisting of "<<vertices.size()<<" vertices and "<<indices.size()<<" indices"<<std::endl;
	}

};

struct DebugSphere
{
	DebugSphere() : va_handle(0), vbo_handle(0), ibo_handle(0) {}

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
				float PI = 3.141592653589793238462643383279502884197169399375105820;

				float lat_sin = sin( (PI/180.0f) * latitude);
				float lon_sin = sin( (PI/180.0f) * longitude);

				float lat_cos = cos( (PI/180.0f) * latitude);
				float lon_cos = cos( (PI/180.0f) * longitude);

				float r = 1.0; //6378137.0;

				vertices.push_back(lon_sin * lat_cos * r);
				vertices.push_back(lat_sin * r);
				vertices.push_back(lat_cos * lon_cos * r);

				longitude += (1/100.0) * 360.0;
			}
			latitude += (1/50.0) * 180.0;
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
	for(int i=0; i<attributes.size(); i++)
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
	 */

	/* Create GLSL program */
	GLuint shader_prgm_handle = createShaderProgram("../src/edge_v.glsl","../src/edge_f.glsl",{"v_geoCoords","v_color"});
	GLuint debug_prgm_handle = createShaderProgram("../src/debug_v.glsl","../src/debug_f.glsl",{"v_position"});

	//GLenum glerror = glGetError();
	//std::cout<<glerror<<std::endl;

	/* Create renderable graph (mesh) */
	GfxGraph lineGraph;
	lineGraph.addSubgraph(nodes,edges);

	/* Create a orbital camera */
	OrbitalCamera camera;
	camera.longitude = 0.0;
	camera.latitude = 0.0;
	camera.orbit = 5.0;
	camera.near = 0.0001;
	camera.far = 10.0;
	camera.fovy = 30.0 * 3.14/180.0f;
	camera.aspect_ratio = 16.0/9.0;

	camera.updateViewMatrix();
	camera.updateProjectionMatrix();

	/* Make camera accessable in window callbacks */
	glfwSetWindowUserPointer(window,&camera);

	/* Create the debug sphere */
	DebugSphere db_sphere;
	db_sphere.create();


	//////////////////////////////////////////////////////
	// Render Loop - Each loop iteration renders one frame
	//////////////////////////////////////////////////////

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

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
