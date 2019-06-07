#ifndef _TRIANGLEMESH_H_
#define _TRIANGLEMESH_H_

#include <Eigen/Eigen>
#include <cstring>

struct TriangleIdx
{
	int pt_id[3];
	int vt_id[3];
	int vn_id[3];
	TriangleIdx(int* ids)
	{
		for(int i = 0; i < 3; i++)
			pt_id[i] = ids[i];
		for(int i = 0; i < 3; i++)
			vt_id[i] = ids[i+3];
		for(int i = 0; i < 3; i++)
			vn_id[i] = ids[i+6];
	}
};

class TriangleMesh
{
public:
	std::string file_path_;
	int num_vertex_;
	int num_vertex_uv_;
	int num_normal_;
	int num_face_;

	int max_length_;

	std::vector<Eigen::Vector3f> vertices_uv_; 
	std::vector<Eigen::Vector3f> vertices_;
	std::vector<Eigen::Vector3f> normals_;
	std::vector<TriangleIdx> faces_;

	TriangleMesh();
	TriangleMesh(std::string);
	~TriangleMesh();

	void loadObj(std::string);
	void dumpObj(std::string);
	void normalize(int max_length);
	void refine();
};

#endif