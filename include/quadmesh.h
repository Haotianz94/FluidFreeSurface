#ifndef _QUADMESH_H_
#define _QUADMESH_H_

#include <Eigen/Eigen>
#include <cstring>

struct QuadIdx
{
	int pt_id[4];
	int vt_id[4];
	int vn_id[4];
	QuadIdx(int* ids)
	{
		for(int i = 0; i < 4; i++)
			pt_id[i] = ids[i];
		for(int i = 0; i < 4; i++)
			vt_id[i] = ids[i+4];
		for(int i = 0; i < 4; i++)
			vn_id[i] = ids[i+8];
	}
};

class QuadMesh
{
public:
	std::string file_path_;
	int num_vertex_;
	int num_vertex_uv_;
	int num_normal_;
	int num_face_;

	// Volcano
	int max_length_;
	int *height_map_;

	std::vector<Eigen::Vector3f> vertices_uv_; 
	std::vector<Eigen::Vector3f> vertices_;
	std::vector<Eigen::Vector3f> normals_;
	std::vector<QuadIdx> faces_;

	QuadMesh();
	QuadMesh(std::string);
	~QuadMesh();

	void loadObj(std::string);
	void dumpObj(std::string);
	void normalize(int max_length);
	void calculateHeightMap(int max_length);
	int getHeight(int x, int z);
};

#endif