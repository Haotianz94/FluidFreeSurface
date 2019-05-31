#include "quadmesh.h"
#include "Configer.h"
#include "Logger.h"
#include <fstream>
#include <memory.h>


QuadMesh::QuadMesh()
{	
	num_vertex_ = 0;
	num_vertex_uv_ = 0;
	num_normal_ = 0;
	num_face_ = 0;
}

QuadMesh::QuadMesh(std::string file_path)
{
	num_vertex_ = 0;
	num_vertex_uv_ = 0;
	num_normal_ = 0;
	num_face_ = 0;
	file_path_ = file_path;
	loadObj(file_path_);
}

void QuadMesh::loadObj(std::string file_path)
{

	std::ifstream fin(file_path);
	assert(!fin.fail() && "Can not open the obj file!");

	char c1, c2;
	std::string s;
	float x, y, z;
	int ids[12];
	while(fin >> c1)
	{
		if(c1 == 'f')
		{
			fin >> ids[0] >> c1 >> ids[4] >> c2 >> ids[8];
			fin >> ids[1] >> c1 >> ids[5] >> c2 >> ids[9];
			fin >> ids[2] >> c1 >> ids[6] >> c2 >> ids[10];
			fin >> ids[3] >> c1 >> ids[7] >> c2 >> ids[11];
			for(int i = 0; i < 12; i++)
				ids[i] -= 1;
			faces_.push_back(QuadIdx(ids));
			num_face_ ++;
		}
		else if(c1 == 'v')
		{
			fin.get(c2);
			fin >> x >> y >> z;
			
			if(c2 == ' ')
			{
				vertices_.push_back(Eigen::Vector3f(x, y, z));
				num_vertex_ ++;
			}
			else if(c2 == 't')
			{
				vertices_uv_.push_back(Eigen::Vector3f(x, y, z));
				num_vertex_uv_ ++;
			}
			else if(c2 == 'n')
			{
				normals_.push_back(Eigen::Vector3f(x, y, z));
				num_normal_ ++;
			}
		}
		else //(c1 == '#')
		{	
			getline(fin, s);
			continue;
		}
	}
	fin.close();

	LOG << "Load " << file_path.c_str() << "\n";
	REPORT(num_vertex_);
	REPORT(num_vertex_uv_);
	REPORT(num_normal_);
	REPORT(num_face_);
}

void QuadMesh::dumpObj(std::string obj_path)
{
	std::ofstream fout(obj_path);

	fout << "mtllib volcano_01.mtl\n";

	// Dump all vertices
	fout << "# vertices " << num_vertex_ << "\n";
	for(auto& v : vertices_)
		fout << "v " << v[0] << ' ' << v[1] << ' ' << v[2] << "\n";

	// Dump UV vertices
	fout << "# uv vertices " << num_vertex_uv_ << "\n";
	for(auto& vt : vertices_uv_)
		fout << "vt " << vt[0] << ' ' << vt[1] << ' ' << vt[2] << "\n";
	
	
	// Dump normal
	fout << "# normals " << num_normal_ << "\n";
	for(auto& vn : normals_)
		fout << "vn " << vn[0] << ' ' << vn[1] << ' ' << vn[2] << "\n";

	fout << "g VOLCANO\n";
	fout << "s 1\n";
	fout << "usemtl Material__32\n";

	// Dump all faces_
	fout << "# faces_ " << num_face_ << "\n";
	
	for(auto& f : faces_)
	{
		fout << "f";
		for(int j = 0; j < 4; j++)
			fout << ' ' << f.pt_id[j]+1 << '/' << f.vt_id[j]+1 << '/' << f.vn_id[j]+1;
		fout << "\n";
	}
	fout.close();

	LOG << "Dump obj file "<< obj_path.c_str() << "\n";
}

QuadMesh::~QuadMesh()
{
	delete[] height_map_;
}

void QuadMesh::normalize(int max_length)
{
	max_length_ = max_length;
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;
	for(auto& v : vertices_)
	{
		min_x = std::min(min_x, v[0]);
		max_x = std::max(max_x, v[0]);
		min_y = std::min(min_y, v[1]);
		max_y = std::max(max_y, v[1]);
		min_z = std::min(min_z, v[2]);
		max_z = std::max(max_z, v[2]);
		REPORT(v[0]);
	}
	REPORT(min_x);
	float scale = 1.0f * max_length / std::max(std::max(max_x - min_x, max_y - min_y), max_z - min_z);
	Eigen::Vector3f shift(min_x * scale, min_y * scale*3, min_z * scale); 
	REPORT(scale);
	std::cout << shift << std::endl;
	for(auto& v : vertices_)
	{
		v = Eigen::Vector3f(v[0] * scale, v[1] * scale*3, v[2] * scale) - shift;
	}
}

void QuadMesh::calculateHeightMap()
{
	height_map_ = new int[max_length_ * max_length_];
	memset(height_map_, 0, sizeof(int) * max_length_ * max_length_);
	for(auto& v : vertices_)
	{
		int x = int(v[0]);
		int y = int(v[1]);
		int z = int(v[2]);
		height_map_[z * max_length_ + x] = std::max(height_map_[z * max_length_ + x], y);
	}
}

int QuadMesh::getHeight(int x, int z)
{
	return height_map_[z * max_length_ + x];
}