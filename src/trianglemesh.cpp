#include "trianglemesh.h"
#include "Configer.h"
#include "Logger.h"
#include <fstream>
#include <memory.h>
#include <math.h>


TriangleMesh::TriangleMesh()
{	
	num_vertex_ = 0;
	num_vertex_uv_ = 0;
	num_normal_ = 0;
	num_face_ = 0;
}

TriangleMesh::TriangleMesh(std::string file_path)
{
	num_vertex_ = 0;
	num_vertex_uv_ = 0;
	num_normal_ = 0;
	num_face_ = 0;
	file_path_ = file_path;
	loadObj(file_path_);
}

void TriangleMesh::loadObj(std::string file_path)
{

	std::ifstream fin(file_path);
	assert(!fin.fail() && "Can not open the obj file!");

	char c1, c2;
	std::string s;
	float x, y, z;
	int ids[9];
	while(fin >> c1)
	{
		if(c1 == 'f')
		{
			//fin >> ids[0] >> c1 >> ids[3] >> c2 >> ids[6];
			//fin >> ids[1] >> c1 >> ids[4] >> c2 >> ids[7];
			//fin >> ids[2] >> c1 >> ids[5] >> c2 >> ids[8];
			 fin >> ids[0] >> c1 >> c2 >> ids[6];
			 fin >> ids[1] >> c1 >> c2 >> ids[7];
			 fin >> ids[2] >> c1 >> c2 >> ids[8];
			for(int i = 0; i < 9; i++)
				ids[i] -= 1;
			faces_.push_back(TriangleIdx(ids));
			num_face_ ++;
		}
		else if(c1 == 'v')
		{
			fin.get(c2);			
			if(c2 == ' ')
			{
				fin >> x >> y >> z;
				vertices_.push_back(Eigen::Vector3f(x, y, z));
				num_vertex_ ++;
			}
			else if(c2 == 't')
			{
				fin >> x >> y;
				vertices_uv_.push_back(Eigen::Vector3f(x, y, 0));
				num_vertex_uv_ ++;
			}
			else if(c2 == 'n')
			{
				fin >> x >> y >> z;
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

void TriangleMesh::dumpObj(std::string obj_path)
{
	std::ofstream fout(obj_path);

	fout << "mtllib volcano.mtl\n";

	// Dump all vertices
	fout << "# vertices " << num_vertex_ << "\n";
	for(auto& v : vertices_)
		fout << "v " << v[0] << ' ' << v[1] << ' ' << v[2] << "\n";

	// Dump UV vertices
	// fout << "# uv vertices " << num_vertex_uv_ << "\n";
	// for(auto& vt : vertices_uv_)
	// 	fout << "vt " << vt[0] << ' ' << vt[1] << ' ' << vt[2] << "\n";
	
	
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
		for(int j = 0; j < 3; j++)
			// fout << ' ' << f.pt_id[j]+1 << '/' << f.vt_id[j]+1 << '/' << f.vn_id[j]+1;
			fout << ' ' << f.pt_id[j]+1 << '/' << '/' << f.vn_id[j]+1;
		fout << "\n";
	}
	fout.close();

	LOG << "Dump obj file "<< obj_path.c_str() << "\n";
}

TriangleMesh::~TriangleMesh()
{
}

void TriangleMesh::normalize(int max_length)
{
	max_length_ = max_length;
	float min_x = 1e9, min_y = 1e9, min_z = 1e9;
	float max_x = -1e9, max_y = -1e9, max_z = -1e9;
	for(auto& v : vertices_)
	{
		min_x = std::min(min_x, v[0]);
		max_x = std::max(max_x, v[0]);
		min_y = std::min(min_y, v[1]);
		max_y = std::max(max_y, v[1]);
		min_z = std::min(min_z, v[2]);
		max_z = std::max(max_z, v[2]);
	}
	float scale = 1.0f * max_length / std::max(std::max(max_x - min_x, max_y - min_y), max_z - min_z);
	Eigen::Vector3f shift(min_x * scale, min_y * scale, min_z * scale); 
	REPORT(scale);
	std::cout << shift << std::endl;
	for(auto& v : vertices_)
	{
		v = Eigen::Vector3f(v[0] * scale, v[1] * scale, v[2] * scale) - shift;
	}
}

void TriangleMesh::refine()
{
	int n = 0;
	std::vector<TriangleIdx> faces_refine;
	for(auto& f : faces_)
	{
   	 	// it.n = Normalize(Normal3f(Cross(p1 - p0, p2 - p0)));
		Eigen::Vector3f u = vertices_[f.pt_id[1]] - vertices_[f.pt_id[0]]; 
		Eigen::Vector3f v = vertices_[f.pt_id[2]] - vertices_[f.pt_id[0]]; 
		Eigen::Vector3f cross = u.cross(v);
		if(abs(cross.dot(cross)) < 1e-3)
		{
			// LOG << cross[0] << " " << cross[1] << " " << cross[2] << '\n';
			n += 1;
		}
		else
		{
			faces_refine.push_back(f);
		}
	}
	REPORT(n);
	faces_ = faces_refine;
}
