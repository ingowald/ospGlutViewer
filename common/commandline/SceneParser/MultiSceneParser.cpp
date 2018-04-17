// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "MultiSceneParser.h"

#include "particle/ParticleSceneParser.h"
#include "streamlines/StreamLineSceneParser.h"
#ifdef OSPRAY_TACHYON_SUPPORT
#  include "tachyon/TachyonSceneParser.h"
#endif
#include "trianglemesh/TriangleMeshSceneParser.h"
#ifndef _WIN32
// # include "volume/VolumeSceneParser.h"
#endif

#include "biff/biff.h"
// ptex:
#include "ptex/Ptexture.h"
// #include "ospcommon/tasking/parallel_for.h"

using namespace ospray;
using namespace ospcommon;

template<typename T>
void microBench(const char *name, const T &lambda) 
{
  static size_t numCalls = 0;
  static size_t nextPing = 16;
  static double totalTime = 0.f;

  double before = getSysTime();
  lambda();
  double after = getSysTime();
  totalTime += (after-before);
  numCalls++;
  if (numCalls == nextPing) {
    std::cout << "avg time for '" << name << "' is " << (totalTime / numCalls) << std::endl;
    nextPing *= 2;
  }
}

// std::map<const biff::Texture *, void *> filterForTexture;

// void *getPTexFilterFor(const biff::Texture *btex)
// {
//   if (filterForTexture.find(btex) == filterForTexture.end()) {
//     uint64_t base = (uint64_t)&btex->rawData[0];
//     uint64_t size = (uint64_t)btex->rawData.size();
//     char ptexName[200];
//     sprintf(ptexName,"mem://%lx/%lx",base,size);
    
//     Ptex::String error;
//     PtexTexture *ptex_texture = PtexTexture::open(ptexName,error);

//     PtexFilter::Options opts;
//     opts.filter = PtexFilter::f_point; //f_bilinear;
//     opts.lerp = 0;
//     opts.sharpness = 1.f;
//     opts.noedgeblend = false;
//     PtexFilter *ptex_filter = PtexFilter::getFilter(ptex_texture,opts);

//     filterForTexture[btex] = (void*)ptex_filter;
//   }
//   return filterForTexture[btex];
// }



struct BiffParser : public SceneParser {
  ospray::cpp::Model rootModel;
  ospcommon::box3f   rootBounds;
  ospray::cpp::Renderer ospRenderer;

  BiffParser(ospray::cpp::Renderer ospRenderer) : ospRenderer(ospRenderer) {}

  virtual std::deque<ospray::cpp::Model> model() const 
  {
    std::deque<ospray::cpp::Model> d; d.push_back(rootModel); return d;
  }
  virtual std::deque<ospcommon::box3f>   bbox()  const 
  {
    std::deque<ospcommon::box3f> d; d.push_back(rootBounds); return d;
  }


  OSPData createOspDataFromTexture(std::shared_ptr<biff::Texture> btex)
  {
    static std::map<std::shared_ptr<biff::Texture>,OSPData> alreadyCreated;
    if (alreadyCreated.find(btex) == alreadyCreated.end()) {
      alreadyCreated[btex] = ospNewData(btex->rawData.size(),
                                        OSP_RAW,
                                        btex->rawData.data(),
                                        OSP_DATA_SHARED_BUFFER);
    }
    return alreadyCreated[btex];
  }

  cpp::Material makeMaterial(std::shared_ptr<biff::Material> material,
                             std::shared_ptr<biff::Scene> scene)
  {
    cpp::Material mat = ospray::cpp::Material("scivis","OBJMaterial");
    mat.set("Kd", .6f, 0.6f, 0.6f);

    for (auto p : material->param_string) 
      ospSetString(mat.handle(),p.first.c_str(),p.second.c_str());
    for (auto p : material->param_int) 
      ospSet1i(mat.handle(),p.first.c_str(),p.second);
    for (auto p : material->param_float) 
      ospSet1f(mat.handle(),p.first.c_str(),p.second);
    for (auto p : material->param_vec3f) 
      ospSet3f(mat.handle(),p.first.c_str(),p.second.x,p.second.y,p.second.z);
    for (auto p : material->param_texture) {
      std::shared_ptr<biff::Texture> btex = scene->getTexture(p.second);
      if (!btex->rawDataSize)
        throw std::runtime_error("can currently do ONLY embedded ptex");
      // ospSetString(mat.handle(),p.first.c_str(),btex->param_string["filename"].c_str());
      ospSetData(mat.handle(),
                 ("ptex://"+p.first).c_str(),
                 createOspDataFromTexture(btex));
      // PRINT(p.first);
      // PRINT(btex->param_string["filename"]);
    }

    mat.commit();
    return mat;
  }

  cpp::Geometry buildTriMesh(std::shared_ptr<biff::Scene> scene, 
                          biff::GeomHandle geom) 
  {
    cpp::Geometry g("triangles");
    std::shared_ptr<biff::TriMesh> mesh 
      = scene->triMeshes[geom.geomID];
    std::shared_ptr<biff::Material> bmat
      = scene->getMaterial(mesh->materialID);
    cpp::Material mat = makeMaterial(bmat// ,btex
                                     ,scene
                                     );
    g.setMaterial(mat);
      
    // padding!
    mesh->vtx.push_back(vec3f(0));
    OSPData vtx = ospNewData(mesh->vtx.size()-1,OSP_FLOAT3,&mesh->vtx[0],OSP_DATA_SHARED_BUFFER);
    g.set("position", vtx);
    OSPData idx = ospNewData(mesh->idx.size(),OSP_INT3,&mesh->idx[0],OSP_DATA_SHARED_BUFFER);
    g.set("index", idx);
    if (mesh->txt.size()) {
      OSPData txt = ospNewData(mesh->txt.size(),OSP_FLOAT2,&mesh->txt[0],OSP_DATA_SHARED_BUFFER);
      g.set("texcoord", txt);
    }
    g.commit();
    return g;
  }

  cpp::Geometry buildQuadMesh(std::shared_ptr<biff::Scene> scene, 
                              biff::GeomHandle geom) 
  {
    // for now, let's split everything into tris - fix as soon as we
    // have 'real' quads in ospray
    cpp::Geometry g("triangles");
    std::shared_ptr<biff::QuadMesh> mesh 
      = scene->quadMeshes[geom.geomID];
    std::shared_ptr<biff::Material> bmat
      = scene->getMaterial(mesh->materialID);
    cpp::Material mat = makeMaterial(bmat// ,btex
                                     ,scene
                                     );
    g.setMaterial(mat);
      
    // padding!
    mesh->vtx.push_back(vec3f(0));
    OSPData vtx = ospNewData(mesh->vtx.size()-1,OSP_FLOAT3,&mesh->vtx[0],OSP_DATA_SHARED_BUFFER);
    g.set("position", vtx);
#if 1
    std::vector<vec3i> asTris;
    for (vec4i quad : mesh->idx) {
      asTris.push_back(vec3i(quad.x,quad.y,quad.z));
      asTris.push_back(vec3i(quad.x,quad.z,quad.w));
    }
    OSPData idx = ospNewData(asTris.size(),OSP_INT3,&asTris[0],0);
#else
    OSPData idx = ospNewData(mesh->idx.size(),OSP_INT3,&mesh->idx[0],OSP_DATA_SHARED_BUFFER);
#endif
    g.set("index", idx);
    if (mesh->txt.size()) {
      OSPData txt = ospNewData(mesh->txt.size(),OSP_FLOAT2,&mesh->txt[0],OSP_DATA_SHARED_BUFFER);
      g.set("texcoord", txt);
    }
    g.commit();
    return g;
  }

  cpp::Geometry buildGeometry(std::shared_ptr<biff::Scene> scene, 
                          biff::GeomHandle geom) 
  {
    switch (geom.type) {
    case biff::TRI_MESH:
      return buildTriMesh(scene,geom);
    case biff::QUAD_MESH:
      return buildQuadMesh(scene,geom);
    default: NOTIMPLEMENTED;
    }
  }


  bool parse(int ac, const char **&av)
  {
    for (int i=0;i<ac;i++) {
      std::string arg = av[i];
      if (arg.substr(arg.size()-5) == ".biff"
          ||
          arg.substr(arg.size()-6) == ".biff/") {
        // HACKS ALL OVER THE PLACE - JUST FOR TESTING!!!!!
        std::cout << "Loading BIFF model" << std::endl;
        static std::shared_ptr<biff::Scene> scene = biff::Scene::read(av[1]);
        std::cout << "done loading input..." << std::endl;
    
        rootBounds = ospcommon::box3f(vec3f(-100),vec3f(+100));
    
        // first - create one model (with shapes) per biff group
        std::vector<cpp::Model> groupModels;
        std::cout << "creating " << scene->groups.size() << " groups..." << std::endl;
        for (auto &group : scene->groups) {
          cpp::Model model;
          std::cout << " - " << group->geometries.size() << " geometries              " << std::endl;
          for (int shapeID=0;shapeID<(int)group->geometries.size();shapeID++) {
            std::cout << "\r[" << shapeID << "/" << group->geometries.size() << "]    ";
            try {
              cpp::Geometry g = buildGeometry(scene,group->geometries[shapeID]);
              model.addGeometry(g);
            } catch (...) {
              static bool warned = false;
              if (!warned) {
                std::cout << "WARNING - AT LEAST ONE UNSUPPORTED GEOMETRY" << std::endl;
                warned = true;
              }
            }
          }
          std::cout << "\r";
          model.commit();
          groupModels.push_back(model);
        }

        // second: create the root model, by instantiating all the groups
        std::cout << "creating max " << scene->instances.size() << " instances" << std::endl;
        for (auto &inst : scene->instances) {
          OSPGeometry ginst =
            ospNewInstance(groupModels[inst.groupID].handle(),
                           (osp::affine3f&)inst.xfm);
          rootModel.addGeometry(ginst);
        }
        std::cout << std::endl;
        std::cout << "done building instances - committing toplevel scene..." << std::endl;
        rootModel.commit();
    
        return true;
      }
    }
  return false;
  }
};








MultiSceneParser::MultiSceneParser(cpp::Renderer renderer) 
  : renderer(renderer)
{
}

bool MultiSceneParser::parse(int ac, const char **&av)
{
  BiffParser biffParser(renderer);
  TriangleMeshSceneParser triangleMeshParser(renderer);
#ifdef OSPRAY_TACHYON_SUPPORT
  TachyonSceneParser      tachyonParser(renderer);
#endif
  ParticleSceneParser     particleParser(renderer);
  StreamLineSceneParser   streamlineParser(renderer);
#ifndef _WIN32
  // VolumeSceneParser       volumeParser(renderer);
#endif

  bool gotBiff = biffParser.parse(ac, av);
  bool gotTriangleMeshScene = triangleMeshParser.parse(ac, av);
#ifdef OSPRAY_TACHYON_SUPPORT
  bool gotTachyonScene      = tachyonParser.parse(ac, av);
#endif
  bool gotParticleScene      = particleParser.parse(ac, av);
  bool gotStreamLineScene   = streamlineParser.parse(ac, av);
#ifndef _WIN32
  // bool gotVolumeScene       = volumeParser.parse(ac, av);
#endif

  SceneParser *parser = nullptr;

  if (gotBiff)
    parser = &biffParser;
  else if (gotTriangleMeshScene)
    parser = &triangleMeshParser;
#ifdef OSPRAY_TACHYON_SUPPORT
  else if (gotTachyonScene)
    parser = &tachyonParser;
#endif
  else if (gotParticleScene)
    parser = &particleParser;
  else if (gotStreamLineScene)
    parser = &streamlineParser;
#ifndef _WIN32
  // else if (gotVolumeScene)
  //   parser = &volumeParser;
#endif

  if (parser) {
    sceneModels = parser->model();
    sceneBboxes = parser->bbox();
  } else {
    sceneModels.push_back(cpp::Model());
    sceneModels[0].commit();
  }

  return parser != nullptr;
}

std::deque<cpp::Model> MultiSceneParser::model() const
{
  return sceneModels;
}

std::deque<box3f> MultiSceneParser::bbox() const
{
  return sceneBboxes;
}
