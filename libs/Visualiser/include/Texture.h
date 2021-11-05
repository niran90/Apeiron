#pragma once

#include <GL/glew.h>
#include "../../../include/Global.h"
#include "GLDebug.h"
#include "GLTypes.h"

namespace Apeiron {

class Texture
{
private:
  UInt RendererID;
  std::string FilePath;
  UChar* LocalBuffer;
  int Width, Height;
  int BitsPerPixel;

public:
  Texture(const std::string& _file_path);
  ~Texture();

  void Bind(UInt _slot = 0) const;
  void Unbind() const;
};

}