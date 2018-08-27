#include <stdio.h>
#include <stdlib.h>

void load_bmp(const char* filename, std::vector<unsigned char> &tex, int &W, int& H) {
  FILE* f;
  f = fopen(filename, "rb");
  unsigned char info[54];
  fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
  
  W = *(int*)&info[18]; // extract image height and width from header
  H = *(int*)&info[22];
  
  int size = 3 * W * H;
  tex.resize(size); // allocate 3 bytes per pixel
  fread(&tex[0], sizeof(unsigned char), size, f); // read the rest of the data at once
  fclose(f);
  
  for (int i = 0; i < size; i += 3) {
    std::swap(tex[i], tex[i + 2]);
  }
}

