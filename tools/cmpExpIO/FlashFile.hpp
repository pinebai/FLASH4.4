#ifndef FLASHFILE_H
#define FLASHFILE_H

#include <string>
#include <vector>

class FlashFile {
public:
  virtual std::vector<std::string> GetAllVariableNames() const = 0;
  virtual size_t GetNumberDataElements() const = 0;
  virtual const double * GetVariableFromFile(const std::string &varName) const = 0;
  virtual ~FlashFile() {};
};

#endif
